use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles::sam;

const MAGIC_NUMBER: &[u8; 4] = b"BAM\x01";

/// Read a BAM header tolerantly: same shape as `noodles::bam::io::Reader::read_header`,
/// but sanitizes any `@PG` records missing the `VN:` tag (which noodles 0.83 rejects
/// even though the SAM spec marks `VN:` as recommended, not required).
///
/// `inner` is the *decompressed* BAM byte stream (typically a `&mut bgzf::Reader<File>`).
/// The function consumes magic + `l_text` + raw header text. It does not read the
/// binary reference section: when `@SQ` lines are present in the SAM text (the case
/// for every real-world BAM produced by samtools/bowtie2/etc.), reference sequences
/// come from there. Callers that need to walk records still use `bam::io::IndexedReader`
/// on the same file — only header parsing is bypassed.
pub fn read_header_lenient<R: Read>(reader: &mut R) -> io::Result<sam::Header> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC_NUMBER {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a BAM file (bad magic)",
        ));
    }

    let l_text = reader.read_u32::<LittleEndian>()? as usize;
    let mut text_bytes = vec![0u8; l_text];
    reader.read_exact(&mut text_bytes)?;

    // SAM header text is NUL-padded inside the l_text region. Trim at the first NUL.
    let end = text_bytes
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(text_bytes.len());
    let text = std::str::from_utf8(&text_bytes[..end])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let sanitized = sanitize_sam_header_text(text);
    sanitized
        .parse::<sam::Header>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Sanitize SAM header text by inserting `VN:unknown` into any `@PG` record that
/// lacks a `VN:` tag. All other lines pass through untouched — including malformed
/// ones, so unrelated header errors still surface from the noodles parser.
fn sanitize_sam_header_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if !line.starts_with("@PG\t") {
            out.push_str(line);
            continue;
        }

        let nl_len = if line.ends_with("\r\n") {
            2
        } else if line.ends_with('\n') {
            1
        } else {
            0
        };
        let body = &line[..line.len() - nl_len];
        let nl = &line[line.len() - nl_len..];

        let has_vn = body.split('\t').any(|f| f.starts_with("VN:"));
        if has_vn {
            out.push_str(line);
        } else {
            out.push_str(body);
            out.push_str("\tVN:unknown");
            out.push_str(nl);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_header_bytes(text: &[u8]) -> Vec<u8> {
        let mut data = Vec::new();
        data.extend_from_slice(MAGIC_NUMBER);
        data.extend_from_slice(&(text.len() as u32).to_le_bytes());
        data.extend_from_slice(text);
        // n_ref = 0; helper does not read the binary reference section, but
        // round-tripping through a real reader would require it.
        data.extend_from_slice(&0u32.to_le_bytes());
        data
    }

    #[test]
    fn sanitize_inserts_vn_when_missing() {
        let s = "@PG\tID:bowtie2\tPN:bowtie2\n";
        assert_eq!(
            sanitize_sam_header_text(s),
            "@PG\tID:bowtie2\tPN:bowtie2\tVN:unknown\n"
        );
    }

    #[test]
    fn sanitize_leaves_pg_with_vn() {
        let s = "@PG\tID:x\tVN:1.0\n";
        assert_eq!(sanitize_sam_header_text(s), s);
    }

    #[test]
    fn sanitize_leaves_non_pg_lines() {
        let s = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
        assert_eq!(sanitize_sam_header_text(s), s);
    }

    #[test]
    fn sanitize_handles_empty() {
        assert_eq!(sanitize_sam_header_text(""), "");
    }

    #[test]
    fn sanitize_handles_mixed_pg_records() {
        let s = "@HD\tVN:1.6\n@PG\tID:a\tPN:a\n@PG\tID:b\tVN:2.0\n";
        let expected = "@HD\tVN:1.6\n@PG\tID:a\tPN:a\tVN:unknown\n@PG\tID:b\tVN:2.0\n";
        assert_eq!(sanitize_sam_header_text(s), expected);
    }

    #[test]
    fn sanitize_preserves_crlf() {
        let s = "@PG\tID:x\tPN:x\r\n";
        assert_eq!(sanitize_sam_header_text(s), "@PG\tID:x\tPN:x\tVN:unknown\r\n");
    }

    #[test]
    fn sanitize_leaves_pg_substring_safely() {
        // A field containing the substring "VN:" but not at start should still trigger
        // VN tag detection only when the field actually starts with "VN:".
        let s = "@PG\tID:x\tPN:somethingVN:weird\n";
        assert_eq!(
            sanitize_sam_header_text(s),
            "@PG\tID:x\tPN:somethingVN:weird\tVN:unknown\n"
        );
    }

    #[test]
    fn read_header_accepts_pg_without_vn() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:bowtie2\tPN:bowtie2\n";
        let bytes = build_header_bytes(text);
        let mut reader = &bytes[..];
        let header = read_header_lenient(&mut reader)?;
        assert_eq!(header.reference_sequences().len(), 1);
        assert!(
            !header.programs().as_ref().is_empty(),
            "expected @PG record to be parsed"
        );
        Ok(())
    }

    #[test]
    fn read_header_round_trips_valid_pg() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:x\tVN:1.0\n";
        let bytes = build_header_bytes(text);
        let mut reader = &bytes[..];
        let header = read_header_lenient(&mut reader)?;
        assert_eq!(header.reference_sequences().len(), 1);
        assert!(!header.programs().as_ref().is_empty());
        Ok(())
    }

    #[test]
    fn read_header_rejects_bad_magic() {
        let mut data: &[u8] = b"NOTB";
        let err = read_header_lenient(&mut data).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn read_header_errors_on_truncated_text() {
        // l_text declares 16 bytes but only 2 are present after the prefix.
        let mut data: Vec<u8> = Vec::new();
        data.extend_from_slice(MAGIC_NUMBER);
        data.extend_from_slice(&16u32.to_le_bytes());
        data.extend_from_slice(b"@H");
        let mut r = &data[..];
        assert!(read_header_lenient(&mut r).is_err());
    }

    #[test]
    fn read_header_two_pg_records_one_missing_vn() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n\
                     @PG\tID:a\tPN:a\n\
                     @PG\tID:b\tPN:b\tVN:2.0\tPP:a\n";
        let bytes = build_header_bytes(text);
        let mut reader = &bytes[..];
        let header = read_header_lenient(&mut reader)?;
        // Both programs present.
        let programs = header.programs();
        assert_eq!(programs.as_ref().len(), 2);
        Ok(())
    }
}
