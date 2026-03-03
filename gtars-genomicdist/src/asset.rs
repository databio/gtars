//! GenomicDist Annotation (GDA) — binary gene model format.
//!
//! Replaces the old GMDL format. Stores gene model components (genes, exons,
//! optionally UTRs) with strand information in a compact column-oriented layout.
//! ChromSizes are served separately (e.g. via refgenie).

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use gtars_core::models::{Region, RegionSet};

use crate::errors::GtarsGenomicDistError;
use crate::models::{Strand, StrandedRegionSet};
use crate::partitions::GeneModel;

/// Magic bytes: "GDA\0" = 0x00414447 (little-endian).
const GDA_MAGIC: u32 = 0x00414447;
/// Current GDA format version.
const GDA_VERSION: u32 = 1;

/// Flag bit: three_utr component present.
const FLAG_HAS_THREE_UTR: u32 = 1;
/// Flag bit: five_utr component present.
const FLAG_HAS_FIVE_UTR: u32 = 2;

/// A genomic distribution annotation: serializable gene model.
///
/// Wraps a `GeneModel` with binary I/O in the GDA format. ChromSizes are
/// not included — they come from a separate source (refgenie, .chrom.sizes
/// file, etc.).
pub struct GenomicDistAnnotation {
    pub gene_model: GeneModel,
}

impl GenomicDistAnnotation {
    /// Construct from a GTF file.
    pub fn from_gtf<P: AsRef<Path>>(gtf_path: P) -> Result<Self, GtarsGenomicDistError> {
        let gene_model = GeneModel::from_gtf(
            gtf_path.as_ref().to_str().ok_or_else(|| {
                GtarsGenomicDistError::CustomError("Invalid GTF path".to_string())
            })?,
            true,
            true,
        )?;
        Ok(GenomicDistAnnotation { gene_model })
    }

    /// Serialize to the GDA binary format and write to a file.
    pub fn save_bin<P: AsRef<Path>>(&self, path: P) -> Result<(), GtarsGenomicDistError> {
        let buf = self.to_bytes();
        let mut file = File::create(path)?;
        file.write_all(&buf)?;
        Ok(())
    }

    /// Deserialize from a GDA binary file.
    pub fn load_bin<P: AsRef<Path>>(path: P) -> Result<Self, GtarsGenomicDistError> {
        let data = std::fs::read(path.as_ref())?;
        Self::load_bin_from_bytes(&data)
    }

    /// Deserialize from a byte slice (for WASM or in-memory use).
    pub fn load_bin_from_bytes(data: &[u8]) -> Result<Self, GtarsGenomicDistError> {
        let err = |msg: &str| GtarsGenomicDistError::CustomError(msg.to_string());
        let mut pos = 0usize;

        macro_rules! check {
            ($n:expr) => {
                if pos + $n > data.len() {
                    return Err(err("Unexpected end of GDA file"));
                }
            };
        }
        macro_rules! read_u32 {
            () => {{
                check!(4);
                let v = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
                pos += 4;
                v
            }};
        }

        // ── Header (16 bytes) ──────────────────────────────────────────
        let magic = read_u32!();
        if magic != GDA_MAGIC {
            return Err(err(
                "Invalid GDA file format — regenerate with 'gtars prep --gtf ...'",
            ));
        }
        let version = read_u32!();
        if version != GDA_VERSION {
            return Err(err(&format!("Unsupported GDA format version {}", version)));
        }
        let _n_components = read_u32!();
        let flags = read_u32!();
        let has_three_utr = flags & FLAG_HAS_THREE_UTR != 0;
        let has_five_utr = flags & FLAG_HAS_FIVE_UTR != 0;

        // ── String intern table ────────────────────────────────────────
        let n_strings = read_u32!() as usize;
        let mut intern_table: Vec<String> = Vec::with_capacity(n_strings);
        for _ in 0..n_strings {
            let len = read_u32!() as usize;
            check!(len);
            let s = String::from_utf8(data[pos..pos + len].to_vec())
                .map_err(|e| err(&format!("Invalid UTF-8 in intern table: {}", e)))?;
            pos += len;
            intern_table.push(s);
        }

        // ── Gene model components ──────────────────────────────────────
        let read_srs =
            |data: &[u8],
             pos: &mut usize,
             intern_table: &[String]|
             -> Result<StrandedRegionSet, GtarsGenomicDistError> {
                if *pos + 4 > data.len() {
                    return Err(err("Unexpected end of file"));
                }
                let n = u32::from_le_bytes(data[*pos..*pos + 4].try_into().unwrap()) as usize;
                *pos += 4;

                let mut chr_ids = Vec::with_capacity(n);
                for _ in 0..n {
                    if *pos + 2 > data.len() {
                        return Err(err("Unexpected end of file"));
                    }
                    chr_ids
                        .push(u16::from_le_bytes(data[*pos..*pos + 2].try_into().unwrap()) as usize);
                    *pos += 2;
                }

                let starts_size = n * 4;
                if *pos + starts_size > data.len() {
                    return Err(err("Unexpected end of file"));
                }
                let starts_bytes = &data[*pos..*pos + starts_size];
                *pos += starts_size;

                let ends_size = n * 4;
                if *pos + ends_size > data.len() {
                    return Err(err("Unexpected end of file"));
                }
                let ends_bytes = &data[*pos..*pos + ends_size];
                *pos += ends_size;

                if *pos + n > data.len() {
                    return Err(err("Unexpected end of file"));
                }
                let strand_bytes = &data[*pos..*pos + n];
                *pos += n;

                let mut regions = Vec::with_capacity(n);
                let mut strands = Vec::with_capacity(n);
                for i in 0..n {
                    let start = u32::from_le_bytes(
                        starts_bytes[i * 4..(i + 1) * 4].try_into().unwrap(),
                    );
                    let end =
                        u32::from_le_bytes(ends_bytes[i * 4..(i + 1) * 4].try_into().unwrap());
                    let chr = intern_table
                        .get(chr_ids[i])
                        .ok_or_else(|| err("Gene model: intern table index out of bounds"))?
                        .clone();
                    regions.push(Region {
                        chr,
                        start,
                        end,
                        rest: None,
                    });
                    strands.push(match strand_bytes[i] {
                        0 => Strand::Plus,
                        1 => Strand::Minus,
                        _ => Strand::Unstranded,
                    });
                }

                Ok(StrandedRegionSet::new(RegionSet::from(regions), strands))
            };

        let genes = read_srs(data, &mut pos, &intern_table)?;
        let exons = read_srs(data, &mut pos, &intern_table)?;
        let three_utr = if has_three_utr {
            Some(read_srs(data, &mut pos, &intern_table)?)
        } else {
            None
        };
        let five_utr = if has_five_utr {
            Some(read_srs(data, &mut pos, &intern_table)?)
        } else {
            None
        };

        Ok(GenomicDistAnnotation {
            gene_model: GeneModel {
                genes,
                exons,
                three_utr,
                five_utr,
            },
        })
    }

    /// Serialize to GDA bytes.
    fn to_bytes(&self) -> Vec<u8> {
        let mut buf: Vec<u8> = Vec::new();

        // ── Build string intern table ──────────────────────────────────
        let mut intern_map: HashMap<String, u16> = HashMap::new();
        let mut intern_table: Vec<String> = Vec::new();

        let all_components: Vec<&StrandedRegionSet> = [
            Some(&self.gene_model.genes),
            Some(&self.gene_model.exons),
            self.gene_model.three_utr.as_ref(),
            self.gene_model.five_utr.as_ref(),
        ]
        .into_iter()
        .flatten()
        .collect();

        for srs in &all_components {
            for r in &srs.inner.regions {
                if !intern_map.contains_key(&r.chr) {
                    let id = intern_table.len() as u16;
                    intern_map.insert(r.chr.clone(), id);
                    intern_table.push(r.chr.clone());
                }
            }
        }

        // ── Header (16 bytes) ──────────────────────────────────────────
        let mut flags: u32 = 0;
        if self.gene_model.three_utr.is_some() {
            flags |= FLAG_HAS_THREE_UTR;
        }
        if self.gene_model.five_utr.is_some() {
            flags |= FLAG_HAS_FIVE_UTR;
        }
        let n_components: u32 = 2
            + self.gene_model.three_utr.is_some() as u32
            + self.gene_model.five_utr.is_some() as u32;

        buf.extend_from_slice(&GDA_MAGIC.to_le_bytes());
        buf.extend_from_slice(&GDA_VERSION.to_le_bytes());
        buf.extend_from_slice(&n_components.to_le_bytes());
        buf.extend_from_slice(&flags.to_le_bytes());

        // ── String intern table ────────────────────────────────────────
        buf.extend_from_slice(&(intern_table.len() as u32).to_le_bytes());
        for s in &intern_table {
            let bytes = s.as_bytes();
            buf.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            buf.extend_from_slice(bytes);
        }

        // ── Gene model components ──────────────────────────────────────
        let write_srs =
            |srs: &StrandedRegionSet, buf: &mut Vec<u8>, intern_map: &HashMap<String, u16>| {
                let n = srs.inner.regions.len() as u32;
                buf.extend_from_slice(&n.to_le_bytes());
                for r in &srs.inner.regions {
                    buf.extend_from_slice(&intern_map[&r.chr].to_le_bytes());
                }
                for r in &srs.inner.regions {
                    buf.extend_from_slice(&r.start.to_le_bytes());
                }
                for r in &srs.inner.regions {
                    buf.extend_from_slice(&r.end.to_le_bytes());
                }
                for s in &srs.strands {
                    buf.push(match s {
                        Strand::Plus => 0,
                        Strand::Minus => 1,
                        Strand::Unstranded => 2,
                    });
                }
            };

        write_srs(&self.gene_model.genes, &mut buf, &intern_map);
        write_srs(&self.gene_model.exons, &mut buf, &intern_map);
        if let Some(ref srs) = self.gene_model.three_utr {
            write_srs(srs, &mut buf, &intern_map);
        }
        if let Some(ref srs) = self.gene_model.five_utr {
            write_srs(srs, &mut buf, &intern_map);
        }

        buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_annotation(with_utrs: bool) -> GenomicDistAnnotation {
        let genes = StrandedRegionSet::new(
            RegionSet::from(vec![
                Region { chr: "chr1".into(), start: 1000, end: 5000, rest: None },
                Region { chr: "chr1".into(), start: 8000, end: 12000, rest: None },
                Region { chr: "chr2".into(), start: 500, end: 3000, rest: None },
            ]),
            vec![Strand::Plus, Strand::Minus, Strand::Plus],
        );

        let exons = StrandedRegionSet::new(
            RegionSet::from(vec![
                Region { chr: "chr1".into(), start: 1000, end: 1500, rest: None },
                Region { chr: "chr1".into(), start: 3000, end: 3500, rest: None },
                Region { chr: "chr2".into(), start: 500, end: 800, rest: None },
            ]),
            vec![Strand::Plus, Strand::Plus, Strand::Plus],
        );

        let (three_utr, five_utr) = if with_utrs {
            (
                Some(StrandedRegionSet::new(
                    RegionSet::from(vec![
                        Region { chr: "chr1".into(), start: 4500, end: 5000, rest: None },
                    ]),
                    vec![Strand::Plus],
                )),
                Some(StrandedRegionSet::new(
                    RegionSet::from(vec![
                        Region { chr: "chr1".into(), start: 1000, end: 1200, rest: None },
                    ]),
                    vec![Strand::Plus],
                )),
            )
        } else {
            (None, None)
        };

        GenomicDistAnnotation {
            gene_model: GeneModel {
                genes,
                exons,
                three_utr,
                five_utr,
            },
        }
    }

    #[test]
    fn test_gda_round_trip_with_utrs() {
        let ann = make_test_annotation(true);
        let bytes = ann.to_bytes();
        let loaded = GenomicDistAnnotation::load_bin_from_bytes(&bytes).unwrap();

        assert_eq!(
            loaded.gene_model.genes.inner.regions.len(),
            ann.gene_model.genes.inner.regions.len()
        );
        assert_eq!(
            loaded.gene_model.exons.inner.regions.len(),
            ann.gene_model.exons.inner.regions.len()
        );
        assert_eq!(loaded.gene_model.genes.strands, ann.gene_model.genes.strands);

        assert!(loaded.gene_model.three_utr.is_some());
        assert!(loaded.gene_model.five_utr.is_some());

        for (a, b) in ann
            .gene_model
            .genes
            .inner
            .regions
            .iter()
            .zip(loaded.gene_model.genes.inner.regions.iter())
        {
            assert_eq!(a.chr, b.chr);
            assert_eq!(a.start, b.start);
            assert_eq!(a.end, b.end);
        }
    }

    #[test]
    fn test_gda_round_trip_without_utrs() {
        let ann = make_test_annotation(false);
        let bytes = ann.to_bytes();
        let loaded = GenomicDistAnnotation::load_bin_from_bytes(&bytes).unwrap();

        assert_eq!(
            loaded.gene_model.genes.inner.regions.len(),
            ann.gene_model.genes.inner.regions.len()
        );
        assert!(loaded.gene_model.three_utr.is_none());
        assert!(loaded.gene_model.five_utr.is_none());
    }

    #[test]
    fn test_gda_file_round_trip() {
        let ann = make_test_annotation(true);
        let dir = tempfile::tempdir().unwrap();
        let bin_path = dir.path().join("test.gda.bin");

        ann.save_bin(&bin_path).unwrap();
        let loaded = GenomicDistAnnotation::load_bin(&bin_path).unwrap();

        assert_eq!(
            loaded.gene_model.genes.inner.regions.len(),
            ann.gene_model.genes.inner.regions.len()
        );
    }

    #[test]
    fn test_gda_rejects_invalid_magic() {
        let data = b"not a valid GDA file at all!";
        let result = GenomicDistAnnotation::load_bin_from_bytes(data);
        let msg = match result {
            Err(e) => format!("{}", e),
            Ok(_) => panic!("Expected error loading invalid data"),
        };
        assert!(msg.contains("Invalid GDA"), "Error should mention GDA: {}", msg);
    }
}
