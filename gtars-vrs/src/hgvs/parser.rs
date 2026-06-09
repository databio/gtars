//! Hand-written HGVS parser.
//!
//! Supports the v1 HGVS grammar slice needed for c./n./g./m./r. variants
//! with substitution, del, dup, ins, delins, inv, identity (`=`), and unknown
//! (`?`) edits. Handles single positions and ranges, intronic offsets,
//! 5' UTR (`c.-N`), 3' UTR (`c.*N`), uncertain positions in parentheses,
//! and optional gene-symbol annotation `ACC(GENE):...`.

use super::ast::{
    Datum, Edit, HgvsVariant, LocationRange, PosEdit, Position, ReferenceType,
};
use super::error::HgvsError;

/// Parse an HGVS variant string into an AST.
pub fn parse(input: &str) -> Result<HgvsVariant<'_>, HgvsError> {
    let mut p = Parser::new(input);
    let v = p.parse_variant()?;
    if !p.is_eof() {
        return Err(p.error("trailing characters after variant"));
    }
    Ok(v)
}

struct Parser<'a> {
    input: &'a str,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Self { input, pos: 0 }
    }

    fn error(&self, msg: &str) -> HgvsError {
        HgvsError::Parse {
            input: self.input.to_string(),
            pos: self.pos,
            msg: msg.to_string(),
        }
    }

    fn is_eof(&self) -> bool {
        self.pos >= self.input.len()
    }

    fn peek(&self) -> Option<u8> {
        self.input.as_bytes().get(self.pos).copied()
    }

    fn consume_byte(&mut self) -> Option<u8> {
        let b = self.peek()?;
        self.pos += 1;
        Some(b)
    }

    fn expect(&mut self, b: u8, ctx: &str) -> Result<(), HgvsError> {
        match self.peek() {
            Some(p) if p == b => {
                self.pos += 1;
                Ok(())
            }
            _ => Err(self.error(ctx)),
        }
    }

    fn try_consume(&mut self, b: u8) -> bool {
        if self.peek() == Some(b) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn try_consume_keyword(&mut self, kw: &str) -> bool {
        let bytes = kw.as_bytes();
        if self.input.as_bytes()[self.pos..].starts_with(bytes) {
            self.pos += bytes.len();
            true
        } else {
            false
        }
    }

    fn parse_variant(&mut self) -> Result<HgvsVariant<'a>, HgvsError> {
        // Accession: chars up to `:` or `(`
        let acc_start = self.pos;
        while let Some(b) = self.peek() {
            if b == b':' || b == b'(' {
                break;
            }
            self.pos += 1;
        }
        if self.pos == acc_start {
            return Err(self.error("expected accession"));
        }
        let accession = &self.input[acc_start..self.pos];

        // Optional gene in parentheses.
        let mut gene = None;
        if self.try_consume(b'(') {
            let g_start = self.pos;
            while let Some(b) = self.peek() {
                if b == b')' {
                    break;
                }
                self.pos += 1;
            }
            if self.pos == g_start {
                return Err(self.error("expected gene symbol after `(`"));
            }
            gene = Some(&self.input[g_start..self.pos]);
            self.expect(b')', "expected `)` after gene symbol")?;
        }

        self.expect(b':', "expected `:` after accession")?;

        // Reference type: g., c., n., m., r., p.
        let rt = match self.consume_byte() {
            Some(b'g') => ReferenceType::G,
            Some(b'c') => ReferenceType::C,
            Some(b'n') => ReferenceType::N,
            Some(b'm') => ReferenceType::M,
            Some(b'r') => ReferenceType::R,
            Some(b'p') => ReferenceType::P,
            _ => return Err(self.error("expected reference type (g/c/n/m/r/p)")),
        };
        self.expect(b'.', "expected `.` after reference type")?;

        // Check for whole-sequence edits (e.g., `g.=` or `g.?`)
        let (pos, edit, uncertain) = if self.peek() == Some(b'=') || self.peek() == Some(b'?') {
            let edit = self.parse_edit(rt)?;
            (LocationRange::WholeSequence, edit, false)
        } else {
            // Check for outer uncertain wrapper: `(...)`
            // This is tricky because `(` can also start inner uncertain positions.
            // Peek ahead to distinguish:
            // - `(1799T>A)` - outer uncertain (number followed by edit)
            // - `(4_6)_246del` - inner uncertain (number_number)
            let outer_uncertain = if self.peek() == Some(b'(') {
                // Save position, peek inside
                let saved = self.pos;
                self.pos += 1; // consume `(`
                // Skip `?` for unknown position
                if self.peek() == Some(b'?') {
                    self.pos = saved;
                    false // `(?_...)` is inner uncertain
                } else {
                    // Skip optional sign and number
                    if self.peek() == Some(b'-') || self.peek() == Some(b'*') {
                        self.pos += 1;
                    }
                    while self.peek().map(|b| b.is_ascii_digit()).unwrap_or(false) {
                        self.pos += 1;
                    }
                    // Also skip intronic offset like +5 or -3
                    if self.peek() == Some(b'+') || self.peek() == Some(b'-') {
                        self.pos += 1;
                        while self.peek().map(|b| b.is_ascii_digit()).unwrap_or(false) {
                            self.pos += 1;
                        }
                    }
                    // If next char is `_`, it's inner uncertain; otherwise outer
                    let is_inner = self.peek() == Some(b'_');
                    self.pos = saved; // restore
                    !is_inner
                }
            } else {
                false
            };

            if outer_uncertain {
                self.pos += 1; // consume `(`
                let pos = self.parse_location_range(rt)?;
                let edit = self.parse_edit(rt)?;
                self.expect(b')', "expected `)` to close uncertain posedit")?;
                (pos, edit, true)
            } else {
                // Position(s) - may include uncertain ranges like (4_6)_246
                let pos = self.parse_location_range(rt)?;
                let edit = self.parse_edit(rt)?;
                // Mark as uncertain if the location itself is uncertain
                let uncertain = matches!(
                    pos,
                    LocationRange::UncertainStart { .. }
                        | LocationRange::UncertainEnd { .. }
                        | LocationRange::UncertainBoth { .. }
                );
                (pos, edit, uncertain)
            }
        };

        Ok(HgvsVariant {
            accession,
            gene,
            reference_type: rt,
            posedit: PosEdit { pos, edit, uncertain },
        })
    }

    fn parse_location_range(
        &mut self,
        rt: ReferenceType,
    ) -> Result<LocationRange, HgvsError> {
        // Check for uncertain start: (pos1_pos2) or (?_pos)
        let (start_pos, start_is_uncertain, start_low, start_high) = if self.try_consume(b'(') {
            let (low, high) = self.parse_uncertain_position_pair(rt)?;
            self.expect(b')', "expected `)` after uncertain position")?;
            // Use whichever is available as the "main" position
            let main = high.or(low).ok_or_else(|| self.error("both bounds unknown"))?;
            (main, true, low, high)
        } else {
            let pos = self.parse_position(rt)?;
            (pos, false, None, None)
        };

        if self.try_consume(b'_') {
            // Check for uncertain end: _(pos1_pos2) or _(pos_?)
            let (end_pos, end_is_uncertain, end_low, end_high) = if self.try_consume(b'(') {
                let (low, high) = self.parse_uncertain_position_pair(rt)?;
                self.expect(b')', "expected `)` after uncertain position")?;
                let main = low.or(high).ok_or_else(|| self.error("both bounds unknown"))?;
                (main, true, low, high)
            } else {
                let pos = self.parse_position(rt)?;
                (pos, false, None, None)
            };

            match (start_is_uncertain, end_is_uncertain) {
                (true, true) => Ok(LocationRange::UncertainBoth {
                    start_low,
                    start_high,
                    end_low,
                    end_high,
                }),
                (true, false) => Ok(LocationRange::UncertainStart {
                    start_low,
                    start_high,
                    end: end_pos,
                }),
                (false, true) => Ok(LocationRange::UncertainEnd {
                    start: start_pos,
                    end_low,
                    end_high,
                }),
                (false, false) => Ok(LocationRange::Range {
                    start: start_pos,
                    end: end_pos,
                }),
            }
        } else {
            Ok(LocationRange::Single(start_pos))
        }
    }

    /// Parse `pos1_pos2` or `?_pos` or `pos_?` inside uncertain parens.
    /// Returns (low: Option<Position>, high: Option<Position>).
    fn parse_uncertain_position_pair(
        &mut self,
        rt: ReferenceType,
    ) -> Result<(Option<Position>, Option<Position>), HgvsError> {
        let low = if self.try_consume(b'?') {
            None
        } else {
            Some(self.parse_position(rt)?)
        };
        self.expect(b'_', "expected `_` in uncertain position range")?;
        let high = if self.try_consume(b'?') {
            None
        } else {
            Some(self.parse_position(rt)?)
        };
        Ok((low, high))
    }

    fn parse_position(&mut self, rt: ReferenceType) -> Result<Position, HgvsError> {
        // Protein positions have amino acid prefix: Ala1, Met100, Ter50, etc.
        if rt == ReferenceType::P {
            return self.parse_protein_position();
        }

        // Datum: detect `*` (CDS end) for `c.*N`; otherwise CDS start for `c.`,
        // else SeqStart.
        let mut datum = match rt {
            ReferenceType::C => Datum::CdsStart,
            _ => Datum::SeqStart,
        };

        if rt == ReferenceType::C && self.try_consume(b'*') {
            datum = Datum::CdsEnd;
        }

        // Sign for the base position (allowed for c. coordinates: c.-14).
        let mut neg_base = false;
        if self.peek() == Some(b'-') {
            // Could be a negative base or part of an offset? Inside parse_position
            // for a c. coord, "-" means negative base. We commit to negative base.
            neg_base = true;
            self.pos += 1;
        } else if self.peek() == Some(b'+') {
            // Explicit positive sign isn't standard for HGVS base positions; skip.
            self.pos += 1;
        }

        let base = self.parse_uint()? as i64;
        let base = if neg_base { -base } else { base };

        // Optional intronic offset: `+N` or `-N` after the base.
        let mut offset: i64 = 0;
        match self.peek() {
            Some(b'+') => {
                self.pos += 1;
                let n = self.parse_uint()? as i64;
                offset = n;
            }
            Some(b'-') => {
                self.pos += 1;
                let n = self.parse_uint()? as i64;
                offset = -n;
            }
            _ => {}
        }

        Ok(Position { base, offset, datum })
    }

    /// Parse protein position like `Ala1`, `Met100`, `Ter50`, `*10`.
    fn parse_protein_position(&mut self) -> Result<Position, HgvsError> {
        // Handle Ter/*
        let datum = if self.try_consume(b'*') || self.try_consume_keyword("Ter") {
            Datum::CdsEnd
        } else {
            // 3-letter or 1-letter amino acid code
            let aa_start = self.pos;
            // Try 3-letter first
            if self.peek().map(|b| b.is_ascii_uppercase()).unwrap_or(false) {
                self.pos += 1;
                // Check for lowercase letters (3-letter code like Ala, Met)
                while let Some(b) = self.peek() {
                    if b.is_ascii_lowercase() {
                        self.pos += 1;
                    } else {
                        break;
                    }
                }
            }
            if self.pos == aa_start {
                return Err(self.error("expected amino acid"));
            }
            Datum::SeqStart
        };

        let base = self.parse_uint()? as i64;
        Ok(Position {
            base,
            offset: 0,
            datum,
        })
    }

    /// Parse protein edit like `Ser` (substitution), `del`, `fs`, etc.
    fn parse_protein_edit(&mut self) -> Result<Edit<'a>, HgvsError> {
        // Check for common protein edit keywords
        if self.try_consume_keyword("del") {
            return Ok(Edit::Del { reference: None });
        }
        if self.try_consume_keyword("dup") {
            return Ok(Edit::Dup { reference: None });
        }
        if self.try_consume_keyword("ins") {
            let alt = self.parse_amino_acid_seq()?;
            return Ok(Edit::Ins { alternate: alt });
        }
        if self.try_consume_keyword("fs") {
            // Frameshift - consume optional Ter/*/number
            while let Some(b) = self.peek() {
                if b.is_ascii_alphanumeric() || b == b'*' {
                    self.pos += 1;
                } else {
                    break;
                }
            }
            return Ok(Edit::Unknown); // Represent frameshift as unknown for now
        }

        // Simple substitution: just an amino acid like `Ser`, `Ter`, `*`
        let alt = self.parse_amino_acid_seq()?;
        // For protein substitution, we don't have the ref in the edit string
        Ok(Edit::Sub {
            reference: "",
            alternate: alt,
        })
    }

    /// Parse one or more amino acids (3-letter or 1-letter codes).
    fn parse_amino_acid_seq(&mut self) -> Result<&'a str, HgvsError> {
        let start = self.pos;
        while let Some(b) = self.peek() {
            if b.is_ascii_alphabetic() || b == b'*' {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err(self.error("expected amino acid"));
        }
        Ok(&self.input[start..self.pos])
    }

    fn parse_uint(&mut self) -> Result<u64, HgvsError> {
        let start = self.pos;
        while let Some(b) = self.peek() {
            if b.is_ascii_digit() {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err(self.error("expected integer"));
        }
        self.input[start..self.pos]
            .parse::<u64>()
            .map_err(|_| self.error("invalid integer"))
    }

    fn parse_edit(&mut self, rt: ReferenceType) -> Result<Edit<'a>, HgvsError> {
        // Identity?
        if self.try_consume(b'=') {
            return Ok(Edit::Identity);
        }
        // Unknown?
        if self.try_consume(b'?') {
            return Ok(Edit::Unknown);
        }

        // Protein edits: amino acid substitution like `Ser` after position `Ala1`
        if rt == ReferenceType::P {
            return self.parse_protein_edit();
        }

        // del / delins / dup / ins / inv
        if self.try_consume_keyword("delins") {
            let alt = self.parse_iupac_run()?;
            return Ok(Edit::DelIns {
                reference: None,
                alternate: alt,
            });
        }
        if self.try_consume_keyword("del") {
            // Optional reference allele; may be followed by `ins<seq>`.
            let reference = self.parse_optional_iupac_run();
            if self.try_consume_keyword("ins") {
                let alt = self.parse_iupac_run()?;
                return Ok(Edit::DelIns {
                    reference,
                    alternate: alt,
                });
            }
            // Also allow `delN` where N is a count (legacy/lax HGVS) - just consume it.
            if reference.is_none() {
                while let Some(b) = self.peek() {
                    if b.is_ascii_digit() {
                        self.pos += 1;
                    } else {
                        break;
                    }
                }
            }
            return Ok(Edit::Del { reference });
        }
        if self.try_consume_keyword("dup") {
            let reference = self.parse_optional_iupac_run();
            return Ok(Edit::Dup { reference });
        }
        if self.try_consume_keyword("ins") {
            let alt = self.parse_iupac_run()?;
            return Ok(Edit::Ins { alternate: alt });
        }
        if self.try_consume_keyword("inv") {
            let reference = self.parse_optional_iupac_run();
            return Ok(Edit::Inv { reference });
        }
        // Invitae extension: `copyN`
        if self.try_consume_keyword("copy") {
            let count = self.parse_uint()? as u32;
            return Ok(Edit::Copy { count });
        }

        // Substitution: <REF>><ALT> or identity with explicit ref: <REF>=
        let ref_start = self.pos;
        while let Some(b) = self.peek() {
            if is_iupac(b) {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == ref_start {
            return Err(self.error("expected edit"));
        }
        let reference = &self.input[ref_start..self.pos];

        // Check for identity with explicit ref (e.g., `G=`)
        if self.try_consume(b'=') {
            return Ok(Edit::Identity);
        }

        // Check for repeat notation (e.g., `CA[4]`)
        if self.try_consume(b'[') {
            let count = self.parse_uint()? as u32;
            self.expect(b']', "expected `]` after repeat count")?;
            return Ok(Edit::Repeat {
                sequence: reference,
                count,
            });
        }

        self.expect(b'>', "expected `>` in substitution")?;
        let alt_start = self.pos;
        while let Some(b) = self.peek() {
            if is_iupac(b) {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == alt_start {
            return Err(self.error("expected alternate allele"));
        }
        let alternate = &self.input[alt_start..self.pos];
        Ok(Edit::Sub {
            reference,
            alternate,
        })
    }

    fn parse_iupac_run(&mut self) -> Result<&'a str, HgvsError> {
        let s = self.pos;
        while let Some(b) = self.peek() {
            if is_iupac(b) {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == s {
            return Err(self.error("expected nucleotide sequence"));
        }
        Ok(&self.input[s..self.pos])
    }

    fn parse_optional_iupac_run(&mut self) -> Option<&'a str> {
        let s = self.pos;
        while let Some(b) = self.peek() {
            if is_iupac(b) {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == s {
            None
        } else {
            Some(&self.input[s..self.pos])
        }
    }
}

fn is_iupac(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N'
            | b'a' | b'c' | b'g' | b't' | b'u' | b'n'
            | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
            | b'r' | b'y' | b's' | b'w' | b'k' | b'm' | b'b' | b'd' | b'h' | b'v'
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genomic_sub() {
        let v = parse("NC_000007.14:g.140753336A>T").unwrap();
        assert_eq!(v.accession, "NC_000007.14");
        assert!(matches!(v.reference_type, ReferenceType::G));
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, 140753336);
                assert_eq!(p.offset, 0);
                assert!(matches!(p.datum, Datum::SeqStart));
            }
            _ => panic!("expected single position"),
        }
        match v.posedit.edit {
            Edit::Sub { reference, alternate } => {
                assert_eq!(reference, "A");
                assert_eq!(alternate, "T");
            }
            _ => panic!("expected substitution"),
        }
    }

    #[test]
    fn test_parse_coding_sub() {
        let v = parse("NM_004333.6:c.1799T>A").unwrap();
        assert_eq!(v.accession, "NM_004333.6");
        assert!(matches!(v.reference_type, ReferenceType::C));
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, 1799);
                assert!(matches!(p.datum, Datum::CdsStart));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_five_prime_utr() {
        let v = parse("NM_004333.6:c.-14C>T").unwrap();
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, -14);
                assert!(matches!(p.datum, Datum::CdsStart));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_three_prime_utr() {
        let v = parse("NM_004333.6:c.*37A>G").unwrap();
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, 37);
                assert!(matches!(p.datum, Datum::CdsEnd));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_intronic_positive_offset() {
        let v = parse("NM_004333.6:c.1798+5G>A").unwrap();
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, 1798);
                assert_eq!(p.offset, 5);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_intronic_negative_offset() {
        let v = parse("NM_004333.6:c.1799-3T>A").unwrap();
        match v.posedit.pos {
            LocationRange::Single(p) => {
                assert_eq!(p.base, 1799);
                assert_eq!(p.offset, -3);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_del_with_reference() {
        let v = parse("NM_004333.6:c.100delAGT").unwrap();
        match v.posedit.edit {
            Edit::Del { reference } => {
                assert_eq!(reference, Some("AGT"));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_del_no_reference() {
        let v = parse("NM_004333.6:c.100del").unwrap();
        match v.posedit.edit {
            Edit::Del { reference } => assert_eq!(reference, None),
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_dup() {
        let v = parse("NM_004333.6:c.100_102dupATG").unwrap();
        match v.posedit.pos {
            LocationRange::Range { start, end } => {
                assert_eq!(start.base, 100);
                assert_eq!(end.base, 102);
            }
            _ => panic!(),
        }
        match v.posedit.edit {
            Edit::Dup { reference } => assert_eq!(reference, Some("ATG")),
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_ins() {
        let v = parse("NM_004333.6:c.100_101insATG").unwrap();
        match v.posedit.edit {
            Edit::Ins { alternate } => assert_eq!(alternate, "ATG"),
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_delins() {
        let v = parse("NM_004333.6:c.100_102delinsCT").unwrap();
        match v.posedit.edit {
            Edit::DelIns { reference, alternate } => {
                assert_eq!(reference, None);
                assert_eq!(alternate, "CT");
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_delins_with_ref() {
        let v = parse("NM_004333.6:c.100_102delAGTinsCT").unwrap();
        match v.posedit.edit {
            Edit::DelIns { reference, alternate } => {
                assert_eq!(reference, Some("AGT"));
                assert_eq!(alternate, "CT");
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_identity() {
        let v = parse("NM_004333.6:c.100=").unwrap();
        assert!(matches!(v.posedit.edit, Edit::Identity));
    }

    #[test]
    fn test_parse_gene_in_parens() {
        let v = parse("NM_004333.6(BRAF):c.1799T>A").unwrap();
        assert_eq!(v.accession, "NM_004333.6");
        assert_eq!(v.gene, Some("BRAF"));
    }

    #[test]
    fn test_parse_gene_symbol_accession() {
        let v = parse("BRAF:c.1799T>A").unwrap();
        assert_eq!(v.accession, "BRAF");
    }

    #[test]
    fn test_parse_uncertain() {
        let v = parse("NM_004333.6:c.(1799T>A)").unwrap();
        assert!(v.posedit.uncertain);
        assert!(matches!(v.posedit.edit, Edit::Sub { .. }));
    }

    #[test]
    fn test_parse_inv() {
        let v = parse("NC_000007.14:g.100_200inv").unwrap();
        assert!(matches!(v.posedit.edit, Edit::Inv { reference: None }));
    }

    #[test]
    fn test_parse_range() {
        let v = parse("NM_004333.6:c.1798_1800delAGT").unwrap();
        match v.posedit.pos {
            LocationRange::Range { start, end } => {
                assert_eq!(start.base, 1798);
                assert_eq!(end.base, 1800);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_parse_error_no_colon() {
        let err = parse("NM_004333.6").unwrap_err();
        assert!(matches!(err, HgvsError::Parse { .. }));
    }

    #[test]
    fn test_parse_error_bad_reftype() {
        let err = parse("NM_004333.6:x.100A>T").unwrap_err();
        assert!(matches!(err, HgvsError::Parse { .. }));
    }
}
