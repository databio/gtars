//! FASTA label templates for store export.
//!
//! A template is literal text plus fields. A field names an object; an
//! optional `:NS` swaps that object's digest for its alias in namespace `NS`:
//!
//! | Field      | Gives                                                  |
//! |------------|--------------------------------------------------------|
//! | `{seq}`    | sequence digest as a refget ID (`SQ.…`), the default   |
//! | `{seq:NS}` | the sequence's alias in sequence-alias namespace `NS`  |
//! | `{col}`    | collection digest                                      |
//! | `{col:NS}` | the collection's alias in collection-alias namespace `NS` |
//! | `{name}`   | the sequence's name in its collection                  |
//!
//! Literal text is copied as-is (`{col:ucsc}/{name}`, `{col:pansn}#{name}`);
//! `{{` and `}}` write literal braces. A template that is exactly one field may
//! drop its braces: `name` means `{name}`, `seq:refseq` means `{seq:refseq}`.

use std::collections::HashMap;
use std::fmt;

use anyhow::Result;

use crate::hashkeyable::DigestKey;

/// A problem with a label template or with the labels it produces: a parse
/// error, an unknown namespace, a whitespace label or a label clash. Bindings
/// map it to a "bad value" error (Python `ValueError`).
#[derive(Debug)]
pub struct LabelError(pub String);

impl fmt::Display for LabelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for LabelError {}

/// Build an `anyhow::Error` carrying a [`LabelError`].
macro_rules! label_err {
    ($($arg:tt)*) => { anyhow::Error::new(LabelError(format!($($arg)*))) };
}
pub(crate) use label_err;

const VALID_FIELDS: &str =
    "valid fields are {seq}, {seq:NS}, {col}, {col:NS} and {name}, where NS is an alias namespace";

/// The object a label field refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LabelObject {
    Seq,
    Col,
    Name,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Part {
    Literal(String),
    Field { object: LabelObject, namespace: Option<String> },
}

/// A parsed FASTA label template. See the module docs for the syntax.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LabelTemplate {
    parts: Vec<Part>,
    source: String,
}

impl Default for LabelTemplate {
    fn default() -> Self {
        Self {
            parts: vec![Part::Field { object: LabelObject::Seq, namespace: None }],
            source: "{seq}".to_string(),
        }
    }
}

impl fmt::Display for LabelTemplate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.source)
    }
}

fn parse_field(body: &str, template: &str) -> Result<Part> {
    let (object_str, namespace) = match body.split_once(':') {
        Some((o, ns)) => (o, Some(ns)),
        None => (body, None),
    };
    let object = match object_str {
        "seq" => LabelObject::Seq,
        "col" => LabelObject::Col,
        "name" => LabelObject::Name,
        other => {
            return Err(label_err!(
                "Invalid label template '{}': unknown field '{{{}}}'; {}",
                template,
                other,
                VALID_FIELDS
            ))
        }
    };
    let namespace = match namespace {
        None => None,
        Some(_) if object == LabelObject::Name => {
            return Err(label_err!(
                "Invalid label template '{}': {{name}} takes no namespace (names have no aliases); {}",
                template,
                VALID_FIELDS
            ))
        }
        Some("") => {
            return Err(label_err!(
                "Invalid label template '{}': empty namespace in '{{{}}}'; {}",
                template,
                body,
                VALID_FIELDS
            ))
        }
        Some(ns) if ns.contains(':') || ns.chars().any(char::is_whitespace) => {
            return Err(label_err!(
                "Invalid label template '{}': bad namespace '{}'; {}",
                template,
                ns,
                VALID_FIELDS
            ))
        }
        Some(ns) => Some(ns.to_string()),
    };
    Ok(Part::Field { object, namespace })
}

/// Is `s` a single field written without braces (`name`, `seq:refseq`)?
fn is_shorthand(s: &str) -> bool {
    let object = s.split_once(':').map_or(s, |(o, _)| o);
    matches!(object, "seq" | "col" | "name") && !s.contains(['{', '}'])
}

impl LabelTemplate {
    /// Parse a template string.
    pub fn parse(template: &str) -> Result<Self> {
        if template.is_empty() {
            return Err(label_err!("Invalid label template: it is empty; {}", VALID_FIELDS));
        }
        if is_shorthand(template) {
            let part = parse_field(template, template)?;
            return Ok(Self { parts: vec![part], source: template.to_string() });
        }

        let mut parts: Vec<Part> = Vec::new();
        let mut literal = String::new();
        let mut chars = template.chars().peekable();
        while let Some(c) = chars.next() {
            match c {
                '{' if chars.peek() == Some(&'{') => {
                    chars.next();
                    literal.push('{');
                }
                '}' if chars.peek() == Some(&'}') => {
                    chars.next();
                    literal.push('}');
                }
                '{' => {
                    let mut body = String::new();
                    let mut closed = false;
                    for c in chars.by_ref() {
                        match c {
                            '}' => {
                                closed = true;
                                break;
                            }
                            '{' => break,
                            c => body.push(c),
                        }
                    }
                    if !closed {
                        return Err(label_err!(
                            "Invalid label template '{}': unbalanced braces; {}",
                            template,
                            VALID_FIELDS
                        ));
                    }
                    if !literal.is_empty() {
                        parts.push(Part::Literal(std::mem::take(&mut literal)));
                    }
                    parts.push(parse_field(&body, template)?);
                }
                '}' => {
                    return Err(label_err!(
                        "Invalid label template '{}': unbalanced braces (write '}}}}' for a literal brace); {}",
                        template,
                        VALID_FIELDS
                    ))
                }
                c => literal.push(c),
            }
        }
        if !literal.is_empty() {
            parts.push(Part::Literal(literal));
        }
        if !parts.iter().any(|p| matches!(p, Part::Field { .. })) {
            return Err(label_err!(
                "Invalid label template '{}': it has no fields, so every label would be the same; {}",
                template,
                VALID_FIELDS
            ));
        }
        Ok(Self { parts, source: template.to_string() })
    }

    /// Namespaces used with an object, deduplicated, in order of first use.
    pub(crate) fn namespaces(&self, which: LabelObject) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for p in &self.parts {
            if let Part::Field { object, namespace: Some(ns) } = p {
                if *object == which && !out.contains(ns) {
                    out.push(ns.clone());
                }
            }
        }
        out
    }

    /// Does the template use `{col}`, `{col:NS}` or `{name}`?
    pub(crate) fn needs_collection(&self) -> bool {
        self.parts.iter().any(|p| {
            matches!(p, Part::Field { object: LabelObject::Col | LabelObject::Name, .. })
        })
    }

    /// Fill in the template for one record.
    ///
    /// `collection` is `(collection digest, sequence name)` when the record
    /// was reached through a collection. Returns the label and whether any
    /// fallback happened (a missing alias, or a collection field on a record
    /// with no collection, which makes the whole label `{seq}`).
    pub(crate) fn render(
        &self,
        seq_digest: &str,
        collection: Option<(&str, &str)>,
        maps: &AliasMaps,
    ) -> (String, bool) {
        let seq_id = || format!("SQ.{}", seq_digest);
        if collection.is_none() && self.needs_collection() {
            return (seq_id(), true);
        }
        let mut fallback = false;
        let mut out = String::new();
        for p in &self.parts {
            match p {
                Part::Literal(s) => out.push_str(s),
                Part::Field { object: LabelObject::Seq, namespace: None } => out.push_str(&seq_id()),
                Part::Field { object: LabelObject::Seq, namespace: Some(ns) } => {
                    match maps.sequence(ns, seq_digest) {
                        Some(alias) => out.push_str(alias),
                        None => {
                            fallback = true;
                            out.push_str(&seq_id());
                        }
                    }
                }
                Part::Field { object: LabelObject::Col, namespace } => {
                    let col = collection.map(|c| c.0).unwrap_or_default();
                    match namespace {
                        None => out.push_str(col),
                        Some(ns) => match maps.collection(ns, col) {
                            Some(alias) => out.push_str(alias),
                            None => {
                                fallback = true;
                                out.push_str(col);
                            }
                        },
                    }
                }
                Part::Field { object: LabelObject::Name, .. } => {
                    out.push_str(collection.map(|c| c.1).unwrap_or_default())
                }
            }
        }
        (out, fallback)
    }
}

/// Digest -> chosen alias maps for the namespaces a template uses, built once
/// before planning so labeling a whole store stays linear.
#[derive(Debug, Default)]
pub(crate) struct AliasMaps {
    sequence: HashMap<String, HashMap<DigestKey, String>>,
    collection: HashMap<String, HashMap<DigestKey, String>>,
}

/// Invert `alias -> digest`, keeping the first alias in sorted order per digest.
fn invert(map: HashMap<String, DigestKey>) -> HashMap<DigestKey, String> {
    let mut pairs: Vec<(String, DigestKey)> = map.into_iter().collect();
    pairs.sort();
    let mut out: HashMap<DigestKey, String> = HashMap::new();
    for (alias, key) in pairs {
        out.entry(key).or_insert(alias);
    }
    out
}

impl AliasMaps {
    pub(crate) fn add_sequence_namespace(&mut self, ns: &str, map: HashMap<String, DigestKey>) {
        self.sequence.insert(ns.to_string(), invert(map));
    }

    pub(crate) fn add_collection_namespace(&mut self, ns: &str, map: HashMap<String, DigestKey>) {
        self.collection.insert(ns.to_string(), invert(map));
    }

    fn sequence(&self, ns: &str, digest: &str) -> Option<&str> {
        use crate::hashkeyable::HashKeyable;
        self.sequence.get(ns)?.get(&digest.to_key()).map(String::as_str)
    }

    fn collection(&self, ns: &str, digest: &str) -> Option<&str> {
        use crate::hashkeyable::HashKeyable;
        self.collection.get(ns)?.get(&digest.to_key()).map(String::as_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hashkeyable::HashKeyable;

    fn render(t: &str, coll: Option<(&str, &str)>, maps: &AliasMaps) -> (String, bool) {
        LabelTemplate::parse(t).unwrap().render("SEQDIGEST", coll, maps)
    }

    #[test]
    fn valid_templates_parse() {
        for t in [
            "{seq}", "{seq:refseq}", "{col}", "{col:ucsc}", "{name}",
            "seq", "seq:refseq", "col", "col:ucsc", "name",
            "{col:ucsc}/{name}", "{col:pansn}#{name}", "x{{{seq}}}y", "{name}_{seq}",
        ] {
            LabelTemplate::parse(t).unwrap_or_else(|e| panic!("{t}: {e}"));
        }
    }

    #[test]
    fn invalid_templates_fail_clearly() {
        for (t, needle) in [
            ("", "empty"),
            ("{foo}", "unknown field"),
            ("{name:ucsc}", "no namespace"),
            ("{seq:}", "empty namespace"),
            ("seq:", "empty namespace"),
            ("{seq", "unbalanced"),
            ("seq}", "unbalanced"),
            ("{seq}}", "unbalanced"),
            ("hello", "no fields"),
            ("{{seq}}", "no fields"),
            ("names", "no fields"),
        ] {
            let err = LabelTemplate::parse(t).unwrap_err().to_string();
            assert!(err.contains(needle), "{t:?}: {err}");
            assert!(err.contains("valid fields"), "{t:?}: {err}");
        }
    }

    #[test]
    fn shorthand_equals_braced() {
        let maps = AliasMaps::default();
        let c = Some(("COL", "chr1"));
        assert_eq!(render("name", c, &maps), render("{name}", c, &maps));
        assert_eq!(render("seq", c, &maps), render("{seq}", c, &maps));
        assert_eq!(render("col", c, &maps), render("{col}", c, &maps));
    }

    #[test]
    fn rendering_and_fallbacks() {
        let mut maps = AliasMaps::default();
        let mut seq_map = HashMap::new();
        seq_map.insert("NC_1".to_string(), "SEQDIGEST".to_key());
        maps.add_sequence_namespace("refseq", seq_map);
        let mut col_map = HashMap::new();
        col_map.insert("zz".to_string(), "COL".to_key());
        col_map.insert("hg38".to_string(), "COL".to_key());
        maps.add_collection_namespace("ucsc", col_map);

        let c = Some(("COL", "chr1"));
        assert_eq!(render("{seq}", c, &maps), ("SQ.SEQDIGEST".into(), false));
        assert_eq!(render("seq:refseq", c, &maps), ("NC_1".into(), false));
        // First sorted alias wins.
        assert_eq!(render("{col:ucsc}/{name}", c, &maps), ("hg38/chr1".into(), false));
        assert_eq!(render("{col}/{name}", c, &maps), ("COL/chr1".into(), false));
        assert_eq!(render("x{{{name}}}", c, &maps), ("x{chr1}".into(), false));
        // Missing alias falls back to the digest.
        assert_eq!(render("{col:ucsc}", Some(("OTHER", "a")), &maps), ("OTHER".into(), true));
        assert_eq!(render("{seq:ensembl}", c, &maps), ("SQ.SEQDIGEST".into(), true));
        // No collection: the whole label becomes {seq}.
        assert_eq!(render("{col}/{name}", None, &maps), ("SQ.SEQDIGEST".into(), true));
        assert_eq!(render("{seq:refseq}", None, &maps), ("NC_1".into(), false));
    }

    #[test]
    fn namespaces_and_needs_collection() {
        let t = LabelTemplate::parse("{col:a}/{seq:b}/{col:a}/{seq:c}").unwrap();
        assert_eq!(t.namespaces(LabelObject::Col), vec!["a"]);
        assert_eq!(t.namespaces(LabelObject::Seq), vec!["b", "c"]);
        assert!(t.needs_collection());
        assert!(!LabelTemplate::parse("seq:x").unwrap().needs_collection());
        assert_eq!(LabelTemplate::default().to_string(), "{seq}");
    }
}
