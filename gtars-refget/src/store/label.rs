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
//! Literal text is copied as-is (`{col:ucsc}/{name}`, `{col:pansn}#{name}`).
//! Every field is written with its braces, and a template cannot contain a
//! literal brace.

use std::collections::HashMap;
use std::fmt;

use anyhow::Result;

use super::alias::AliasKind;
use crate::hashkeyable::{DigestKey, HashKeyable};

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

/// One piece of a template: literal text, or a field with its alias
/// namespace, if any.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Part {
    Lit(String),
    Seq(Option<String>),
    Col(Option<String>),
    Name,
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
            parts: vec![Part::Seq(None)],
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
    let (object, namespace) = match body.split_once(':') {
        Some((o, ns)) => (o, Some(ns)),
        None => (body, None),
    };
    if !matches!(object, "seq" | "col" | "name") {
        return Err(label_err!(
            "Invalid label template '{}': unknown field '{{{}}}'; {}",
            template,
            object,
            VALID_FIELDS
        ));
    }
    let namespace = match namespace {
        None => None,
        Some(_) if object == "name" => {
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
    Ok(match object {
        "seq" => Part::Seq(namespace),
        "col" => Part::Col(namespace),
        _ => Part::Name,
    })
}

/// For each alias namespace a template uses: digest -> the alias it labels
/// with. Built once before planning so labeling a whole store stays linear.
pub(crate) type AliasIndex = HashMap<String, HashMap<DigestKey, String>>;

/// Invert `alias -> digest`, keeping the first alias in sorted order per digest.
pub(crate) fn invert(map: HashMap<String, DigestKey>) -> HashMap<DigestKey, String> {
    let mut pairs: Vec<(String, DigestKey)> = map.into_iter().collect();
    pairs.sort();
    let mut out: HashMap<DigestKey, String> = HashMap::new();
    for (alias, key) in pairs {
        out.entry(key).or_insert(alias);
    }
    out
}

fn alias<'a>(index: &'a AliasIndex, ns: &str, digest: &str) -> Option<&'a str> {
    index.get(ns)?.get(&digest.to_key()).map(String::as_str)
}

impl LabelTemplate {
    /// Parse a template string.
    pub fn parse(template: &str) -> Result<Self> {
        if template.is_empty() {
            return Err(label_err!("Invalid label template: it is empty; {}", VALID_FIELDS));
        }
        let unbalanced = || {
            label_err!("Invalid label template '{}': unbalanced braces; {}", template, VALID_FIELDS)
        };
        let mut parts: Vec<Part> = Vec::new();
        let mut rest = template;
        while let Some(open) = rest.find(['{', '}']) {
            if rest[open..].starts_with('}') {
                return Err(unbalanced());
            }
            if open > 0 {
                parts.push(Part::Lit(rest[..open].to_string()));
            }
            let field = &rest[open + 1..];
            let close = match field.find(['{', '}']) {
                Some(i) if field[i..].starts_with('}') => i,
                _ => return Err(unbalanced()),
            };
            parts.push(parse_field(&field[..close], template)?);
            rest = &field[close + 1..];
        }
        if !rest.is_empty() {
            parts.push(Part::Lit(rest.to_string()));
        }
        if parts.iter().all(|p| matches!(p, Part::Lit(_))) {
            return Err(label_err!(
                "Invalid label template '{}': it has no fields, so every label would be the same; {}",
                template,
                VALID_FIELDS
            ));
        }
        Ok(Self { parts, source: template.to_string() })
    }

    /// Alias namespaces of `kind` the template uses, deduplicated, in order
    /// of first use.
    pub(crate) fn namespaces(&self, kind: AliasKind) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for p in &self.parts {
            let ns = match (kind, p) {
                (AliasKind::Sequence, Part::Seq(Some(ns))) => ns,
                (AliasKind::Collection, Part::Col(Some(ns))) => ns,
                _ => continue,
            };
            if !out.contains(ns) {
                out.push(ns.clone());
            }
        }
        out
    }

    /// Does the template use `{col}`, `{col:NS}` or `{name}`?
    pub(crate) fn needs_collection(&self) -> bool {
        self.parts.iter().any(|p| matches!(p, Part::Col(_) | Part::Name))
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
        seq_aliases: &AliasIndex,
        col_aliases: &AliasIndex,
    ) -> (String, bool) {
        let seq_id = || format!("SQ.{}", seq_digest);
        if collection.is_none() && self.needs_collection() {
            return (seq_id(), true);
        }
        let mut fallback = false;
        let mut out = String::new();
        for p in &self.parts {
            match p {
                Part::Lit(s) => out.push_str(s),
                Part::Seq(None) => out.push_str(&seq_id()),
                Part::Seq(Some(ns)) => match alias(seq_aliases, ns, seq_digest) {
                    Some(a) => out.push_str(a),
                    None => {
                        fallback = true;
                        out.push_str(&seq_id());
                    }
                },
                Part::Col(namespace) => {
                    let col = collection.map(|c| c.0).unwrap_or_default();
                    match namespace {
                        None => out.push_str(col),
                        Some(ns) => match alias(col_aliases, ns, col) {
                            Some(a) => out.push_str(a),
                            None => {
                                fallback = true;
                                out.push_str(col);
                            }
                        },
                    }
                }
                Part::Name => out.push_str(collection.map(|c| c.1).unwrap_or_default()),
            }
        }
        (out, fallback)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render(t: &str, coll: Option<(&str, &str)>, maps: &(AliasIndex, AliasIndex)) -> (String, bool) {
        LabelTemplate::parse(t).unwrap().render("SEQDIGEST", coll, &maps.0, &maps.1)
    }

    #[test]
    fn valid_templates_parse() {
        for t in [
            "{seq}", "{seq:refseq}", "{col}", "{col:ucsc}", "{name}",
            "{col:ucsc}/{name}", "{col:pansn}#{name}", "x{seq}y", "{name}_{seq}",
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
            ("{seq", "unbalanced"),
            ("seq}", "unbalanced"),
            ("{seq}}", "unbalanced"),
            ("{{seq}}", "unbalanced"),
            ("hello", "no fields"),
            ("name", "no fields"),
            ("seq:refseq", "no fields"),
            ("names", "no fields"),
        ] {
            let err = LabelTemplate::parse(t).unwrap_err().to_string();
            assert!(err.contains(needle), "{t:?}: {err}");
            assert!(err.contains("valid fields"), "{t:?}: {err}");
        }
    }

    #[test]
    fn rendering_and_fallbacks() {
        let mut seq_map = HashMap::new();
        seq_map.insert("NC_1".to_string(), "SEQDIGEST".to_key());
        let mut col_map = HashMap::new();
        col_map.insert("zz".to_string(), "COL".to_key());
        col_map.insert("hg38".to_string(), "COL".to_key());
        let maps = (
            AliasIndex::from([("refseq".to_string(), invert(seq_map))]),
            AliasIndex::from([("ucsc".to_string(), invert(col_map))]),
        );

        let c = Some(("COL", "chr1"));
        assert_eq!(render("{seq}", c, &maps), ("SQ.SEQDIGEST".into(), false));
        assert_eq!(render("{seq:refseq}", c, &maps), ("NC_1".into(), false));
        // First sorted alias wins.
        assert_eq!(render("{col:ucsc}/{name}", c, &maps), ("hg38/chr1".into(), false));
        assert_eq!(render("{col}/{name}", c, &maps), ("COL/chr1".into(), false));
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
        assert_eq!(t.namespaces(AliasKind::Collection), vec!["a"]);
        assert_eq!(t.namespaces(AliasKind::Sequence), vec!["b", "c"]);
        assert!(t.needs_collection());
        assert!(!LabelTemplate::parse("{seq:x}").unwrap().needs_collection());
        assert_eq!(LabelTemplate::default().to_string(), "{seq}");
    }
}
