use std::io::{self, Write};

use anyhow::{Context, Result};
use serde::Serialize;

/// Output format for CLI results.
pub enum Format {
    Table,
    Json,
    Tsv,
}

impl Format {
    pub fn from_str(s: &str) -> Self {
        match s {
            "json" => Format::Json,
            "tsv" => Format::Tsv,
            _ => Format::Table,
        }
    }
}

/// Write a serializable value as JSON to stdout.
pub fn write_json<T: Serialize>(value: &T) -> Result<()> {
    let json = serde_json::to_string_pretty(value).context("serializing to JSON")?;
    io::stdout().write_all(json.as_bytes())?;
    println!();
    Ok(())
}

/// Write a TSV header + rows to stdout.
pub fn write_tsv(headers: &[&str], rows: &[Vec<String>]) -> Result<()> {
    let out = io::stdout();
    let mut w = out.lock();
    writeln!(w, "{}", headers.join("\t"))?;
    for row in rows {
        writeln!(w, "{}", row.join("\t"))?;
    }
    Ok(())
}

/// Write a human-readable table to stdout.
pub fn write_table(headers: &[&str], rows: &[Vec<String>]) -> Result<()> {
    if rows.is_empty() {
        return Ok(());
    }

    // Compute column widths
    let n_cols = headers.len();
    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if i < n_cols {
                widths[i] = widths[i].max(cell.len());
            }
        }
    }

    let out = io::stdout();
    let mut w = out.lock();

    // Header
    for (i, header) in headers.iter().enumerate() {
        if i > 0 {
            write!(w, "  ")?;
        }
        write!(w, "{:>width$}", header, width = widths[i])?;
    }
    writeln!(w)?;

    // Rows
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if i >= n_cols {
                break;
            }
            if i > 0 {
                write!(w, "  ")?;
            }
            write!(w, "{:>width$}", cell, width = widths[i])?;
        }
        writeln!(w)?;
    }

    Ok(())
}
