use std::io::{self, Write};

use anyhow::{Context, Result};
use serde::Serialize;

/// Output format for CLI results.
pub enum Format {
    Table,
    Json,
    JsonCompact,
    Jsonl,
    Tsv,
    #[cfg(feature = "sc-parquet")]
    Parquet,
}

impl Format {
    pub fn from_str(s: &str) -> Self {
        match s {
            "json" => Format::Json,
            "jsonl" => Format::Jsonl,
            "tsv" => Format::Tsv,
            #[cfg(feature = "sc-parquet")]
            "parquet" => Format::Parquet,
            _ => Format::Table,
        }
    }

    /// Apply --compact flag: converts Json to JsonCompact.
    pub fn with_compact(self, compact: bool) -> Self {
        if compact {
            match self {
                Format::Json => Format::JsonCompact,
                other => other,
            }
        } else {
            self
        }
    }

    pub fn is_json(&self) -> bool {
        matches!(self, Format::Json | Format::JsonCompact | Format::Jsonl)
    }
}

/// Write a serializable value as JSON to stdout (pretty-printed).
pub fn write_json<T: Serialize>(value: &T) -> Result<()> {
    let json = serde_json::to_string_pretty(value).context("serializing to JSON")?;
    io::stdout().write_all(json.as_bytes())?;
    println!();
    Ok(())
}

/// Write a serializable value as compact JSON (single line) to stdout.
pub fn write_json_compact<T: Serialize>(value: &T) -> Result<()> {
    let json = serde_json::to_string(value).context("serializing to JSON")?;
    io::stdout().write_all(json.as_bytes())?;
    println!();
    Ok(())
}

/// Write a sequence of serializable values as JSONL (one JSON object per line).
pub fn write_jsonl<T: Serialize>(values: &[T]) -> Result<()> {
    let out = io::stdout();
    let mut w = out.lock();
    for value in values {
        let line = serde_json::to_string(value).context("serializing to JSONL")?;
        writeln!(w, "{}", line)?;
    }
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

/// Write tabular data as a Parquet file.
#[cfg(feature = "sc-parquet")]
pub fn write_parquet(
    path: &std::path::Path,
    headers: &[&str],
    columns: &[ParquetColumn],
) -> Result<()> {
    use arrow::array::{Float64Array, StringArray, UInt32Array};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;

    let fields: Vec<Field> = headers
        .iter()
        .zip(columns.iter())
        .map(|(name, col)| match col {
            ParquetColumn::Str(_) => Field::new(*name, DataType::Utf8, false),
            ParquetColumn::F64(_) => Field::new(*name, DataType::Float64, false),
            ParquetColumn::U32(_) => Field::new(*name, DataType::UInt32, false),
        })
        .collect();

    let schema = Arc::new(Schema::new(fields));

    let arrays: Vec<Arc<dyn arrow::array::Array>> = columns
        .iter()
        .map(|col| -> Arc<dyn arrow::array::Array> {
            match col {
                ParquetColumn::Str(v) => Arc::new(StringArray::from(v.clone())),
                ParquetColumn::F64(v) => Arc::new(Float64Array::from(v.clone())),
                ParquetColumn::U32(v) => Arc::new(UInt32Array::from(v.clone())),
            }
        })
        .collect();

    let batch = RecordBatch::try_new(schema.clone(), arrays)
        .context("creating Arrow record batch")?;

    let file = std::fs::File::create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    let mut writer = ArrowWriter::try_new(file, schema, None)
        .context("creating Parquet writer")?;

    writer.write(&batch).context("writing Parquet batch")?;
    writer.close().context("closing Parquet writer")?;

    Ok(())
}

/// Column data for Parquet output.
#[cfg(feature = "sc-parquet")]
pub enum ParquetColumn {
    Str(Vec<String>),
    F64(Vec<f64>),
    U32(Vec<u32>),
}
