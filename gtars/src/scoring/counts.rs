use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::{Add, AddAssign};

use anyhow::Result;
use flate2::write::GzEncoder;
use flate2::Compression;

pub struct CountMatrix<T> {
    data: Vec<T>,
    pub rows: usize,
    pub cols: usize,
}

pub struct RowIterator<'a, T> {
    matrix: &'a CountMatrix<T>,
    current_row: usize,
}

impl<T> CountMatrix<T>
where
    T: Copy + Default + Add<Output = T> + AddAssign + From<u8>,
{
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![T::default(); rows * cols],
            rows,
            cols,
        }
    }

    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        self.data.get(row * self.cols + col)
    }

    pub fn set(&mut self, row: usize, col: usize, value: T) -> Result<(), String> {
        if row < self.rows && col < self.cols {
            self.data[row * self.cols + col] = value;
            Ok(())
        } else {
            Err(format!("Index out of bounds: row {}, col {}", row, col))
        }
    }

    pub fn increment(&mut self, row: usize, col: usize) {
        if row < self.rows && col < self.cols {
            let index = row * self.cols + col;
            if let Some(value) = self.data.get_mut(index) {
                *value += 1.into();
            }
        }
    }
}

impl<'a, T> Iterator for RowIterator<'a, T>
where
    T: Copy + Default,
{
    type Item = &'a [T];

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_row < self.matrix.rows {
            let start = self.current_row * self.matrix.cols;
            let end = start + self.matrix.cols;
            self.current_row += 1;
            Some(&self.matrix.data[start..end])
        } else {
            None
        }
    }
}

impl<T> CountMatrix<T>
where
    T: Copy + Default,
{
    pub fn iter_rows(&self) -> RowIterator<T> {
        RowIterator {
            matrix: self,
            current_row: 0,
        }
    }
}

impl<T> CountMatrix<T>
where
    T: Copy + Default + ToString,
{
    pub fn write_to_file(&self, filename: &str) -> std::io::Result<()> {
        let file = File::create(filename)?;
        let mut buf_writer = BufWriter::new(GzEncoder::new(file, Compression::default()));

        for row in self.iter_rows() {
            let row_str: String = row
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<String>>()
                .join(",");
            buf_writer.write_all(row_str.as_bytes())?;
            buf_writer.write_all(b"\n")?; // Add a newline after each row
        }

        Ok(())
    }
}
