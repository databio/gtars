use std::ops::Add;

pub struct CountMatrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

pub struct RowIterator<'a, T> {
    matrix: &'a CountMatrix<T>,
    current_row: usize,
}


impl<T> CountMatrix<T>
where
    T: Copy + Default + Add<Output = T>,
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
                *value = *value + T::default();
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