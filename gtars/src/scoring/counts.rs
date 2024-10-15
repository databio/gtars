use std::ops::Add;

pub struct CountMatrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
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
