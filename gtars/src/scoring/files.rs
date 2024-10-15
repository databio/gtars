use std::path::PathBuf;
use std::vec::IntoIter;

use anyhow::Result;
use glob::glob;

pub struct FragmentFileGlob {
    curr: usize,
    files: Vec<PathBuf>,
}

impl FragmentFileGlob {
    pub fn new(pattern: &str) -> Result<Self> {
        let files = glob(pattern)?;
        let files = files
            .map(|f| match f {
                Ok(path) => Ok(path),
                Err(_) => anyhow::bail!(format!("Error reading file entry: {:?}", f)),
            })
            .collect::<Result<Vec<_>>>()?;
        let curr = 0_usize;
        Ok(FragmentFileGlob { files, curr })
    }
}

impl Iterator for FragmentFileGlob {
    type Item = PathBuf;
    fn next(&mut self) -> Option<Self::Item> {
        let result = self.files.get(self.curr).cloned();
        self.curr +=1;
        result
    }
}