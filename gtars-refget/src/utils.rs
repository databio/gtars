use std::path::{Path, PathBuf};
// ...existing imports...

/// Extension trait to add replace_ext_with method to Path types
pub trait PathExtension {
    fn replace_exts_with(&self, new_extension: &str) -> PathBuf;
}

impl<T: AsRef<Path>> PathExtension for T {
    fn replace_exts_with(&self, new_extension: &str) -> PathBuf {
        let path = self.as_ref();
        if let Some(file_name) = path.file_name() {
            let file_name_str = file_name.to_string_lossy();

            // Strip only known FASTA extensions, preserving dots in the base name
            // Check longest extensions first (.fa.gz, .fasta.gz) before shorter ones
            let known_extensions = [".fa.gz", ".fasta.gz", ".fa", ".fasta"];

            let base_name = known_extensions.iter()
                .find_map(|ext| file_name_str.strip_suffix(ext))
                .unwrap_or(&file_name_str);

            if let Some(parent) = path.parent() {
                parent.join(format!("{}.{}", base_name, new_extension))
            } else {
                PathBuf::from(format!("{}.{}", base_name, new_extension))
            }
        } else {
            // Fallback if no file name exists
            path.with_extension(new_extension)
        }
    }
}
