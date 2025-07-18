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
            // Find the first dot to get the base name without any extensions
            let base_name = if let Some(dot_pos) = file_name_str.find('.') {
                &file_name_str[..dot_pos]
            } else {
                &file_name_str
            };

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
