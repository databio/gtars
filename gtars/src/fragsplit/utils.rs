use std::path::Path;

pub fn remove_all_extensions(path: &Path) -> String {
    let mut stem = path.file_stem().unwrap().to_string_lossy().to_string();

    let mut parent_path = path.with_file_name(stem.clone());
    while let Some(_extension) = parent_path.extension() {
        // Remove the extension by recreating the path without it
        parent_path = parent_path.with_extension("");
        stem = parent_path
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
    }

    stem
}
