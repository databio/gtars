pub enum UniverseFileType {
    BedThree,
    BedFivePlus,
    Unknown,
}

impl From<&String> for UniverseFileType {
    fn from(value: &String) -> Self {
        if value.starts_with("track") {
            return UniverseFileType::Unknown;
        }
        let parts: Vec<&str> = value.split('\t').collect();
        if parts.len() == 3 {
            return UniverseFileType::BedThree;
        } else if parts.len() >= 5 {
            return UniverseFileType::BedFivePlus;
        }
        UniverseFileType::Unknown
    }
}
