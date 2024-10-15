use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use rust_lapper::{Interval, Lapper};

use crate::common::models::{Region, RegionSet};
use crate::common::utils::extract_regions_from_bed_file;

