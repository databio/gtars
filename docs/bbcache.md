# BED file caching and loading from BEDbase

## Introduction

With `gtars.bbcache`, users can download and cache BED files and BED sets from the BEDbase API to their local disk, as well as process, search, and remove cached files.

This document provides tutorials for using `gtars.bbcache` via either:

1. the [Rust interface](#rust-interface), or
2. the [command-line interface](#command-line-interface).

## Rust interface

### Create an instance of the BBClient Class:

```rust
use gtars::bbcache::client::BBClient;
use std::path::PathBuf;

let cache_folder = PathBuf::from("cache");
let mut bbclient = BBClient::new(Some(cache_folder), None).expect("Failed to create BBClient");
```

### Download and cache a remote BED file from BEDbase / Load a BED file from cache

```rust
let bedfile_id = "..."; // find interesting bedfile on bedbase
// download, cache and return a RegionSet object
// or return a RegionSet from a cached BED file
let regionset = bbclient.load_bed(bedfile_id).expect("Failed to load bed file");
```

### Cache a local BED file

```rust
// compute its ID and add it to the cache
let bedfile_id = bbclient.add_local_bed_to_cache(PathBuf::from("path/to/bedfile"), Some(false)).unwrap();
```

### Cache a BED file from within Python memory

You can also provide a URL and it will add to cache for you:

```rust
use gtars::common::models::RegionSet;

let file_url = String::from("https://...");
let regionset = RegionSet::try_from(file_url)?;
let bedfile_id = bbclient.add_regionset_to_cache(regionset, Some(false)).unwrap();
```

### Download and cache a BEDset from BEDbase

```rust
let bedset_identifier = "xyz"; // find some interesting bedset on bedbase.org
let bedset = bbclient.load_bedset(bedset_identifier).unwrap(); //download, cache and return a BedSet object
```

### Cache a local BEDset

```rust
let bedset_id = bbc.add_local_folder_as_bedset(PathBuf::from("path/go/bed/files/folder")).unwrap();
```

### Seek for a cached BED file / BEDset

```rust
let bedset_id = "...";
let bedfile_id = "...";

let cached_bed_path = bbclient.seek(bedfile_id).expect("Failed to seek cached bed file");
let cached_bedset_path = bbclient.seek(bedset_id).expect("Failed to seek cached bed file");
```

### Remove a cached BED file / BEDset

```rust
let bedset_id = "...";
let bedfile_id = "...";

bbc.remove(bedset_id).expect("Failed to remove bedset file and its bed files");
bbc.remove(bedfile_id).expect("Failed to remove bedset file and its bed files");
```

## Command line interface

### Overall

```
Downloads, processes, and caches BED files from the BEDbase API

Usage: bbcache <COMMAND>

Commands:
  cache-bed           Cache a BED file from local file or BEDbase
  cache-bedset        Cache a BED set from local file or BEDbase
  seek                Seek the BED file path by giving identifier
  inspect-bedfiles    Inspect the contents of bedfile cache folder
  inspect-bedsets     Inspect the contents of bedsets cache folder
  rm                  Remove the BED file or BED set from cache with given identifier
  help                Print this message or the help of the given subcommand(s)

Options (for applicable subcommands):
  -i, --identifier <identifier>    BED file or BED set identifier, URL, or file path
  -f, --cache-folder <cache-folder>  Cache folder path
  -h, --help                       Print help
```

### Cache BED file

```bash
gtars bbcache cache-bed -i <BED_file_or_identifier_or_url>
```

The `<BED_file_or_identifier_or_url>` variable can be one of 3 things:

1. a path to a local BED file;
2. a BED record identifier from BEDbase; or,
3. a URL to a BED file hosted anywhere.

### Cache BEDset

```bash
gtars bbcache cache-bedset -i <BED_files_folder_or_identifier>
```

The `<BED_files_folder_or_identifier>` variable may be:

1. local path to a folder containing BED files; or,
2. a BEDbase BEDset identifier

### Seek the path of a BED file or BEDset in cache folder

To retrieve the local file path to a BED file stored locally,

```bash
gtars bbcache seek -i <identifier>
```

Replace <identifier> with the identifier of the BED file or BEDset you want to seek.


### Count the subdirectories and files in `bedfiles` & `bedsets` folder

```bash
gtars bbcache inspect-bedfiles
gtars bbcache inspect-bedsets
```

### Remove a BED file or BEDset from the cache folder 

```bash
gtars bbcache rm  <identifier>
```

Replace <identifier> with the identifier of the BED file or BEDset you want to remove.


### Cache Folder

By default, the downloaded and processed BED files are cached in the bed_cache folder. You can specify a different cache folder using the `--cache-folder` argument, or set the environment variable `$BBCLIENT_CACHE_ENV`.
The cache folder has this structure:
```
cache_folder
  bedfiles
    a/b/ab1234xyz.bed.gz
    ..
  bedsets
    c/d/cd123hij.txt
```
