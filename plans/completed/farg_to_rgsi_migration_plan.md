# Plan: Complete `.farg` to `.rgsi` Migration

## Overview

The gtars-refget codebase underwent a partial migration from `.farg` file extension to `.rgsi`/`.rgci` extensions. The **RefgetStore** code was updated, but the **SequenceCollection** caching layer was not. This results in `.farg` files being auto-generated when loading FASTA files, polluting directories in downstream projects (like the refget Python package).

**Goal**: Complete the migration by eliminating all `.farg` references and replacing them with `.rgsi`.

**Important**: DO NOT MAINTAIN BACKWARDS COMPATIBILITY. This is developmental software. We are eliminating the old format entirely, not keeping it around.

## Current State

### What's Using `.rgsi` (new format):
- `RefgetStore.write_store_to_dir()` - writes `sequences.rgsi`, `collections.rgci`, `rgstore.json`
- `RefgetStore.write_index_files()` - writes `.rgsi` files
- `RefgetStore.write_collection_to_disk_single()` - writes `collections/{digest}.rgsi`

### What's Still Using `.farg` (old format):
- `SequenceCollection::from_fasta()` - creates `.farg` cache files next to FASTA files
- `SequenceCollection::from_path_with_cache()` - reads/writes `.farg` files
- `SequenceCollection::from_farg()` - reads `.farg` files
- `SequenceCollection::write_farg()` - writes `.farg` files
- `SequenceCollection::write_collection_farg()` - writes `.farg` files
- `RefgetStore.write_sequences_farg()` - called by `write_sequences_rgsi()`
- `RefgetStore.load_local()` - backward-compat code checks for `.farg` extension
- Test files in `tests/data/fasta/*.farg`
- `.gitignore` references `sequences.farg`
- Python binding docstrings mention `.farg`

## Files to Read Before Implementation

1. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/gtars-refget/src/collection.rs` - Lines 460-620 (SequenceCollection caching methods)
2. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/gtars-refget/src/store.rs` - Lines 1340-1450, 1500-1600, 1870-1910, 1990-2010, 2287-2320, 3049-3090 (farg references)
3. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/gtars-refget/src/fasta.rs` - Lines 610-620, 786-795 (read_fasta_refget_file, test)
4. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/gtars-python/src/refget/mod.rs` - Lines 1220-1260 (docstrings)
5. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/gtars-python/py_src/gtars/refget/__init__.pyi` - Line 249 (docstring)
6. `/home/nsheff/Dropbox/workspaces/refgenie/repos/gtars/.gitignore` - Line 30

## Implementation Steps

### Step 1: Update `collection.rs` - Rename Methods and Extensions

**File**: `gtars-refget/src/collection.rs`

1. Rename `from_farg()` → `from_rgsi()` and change extension from `"farg"` to `"rgsi"`
2. In `from_path_with_cache()`:
   - Rename `farg_file_path` → `rgsi_file_path`
   - Change `replace_exts_with("farg")` → `replace_exts_with("rgsi")`
   - Update all print statements: "farg" → "rgsi"
   - Change call from `write_farg()` → `write_rgsi()`
3. Rename `write_collection_farg()` → `write_collection_rgsi()` and update print statement
4. Rename `write_farg()` → `write_rgsi()`:
   - Change `replace_exts_with("farg")` → `replace_exts_with("rgsi")`
   - Change call from `write_collection_farg()` → `write_collection_rgsi()`
   - Update error message: "FARG" → "RGSI"
5. Remove the existing `write_rgsi()` method at line ~602 (it just calls `write_collection_farg`, so it's redundant after rename)

### Step 2: Update `store.rs` - Rename and Remove Backward Compatibility

**File**: `gtars-refget/src/store.rs`

1. Rename `write_sequences_farg()` → `write_sequences_rgsi()` (line ~1999):
   - Update print statement: "farg" → "rgsi"
   - Update internal references

2. Update `write_sequences_rgsi()` wrapper (line ~1431):
   - It currently just calls `write_sequences_farg()`, make it the primary implementation
   - Or simply inline the code and remove the wrapper

3. In `load_local()` (line ~1502):
   - Remove comment about "old format"
   - Remove `is_farg` checks
   - Remove the code that checks for `.farg` extension (only check for `.rgsi`)

4. In `load_remote()` (around line ~1830):
   - Remove `is_farg` checks
   - Remove the code that checks for `.farg` extension

5. In `get_collection()` (around line ~1883):
   - Remove `relative_path_old` (the `.farg` fallback)
   - Only try `relative_path_new` (`.rgsi`)

6. Update test `store_fa_to_farg` (line ~2287):
   - Rename to `store_fa_to_rgsi`
   - Change `temp_farg` → `temp_rgsi`
   - Change `"base.farg"` → `"base.rgsi"`
   - Update call from `write_farg()` → `write_rgsi()`

7. Update test `test_farg_filename_with_dots` (line ~3049):
   - Rename to `test_rgsi_filename_with_dots`
   - Change all `.farg` references to `.rgsi`
   - Update expected filenames and assertions
   - Update print statement

8. Remove/update comments that say "renamed from .farg" - they're no longer relevant

### Step 3: Update `fasta.rs` - Rename Function and Test

**File**: `gtars-refget/src/fasta.rs`

1. Rename function `read_fasta_refget_file()` → `read_rgsi_file()` (line ~615)
   - Consider keeping the generic name since it reads the TSV format, not specifically "fasta refget"

2. Update test `digests_fa_to_farg` (line ~786):
   - Rename to `digests_fa_to_rgsi`
   - Change `write_farg()` → `write_rgsi()`
   - Change file path from `"../tests/data/fasta/base.farg"` → `"../tests/data/fasta/base.rgsi"`

### Step 4: Update All Imports and References

**File**: `gtars-refget/src/collection.rs` (line 3)
- If renaming `read_fasta_refget_file` → `read_rgsi_file`, update import

**File**: `gtars-refget/src/store.rs` (line 37)
- If renaming `read_fasta_refget_file` → `read_rgsi_file`, update import

### Step 5: Update Python Bindings

**File**: `gtars-python/src/refget/mod.rs`

1. Update docstring at line ~1227: change `sequences.farg` → `sequences.rgsi`
2. Update docstring at line ~1250: change `sequences.farg` → `sequences.rgsi`

**File**: `gtars-python/py_src/gtars/refget/__init__.pyi`

1. Update docstring at line ~249: remove mention of "old format (index.json, sequences.farg)"

### Step 6: Update Test Data Files

**Directory**: `tests/data/fasta/`

1. Delete all `.farg` files:
   - `base.farg`
   - `crlf_endings.farg`
   - `HG002.alt.pat.f1_v2.unmasked.farg`
   - `unwrapped.farg`
   - `wrapped_20.farg`
   - `wrapped_40.farg`

2. The tests will regenerate them as `.rgsi` files when run

### Step 7: Update .gitignore

**File**: `.gitignore`

1. Change line 30 from `/gtars-refget/tests/store_test/sequences.farg` to `/gtars-refget/tests/store_test/sequences.rgsi`

### Step 8: Run Tests and Fix Any Remaining References

1. Run `cargo test` in `gtars-refget/`
2. Search for any remaining "farg" strings: `grep -r "farg" gtars-refget/src/`
3. Fix any compilation errors or test failures
4. Run Python tests if applicable

## Summary of Changes

After implementation, verify:

| Before | After |
|--------|-------|
| `SequenceCollection::from_farg()` | `SequenceCollection::from_rgsi()` |
| `SequenceCollection::write_farg()` | `SequenceCollection::write_rgsi()` |
| `SequenceCollection::write_collection_farg()` | `SequenceCollection::write_collection_rgsi()` |
| `RefgetStore::write_sequences_farg()` | `RefgetStore::write_sequences_rgsi()` |
| `read_fasta_refget_file()` | `read_rgsi_file()` (optional rename) |
| `.farg` extension | `.rgsi` extension |
| Backward-compat `.farg` checks | Removed |
| 6 `.farg` test data files | Deleted (regenerate as `.rgsi`) |

## Final Verification

1. `grep -r "farg" gtars-refget/` should return no results
2. `grep -r "farg" gtars-python/` should return no results
3. All 62+ Rust tests pass
4. All Python tests pass
5. No `.farg` files exist in `tests/data/fasta/`

## IMPORTANT: NO BACKWARDS COMPATIBILITY

This is developmental software. We are NOT:
- Keeping fallback code that reads old `.farg` files
- Adding migration utilities
- Maintaining dual-format support
- Preserving the old method names as aliases

We ARE:
- Completely eliminating `.farg` as a concept
- Renaming all methods and variables
- Removing all backward-compatibility code paths
- Deleting old test data files
