# `gtars.lola`

Python bindings for [`gtars-lola`](../lola.md) — LOLA (Locus Overlap Analysis) enrichment testing. Given a user set of regions, a universe (background), and a region database, computes Fisher's exact test enrichment of the user set against every database region set.

See the [Rust gtars-lola page](../lola.md) for the full statistical reference (Fisher's exact test, CMLE odds ratio, contingency table semantics, R LOLA compatibility notes). This page is a Python-focused reference.

!!! warning "LOLA uses tuples, not `RegionSet`"
    The Python `gtars.lola` API takes **`list[tuple[str, int, int]]`** — lists of `(chr, start, end)` tuples — rather than `RegionSet` objects. This is a short-term quirk of the binding. If you have a `RegionSet` in hand, convert with:

    ```python
    def to_tuples(rs):
        return [(r.chr, r.start, r.end) for r in rs]
    ```

## End-to-end workflow

```python
from gtars.models import RegionSet
from gtars.lola import RegionDB, run_lola

# 1. Load the region database from a LOLA-format folder
db = RegionDB.from_folder(
    "LOLACore/hg38",
    collections=None,   # None = all; or pass a list of collection names
    limit=None,         # None = all files; or pass a per-collection cap
)

# 2. Prepare user sets and universe as lists of (chr, start, end) tuples
peaks = RegionSet("peaks.bed")
universe = RegionSet("universe.bed")

def to_tuples(rs):
    return [(r.chr, r.start, r.end) for r in rs]

user_sets = [to_tuples(peaks)]          # list of user sets; each is a list of tuples
universe_tuples = to_tuples(universe)

# 3. Run the enrichment
results = run_lola(
    user_sets,
    universe_tuples,
    db,
    min_overlap=1,
    direction="enrichment",   # or "depletion"
)

# results is a column-oriented dict — DataFrame-friendly
import pandas as pd
df = pd.DataFrame(results)
df = df.sort_values("pValueLog", ascending=False)
print(df.head(20))
```

The 4-step workflow mirrors R LOLA: load database → convert inputs → run enrichment → pandas-ify.

## `RegionDB`

Wraps a `gtars-lola` `RegionDB` — an IGD overlap index plus original region sets plus per-file annotations, loaded from the standard LOLA directory layout:

```text
db_path/
├── collection1/
│   ├── collection.txt      # collector, date, source, description (TSV/CSV)
│   ├── index.txt           # per-file annotations (TSV/CSV)
│   └── regions/
│       ├── file1.bed
│       └── file2.bed
├── collection2/
│   └── …
```

### `RegionDB.from_folder`

```python
RegionDB.from_folder(db_path, collections=None, limit=None) -> RegionDB
```

- `db_path` — path to the LOLA database root.
- `collections` — optional list of collection names to include (others are skipped).
- `limit` — optional per-collection file cap (useful for smoke tests).

`collection.txt` and `index.txt` are auto-detected TSV or CSV. Recognized `index.txt` columns: `filename`, `cellType`, `description`, `tissue`, `dataSource`, `antibody`, `treatment`. Unknown columns are ignored; missing optional columns become `None`.

### `RegionDB.from_bed_files`

Build a database from a list of BED file paths without the full LOLA layout:

```python
db = RegionDB.from_bed_files(
    ["file1.bed", "file2.bed", "file3.bed"],
    filenames=None,   # optional; defaults to basenames
)
```

Useful for ad-hoc analyses where you have a list of peak files but don't want to set up a collection directory.

### Accessors

```python
db.num_region_sets         # int — total files in the database

db.list_region_sets()                       # list[str] — all filenames
db.list_region_sets(collections=["enc3"])   # filtered

# Extract region sets as a gtars.models.RegionSetList (names populated from filenames)
rsl = db.get_region_sets()                  # all sets
rsl = db.get_region_sets(indices=[0, 5, 12])

# Annotations as lists of dicts
db.region_anno       # per-file dicts with filename/cellType/description/tissue/...
db.collection_anno   # per-collection dicts with collectionname/collector/date/source/description
```

`get_region_sets()` returns a [`gtars.models.RegionSetList`](models.md#regionsetlist), which has `names` populated from the database filenames. This is the one path in Python where `RegionSetList.names` is non-None — normal `RegionSetList(...)` construction can't set names.

## Running enrichment

### `run_lola`

```python
run_lola(
    user_sets: list[list[tuple[str, int, int]]],
    universe: list[tuple[str, int, int]],
    region_db: RegionDB,
    min_overlap: int = 1,
    direction: str = "enrichment",
) -> dict
```

**Parameters:**

- `user_sets` — a **list of user sets**, where each user set is a list of `(chr, start, end)` tuples. Pass one-element list `[peaks_tuples]` for a single-user-set analysis.
- `universe` — a single list of `(chr, start, end)` tuples representing the background.
- `region_db` — a `RegionDB`.
- `min_overlap` — minimum bp overlap to count as overlapping (default 1).
- `direction` — `"enrichment"` (default, P(X ≥ a), alternative "greater") or `"depletion"` (P(X ≤ a), alternative "less"). The strings `"greater"` / `"less"` are accepted as aliases.

**Returns** a column-oriented dict with these parallel lists (one entry per `(user_set, db_set)` pair):

| key | type | meaning |
|---|---|---|
| `userSet` | `list[int]` | 0-based user-set index |
| `dbSet` | `list[int]` | 0-based db-set index |
| `collection` | `list[str | None]` | from `index.txt` |
| `pValueLog` | `list[float]` | `-log10(p)` from Fisher's exact test, capped at ~322 |
| `oddsRatio` | `list[float]` | CMLE odds ratio (matches R `fisher.test()$estimate`) |
| `support` | `list[int]` | overlap count between user set and db set — the contingency `a` |
| `rnkPV`, `rnkOR`, `rnkSup` | `list[int]` | 1-based per-metric ranks within the user set |
| `maxRnk` | `list[int]` | max of the three ranks |
| `meanRnk` | `list[float]` | mean of the three ranks |
| `b`, `c`, `d` | `list[int]` | signed contingency values (can be negative if user set extends outside universe) |
| `qValue` | `list[float | None]` | BH-adjusted p-value (applied automatically inside `run_lola`) |
| `description`, `cellType`, `tissue`, `antibody`, `treatment`, `dataSource` | `list[str | None]` | from `index.txt` |
| `filename` | `list[str]` | db set file name |
| `size` | `list[int]` | number of regions in the db set |

The Python `run_lola` already calls `annotate_results` + `apply_fdr_correction` internally, so `qValue` and the annotation columns come back populated — you don't need to post-process.

Rows are sorted by descending `pValueLog`, then ascending `meanRnk` (matching R LOLA output order).

!!! note "Negative contingency values"
    If your user set contains regions *outside* the universe, the contingency table produces negative `b`/`c`/`d`. The binding prints a warning to stderr, stores the signed values, and returns `pValueLog = 0.0`, `oddsRatio = NaN` for that row. To avoid this, preprocess with `redefine_user_sets` (below).

## Universe preparation

LOLA results are only meaningful when your universe contains your user sets. Three helpers for checking and shaping the universe:

### `check_universe`

```python
from gtars.lola import check_universe

report = check_universe(
    user_sets,           # list[list[tuple[str, int, int]]]
    universe_tuples,     # list[tuple[str, int, int]]
)

# report is a dict with parallel lists:
report["userSet"]           # list[int]
report["totalRegions"]      # list[int]
report["regionsInUniverse"] # list[int]
report["coverage"]          # list[float]  — 0.0 to 1.0
report["manyToMany"]        # list[int]   — number of user regions overlapping >1 universe region
report["warnings"]          # list[str]   — free-form warning messages
```

Warnings fire when coverage is under 50% (severe) or under 90% (moderate), and when any user region overlaps more than one universe region.

### `redefine_user_sets`

Rewrite each user set in terms of universe regions — for every user region, find the universe regions it overlaps and emit *those* as the new user set. Eliminates many-to-many mapping artifacts. This is the Python equivalent of R LOLA's `redefineUserSets()`.

```python
from gtars.lola import redefine_user_sets

redefined = redefine_user_sets(user_sets, universe_tuples)
# redefined is list[list[tuple[str, int, int]]] — same shape as input, but rewritten

# Now run enrichment on the redefined version
results = run_lola(redefined, universe_tuples, db)
```

### `build_restricted_universe`

For differential enrichment analysis — build a universe that is exactly the union of all user sets, disjoined into non-overlapping pieces. This is R LOLA's `disjoin(unlist(userSets))`.

```python
from gtars.lola import build_restricted_universe

restricted = build_restricted_universe(user_sets)
# restricted: list[tuple[str, int, int]]

# Use it as the universe in a differential run
results = run_lola(user_sets, restricted, db)
```

## Complete example with pandas

```python
import pandas as pd
from gtars.models import RegionSet
from gtars.lola import RegionDB, run_lola, check_universe, redefine_user_sets

# Convert a RegionSet to the tuple form lola expects
def to_tuples(rs):
    return [(r.chr, r.start, r.end) for r in rs]

# Load inputs
db = RegionDB.from_folder("LOLACore/hg38")
peaks_a = to_tuples(RegionSet("condition_a.bed"))
peaks_b = to_tuples(RegionSet("condition_b.bed"))
universe = to_tuples(RegionSet("universe.bed"))

user_sets = [peaks_a, peaks_b]

# Universe sanity check
report = check_universe(user_sets, universe)
for us, cov, mm in zip(
    report["userSet"], report["coverage"], report["manyToMany"]
):
    print(f"user set {us}: {cov:.1%} coverage, {mm} many-to-many")
for w in report["warnings"]:
    print(f"  ⚠ {w}")

# Optional: eliminate many-to-many artifacts
user_sets = redefine_user_sets(user_sets, universe)

# Run enrichment
results = run_lola(user_sets, universe, db, direction="enrichment")

df = pd.DataFrame(results)
print(df.shape)   # (n_user_sets * n_db_sets, 23)

# Top 10 per user set
for us in df["userSet"].unique():
    top = (
        df[df["userSet"] == us]
        .sort_values("pValueLog", ascending=False)
        .head(10)
    )
    print(f"\n=== User set {us} — top 10 ===")
    print(top[["filename", "cellType", "pValueLog", "oddsRatio", "support", "qValue"]])
```

## See also

- **[`gtars.models`](models.md)** — `RegionSet`, `RegionSetList`, and the types returned by `RegionDB.get_region_sets()`.
- **[gtars-lola](../lola.md)** — full Rust API reference, including contingency table math, p-value computation via hypergeometric survival function, CMLE odds ratio, and the R LOLA TSV output format.
