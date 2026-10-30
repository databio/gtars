# Models and Region Set objects in Gtars

Gtars has multiple objects (structs/models) for representation of genomic regions and other related data. 

### 🟢 Region

Region is Python representation of a genomic region. e.g. `chr1:100-200` + additional information.

#### Example

=== "Python"
    ```python
    from gtars.models import Region

    # Create a Region
    genomic_region = Region(chr="chr1", 
                             start=100, 
                             end=200, 
                             rest="peak1")
    print(genomic_region)

    ```

=== "Rust"
    ```rust
    use gtars_core::models::Region;
    
    // Create a Region
    let genomic_region: Region = Region { chr: "chr1".to_string(), 
                                          start: 100, 
                                          end: 200, 
                                          rest: Some("peak1".to_string()) 
                                        };
    let identifier = genomic_region.digest();

    println!("{:?}", identifier);
    
    ```


### 🟢 RegionSet

RegionSet is Python representation of a genomic region set, commonly named as BED file.


#### Quick example
Open BED file from URL and get its identifier.

=== "Python"
    ```python
    
    from gtars.models import RegionSet
    
    # Create a RegionSet from a url, or lcoal BED file.
    rs = RegionSet("https://data2.bedbase.org/files/d/a/dafd661aa70590999e0ff9e1980217db.bed.gz")
    
    # Get identifier for the RegionSet
    rs.identifier
    
    print(rs)

    ```
=== "Rust"
    ```rust
    use gtars_core::models::RegionSet;
    
    // Create a RegionSet from a url (requires the `http` feature), or lcoal BED file.
    let rs = RegionSet::try_from("https://data2.bedbase.org/files/d/a/dafd661aa70590999e0ff9e1980217db.bed.gz").unwrap();
    
    // Get identifier for the RegionSet
    let id = rs.identifier();

    println!("{:?}", rs);
    ```

=== "TypeScript"
    
    ❗ Note: This is test example and may require additional setup to run.
    
    ```typescript
    import init from '@databio/gtars';
    import { RegionSet } from '@databio/gtars';

    init();

    export type BedEntry1 = [string, number, number, string];
    
    // Define entries (regions)
    export const entries1: BedEntry1[] = [
      ['chr1', 100, 200, 'peak1'],
      ['chr2', 150, 250, 'peak2'],
      ['chr3', 300, 400, 'peak3'],
    ];

    // Create a Region
    const rs = new RegionSet(entries1);

    console.log(rs);

    ```

❗ Note: RegionSet can be created from a local file path, URL, or from a list (vector) of Region objects (in Python, `RegionSet.from_regions([...])`).

#### Main commands in Python

- Load a BED file from local path or URL
```python
rs = RegionSet("path/to/bedfile.bed")
```
- Get number of regions
```python
len(rs)
```
- Calculate mean reagion width
```python
rs.mean_region_width()
```
- Get last base pair location for each chromosome
```python
rs.get_max_end_per_chr()
```
- Get number of base pairs in the region set
```python
rs.get_nucleotide_length()
```
- Save the regionSet as a BED file
```python
rs.to_bed("path/to/save/bedfile.bed")
rs.to_bed_gz("path/to/save/bedfile.bed.gz")  # gzipped
```
- Save the regionSet as a bigBed file
```python
rs.to_bigbed("path/to/save/bedfile.bb", "path/to/chrom.sizes")
```

!!! info 
    - Detailed documentation for RegionSet is available in the [API reference](https://docs.rs/gtars-core/latest/gtars_core/models/region_set/).

### 🟢 RegionSetList

`RegionSetList` is the gtars equivalent of Bioconductor's `GRangesList` — an ordered collection of `RegionSet`s with optional names. It's the type downstream crates (genomicdist, lola) use to pass multiple region sets across FFI boundaries without paying N×clone costs.

#### Example

=== "Python"
    ```python
    from gtars.models import RegionSet, RegionSetList

    rs1 = RegionSet("rep1.bed")
    rs2 = RegionSet("rep2.bed")
    rs3 = RegionSet("rep3.bed")

    rsl = RegionSetList([rs1, rs2, rs3])

    print(len(rsl))                # number of sets
    first = rsl[0]                 # index like a list
    combined = rsl.concat()        # flatten into a single RegionSet (no merge)
    jac = rsl.pairwise_jaccard()   # N x N Jaccard matrix (list of lists)
    ```

=== "Rust"
    ```rust
    use gtars_core::models::{RegionSet, RegionSetList};

    let rs1 = RegionSet::try_from("rep1.bed")?;
    let rs2 = RegionSet::try_from("rep2.bed")?;
    let rs3 = RegionSet::try_from("rep3.bed")?;

    let rsl = RegionSetList::with_names(
        vec![rs1, rs2, rs3],
        vec!["rep1".into(), "rep2".into(), "rep3".into()],
    );

    // Iterate
    for rs in &rsl {
        println!("{} regions", rs.len());
    }

    // Flatten all regions into a single RegionSet (no merge/dedup)
    let combined = rsl.concat();
    let id = rsl.identifier();
    # Ok::<(), gtars_core::errors::RegionSetError>(())
    ```

`RegionSetList::try_from` in Rust also accepts a **bedset manifest file** (text file listing one BED path per line) or a `Vec<&Path>` / `Vec<&str>` / `Vec<String>` / `Vec<PathBuf>`.

`concat()` flattens without merging; if you need a reduced union, call `.reduce()` on the result. In Rust, `reduce` and the other interval set operations are methods on `RegionSet` in [gtars-core](core.md).

## See also

- **[gtars-core](core.md)** — the canonical Rust API reference for `Region`, `RegionSet`, `RegionSetList`, `Interval`, `Fragment`, `CoordinateMode`, `RegionSetError`, and the interval set operations.
- **[gtars-genomicdist](genomicdist.md)** — the `GenomicIntervalSetStatistics` trait that extends `RegionSet` with summary stats.
- **[gtars-lola](lola.md)** — LOLA enrichment, which consumes `RegionSetList` for user-set and database-set inputs.
