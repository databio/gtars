# Benchmark results

Measured on a 16-core, 31 GB machine shared with other jobs (the load is recorded in each row).
Every `*.jsonl` file holds one JSON object per run, written by `scripts/bench/bench.py`
(`label`, `rep`, `wall`, `maxrss_kb`, `load`, `rc`); `scripts/bench/summarize.py FILE` prints
medians per label. Docker runs use `ncbi/blast:2.17.0`.

| File | What |
|---|---|
| `bench_hg38.jsonl` | hg38 (jungle collection `EiFob05aCWgVU_B_Ae0cypnQut3cxUP1`, 195 sequences, 3.1 Gbp): makeblastdb against the database builder at 1 and 16 threads (`scripts/bench/bench_hg38.sh`) |
| `bench_hg38_rgmakeblastdb.jsonl` | hg38 with `rgmakeblastdb` and the workspace release profile, 1 and 16 threads |
| `bench_jungle.jsonl` | the whole jungle store (79 assemblies, 6,359 unique sequences, 41.8 Gbp; `scripts/bench/bench_jungle.sh`) |
| `bench_txome.jsonl` | GENCODE v47 transcripts (384,354 sequences, 560 Mbp) as a local store (`scripts/bench/txome_store.sh`) |
| `remote_fetch_hg38.jsonl` | hg38 built from the S3 store through gtars-refget's fetcher, cold and warm cache |
| `bench_rgblastn_search.jsonl` | rgblastn against blastn on hg38 and chr1 |
| `bench_jungle_search.jsonl` | rgblastn against blastn on the whole jungle store |
| `q_hg38.fa`, `q_jungle.fa` | the search queries (10 and 20) |

## Building databases

| Input | makeblastdb | rgmakeblastdb, 16 threads |
|---|---|---|
| hg38 (3.1 Gbp) | 62-87 s, plus 13-25 s to export FASTA | 0.55-2.2 s |
| jungle (41.8 Gbp) | 770-931 s, plus 104-195 s to export FASTA | 15-26 s |
| GENCODE transcripts (384,354 sequences) | 12.3-12.8 s (0.43 GB) | 6.9-7.8 s (1.55 GB) |

Every database is byte-identical to makeblastdb's on the same sequences (`tests/golden.rs`,
`scripts/golden/remote_verify.sh`). hg38 from S3 (1.36 GB packed, 32 parallel requests) took
177 s on a 7.7 MB/s link; the same command again reads only the cache and takes 1.25 s.

## Searching with rgblastn

Output is byte-identical to blastn 2.17.0 on every row below.

| Search | blastn | rgblastn |
|---|---|---|
| hg38 megablast, 10 queries, 1 thread | 42.5-45.5 s | 40.3-47.2 s (1.63 GB) |
| same, 4 threads | 24.6-27.2 s | 15.1-17.0 s (2.4 GB) |
| same, 8 threads | 21.8-22.4 s | 10.4-11.1 s (3.2 GB) |
| whole jungle megablast, 20 queries, 16 threads (809,785 rows) | 60.2-60.8 s | 45.3-45.6 s (7.65 GB) |
| GENCODE transcripts megablast, 1 thread | 2.2-2.3 s | 7.4-7.6 s (1.39 GB) |
| GENCODE transcripts `-task blastn`, 4 threads | 9.6-10.0 s | 14.3-14.8 s |

rgblastn is faster than blastn on large subjects and slower on stores of many small sequences,
where it pays about 18 µs per subject to plan and open the store's per-sequence files.
`-task blastn` on hg38 with the 10 queries together (11.1 million rows) follows blastn's
adaptive query batching: the queries run in batches of 6 + 1 + 1 + 1 + 1, as in blastn.
