# BLAST on a RefgetStore: design

This document explains the BLAST support added to gtars (the `gtars-rgblast` crate) and to the
`refget` command line: what is new, how it relates to NCBI BLAST+ and to RefgetStore, and what
changed on the build side and the query side so that BLAST works on a store. It describes branch
`feat/gtars-rgblast` as of 2026-11-09. Commands and flags are in the user guide,
[docs/blast.md](../../docs/blast.md).

## Terms

| Term | Meaning |
|---|---|
| RefgetStore | A content-addressed store of reference sequences ([gtars-refget](../../gtars-refget)). Each sequence is stored once, under its refget digest. |
| refget digest | A sequence's identifier from the GA4GH [refget sequences](https://ga4gh.github.io/refget/sequences/) standard: `sha512t24u`, the first 24 bytes of the sequence's SHA-512 in base64url (32 characters). Written `SQ.<digest>` in ids. |
| collection | An ordered set of named sequences, usually one genome assembly, with its own digest under the GA4GH [sequence collections](https://ga4gh.github.io/refget/seqcols/) (seqcol) standard. Collections can share sequences. |
| alias | A readable name for a collection or sequence within a namespace, such as `name:hg38-p14-ucsc` or `refseq:GCF_000001405.39`. |
| FHR | [FAIR Headers Reference genome](https://github.com/FAIR-bioHeaders/FHR-Specification) metadata per collection: species, taxid, assembly accession. |
| OID | Ordinal id: BLAST's index of a sequence inside a database (0, 1, 2, ...). |
| defline | One description of a database sequence: a list of ids, a title and an optional taxid. In a BLAST version 5 database one sequence can carry several deflines. |
| SeqDB | The BLAST+ library that opens database files and hands sequences to the search engine. |
| alias file (`.nal`) | A small text file that presents one or more databases, or a subset of one, as a single database. It can also set the statistics BLAST uses (`STATS_NSEQ`, `STATS_TOTLEN`). |
| seqidlist (`.bsl`) | A binary list of sequence ids that limits a search or an alias file to those sequences. |
| ncbi2na, blastna | BLAST's 2-bit code (A=0, C=1, G=2, T=3) and its one-letter-per-byte code that includes ambiguity letters. |

## 1. What existed before

### BLAST+

In [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/), `makeblastdb` reads a FASTA file and
writes a version 5 database; `blastn` reads it through SeqDB.

```
FASTA --makeblastdb--> .nin .nhr .nsq  .ndb .ntf  .nos .not .nto  .njs  [.nal]  --SeqDB--> blastn
```

`.nin` is the index, `.nhr` the deflines, `.nsq` the sequences (ncbi2na plus a table of ambiguous
bases); `.ndb` and `.ntf` are LMDB lookups (id to OID, taxid to offset); `.nos`, `.not` and `.nto`
map OIDs to ids and taxids; `.njs` is JSON metadata. Alias files and seqidlists define subsets.

### RefgetStore

A RefgetStore keeps each unique sequence once, in a file named by its digest
(`sequences/<2 chars>/<digest>.seq`). Index files list sequences, collections and aliases; FHR
metadata sits in `fhr/<digest>.fhr.json`; `rgstore.json` holds digests of the store's content.
Sequences are stored as packed codes (2-bit for pure ACGT, in UCSC order T, C, A, G; 3-bit for
ACGTN and a few more letters; 4-bit for full IUPAC), as raw ASCII, or as zstd-compressed ASCII. A
store is a local directory or any plain HTTP prefix, such as the public stores on S3. The "jungle"
store used below holds 79 human and mouse assemblies: 6,359 unique sequences, 41.8 Gbp.

### The old path, and what it cost

You exported FASTA (`gtars refget export`) and ran `makeblastdb`. The database duplicated the store
(about 11 GB for jungle) and had to be rebuilt to stay current. Measured on 16 cores on 2026-11-08
against the prototype that became `rgmakeblastdb`
([benchmarks/results](../benchmarks/results/README.md)):

| Task | Export + makeblastdb | Native build (16 threads) |
|---|---|---|
| hg38 (195 sequences, 3.1 Gbp) | 13–25 s + 63–85 s | 1.1–2.2 s (median 1.3 s) |
| jungle, every assembly | 104–195 s + 770–932 s (about 17 min) | 16.5–26 s (median 22 s) |
| hg38 from the S3 store | 657 s (export piped into makeblastdb) | 49.6 s (32 parallel requests) |

## 2. What is new

```
                       RefgetStore (local or http)
                                 |
                select (collections, taxa, N filter)
                                 |
              plan::SubjectPlan: OIDs, deflines, taxids, totals
                 /                                   \
   rgmakeblastdb (database writer)             rgblastn engine
   same bytes as makeblastdb                   reads the store directly
                 |                                   |
   stock blastn / blastdbcmd (SeqDB)           same tabular bytes as blastn
                 \                                   /
                  rgblast: search, report, origin, reps
```

The `gtars-rgblast` crate ships three binaries:

| Binary | Verbs | Role |
|---|---|---|
| `rgmakeblastdb` | `build`, `verify`, `view`, `nstats` | Writes BLAST version 5 databases from a store, checks them against it, writes subset views, counts N bases |
| `rgblast` | `search`, `report`, `origin`, `reps` | Searches a store selection with rgblastn (or a database you built with blastn); says which species and genomes a query comes from |
| `rgblastn` | blastn's flags | A Rust port of blastn's nucleotide engine that searches a store with no database files |

The same verbs and options are available as `gtars blast <verb>`, Python `gtars.blast`, R
`blast_*()` and `refget blast <verb>` ([table](../../docs/blast.md#one-command-set-five-ways-to-call-it);
tutorials for [Python](../../docs/python/blast.md), [R](../../docs/r/blast.md) and
[refget](https://docs.refgenie.org/refget/using-services/blast-store/)). For remote stores,
gtars-refget gained a parallel fetcher and a shared download cache (section 5).

## 3. Build side: from a store to a BLAST database

`rgmakeblastdb build` resolves a selection of collections, builds a `SubjectPlan`, then reads each
selected sequence once, converts it and writes the volumes.

| Topic | FASTA + makeblastdb | rgmakeblastdb |
|---|---|---|
| Sequences shared by assemblies | one OID per FASTA record | one OID per unique sequence, one defline per assembly |
| Ids | whatever the FASTA headers say | fixed scheme: `SQ.<digest>`, plus `gnl\|<db_tag>\|<name>` per assembly |
| Taxids | `-taxid_map`, keyed by id | from FHR, per collection |
| Subsets | rebuild, or `blastdb_aliastool` with estimated statistics | views (`.nal` + `.bsl`) with exact statistics |

### The SubjectPlan

`plan::SubjectPlan` ([src/plan/](../src/plan/)) decides, before any write, which sequences become
OIDs, their order, deflines, taxids and the database totals. The database writer and the rgblastn
engine both use it, so a database search and a store search agree.

- **One OID per unique sequence.** A mitochondrial genome shared by 37 assemblies is stored and
  searched once.
- **Canonical OID order.** Collections sorted by first alias (else digest), each collection's
  sequences in stored order, each sequence taking the next OID where it first appears; sequences in
  no collection last. Order matters: on an exact tie blastn keeps the hit with the higher OID, so
  matching blastn's output means matching its OID order.

### Ids, deflines and titles

| Defline | Ids | Title (example) |
|---|---|---|
| one for the sequence | `SQ.<digest>`, then the sequence's own aliases (INSDC, RefSeq) | `Homo sapiens chr1; 12 assemblies; length:248956422; md5:...` |
| one per (collection, name) | `gnl\|<db_tag>\|<name>`, `gnl\|<other alias>\|<name>`, `gnl\|<collection digest>\|<name>` | `Homo sapiens GRCh38.p13-fasta-genomic (GCF_000001405.39)` |

`db_tag` is the collection's first alias that no other collection in the store uses, else its
digest, so ids do not depend on the selection. The scheme follows from how blastn prints ids:
`sseqid` prints the first id of the first defline, so every hit names its sequence by digest;
`sallseqid` prints the first id of *each* defline, so a hit reads
`SQ.<digest>;gnl|hg38-p14-ucsc|chr1;gnl|...`, one `gnl` id per assembly holding the sequence. That
column alone tells `rgblast report` every genome a hit belongs to. Characters BLAST ids cannot hold
are replaced and listed in `<db>.sanitized.tsv`; integer-looking ids are stored as integers, as
makeblastdb does ([README](../README.md#id-scheme)).

**Taxids** come from each collection's FHR file. Assembly deflines carry their collection's taxid;
the sequence defline carries one only when all its collections agree, so `staxids` lists each
taxon of a shared sequence (`9606;10090`) without a stray `0`. The taxonomy files index every
taxid, so `blastn -taxids` works. `--taxid-map` overrides FHR per collection; makeblastdb's
`-taxid_map` keys on bare ids, so it cannot give human `chr1` and mouse `chr1` different taxids.

### Sequence bytes

Packed codes convert straight to ncbi2na with no ASCII step ([src/convert/](../src/convert/)):

- **2-bit order.** The store's 2-bit code is T=0, C=1, A=2, G=3; ncbi2na is A=0, C=1, G=2, T=3.
  Each 2-bit field `(h, l)` maps to `(!(h^l), !h)`, computed on 64-bit words (32 bases at a time).
- **Ambiguous bases.** makeblastdb writes a random base in place of each N (or other ambiguity
  code) and keeps the real codes in a separate table. The random bases come from NCBI's `CRandom`,
  a lagged Fibonacci generator seeded with the sequence length, one draw per ambiguous base.
  `src/convert/random.rs` ports it, so the bytes match.
- Code tables are derived at compile time from gtars' alphabet definitions, so a change there
  breaks the build rather than the output.

### Byte identity with makeblastdb, and how it is tested

For the same sequences, ids and taxids, `.nhr`, `.nsq`, `.nos`, `.not` and `.nto` are byte-identical
to makeblastdb 2.17's; `.nin`, `.nal` and `.njs` differ only in the date; `.ndb` and `.ntf` hold the
same LMDB keys and values.

- `tests/golden.rs` (feature `golden`) writes a FASTA that `makeblastdb -parse_seqids` turns into
  the same database (deflines joined by Ctrl-A, `src/reference_fasta.rs`), runs makeblastdb in
  docker `ncbi/blast:2.17.0`, and compares file by file.
- The default `cargo test` needs no BLAST+: every id resolves, every sequence read back has the
  store's md5, totals agree. `rgmakeblastdb verify` runs the same checks on a database and its store.
- At scale, every one of 209,486 jungle ids resolved with `blastdbcmd`, and the md5 of all 6,359
  sequences matched the store.

Byte identity is pinned to BLAST+ 2.17.0; a newer BLAST+ that changes bytes but passes the
functional tests is a signal to investigate.

### Subsets and masked assemblies

- **Views** ([src/views.rs](../src/views.rs)): one per assembly, or one for any selection, a `.bsl`
  of the subset's ids plus a `.nal` naming the database, with the subset's exact sequence count and length, so
  E-values equal a fresh build of the subset. (`blastdb_aliastool -seqidlist` writes SeqDB's
  estimated length; see 4a.) No BLAST+ install is needed.
- **Hard-masked assemblies** replace repeats with N (46–56% N in jungle). `nstats` counts N per
  sequence once and caches it by digest; `--exclude-masked` (N fraction 0.25 or more) leaves those
  assemblies out and keeps soft-masked ones, which match their unmasked twins because the store
  ignores case.

### Many small sequences, and remote stores

Transcriptome and microbiome stores hold up to millions of short sequences. The builder gives
workers batches of consecutive OIDs (up to 4,096 sequences or 8 MiB), computes exact header sizes
up front instead of holding headers in memory, and keeps LMDB keys in one arena. GENCODE v47
(384,354 transcripts) builds in about 7 s and 1.6 GB; makeblastdb takes 12.3–12.8 s.

A remote store's selected sequences are first pulled into the shared cache with gtars-refget's
`cache_sequences`, the same steps as `refget store pull` (parallel requests, digest checks). Every
command then reads sequence files, the store's own or the cached ones
([src/store/fetch.rs](../src/store/fetch.rs)). Pulled sequences stay in the cache, so a rerun
downloads only what is missing. Remote builds are byte-identical to local ones. Link speeds vary,
so compare rates: on 2026-11-08 a cold remote hg38 build ran at 7.7 MB/s (177 s), and the same
command again took 1.25 s.

## 4. Query side

### 4a. Stock blastn on a built database

`rgblast search --engine blastn --db BASE` runs `blastn` (from `PATH` or a given path) on a
database you built, as given, appends the columns the report needs (`qlen score sallseqid
salltitles staxids`) and records the run in a `.meta.json`. Three blastn behaviours turned out to
change results silently:

| Behaviour (BLAST+ 2.17) | Effect | Handling |
|---|---|---|
| `-taxids` expands a taxid to its descendants through NCBI's `taxonomy4blast.sqlite3` | Without that file, blastn prints a notice, exits 0 and **searches the whole database** ([NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569839/)) | `rgblast search` refuses `-taxids` unless `-no_taxid_expansion` is given or the file is installed |
| `-dbsize N` | Replaces the database length but not the sequence count: exact for megablast, last-digit E-value differences for some `-task blastn` hits (9 of 115 in plan 5's test) | A wrapper `.nal` with `STATS_NSEQ` and `STATS_TOTLEN` is exact for both tasks; staged origin gives rgblastn the same override |
| Search-time `-seqidlist` or `-taxids` | Statistics use SeqDB's estimated subset length, Σ(4·⌊L/4⌋ + OID mod 4) | Views write the exact length |

### 4b. Where did this sequence come from?

Details are in [docs/origin.md](origin.md).

- **`report`** turns blastn hits into one line per query: the best hit per genome, grouped by
  species, with ties and an optional lowest common ancestor from an NCBI taxdump. It maps each
  `gnl` id in `sallseqid` to its genome through the `collections.tsv` and `members.tsv` sidecars.
- **`origin`** searches and reports in one step, Kraken-style but with alignments. Staged mode
  first searches one representative assembly per species (its deflines still name every genome
  holding each sequence), then the rest of only those species that scored within 20% of the best.
  Every stage uses the full selection's `STATS_NSEQ`/`STATS_TOTLEN`, so E-values equal a full
  search. Staged mode runs rgblastn; with blastn on a database, origin runs one exhaustive search.
  On jungle with blastn megablast: 41.1 s staged, 49.9 s exhaustive, one query answered
  differently.
- **`reps`** shows each species' representative: the one assembly with a `reference` alias, else
  the one sharing the most bases with the rest of its species.

### 4c. rgblastn: searching the store directly

`rgblastn` ports BLAST+ 2.17.0's nucleotide search core (`algo/blast/core`) and the parts of the
C++ layer that decide output: task defaults, the symmetric DUST filter, `CRandom`, query batching
and the tabular formatter ([src/engine/](../src/engine/), about 10,600 lines; each module header
names the C file it ports). It writes no database files. `rgblast search` (its default engine) and
the Python and R bindings call it in process (`api::rgblastn`), with no BLAST+ install.

**What replaces SeqDB.** BLAST+'s engine asks a "sequence source" for subjects; for databases that
is SeqDB, which memory-maps `.nsq` volumes. rgblastn splits the role in three:

| Layer | File | Job |
|---|---|---|
| `SubjectSource` | `engine/subject.rs` | What the engine asks: counts, lengths, statistics overrides, a batch of subjects. The port of BLAST+'s sequence source. |
| `StoreSubjects` | `engine/subject.rs` | Implements it from a `SubjectPlan` and a `SubjectAccess`. |
| `SubjectAccess` | `store/access.rs` | Stored bytes of one OID. `MappedFiles` memory-maps `.seq` files on demand (files under 1 MB are read instead); the search reads a long file piece by piece, as on GPFS every 4 KB page of a mapping costs a page fault. |

For the first search stage, each subject is decoded to ncbi2na with makeblastdb's random fill,
so seeds inside N runs see the bases blastn sees, one 5 Mbp chunk at a time just before the
chunk is searched, so a thread holds one chunk rather than a chromosome. For traceback (the final
gapped alignment) subjects are decoded to blastna with real ambiguity codes, as blastn does; for
subjects of 90,000 bases or more, only the letters around each hit, in a whole-subject buffer of
memory that is never touched elsewhere (SeqDB's partial fetch works the same way).

**Byte-identical to blastn 2.17.0** (tabular output, against a database of the same selection):

| Area | Covered |
|---|---|
| Tasks | `megablast`, `blastn`, `blastn-short` |
| Output | `-outfmt 6` with every nucleotide field (`qseq`, `sseq`, `btop`, `sallseqid`, `staxids`, `qcovs`, ...) |
| Filters | DUST |
| Statistics | the `STATS_NSEQ`/`STATS_TOTLEN` override, `--sequences-file` subsets |

**blastn behaviours reproduced exactly**, each found when output differed:

- **Adaptive query batching.** blastn searches queries in batches sized by how many hits the last
  batch saved (`CBatchSizeMixer`); each batch is a separate search with its own lookup table, and
  under `-task blastn` a query that is not first in its batch can gain or lose hits at its first
  base. The port copies the mixer member for member, including a 32-bit cast that overflows when a
  batch saves no hits. On chr1 with 10 queries, all 1,037,968 rows match blastn.
- **The megablast x-drop default.** The C core's default ungapped x-drop is 20, but blastn's
  megablast setup never sets it, so it is 0 and the cutoff score is used. With 20, 497 hg38 HSPs took
  a different traceback seed.
- **The purge order.** When two HSPs share an endpoint, blastn's traceback cuts one back and leaves
  the cut HSPs in reverse order of moving; later stable sorts depend on that order
  (`engine/hsp.rs`).
- Also: the higher-OID tie-break, 5 Mbp subject chunks with 100-base overlap, and floating-point
  operations in the C's order.

**Threads** split the database, as blastn's default `-mt_mode 0` does. Work items (runs of
consecutive OIDs; a subject over 5 Mbp alone, as one blastn thread searches it) come from a shared
counter, longest first; the hit collector does not depend on arrival order, and traceback results
are placed by index. `tests/threads.rs` shuffles work items and varies threads and batch sizes
without changing a byte.

**Remote stores.** rgblastn pulls the selected sequences into the shared cache, then maps the
cached files as for a local store; a rerun downloads only what is missing.
`remote_rgblastn_equals_local` checks that a remote search gives the local bytes.

**Speed** (docker `ncbi/blast:2.17.0` against rgblastn, same machine, shared with other jobs;
[benchmarks/results](../benchmarks/results/README.md)):

| Search (megablast unless noted) | blastn | rgblastn |
|---|---|---|
| hg38, 10 queries, 1 thread | 42.5–45.5 s | 41.9–47.2 s |
| hg38, 10 queries, 4 threads | 24.6–27.2 s | 15.1–17.0 s |
| jungle (41.8 Gbp), 20 queries, 16 threads | 60.2–60.8 s | 45.3–45.6 s |
| GENCODE (384,354 subjects), 1 thread | 2.2–2.3 s | 7.4–7.6 s |
| GENCODE, `-task blastn`, 4 threads | 9.6–10.0 s | 14.3–14.8 s |

**Checks.** `tests/golden_search*.rs` (feature `golden`) compare full output with docker blastn.
At scale, hg38 megablast (495,970 rows) and jungle megablast (809,785 rows) are byte-identical.
hg38 `-task blastn` with all 10 queries in one run has no blastn reference: blastn ran out of
memory at 18 GB.

## 5. Changes to existing code

**gtars-refget** (plan 3; new [`store/fetch.rs`](../../gtars-refget/src/store/fetch.rs),
[`store/cache.rs`](../../gtars-refget/src/store/cache.rs),
[`store/prefetch.rs`](../../gtars-refget/src/store/prefetch.rs)):

- **Parallel fetch:** keep-alive connections (32 by default), sequence files over 8 MiB as
  parallel 8 MiB ranges (the server must support HTTP Range requests), retries with jittered
  backoff, every download checked against its digest, metadata fetched in parallel batches.
- **Shared, content-addressed cache:** one cache root (`--cache-dir`, else `$REFGET_CACHE_DIR`,
  else `$XDG_CACHE_HOME/refget`, else `~/.cache/refget`) replaces each store's `.remote_cache/`.
  Sequences live under `sequences/<alphabet>/<d0d1>/<digest>.seq`, so a sequence held by two
  stores downloads once. Cached metadata is refetched when the store's manifest digests change.
- **Bulk pull:** `prefetch` streams many sequences in parallel and `cache_sequences` downloads
  them into the cache; gtars-rgblast pulls a remote store's selection with `cache_sequences`. A
  sequence is published to the cache only when complete, so a rerun downloads only what is
  missing.
- **Three fetch flags**, the same on every BLAST verb that reads a store and on FASTA export:
  `--cache-dir`, `--resume` (always on for BLAST: pulled sequences stay in the cache) and
  `--fetch-jobs`.
- `sequence_relative_path` (`store/readonly.rs`) is now public and refuses absolute paths and `..`,
  so it is safe to join onto a local root or a remote URL.

**gtars CLI and bindings:** `gtars blast <verb>` (`gtars-cli/src/blast/`, feature `blast`, on by
default) mounts `gtars_rgblast::cli::BlastCommand`, the clap parser the `rgmakeblastdb` and
`rgblast` binaries use, so flags, errors and exit codes match. Python `gtars.blast`
(`gtars-python/py_src/gtars/blast/`) and R `blast_*()` (`gtars-r/R/blast.R`) turn keyword options
into the verb's flags and call `gtars_rgblast::cli::run_json`, so they parse with the same clap
structs; they add progress callbacks and Ctrl-C cancel.
`gtars refget export` gained `--cache-dir`, `--resume` and `--fetch-jobs`, and the Python and R
`open_remote` gained `fetch_jobs`.

**refget CLI:** `refget blast <verb>` (`refget/cli/blast.py`) passes its arguments to the same
parser (`gtars.blast._main`), mapping `-p/--path` or `-r/--remote` to `--store`; its help and exit
codes are clap's. In `refget/cli/store.py`, `export` and `pull` gained `--cache-dir` and
`--fetch-jobs` (and `export` `--resume`), and `refget store cache info|clear` is new.

## 6. Differences from BLAST+

| Topic | BLAST+ | gtars-rgblast |
|---|---|---|
| Identical sequences | separate OIDs | one OID with several deflines; E-values equal a makeblastdb database when its FASTA has no duplicate sequences |
| `.nog` (OID to GI) | written | not written; BLAST+ does not need it |
| Taxonomy names (`sscinames`) | NCBI's taxdb | install NCBI's taxdb; rgblastn refuses those fields |
| Sequences of 2^31 bp or more | makeblastdb crashes, refuses or wraps the length, depending on size | skipped and listed in `<db>.skipped.tsv` |
| Proteins | protein databases | skipped; nucleotide only |
| `U` | n/a | stored as `T` |

rgblastn covers the common searches: tasks `megablast`, `blastn` and `blastn-short` with
`-outfmt 6`. It refuses any other blastn option by name and points to building a database with
`rgmakeblastdb` and running `blastn` on it. It searches each query in one piece, up to the task's
query chunk (1 Mbp for `blastn` and `blastn-short`, 5 Mbp for `megablast`), and refuses a longer
query; blastn would split it (`split_query.c`).

## 7. Open work

- **Packed stores** (plan 8): millions of small sequences packed into a few large files instead of
  one file each. The build-side batching is done; the packed layout itself waits for gtars-refget's
  packed-layout work to reach `dev`. Until then `store::access::open` refuses packed stores.
- **Transcriptome search speed.** On GENCODE, rgblastn is about 3x slower than blastn for megablast
  and 1.3–1.5x slower for `-task blastn`. Most of it is per-subject cost: 4.7 s of plan setup
  (mostly freeing per-subject strings and maps) and one file open per subject. Plan 8 replaces
  both. The rest of the `-task blastn` gap is the dynamic-programming traceback.
- **CI.** The BLAST+ steps and the `blast` job in `.github/workflows/tests.yml` and `codecov.yml`
  have not been run in CI.
- **Not yet measured:** jungle with `-task blastn`; remote rgblastn speed.
- **Status.** In the metaplan, plans 1–7 are marked implemented and plan 8 (packed stores) is
  open; none is marked verified yet.

## 8. Credit

The search algorithms, DUST masker, `CRandom` and database format come from NCBI BLAST+ 2.17.0
([source](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/)), a United States
Government work in the public domain. [NOTICE](../NOTICE) reproduces NCBI's notice and lists the
ported parts. Please cite:

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture
and applications. *BMC Bioinformatics* 10, 421 (2009). <https://doi.org/10.1186/1471-2105-10-421>
