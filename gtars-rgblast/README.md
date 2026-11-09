# gtars-rgblast

BLAST nucleotide databases built straight from a [RefgetStore](../gtars-refget), with no FASTA
in between. `rgmakeblastdb` writes BLAST+ version 5 databases that `blastn` and `blastdbcmd`
read like any makeblastdb database: hg38 in about 2 s and a 79-assembly store (41.8 Gbp) in
about 20 s on 16 cores, against minutes for FASTA export plus `makeblastdb`.

## Quick start

```bash
# one assembly, by alias or digest
rgmakeblastdb build --store /data/stores/jungle --collection hg38-primary-refgenie --out db/hg38
# every assembly in a store (each unique sequence once)
rgmakeblastdb build --store /data/stores/jungle --all --out db/jungle
# all unmasked mouse assemblies, then one view per assembly
rgmakeblastdb build --store /data/stores/jungle --taxon "Mus musculus" --exclude-masked --out db/mouse
rgmakeblastdb view --db db/mouse --per-collection
# a remote store (any plain HTTP prefix, such as the public S3 stores)
rgmakeblastdb build --store https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle \
  --collection EiFob05aCWgVU_B_Ae0cypnQut3cxUP1 --out db/hg38
# see what a build would do, then search
rgmakeblastdb build --store /data/stores/jungle --all --out db/jungle --dry-run
blastn -db db/jungle -query query.fa -outfmt "6 qseqid sseqid pident length evalue bitscore sallseqid"
# check a database against its store (no BLAST+ needed)
rgmakeblastdb verify --db db/jungle --store /data/stores/jungle
```

`rgmakeblastdb build --help` lists every option. Existing databases are never replaced without
`--force`; a failed or cancelled build (Ctrl-C) removes everything it wrote.

Every `--store` also takes an `http(s)` URL. When a command needs the sequences, it first pulls
the selected ones into the shared cache with gtars-refget, the same steps as `refget store pull` (progress line
`pulling N sequences (X bp) into the cache`), then reads the cached files as it reads a local
store. Each sequence is checked against its digest. `--fetch-jobs N` sets the parallel requests
(default 32). The cache root is `--cache-dir DIR`, else `$REFGET_CACHE_DIR`, else
`$XDG_CACHE_HOME/refget`, else `~/.cache/refget` (the same root the gtars and refget CLIs use).
Pulled sequences stay in the cache, so a rerun downloads only what is missing (`--resume` is
accepted and always on). `--verbose` prints the pull's network line (bytes, requests, retries,
MB/s, cache hits).

From Rust, `gtars_rgblast::api::{preflight, build, verify, view, nstats}` take plain option
structs and report progress through a callback; `gtars_rgblast::progress::describe` gives the
lines the binary prints. The clap argument structs (`gtars_rgblast::cli::{BuildArgs,
VerifyArgs, ViewArgs, NstatsArgs, SelectionArgs}`) convert to those requests with `TryFrom`,
which re-checks every rule. `gtars_rgblast::cli::run_json(argv, progress)` runs any
verb from an argument list and returns its result as JSON; the Python and R bindings call it.

## Searching: `rgblast`

`rgblast` searches a store without a hand-built database and says where each query comes from.
See [docs/origin.md](docs/origin.md) for every option, output column and rule.

```bash
# search the selection with rgblastn, in process (no database)
rgblast search --store /data/stores/jungle --all -- -query q.fa -outfmt 6 -out hits.tsv
# best hit per genome, grouped by species, one line per query (reads hits.tsv.meta.json)
rgblast report hits.tsv --query q.fa
# both in one step, searching one representative per species first
rgblast origin --store /data/stores/jungle --all --exclude-masked --query q.fa
# stock blastn on a database you built (rgmakeblastdb build or view), searched as built
rgblast search --engine blastn --db db/jungle -- -query q.fa -outfmt 6 -out hits.tsv
```

The default engine is rgblastn (below), on the store and selection. `--engine blastn --db BASE`
runs stock blastn on a database you built; select collections when building it (or build a
view with `rgmakeblastdb view`), since it is searched as given. `origin` with blastn runs one
exhaustive search; staged origin needs rgblastn. blastn comes from `--blastn PATH`,
`$GTARS_RGBLAST_BLASTN`, or `PATH`. From Rust, `gtars_rgblast::api::{search, report, origin,
reps}` take the same requests the `rgblast` argument structs convert to.

## Searching without any database: `rgblastn`

`rgblastn` is a Rust port of NCBI BLASTN 2.17.0's nucleotide engine that searches the store
directly, local or remote, for common searches. It takes blastn's flags (single-dash spelling)
plus the store and selection flags, and its `-outfmt 6` output is byte-identical to `blastn`
against a database built from the same selection. It supports:

- `-task megablast` (the default), `blastn` and `blastn-short`;
- `-outfmt 6` with any field list (`std`, `qseq`, `sseq`, `btop`, `sallseqid`, `salltitles`,
  `sallacc`, `staxid(s)`, `qcovs`, ...; the taxonomy *name* fields need NCBI's taxdb and are
  refused);
- `-evalue`, `-max_target_seqs`, `-max_hsps`, `-perc_identity`, `-strand`, `-word_size`,
  `-dust`, `-soft_masking` and `-num_threads` (the output does not depend on the thread count);
- `--sequences-file` to search a subset with the statistics of a database of just those
  subjects, and `--stats-nseq N --stats-totlen L` to search it under another database's
  statistics (what a wrapper alias with `STATS_NSEQ`/`STATS_TOTLEN` does).

rgblastn searches each query in one piece, up to blastn's query chunk (1 Mbp for `blastn` and
`blastn-short`, 5 Mbp for `megablast`). Any other blastn option is refused by name. For a
longer query or another option, build a database with `rgmakeblastdb` and run `blastn` on it (or
`rgblast search --engine blastn --db BASE`). `rgblast search` runs the same engine in process.

```bash
rgblastn --store /data/stores/jungle --collection hg38-primary-refgenie \
  -query q.fa -task megablast -evalue 1e-5 -outfmt "6 std stitle" -out hits.tsv
```

From Rust: `gtars_rgblast::api::rgblastn::rgblastn(&RgblastnRequest, progress)` runs the search
in process and `api::rgblastn::write_tabular` (over `gtars_rgblast::format::tabular::write`)
renders it; `engine::search::search` is the engine itself. The algorithms are NCBI BLAST+
2.17.0's, a US Government work in the public domain (see `NOTICE` for NCBI's notice and the
citation to use; each module header under `src/engine/` names the C file it ports);
`tests/golden_search.rs` and `tests/golden_search_blastn.rs` (feature `golden`, docker
`ncbi/blast:2.17.0`) pin the output bytes, and
`scripts/golden/search_matrix.py` runs the flag matrix on a real store.

## Choosing collections

Every command that reads a store takes the same selection flags. Choose collections in one of
three ways:

- **A list:** `--collection SEL` (repeatable) and `--collections-file FILE` (one selector per
  line; `-` is stdin; `#` comments and blank lines are skipped; only the first field counts).
  A selector is a collection digest; `NS:ALIAS` when `NS` is one of the store's collection
  alias namespaces (`name:hg38-p14-ucsc`, `refseq:GCF_000001405.39`); or a bare alias that only
  one collection uses. A collection listed twice is kept once.
- **By attribute:** `--namespace NS` (collections with an alias in `NS`) and `--taxon TAXID|NAME`
  (an NCBI taxid, or a species name from FHR, matched ignoring case, with `_` as a space). Values
  of one flag are OR'd; `--namespace` and `--taxon` together are AND'd.
- **Everything:** `--all`, every collection plus the sequences in no collection.

Two filters work with any mode:

- `--max-n-fraction F` (0 < F ≤ 1) leaves out collections whose N fraction (N bases over all
  bases of the collection's sequences) is F or more. `--exclude-masked` means exactly
  `--max-n-fraction 0.25`: on the jungle store it leaves out the 12 hard-masked assemblies
  (46–56% N) and keeps the soft-masked ones, which are identical to their unmasked twins
  because the store ignores case. In a list, a collection that fails the filter is an error, not
  a silent drop. With `--all`, sequences that belong to no collection (orphans) are kept, but not
  when the filter excluded something, because an excluded collection's sequences would come back
  as orphans.
- `--sequences-file FILE` keeps only the listed sequences (digests, `SQ.` optional), each with
  the deflines of every selected collection that holds it. A digest that no selected collection
  holds is an error (under `--all`, sequences in no collection count as selected).

| Flags | Allowed with | Error with |
|---|---|---|
| `--collection`, `--collections-file` | each other, an N filter, `--sequences-file` | `--all`, `--namespace`, `--taxon` |
| `--namespace`, `--taxon` | each other, an N filter, `--sequences-file` | `--all`, a list |
| `--all` | an N filter, `--sequences-file` | a list, `--namespace`, `--taxon` |
| `--max-n-fraction` | any mode | `--exclude-masked` |
| `--exclude-masked` | any mode | `--max-n-fraction` |

Choosing no mode is an error too. Every violation exits with code 2 and clap's wording.

N counts are computed once per sequence and cached by digest in `<cache root>/blast/nstats/`
(the shared cache root: `--cache-dir`, else `$REFGET_CACHE_DIR`, else `~/.cache/refget`), so they
are valid for any store. Every build adds the
counts of the sequences it converts, so after one full build, filtered builds read nothing extra.
`rgmakeblastdb nstats --store S [selection]` fills the cache and prints each collection's N
fraction, highest first, to help choose a threshold.

## Id scheme

Every database uses one scheme, whatever the selection:

- **One OID per unique sequence.** A sequence shared by many assemblies is stored and searched
  once. OIDs follow the canonical order: collections sorted by their first alias (else digest),
  each collection's sequences in their stored order, a sequence taking the next OID where it
  first appears; with `--all`, sequences in no collection come last.
- **Alias order.** A collection's aliases are ranked by namespace: `name` first, then the other
  namespaces in `rgstore.json` order, then the accession namespaces `insdc` and `refseq`; within a
  namespace, in byte order. The first alias is the collection's display name.
- **Defline 0** of each OID has the id `SQ.<sha512t24u digest>`, then every sequence alias
  (INSDC, RefSeq) as a local id. Its title is `<species> <name>; <n> assemblies; length:<len>;
  md5:<md5>`, for example `Homo sapiens chr1; 12 assemblies; length:248956422; md5:…`.
- **One more defline per (collection, name)**, with ids `gnl|<db_tag>|<name>`, then
  `gnl|<alias>|<name>` for the collection's other aliases, then `gnl|<collection digest>|<name>`.
  Its title is `<species> <display name>`, then ` (<assembly accession>)` when known, then the
  FASTA description if there is one: `Homo sapiens GRCh38.p13-fasta-genomic (GCF_000001405.39)`.
  The name is not repeated (it is already the id), so `-outfmt 0` prints
  `hg38-p14-ucsc:chr1 Homo sapiens hg38-p14-ucsc`, and `stitle`/`salltitles` still name the
  assembly.
- **`db_tag`** is the collection's first alias that no other collection in the store uses
  (compared ignoring case), else its digest. It depends only on the store, not on the selection.
  An alias used by more than one collection is left out of the ids (and reported).

blastn's `sallseqid` prints only the first id of each defline, so a hit shows as
`SQ.<digest>;gnl|<db_tag>|<name>;…`, one `gnl` id per assembly that contains the sequence.

### Id rules

- **Sanitizing.** Characters BLAST ids can't hold (space `>` `[` `]` `|` `"` `;`, tab and control
  characters) become `_` in ids. The original name stays in the title and in `<db>.members.tsv`;
  every change is listed in `<db>.sanitized.tsv`. If two names of one collection become the same
  id (for example `a|b` and `a_b`, or names differing only in case), the build stops before
  writing anything and names both.
- **Integer ids.** makeblastdb stores an id that is a canonical positive decimal up to
  2,147,483,647 (no leading zero, no sign) as an integer. `rgmakeblastdb` does the same, for local
  ids and `gnl` tags alike, so the bytes match. Checked against makeblastdb 2.17: `1`, `22`,
  `123456` and `2147483647` are integers; `0`, `00`, `007`, `0099`, `2147483648`, `+1` and `-1`
  are strings. Lookups work the same either way.
- **Digests starting with `-`.** A `gnl|-…|name` id works everywhere. Its other spelling,
  `-…:name`, resolves inside BLAST but looks like a flag on a command line, so `<db>.ids.tsv`
  lists only the `gnl|` form.

## Looking up ids

```bash
blastdbcmd -db db/jungle -entry 'SQ.IefZqMQcwYR9aHjkjMlErB-KyzfUV9ae'
blastdbcmd -db db/jungle -entry 'gnl|hg38-primary-refgenie|chrM'
blastdbcmd -db db/jungle -entry_batch ids.txt -outfmt '%o %a %l'
```

Every id is also stored upper-cased, so lookups ignore case.

## Sidecars

Written next to the database (header-only when empty), not listed in `.njs`. `ids.tsv` and
`members.tsv` take about 370 bytes per sequence (140 MB for the 384,354 GENCODE transcripts).

| File | Rows |
|---|---|
| `<db>.ids.tsv` | `id`, `oid`, `sequence_digest`: every id that is safe to paste into `blastdbcmd -entry` |
| `<db>.collections.tsv` | one per selected collection: digest, `db_tag`, display name, aliases, taxid, species, assembly accession, sequence count, length and N fraction (an empty field when unknown) |
| `<db>.members.tsv` | one per (collection, name) in the database: `oid`, digest, collection, original name, id name |
| `<db>.oids.tsv` | one per OID, in order: digest, length, N count, ambiguous count |
| `<db>.manifest.json` | the store and its digests, the selection, the taxid map, the volumes and files, totals |
| `<db>.sanitized.tsv` | every name or alias that was changed for use in ids |
| `<db>.skipped.tsv` | sequences left out: `oversize`, `non_nucleotide` (proteins) or `empty` |

`gtars_rgblast::blastdb::sidecars::{read_collections, read_members, read_ids, read_oids,
read_manifest}` read them back.

## Taxonomy

Each collection's taxid and species come from its FHR metadata (`fhr/<digest>.fhr.json`:
`taxon.uri`, such as `https://identifiers.org/taxonomy:9606`, and `taxon.name`). Its assembly
accession is FHR `accessionID.name`, else its first `refseq` alias, else its first `insdc` alias.
A collection without an FHR taxon gets taxid 0, BLAST's "unknown", and the build says which.

- **`--taxid-map FILE`** overrides FHR, one collection per line:
  `<collection selector><TAB><taxid>[<TAB><species name>]`. Lines for collections outside the
  selection are ignored, so one map serves many builds, and the map is recorded in the manifest.
  This is a per-collection map: makeblastdb's `-taxid_map` looks ids up by their bare
  tag (`chr1`), so it cannot give human `chr1` and mouse `chr1` different taxids.
- **Every collection defline** carries its collection's taxid. The `SQ.` defline carries the
  taxid all of its collections agree on; when they disagree it carries none, so `staxids` lists
  each taxon of a shared sequence (`9606;10090`) without a stray `0`.
- **The taxonomy files** (`.not`, `.nto`, `.ntf`) index every taxid, so `-taxids` filters.
- **The `-taxids` trap.** `blastn -taxids T` expands T to its descendant taxa through NCBI's
  `taxonomy4blast.sqlite3`. When that file is in neither the working directory nor a `BLASTDB`
  directory, blastn 2.17 prints a notice, exits 0, and **searches the whole database**. Always
  add `-no_taxid_expansion` unless the file is installed:
  `blastn -db db/jungle -query q.fa -taxids 10090 -no_taxid_expansion`. Filtering keeps hits on
  sequences with at least one defline of the taxon and shows only those deflines;
  `-negative_taxids` drops them. `-taxids` cannot be combined with `-seqidlist`.
- **Species names in output** (`sscinames`, `scomnames`, `sskingdoms`) and descendant taxa need
  NCBI's taxdb: download https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz (66 MB:
  `taxdb.btd`, `taxdb.bti`, `taxonomy4blast.sqlite3`) into a directory on `BLASTDB` (or run
  `update_blastdb.pl --decompress taxdb` there). BLAST does not look for these files next to the
  database, so rgmakeblastdb does not write its own.
- **Lineage and LCA** come from an NCBI taxdump in `rgblast report`, not from FHR, which has only
  the taxid and species name.

## Views

A view searches part of a database with no rebuild:

```bash
rgmakeblastdb view --db db/jungle --per-collection             # one per assembly
rgmakeblastdb view --db db/jungle --name mouse --taxon 10090   # one for a species
rgmakeblastdb view --db db/jungle --name unmasked --all --exclude-masked
blastn -db db/jungle.views/hg38-p14-ucsc -query q.fa
blastn -db db/jungle.views/mouse -query q.fa
```

Each view is a `.nal` naming the database plus a `.bsl` seqid list of the member rows' collection
ids, written to `<base>.views/` (`--out-dir` to change), with an index in `views.tsv`. Selection
flags limit which collections get views, with N fractions taken from `collections.tsv`; no store
is read. `--per-collection` writes one view per selected collection; `--name NAME` writes one
view of the whole selection.

- **Exact statistics.** A view's `.nal` holds the exact sequence count and length of the subset,
  so `blastdbcmd -info`, the blastn footer and every E-value equal a fresh build of the same
  selection (and makeblastdb over the assembly's FASTA when it holds no duplicate sequence).
  Hits show the collection's own ids. The order of tied hits can differ from a fresh build.
- **`--full-db-stats`** also writes the whole database's sequence count and length as the
  statistics, so E-values are on the scale of the whole database while `-info` still describes
  the view. Use it, not `-dbsize`, which is exact only for megablast.
- **Why not search-time `-seqidlist` or `-taxids`?** They work, but BLAST computes their
  statistics from an approximate subset length (`4*floor(L/4) + (OID mod 4)` per sequence), so
  E-values differ slightly from a database of the subset; `blastdb_aliastool -seqidlist` writes
  that approximate length into its `.nal` too. A view of the taxon (`--name NAME --taxon T`)
  gives exact statistics; don't combine `-taxids` with a view.

## Many small sequences

Transcriptomes and microbiome stores hold hundreds of thousands to millions of short sequences.
`rgmakeblastdb` hands its workers batches of consecutive OIDs (up to 4,096 sequences or 8 MiB
of stored bytes), each decoded into one reused buffer, so scheduling
costs nothing per sequence; the output does not depend on the batch size. Headers and lookup keys
are not held in memory before writing: the preflight computes their exact sizes, and the writer
encodes each header as it writes it and keeps the lookup keys in one byte arena. GENCODE v47
(384,354 transcripts, a standalone store) builds in about 7 s with 16 threads and 1.6 GB of
memory; makeblastdb takes 12.5 s on the same sequences
(`benchmarks/results/bench_txome.jsonl`, `scripts/bench/txome_store.sh`).

## Limits

- **Sequences of 2^31 bp or more** can't be stored in a BLAST database: lengths are 32-bit
  signed integers throughout BLAST+ (the writer, the `.nin` index and the search engine).
  `rgmakeblastdb` skips them by default and lists them in `<db>.skipped.tsv`;
  `--on-oversize fail` stops the build instead. makeblastdb 2.17 can't hold them either
  (measured 2026-11-08): at exactly 2^31 − 1 bp it
  crashes (`basic_string::_M_create`, exit 139, no database); at 2^31 + 10 bp it stops with
  "Sequences longer than 2,147,483,647 bases are not supported"; and at 2^32 + 10 bp it silently
  writes a database holding a 10-base sequence (the length wrapped modulo 2^32).
- **Proteins** and other non-nucleotide sequences are skipped (`non_nucleotide`).
- **U** is stored as T (BLAST has no U); `verify` accounts for that.
- makeblastdb also writes a `.nog` (OID to GI) file per volume; `rgmakeblastdb` does not, and
  BLAST+ does not need it.

## Tests

- `cargo test -p gtars-rgblast` needs no network, docker or BLAST+. On small committed fixtures
  (`tests/data/fixtures`) it checks the md5 of every entry, that every id resolves, the totals,
  identical output across storage modes and local/remote stores, and the error paths. For
  `rgblast` it checks the report on recorded jungle hits, `origin`'s stages by replaying recorded
  blastn output (`origin_store_hits/`) through a fake runner, the representative rule on a
  jungle-shaped membership table, and the heartbeat and cancel of a fake blastn.
- `cargo test -p gtars-rgblast --features blastplus` runs the same checks through real BLAST+
  tools, and compares blastn tabular output with a makeblastdb database built from the same
  sequences. `blastplus_taxonomy` checks `-taxids` (and the trap), `staxids`, `-outfmt 0` and
  view statistics against fresh builds; its `sscinames_with_taxdb` test runs when
  `GTARS_RGBLAST_TAXDB` names an unpacked NCBI taxdb.
- `cargo test -p gtars-rgblast --features golden` adds the byte-identity canary, pinned to BLAST+
  2.17.0. A newer BLAST+ that breaks byte identity but passes `blastplus` is a canary to
  investigate, not a broken product. `-- --ignored golden_jungle` runs it on a local jungle
  mirror (`GTARS_RGBLAST_STORE_JUNGLE`, about 20 min), and `-- --ignored golden_txome_gencode` on
  the GENCODE transcript store of `scripts/bench/txome_store.sh` (`GTARS_RGBLAST_STORE_TXOME`).
  The `txome` fixture (2,826 short sequences in 42 collections, generated at test time by
  `tests/common/txome.rs` and pinned by `tests/data/fixtures/txome.expected.json`) is in the
  default suite and the canary. `golden_origin` runs staged `origin` with rgblastn, and
  `search`, `origin` and `report` with real blastn on a database built from a synthetic store,
  against the recorded `origin_store_hits/`, and checks the BLAST behaviors `rgblast` relies on;
  `GOLDEN_UPDATE=1` re-records `origin_store_hits/`.

BLAST+ comes from `GTARS_RGBLAST_BIN` (a directory), `PATH`, or the docker image
`GTARS_RGBLAST_DOCKER` (default `ncbi/blast:2.17.0`; tests that run blastn through `rgblast`
pass it a test-only shim script that calls docker). With a feature on and no BLAST+, the tests
fail and say how to provide it.

## Scripts

- `scripts/golden/`: manual checks at scale. `remote_verify.sh` builds hg38 from the S3 store,
  compares it byte for byte with a local build and checks every md5; `search_matrix.py` runs the
  rgblastn flag matrix against blastn on a real store. Helpers: `compare.py`, `md5check.py` (the
  Python scripts run with `uv run --script`).
- `scripts/bench/`: benchmarks (`bench_hg38.sh`, `bench_jungle.sh`, and `txome_store.sh` for
  GENCODE transcripts), run through `bench.py` and summarized by `summarize.py`; results are in
  `benchmarks/results/`.
- `scripts/make_origin_fixtures.sh`: regenerates the `rgblast report` fixtures
  (`tests/data/fixtures/origin/`) from a jungle mirror; `scripts/jungle_queries.sh` writes the
  jungle query set.

Paths come from `WORK` (default `~/work/gtars-rgblast`) and `STORES` (default
`~/work/refget-blast/stores`); the binary from `B` (default `$WORK/target/release/rgmakeblastdb`).
