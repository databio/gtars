# Where did this sequence come from? `rgblast search`, `report` and `origin`

`rgblast` answers "which species, and which genomes, could this DNA sequence have come from?"
against a [RefgetStore](../../gtars-refget) such as jungle (79 human and mouse assemblies) or VGP
(605 vertebrate assemblies), with the in-process `rgblastn` engine on the store, or stock
[BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) `blastn` on a database you built.

- `rgblast search` searches a selection with rgblastn (or a database with blastn), adding the
  output columns the report needs.
- `rgblast report` turns tabular hits into the best hit per genome, grouped by species, one line
  per query.
- `rgblast origin` does both, searching one representative assembly per species first and then
  the rest of the species that matched.
- `rgblast reps` shows the representatives and the predicted cost.

## `rgblast search`

```bash
rgblast search --store /data/stores/jungle --all --exclude-masked \
    -- -query q.fa -outfmt 6 -evalue 1e-20 -out hits.tsv
```

Selection flags are the same as `rgmakeblastdb` (see the [README](../README.md#choosing-collections)).
Everything after `--`, or from the first single-dash flag, goes to the engine.

- **Engine.** The default, rgblastn, searches the store and selection in process with no
  database (megablast, blastn and blastn-short; see the [README](../README.md)); a flag it does
  not implement is refused by name. `--engine blastn --db BASE` runs stock blastn on a database
  you built with `rgmakeblastdb build` or `view`, searched as given: select collections when
  building it, so `--store` and selection flags (other than `--all`) are refused. blastn comes
  from `--blastn PATH`, `$GTARS_RGBLAST_BLASTN`, then `blastn` on `PATH`; BLAST+ 2.10 or newer is
  required, and tests pin 2.17.0. With either engine, `-num_threads` defaults to all cores.
- **Output columns.** For `-outfmt 6` or `7`, `qlen score sallseqid salltitles staxids` are
  appended when missing (`-outfmt 6` becomes `6 std qlen score sallseqid salltitles staxids`), so
  your own columns keep their positions. `--no-extra-fields` turns this off.
- **Refused arguments.** `-db` (give `--db BASE` with `--engine blastn`), `-subject`,
  `-subject_loc` and `-remote`. With blastn, `-taxids`/`-taxidlist` need `-no_taxid_expansion`
  unless NCBI's `taxonomy4blast.sqlite3` is installed, because without it blastn 2.17 searches
  the whole database and only prints a notice.
- **Metadata.** `--meta FILE`, or `FILE.meta.json` next to `-out FILE`, records the engine, the
  store, its digests and the selection's collections (rgblastn) or the database (blastn), and the
  final arguments. `rgblast report` reads it. The exit status is the engine's.

## `rgblast report`

```bash
rgblast report hits.tsv --query q.fa          # uses hits.tsv.meta.json
rgblast report hits.tsv --db db/jungle        # any rgmakeblastdb database
```

Without `--db`, the report reads the database the meta file names (blastn), else the store and
collections it names (rgblastn).

Hits are `-outfmt 7` (columns from its `# Fields:` lines) or `-outfmt 6` (columns from
`--fields`, default `std qlen score sallseqid salltitles staxids`). The report needs a query id,
`pident length qstart qend sstart send evalue bitscore qlen sallseqid`; a missing column is an
error that names it.

**Fan-out.** blastn's `sallseqid` prints the first id of every defline: `SQ.<digest>` for the
sequence, then `gnl|<db_tag>|<name>` for each genome holding it (rgblastn prints the same ids).
The report maps each `gnl` id to its genome and original name through the database's
`collections.tsv` and `members.tsv` (or the same tables built from the store). When a
seqid list kept only `SQ.` ids, the genomes come from the digest instead.

### Output

`--format tsv` (default), one line per query in input order. Real lines from the jungle
sequences up to 200 kb:

```
#query	status	taxid	taxon	rank	top_genome	top_seq	pident	qcov	bits	evalue	tied_genomes	tied_taxa	searched	flags
q_human_chrM_2001	taxon	9606	Homo sapiens	-	GRCh37-igenomes-ensembl	MT	100.000	100.0	1847	0.0	37	1	79/79	-
q_mouse_chrM_5001	taxon	10090	Mus musculus	-	GRCm38-igenomes-ensembl	MT	100.000	100.0	1847	0.0	27	1	79/79	-
q_ebv_100001	taxon	9606	Homo sapiens	-	GRCh37.p13-fasta-full-analysis	chrEBV	100.000	100.0	1847	0.0	13	1	79/79	-
q_random	none	-	-	-	-	-	-	-	-	-	0	0	79/79	-
```

| Column | Meaning |
|---|---|
| `status` | `taxon` (all tied genomes in one group), `lca` (several groups, called as their lowest common ancestor), `multi` (several groups, no taxonomy loaded), `none` (no hit survives the filters) |
| `taxid`, `taxon`, `rank` | the call; `rank` is filled only with a taxonomy. A genome without an FHR taxid is its own group, `unknown:<name>`; `multi` lists the groups, most tied genomes first |
| `top_genome`, `top_seq` | the first tied genome by identity, then alignment length, then name, and its sequence name |
| `pident`, `qcov`, `bits`, `evalue` | that genome's best HSP, printed as blastn printed it; `qcov` is the HSP's query coverage |
| `tied_genomes`, `tied_taxa` | genomes within `--tie-percent` of the top score, and their groups |
| `searched` | genomes whose every sequence the query was compared with, over the selection size |
| `flags` | `truncated` (the query hit `-max_target_seqs` subjects, so more genomes may exist), `reps-only` (origin, below) |

`--format json` is one document with every query, its call, and every hit genome grouped by
taxon. `--detail FILE` writes one row per (query, genome) with the sequence, coordinates,
strand, identity, E-value, whether it ties, and which stage found it (`--detail-top N` keeps
the best N genomes per taxon). `--output FILE` writes the summary to a file.

### Ties and calls

- HSPs failing `--min-pident`, `--min-qcov` or `--max-evalue` are dropped first (none by default).
- Each genome keeps its best HSP: highest bit score, then lowest E-value, highest identity,
  longest alignment, smallest subject start.
- The top score is the raw `score` when every hit has it (exact ties for any scoring scheme),
  else the bit score. Genomes scoring at least `top x (1 - P/100)` tie, with P =
  `--tie-percent` (default 0, exact ties only). P plays the role of
  [MEGAN](https://software-ab.cs.uni-tuebingen.de/download/megan6/manual.pdf)'s top percent.
- Groups are taxids. With `--group-rank RANK` and a taxonomy, genomes are grouped by their
  ancestor at that rank (`--group-rank species` folds subspecies into their species).

### Taxonomy

FHR metadata names each genome's species and taxid but no lineage, so LCA calls use the NCBI
taxdump, the same source as Kraken2, mmseqs and sourmash lineages. It is optional and never
downloaded:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz nodes.dmp names.dmp merged.dmp
rgblast report hits.tsv --taxdump .        # or export GTARS_RGBLAST_TAXDUMP=$PWD
```

Merged taxids are followed; an unknown taxid gets a warning and becomes a leaf under root. Only
the names that appear in the output are loaded. Tied human and mouse genomes are called
`lca 314146 Euarchontoglires superorder`.

## `rgblast origin`

```bash
rgblast origin --store /data/stores/jungle --all --exclude-masked --query q.fa
```

Origin reads the queries (or stdin), runs the searches and prints the report, taking the
report flags (`--format`, `--detail`, `--tie-percent`, `--taxdump`, ...) and the engine flags
of `search`. With `--engine blastn --db BASE` it runs one exhaustive search of that database;
staged mode needs the default engine, rgblastn. Arguments after `--` go to every stage, except
the ones origin owns (`-outfmt`, `-out`, `-dbsize`, `-seqidlist`, `-negative_seqidlist` and the
`-taxids` family).

**Staged mode** (the default when it pays off):

1. **Representatives.** Stage 1 searches one assembly per species. Its deflines still list every
   genome that holds each sequence, so a hit on the human representative's chrM already names
   all 37 genomes with that identical chrM.
2. **Fan-out.** A query's candidate species are those whose best stage 1 genome scores within
   `--fanout-percent` (default 20) of its best. Stage 2 searches, for the queries with candidates,
   the candidates' sequences that stage 1 did not cover. A query's stage 2 hits are reported only
   for its own candidate species, so its result never depends on the other queries in the batch.
3. **Statistics.** Every stage searches under the full selection's sequence count and length
   (rgblastn's `STATS_NSEQ`/`STATS_TOTLEN` override, what a wrapper alias does for blastn), so
   E-values and `-evalue` cutoffs are those of a full search; the golden test checks that each
   stage's hits equal blastn's under a wrapper alias.

`--exhaustive` runs one search of the whole selection; `--staged` forces staged mode. Without
either, origin runs staged when the representatives plus the largest species remainder are at
most 0.8 of the selection, else exhaustively, saying so. `--keep-hits DIR` keeps each stage's
hits, the query files and `origin.meta.json` (mode, statistics, representatives, each query's
candidates and scope).

### Representatives (`rgblast reps`)

Within each species (taxid, or `--group-rank` ancestor), assemblies with an N fraction below
0.25 are candidates (all of them, when none qualifies). The representative is the first that
applies of:

1. `--rep COLLECTION` (repeatable);
2. the only candidate with an alias in `--rep-namespace` (default `reference`);
3. **most-shared**: the largest sum, over its sequences, of length x (number of other candidates
   holding the sequence), so the assembly sharing the most bases with the rest of its species;
   ties go to the larger assembly, then the smaller digest;
4. `largest` when nothing is shared, `only` for a species with one assembly.

```
$ rgblast reps --store jungle --all --exclude-masked
#group	taxon	collections	rep	rep_collection	rule	rep_gbp	group_gbp	not_in_rep_seqs	not_in_rep_gbp
9606	Homo sapiens	45	GRCh38.p14-fasta-genomic	u1HyLgIlq8M_XvEwy0oGqAvKGHJMGtxH	most-shared	3.30	12.75	3225	9.45
10090	Mus musculus	22	GRCm38.p6-fasta-genomic-GCA	JL56x8L1q1Fs_-jHvZxBG01Vitac-CmO	most-shared	2.82	8.33	219	5.51
# representatives: 2 of 67 collections, 6.12 of 21.08 Gbp (3.4x smaller)
# worst-case staged query (genome level): 15.57 Gbp = 0.74 of a full search -> staged mode on
```

### When staged and exhaustive results differ

1. **The best sequence is in no candidate species' representative.** With no stage 1 hit at all
   the query is `none` with `reps-only` (rerun it with `--exhaustive`). On jungle this is
   the EBV decoy (`SQ.uaSHw…`, in 13 human genomes, not in GRCh38.p14-fasta-genomic). Decoys, alt
   and patch contigs, HLA alleles and regions only another build has behave the same way. If the
   only stage 1 hits are in other species, the wrong species are fanned out. If a weaker copy is
   in the right representative, stage 2 still finds the exact genomes (the hs38d1 decoy below).
2. **The true species' representative scores below the fan-out cut** while another of its
   genomes would have scored best.
3. **`-max_target_seqs` truncation** sees different subject pools (flag `truncated`).
4. **A query whose stage 2 hit is shared with a species that is not its candidate** reports that
   hit only for its own candidates.

E-values are the same in both modes. The `searched` column and the `reps-only` flag show each
case per query.

### Cost, and the honest scope

Bases each query is compared with on jungle (stage 1 = representatives; stage 2 = the matched
species minus what stage 1 covered):

| Query species | All 79: staged / full | | Masked excluded (67): staged / full | |
|---|---|---|---|---|
| Homo sapiens | 25.31 / 41.81 Gbp | 1.65x | 15.57 / 21.08 Gbp | 1.35x |
| Mus musculus | 22.62 / 41.81 Gbp | 1.85x | 11.63 / 21.08 Gbp | 1.81x |
| Species answer only (stage 1) | 6.12 / 41.81 Gbp | 6.8x | 6.12 / 21.08 Gbp | 3.4x |

Timings on the jungle mirror (masked assemblies excluded, 16 cores, 11 queries from
`scripts/jungle_queries.sh`, megablast), with each stage searched by BLAST+ 2.17.0 blastn on a
prebuilt database:

| staged | exhaustive | queries whose summary differs |
|---|---|---|
| 41.1 s | 49.9 s | `q_ebv_100001` only (staged `none reps-only`; exhaustive 13 human genomes) |

Wall time follows the scanned bases only roughly: each stage is a separate search with its own
startup, and a query's cost depends on how repetitive its region is (the MHC alt-contig query
dominates the time above).

Staging saves the most for stores with many species and several
assemblies per species. At genome level on jungle it saves little, and on VGP-like stores (585
species for 605 assemblies) the representatives are nearly the whole store, so origin runs one
full search.

### When BLAST is the right tool

BLAST alignments suit queries of hundreds to thousands of bases, divergent queries (80–90%
identity) and answers that need an alignment people can inspect. For classifying billions of
short reads by exact k-mer matches,
[Kraken2](https://github.com/DerrickWood/kraken2) or [sourmash](https://sourmash.readthedocs.io/)
fit better. A k-mer prefilter would be about as large as the store at mammal scale (a dense
minimizer index of jungle is 16–50 GB against a 15.9 GB store), so `origin` uses representatives
and fan-out instead.

## From Rust

`gtars_rgblast::api::{search, report, origin, reps}` print nothing; progress, engine stderr
lines and a heartbeat (at least once a second while blastn runs) arrive through the progress
callback, and returning `false` from it stops the engine and returns `Error::Cancelled`. The
`rgblast` argument structs (`gtars_rgblast::cli::{SearchArgs, ReportArgs, OriginArgs,
RepsArgs}`) convert to
these requests with `TryFrom`, which re-checks every command-line rule.
