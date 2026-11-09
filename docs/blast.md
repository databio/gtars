# gtars-rgblast: BLAST over a RefgetStore

`gtars-rgblast` turns a [RefgetStore](refget.md) into [BLAST](https://blast.ncbi.nlm.nih.gov/) databases, searches a store, and reports where a sequence comes from. A store keeps every unique sequence once, so building a database from it reads each sequence once, however many assemblies share it. Building hg38 takes about 1.3 s on 16 cores, against 64 to 87 s for `makeblastdb` on a FASTA file.

What you can do:

- **Build** a BLAST version 5 nucleotide database from a store (`build`). `blastn` and `blastdbcmd` read it like any `makeblastdb` database. Check it against the store with `verify`, write subset views with `view`, and count N bases with `nstats`.
- **Search** a store directly, with no database and no BLAST+ install (`search`), or run BLAST+ `blastn` on a database you built.
- **Find where a sequence comes from** (`origin`): which species and which genomes could it have come from. `report` turns a hit file from `search` into the same answer, and `reps` shows the representative assembly of each species.

## One command set, five ways to call it

Every verb exists with the same name and the same options on every surface.

| Verb | Binary | `gtars` CLI | Python | R | `refget` CLI |
|---|---|---|---|---|---|
| build | `rgmakeblastdb build` | `gtars blast build` | `gtars.blast.build()` | `blast_build()` | `refget blast build` |
| verify | `rgmakeblastdb verify` | `gtars blast verify` | `gtars.blast.verify()` | `blast_verify()` | `refget blast verify` |
| view | `rgmakeblastdb view` | `gtars blast view` | `gtars.blast.view()` | `blast_view()` | `refget blast view` |
| nstats | `rgmakeblastdb nstats` | `gtars blast nstats` | `gtars.blast.nstats()` | `blast_nstats()` | `refget blast nstats` |
| search | `rgblast search` | `gtars blast search` | `gtars.blast.search()` | `blast_search()` | `refget blast search` |
| report | `rgblast report` | `gtars blast report` | `gtars.blast.report()` | `blast_report()` | `refget blast report` |
| origin | `rgblast origin` | `gtars blast origin` | `gtars.blast.origin()` | `blast_origin()` | `refget blast origin` |
| reps | `rgblast reps` | `gtars blast reps` | `gtars.blast.reps()` | `blast_reps()` | `refget blast reps` |

The [`refget` command line](https://docs.refgenie.org/refget/reference/cli/#blast-commands) also takes the store from `-p/--path` or `-r/--remote`.

The Python and R names are the long flag names with `-` written as `_` (`--max-n-fraction` is `max_n_fraction`). One table of options therefore covers all of them, and `gtars blast <verb> --help` lists every flag of a verb:

| Group | Options | Used by |
|---|---|---|
| Store | `--store PATH\|URL` (a local folder or an `http(s)` prefix) | build, verify, nstats, reps; search and origin with rgblastn |
| Selection | `--all`, `--collection`, `--collections-file`, `--namespace`, `--taxon`, `--exclude-masked`, `--max-n-fraction`, `--sequences-file` | build, view, nstats, reps; search and origin with rgblastn |
| Fetch (remote stores) | `--cache-dir`, `--resume`, `--fetch-jobs` | every verb that takes `--store` |
| Build | `--out`, `--title`, `--taxid-map`, `--max-file-size`, `--threads`, `--budget-mb`, `--on-oversize`, `--dry-run`, `--force`, `--no-space-check` | build |
| Database | `--db BASE` | verify, view, report; search and origin with blastn |
| Views | `--per-collection`, `--name`, `--out-dir`, `--title`, `--full-db-stats` | view |
| Engine | `--engine rgblastn\|blastn`, `--blastn PATH` | search, origin |
| Search | `--no-extra-fields`, `--meta FILE` | search |
| Report | `--format tsv\|json`, `--output`, `--detail`, `--detail-top`, `--taxdump`, `--group-rank`, `--tie-percent`, `--min-pident`, `--min-qcov`, `--max-evalue` | report, origin |
| Report input | `--meta`, `--query`, `--fields`, `--max-target-seqs` | report |
| Origin | `--query`, `--exhaustive`, `--staged`, `--fanout-percent`, `--rep`, `--rep-namespace`, `--keep-hits` | origin (`--rep` and `--rep-namespace` also on reps) |

Collections are chosen by digest, `NAMESPACE:ALIAS` (for example `ucsc:hg38`) or a bare alias. The engine arguments of `search` and `origin` (blastn's own flags, such as `-query` and `-evalue`) go after `--` on the command line, and in `args=` in Python and R.

## Quick start

```bash
# build a database from one assembly, check it, and search it with blastn
gtars blast build --store ./store --collection ucsc:hg38 --out db/hg38
gtars blast verify --db db/hg38 --store ./store
blastn -db db/hg38 -query q.fa -outfmt 6

# or search the store directly, with no database
gtars blast search --store ./store --all -- -query q.fa -outfmt 6

# which genomes could q.fa come from?
gtars blast origin --store ./store --all --query q.fa
```

## Engines

`rgblastn` (the default) is a Rust port of `blastn` inside gtars. It reads the store and selection directly, needs no database and no BLAST+, and covers the common searches: tasks `megablast`, `blastn` and `blastn-short`, with `-outfmt 6`. It refuses any other blastn option by name. The standalone `rgblastn` binary is the same engine.

`--engine blastn --db BASE` runs [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) `blastn` on a database you built with `build` or `view`, for every other blastn option. Choose collections when you build the database; `search` then takes no `--store` or selection flags. blastn comes from `--blastn PATH`, else `$GTARS_RGBLAST_BLASTN`, else the `PATH`. If none is found, the command exits 2 and says how to install BLAST+. `origin` with blastn runs one search of the whole database; its staged search (representatives first) runs with rgblastn.

## The ids in the hits

Each unique sequence is one BLAST subject, named by its sequence digest: `SQ.<digest>`. The hit's `sallseqid` column lists, for every assembly that holds the sequence, a second id `gnl|<collection alias>|<sequence name>`, and `salltitles` carries the species and assembly count. `report` and `origin` read these columns to name the genomes. For `-outfmt 6` and `7`, `search` appends the columns they need (`qlen score sallseqid salltitles staxids`) unless you pass `--no-extra-fields`.

## Remote stores

`--store` also takes an `http(s)` URL. A command that needs the sequences first pulls the selected ones into a local cache (the same steps as `refget store pull`), then reads them from there. Pulled sequences stay in the cache, so running the command again downloads only what is missing. The cache root is `--cache-dir`, else `$REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`, else `~/.cache/refget`. `--fetch-jobs N` sets the number of parallel downloads (default 32).

## Exit codes

0 success, 1 runtime error, 2 usage error or no BLAST engine, 130 cancelled with Ctrl-C (a database being built is removed). For `search`, the engine's own exit status passes through. Closing the output pipe (`| head`) exits 0. `refget blast` uses the same codes.

## More

- Python: [BLAST tutorial](python/blast.md) and [API reference](python/blast-api.md)
- R: [BLAST tutorial](r/blast.md)
- Command line: [gtars blast](cli.md#blast)
- `refget`: [Search a RefgetStore with BLAST](https://docs.refgenie.org/refget/using-services/blast-store/)
- Crate: [README](https://github.com/databio/gtars/tree/master/gtars-rgblast), [`origin` details](https://github.com/databio/gtars/blob/master/gtars-rgblast/docs/origin.md), and [docs.rs](https://docs.rs/gtars-rgblast/latest/gtars_rgblast/)
