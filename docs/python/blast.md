# BLAST over a RefgetStore in Python

`gtars.blast` builds [BLAST](https://blast.ncbi.nlm.nih.gov/) databases from a [RefgetStore](refget-store.md), searches a store, and tells you which genomes a sequence could come from. The functions have the same names and options as [`gtars blast`](../blast.md) on the command line: each keyword is a long flag with `-` written as `_`, and `gtars blast <verb> --help` lists them all. Results are plain dicts, so they work with `pandas`, `json` and everything else. See the [API reference](blast-api.md) for the functions.

## Make a store

Any store works, local or an `http(s)` URL. This one has two small genomes that share a chromosome:

```python
import random
from pathlib import Path
from gtars.refget import RefgetStore

rng = random.Random(7)
chroms = ["".join(rng.choice("ACGT") for _ in range(2000)) for _ in range(3)]

work = Path("blast_demo"); work.mkdir(exist_ok=True)
store = RefgetStore.on_disk(str(work / "store"))
for name, recs in {"genome_a": chroms[:2], "genome_b": chroms[1:]}.items():
    fa = work / f"{name}.fa"
    fa.write_text("".join(f">{name}_{i}\n{s}\n" for i, s in enumerate(recs)))
    store.add_sequence_collection_from_fasta(str(fa))
del store                       # the on-disk store is complete once the collections are added

(work / "query.fa").write_text(">q1\n" + chroms[0][500:900] + "\n")
```

## Build a database

```python
from gtars import blast

summary = blast.build(work / "store", work / "db", all=True)
print(summary["preflight"]["subjects"], "unique sequences")
print(blast.verify(work / "db", work / "store")["ok"])
```

`store` and `out` can be strings or paths. Choose collections with `all=True`, or `collection="ucsc:hg38"`, `namespace=`, `taxon=` and the other selection options (a value or a list). `dry_run=True` returns what a build would do without writing anything.

### Progress and Ctrl-C

`progress=True` prints a line per step to stderr. Pass a function to get each event as a dict, for example to drive a [tqdm](https://tqdm.github.io/) bar. Every event has `event` and `message`; `sequences` events carry `done` and `total`:

```python
from tqdm import tqdm

bar = None
def show(ev):
    global bar
    if ev["event"] == "sequences":
        if bar is None:
            bar = tqdm(total=ev["total"], unit="seq")
        bar.n = ev["done"]; bar.refresh()

blast.build(work / "store", work / "db2", all=True, progress=show)
```

The callback runs on the calling thread. Ctrl-C cancels the build at the next event, raises `KeyboardInterrupt`, and removes partial output. An exception raised inside the callback cancels the build too, and propagates unchanged.

## Search

`search` runs rgblastn, a Rust port of blastn, in process on your selection of the store. It needs no database and no BLAST+ install, and covers tasks `megablast`, `blastn` and `blastn-short` with `-outfmt 6`.

```python
r = blast.search(work / "store", all=True, query=work / "query.fa", outfmt="6")
print(r["output"])
```

With `out=None` the hits come back in `r["output"]`. Give `out="hits.tsv"` to write a file (the run metadata goes next to it in `hits.tsv.meta.json`). Extra blastn options go in `args=["-evalue", "1e-20"]`.

For any other blastn option, run BLAST+ on a database you built: `blast.search(query=..., engine="blastn", db=work / "db")`. blastn comes from `blastn=`, else `$GTARS_RGBLAST_BLASTN`, else the `PATH`.

### Errors

```python
try:
    blast.search(query=work / "query.fa", engine="blastn", db=work / "db",
                 blastn="/nonexistent/blastn")
except blast.EngineNotFoundError as e:      # also a FileNotFoundError
    print("no BLAST+:", e)                  # the message says how to install it
```

A bad or conflicting option raises `ValueError` (with the command line's message), and a missing store, database or collection raises `FileNotFoundError`. A blastn run that fails raises `blast.EngineError` with `returncode`. Every exception here has a `kind` (`"not_found"`, `"invalid_input"`, `"conflict"`, `"io"`, `"network"`, `"data"`) that matches the R condition classes.

## Where did this sequence come from?

`origin` searches one representative assembly per species, then the rest of each species that matched, and returns the best hit per genome:

```python
r = blast.origin(work / "store", all=True, query=work / "query.fa")
for row in r["summary"]:
    print(row["query"], row["status"], row["taxon"], row["top_genome"], row["pident"])
```

`r["summary"]` is a list of dicts, one per query, with the columns of the TSV report. `r["report"]` is the full JSON document and `r["meta"]` describes the run. To start from hits you already have, use `blast.report("hits.tsv")`, which reads `hits.tsv.meta.json` from `search`, or `blast.report("hits.tsv", db="db/hg38")` for hits from a database. Load the summary into pandas with `pandas.DataFrame(r["summary"])`.

`blast.reps(store, all=True)` shows the representative of each species and the predicted cost of a staged search.

## Views and N counts

`blast.view(db, per_collection=True)` writes one searchable view per assembly with exact statistics, and `blast.nstats(store, all=True)` counts N bases per collection, to choose `max_n_fraction`.
