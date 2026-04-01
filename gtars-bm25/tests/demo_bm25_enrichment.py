"""
BM25 vs LOLA Enrichment Demo

Compares enrichment analysis using:
  1. BM25 sparse embeddings (new approach)
  2. Traditional LOLA Fisher's exact test (gtars-lola)

Uses the lola_hg38_ucsc_features bedset from BEDbase as the reference database.

Usage:
    python demo_bm25_enrichment.py <query.bed>
    python demo_bm25_enrichment.py  # uses a sample BED from BEDbase
"""

import sys
import gzip
import json
import time
import urllib.request
from pathlib import Path
from collections import defaultdict

from gtars.bm25 import Bm25
from gtars.lola import RegionDB, run_lola

BEDSET_ID = "lola_hg38_ucsc_features"
API_BASE = "https://api.bedbase.org/v1"
CACHE_DIR = Path(__file__).parent / "demo_cache"
RESULTS_PATH = Path(__file__).parent / "demo_cache" / "results.json"


# ── Helpers ──────────────────────────────────────────────────────────────

def _request(url: str) -> urllib.request.Request:
    return urllib.request.Request(url, headers={"User-Agent": "gtars-bm25-demo/0.1"})


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(_request(url)) as r:
        return json.loads(r.read())


def download_file(url: str, dest: Path):
    if dest.exists():
        return
    print(f"    {dest.name}")
    with urllib.request.urlopen(_request(url)) as r, open(dest, "wb") as f:
        f.write(r.read())


def fetch_bedset_files() -> list[dict]:
    """Fetch metadata for all BED files in the LOLA UCSC bedset."""
    data = fetch_json(f"{API_BASE}/bedset/{BEDSET_ID}/bedfiles")
    return data["results"]


def download_bed(bed_id: str, name: str, cache: Path) -> Path:
    """Download a BED file from BEDbase, return path to decompressed file."""
    gz_path = cache / f"{name}.bed.gz"
    bed_path = cache / f"{name}.bed"

    if bed_path.exists():
        return bed_path

    files_meta = fetch_json(f"{API_BASE}/bed/{bed_id}/metadata/files")
    url = files_meta["bed_file"]["access_methods"][0]["access_url"]["url"]
    download_file(url, gz_path)

    with gzip.open(gz_path, "rb") as f_in, open(bed_path, "wb") as f_out:
        f_out.write(f_in.read())

    return bed_path


def read_bed_tuples(path: Path) -> list[tuple[str, int, int]]:
    """Read a BED file into a list of (chr, start, end) tuples."""
    regions = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def get_sample_query(cache: Path) -> str:
    """Download a sample BED file to use as query."""
    sample_path = cache / "sample_query.bed"
    if sample_path.exists():
        return str(sample_path)

    print("  No query provided — fetching a sample from BEDbase...")
    beds = fetch_json(f"{API_BASE}/bed/list?genome=hg38&limit=1&offset=50")
    bed = beds["results"][0]
    bed_id = bed["id"]
    name = bed.get("name", bed_id)
    print(f"  Sample: {name} ({bed_id})")

    files_meta = fetch_json(f"{API_BASE}/bed/{bed_id}/metadata/files")
    url = files_meta["bed_file"]["access_methods"][0]["access_url"]["url"]
    download_file(url, cache / "sample_query.bed.gz")

    with gzip.open(cache / "sample_query.bed.gz", "rb") as f_in:
        with open(sample_path, "wb") as f_out:
            f_out.write(f_in.read())

    return str(sample_path)


# ── Data setup ───────────────────────────────────────────────────────────

def download_reference_beds(bed_files: list[dict], cache: Path) -> dict[str, Path]:
    """Download all reference BED files. Returns {name: path}."""
    ref_beds = {}
    for meta in bed_files:
        path = download_bed(meta["id"], meta["name"], cache)
        ref_beds[meta["name"]] = path
    return ref_beds


def build_vocab_and_mapping(
    bed_files: list[dict], ref_beds: dict[str, Path], cache: Path
) -> tuple[Path, dict[int, str], dict[str, str]]:
    """Concatenate reference BEDs into one vocab, tracking annotation per line."""
    vocab_path = cache / "vocab.bed"
    token_to_anno: dict[int, str] = {}
    anno_descriptions: dict[str, str] = {}
    anno_region_counts: dict[str, int] = {}
    line_num = 0

    with open(vocab_path, "w") as vocab_f:
        for meta in bed_files:
            name = meta["name"]
            anno_descriptions[name] = meta["description"]
            count = 0

            for chr, start, end in read_bed_tuples(ref_beds[name]):
                vocab_f.write(f"{chr}\t{start}\t{end}\n")
                token_to_anno[line_num] = name
                line_num += 1
                count += 1

            anno_region_counts[name] = count

    return vocab_path, token_to_anno, anno_descriptions


# ── IDF computation ──────────────────────────────────────────────────────

def compute_idf(
    bm25: Bm25,
    ref_beds: dict[str, Path],
) -> dict[int, float]:
    """Compute IDF for each token across the reference corpus.

    IDF(token) = ln(1 + (N - df + 0.5) / (df + 0.5))
    where N = number of documents, df = documents containing token.
    """
    import math

    # Embed each reference file, collect which tokens appear in each
    doc_freq: dict[int, int] = defaultdict(int)
    n_docs = len(ref_beds)

    for name, path in ref_beds.items():
        sv = bm25.embed(str(path))
        seen_tokens = set(sv.indices)
        for token_id in seen_tokens:
            doc_freq[token_id] += 1

    # Compute IDF
    idf = {}
    for token_id, df in doc_freq.items():
        idf[token_id] = math.log(1 + (n_docs - df + 0.5) / (df + 0.5))

    return idf


# ── BM25 enrichment ─────────────────────────────────────────────────────

def run_bm25_enrichment(
    query_path: str,
    bm25: Bm25,
    token_to_anno: dict[int, str],
    idf: dict[int, float],
) -> tuple[dict[str, dict], float]:
    """Run BM25 enrichment with IDF. Returns ({name: {score, hits}}, elapsed_seconds)."""
    t0 = time.perf_counter()
    sv = bm25.embed(query_path)
    elapsed = time.perf_counter() - t0

    if len(sv) == 0:
        return {}, elapsed

    anno_scores: dict[str, float] = defaultdict(float)
    anno_hits: dict[str, int] = defaultdict(int)

    for idx, val in zip(sv.indices, sv.values):
        anno = token_to_anno.get(idx)
        if anno is None:
            for offset in range(1, 10):
                anno = token_to_anno.get(idx - offset)
                if anno:
                    break
        if anno is None:
            continue
        # Apply IDF: full BM25 = IDF * TF_score
        token_idf = idf.get(idx, 0.0)
        anno_scores[anno] += val * token_idf
        anno_hits[anno] += 1

    results = {}
    for name in anno_scores:
        results[name] = {
            "score": round(anno_scores[name], 4),
            "hits": anno_hits[name],
        }

    return results, elapsed


# ── LOLA enrichment ──────────────────────────────────────────────────────

def run_lola_enrichment(
    query_path: str,
    ref_beds: dict[str, Path],
) -> tuple[dict[str, dict], float]:
    """Run traditional LOLA. Returns ({name: {support, fraction, ...}}, elapsed)."""
    bed_paths = list(ref_beds.values())
    bed_names = list(ref_beds.keys())

    query_tuples = read_bed_tuples(Path(query_path))
    n_query = len(query_tuples)

    # Universe = union of query + all reference regions
    universe_tuples = list(query_tuples)
    for path in bed_paths:
        universe_tuples.extend(read_bed_tuples(path))

    t0 = time.perf_counter()
    region_db = RegionDB.from_bed_files(
        [str(p) for p in bed_paths],
        filenames=bed_names,
    )

    raw = run_lola(
        user_sets=[query_tuples],
        universe=universe_tuples,
        region_db=region_db,
        min_overlap=1,
        direction="enrichment",
    )
    elapsed = time.perf_counter() - t0

    results = {}
    n = len(raw["filename"])
    for i in range(n):
        name = raw["filename"][i]
        support = raw["support"][i]
        results[name] = {
            "support": support,
            "fraction": round(support / n_query, 4) if n_query > 0 else 0,
        }

    return results, elapsed


# ── Display ──────────────────────────────────────────────────────────────

def print_results(
    bm25_results: dict[str, dict],
    bm25_time: float,
    lola_results: dict[str, dict],
    lola_time: float,
    anno_descriptions: dict[str, str],
):
    """Pretty-print side-by-side comparison."""

    # Rank by BM25 score and by LOLA support
    bm25_ranked = sorted(bm25_results.keys(), key=lambda n: bm25_results[n]["score"], reverse=True)
    sup_ranked = sorted(lola_results.keys(), key=lambda n: lola_results[n]["support"], reverse=True)

    bm25_rank = {name: i + 1 for i, name in enumerate(bm25_ranked)}
    sup_rank = {name: i + 1 for i, name in enumerate(sup_ranked)}

    all_annos = list(dict.fromkeys(bm25_ranked + sup_ranked))

    # Header
    print()
    print("┌────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print("│                          BM25 vs Overlap Enrichment Results                                   │")
    print("├──────────────────────┬────────┬──────────┬────────┬──────────┬──────────┬─────────────────────┤")
    print("│ Annotation           │ BM25 # │ BM25 Scr │ Sup #  │ Support  │ Fraction │ Description         │")
    print("├──────────────────────┼────────┼──────────┼────────┼──────────┼──────────┼─────────────────────┤")

    for anno in sorted(all_annos, key=lambda a: bm25_rank.get(a, 99)):
        br = str(bm25_rank.get(anno, "-"))
        bs = bm25_results.get(anno, {}).get("score", 0)
        sr = str(sup_rank.get(anno, "-"))
        sup = lola_results.get(anno, {}).get("support", 0)
        frac = lola_results.get(anno, {}).get("fraction", 0)
        desc = anno_descriptions.get(anno, "")[:19]

        print(
            f"│ {anno:<20} │ {br:>6} │ {bs:>8.2f} │ {sr:>6} │ {sup:>8} │ {frac:>8.1%} │ {desc:<19} │"
        )

    print("└──────────────────────┴────────┴──────────┴────────┴──────────┴──────────┴─────────────────────┘")

    # Timing
    print(f"\n  BM25 time:  {bm25_time:.3f}s")
    print(f"  LOLA time:  {lola_time:.3f}s")

    # Rank correlation (BM25 vs support rank)
    shared = [a for a in all_annos if a in bm25_rank and a in sup_rank]
    if len(shared) >= 2:
        bm25_r = [bm25_rank[a] for a in shared]
        sup_r = [sup_rank[a] for a in shared]
        n = len(shared)
        d_sq = sum((b - s) ** 2 for b, s in zip(bm25_r, sup_r))
        rho = 1 - (6 * d_sq) / (n * (n ** 2 - 1))
        print(f"  Rank corr:  {rho:.3f} (Spearman, BM25 vs support, n={n})")


def save_json(
    query_path: str,
    bm25_results: dict[str, dict],
    bm25_time: float,
    lola_results: dict[str, dict],
    lola_time: float,
    anno_descriptions: dict[str, str],
    output_path: Path,
):
    """Save comprehensive JSON output."""
    # Add ranks
    bm25_ranked = sorted(bm25_results.keys(), key=lambda n: bm25_results[n]["score"], reverse=True)
    lola_ranked = sorted(lola_results.keys(), key=lambda n: lola_results[n]["support"], reverse=True)

    combined = []
    all_annos = list(dict.fromkeys(bm25_ranked + lola_ranked))

    for anno in all_annos:
        entry = {
            "annotation": anno,
            "description": anno_descriptions.get(anno, ""),
            "bm25": {
                "rank": bm25_ranked.index(anno) + 1 if anno in bm25_ranked else None,
                **(bm25_results.get(anno, {})),
            },
            "lola": {
                "rank": lola_ranked.index(anno) + 1 if anno in lola_ranked else None,
                **(lola_results.get(anno, {})),
            },
        }
        combined.append(entry)

    output = {
        "query": str(query_path),
        "reference_bedset": BEDSET_ID,
        "num_annotations": len(all_annos),
        "timing": {
            "bm25_seconds": round(bm25_time, 4),
            "lola_seconds": round(lola_time, 4),
        },
        "results": combined,
    }

    output_path.parent.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n  JSON saved: {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    cache = CACHE_DIR
    cache.mkdir(exist_ok=True)

    print()
    print("  BM25 vs LOLA Enrichment Comparison")
    print(f"  Reference: {BEDSET_ID}")
    print()

    # Fetch and download
    print("  Fetching reference files...")
    bed_files = fetch_bedset_files()
    print(f"  {len(bed_files)} annotations found\n")

    print("  Downloading BED files...")
    ref_beds = download_reference_beds(bed_files, cache)

    # Build BM25 vocab
    print("\n  Building BM25 vocabulary...")
    vocab_path, token_to_anno, anno_descriptions = build_vocab_and_mapping(
        bed_files, ref_beds, cache
    )
    bm25 = Bm25(tokenizer=str(vocab_path), k=1.5, b=0.75, avg_doc_length=1000.0)
    print(f"  Vocab: {bm25.vocab_size} tokens")

    # Compute IDF across reference corpus
    print("\n  Computing IDF across reference files...")
    idf = compute_idf(bm25, ref_beds)
    print(f"  {len(idf)} tokens with IDF scores")

    # Query
    if len(sys.argv) > 1:
        query_path = sys.argv[1]
    else:
        query_path = get_sample_query(cache)

    query_regions = read_bed_tuples(Path(query_path))
    print(f"\n  Query: {Path(query_path).name} ({len(query_regions)} regions)")

    # Run both
    print("\n  Running BM25 (TF * IDF)...")
    bm25_results, bm25_time = run_bm25_enrichment(query_path, bm25, token_to_anno, idf)
    print(f"  done ({bm25_time:.3f}s)")

    print("  Running LOLA...")
    lola_results, lola_time = run_lola_enrichment(query_path, ref_beds)
    print(f"  done ({lola_time:.3f}s)")

    # Output
    print_results(bm25_results, bm25_time, lola_results, lola_time, anno_descriptions)
    save_json(query_path, bm25_results, bm25_time, lola_results, lola_time, anno_descriptions, RESULTS_PATH)


if __name__ == "__main__":
    main()
