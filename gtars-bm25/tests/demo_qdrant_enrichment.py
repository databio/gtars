"""
BM25 Enrichment via Qdrant Sparse Vector Search

Indexes a LOLA reference database into Qdrant (in-memory) as sparse vectors,
then searches with a query BED file. Qdrant applies IDF at search time.

This demonstrates the production-ready approach: embed once, search many times,
IDF computed automatically from the corpus.

Usage:
    python demo_qdrant_enrichment.py [bedset_id] [query.bed]
    python demo_qdrant_enrichment.py                                    # UCSC features, sample query
    python demo_qdrant_enrichment.py lola_hg38_encode_tfbs              # ENCODE TFBS, sample query
    python demo_qdrant_enrichment.py lola_hg38_encode_tfbs query.bed    # ENCODE TFBS, custom query
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
from qdrant_client import QdrantClient, models

API_BASE = "https://api.bedbase.org/v1"
CACHE_DIR = Path(__file__).parent / "demo_cache"


# ── Helpers ──────────────────────────────────────────────────────────────

def _request(url: str) -> urllib.request.Request:
    return urllib.request.Request(url, headers={"User-Agent": "gtars-bm25-demo/0.1"})


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(_request(url)) as r:
        return json.loads(r.read())


def download_file(url: str, dest: Path):
    if dest.exists():
        return
    with urllib.request.urlopen(_request(url)) as r, open(dest, "wb") as f:
        f.write(r.read())


def fetch_bedset_files(bedset_id: str) -> list[dict]:
    """Fetch all BED file metadata from a bedset, handling pagination."""
    results = []
    data = fetch_json(f"{API_BASE}/bedset/{bedset_id}/bedfiles")
    results.extend(data["results"])
    return results


def download_bed(bed_id: str, name: str, cache: Path) -> Path:
    """Download and decompress a BED file from BEDbase."""
    bed_path = cache / f"{name}.bed"
    if bed_path.exists():
        return bed_path

    gz_path = cache / f"{name}.bed.gz"
    files_meta = fetch_json(f"{API_BASE}/bed/{bed_id}/metadata/files")
    url = files_meta["bed_file"]["access_methods"][0]["access_url"]["url"]
    download_file(url, gz_path)

    with gzip.open(gz_path, "rb") as f_in, open(bed_path, "wb") as f_out:
        f_out.write(f_in.read())

    return bed_path


def read_bed_tuples(path: Path) -> list[tuple[str, int, int]]:
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


# ── Qdrant setup ─────────────────────────────────────────────────────────

def create_collection(client: QdrantClient, name: str):
    """Create a Qdrant collection with sparse vectors + IDF modifier."""
    client.create_collection(
        collection_name=name,
        vectors_config={},
        sparse_vectors_config={
            "bm25": models.SparseVectorParams(
                modifier=models.Modifier.IDF,
            )
        },
    )


def index_reference_files(
    client: QdrantClient,
    collection: str,
    bm25: Bm25,
    bed_files: list[dict],
    cache: Path,
) -> dict[str, str]:
    """Download, embed, and index all reference BED files into Qdrant.

    Returns {name: description} mapping.
    """
    anno_descriptions = {}
    points = []
    skipped = 0

    for i, meta in enumerate(bed_files):
        name = meta["name"]
        desc = meta.get("description", "")
        bed_id = meta["id"]
        anno_descriptions[name] = desc

        try:
            bed_path = download_bed(bed_id, name, cache)
        except Exception as e:
            print(f"    SKIP {name}: {e}")
            skipped += 1
            continue

        sv = bm25.embed(str(bed_path))

        if len(sv) == 0:
            skipped += 1
            continue

        points.append(
            models.PointStruct(
                id=i,
                vector={
                    "bm25": models.SparseVector(
                        indices=list(sv.indices),
                        values=list(sv.values),
                    )
                },
                payload={
                    "name": name,
                    "description": desc,
                    "bed_id": bed_id,
                    "annotation": meta.get("annotation", {}),
                },
            )
        )

        if len(points) % 50 == 0:
            # Batch upsert
            client.upsert(collection_name=collection, points=points)
            print(f"    indexed {i + 1}/{len(bed_files)}...")
            points = []

    # Final batch
    if points:
        client.upsert(collection_name=collection, points=points)

    print(f"    indexed {len(bed_files) - skipped}/{len(bed_files)} files ({skipped} skipped)")
    return anno_descriptions


def search_enrichment(
    client: QdrantClient,
    collection: str,
    bm25: Bm25,
    query_path: str,
    limit: int = 20,
) -> tuple[list[dict], float]:
    """Embed query and search Qdrant. Returns (results, elapsed)."""
    t0 = time.perf_counter()

    sv = bm25.embed(query_path)
    if len(sv) == 0:
        return [], time.perf_counter() - t0

    results = client.query_points(
        collection_name=collection,
        query=models.SparseVector(
            indices=list(sv.indices),
            values=list(sv.values),
        ),
        using="bm25",
        limit=limit,
        with_payload=True,
    )

    elapsed = time.perf_counter() - t0

    ranked = []
    for point in results.points:
        ranked.append({
            "name": point.payload["name"],
            "description": point.payload["description"],
            "bed_id": point.payload["bed_id"],
            "score": round(point.score, 4),
            "annotation": point.payload.get("annotation", {}),
        })

    return ranked, elapsed


# ── LOLA support ─────────────────────────────────────────────────────────

def compute_lola_support(
    query_path: str,
    ref_beds: dict[str, Path],
) -> tuple[dict[str, dict], float]:
    """Compute LOLA support/fraction for comparison. Returns ({name: {support, fraction}}, elapsed)."""
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
    for i in range(len(raw["filename"])):
        name = raw["filename"][i]
        support = raw["support"][i]
        results[name] = {
            "support": support,
            "fraction": round(support / n_query, 4) if n_query > 0 else 0,
        }

    return results, elapsed


# ── Display ──────────────────────────────────────────────────────────────

def print_results(
    qdrant_results: list[dict],
    search_time: float,
    lola_results: dict[str, dict],
    lola_time: float,
    bedset_id: str,
):
    # Build LOLA support rank lookup
    sup_ranked = sorted(lola_results.keys(), key=lambda n: lola_results[n]["support"], reverse=True)
    sup_rank = {name: i + 1 for i, name in enumerate(sup_ranked)}

    print()
    print(f"┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print(f"│  Qdrant BM25+IDF vs LOLA Support          Bedset: {bedset_id:<40}          │")
    print(f"├──────┬──────────────────────────────┬──────────┬────────┬──────────┬──────────┬────────────────────────────┤")
    print(f"│ BM25 │ Name                         │ BM25 Scr │ Sup #  │ Support  │ Fraction │ Description                │")
    print(f"├──────┼──────────────────────────────┼──────────┼────────┼──────────┼──────────┼────────────────────────────┤")

    for i, r in enumerate(qdrant_results, 1):
        name = r["name"]
        score = r["score"]
        desc = r["description"][:26]
        sr = str(sup_rank.get(name, "-"))
        sup = lola_results.get(name, {}).get("support", 0)
        frac = lola_results.get(name, {}).get("fraction", 0)

        print(
            f"│ {i:>4} │ {name[:28]:<28} │ {score:>8.2f} │ {sr:>6} │ {sup:>8} │ {frac:>8.1%} │ {desc:<26} │"
        )

    print(f"└──────┴──────────────────────────────┴──────────┴────────┴──────────┴──────────┴────────────────────────────┘")

    print(f"\n  BM25 search: {search_time:.3f}s (embed + Qdrant)")
    print(f"  LOLA time:   {lola_time:.3f}s")

    # Rank correlation
    shared = [r["name"] for r in qdrant_results if r["name"] in sup_rank]
    if len(shared) >= 2:
        bm25_r = [i + 1 for i, r in enumerate(qdrant_results) if r["name"] in sup_rank]
        sup_r = [sup_rank[r["name"]] for r in qdrant_results if r["name"] in sup_rank]
        n = len(shared)
        d_sq = sum((b - s) ** 2 for b, s in zip(bm25_r, sup_r))
        rho = 1 - (6 * d_sq) / (n * (n ** 2 - 1))
        print(f"  Rank corr:   {rho:.3f} (Spearman, BM25 vs support, n={n})")


def save_json(
    query_path: str,
    bedset_id: str,
    qdrant_results: list[dict],
    lola_results: dict[str, dict],
    search_time: float,
    lola_time: float,
    index_time: float,
    n_indexed: int,
    output_path: Path,
):
    # Merge BM25 and LOLA results
    combined = []
    for i, r in enumerate(qdrant_results, 1):
        name = r["name"]
        lola = lola_results.get(name, {})
        combined.append({
            "bm25_rank": i,
            "name": name,
            "description": r["description"],
            "bm25_score": r["score"],
            "lola_support": lola.get("support"),
            "lola_fraction": lola.get("fraction"),
        })

    output = {
        "query": str(query_path),
        "reference_bedset": bedset_id,
        "num_indexed": n_indexed,
        "timing": {
            "index_seconds": round(index_time, 4),
            "bm25_search_seconds": round(search_time, 4),
            "lola_seconds": round(lola_time, 4),
        },
        "results": combined,
    }

    output_path.parent.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"  JSON saved: {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    # Parse args
    bedset_id = sys.argv[1] if len(sys.argv) > 1 else "lola_hg38_ucsc_features"
    query_arg = sys.argv[2] if len(sys.argv) > 2 else None

    cache = CACHE_DIR / bedset_id
    cache.mkdir(parents=True, exist_ok=True)
    collection = f"lola_{bedset_id}"
    results_path = CACHE_DIR / f"qdrant_results_{bedset_id}.json"

    print()
    print("  Qdrant Sparse Vector Enrichment Demo")
    print(f"  Bedset: {bedset_id}")
    print()

    # Fetch bedset metadata
    print("  Fetching bedset metadata...")
    bed_files = fetch_bedset_files(bedset_id)
    print(f"  {len(bed_files)} files in bedset")

    # Build vocab from all reference files
    print("\n  Downloading reference files & building vocabulary...")
    t0 = time.perf_counter()

    # Download all first
    ref_beds = {}
    for meta in bed_files:
        try:
            path = download_bed(meta["id"], meta["name"], cache)
            ref_beds[meta["name"]] = path
        except Exception as e:
            print(f"    SKIP {meta['name']}: {e}")

    # Build concatenated vocab
    vocab_path = cache / "vocab.bed"
    with open(vocab_path, "w") as vocab_f:
        for name, path in ref_beds.items():
            for chr, start, end in read_bed_tuples(path):
                vocab_f.write(f"{chr}\t{start}\t{end}\n")

    download_time = time.perf_counter() - t0
    print(f"  Downloaded {len(ref_beds)} files in {download_time:.1f}s")

    # Build BM25 model
    print("\n  Building BM25 model...")
    bm25 = Bm25(tokenizer=str(vocab_path), k=1.5, b=0.75, avg_doc_length=1000.0)
    print(f"  Vocab: {bm25.vocab_size} tokens")

    # Create Qdrant collection and index
    print("\n  Indexing into Qdrant (in-memory)...")
    client = QdrantClient(":memory:")
    create_collection(client, collection)

    t0 = time.perf_counter()
    anno_descriptions = index_reference_files(client, collection, bm25, bed_files, cache)
    index_time = time.perf_counter() - t0
    print(f"  Index time: {index_time:.1f}s")

    # Get collection info
    info = client.get_collection(collection)
    print(f"  Points in collection: {info.points_count}")

    # Query
    if query_arg:
        query_path = query_arg
    else:
        query_path = get_sample_query(CACHE_DIR)

    query_regions = read_bed_tuples(Path(query_path))
    print(f"\n  Query: {Path(query_path).name} ({len(query_regions)} regions)")

    # BM25 search
    print("\n  Searching Qdrant (BM25+IDF)...")
    qdrant_results, search_time = search_enrichment(client, collection, bm25, query_path, limit=20)

    # LOLA support
    print("  Computing LOLA support...")
    lola_results, lola_time = compute_lola_support(query_path, ref_beds)
    print(f"  done ({lola_time:.1f}s)")

    # Output
    n_indexed = info.points_count
    print_results(qdrant_results, search_time, lola_results, lola_time, bedset_id)
    save_json(query_path, bedset_id, qdrant_results, lola_results, search_time, lola_time, index_time, n_indexed, results_path)


if __name__ == "__main__":
    main()
