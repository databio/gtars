//! KNN and SNN graph construction for single-cell analysis.
//!
//! Builds a k-nearest-neighbor graph from PCA/LSI embeddings using approximate
//! nearest neighbors (NN-descent algorithm), then derives a Shared Nearest
//! Neighbor (SNN) graph with Jaccard-weighted edges for downstream clustering.
//! Matches Seurat's `FindNeighbors` behavior.
//!
//! NN-descent is the algorithm used by pynndescent/UMAP. It starts with a random
//! KNN graph and iteratively improves it by checking "neighbors of neighbors" until
//! convergence.

use std::collections::HashSet;

use anyhow::{Result, bail};
use ndarray::Array2;
use rand::Rng;
use rand::rngs::SmallRng;
use rand::SeedableRng;

use crate::types::{KnnGraph, SnnGraph};

// ---------------------------------------------------------------------------
// NN-Descent (approximate nearest neighbors)
// ---------------------------------------------------------------------------

/// Compute squared Euclidean distance between two rows of a matrix.
#[inline]
fn dist_sq(data: &Array2<f64>, i: usize, j: usize) -> f64 {
    let n_dims = data.ncols();
    let mut d2 = 0.0f64;
    for d in 0..n_dims {
        let diff = data[[i, d]] - data[[j, d]];
        d2 += diff * diff;
    }
    d2
}

/// A single entry in a KNN heap: (distance, index).
/// We use a max-heap (largest distance at top) so we can efficiently evict the
/// farthest current neighbor when a closer one is found.
#[derive(Clone, Copy)]
struct HeapEntry {
    dist: f64,
    idx: usize,
    is_new: bool, // NN-descent flag: was this entry added in the current iteration?
}

/// Fixed-size max-heap for maintaining k nearest neighbors.
struct KnnHeap {
    entries: Vec<HeapEntry>,
    k: usize,
}

impl KnnHeap {
    fn new(k: usize) -> Self {
        KnnHeap {
            entries: Vec::with_capacity(k),
            k,
        }
    }

    /// Try to insert a new neighbor. Returns true if it was inserted (improved the heap).
    fn insert(&mut self, idx: usize, dist: f64) -> bool {
        // Don't insert duplicates
        if self.entries.iter().any(|e| e.idx == idx) {
            return false;
        }

        if self.entries.len() < self.k {
            self.entries.push(HeapEntry {
                dist,
                idx,
                is_new: true,
            });
            return true;
        }

        // Find the farthest current neighbor
        let (max_pos, max_dist) = self
            .entries
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.dist.partial_cmp(&b.1.dist).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, e)| (i, e.dist))
            .unwrap();

        if dist < max_dist {
            self.entries[max_pos] = HeapEntry {
                dist,
                idx,
                is_new: true,
            };
            return true;
        }

        false
    }

    /// Get sorted indices and distances.
    fn to_sorted(&self) -> (Vec<usize>, Vec<f64>) {
        let mut sorted: Vec<(f64, usize)> = self.entries.iter().map(|e| (e.dist, e.idx)).collect();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let indices: Vec<usize> = sorted.iter().map(|&(_, i)| i).collect();
        let distances: Vec<f64> = sorted.iter().map(|&(d, _)| d.sqrt()).collect();
        (indices, distances)
    }

    /// Get indices marked as "new" (added in the current iteration).
    fn new_indices(&self) -> Vec<usize> {
        self.entries
            .iter()
            .filter(|e| e.is_new)
            .map(|e| e.idx)
            .collect()
    }

    /// Mark all entries as "old" (not new).
    fn mark_all_old(&mut self) {
        for e in &mut self.entries {
            e.is_new = false;
        }
    }

    /// Get all neighbor indices.
    fn all_indices(&self) -> Vec<usize> {
        self.entries.iter().map(|e| e.idx).collect()
    }
}

/// Build an approximate KNN graph using NN-descent.
///
/// NN-descent starts with a random KNN and iteratively refines it by exploring
/// "neighbors of neighbors". This is the same algorithm used by pynndescent (UMAP).
/// It produces high-quality approximate nearest neighbors without tree structures.
///
/// Parameters:
/// - `max_iterations`: Maximum NN-descent iterations (default: 10)
/// - `delta`: Convergence threshold — stop when fewer than delta*n*k updates occur (default: 0.001)
/// - `seed`: Random seed for reproducibility
pub fn build_knn(embedding: &Array2<f64>, k: usize) -> Result<KnnGraph> {
    let n_cells = embedding.nrows();
    // Adaptive iteration count following pynndescent's formula:
    // max_iterations = max(5, round(log2(n_cells)))
    // This scales naturally: 2K cells → 11 iters, 10K → 14, 100K → 17
    let max_iterations = 5.max((n_cells as f64).log2().round() as usize);
    build_knn_nndescent(embedding, k, max_iterations, 0.001, 42)
}

/// Build an approximate KNN graph using NN-descent with full parameter control.
pub fn build_knn_nndescent(
    embedding: &Array2<f64>,
    k: usize,
    max_iterations: usize,
    delta: f64,
    seed: u64,
) -> Result<KnnGraph> {
    let n_cells = embedding.nrows();
    if k == 0 {
        bail!("k must be at least 1");
    }
    if k >= n_cells {
        bail!("k ({k}) must be less than n_cells ({n_cells})");
    }

    let mut rng = SmallRng::seed_from_u64(seed);

    // Initialize with random neighbors
    let mut heaps: Vec<KnnHeap> = Vec::with_capacity(n_cells);
    for i in 0..n_cells {
        let mut heap = KnnHeap::new(k);
        // Pick k random distinct neighbors
        let mut used: HashSet<usize> = HashSet::new();
        used.insert(i);
        while heap.entries.len() < k {
            let j = rng.random_range(0..n_cells);
            if used.insert(j) {
                let d = dist_sq(embedding, i, j);
                heap.insert(j, d);
            }
        }
        heaps.push(heap);
    }

    // NN-descent iterations: only compare pairs of cells that share a
    // neighbor (the standard NN-descent "local join"). This is more conservative
    // than full neighbor-of-neighbor exploration.
    let threshold = (delta * n_cells as f64 * k as f64) as usize;
    // Sample rate for local join. Following pynndescent, use full exploration
    // (sample_rate=1.0) and rely on convergence-based early stopping via delta.
    let sample_rate = 1.0;

    for _iter in 0..max_iterations {
        let mut n_updates = 0usize;

        // Collect new neighbors for each cell
        let new_neighbors: Vec<Vec<usize>> = heaps.iter().map(|h| h.new_indices()).collect();

        // Mark all as old for next iteration
        for heap in &mut heaps {
            heap.mark_all_old();
        }

        // Local join: for each cell i with new neighbors, compare pairs of
        // i's neighbors with each other. If neighbor u and neighbor v of cell i
        // are close to each other, they might be each other's KNN.
        for i in 0..n_cells {
            if new_neighbors[i].is_empty() {
                continue;
            }

            let all_nn: Vec<usize> = heaps[i].all_indices();

            // Compare new neighbors with all neighbors
            for &u in &new_neighbors[i] {
                for &v in &all_nn {
                    if u == v {
                        continue;
                    }
                    // Subsample to control exploration rate
                    if rng.random::<f64>() > sample_rate {
                        continue;
                    }
                    let d = dist_sq(embedding, u, v);
                    if heaps[u].insert(v, d) {
                        n_updates += 1;
                    }
                    if heaps[v].insert(u, d) {
                        n_updates += 1;
                    }
                }
            }
        }

        if n_updates <= threshold {
            break;
        }
    }

    // Extract results
    let mut indices = Vec::with_capacity(n_cells);
    let mut distances = Vec::with_capacity(n_cells);
    for heap in &heaps {
        let (idx, dist) = heap.to_sorted();
        indices.push(idx);
        distances.push(dist);
    }

    Ok(KnnGraph {
        n_cells,
        k,
        indices,
        distances,
    })
}

/// Build an exact brute-force KNN graph.
///
/// Computes exact Euclidean distances between all pairs of cells. Only use for
/// small datasets or when exact results are required.
pub fn build_knn_exact(embedding: &Array2<f64>, k: usize) -> Result<KnnGraph> {
    let n_cells = embedding.nrows();

    if k == 0 {
        bail!("k must be at least 1");
    }
    if k >= n_cells {
        bail!("k ({k}) must be less than n_cells ({n_cells})");
    }

    #[cfg(feature = "parallel")]
    let (indices, distances) = {
        use rayon::prelude::*;
        let results: Vec<(Vec<usize>, Vec<f64>)> = (0..n_cells)
            .into_par_iter()
            .map(|i| knn_for_cell(embedding, i, k))
            .collect();
        results.into_iter().unzip()
    };

    #[cfg(not(feature = "parallel"))]
    let (indices, distances) = {
        let mut all_indices = Vec::with_capacity(n_cells);
        let mut all_distances = Vec::with_capacity(n_cells);
        for i in 0..n_cells {
            let (idx, dist) = knn_for_cell(embedding, i, k);
            all_indices.push(idx);
            all_distances.push(dist);
        }
        (all_indices, all_distances)
    };

    Ok(KnnGraph {
        n_cells,
        k,
        indices,
        distances,
    })
}

/// Find k nearest neighbors for a single cell (brute-force).
fn knn_for_cell(
    embedding: &Array2<f64>,
    i: usize,
    k: usize,
) -> (Vec<usize>, Vec<f64>) {
    let n_cells = embedding.nrows();
    let n_dims = embedding.ncols();

    let mut dists: Vec<(usize, f64)> = Vec::with_capacity(n_cells - 1);
    for j in 0..n_cells {
        if j == i {
            continue;
        }
        let mut d2 = 0.0f64;
        for dim in 0..n_dims {
            let diff = embedding[[i, dim]] - embedding[[j, dim]];
            d2 += diff * diff;
        }
        dists.push((j, d2.sqrt()));
    }

    // O(n) partial sort to find k smallest
    dists.select_nth_unstable_by(k - 1, |a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    dists.truncate(k);
    dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let indices: Vec<usize> = dists.iter().map(|&(j, _)| j).collect();
    let distances: Vec<f64> = dists.iter().map(|&(_, d)| d).collect();
    (indices, distances)
}

// ---------------------------------------------------------------------------
// SNN graph construction
// ---------------------------------------------------------------------------

/// Build an SNN (Shared Nearest Neighbor) graph from a KNN graph.
///
/// Matches Seurat's `FindNeighbors` SNN construction: computes the Jaccard
/// index for ALL pairs of cells that share at least one KNN neighbor (not just
/// direct KNN pairs). This is equivalent to Seurat's matrix-multiply approach:
/// `SNN = nn_indicator %*% t(nn_indicator)`, then normalizing to Jaccard.
///
/// Edges with Jaccard below `prune_threshold` are dropped.
pub fn build_snn(knn: &KnnGraph, prune_threshold: f64) -> SnnGraph {
    let n = knn.n_cells;
    let k = knn.k;

    // Build neighbor sets: self + first k-1 true neighbors = k elements.
    // Matches Seurat's convention where k.param=20 returns 20 items including
    // self (19 true neighbors + self).
    let neighbor_sets: Vec<HashSet<usize>> = knn
        .indices
        .iter()
        .enumerate()
        .map(|(i, nn)| {
            let mut set: HashSet<usize> = nn.iter().take(k - 1).copied().collect();
            set.insert(i);
            set
        })
        .collect();

    // Build inverted index: for each cell c, which cells have c in their
    // neighbor set? This lets us find ALL pairs sharing at least one neighbor.
    let mut inv_index: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, set) in neighbor_sets.iter().enumerate() {
        for &neighbor in set {
            inv_index[neighbor].push(i);
        }
    }

    // For each cell i, find all cells j that share at least one neighbor.
    // Count shared neighbors via the inverted index, then compute Jaccard.
    let mut edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    for i in 0..n {
        // Count how many neighbors cell i shares with each other cell
        let mut shared_counts: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for &neighbor in &neighbor_sets[i] {
            for &j in &inv_index[neighbor] {
                if j > i {
                    // Only count upper triangle to avoid duplicates
                    *shared_counts.entry(j).or_insert(0) += 1;
                }
            }
        }

        // Compute Jaccard for each pair with shared neighbors
        for (j, intersection) in shared_counts {
            // Each set has k elements, |A ∪ B| = |A| + |B| - |A ∩ B|
            let union_size = 2 * k - intersection;
            let jaccard = intersection as f64 / union_size as f64;

            if jaccard >= prune_threshold {
                edges[i].push((j, jaccard));
                edges[j].push((i, jaccard));
            }
        }
    }

    SnnGraph { n_cells: n, edges }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_knn_exact_basic() {
        // 4 points in 2D: two pairs close together
        let embedding = array![
            [0.0, 0.0], // 0
            [0.1, 0.0], // 1 — close to 0
            [5.0, 0.0], // 2
            [5.1, 0.0], // 3 — close to 2
        ];
        let knn = build_knn_exact(&embedding, 2).unwrap();
        assert_eq!(knn.n_cells, 4);
        assert_eq!(knn.k, 2);

        // Cell 0's nearest neighbor should be cell 1
        assert_eq!(knn.indices[0][0], 1);
        // Cell 3's nearest neighbor should be cell 2
        assert_eq!(knn.indices[3][0], 2);
    }

    #[test]
    fn test_knn_nndescent_basic() {
        // Well-separated clusters — NN-descent should find correct neighbors
        let embedding = array![
            [0.0, 0.0],  // cluster A
            [0.1, 0.0],
            [0.0, 0.1],
            [0.1, 0.1],
            [10.0, 0.0], // cluster B
            [10.1, 0.0],
            [10.0, 0.1],
            [10.1, 0.1],
        ];
        let knn = build_knn_nndescent(&embedding, 3, 20, 0.001, 42).unwrap();
        assert_eq!(knn.k, 3);

        // Cell 0's neighbors should all be in cluster A
        for &neighbor in &knn.indices[0] {
            assert!(
                neighbor <= 3,
                "cell 0 neighbor {neighbor} should be in cluster A"
            );
        }
        // Cell 4's neighbors should all be in cluster B
        for &neighbor in &knn.indices[4] {
            assert!(
                neighbor >= 4,
                "cell 4 neighbor {neighbor} should be in cluster B"
            );
        }
    }

    #[test]
    fn test_knn_nndescent_recall() {
        // Generate a larger random dataset and verify high recall
        let mut rng = SmallRng::seed_from_u64(123);
        let n = 200;
        let d = 10;
        let mut data = Array2::zeros((n, d));
        for i in 0..n {
            for j in 0..d {
                data[[i, j]] = rng.random::<f64>();
            }
        }

        let k = 10;
        let exact = build_knn_exact(&data, k).unwrap();
        let approx = build_knn_nndescent(&data, k, 10, 0.001, 42).unwrap();

        // Compute recall
        let mut total_found = 0;
        let total = n * k;
        for i in 0..n {
            let exact_set: HashSet<usize> = exact.indices[i].iter().copied().collect();
            for &j in &approx.indices[i] {
                if exact_set.contains(&j) {
                    total_found += 1;
                }
            }
        }
        let recall = total_found as f64 / total as f64;
        assert!(
            recall > 0.80,
            "NN-descent recall {recall:.3} is too low (expected > 0.80)"
        );
    }

    #[test]
    fn test_knn_distances_sorted() {
        let embedding = array![
            [0.0, 0.0],
            [1.0, 0.0],
            [3.0, 0.0],
            [6.0, 0.0],
        ];
        let knn = build_knn_exact(&embedding, 3).unwrap();
        for cell in 0..4 {
            for w in knn.distances[cell].windows(2) {
                assert!(w[0] <= w[1], "distances not sorted for cell {cell}");
            }
        }
    }

    #[test]
    fn test_knn_k_too_large() {
        let embedding = array![[0.0], [1.0], [2.0]];
        assert!(build_knn(&embedding, 3).is_err());
    }

    #[test]
    fn test_snn_jaccard() {
        // 8 cells: 4 in each cluster, well separated. k=3.
        let embedding = array![
            [0.0, 0.0],  // cluster A
            [0.1, 0.0],
            [0.0, 0.1],
            [0.1, 0.1],
            [10.0, 0.0], // cluster B
            [10.1, 0.0],
            [10.0, 0.1],
            [10.1, 0.1],
        ];
        let knn = build_knn_exact(&embedding, 3).unwrap();

        // All of cell 0's neighbors should be in cluster A (cells 1,2,3)
        for &neighbor in &knn.indices[0] {
            assert!(neighbor <= 3, "cell 0 should only have cluster A neighbors");
        }

        let snn = build_snn(&knn, 0.0); // no pruning

        // Intra-cluster edges should exist with positive weights
        let cell0_edges: Vec<(usize, f64)> = snn.edges[0].clone();
        assert!(!cell0_edges.is_empty(), "cell 0 should have SNN edges");

        // All of cell 0's SNN edges should be to cluster A
        for &(j, w) in &cell0_edges {
            assert!(j <= 3, "cell 0 SNN edge to {j} crosses clusters");
            assert!(w > 0.0, "SNN weight should be positive");
        }
    }

    #[test]
    fn test_snn_pruning() {
        let embedding = array![
            [0.0, 0.0],
            [0.1, 0.0],
            [0.2, 0.0],
            [10.0, 0.0],
            [10.1, 0.0],
        ];
        let knn = build_knn_exact(&embedding, 3).unwrap();

        // High threshold should remove weak edges
        let snn_strict = build_snn(&knn, 0.5);
        // Low threshold should keep more edges
        let snn_loose = build_snn(&knn, 0.0);

        let strict_edges: usize = snn_strict.edges.iter().map(|e| e.len()).sum();
        let loose_edges: usize = snn_loose.edges.iter().map(|e| e.len()).sum();
        assert!(strict_edges <= loose_edges);
    }
}
