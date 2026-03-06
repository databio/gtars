//! Brute-force KNN and SNN graph construction.
//!
//! Builds a k-nearest-neighbor graph from PCA/LSI embeddings, then derives
//! a Shared Nearest Neighbor (SNN) graph with Jaccard-weighted edges for
//! downstream clustering. Matches Seurat's `FindNeighbors` behavior.

use std::collections::HashSet;

use anyhow::{Result, bail};
use ndarray::Array2;

use crate::types::{KnnGraph, SnnGraph};

/// Build a brute-force KNN graph from a cells × components embedding matrix.
///
/// For each cell, computes Euclidean distance to all other cells and retains
/// the k nearest neighbors. Uses `select_nth_unstable_by` for O(n) partial
/// sort per cell.
///
/// Parallelizes over cells with rayon when `parallel` feature is enabled.
pub fn build_knn(embedding: &Array2<f64>, k: usize) -> Result<KnnGraph> {
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

/// Find k nearest neighbors for a single cell.
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

/// Build an SNN (Shared Nearest Neighbor) graph from a KNN graph.
///
/// For each pair of cells (i, j) connected by KNN in either direction,
/// computes the Jaccard index: `|KNN(i) ∩ KNN(j)| / |KNN(i) ∪ KNN(j)|`.
/// Edges with Jaccard below `prune_threshold` are dropped.
///
/// Matches Seurat's `FindNeighbors` SNN construction with `prune.SNN` parameter.
pub fn build_snn(knn: &KnnGraph, prune_threshold: f64) -> SnnGraph {
    let n = knn.n_cells;
    let k = knn.k;

    // Build HashSets for O(1) neighbor lookup
    let neighbor_sets: Vec<HashSet<usize>> = knn
        .indices
        .iter()
        .map(|nn| nn.iter().copied().collect())
        .collect();

    // Collect all unique unordered pairs connected by KNN in either direction.
    // KNN is asymmetric: j ∈ KNN(i) does not imply i ∈ KNN(j).
    // Seurat symmetrizes, so we consider the union of both directions.
    let mut pairs: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..n {
        for &j in &knn.indices[i] {
            if i != j {
                pairs.insert((i.min(j), i.max(j)));
            }
        }
    }

    let mut edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    for (i, j) in pairs {
        let intersection = neighbor_sets[i].intersection(&neighbor_sets[j]).count();
        let union_size = 2 * k - intersection; // |A ∪ B| = |A| + |B| - |A ∩ B|
        let jaccard = intersection as f64 / union_size as f64;

        if jaccard > prune_threshold {
            edges[i].push((j, jaccard));
            edges[j].push((i, jaccard));
        }
    }

    SnnGraph { n_cells: n, edges }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_knn_basic() {
        // 4 points in 2D: two pairs close together
        let embedding = array![
            [0.0, 0.0], // 0
            [0.1, 0.0], // 1 — close to 0
            [5.0, 0.0], // 2
            [5.1, 0.0], // 3 — close to 2
        ];
        let knn = build_knn(&embedding, 2).unwrap();
        assert_eq!(knn.n_cells, 4);
        assert_eq!(knn.k, 2);

        // Cell 0's nearest neighbor should be cell 1
        assert_eq!(knn.indices[0][0], 1);
        // Cell 3's nearest neighbor should be cell 2
        assert_eq!(knn.indices[3][0], 2);
    }

    #[test]
    fn test_knn_distances_sorted() {
        let embedding = array![
            [0.0, 0.0],
            [1.0, 0.0],
            [3.0, 0.0],
            [6.0, 0.0],
        ];
        let knn = build_knn(&embedding, 3).unwrap();
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
        // With k=3 out of 7 other cells, intra-cluster neighbors dominate.
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
        let knn = build_knn(&embedding, 3).unwrap();

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
        let knn = build_knn(&embedding, 3).unwrap();

        // High threshold should remove weak edges
        let snn_strict = build_snn(&knn, 0.5);
        // Low threshold should keep more edges
        let snn_loose = build_snn(&knn, 0.0);

        let strict_edges: usize = snn_strict.edges.iter().map(|e| e.len()).sum();
        let loose_edges: usize = snn_loose.edges.iter().map(|e| e.len()).sum();
        assert!(strict_edges <= loose_edges);
    }
}
