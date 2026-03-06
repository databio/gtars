//! Silhouette score computation for cluster quality assessment.
//!
//! For each cell i with cluster label c:
//!   a(i) = mean Euclidean distance to other cells in cluster c
//!   b(i) = min over other clusters c' of mean distance to cells in c'
//!   s(i) = (b - a) / max(a, b)
//!
//! Scores range from -1 (misclassified) to +1 (well-clustered).

use ndarray::Array2;

/// Compute per-cell silhouette scores from embeddings and cluster assignments.
///
/// Cells in singleton clusters receive a score of 0.0.
/// Returns a vector of length `n_cells`.
pub fn silhouette_scores(embedding: &Array2<f64>, clusters: &[u32]) -> Vec<f64> {
    let n_cells = embedding.nrows();
    let n_dims = embedding.ncols();

    if n_cells == 0 {
        return Vec::new();
    }

    let n_clusters = clusters.iter().max().map(|&m| m as usize + 1).unwrap_or(0);
    let mut cluster_members: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (i, &c) in clusters.iter().enumerate() {
        cluster_members[c as usize].push(i);
    }

    let mut scores = vec![0.0f64; n_cells];

    for i in 0..n_cells {
        let my_cluster = clusters[i] as usize;
        if cluster_members[my_cluster].len() <= 1 {
            continue; // singleton → 0.0
        }

        // a(i): mean distance to same-cluster cells
        let a = mean_dist_to_group(embedding, i, &cluster_members[my_cluster], n_dims);

        // b(i): min mean distance to any other cluster
        let mut b = f64::INFINITY;
        for c in 0..n_clusters {
            if c == my_cluster || cluster_members[c].is_empty() {
                continue;
            }
            let d = mean_dist_to_group(embedding, i, &cluster_members[c], n_dims);
            b = b.min(d);
        }

        let max_ab = a.max(b);
        scores[i] = if max_ab > 0.0 { (b - a) / max_ab } else { 0.0 };
    }

    scores
}

/// Mean silhouette score across all cells.
pub fn mean_silhouette(scores: &[f64]) -> f64 {
    if scores.is_empty() {
        return 0.0;
    }
    scores.iter().sum::<f64>() / scores.len() as f64
}

/// Mean Euclidean distance from cell `i` to all cells in `members` (excluding self).
fn mean_dist_to_group(
    embedding: &Array2<f64>,
    i: usize,
    members: &[usize],
    n_dims: usize,
) -> f64 {
    let mut total = 0.0f64;
    let mut count = 0usize;
    for &j in members {
        if j == i {
            continue;
        }
        let mut d2 = 0.0;
        for d in 0..n_dims {
            let diff = embedding[[i, d]] - embedding[[j, d]];
            d2 += diff * diff;
        }
        total += d2.sqrt();
        count += 1;
    }
    if count > 0 {
        total / count as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_perfectly_separated() {
        // Two well-separated clusters
        let embedding = array![
            [0.0, 0.0],
            [0.1, 0.0],
            [0.0, 0.1],
            [10.0, 10.0],
            [10.1, 10.0],
            [10.0, 10.1],
        ];
        let clusters = vec![0, 0, 0, 1, 1, 1];
        let scores = silhouette_scores(&embedding, &clusters);
        let avg = mean_silhouette(&scores);

        // All scores should be close to 1.0 for well-separated data
        for &s in &scores {
            assert!(s > 0.9, "expected high silhouette, got {s}");
        }
        assert!(avg > 0.9);
    }

    #[test]
    fn test_singleton_cluster() {
        let embedding = array![[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]];
        let clusters = vec![0, 1, 2]; // each cell alone
        let scores = silhouette_scores(&embedding, &clusters);
        for &s in &scores {
            assert_eq!(s, 0.0);
        }
    }

    #[test]
    fn test_empty() {
        let embedding = Array2::<f64>::zeros((0, 2));
        let scores = silhouette_scores(&embedding, &[]);
        assert!(scores.is_empty());
    }
}
