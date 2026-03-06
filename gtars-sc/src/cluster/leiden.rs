//! Leiden community detection with resolution parameter.
//!
//! Implements a multi-level Louvain/Leiden algorithm (Traag, Waltman, van Eck 2019)
//! with the modularity quality function:
//!
//!   Q = Σ_c [ L_c/m − γ (K_c / 2m)² ]
//!
//! where L_c = sum of edge weights within community c,
//!       K_c = sum of node strengths in community c,
//!       m = total edge weight / 2,
//!       γ = resolution parameter.
//!
//! Higher resolution → more clusters; lower resolution → fewer clusters.
//!
//! The multi-level approach:
//! 1. Run local moves on the current graph
//! 2. Coarsen the graph: communities become super-nodes
//! 3. Run local moves on the coarsened graph
//! 4. Repeat until no further reduction

use anyhow::Result;
use rand::seq::SliceRandom;
use rand::SeedableRng;

use crate::types::{ClusterResult, SnnGraph};

/// Run Leiden community detection on an SNN graph.
///
/// # Arguments
/// * `graph` — Weighted SNN graph (Jaccard-weighted edges)
/// * `resolution` — Controls granularity (Seurat default 0.8)
/// * `max_iterations` — Maximum number of local move iterations per level (default 10)
pub fn leiden_clustering(
    graph: &SnnGraph,
    resolution: f64,
    max_iterations: usize,
) -> Result<ClusterResult> {
    let n = graph.n_cells;
    if n == 0 {
        return Ok(ClusterResult {
            assignments: Vec::new(),
            n_clusters: 0,
            quality: Some(0.0),
        });
    }

    // Precompute node strengths on the original graph (for final modularity)
    let orig_node_strength: Vec<f64> = (0..n)
        .map(|i| graph.edges[i].iter().map(|&(_, w)| w).sum())
        .collect();
    let orig_total_weight: f64 = orig_node_strength.iter().sum::<f64>() / 2.0;

    if orig_total_weight == 0.0 {
        return Ok(ClusterResult {
            assignments: (0..n as u32).collect(),
            n_clusters: n as u32,
            quality: Some(0.0),
        });
    }

    let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

    // Global assignment: original node → current community
    let mut global_community: Vec<u32> = (0..n as u32).collect();

    // Start with the original graph
    let mut current_graph = SnnGraph {
        n_cells: graph.n_cells,
        edges: graph.edges.clone(),
    };

    // Multi-level loop: local moves → coarsen → repeat
    for _level in 0..20 {
        let n_nodes = current_graph.n_cells;

        // Run local move phase on current graph
        let level_community = local_move_phase(
            &current_graph,
            resolution,
            max_iterations,
            &mut rng,
        );

        // Count communities
        let n_comm = level_community.iter().max().map(|&m| m as usize + 1).unwrap_or(0);

        // If no reduction (each node stayed in its own community), we're done
        if n_comm >= n_nodes {
            break;
        }

        // Update global assignment: map through the level assignment
        for g in global_community.iter_mut() {
            *g = level_community[*g as usize];
        }

        // Coarsen the graph: communities become super-nodes
        current_graph = coarsen_graph(&current_graph, &level_community, n_comm);
    }

    // Final renumbering
    renumber_communities(&mut global_community);
    let n_clusters = global_community.iter().max().map(|&m| m + 1).unwrap_or(0);

    // Compute final modularity on the original graph
    let quality = compute_modularity(
        graph,
        &global_community,
        &orig_node_strength,
        orig_total_weight,
        resolution,
    );

    Ok(ClusterResult {
        assignments: global_community,
        n_clusters,
        quality: Some(quality),
    })
}

/// Run local move phase: iterate nodes, move each to the neighboring community
/// that maximizes ΔQ. Returns a community assignment vector (renumbered 0..k-1).
fn local_move_phase(
    graph: &SnnGraph,
    resolution: f64,
    max_iterations: usize,
    rng: &mut rand::rngs::SmallRng,
) -> Vec<u32> {
    let n = graph.n_cells;
    let mut community: Vec<u32> = (0..n as u32).collect();

    // Precompute node strengths for this graph
    let node_strength: Vec<f64> = (0..n)
        .map(|i| graph.edges[i].iter().map(|&(_, w)| w).sum())
        .collect();

    let total_weight: f64 = node_strength.iter().sum::<f64>() / 2.0;
    if total_weight == 0.0 {
        return community;
    }

    // Community-level strength aggregate
    let mut comm_strength: Vec<f64> = node_strength.clone();

    for _iter in 0..max_iterations {
        let mut improved = false;

        // Iterate nodes in random order
        let mut order: Vec<usize> = (0..n).collect();
        order.shuffle(rng);

        for &node in &order {
            let current_comm = community[node] as usize;

            // Compute edges from this node to each neighboring community.
            // Skip self-loops: they represent internal edges from coarsening
            // and should not influence community move decisions.
            let mut neighbor_comm_weight: Vec<(u32, f64)> = Vec::new();
            for &(j, w) in &graph.edges[node] {
                if j == node {
                    continue;
                }
                let j_comm = community[j];
                if let Some(entry) = neighbor_comm_weight.iter_mut().find(|e| e.0 == j_comm) {
                    entry.1 += w;
                } else {
                    neighbor_comm_weight.push((j_comm, w));
                }
            }

            // Weight of edges from node to its own community
            let e_to_own = neighbor_comm_weight
                .iter()
                .find(|e| e.0 == current_comm as u32)
                .map(|e| e.1)
                .unwrap_or(0.0);

            let s_i = node_strength[node];

            // Try each neighboring community
            let mut best_comm = current_comm as u32;
            let mut best_delta = 0.0f64;

            for &(target_comm, e_to_target) in &neighbor_comm_weight {
                if target_comm as usize == current_comm {
                    continue;
                }

                // ΔQ = (e_to_target - e_to_own)/m - γ * s_i * (K_target - K_own + s_i) / (2m²)
                let k_target = comm_strength[target_comm as usize];
                let k_own = comm_strength[current_comm] - s_i;

                let delta = (e_to_target - e_to_own) / total_weight
                    - resolution * s_i * (k_target - k_own)
                        / (2.0 * total_weight * total_weight);

                if delta > best_delta {
                    best_delta = delta;
                    best_comm = target_comm;
                }
            }

            // Move node if improvement found
            if best_comm != current_comm as u32 {
                comm_strength[current_comm] -= s_i;
                comm_strength[best_comm as usize] += s_i;
                community[node] = best_comm;
                improved = true;
            }
        }

        if !improved {
            break;
        }

        // Renumber communities to be contiguous 0..k
        renumber_communities(&mut community);

        // Recompute comm_strength after renumbering
        let max_c = *community.iter().max().unwrap_or(&0) as usize + 1;
        comm_strength = vec![0.0; max_c];
        for i in 0..n {
            comm_strength[community[i] as usize] += node_strength[i];
        }
    }

    // Final renumber for this level
    renumber_communities(&mut community);
    community
}

/// Build a coarsened graph where each community becomes a super-node.
///
/// Edge weights between super-nodes are the sum of original edge weights
/// between nodes in those communities. Self-loops represent internal edges
/// and are included to preserve node strength invariants.
fn coarsen_graph(graph: &SnnGraph, community: &[u32], n_comm: usize) -> SnnGraph {
    let mut edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_comm];

    // Iterate all edges in the adjacency list (bidirectional).
    // Each original edge (i,j) with weight w appears as (j,w) in adj[i]
    // and (i,w) in adj[j]. By processing both, we naturally get the correct
    // bidirectional representation and self-loop weights in the coarsened graph.
    for i in 0..graph.n_cells {
        let ci = community[i] as usize;
        for &(j, w) in &graph.edges[i] {
            let cj = community[j] as usize;
            if let Some(entry) = edges[ci].iter_mut().find(|e| e.0 == cj) {
                entry.1 += w;
            } else {
                edges[ci].push((cj, w));
            }
        }
    }

    SnnGraph {
        n_cells: n_comm,
        edges,
    }
}

/// Renumber communities to contiguous 0..k.
fn renumber_communities(community: &mut [u32]) {
    let mut seen = std::collections::HashMap::new();
    let mut next_id = 0u32;
    for c in community.iter_mut() {
        let new_id = *seen.entry(*c).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
        *c = new_id;
    }
}

/// Compute modularity Q = Σ_c [L_c/m − γ (K_c / 2m)²].
fn compute_modularity(
    graph: &SnnGraph,
    community: &[u32],
    node_strength: &[f64],
    total_weight: f64,
    resolution: f64,
) -> f64 {
    let n = graph.n_cells;
    let n_comm = community.iter().max().map(|&m| m as usize + 1).unwrap_or(0);
    let two_m = 2.0 * total_weight;

    let mut l_c = vec![0.0f64; n_comm]; // internal weight
    let mut k_c = vec![0.0f64; n_comm]; // total strength

    for i in 0..n {
        let ci = community[i] as usize;
        k_c[ci] += node_strength[i];
        for &(j, w) in &graph.edges[i] {
            if j > i && community[j] as usize == ci {
                l_c[ci] += w;
            }
        }
    }

    let mut q = 0.0f64;
    for c in 0..n_comm {
        q += l_c[c] / total_weight - resolution * (k_c[c] / two_m).powi(2);
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an SNN graph from an edge list.
    fn graph_from_edges(n: usize, edges: &[(usize, usize, f64)]) -> SnnGraph {
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for &(i, j, w) in edges {
            adj[i].push((j, w));
            adj[j].push((i, w));
        }
        SnnGraph {
            n_cells: n,
            edges: adj,
        }
    }

    #[test]
    fn test_two_cliques() {
        // Two fully-connected triangles connected by a weak edge
        let graph = graph_from_edges(
            6,
            &[
                // Clique A: 0-1-2
                (0, 1, 1.0),
                (0, 2, 1.0),
                (1, 2, 1.0),
                // Clique B: 3-4-5
                (3, 4, 1.0),
                (3, 5, 1.0),
                (4, 5, 1.0),
                // Weak bridge
                (2, 3, 0.1),
            ],
        );

        let result = leiden_clustering(&graph, 1.0, 10).unwrap();
        assert_eq!(result.n_clusters, 2, "expected 2 clusters, got {}", result.n_clusters);

        // Nodes 0,1,2 should be in the same cluster
        assert_eq!(result.assignments[0], result.assignments[1]);
        assert_eq!(result.assignments[1], result.assignments[2]);

        // Nodes 3,4,5 should be in the same cluster
        assert_eq!(result.assignments[3], result.assignments[4]);
        assert_eq!(result.assignments[4], result.assignments[5]);

        // The two groups should be in different clusters
        assert_ne!(result.assignments[0], result.assignments[3]);

        assert!(result.quality.unwrap() > 0.0, "modularity should be positive");
    }

    #[test]
    fn test_all_disconnected() {
        let graph = SnnGraph {
            n_cells: 4,
            edges: vec![Vec::new(); 4],
        };
        let result = leiden_clustering(&graph, 1.0, 10).unwrap();
        assert_eq!(result.n_clusters, 4); // each node its own community
    }

    #[test]
    fn test_single_clique() {
        // Fully connected 4-node graph at low resolution → should be 1 cluster
        let graph = graph_from_edges(
            4,
            &[
                (0, 1, 1.0),
                (0, 2, 1.0),
                (0, 3, 1.0),
                (1, 2, 1.0),
                (1, 3, 1.0),
                (2, 3, 1.0),
            ],
        );
        let result = leiden_clustering(&graph, 0.5, 10).unwrap();
        assert_eq!(result.n_clusters, 1);
    }

    #[test]
    fn test_empty_graph() {
        let graph = SnnGraph {
            n_cells: 0,
            edges: Vec::new(),
        };
        let result = leiden_clustering(&graph, 1.0, 10).unwrap();
        assert_eq!(result.n_clusters, 0);
    }

    #[test]
    fn test_resolution_affects_clusters() {
        // Ring of 6 nodes: should split differently at different resolutions
        let graph = graph_from_edges(
            6,
            &[
                (0, 1, 1.0),
                (1, 2, 1.0),
                (2, 3, 1.0),
                (3, 4, 1.0),
                (4, 5, 1.0),
                (5, 0, 1.0),
            ],
        );
        let low_res = leiden_clustering(&graph, 0.1, 10).unwrap();
        let high_res = leiden_clustering(&graph, 5.0, 10).unwrap();
        assert!(
            high_res.n_clusters >= low_res.n_clusters,
            "higher resolution should produce more or equal clusters"
        );
    }

    #[test]
    fn test_all_nodes_assigned() {
        let graph = graph_from_edges(
            5,
            &[
                (0, 1, 1.0),
                (1, 2, 0.5),
                (2, 3, 1.0),
                (3, 4, 0.5),
            ],
        );
        let result = leiden_clustering(&graph, 1.0, 10).unwrap();
        assert_eq!(result.assignments.len(), 5);
        // All cluster IDs should be < n_clusters
        for &c in &result.assignments {
            assert!(c < result.n_clusters);
        }
    }

    #[test]
    fn test_multi_level_two_large_cliques() {
        // Two cliques of 20 nodes each, weakly connected
        // This requires multi-level aggregation to resolve correctly
        let mut edges = Vec::new();
        // Clique A: nodes 0-19, fully connected
        for i in 0..20usize {
            for j in (i + 1)..20 {
                edges.push((i, j, 1.0));
            }
        }
        // Clique B: nodes 20-39, fully connected
        for i in 20..40usize {
            for j in (i + 1)..40 {
                edges.push((i, j, 1.0));
            }
        }
        // Weak bridge
        edges.push((19, 20, 0.01));

        let graph = graph_from_edges(40, &edges);
        let result = leiden_clustering(&graph, 1.0, 10).unwrap();

        assert_eq!(result.n_clusters, 2, "expected 2 clusters for two large cliques, got {}", result.n_clusters);

        // All of clique A should be in the same cluster
        for i in 1..20 {
            assert_eq!(result.assignments[0], result.assignments[i]);
        }
        // All of clique B should be in the same cluster
        for i in 21..40 {
            assert_eq!(result.assignments[20], result.assignments[i]);
        }
        // Different clusters
        assert_ne!(result.assignments[0], result.assignments[20]);
    }
}
