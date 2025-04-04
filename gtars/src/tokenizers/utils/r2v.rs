#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Algorithm {
    Cbow,
    Sg,
}

#[derive(Debug, PartialEq)]
pub enum Instance {
    Cbow {
        context_ids: Vec<u32>,
        target_id: u32,
    },
    Sg {
        center_id: u32,
        context_ids: Vec<u32>,
    },
}

///
/// Create instances for the CBOW or Skip-gram algorithm.
///
/// # Arguments:
/// - `sequence`: The input sequence of integers.
/// - `window_size`: The size of the context window.
/// - `algorithm`: The algorithm to use (CBOW or Skip-gram).
///
/// # Returns:
/// - A vector of instances, each representing a context-target pair.
///
pub fn create_instances(
    sequence: &[u32],
    window_size: usize,
    algorithm: Algorithm,
) -> Vec<Instance> {
    let mut instances = Vec::new();
    if sequence.len() < window_size * 2 + 1 {
        return instances;
    }
    for (i, &target) in sequence.iter().enumerate() {
        let start = i.saturating_sub(window_size);
        let end = (i + window_size + 1).min(sequence.len());
        let mut context = Vec::new();
        for &value in &sequence[start..end] {
            if value != target {
                context.push(value);
            }
        }
        if context.is_empty() {
            continue;
        }
        match algorithm {
            Algorithm::Cbow => instances.push(Instance::Cbow {
                context_ids: context,
                target_id: target,
            }),
            Algorithm::Sg => instances.push(Instance::Sg {
                center_id: target,
                context_ids: context,
            }),
        }
    }
    instances
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_sequence() {
        let seq = vec![];
        let res = create_instances(&seq, 2, Algorithm::Cbow);
        assert!(res.is_empty());
    }

    #[test]
    fn test_short_sequence() {
        let seq = vec![1, 2];
        let res = create_instances(&seq, 2, Algorithm::Cbow);
        assert!(res.is_empty());
    }

    #[test]
    fn test_cbow_instances() {
        let seq = vec![1, 2, 3, 4, 5];
        let res = create_instances(&seq, 1, Algorithm::Cbow);
        assert_eq!(
            res,
            vec![
                Instance::Cbow {
                    context_ids: vec![2],
                    target_id: 1
                },
                Instance::Cbow {
                    context_ids: vec![1, 3],
                    target_id: 2
                },
                Instance::Cbow {
                    context_ids: vec![2, 4],
                    target_id: 3
                },
                Instance::Cbow {
                    context_ids: vec![3, 5],
                    target_id: 4
                },
                Instance::Cbow {
                    context_ids: vec![4],
                    target_id: 5
                },
            ]
        );
    }

    #[test]
    fn test_sg_instances() {
        let seq = vec![1, 2, 3];
        let res = create_instances(&seq, 1, Algorithm::Sg);
        assert_eq!(
            res,
            vec![
                Instance::Sg {
                    center_id: 1,
                    context_ids: vec![2]
                },
                Instance::Sg {
                    center_id: 2,
                    context_ids: vec![1, 3]
                },
                Instance::Sg {
                    center_id: 3,
                    context_ids: vec![2]
                },
            ]
        );
    }

    #[test]
    fn test_cbow_larger_window() {
        let seq = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let res = create_instances(&seq, 2, Algorithm::Cbow);
        let expected = vec![
            Instance::Cbow {
                context_ids: vec![2, 3],
                target_id: 1,
            },
            Instance::Cbow {
                context_ids: vec![1, 3, 4],
                target_id: 2,
            },
            Instance::Cbow {
                context_ids: vec![1, 2, 4, 5],
                target_id: 3,
            },
            Instance::Cbow {
                context_ids: vec![2, 3, 5, 6],
                target_id: 4,
            },
            Instance::Cbow {
                context_ids: vec![3, 4, 6, 7],
                target_id: 5,
            },
            Instance::Cbow {
                context_ids: vec![4, 5, 7, 8],
                target_id: 6,
            },
            Instance::Cbow {
                context_ids: vec![5, 6, 8, 9],
                target_id: 7,
            },
            Instance::Cbow {
                context_ids: vec![6, 7, 9],
                target_id: 8,
            },
            Instance::Cbow {
                context_ids: vec![7, 8],
                target_id: 9,
            },
        ];
        assert_eq!(res, expected);
    }
}
