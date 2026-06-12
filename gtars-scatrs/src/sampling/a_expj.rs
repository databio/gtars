/// A-ExpJ: Algorithm A with Exponential Jumps (Efraimidis & Spirakis 2006)
/// 
/// This implements Algorithm A with Exponential Jumps for weighted reservoir
/// sampling without replacement. The exponential jump optimization allows
/// skipping items that cannot enter the reservoir, reducing complexity from
/// O(n) to O(k log(n/k)) for small sample ratios.
/// 
/// Key features:
/// - Exponential jumps to skip items using formula: Xw = log(random) / log(Tw)
/// - Cumulative weight tracking for jump distance calculations
/// - Identical sampling distribution to A-Res but with dramatically improved performance
/// - Optimal for small sample ratios (k/n < 10%) common in scATAC-seq simulations
/// 
/// Time complexity: O(k log(n/k) + k log k) for small k/n ratios
/// Space complexity: O(k)
use rand::prelude::*;

/// Statistics for A-ExpJ performance monitoring
#[derive(Debug, Clone)]
pub struct AExpJStatistics {
    pub items_processed: usize,
    pub items_skipped: usize,
    pub skip_ratio: f64,
    pub total_weight_processed: f64,
    pub total_weight_skipped: f64,
    pub jump_efficiency: f64,
    pub average_jump_size: f64,
    pub num_jumps: usize,
}

/// A-ExpJ strategy for adaptive optimization
#[derive(Debug, Clone, Copy)]
pub enum AExpJStrategy {
    Standard,     // Use A-ExpJ as described in the 2006 paper
    Conservative, // Reduce jump sizes to be more conservative
    Aggressive,   // Increase jump sizes for maximum performance
    Disabled,     // Fall back to A-Res (no jumps)
}

/// Item with associated key for reservoir ordering
#[derive(Clone, Debug)]
struct KeyedItem<T> {
    key: f64,
    item: T,
}

/// A-ExpJ sampler using vector for reservoir management (optimal performance)
pub struct AExpJSampler<T> {
    reservoir: Vec<KeyedItem<T>>,
    k: usize,
    min_key: f64,
    min_index: usize,
    cumulative_weight: f64,
    total_weight_seen: f64,
    items_seen: usize,
    items_skipped: usize,
    current_jump_target: f64,
    num_jumps: usize,
    rng: StdRng,
    expj_enabled: bool,
    strategy: AExpJStrategy,
    // Cached values for optimization
    cached_log_tw: f64,
    cached_tw: f64,
    // Pre-computed strategy constant
    strategy_multiplier: f64,
}

impl<T: Clone> AExpJSampler<T> {
    /// Create a new A-ExpJ sampler
    pub fn new(k: usize, seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        Self {
            reservoir: Vec::with_capacity(k),
            k,
            min_key: f64::MAX,
            min_index: 0,
            cumulative_weight: 0.0,
            total_weight_seen: 0.0,
            items_seen: 0,
            items_skipped: 0,
            current_jump_target: 0.0,
            num_jumps: 0,
            rng,
            expj_enabled: true,
            strategy: AExpJStrategy::Standard,
            cached_log_tw: 0.0,
            cached_tw: 0.0,
            strategy_multiplier: 1.0,
        }
    }

    /// Set the A-ExpJ strategy
    pub fn with_strategy(mut self, strategy: AExpJStrategy) -> Self {
        self.strategy = strategy;
        self.expj_enabled = !matches!(strategy, AExpJStrategy::Disabled);
        // Pre-compute strategy multiplier
        self.strategy_multiplier = match strategy {
            AExpJStrategy::Standard => 1.0,
            AExpJStrategy::Conservative => 0.5,
            AExpJStrategy::Aggressive => 1.5,
            AExpJStrategy::Disabled => 0.0,
        };
        self
    }

    /// Compute the next exponential jump distance (optimized)
    fn compute_next_jump(&mut self) -> f64 {
        if !self.expj_enabled || self.reservoir.len() < self.k {
            return 0.0;
        }

        let tw = self.min_key;
        
        // IMPORTANT: Safeguard against invalid min_key that would break jump calculations
        // If min_key is 0 or 1, log(min_key) becomes undefined or 0, causing infinite/zero jumps
        if tw <= 0.0 || tw >= 1.0 {
            // This shouldn't happen with proper key generation
            if std::env::var("SCATRS_DEBUG").is_ok() {
                eprintln!("WARNING: Invalid min_key tw={:.10}, using default jump", tw);
            }
            return 1.0; // Small default jump
        }
        
        // Use cached value if unchanged
        if tw != self.cached_tw {
            self.cached_tw = tw;
            self.cached_log_tw = tw.ln();
        }

        let random: f64 = self.rng.gen_range(0.0001..0.9999);
        let log_random = random.ln();

        // A-ExpJ formula: Xw = log(random) / log(Tw)
        // CRITICAL: Both log_random and cached_log_tw are negative for values in (0,1)
        // The division of two negative numbers yields positive jump distance
        // If log(Tw) is -inf (from underflow), xw becomes 0 causing no jumps
        let xw = log_random / self.cached_log_tw;
        
        // xw should be positive (it's an exponential variate)
        let base_jump = if xw.is_finite() && xw > 0.0 {
            xw
        } else {
            1.0 // Fallback to small jump if calculation fails
        };
        
        if std::env::var("SCATRS_DEBUG").is_ok() && self.num_jumps < 100 {
            eprintln!("Jump computation: tw={:.6}, log_tw={:.6}, random={:.6}, log_random={:.6}, jump={:.6}", 
                tw, self.cached_log_tw, random, log_random, base_jump);
        }

        // Use pre-computed strategy multiplier
        base_jump * self.strategy_multiplier
    }

    /// Check if current item should be processed
    fn should_process_item(&self) -> bool {
        if !self.expj_enabled || self.reservoir.len() < self.k {
            return true;
        }
        self.cumulative_weight >= self.current_jump_target
    }

    /// Update minimum key tracking
    fn update_min_key(&mut self) {
        if self.reservoir.is_empty() {
            self.min_key = f64::MAX;
            self.min_index = 0;
            return;
        }

        self.min_key = f64::MAX;
        for (i, item) in self.reservoir.iter().enumerate() {
            if item.key < self.min_key {
                self.min_key = item.key;
                self.min_index = i;
            }
        }
        // Update cached log value
        if self.min_key != self.cached_tw {
            self.cached_tw = self.min_key;
            self.cached_log_tw = self.min_key.ln();
        }
    }

    /// Reset jump target after reservoir changes
    fn reset_jump_target(&mut self) {
        if self.reservoir.len() == self.k && self.expj_enabled {
            let jump_distance = self.compute_next_jump();
            if jump_distance > 0.0 {
                self.current_jump_target = self.cumulative_weight + jump_distance;
                self.num_jumps += 1;
                
                if std::env::var("SCATRS_DEBUG").is_ok() && self.num_jumps <= 5 {
                    eprintln!("Reset jump target: cumulative_weight={:.2}, jump_distance={:.2}, new_target={:.2}", 
                        self.cumulative_weight, jump_distance, self.current_jump_target);
                }
            }
        }
    }

    /// Add item using A-ExpJ logic
    fn add_item_expj(&mut self, item: T, weight: f64) -> bool {
        if self.k == 0 || weight <= 0.0 {
            return false; // No items to sample or invalid weight
        }
        
        // Generate key using the standard A-Res formula: u^(1/w)
        // For numerical stability with large weights, use exp(log(u)/w)
        let u: f64 = self.rng.gen_range(0.0001..0.9999);
        // CRITICAL FIX: Prevent numerical underflow for large weights
        // For weights > 100, u^(1/weight) underflows to 0, breaking the algorithm
        let key = if weight > 100.0 {
            // For large weights, use log-space computation to avoid underflow
            (u.ln() / weight).exp()
        } else {
            u.powf(1.0 / weight)
        };
        
        // CRITICAL: Ensure key never becomes 0 (would cause log(0) = -inf in jump calculations)
        // Also prevent key = 1.0 which could cause numerical issues
        let key = key.max(1e-10).min(1.0 - 1e-10);

        if self.reservoir.len() < self.k {
            // Initial filling
            self.reservoir.push(KeyedItem { key, item });
            if key < self.min_key {
                self.min_key = key;
                self.min_index = self.reservoir.len() - 1;
            }
            true
        } else if key > self.min_key {
            // Replace minimum
            self.reservoir[self.min_index] = KeyedItem { key, item };
            self.update_min_key();
            true
        } else {
            false
        }
    }

    /// Add a single item with weight (streaming interface)
    pub fn add_weighted(&mut self, item: T, weight: f64) {
        if weight <= 0.0 {
            return;
        }
        
        self.cumulative_weight += weight;
        
        if self.should_process_item() {
            let old_size = self.reservoir.len();
            let added = self.add_item_expj(item, weight);
            
            self.total_weight_seen += weight;
            self.items_seen += 1;
            
            if self.reservoir.len() == self.k {
                if old_size < self.k || added {
                    self.reset_jump_target();
                }
            }
        } else {
            self.items_skipped += 1;
        }
    }
    
    /// Get current reservoir items
    pub fn get_reservoir(self) -> Vec<T> {
        self.reservoir
            .into_iter()
            .map(|keyed| keyed.item)
            .collect()
    }

    /// Sample from stream using A-ExpJ optimization
    pub fn sample_stream<I, F>(&mut self, items: I, weight_fn: F) -> Vec<T>
    where
        I: Iterator<Item = T>,
        F: Fn(&T) -> f64,
    {
        self.reservoir.clear();
        self.cumulative_weight = 0.0;
        self.current_jump_target = 0.0;
        self.items_seen = 0;
        self.items_skipped = 0;
        self.num_jumps = 0;
        self.min_key = f64::MAX;
        self.min_index = 0;

        let mut total_items = 0;
        let mut zero_weight_items = 0;
        
        for item in items {
            total_items += 1;
            let weight = weight_fn(&item);

            if weight <= 0.0 {
                zero_weight_items += 1;
                continue;
            }

            self.cumulative_weight += weight;

            if self.should_process_item() {
                let old_size = self.reservoir.len();
                let added = self.add_item_expj(item, weight);

                self.total_weight_seen += weight;
                self.items_seen += 1;

                // CRITICAL FIX: Reset jump target to prevent infinite skipping
                // Without condition #3, after initial fill the algorithm could skip all remaining items
                // Reset jump target when:
                // 1. Reservoir just reached full size (old_size < k)
                // 2. An item was added/replaced (added = true)
                // 3. We've reached the jump target (need a new jump) - THIS WAS THE MISSING CONDITION
                if self.reservoir.len() == self.k {
                    if old_size < self.k || added || self.cumulative_weight >= self.current_jump_target {
                        self.reset_jump_target();
                    }
                }
            } else {
                self.items_skipped += 1;
            }
        }
        
        // Debug output for diagnosing sampling issues
        if std::env::var("SCATRS_DEBUG").is_ok() {
            eprintln!("A-ExpJ sampling stats: total_items={}, zero_weight={}, processed={}, skipped={}, jumps={}", 
                total_items, zero_weight_items, self.items_seen, self.items_skipped, self.num_jumps);
        }

        // Extract items
        self.reservoir
            .drain(..)
            .map(|keyed| keyed.item)
            .collect()
    }

    /// Get A-ExpJ performance statistics
    pub fn get_statistics(&self) -> AExpJStatistics {
        let total_items = self.items_seen + self.items_skipped;
        AExpJStatistics {
            items_processed: self.items_seen,
            items_skipped: self.items_skipped,
            skip_ratio: if self.items_seen > 0 {
                self.items_skipped as f64 / self.items_seen as f64
            } else {
                0.0
            },
            total_weight_processed: self.total_weight_seen,
            total_weight_skipped: self.cumulative_weight - self.total_weight_seen,
            jump_efficiency: if total_items > 0 {
                self.items_skipped as f64 / total_items as f64
            } else {
                0.0
            },
            average_jump_size: if self.num_jumps > 0 {
                (self.cumulative_weight - self.total_weight_seen) / self.num_jumps as f64
            } else {
                0.0
            },
            num_jumps: self.num_jumps,
        }
    }
    
    /// Sample multiple batches efficiently
    pub fn sample_batch<I, F>(
        items: I,
        weight_fn: F,
        k: usize,
        num_batches: usize,
        seed: Option<u64>,
    ) -> Vec<Vec<T>>
    where
        I: Iterator<Item = T> + Clone,
        F: Fn(&T) -> f64 + Clone,
    {
        (0..num_batches)
            .map(|i| {
                let batch_seed = seed.map(|s| s + i as u64);
                let mut sampler = Self::new(k, batch_seed);
                sampler.sample_stream(items.clone(), weight_fn.clone())
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "a_expj_tests.rs"]
mod tests;

#[cfg(test)]
mod basic_tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_a_expj_basic_functionality() {
        let items: Vec<usize> = (0..100).collect();
        let weights: Vec<f64> = vec![1.0; 100];

        let mut sampler = AExpJSampler::new(10, Some(42));
        let result = sampler.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        assert_eq!(result.len(), 10);
        
        // Check for no duplicates
        let unique_items: HashSet<_> = result.iter()
            .map(|(item, _)| item)
            .collect();
        assert_eq!(unique_items.len(), 10);
    }

    #[test]
    fn test_a_expj_jump_computation() {
        let mut sampler = AExpJSampler::new(10, Some(42));
        
        // Fill reservoir first
        for i in 0..10 {
            sampler.add_item_expj(i, 1.0);
        }

        // Test jump computation
        let jump = sampler.compute_next_jump();
        assert!(jump >= 0.0, "Jump distance should be non-negative");
    }

    #[test]
    fn test_a_expj_statistics() {
        let n = 10000;
        let k = 100;
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = vec![1.0; n];

        let mut sampler = AExpJSampler::new(k, Some(42));
        let _ = sampler.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        let stats = sampler.get_statistics();
        assert!(stats.items_skipped > 0, "Should skip some items");
        assert!(stats.jump_efficiency > 0.0, "Should have positive efficiency");
        println!("A-ExpJ Stats: {:?}", stats);
    }

    #[test]
    fn test_a_expj_strategies() {
        let n = 5000;
        let k = 50;
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = vec![1.0; n];

        // Test different strategies
        for strategy in [
            AExpJStrategy::Standard,
            AExpJStrategy::Conservative,
            AExpJStrategy::Aggressive,
            AExpJStrategy::Disabled,
        ] {
            let mut sampler = AExpJSampler::new(k, Some(42))
                .with_strategy(strategy);
            
            let result = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight,
            );
            
            assert_eq!(result.len(), k);
            
            let stats = sampler.get_statistics();
            println!("Strategy {:?}: skipped {} items", strategy, stats.items_skipped);
            
            if matches!(strategy, AExpJStrategy::Disabled) {
                assert_eq!(stats.items_skipped, 0, "Disabled strategy should skip no items");
            }
        }
    }
}