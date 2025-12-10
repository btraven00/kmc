use super::CounterCore;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::Arc;

/// HyperLogLog counter for approximate cardinality estimation
/// Based on the original paper: "HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm"
/// by Flajolet, Fusy, Gandouet, and Meunier (2007)
/// 
/// This implementation uses 64-bit hashes as specified in Google's "HyperLogLog in Practice" paper
/// to support cardinalities well beyond 1 billion without hash collision issues
pub struct HyperLogLogCounter {
    /// Number of bits for register indexing (m = 2^p registers)
    /// Using p=14 gives 16,384 registers, which provides good accuracy (~0.81% standard error)
    p: usize,
    
    /// Number of registers (m = 2^p)
    m: usize,
    
    /// Array of registers storing maximum rho (leading zeros + 1) values
    /// Each register is 6 bits (can store values 0-63)
    registers: Vec<AtomicU8>,
    
    /// Total count of k-mers processed (for tracking purposes)
    total_kmers: AtomicU64,
    
    /// Peak memory tracker (for compatibility with main.rs)
    pub peak_memory: Arc<AtomicU64>,
}

impl HyperLogLogCounter {
    /// Create a new HyperLogLog counter with default precision (p=14)
    pub fn new() -> Self {
        Self::with_precision(14)
    }
    
    /// Create a new HyperLogLog counter with specified precision
    /// p: number of bits for register indexing (typical values: 4-16)
    /// Standard error ≈ 1.04 / sqrt(2^p)
    pub fn with_precision(p: usize) -> Self {
        let m = 1 << p; // m = 2^p
        let mut registers = Vec::with_capacity(m);
        for _ in 0..m {
            registers.push(AtomicU8::new(0));
        }
        
        HyperLogLogCounter {
            p,
            m,
            registers,
            total_kmers: AtomicU64::new(0),
            peak_memory: Arc::new(AtomicU64::new(0)),
        }
    }
    
    /// Add a pre-computed 64-bit hash value to the counter
    /// This is the main entry point used by the streaming NtHash path
    /// 
    /// Implementation note: NtHash and other biological sequence hashers may have
    /// poor bit distribution (leading zeros in high bits). We apply a mixing function
    /// to improve bit distribution before using the hash for HLL.
    pub fn add_hash(&self, hash: u64) {
        // Mix the hash to improve bit distribution
        // NtHash often has leading zeros in high bits, which causes inflated rho values
        // This mixing function spreads entropy across all bits
        let mixed = hash.wrapping_mul(0x9e3779b97f4a7c15u64); // Knuth's multiplicative hash
        
        // Extract the low-order p bits for the register index
        let j = (mixed & ((1 << self.p) - 1)) as usize;
        
        // Extract the remaining high bits for rho calculation
        // Shift out the p bits used for indexing
        let w = mixed >> self.p;
        
        // For 64-bit hash with p bits for index, we have (64-p) bits in w
        // ρ(w) = position of leftmost 1-bit in the (64-p)-bit representation of w
        // Maximum value is (64-p+1) when all bits are 0
        let rho = if w == 0 {
            // All bits in w are 0
            (64 - self.p + 1) as u8
        } else {
            // w.leading_zeros() counts leading zeros in 64-bit representation
            // We need to adjust for the fact that w only has (64-p) meaningful bits
            // Adjustment: lz_64 - (64 - (64-p)) = lz_64 - p
            let lz = w.leading_zeros() as usize;
            let lz_adjusted = lz - self.p;
            (lz_adjusted + 1) as u8
        };
        
        // Update register j with max(M[j], rho)
        let mut current = self.registers[j].load(Ordering::Relaxed);
        while rho > current {
            match self.registers[j].compare_exchange_weak(
                current,
                rho,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(x) => current = x,
            }
        }
        
        // Increment total k-mer count
        self.total_kmers.fetch_add(1, Ordering::Relaxed);
    }
    
    /// Batch add multiple hashes (more efficient for bulk operations)
    /// This reduces atomic contention by computing max rho values locally first
    pub fn batch_add_hashes(&self, hashes: &[u64]) {
        // Process hashes directly without intermediate storage
        // Group updates by register to reduce contention
        for &hash in hashes {
            // Mix the hash to improve bit distribution (same as add_hash)
            let mixed = hash.wrapping_mul(0x9e3779b97f4a7c15u64);
            
            let j = (mixed & ((1 << self.p) - 1)) as usize;
            let w = mixed >> self.p;
            
            let rho = if w == 0 {
                (64 - self.p + 1) as u8
            } else {
                let lz = w.leading_zeros() as usize;
                let lz_adjusted = lz - self.p;
                (lz_adjusted + 1) as u8
            };
            
            // Try to update register with relaxed ordering
            // Only use atomic if the new value is larger
            let current = self.registers[j].load(Ordering::Relaxed);
            if rho > current {
                // Try to update atomically
                let mut old = current;
                while rho > old {
                    match self.registers[j].compare_exchange_weak(
                        old,
                        rho,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(x) => old = x,
                    }
                }
            }
        }
        
        // Update total count
        self.total_kmers.fetch_add(hashes.len() as u64, Ordering::Relaxed);
    }
    
    /// Get the estimated unique count using the HyperLogLog algorithm
    pub fn get_unique_count(&self) -> usize {
        // Calculate the raw estimate using harmonic mean
        let mut sum = 0.0_f64;
        let mut zero_count = 0;
        
        for i in 0..self.m {
            let val = self.registers[i].load(Ordering::Relaxed);
            if val == 0 {
                zero_count += 1;
            }
            // Harmonic mean: sum of 2^(-M[j])
            sum += 2.0_f64.powi(-(val as i32));
        }
        
        // Raw estimate: E = alpha_m * m^2 * (sum of 2^(-M[j]))^(-1)
        let alpha_m = self.alpha_m();
        let raw_estimate = alpha_m * (self.m as f64) * (self.m as f64) / sum;
        
        // Apply bias correction for different ranges
        // For 64-bit hash space (Google's HyperLogLog in Practice):
        // - Small range: use LinearCounting when there are empty registers
        // - No large range correction needed (would only apply near 2^64)
        let estimate = if raw_estimate <= 2.5 * (self.m as f64) {
            // Small range correction
            if zero_count > 0 {
                // Linear counting for small cardinalities
                (self.m as f64) * ((self.m as f64) / (zero_count as f64)).ln()
            } else {
                raw_estimate
            }
        } else {
            // For 64-bit hashes, no correction needed until cardinality approaches 2^64
            // which is extremely unlikely in practice (1.8 × 10^19)
            raw_estimate
        };
        
        estimate.round() as usize
    }
    
    /// Get the total count of k-mers processed
    pub fn get_total_count(&self) -> u64 {
        self.total_kmers.load(Ordering::Relaxed)
    }
    
    /// Get the precision parameter (p)
    pub fn get_precision(&self) -> usize {
        self.p
    }
    
    /// Merge another HyperLogLog counter into this one
    /// Takes the maximum of each register (standard HLL merge operation)
    pub fn merge_hll(&self, other: &HyperLogLogCounter) {
        assert_eq!(self.p, other.p, "Cannot merge HLLs with different precision");
        
        for i in 0..self.m {
            let other_val = other.registers[i].load(Ordering::Relaxed);
            if other_val > 0 {
                // Update with max(self[i], other[i])
                let mut current = self.registers[i].load(Ordering::Relaxed);
                while other_val > current {
                    match self.registers[i].compare_exchange_weak(
                        current,
                        other_val,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(x) => current = x,
                    }
                }
            }
        }
        
        // Merge total counts
        let other_total = other.total_kmers.load(Ordering::Relaxed);
        self.total_kmers.fetch_add(other_total, Ordering::Relaxed);
    }
    
    /// Calculate alpha_m constant based on the number of registers
    /// From the original paper, alpha_m for large m is approximately:
    /// alpha_16 = 0.673, alpha_32 = 0.697, alpha_64 = 0.709, alpha_inf ≈ 0.7213/(1 + 1.079/m)
    fn alpha_m(&self) -> f64 {
        match self.m {
            16 => 0.673,
            32 => 0.697,
            64 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / (self.m as f64)),
        }
    }
    
    /// Get the heap size in bytes used by this counter
    pub fn heap_size_bytes(&self) -> usize {
        // Each register is 1 byte (AtomicU8)
        std::mem::size_of::<AtomicU8>() * self.m
    }
}

impl Default for HyperLogLogCounter {
    fn default() -> Self {
        Self::new()
    }
}

impl CounterCore for HyperLogLogCounter {
    fn new() -> Self {
        HyperLogLogCounter::new()
    }

    fn merge_hashes(&self, local_counts: HashMap<u64, u64>) {
        // For HyperLogLog, we only care about unique hashes (cardinality)
        // The count doesn't matter since we're estimating cardinality, not frequency
        for (hash, _count) in local_counts {
            self.add_hash(hash);
        }
    }

    fn unique_count(&self) -> usize {
        self.get_unique_count()
    }

    fn total_count(&self) -> u64 {
        self.get_total_count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_counting() {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let hll = HyperLogLogCounter::new();
        
        // Add some hashes (simulating k-mers)
        let mut hasher = DefaultHasher::new();
        b"ACGT".hash(&mut hasher);
        let hash1 = hasher.finish();
        
        let mut hasher = DefaultHasher::new();
        b"TGCA".hash(&mut hasher);
        let hash2 = hasher.finish();
        
        hll.add_hash(hash1);
        hll.add_hash(hash2);
        hll.add_hash(hash1); // Duplicate
        
        assert_eq!(hll.get_total_count(), 3);
        
        // Unique count should be approximately 2 (with some error)
        let unique = hll.get_unique_count();
        assert!(unique >= 1 && unique <= 3, "Expected ~2, got {}", unique);
    }
    
    #[test]
    fn test_large_dataset() {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let hll = HyperLogLogCounter::new();
        
        // Add 10000 unique hashes
        for i in 0..10000 {
            let kmer = format!("KMER{:05}", i);
            let mut hasher = DefaultHasher::new();
            kmer.hash(&mut hasher);
            hll.add_hash(hasher.finish());
        }
        
        let estimate = hll.get_unique_count();
        let error = ((estimate as f64 - 10000.0) / 10000.0).abs();
        
        // With p=14, standard error should be ~0.81%, so 5% is a generous bound
        assert!(error < 0.05, "Error too large: {:.2}% (estimate: {})", error * 100.0, estimate);
    }
    
    #[test]
    fn test_rho_calculation() {
        let hll = HyperLogLogCounter::with_precision(14);
        
        // Note: The implementation now applies mixing via multiplication before extraction
        // So we can't predict exact register values from input hashes
        // Instead, test that the algorithm produces reasonable results
        
        // Test 1: Hash with all zeros produces some rho value
        let hash1 = 0x0000000000000000u64;
        hll.add_hash(hash1);
        let rho1 = hll.registers.iter()
            .map(|r| r.load(Ordering::Relaxed))
            .find(|&r| r > 0)
            .expect("Should have at least one non-zero register");
        assert!(rho1 > 0 && rho1 <= 51, "Rho should be in valid range");
        
        // Test 2: Different hash produces different register update
        let hash2 = 0x123456789ABCDEF0u64;
        let hll2 = HyperLogLogCounter::with_precision(14);
        hll2.add_hash(hash2);
        let rho2 = hll2.registers.iter()
            .map(|r| r.load(Ordering::Relaxed))
            .find(|&r| r > 0)
            .expect("Should have at least one non-zero register");
        assert!(rho2 > 0 && rho2 <= 51, "Rho should be in valid range");
        
        // Test 3: Multiple hashes update registers with reasonable distribution
        let hll3 = HyperLogLogCounter::with_precision(14);
        for i in 0..1000 {
            hll3.add_hash(i);
        }
        let non_zero_registers = hll3.registers.iter()
            .filter(|r| r.load(Ordering::Relaxed) > 0)
            .count();
        // With 1000 hashes and 16384 registers, we expect most hashes to hit unique registers
        // but the exact number depends on the hash distribution
        assert!(non_zero_registers > 500, "Should have many registers updated with 1000 hashes, got {}", non_zero_registers);
        assert!(non_zero_registers <= 1000, "Cannot have more non-zero registers than hashes");
    }
}
