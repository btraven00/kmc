use super::HashCounterCore;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::Arc;

/// HyperLogLog counter for approximate cardinality estimation
/// Based on the original paper: "HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm"
/// by Flajolet, Fusy, Gandouet, and Meunier (2007)
/// 
/// This implementation uses 32-bit hashes as specified in the original paper
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
    
    /// Add a pre-computed 64-bit hash value to the counter (converts to 32-bit)
    /// This is the main entry point used by the streaming NtHash path
    pub fn add_hash(&self, hash: u64) {
        let hash32 = hash as u32; // Use lower 32 bits
        self.add_hash32(hash32);
    }
    
    /// Add a 32-bit hash value to the counter (as per original HyperLogLog paper)
    fn add_hash32(&self, hash: u32) {
        // Extract the first p bits for the register index
        let j = (hash & ((1 << self.p) - 1)) as usize;
        
        // Calculate rho(w): position of first 1-bit in the remaining bits
        // Shift out the p bits used for indexing
        let w = hash >> self.p;
        
        // For 32-bit hash with p bits used for index, we have (32-p) bits remaining
        // We need to count leading zeros in just those (32-p) bits
        let rho = if w == 0 {
            // All remaining bits are 0
            (32 - self.p + 1) as u8
        } else {
            // Count leading zeros, but adjust because w only has (32-p) bits meaningful
            // w.leading_zeros() counts zeros in full 32-bit representation
            // We need to subtract the p high-order bits that don't exist in w
            let lz = w.leading_zeros();
            // Since w was shifted right by p, the leading zeros include p extra bits
            (lz - self.p as u32 + 1).min(32) as u8
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
            let hash32 = hash as u32;
            let j = (hash32 & ((1 << self.p) - 1)) as usize;
            let w = hash32 >> self.p;
            
            let rho = if w == 0 {
                (32 - self.p + 1) as u8
            } else {
                let lz = w.leading_zeros();
                (lz - self.p as u32 + 1).min(32) as u8
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
        
        // Apply bias correction for different ranges (from original paper)
        // Using 32-bit hash space as per the original HyperLogLog paper
        let estimate = if raw_estimate <= 2.5 * (self.m as f64) {
            // Small range correction
            if zero_count > 0 {
                // Linear counting for small cardinalities
                (self.m as f64) * ((self.m as f64) / (zero_count as f64)).ln()
            } else {
                raw_estimate
            }
        } else if raw_estimate <= (1.0 / 30.0) * (1u64 << 32) as f64 {
            // Intermediate range: no correction
            raw_estimate
        } else {
            // Large range correction (for 32-bit hash space)
            -((1u64 << 32) as f64) * (1.0 - raw_estimate / ((1u64 << 32) as f64)).ln()
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

impl HashCounterCore for HyperLogLogCounter {
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
        
        // Test with a known hash value
        // Using a 32-bit hash: the first 14 bits are for register index,
        // remaining 18 bits are for rho calculation
        
        // Hash with all zeros: index=0, w=0, rho should be 19 (18 zeros + 1)
        let hash1 = 0x00000000u32;
        hll.add_hash32(hash1);
        assert_eq!(hll.registers[0].load(Ordering::Relaxed), 19);
        
        // Hash with immediate 1-bit in MSB of w: index=1, w has MSB set
        // w = 0x0003_FFFF >> 0 = 0x0003_FFFF, but we need MSB of remaining 18 bits set
        // 32-bit hash structure: [18 bits for w][14 bits for index]
        // We want index=1, w=0x20000 (bit 17 set, which is MSB of 18-bit space)
        let hash2 = 0x00004001u32; // index=1, w=1 (after shift by 14)
        hll.add_hash32(hash2);
        // w = hash >> 14 = 0x00004001 >> 14 = 1
        // leading_zeros(1) = 31, rho = 31 - 14 + 1 = 18
        assert_eq!(hll.registers[1].load(Ordering::Relaxed), 18);
        
        // Hash where w has bit 17 set (MSB of 18-bit range)
        // We want w = 0x20000 after shifting by 14
        let hash3 = (0x20000u32 << 14) | 2; // index=2, w=0x20000
        hll.add_hash32(hash3);
        // w = 0x20000, leading_zeros(0x20000) = 14, rho = 14 - 14 + 1 = 1
        assert_eq!(hll.registers[2].load(Ordering::Relaxed), 1);
    }
}
