use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

pub mod hashmap_counter;
pub mod dashmap_counter;
pub mod hyperloglog_counter;
pub mod hashmap_fxhash;
pub mod hashmap_ahash;
pub mod hashmap_generic;
pub mod dashmap_generic;
pub mod hashmap_hash_generic;

pub use hashmap_counter::HashMapCounter;
pub use dashmap_counter::DashMapCounter;
pub use hyperloglog_counter::HyperLogLogCounter;
pub use hashmap_fxhash::HashMapFxHashCounter;
pub use hashmap_ahash::HashMapAHashCounter;
pub use hashmap_generic::GenericHashMapCounter;
pub use dashmap_generic::GenericDashMapCounter;
pub use hashmap_hash_generic::GenericHashMapHashCounter;

// Re-export hashers for convenience
pub use std::collections::hash_map::RandomState as DefaultHasher;
pub use rustc_hash::FxHasher;
pub use ahash::RandomState as AHasher;
pub use twox_hash::xxhash64::RandomState as XxHasher;

// NtHash wrapper for BuildHasher trait
use std::hash::{BuildHasher, Hasher};
use nthash::NtHashIterator;

#[derive(Clone)]
pub struct NtHashBuilder;

impl BuildHasher for NtHashBuilder {
    type Hasher = NtHashWrapper;

    fn build_hasher(&self) -> Self::Hasher {
        NtHashWrapper::new()
    }
}

impl Default for NtHashBuilder {
    fn default() -> Self {
        NtHashBuilder
    }
}

pub struct NtHashWrapper {
    hash: u64,
}

impl NtHashWrapper {
    fn new() -> Self {
        NtHashWrapper { hash: 0 }
    }
}

impl Hasher for NtHashWrapper {
    fn write(&mut self, bytes: &[u8]) {
        // Fast path: try NtHash directly if sequence looks like DNA
        // Single pass validation and normalization
        let mut normalized = [0u8; 64]; // Stack buffer for typical k-mer sizes
        let mut is_dna = true;

        if bytes.len() <= normalized.len() {
            for (i, &b) in bytes.iter().enumerate() {
                let upper_b = b.to_ascii_uppercase();
                if matches!(upper_b, b'A' | b'C' | b'G' | b'T' | b'N') {
                    normalized[i] = upper_b;
                } else {
                    is_dna = false;
                    break;
                }
            }

            if is_dna {
                // Try NtHash with stack-allocated buffer (no heap allocation)
                if let Ok(mut nt_iter) = NtHashIterator::new(&normalized[..bytes.len()], bytes.len()) {
                    if let Some(hash_val) = nt_iter.next() {
                        self.hash = hash_val;
                        return;
                    }
                }
            }
        }

        // Fallback: FNV-1a hash for non-DNA sequences, long sequences, or nthash failures
        let mut hash = 14695981039346656037u64; // FNV offset basis
        for &byte in bytes {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(1099511628211u64); // FNV prime
        }
        self.hash = hash;
    }

    fn finish(&self) -> u64 {
        self.hash
    }
}

pub use NtHashBuilder as NtHasher;

#[cfg(target_os = "linux")]
use procfs::process::Process;

fn get_memory_usage_bytes() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        Process::myself()
            .ok()
            .and_then(|p| p.status().ok())
            .map(|s| s.vmrss.unwrap_or(0) * 1024) // VmRSS is in KB
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

// Core trait for k-mer counting strategies (without memory tracking)
pub trait KmerCounterCore: Send + Sync {
    fn new() -> Self where Self: Sized;
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>);
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
}

// Hash-based counter trait for streaming approach
pub trait HashCounterCore: Send + Sync {
    fn new() -> Self where Self: Sized;
    fn merge_hashes(&self, local_counts: HashMap<u64, u64>);
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
}

// Wrapper that adds memory tracking to any counter implementation
pub struct MemoryTrackedCounter<C: KmerCounterCore> {
    inner: C,
    peak_memory: Arc<AtomicU64>,
}

impl<C: KmerCounterCore> MemoryTrackedCounter<C> {
    pub fn new() -> Self {
        let counter = MemoryTrackedCounter {
            inner: C::new(),
            peak_memory: Arc::new(AtomicU64::new(0)),
        };
        counter.do_update_peak_memory();
        counter
    }

    fn do_update_peak_memory(&self) {
        if let Some(current_mem) = get_memory_usage_bytes() {
            let mut peak = self.peak_memory.load(Ordering::Relaxed);
            while current_mem > peak {
                match self.peak_memory.compare_exchange_weak(
                    peak,
                    current_mem,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(x) => peak = x,
                }
            }
        }
    }
}

// Hash-based memory tracked counter
pub struct MemoryTrackedHashCounter<C: HashCounterCore> {
    inner: C,
    peak_memory: Arc<AtomicU64>,
}

impl<C: HashCounterCore> MemoryTrackedHashCounter<C> {
    pub fn new() -> Self {
        let counter = MemoryTrackedHashCounter {
            inner: C::new(),
            peak_memory: Arc::new(AtomicU64::new(0)),
        };
        counter.do_update_peak_memory();
        counter
    }

    fn do_update_peak_memory(&self) {
        if let Some(current_mem) = get_memory_usage_bytes() {
            let mut peak = self.peak_memory.load(Ordering::Relaxed);
            while current_mem > peak {
                match self.peak_memory.compare_exchange_weak(
                    peak,
                    current_mem,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(x) => peak = x,
                }
            }
        }
    }
}

// Public trait that includes memory tracking
pub trait KmerCounter: Send + Sync {
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>);
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
    fn update_peak_memory(&self);
    fn peak_memory_bytes(&self) -> Option<u64>;
}

// Hash-based counter trait with memory tracking
pub trait HashCounter: Send + Sync {
    fn merge_hashes(&self, local_counts: HashMap<u64, u64>);
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
    fn update_peak_memory(&self);
    fn peak_memory_bytes(&self) -> Option<u64>;
}

// Implement KmerCounter for the wrapper
impl<C: KmerCounterCore> KmerCounter for MemoryTrackedCounter<C> {
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>) {
        self.inner.merge(local_counts);
        self.do_update_peak_memory();
    }

    fn unique_count(&self) -> usize {
        self.inner.unique_count()
    }

    fn total_count(&self) -> u64 {
        self.inner.total_count()
    }

    fn update_peak_memory(&self) {
        self.do_update_peak_memory();
    }

    fn peak_memory_bytes(&self) -> Option<u64> {
        let peak = self.peak_memory.load(Ordering::Relaxed);
        if peak > 0 {
            Some(peak)
        } else {
            None
        }
    }
}

// Implement HashCounter for the hash-based wrapper
impl<C: HashCounterCore> HashCounter for MemoryTrackedHashCounter<C> {
    fn merge_hashes(&self, local_counts: HashMap<u64, u64>) {
        self.inner.merge_hashes(local_counts);
        self.do_update_peak_memory();
    }

    fn unique_count(&self) -> usize {
        self.inner.unique_count()
    }

    fn total_count(&self) -> u64 {
        self.inner.total_count()
    }

    fn update_peak_memory(&self) {
        self.do_update_peak_memory();
    }

    fn peak_memory_bytes(&self) -> Option<u64> {
        let peak = self.peak_memory.load(Ordering::Relaxed);
        if peak > 0 {
            Some(peak)
        } else {
            None
        }
    }
}