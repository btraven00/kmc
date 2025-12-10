use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

pub mod hyperloglog_counter;
pub mod hashmap_generic;
pub mod dashmap_generic;
pub mod hashmap_hash_generic;

pub use hyperloglog_counter::HyperLogLogCounter;
pub use hashmap_generic::GenericHashMapCounter;
pub use dashmap_generic::GenericDashMapCounter;
pub use hashmap_hash_generic::GenericHashMapHashCounter;

// Re-export hashers for convenience
pub use std::collections::hash_map::RandomState as DefaultHasher;
pub use rustc_hash::FxBuildHasher;

#[cfg(target_os = "linux")]
use procfs::process::Process;

pub fn get_memory_usage_bytes() -> Option<u64> {
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

// Core trait for counter implementations (without memory tracking)
// Counters can support k-mer-based counting, hash-based counting, or both
pub trait CounterCore: Send + Sync {
    fn new() -> Self where Self: Sized;
    
    // K-mer-based interface (optional)
    fn merge(&self, _local_counts: HashMap<Vec<u8>, u64>) {
        panic!("K-mer merge not supported by this counter");
    }
    
    // Hash-based interface (optional)
    fn merge_hashes(&self, _local_counts: HashMap<u64, u64>) {
        panic!("Hash merge not supported by this counter");
    }
    
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
}

// Unified wrapper that adds memory tracking to any counter implementation
#[derive(Clone)]
pub struct MemoryTrackedCounter<C: CounterCore> {
    inner: C,
    peak_memory: Arc<AtomicU64>,
}

impl<C: CounterCore> MemoryTrackedCounter<C> {
    pub fn new() -> Self {
        let counter = MemoryTrackedCounter {
            inner: C::new(),
            peak_memory: Arc::new(AtomicU64::new(0)),
        };
        counter.update_peak_memory();
        counter
    }

    fn update_peak_memory(&self) {
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
    
    /// Get reference to the inner counter for accessing type-specific methods
    pub fn inner(&self) -> &C {
        &self.inner
    }
    
    pub fn peak_memory_bytes(&self) -> Option<u64> {
        let peak = self.peak_memory.load(Ordering::Relaxed);
        if peak > 0 { Some(peak) } else { None }
    }
}

// Unified public trait for all counters with memory tracking
pub trait Counter: Send + Sync {
    // K-mer-based counting
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>);
    
    // Hash-based counting  
    fn merge_hashes(&self, local_counts: HashMap<u64, u64>);
    
    // Common interface
    fn unique_count(&self) -> usize;
    fn total_count(&self) -> u64;
    fn update_peak_memory(&self);
    fn peak_memory_bytes(&self) -> Option<u64>;
}

// Implement unified Counter trait for the wrapper
impl<C: CounterCore> Counter for MemoryTrackedCounter<C> {
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>) {
        self.inner.merge(local_counts);
        self.update_peak_memory();
    }

    fn merge_hashes(&self, local_counts: HashMap<u64, u64>) {
        self.inner.merge_hashes(local_counts);
        self.update_peak_memory();
    }

    fn unique_count(&self) -> usize {
        self.inner.unique_count()
    }

    fn total_count(&self) -> u64 {
        self.inner.total_count()
    }

    fn update_peak_memory(&self) {
        self.update_peak_memory();
    }

    fn peak_memory_bytes(&self) -> Option<u64> {
        self.peak_memory_bytes()
    }
}

// Deprecated aliases for backward compatibility
pub trait KmerCounter: Counter {}
impl<T: Counter> KmerCounter for T {}

pub trait HashCounter: Counter {}
impl<T: Counter> HashCounter for T {}

pub type MemoryTrackedHashCounter<C> = MemoryTrackedCounter<C>;