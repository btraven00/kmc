use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

pub mod hashmap_counter;
pub mod dashmap_counter;
pub mod hyperloglog_counter;

pub use hashmap_counter::HashMapCounter;
pub use dashmap_counter::DashMapCounter;
pub use hyperloglog_counter::HyperLogLogCounter;

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

// Public trait that includes memory tracking
pub trait KmerCounter: Send + Sync {
    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>);
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
