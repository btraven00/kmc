use super::CounterCore;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;
use std::sync::Mutex;

// Generic HashMap-based counter that works with hash values (u64) directly
#[derive(Clone)]
pub struct GenericHashMapHashCounter<H: BuildHasher + Default> {
    counts: Arc<Mutex<HashMap<u64, u64, H>>>,
}

impl<H: BuildHasher + Default + Send + Sync> CounterCore for GenericHashMapHashCounter<H> {
    fn new() -> Self {
        GenericHashMapHashCounter {
            counts: Arc::new(Mutex::new(HashMap::with_hasher(H::default()))),
        }
    }

    fn merge_hashes(&self, local_counts: HashMap<u64, u64>) {
        let mut counts = self.counts.lock().unwrap();
        for (hash, count) in local_counts {
            *counts.entry(hash).or_insert(0) += count;
        }
    }

    fn unique_count(&self) -> usize {
        self.counts.lock().unwrap().len()
    }

    fn total_count(&self) -> u64 {
        self.counts.lock().unwrap().values().sum()
    }
}

impl<H: BuildHasher + Default + Send + Sync> GenericHashMapHashCounter<H> {
    /// Returns the approximate heap memory used by the HashMap data structure in bytes.
    /// This includes the memory for keys (u64 hashes), values (counts), and hash table overhead.
    pub fn heap_size_bytes(&self) -> usize {
        let counts = self.counts.lock().unwrap();
        
        // HashMap capacity-based memory: entries + internal overhead
        let capacity = counts.capacity();
        let entry_size = std::mem::size_of::<(u64, u64)>();
        
        capacity * entry_size
    }
}