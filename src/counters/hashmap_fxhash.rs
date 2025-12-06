use super::KmerCounterCore;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;

// HashMap-based counter with FxHash (fast, non-cryptographic)
pub struct HashMapFxHashCounter {
    counts: Arc<Mutex<FxHashMap<Vec<u8>, u64>>>,
}

impl KmerCounterCore for HashMapFxHashCounter {
    fn new() -> Self {
        HashMapFxHashCounter {
            counts: Arc::new(Mutex::new(FxHashMap::default())),
        }
    }

    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>) {
        let mut counts = self.counts.lock().unwrap();
        for (kmer, count) in local_counts {
            *counts.entry(kmer).or_insert(0) += count;
        }
    }

    fn unique_count(&self) -> usize {
        self.counts.lock().unwrap().len()
    }

    fn total_count(&self) -> u64 {
        self.counts.lock().unwrap().values().sum()
    }
}
