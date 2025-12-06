use super::KmerCounterCore;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;

// HashMap-based counter with mutex (simple thread-safe approach)
pub struct HashMapCounter {
    counts: Arc<Mutex<HashMap<Vec<u8>, u64>>>,
}

impl KmerCounterCore for HashMapCounter {
    fn new() -> Self {
        HashMapCounter {
            counts: Arc::new(Mutex::new(HashMap::new())),
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
