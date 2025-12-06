use super::KmerCounterCore;
use std::collections::HashMap;
use std::sync::Arc;
use dashmap::DashMap;

// DashMap-based counter (thread-safe concurrent HashMap)
pub struct DashMapCounter {
    counts: Arc<DashMap<Vec<u8>, u64>>,
}

impl KmerCounterCore for DashMapCounter {
    fn new() -> Self {
        DashMapCounter {
            counts: Arc::new(DashMap::new()),
        }
    }

    fn merge(&self, local_counts: HashMap<Vec<u8>, u64>) {
        for (kmer, count) in local_counts {
            self.counts.entry(kmer).and_modify(|e| *e += count).or_insert(count);
        }
    }

    fn unique_count(&self) -> usize {
        self.counts.len()
    }

    fn total_count(&self) -> u64 {
        self.counts.iter().map(|e| *e.value()).sum()
    }
}
