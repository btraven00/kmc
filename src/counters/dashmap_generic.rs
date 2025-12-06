use super::KmerCounterCore;
use dashmap::DashMap;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;

// Generic DashMap-based counter that accepts any hasher
pub struct GenericDashMapCounter<H: BuildHasher + Clone> {
    counts: Arc<DashMap<Vec<u8>, u64, H>>,
}

impl<H: BuildHasher + Clone + Default + Send + Sync> KmerCounterCore for GenericDashMapCounter<H> {
    fn new() -> Self {
        GenericDashMapCounter {
            counts: Arc::new(DashMap::with_hasher(H::default())),
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
