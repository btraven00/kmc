use super::KmerCounterCore;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;
use std::sync::Mutex;

// Generic HashMap-based counter that accepts any hasher
pub struct GenericHashMapCounter<H: BuildHasher + Default> {
    counts: Arc<Mutex<HashMap<Vec<u8>, u64, H>>>,
}

impl<H: BuildHasher + Default + Send + Sync> KmerCounterCore for GenericHashMapCounter<H> {
    fn new() -> Self {
        GenericHashMapCounter {
            counts: Arc::new(Mutex::new(HashMap::with_hasher(H::default()))),
        }
    }

    fn merge(&self, local_counts: std::collections::HashMap<Vec<u8>, u64>) {
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
