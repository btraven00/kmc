use super::KmerCounterCore;
use dashmap::DashMap;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;

// Generic DashMap-based counter that accepts any hasher
#[derive(Clone)]
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

impl<H: BuildHasher + Clone + Default + Send + Sync> GenericDashMapCounter<H> {
    /// Returns the approximate heap memory used by the DashMap data structure in bytes.
    /// This includes the memory for keys (k-mer sequences), values (counts), and hash table overhead.
    pub fn heap_size_bytes(&self) -> usize {
        // DashMap capacity-based memory: entries + internal overhead
        let capacity = self.counts.capacity();
        let entry_size = std::mem::size_of::<(Vec<u8>, u64)>();
        let capacity_memory = capacity * entry_size;
        
        // Actual k-mer data: each Vec<u8> has its own heap allocation
        let kmer_memory: usize = self.counts.iter()
            .map(|entry| entry.key().capacity() * std::mem::size_of::<u8>())
            .sum();
        
        capacity_memory + kmer_memory
    }
}
