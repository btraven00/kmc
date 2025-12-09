use super::KmerCounterCore;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;
use std::sync::Mutex;

// Generic HashMap-based counter that accepts any hasher
#[derive(Clone)]
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

impl<H: BuildHasher + Default + Send + Sync> GenericHashMapCounter<H> {
    /// Returns the approximate heap memory used by the HashMap data structure in bytes.
    /// This includes the memory for keys (k-mer sequences), values (counts), and hash table overhead.
    pub fn heap_size_bytes(&self) -> usize {
        let counts = self.counts.lock().unwrap();
        
        // HashMap capacity-based memory: entries + internal overhead
        let capacity = counts.capacity();
        let entry_size = std::mem::size_of::<(Vec<u8>, u64)>();
        let capacity_memory = capacity * entry_size;
        
        // Actual k-mer data: each Vec<u8> has its own heap allocation
        let kmer_memory: usize = counts.keys()
            .map(|k| k.capacity() * std::mem::size_of::<u8>())
            .sum();
        
        capacity_memory + kmer_memory
    }
}
