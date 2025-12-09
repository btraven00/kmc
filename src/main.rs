mod counters;

use clap::Parser;
use counters::{
    DefaultHasher, FxBuildHasher, GenericDashMapCounter, GenericHashMapCounter,
    GenericHashMapHashCounter, HashCounter, HyperLogLogCounter, KmerCounter, 
    MemoryTrackedCounter, MemoryTrackedHashCounter,
};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::time::Instant;

fn format_bytes(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    
    if bytes == 0 {
        return "0 B".to_string();
    }
    
    let base = 1024_f64;
    let exp = (bytes as f64).log(base).floor() as usize;
    let exp = exp.min(UNITS.len() - 1);
    
    let value = bytes as f64 / base.powi(exp as i32);
    
    if exp == 0 {
        format!("{} {}", bytes, UNITS[exp])
    } else {
        format!("{:.2} {}", value, UNITS[exp])
    }
}

#[derive(Debug, Clone, Copy)]
enum CounterType {
    HashMap,
    DashMap,
    HyperLogLog,
}

impl std::str::FromStr for CounterType {
    type Err = String;
    
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "hashmap" => Ok(CounterType::HashMap),
            "dashmap" => Ok(CounterType::DashMap),
            "hyperloglog" | "hll" => Ok(CounterType::HyperLogLog),
            _ => Err(format!("Unknown counter type: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum HashType {
    FxHash,   // Fast non-cryptographic hash
    NtHash,   // DNA k-mer specialized rolling hash
}

impl std::str::FromStr for HashType {
    type Err = String;
    
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "fxhash" | "fx" => Ok(HashType::FxHash),
            "nthash" | "nt" => Ok(HashType::NtHash),
            _ => Err(format!("Unknown hash type: {} (available: fxhash/fx, nthash/nt)", s)),
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "kmc")]
#[command(about = "K-mer counter for FASTA files", long_about = None)]
struct Args {
    /// K-mer size
    #[arg(short = 'k', long, default_value_t = 31)]
    kmer_size: usize,

    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Counter implementation to use [default: hashmap] [possible: hashmap, dashmap, hyperloglog/hll]
    #[arg(short = 'c', long, default_value = "hashmap")]
    counter: CounterType,

    /// Hash function to use [default: fxhash] [possible: fxhash/fx, nthash/nt]
    #[arg(short = 'H', long, default_value = "fxhash")]
    hash: HashType,

    /// HyperLogLog precision (p): number of bits for register indexing (4-16)
    /// Higher values = more accurate but more memory. m=2^p registers, memory=2^p bytes
    /// Standard error ≈ 1.04/sqrt(2^p). Default p=14 gives ~0.81% error with 16KB memory
    #[arg(long, default_value_t = 14)]
    hll_precision: usize,

    /// Pin threads to CPU cores for better cache locality and NUMA performance
    #[arg(long)]
    pin_threads: bool,

    /// Output results as JSON
    #[arg(long)]
    json: bool,

    /// Input FASTA/FASTQ file
    filename: String,
}

#[derive(Clone)]
struct Chunk {
    start: u64,
    end: u64,
}

#[derive(Serialize)]
struct Config {
    k: usize,
    threads: usize,
    counter: String,
    filename: String,
}

fn round_to_3_decimals(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}

#[derive(Serialize)]
struct Results {
    unique_kmers: usize,
    total_kmers: u64,
    #[serde(serialize_with = "serialize_f64_3dp")]
    complexity: f64,
    setup_time_ms: f64,
    count_time_ms: f64,
    total_time_ms: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    peak_memory_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    peak_memory_human: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    counter_heap_bytes: Option<usize>,
}

fn serialize_f64_3dp<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_f64(round_to_3_decimals(*value))
}

#[derive(Serialize)]
struct OutputData {
    config: Config,
    results: Results,
}

fn find_next_record_boundary(file: &mut File, pos: u64) -> std::io::Result<u64> {
    file.seek(SeekFrom::Start(pos))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    
    // Read until we find a line starting with '>' (FASTA) or '@' (FASTQ)
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // End of file
            return Ok(pos);
        }
        
        if line.starts_with('>') || line.starts_with('@') {
            // Found a record boundary
            let current_pos = reader.stream_position()?;
            return Ok(current_pos - bytes_read as u64);
        }
    }
}

// New gzip processing with hash streaming
fn process_chunk_gzip_hashes(
    filename: &str,
    k: usize,
    show_progress: bool,
) -> std::io::Result<HashMap<u64, u64>> {
    use needletail::parse_fastx_file;
    
    let mut hash_counts = HashMap::new();
    let mut reader = parse_fastx_file(filename)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    
    let mut record_count = 0;
    let mut total_bases = 0;
    
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        let seq = seqrec.seq();
        total_bases += seq.len() as u64;
        count_hashes_in_sequence(&seq, k, &mut hash_counts);
        
        record_count += 1;
        if show_progress && record_count % 1000 == 0 {
            eprintln!("  Processed {} sequences, {} bases...", record_count, total_bases);
        }
    }
    
    if show_progress {
        eprintln!("  Final: {} sequences, {} bases", record_count, total_bases);
    }
    
    Ok(hash_counts)
}

fn process_chunk_gzip(
    filename: &str,
    k: usize,
    show_progress: bool,
) -> std::io::Result<HashMap<Vec<u8>, u64>> {
    use needletail::parse_fastx_file;
    
    let mut kmer_counts = HashMap::new();
    let mut reader = parse_fastx_file(filename)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    
    if show_progress {
        eprintln!("  Starting to process gzipped file...");
    }
    
    let mut record_count = 0;
    let mut total_bases = 0u64;
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        let seq = seqrec.seq();
        total_bases += seq.len() as u64;
        count_kmers_in_sequence(&seq, k, &mut kmer_counts);
        
        record_count += 1;
        if show_progress && record_count % 1000 == 0 {
            eprintln!("  Processed {} sequences, {} bases...", record_count, total_bases);
        }
    }
    
    if show_progress {
        eprintln!("  Completed: {} sequences, {} bases total", record_count, total_bases);
    }
    
    Ok(kmer_counts)
}

// New hash-based processing for NtHash optimization
fn process_chunk_hashes(
    filename: &str,
    chunk: Chunk,
    k: usize,
) -> std::io::Result<HashMap<u64, u64>> {
    let mut file = File::open(filename)?;
    file.seek(SeekFrom::Start(chunk.start))?;
    
    let mut reader = BufReader::new(file);
    let mut hash_counts = HashMap::new();
    let mut current_pos = chunk.start;
    let mut line = String::new();
    let mut current_seq = Vec::new();
    
    while current_pos < chunk.end {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }
        
        current_pos += bytes_read as u64;
        
        // Check if this is a header line
        if line.starts_with('>') || line.starts_with('@') {
            // Process previous sequence if any
            if !current_seq.is_empty() {
                count_hashes_in_sequence(&current_seq, k, &mut hash_counts);
                current_seq.clear();
            }
        } else if line.starts_with('+') {
            // FASTQ quality header - skip quality scores
            line.clear();
            let qual_bytes = reader.read_line(&mut line)?;
            current_pos += qual_bytes as u64;
        } else {
            // Sequence line - accumulate it
            let seq_line = line.trim().as_bytes();
            current_seq.extend_from_slice(seq_line);
        }
    }
    
    // Process last sequence
    if !current_seq.is_empty() {
        count_hashes_in_sequence(&current_seq, k, &mut hash_counts);
    }
    
    Ok(hash_counts)
}

// New function for hash-based counter processing
fn run_with_hash_counter<C: HashCounter>(
    args: &Args,
    _num_threads: usize,
    chunks: Vec<Chunk>,
    setup_time: std::time::Duration,
    counter: C,
) -> Result<(usize, u64, std::time::Duration, std::time::Duration, Option<u64>), Box<dyn std::error::Error>> {
    let k = args.kmer_size;
    let filename = &args.filename;
    let json_mode = args.json;
    
    let count_start = std::time::Instant::now();
    
    // Check if file is gzipped
    let is_gzip = filename.ends_with(".gz");
    
    if is_gzip {
        // For gzipped files, use needletail (single-threaded)
        match process_chunk_gzip_hashes(&filename, k, !json_mode) {
            Ok(local_counts) => {
                counter.merge_hashes(local_counts);
            }
            Err(e) => eprintln!("Error processing gzip file: {}", e),
        }
    } else {
        // Process chunks in parallel for regular files
        chunks.par_iter().enumerate().for_each(|(i, chunk)| {
            if !json_mode {
                let thread_idx = rayon::current_thread_index().unwrap_or(0);
                println!("  Thread {} processing chunk {}...", thread_idx, i);
            }
            match process_chunk_hashes(&filename, Chunk { start: chunk.start, end: chunk.end }, k) {
                Ok(local_counts) => {
                    counter.merge_hashes(local_counts);
                }
                Err(e) => eprintln!("Error processing chunk {}: {}", i, e),
            }
        });
    }
    
    let count_time = count_start.elapsed();
    let total_time = setup_time + count_time;
    
    // Update peak memory after processing
    counter.update_peak_memory();
    
    let unique_kmers = counter.unique_count();
    let total_kmers = counter.total_count();
    let peak_memory = counter.peak_memory_bytes();
    
    Ok((unique_kmers, total_kmers, count_time, total_time, peak_memory))
}

fn process_chunk(
    filename: &str,
    chunk: Chunk,
    k: usize,
) -> std::io::Result<HashMap<Vec<u8>, u64>> {
    let mut file = File::open(filename)?;
    file.seek(SeekFrom::Start(chunk.start))?;
    
    let mut reader = BufReader::new(file);
    let mut kmer_counts = HashMap::new();
    let mut current_pos = chunk.start;
    let mut line = String::new();
    let mut current_seq = Vec::new();
    
    while current_pos < chunk.end {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }
        
        current_pos += bytes_read as u64;
        
        // Check if this is a header line
        if line.starts_with('>') || line.starts_with('@') {
            // Process previous sequence if any
            if !current_seq.is_empty() {
                count_kmers_in_sequence(&current_seq, k, &mut kmer_counts);
                current_seq.clear();
            }
        } else if line.starts_with('+') {
            // FASTQ quality header - skip quality scores
            line.clear();
            let qual_bytes = reader.read_line(&mut line)?;
            current_pos += qual_bytes as u64;
        } else {
            // Sequence line - accumulate it
            let seq_line = line.trim().as_bytes();
            current_seq.extend_from_slice(seq_line);
        }
    }
    
    // Process last sequence
    if !current_seq.is_empty() {
        count_kmers_in_sequence(&current_seq, k, &mut kmer_counts);
    }
    
    Ok(kmer_counts)
}

fn count_kmers_in_sequence(seq: &[u8], k: usize, kmer_counts: &mut HashMap<Vec<u8>, u64>) {
    if seq.len() < k {
        return;
    }

    // Slide through the sequence and count k-mers
    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        *kmer_counts.entry(kmer.to_vec()).or_insert(0) += 1;
    }
}

// Direct k-mer processing for HyperLogLog using NtHash rolling hash (with batching)
// Collects hashes in a local buffer and flushes in batches to reduce lock contention
fn process_kmers_streaming_nthash(seq: &[u8], k: usize, counter: &HyperLogLogCounter) {
    if seq.len() < k {
        return;
    }

    const BATCH_SIZE: usize = 10000; // Batch hashes before updating shared counter
    let mut hash_buffer = Vec::with_capacity(BATCH_SIZE);

    // Normalize sequence to uppercase and validate for standard DNA bases
    let normalized: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
    
    // Check if sequence contains only standard DNA bases (A, C, G, T, N)
    // NtHash panics on non-standard bases like Y, R, etc.
    let is_valid_dna = normalized.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'));
    
    if !is_valid_dna {
        // Skip sequences with non-standard nucleotides
        return;
    }
    
    // Use NtHash rolling hash for efficient DNA k-mer hashing
    if let Ok(nthash_iter) = nthash::NtHashIterator::new(&normalized, k) {
        for hash_value in nthash_iter {
            hash_buffer.push(hash_value);
            
            // Flush batch when buffer is full
            if hash_buffer.len() >= BATCH_SIZE {
                counter.batch_add_hashes(&hash_buffer);
                hash_buffer.clear();
            }
        }
        
        // Flush remaining hashes
        if !hash_buffer.is_empty() {
            counter.batch_add_hashes(&hash_buffer);
        }
    }
    // Note: Non-DNA sequences or NtHash failures are silently skipped
    // This is appropriate for DNA k-mer counting applications
}

// Generic k-mer processing for HyperLogLog using sliding window hashing
fn process_kmers_streaming_generic(seq: &[u8], k: usize, counter: &HyperLogLogCounter) {
    use std::hash::{Hash, Hasher};
    use std::collections::hash_map::DefaultHasher;
    
    if seq.len() < k {
        return;
    }

    const BATCH_SIZE: usize = 10000;
    let mut hash_buffer = Vec::with_capacity(BATCH_SIZE);

    // Slide through the sequence and hash each k-mer
    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        let hash_value = hasher.finish();
        
        hash_buffer.push(hash_value);
        
        if hash_buffer.len() >= BATCH_SIZE {
            counter.batch_add_hashes(&hash_buffer);
            hash_buffer.clear();
        }
    }
    
    if !hash_buffer.is_empty() {
        counter.batch_add_hashes(&hash_buffer);
    }
}

// HLL-specific file chunk processing with streaming (no intermediate HashMap)
// Processes one file chunk (for parallelization) by reading sequences and updating
// the shared HLL counter directly. Uses sequence windows to limit memory.
fn process_file_chunk_hll_streaming(
    filename: &str,
    chunk: Chunk,
    k: usize,
    counter: &HyperLogLogCounter,
    use_nthash: bool,
) -> std::io::Result<()> {
    const MAX_SEQ_WINDOW: usize = 10_000_000; // Process long sequences in 10MB windows to limit memory
    
    let mut file = File::open(filename)?;
    file.seek(SeekFrom::Start(chunk.start))?;
    
    let mut reader = BufReader::new(file);
    let mut current_pos = chunk.start;
    let mut line = String::new();
    let mut current_seq = Vec::new();
    
    while current_pos < chunk.end {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }
        
        current_pos += bytes_read as u64;
        
        // Check if this is a header line
        if line.starts_with('>') || line.starts_with('@') {
            // Process previous sequence if any
            if !current_seq.is_empty() {
                if use_nthash {
                    process_kmers_streaming_nthash(&current_seq, k, counter);
                } else {
                    process_kmers_streaming_generic(&current_seq, k, counter);
                }
                current_seq.clear();
            }
        } else if line.starts_with('+') {
            // FASTQ quality header - skip quality scores
            line.clear();
            let qual_bytes = reader.read_line(&mut line)?;
            current_pos += qual_bytes as u64;
        } else {
            // Sequence line - accumulate it
            let seq_line = line.trim().as_bytes();
            current_seq.extend_from_slice(seq_line);
            
            // Process in windows if sequence gets too large (for very long chromosomes)
            // This limits memory usage while maintaining k-mer continuity across windows
            if current_seq.len() >= MAX_SEQ_WINDOW {
                // Process what we have so far
                if use_nthash {
                    process_kmers_streaming_nthash(&current_seq, k, counter);
                } else {
                    process_kmers_streaming_generic(&current_seq, k, counter);
                }
                
                // Keep last k-1 bases for continuity (to handle k-mers spanning the window boundary)
                if current_seq.len() >= k {
                    current_seq.drain(0..current_seq.len() - (k - 1));
                } else {
                    current_seq.clear();
                }
            }
        }
    }
    
    // Process last sequence
    if !current_seq.is_empty() {
        if use_nthash {
            process_kmers_streaming_nthash(&current_seq, k, counter);
        } else {
            process_kmers_streaming_generic(&current_seq, k, counter);
        }
    }
    
    Ok(())
}

// HLL-specific gzip processing with streaming
fn process_chunk_gzip_hll_streaming(
    filename: &str,
    k: usize,
    counter: &HyperLogLogCounter,
    show_progress: bool,
    use_nthash: bool,
) -> std::io::Result<()> {
    use needletail::parse_fastx_file;
    
    let mut reader = parse_fastx_file(filename)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    
    if show_progress {
        eprintln!("  Starting to process gzipped file...");
    }
    
    let mut record_count = 0;
    let mut total_bases = 0u64;
    
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        let seq = seqrec.seq();
        total_bases += seq.len() as u64;
        if use_nthash {
            process_kmers_streaming_nthash(&seq, k, counter);
        } else {
            process_kmers_streaming_generic(&seq, k, counter);
        }
        
        record_count += 1;
        if show_progress && record_count % 1000 == 0 {
            eprintln!("  Processed {} sequences, {} bases...", record_count, total_bases);
        }
    }
    
    if show_progress {
        eprintln!("  Completed: {} sequences, {} bases total", record_count, total_bases);
    }
    
    Ok(())
}

fn count_hashes_in_sequence(seq: &[u8], k: usize, hash_counts: &mut HashMap<u64, u64>) {
    if seq.len() < k {
        return;
    }

    // Check if sequence is valid DNA (only A,C,G,T,N after uppercase conversion)
    let is_valid_dna = seq.iter().all(|&b| {
        let upper_b = b.to_ascii_uppercase();
        matches!(upper_b, b'A' | b'C' | b'G' | b'T' | b'N')
    });

    if is_valid_dna {
        // Normalize sequence to uppercase
        let normalized: Vec<u8> = seq
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .collect();
        
        // Try NtHash rolling hash for valid DNA sequences
        if let Ok(nthash_iter) = nthash::NtHashIterator::new(&normalized, k) {
            // Use rolling hash - much faster for DNA
            for hash_value in nthash_iter {
                *hash_counts.entry(hash_value).or_insert(0) += 1;
            }
            return;
        }
    }
    
    // Fallback to sliding window with FNV hash for non-DNA or NtHash failures
    use std::hash::{Hash, Hasher};
    use std::collections::hash_map::DefaultHasher;
    
    for i in 0..=seq.len() - k {
        let mut hasher = DefaultHasher::new();
        seq[i..i + k].hash(&mut hasher);
        let hash_value = hasher.finish();
        *hash_counts.entry(hash_value).or_insert(0) += 1;
    }
}

// HLL-specific runner with streaming (minimal memory usage)
fn run_with_hll_streaming(
    args: &Args,
    _num_threads: usize,
    chunks: Vec<Chunk>,
    setup_time: std::time::Duration,
    counter: HyperLogLogCounter,
    use_nthash: bool,
) -> Result<(usize, u64, std::time::Duration, std::time::Duration, Option<u64>), Box<dyn std::error::Error>> {
    
    let k = args.kmer_size;
    let filename = args.filename.clone();
    let json_mode = args.json;
    let is_gzip = args.filename.ends_with(".gz");
    let precision = counter.get_precision();
    
    // Wrap in Arc for peak memory tracking
    let peak_memory = counter.peak_memory.clone();
    
    // Initial memory reading
    if let Some(mem) = counters::get_memory_usage_bytes() {
        peak_memory.store(mem, std::sync::atomic::Ordering::Relaxed);
    }
    
    let count_start = std::time::Instant::now();
    
    if is_gzip {
        // For gzipped files, use needletail (single-threaded) - process directly into main counter
        if let Err(e) = process_chunk_gzip_hll_streaming(&filename, k, &counter, !json_mode, use_nthash) {
            eprintln!("Error processing gzip file: {}", e);
        }
    } else {
        // Thread-local approach: each thread gets its own HLL, merge at end
        use rayon::prelude::*;
        
        let local_hlls: Vec<HyperLogLogCounter> = chunks.par_iter().enumerate().map(|(i, chunk)| {
            if !json_mode {
                let thread_idx = rayon::current_thread_index().unwrap_or(0);
                println!("  Thread {} processing chunk {}...", thread_idx, i);
            }
            
            // Create thread-local HLL with same precision
            let local_hll = HyperLogLogCounter::with_precision(precision);
            
            // Process chunk into thread-local HLL
            if let Err(e) = process_file_chunk_hll_streaming(
                &filename, 
                Chunk { start: chunk.start, end: chunk.end }, 
                k, 
                &local_hll,
                use_nthash
            ) {
                eprintln!("Error processing chunk {}: {}", i, e);
            }
            
            local_hll
        }).collect();
        
        // Merge all thread-local HLLs into the main counter
        for local_hll in local_hlls {
            counter.merge_hll(&local_hll);
        }
    }
    
    let count_time = count_start.elapsed();
    let total_time = setup_time + count_time;
    
    // Final memory reading
    if let Some(mem) = counters::get_memory_usage_bytes() {
        let mut peak = peak_memory.load(std::sync::atomic::Ordering::Relaxed);
        while mem > peak {
            match peak_memory.compare_exchange_weak(
                peak,
                mem,
                std::sync::atomic::Ordering::Relaxed,
                std::sync::atomic::Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(x) => peak = x,
            }
        }
    }
    
    let unique_kmers = counter.get_unique_count();
    let total_kmers = counter.get_total_count();
    
    // Report actual counter heap size instead of process memory
    if !json_mode {
        let counter_heap_size = counter.heap_size_bytes() as u64;
        eprintln!("HLL counter heap size: {} bytes ({:.2} KB)", 
                  counter_heap_size, counter_heap_size as f64 / 1024.0);
    }
    
    let peak_memory_val = {
        let peak = peak_memory.load(std::sync::atomic::Ordering::Relaxed);
        if peak > 0 { Some(peak) } else { None }
    };
    
    Ok((unique_kmers, total_kmers, count_time, total_time, peak_memory_val))
}

fn run_with_counter<C: KmerCounter + 'static>(
    args: &Args,
    _num_threads: usize,
    chunks: Vec<Chunk>,
    setup_time: std::time::Duration,
    counter: C,
) -> Result<(usize, u64, std::time::Duration, std::time::Duration, Option<u64>), Box<dyn std::error::Error>> {
    
    // Initial memory reading
    counter.update_peak_memory();
    
    // Process chunks in parallel
    let count_start = Instant::now();
    let filename = args.filename.clone();
    let k = args.kmer_size;
    let json_mode = args.json;
    let is_gzip = args.filename.ends_with(".gz");

    if is_gzip {
        // For gzipped files, use needletail (single-threaded)
        match process_chunk_gzip(&filename, k, !json_mode) {
            Ok(local_counts) => {
                counter.merge(local_counts);
            }
            Err(e) => eprintln!("Error processing gzip file: {}", e),
        }
    } else {
        // For uncompressed files, use chunked processing
        chunks.par_iter().enumerate().for_each(|(i, chunk)| {
            if !json_mode {
                let thread_idx = rayon::current_thread_index().unwrap_or(0);
                println!("  Thread {} processing chunk {}...", thread_idx, i);
            }
            match process_chunk(&filename, Chunk { start: chunk.start, end: chunk.end }, k) {
                Ok(local_counts) => {
                    counter.merge(local_counts);
                }
                Err(e) => eprintln!("Error processing chunk {}: {}", i, e),
            }
        });
    }

    // Final memory reading
    counter.update_peak_memory();

    let count_time = count_start.elapsed();
    let total_time = setup_time + count_time;

    let unique_kmers = counter.unique_count();
    let total_kmers = counter.total_count();
    let peak_memory = counter.peak_memory_bytes();
    
    // Note: heap_size_bytes() is type-specific and can't be called here generically
    // It's reported separately for HLL in run_with_hll_streaming
    // For HashMap/DashMap, we rely on peak_memory from the OS

    Ok((unique_kmers, total_kmers, count_time, total_time, peak_memory))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set number of threads with optional core pinning
    if args.pin_threads {
        // Get available CPU cores
        let core_ids = core_affinity::get_core_ids().unwrap_or_default();
        
        if !args.json {
            println!("Core pinning enabled. Available cores: {}", core_ids.len());
        }
        
        // Build thread pool with core pinning
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .start_handler(move |thread_index| {
                if thread_index < core_ids.len() {
                    let core_id = core_ids[thread_index];
                    if !core_affinity::set_for_current(core_id) {
                        eprintln!("Warning: Failed to pin thread {} to core {:?}", thread_index, core_id);
                    }
                } else {
                    eprintln!("Warning: Thread {} has no core to pin to (available: {})", 
                              thread_index, core_ids.len());
                }
            })
            .build_global()
            .unwrap();
    } else {
        // Build thread pool without core pinning
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    let num_threads = rayon::current_num_threads();

    if !args.json {
        println!("K-mer size: {}", args.kmer_size);
        println!("Input file: {}", args.filename);
        println!("Threads: {}", num_threads);
        println!();
    }

    // Start timing
    let start = Instant::now();

    // Check if file is gzip-compressed
    let is_gzip = args.filename.ends_with(".gz");
    
    // Get file size
    let file_size = std::fs::metadata(&args.filename)?.len();
    if !args.json {
        println!("File size: {} bytes{}", file_size, if is_gzip { " (compressed)" } else { "" });
    }

    // For compressed files, we can't chunk efficiently, so use single chunk
    let chunks = if is_gzip {
        if !args.json && num_threads > 1 {
            println!("Note: Compressed files cannot be chunked; using single-threaded processing");
        }
        vec![Chunk {
            start: 0,
            end: file_size,
        }]
    } else {
        // Calculate chunk boundaries for uncompressed files
        let chunk_size = file_size / num_threads as u64;
        let mut chunks = Vec::new();

        if !args.json {
            println!("Calculating chunk boundaries...");
        }
        
        for i in 0..num_threads {
        let theoretical_start = i as u64 * chunk_size;
        let theoretical_end = if i == num_threads - 1 {
            file_size
        } else {
            (i + 1) as u64 * chunk_size
        };

        // For first chunk, start at 0; otherwise find next record boundary
        let actual_start = if i == 0 {
            0
        } else {
            let mut file = File::open(&args.filename)?;
            find_next_record_boundary(&mut file, theoretical_start)?
        };

        // For last chunk, end at file_size; otherwise find next record boundary
        let actual_end = if i == num_threads - 1 {
            file_size
        } else {
            let mut file = File::open(&args.filename)?;
            find_next_record_boundary(&mut file, theoretical_end)?
        };

        chunks.push(Chunk {
            start: actual_start,
            end: actual_end,
        });

        if !args.json {
            println!(
                "  Chunk {}: {} - {} ({} bytes)",
                i,
                actual_start,
                actual_end,
                actual_end - actual_start
            );
        }
        }
        
        chunks
    };

    let setup_time = start.elapsed();
    if !args.json {
        println!("Setup time: {:.2?}", setup_time);
        println!();
        println!("Processing chunks in parallel...");
    }

    // Select counter implementation and run
    let (unique_kmers, total_kmers, count_time, total_time, peak_memory, counter_heap_bytes) = match (args.counter, args.hash) {
        (CounterType::HashMap, HashType::FxHash) => {
            let counter = MemoryTrackedCounter::<GenericHashMapCounter<FxBuildHasher>>::new();
            let counter_clone = counter.clone();
            let (unique, total, ctime, ttime, pmem) = run_with_counter(&args, num_threads, chunks.clone(), setup_time, counter)?;
            let heap_size = counter_clone.inner().heap_size_bytes();
            if !args.json {
                eprintln!("HashMap counter heap size: {} bytes ({:.2} MB)", 
                          heap_size, heap_size as f64 / 1_048_576.0);
            }
            (unique, total, ctime, ttime, pmem, Some(heap_size))
        },
        (CounterType::HashMap, HashType::NtHash) => {
            let counter = MemoryTrackedHashCounter::<GenericHashMapHashCounter<DefaultHasher>>::new();
            let counter_clone = counter.clone();
            let (unique, total, ctime, ttime, pmem) = run_with_hash_counter(&args, num_threads, chunks.clone(), setup_time, counter)?;
            let heap_size = counter_clone.inner().heap_size_bytes();
            if !args.json {
                eprintln!("HashMap counter heap size: {} bytes ({:.2} MB)", 
                          heap_size, heap_size as f64 / 1_048_576.0);
            }
            (unique, total, ctime, ttime, pmem, Some(heap_size))
        },
        (CounterType::DashMap, HashType::FxHash) => {
            let counter = MemoryTrackedCounter::<GenericDashMapCounter<FxBuildHasher>>::new();
            let counter_clone = counter.clone();
            let (unique, total, ctime, ttime, pmem) = run_with_counter(&args, num_threads, chunks.clone(), setup_time, counter)?;
            let heap_size = counter_clone.inner().heap_size_bytes();
            if !args.json {
                eprintln!("DashMap counter heap size: {} bytes ({:.2} MB)", 
                          heap_size, heap_size as f64 / 1_048_576.0);
            }
            (unique, total, ctime, ttime, pmem, Some(heap_size))
        },
        (CounterType::DashMap, HashType::NtHash) => {
            let counter = MemoryTrackedHashCounter::<GenericHashMapHashCounter<DefaultHasher>>::new();
            let counter_clone = counter.clone();
            let (unique, total, ctime, ttime, pmem) = run_with_hash_counter(&args, num_threads, chunks.clone(), setup_time, counter)?;
            let heap_size = counter_clone.inner().heap_size_bytes();
            if !args.json {
                eprintln!("HashMap counter heap size: {} bytes ({:.2} MB)", 
                          heap_size, heap_size as f64 / 1_048_576.0);
            }
            (unique, total, ctime, ttime, pmem, Some(heap_size))
        },
        (CounterType::HyperLogLog, hash_type) => {
            // Validate precision parameter
            if args.hll_precision < 4 || args.hll_precision > 16 {
                eprintln!("Error: HyperLogLog precision must be between 4 and 16");
                std::process::exit(1);
            }
            
            let counter = HyperLogLogCounter::with_precision(args.hll_precision);
            let heap_size = counter.heap_size_bytes();
            
            // Determine if we should use NtHash
            let use_nthash = matches!(hash_type, HashType::NtHash);
            
            if !args.json {
                let std_error = 1.04 / (1_u64 << (args.hll_precision / 2)) as f64;
                let hash_name = match hash_type {
                    HashType::NtHash => "NtHash (rolling hash)",
                    HashType::FxHash => "FxHash",
                };
                eprintln!("HyperLogLog configuration:");
                eprintln!("  Precision (p): {}", args.hll_precision);
                eprintln!("  Registers (m): {}", 1 << args.hll_precision);
                eprintln!("  Memory: {} bytes ({:.2} KB)", heap_size, heap_size as f64 / 1024.0);
                eprintln!("  Expected std error: ~{:.2}%", std_error * 100.0);
                eprintln!("  Hash function: {}", hash_name);
            }
            
            let (unique, total, ctime, ttime, pmem) = run_with_hll_streaming(&args, num_threads, chunks, setup_time, counter, use_nthash)?;
            (unique, total, ctime, ttime, pmem, Some(heap_size))
        },
    };

    // Calculate complexity
    let complexity = if total_kmers > 0 {
        unique_kmers as f64 / total_kmers as f64
    } else {
        0.0
    };

    if args.json {
        // Output as JSON
        let counter_name = match args.counter {
            CounterType::HashMap => "hashmap",
            CounterType::DashMap => "dashmap",
            CounterType::HyperLogLog => "hyperloglog",
        };

        let output = OutputData {
            config: Config {
                k: args.kmer_size,
                threads: num_threads,
                counter: counter_name.to_string(),
                filename: args.filename.clone(),
            },
            results: Results {
                unique_kmers,
                total_kmers,
                complexity,
                setup_time_ms: setup_time.as_secs_f64() * 1000.0,
                count_time_ms: count_time.as_secs_f64() * 1000.0,
                total_time_ms: total_time.as_secs_f64() * 1000.0,
                peak_memory_bytes: peak_memory,
                peak_memory_human: peak_memory.map(format_bytes),
                counter_heap_bytes,
            },
        };

        println!("{}", serde_json::to_string_pretty(&output)?);
    } else {
        // Standard output
        println!();
        println!("Results:");
        println!("  Unique k-mers: {}", unique_kmers);
        println!("  Total k-mers: {}", total_kmers);
        println!("  Complexity: {:.3}", complexity);
        if let Some(peak_mem) = peak_memory {
            println!("  Peak memory: {}", format_bytes(peak_mem));
        }
        println!();
        println!("Timing:");
        println!("  Setup time: {:.2?}", setup_time);
        println!("  Count time: {:.2?}", count_time);
        println!("  Total time: {:.2?}", total_time);
    }

    Ok(())
}
