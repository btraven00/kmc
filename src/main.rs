mod counters;

use clap::Parser;
use counters::{
    DashMapCounter, HashMapCounter, HashMapFxHashCounter, HashMapAHashCounter,
    HyperLogLogCounter, KmerCounter, MemoryTrackedCounter,
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
    Default,  // SipHash (secure but slow)
    FxHash,   // Fast non-cryptographic hash
    AHash,    // Very fast, high-quality hash
    NtHash,   // DNA k-mer specialized rolling hash
}

impl std::str::FromStr for HashType {
    type Err = String;
    
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "default" | "siphash" => Ok(HashType::Default),
            "fxhash" | "fx" => Ok(HashType::FxHash),
            "ahash" => Ok(HashType::AHash),
            "nthash" | "nt" => Ok(HashType::NtHash),
            _ => Err(format!("Unknown hash type: {}", s)),
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

    /// Counter implementation to use (hashmap, dashmap, hyperloglog/hll)
    #[arg(short = 'c', long, default_value = "hashmap")]
    counter: CounterType,

    /// Hash function to use (default/siphash, fxhash/fx, ahash, nthash/nt)
    #[arg(short = 'H', long, default_value = "fxhash")]
    hash: HashType,

    /// Output results as JSON
    #[arg(long)]
    json: bool,

    /// Input FASTA/FASTQ file
    filename: String,
}

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

    Ok((unique_kmers, total_kmers, count_time, total_time, peak_memory))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set number of threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

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
    let (unique_kmers, total_kmers, count_time, total_time, peak_memory) = match (args.counter, args.hash) {
        (CounterType::HashMap, HashType::Default) => {
            let counter = MemoryTrackedCounter::<HashMapCounter>::new();
            run_with_counter(&args, num_threads, chunks, setup_time, counter)?
        },
        (CounterType::HashMap, HashType::FxHash) => {
            let counter = MemoryTrackedCounter::<HashMapFxHashCounter>::new();
            run_with_counter(&args, num_threads, chunks, setup_time, counter)?
        },
        (CounterType::HashMap, HashType::AHash) => {
            let counter = MemoryTrackedCounter::<HashMapAHashCounter>::new();
            run_with_counter(&args, num_threads, chunks, setup_time, counter)?
        },
        (CounterType::HashMap, HashType::NtHash) => {
            eprintln!("NtHash not yet implemented");
            std::process::exit(1);
        },
        (CounterType::DashMap, _) => {
            let counter = MemoryTrackedCounter::<DashMapCounter>::new();
            run_with_counter(&args, num_threads, chunks, setup_time, counter)?
        },
        (CounterType::HyperLogLog, _) => {
            let counter = MemoryTrackedCounter::<HyperLogLogCounter>::new();
            run_with_counter(&args, num_threads, chunks, setup_time, counter)?
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
