#!/usr/bin/env python3
# /// script
# dependencies = [
#   "matplotlib",
#   "pandas",
#   "numpy",
# ]
# ///
"""
Benchmark HyperLogLog counter with different thread counts.
Plots execution time vs number of threads.

Usage: 
  uv run benchmark_threads.py                    # Use drosophila (dm6) with nthash
  uv run benchmark_threads.py --human            # Use human genome (GRCh38)
  uv run benchmark_threads.py --hash=xxhash      # Use xxhash instead of nthash
  uv run benchmark_threads.py --human --ahash    # Human genome with ahash
  
Hash options: --nthash/--nt (default), --fxhash/--fx, --ahash, --xxhash/--xx
"""

import subprocess
import json
import time
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
COUNTER = "hll"
HLL_PRECISION = 14
THREAD_COUNTS = [1, 2, 4, 8, 12, 16]
RUNS_PER_CONFIG = 3  # Number of runs to average
PIN_THREADS = False # Set to True to enable core pinning

# Parse command line arguments
USE_HUMAN = "--human" in sys.argv or "-h" in sys.argv

# Hash type selection (default: nthash)
HASH_TYPE = "nt"  # default
for arg in sys.argv[1:]:
    if arg.startswith("--hash="):
        HASH_TYPE = arg.split("=")[1]
    elif arg in ["--fxhash", "--fx"]:
        HASH_TYPE = "fx"
    elif arg in ["--ahash"]:
        HASH_TYPE = "ahash"
    elif arg in ["--xxhash", "--xx"]:
        HASH_TYPE = "xx"
    elif arg in ["--nthash", "--nt"]:
        HASH_TYPE = "nt"

# Dataset selection
if USE_HUMAN:
    DATA_FILE = "../data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    DATASET_NAME = "Human (GRCh38)"
else:
    DATA_FILE = "../data/dm6.fa"
    DATASET_NAME = "Drosophila (dm6)"

def run_benchmark(threads: int) -> dict:
    """Run kmc with specified thread count and return timing results."""
    cmd = [
        "cargo", "run", "--release", "--",
        "--threads", str(threads),
        "--json",
        "-c", COUNTER,
        "-H", HASH_TYPE,
        "--hll-precision", str(HLL_PRECISION),
    ]
    
    if PIN_THREADS:
        cmd.append("--pin-threads")
    
    cmd.append(DATA_FILE)
    
    print(f"Running with {threads} threads...", end=" ", flush=True)
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    
    # Parse JSON output
    data = json.loads(result.stdout)
    print(f"✓ {data['results']['count_time_ms']:.2f}ms")
    
    return data

def main():
    print("=" * 70)
    print("Benchmarking HyperLogLog with NtHash")
    print("=" * 70)
    print(f"Dataset: {DATASET_NAME}")
    print(f"Data file: {DATA_FILE}")
    print(f"Counter: {COUNTER}")
    print(f"Hash: {HASH_TYPE}")
    print(f"HLL Precision: {HLL_PRECISION}")
    print(f"Core pinning: {'ENABLED' if PIN_THREADS else 'DISABLED'}")
    print(f"Thread counts: {THREAD_COUNTS}")
    print(f"Runs per config: {RUNS_PER_CONFIG}")
    print("=" * 70)
    print()
    
    # Check if data file exists
    if not Path(DATA_FILE).exists():
        print(f"Error: Data file not found: {DATA_FILE}")
        return
    
    # Collect results
    results = []
    
    for threads in THREAD_COUNTS:
        thread_times = []
        
        for run in range(RUNS_PER_CONFIG):
            print(f"  Run {run + 1}/{RUNS_PER_CONFIG}: ", end="")
            data = run_benchmark(threads)
            thread_times.append(data['results']['count_time_ms'])
        
        # Calculate statistics
        avg_time = np.mean(thread_times)
        std_time = np.std(thread_times)
        min_time = np.min(thread_times)
        max_time = np.max(thread_times)
        
        results.append({
            'threads': threads,
            'avg_time_ms': avg_time,
            'std_time_ms': std_time,
            'min_time_ms': min_time,
            'max_time_ms': max_time,
        })
        
        print(f"  Average: {avg_time:.2f}ms (±{std_time:.2f}ms)")
        print()
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Calculate speedup (relative to 1 thread)
    baseline_time = df[df['threads'] == 1]['avg_time_ms'].values[0]
    df['speedup'] = baseline_time / df['avg_time_ms']
    df['efficiency'] = df['speedup'] / df['threads'] * 100
    
    # Print results table
    print("=" * 70)
    print("Results Summary")
    print("=" * 70)
    print(df.to_string(index=False, float_format=lambda x: f'{x:.2f}'))
    print()
    
    # Create plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    hash_name = {"nt": "NtHash", "fx": "FxHash", "ahash": "AHash", "xx": "XxHash"}.get(HASH_TYPE, HASH_TYPE.upper())
    fig.suptitle(f'HyperLogLog Performance with {hash_name} - {DATASET_NAME}', fontsize=16, fontweight='bold')
    
    # Plot 1: Execution Time vs Threads
    ax1 = axes[0]
    ax1.errorbar(df['threads'], df['avg_time_ms'], yerr=df['std_time_ms'], 
                 marker='o', linewidth=2, markersize=8, capsize=5)
    ax1.set_xlabel('Number of Threads', fontsize=12)
    ax1.set_ylabel('Execution Time (ms)', fontsize=12)
    ax1.set_title('Execution Time vs Thread Count', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(THREAD_COUNTS)
    
    # Plot 2: Speedup vs Threads
    ax2 = axes[1]
    ax2.plot(df['threads'], df['speedup'], marker='o', linewidth=2, markersize=8, label='Actual')
    ax2.plot(df['threads'], df['threads'], '--', color='red', alpha=0.5, label='Ideal (linear)')
    ax2.set_xlabel('Number of Threads', fontsize=12)
    ax2.set_ylabel('Speedup', fontsize=12)
    ax2.set_title('Speedup vs Thread Count', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(THREAD_COUNTS)
    ax2.legend()
    
    # Plot 3: Throughput (inversely proportional to time)
    ax3 = axes[2]
    throughput = 1000.0 / df['avg_time_ms']  # runs per second
    ax3.plot(df['threads'], throughput, marker='o', linewidth=2, markersize=8, color='purple')
    ax3.set_xlabel('Number of Threads', fontsize=12)
    ax3.set_ylabel('Throughput (runs/sec)', fontsize=12)
    ax3.set_title('Throughput vs Thread Count', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.set_xticks(THREAD_COUNTS)
    
    plt.tight_layout()
    
    # Save plot with dataset and hash-specific filename
    dataset_suffix = "human" if USE_HUMAN else "dm6"
    output_file = f'benchmark_threads_hll_{dataset_suffix}_{HASH_TYPE}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    
    # Save data to CSV
    csv_file = f'benchmark_threads_hll_{dataset_suffix}_{HASH_TYPE}.csv'
    df.to_csv(csv_file, index=False)
    print(f"Data saved to: {csv_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
