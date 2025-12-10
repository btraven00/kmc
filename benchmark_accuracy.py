#!/usr/bin/env -S uv run --quiet --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "matplotlib",
#     "numpy",
# ]
# ///

"""
K-mer Counting Accuracy Benchmark
Compares: Exact vs FxHash vs NtHash vs ntCard
Dataset: Drosophila melanogaster (dm6)
"""

import subprocess
import json
import time
import sys
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration
DATASET = "../data/dm6.fa"
K = 31
THREADS = 8
OUTPUT_CSV = "benchmark_accuracy_dm6.csv"
NTCARD_BIN = "../ntCard/ntcard"

def run_kmc(counter_type, hash_type=None):
    """Run kmc with specified counter and hash type"""
    cmd = [
        "cargo", "run", "--release", "--",
        "--threads", str(THREADS),
        "--json",
        "-c", counter_type
    ]
    if hash_type:
        cmd.extend(["-H", hash_type])
    cmd.append(DATASET)
    
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return json.loads(result.stdout)

def run_ntcard():
    """Run ntCard"""
    start_time = time.time()
    result = subprocess.run(
        [NTCARD_BIN, "-k", str(K), "-t", str(THREADS), "-o", "/tmp/ntcard_dm6", DATASET],
        capture_output=True, text=True, check=True
    )
    elapsed = time.time() - start_time
    
    # Parse output from stderr (ntCard outputs to stderr, not stdout)
    # Format: k=31\tF0\t123456789
    lines = result.stderr.strip().split('\n')
    f0 = f1 = None
    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 3 and parts[0] == f"k={K}":
            if parts[1] == "F0":
                f0 = int(parts[2])
            elif parts[1] == "F1":
                f1 = int(parts[2])
    
    if f0 is None or f1 is None:
        raise ValueError(f"Failed to parse ntCard output. Got F0={f0}, F1={f1}\nOutput: {result.stderr}")
    
    return {
        "unique_kmers": f0,
        "total_kmers": f1,
        "complexity": f0 / f1 if f1 else 0,
        "time_seconds": elapsed,
        "memory_mb": 0
    }

def main():
    print("=== K-mer Counting Accuracy Benchmark ===")
    print(f"Dataset: {DATASET}")
    print(f"K-mer size: {K}")
    print(f"Threads: {THREADS}")
    print()
    
    results = []
    
    # 0. Exact counting with DashMap + NtHash (ground truth)
    print("[0/3] Running exact counting (DashMap + NtHash) for ground truth...")
    exact_data = run_kmc("dashmap", "nt")
    exact_result = {
        "method": "exact_dashmap",
        "unique_kmers": exact_data["results"]["unique_kmers"],
        "total_kmers": exact_data["results"]["total_kmers"],
        "complexity": exact_data["results"]["complexity"],
        "time_seconds": exact_data["results"]["total_time_ms"] / 1000,
        "memory_mb": exact_data["results"].get("peak_memory_bytes", 0) / 1048576
    }
    results.append(exact_result)
    print(f"  Exact unique k-mers: {exact_result['unique_kmers']:,}")
    print(f"  Exact total k-mers: {exact_result['total_kmers']:,}")
    print()
    
    # 1. ntCard
    print("[1/3] Running ntCard...")
    ntcard_result = run_ntcard()
    ntcard_result["method"] = "ntcard"
    results.append(ntcard_result)
    print(f"  ntCard unique k-mers: {ntcard_result['unique_kmers']:,}")
    print(f"  ntCard total k-mers: {ntcard_result['total_kmers']:,}")
    print()
    
    # 2. HLL + NtHash
    print("[2/3] Running HyperLogLog + NtHash...")
    nt_data = run_kmc("hll", "nt")
    nt_result = {
        "method": "hll_nthash",
        "unique_kmers": nt_data["results"]["unique_kmers"],
        "total_kmers": nt_data["results"]["total_kmers"],
        "complexity": nt_data["results"]["complexity"],
        "time_seconds": nt_data["results"]["total_time_ms"] / 1000,
        "memory_mb": nt_data["results"].get("peak_memory_bytes", 0) / 1048576
    }
    results.append(nt_result)
    print(f"  HLL+NtHash unique k-mers: {nt_result['unique_kmers']:,}")
    print(f"  HLL+NtHash total k-mers: {nt_result['total_kmers']:,}")
    print()
    
    # 3. HLL + FxHash  
    print("[3/3] Running HyperLogLog + FxHash...")
    fx_data = run_kmc("hll", "fx")
    fx_result = {
        "method": "hll_fxhash",
        "unique_kmers": fx_data["results"]["unique_kmers"],
        "total_kmers": fx_data["results"]["total_kmers"],
        "complexity": fx_data["results"]["complexity"],
        "time_seconds": fx_data["results"]["total_time_ms"] / 1000,
        "memory_mb": fx_data["results"].get("peak_memory_bytes", 0) / 1048576
    }
    results.append(fx_result)
    print(f"  HLL+FxHash unique k-mers: {fx_result['unique_kmers']:,}")
    print(f"  HLL+FxHash total k-mers: {fx_result['total_kmers']:,}")
    print()
    
    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    
    # Calculate errors and speedups using DashMap as ground truth
    ground_truth = results[0]
    exact_unique = ground_truth["unique_kmers"]
    exact_total = ground_truth["total_kmers"]
    exact_time = ground_truth["time_seconds"]
    
    df['error_unique_pct'] = ((df['unique_kmers'] - exact_unique) / exact_unique) * 100
    df['error_total_pct'] = ((df['total_kmers'] - exact_total) / exact_total) * 100
    df['speedup'] = exact_time / df['time_seconds']
    
    # Print summary
    print("=== Results Summary ===")
    print()
    print(df.to_string(index=False))
    print()
    
    print("=== Accuracy Analysis (vs Exact DashMap) ===")
    print()
    print(f"{'Method':<16} | {'Unique K-mers':>13} | {'Total K-mers':>13} | {'Unique Err %':>12} | {'Total Err %':>11} | {'Time (s)':>8} | {'Speedup':>7}")
    print("-" * 115)
    
    for idx, row in df.iterrows():
        method = row['method']
        unique = int(row['unique_kmers'])
        total = int(row['total_kmers'])
        err_unique = row['error_unique_pct']
        err_total = row['error_total_pct']
        time_s = row['time_seconds']
        speedup = row['speedup']
        
        if idx == 0:  # First row is ground truth
            print(f"{method:<16} | {unique:>13,} | {total:>13,} | {'0.00%':>12} | {'0.00%':>11} | {time_s:>8.2f} | {'1.0x':>7}")
        else:
            print(f"{method:<16} | {unique:>13,} | {total:>13,} | {err_unique:>11.2f}% | {err_total:>10.2f}% | {time_s:>8.2f} | {speedup:>6.1f}x")
    
    print()
    print(f"Results saved to: {OUTPUT_CSV}")
    
    # Visualization
    print()
    print("Generating visualization...")
    create_visualization(df, exact_unique)

def create_visualization(df, exact_count):
    """Create visualization plots"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('K-mer Counting Accuracy & Performance (Drosophila dm6, k=31)', 
                 fontsize=16, fontweight='bold')
    
    methods = df['method'].values
    colors = ['#2ecc71', '#e74c3c', '#3498db', '#f39c12']  # green, red, blue, orange
    
    # 1. Unique K-mers comparison
    ax1.bar(methods, df['unique_kmers'], color=colors, alpha=0.8, edgecolor='black')
    ax1.axhline(y=exact_count, color='black', linestyle='--', linewidth=2, label='Ground Truth')
    ax1.set_ylabel('Unique K-mers', fontsize=12, fontweight='bold')
    ax1.set_title('Unique K-mer Counts', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    for i, v in enumerate(df['unique_kmers']):
        ax1.text(i, v + 0.02*max(df['unique_kmers']), f'{int(v):,}', 
                ha='center', va='bottom', fontsize=9)
    
    # 2. Error percentage (unique k-mers)
    df_no_exact = df[df['method'] != 'exact_dashmap']
    errors = df_no_exact['error_unique_pct'].values
    error_methods = df_no_exact['method'].values
    error_colors = [colors[i+1] for i in range(len(error_methods))]
    
    bars = ax2.bar(error_methods, errors, color=error_colors, alpha=0.8, edgecolor='black')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.set_ylabel('Error (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Accuracy Error vs Exact Counting', fontsize=13, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    for i, v in enumerate(errors):
        ax2.text(i, v + (0.1 if v > 0 else -0.1), f'{v:.2f}%', 
                ha='center', va='bottom' if v > 0 else 'top', 
                fontsize=10, fontweight='bold')
    
    # 3. Runtime comparison
    ax3.bar(methods, df['time_seconds'], color=colors, alpha=0.8, edgecolor='black')
    ax3.set_ylabel('Time (seconds)', fontsize=12, fontweight='bold')
    ax3.set_title('Runtime Comparison', fontsize=13, fontweight='bold')
    ax3.set_yscale('log')
    ax3.grid(axis='y', alpha=0.3)
    for i, v in enumerate(df['time_seconds']):
        ax3.text(i, v * 1.2, f'{v:.2f}s', ha='center', va='bottom', 
                fontsize=9, fontweight='bold')
    
    # 4. Speedup vs Accuracy tradeoff
    speedups = df_no_exact['speedup'].values
    abs_errors = np.abs(df_no_exact['error_unique_pct'].values)
    
    for i, method in enumerate(error_methods):
        ax4.scatter(abs_errors[i], speedups[i], s=300, color=error_colors[i], 
                    alpha=0.7, edgecolor='black', linewidth=2, zorder=3)
        ax4.annotate(method, (abs_errors[i], speedups[i]), 
                     xytext=(10, 10), textcoords='offset points', 
                     fontsize=10, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.5', facecolor=error_colors[i], alpha=0.3))
    
    ax4.set_xlabel('Absolute Error (%)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Speedup vs Exact', fontsize=12, fontweight='bold')
    ax4.set_title('Speed vs Accuracy Tradeoff', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(left=-0.5)
    
    # Add annotation for best choice
    best_speed_idx = np.argmax(speedups)
    if best_speed_idx is not None:
        ax4.annotate('Fastest!', 
                     xy=(abs_errors[best_speed_idx], speedups[best_speed_idx]),
                     xytext=(-50, -40), textcoords='offset points',
                     fontsize=11, fontweight='bold', color='green',
                     bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                     arrowprops=dict(arrowstyle='->', color='green', lw=2))
    
    plt.tight_layout()
    output_png = 'benchmark_accuracy_dm6.png'
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to: {output_png}")
    
    # Print recommendation
    print()
    print("=== Recommendation ===")
    best_method = error_methods[best_speed_idx]
    print(f"Best method for speed: {best_method}")
    print(f"  - Error: {abs_errors[best_speed_idx]:.2f}%")
    print(f"  - Speedup: {speedups[best_speed_idx]:.1f}x faster than exact")
    print(f"  - Runtime: {df[df['method'] == best_method]['time_seconds'].values[0]:.2f}s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="K-mer counting accuracy benchmark")
    parser.add_argument("--plot", action="store_true", help="Only generate visualization from existing CSV")
    args = parser.parse_args()
    
    if args.plot:
        # Just visualize from existing CSV
        if not Path(OUTPUT_CSV).exists():
            print(f"Error: {OUTPUT_CSV} not found. Run without --plot first to generate data.")
            sys.exit(1)
        
        df = pd.read_csv(OUTPUT_CSV)
        
        # Calculate errors
        exact_unique = df.iloc[0]['unique_kmers']
        exact_total = df.iloc[0]['total_kmers']
        exact_time = df.iloc[0]['time_seconds']
        
        df['error_unique_pct'] = ((df['unique_kmers'] - exact_unique) / exact_unique) * 100
        df['error_total_pct'] = ((df['total_kmers'] - exact_total) / exact_total) * 100
        df['speedup'] = exact_time / df['time_seconds']
        
        print("Generating visualization from existing data...")
        create_visualization(df, exact_unique)
    else:
        # Run full benchmark
        main()
