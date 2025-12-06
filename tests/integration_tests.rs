use std::process::Command;

#[test]
fn test_kmer_count_non_zero() {
    let output = Command::new("cargo")
        .args(&["run", "--release", "--", "--json", "tests/data/28S.fasta"])
        .output()
        .expect("Failed to execute kmc");

    assert!(output.status.success(), "kmc should exit successfully");

    let stdout = String::from_utf8(output.stdout).expect("Invalid UTF-8 output");
    let json: serde_json::Value = serde_json::from_str(&stdout)
        .expect("Failed to parse JSON output");

    let unique_kmers = json["results"]["unique_kmers"]
        .as_u64()
        .expect("unique_kmers should be a number");
    
    let total_kmers = json["results"]["total_kmers"]
        .as_u64()
        .expect("total_kmers should be a number");

    assert!(unique_kmers > 0, "unique_kmers should be greater than 0, got {}", unique_kmers);
    assert!(total_kmers > 0, "total_kmers should be greater than 0, got {}", total_kmers);
    assert!(unique_kmers <= total_kmers, "unique_kmers should be <= total_kmers");
}

#[test]
fn test_different_hash_functions() {
    let hash_types = vec!["default", "fxhash", "ahash"];
    
    for hash in &hash_types {
        let output = Command::new("cargo")
            .args(&[
                "run", "--release", "--",
                "--json",
                "-H", hash,
                "tests/data/28S.fasta"
            ])
            .output()
            .expect(&format!("Failed to execute kmc with hash {}", hash));

        assert!(output.status.success(), "kmc should exit successfully with hash {}", hash);

        let stdout = String::from_utf8(output.stdout).expect("Invalid UTF-8 output");
        let json: serde_json::Value = serde_json::from_str(&stdout)
            .expect(&format!("Failed to parse JSON output for hash {}", hash));

        let unique_kmers = json["results"]["unique_kmers"]
            .as_u64()
            .expect("unique_kmers should be a number");
        
        assert!(unique_kmers > 0, "Hash {} should produce non-zero k-mers, got {}", hash, unique_kmers);
    }
}

#[test]
fn test_complexity_calculation() {
    let output = Command::new("cargo")
        .args(&["run", "--release", "--", "--json", "tests/data/28S.fasta"])
        .output()
        .expect("Failed to execute kmc");

    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("Invalid UTF-8 output");
    let json: serde_json::Value = serde_json::from_str(&stdout)
        .expect("Failed to parse JSON output");

    let complexity = json["results"]["complexity"]
        .as_f64()
        .expect("complexity should be a number");

    assert!(complexity > 0.0 && complexity <= 1.0, 
        "complexity should be between 0 and 1, got {}", complexity);
}
