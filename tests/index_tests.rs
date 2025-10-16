use assert_cmd::Command;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

// Create one of several test fastas with sequences long enough for k=31
fn create_test_fasta(path: &Path, variant: usize) {
    let fasta_content = match variant {
        1 => {
            ">seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>seq2\nCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA\n"
        }
        2 => {
            ">seq1\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n>seq2\nGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT\n"
        }
        _ => {
            ">seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>seq2\nGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n"
        }
    };
    fs::write(path, fasta_content).unwrap();
}

// Index builder helper
fn build_index(fasta_path: &Path, bin_path: &Path) {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg(fasta_path)
        .arg("-o")
        .arg(bin_path)
        .assert()
        .success();

    // Check file exists and isn't empty
    assert!(
        bin_path.exists(),
        "Index file wasn't created at {:?}",
        bin_path
    );
    assert!(
        fs::metadata(bin_path).unwrap().len() > 0,
        "Index file is empty"
    );
}

#[test]
fn test_index_build() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");

    create_test_fasta(&fasta_path, 1);

    // Build index and save to file using -o
    build_index(&fasta_path, &bin_path);
}

#[test]
fn test_index_build_with_custom_kmer_window() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");

    create_test_fasta(&fasta_path, 1);

    // Build index with custom k-mer length and window size using -o
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg(fasta_path)
        .arg("-k")
        .arg("15")
        .arg("-w")
        .arg("11")
        .arg("-o")
        .arg(&bin_path)
        .assert()
        .success();

    // Check file exists and isn't empty
    assert!(bin_path.exists());
    assert!(fs::metadata(&bin_path).unwrap().len() > 0);
}

#[test]
fn test_index_union() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let combined_path = temp_dir.path().join("combined.bin");

    // Create different test FASTA files
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build indexes
    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    // Combine indexes
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("union")
        .arg("-o")
        .arg(&combined_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .assert()
        .success();

    // Check combined file exists
    assert!(combined_path.exists());

    // The combined size, should be larger than either index
    let combined_size = fs::metadata(&combined_path).unwrap().len();
    let bin1_size = fs::metadata(&bin1_path).unwrap().len();
    let bin2_size = fs::metadata(&bin2_path).unwrap().len();

    let max_individual_size = std::cmp::max(bin1_size, bin2_size);
    assert!(
        combined_size >= max_individual_size,
        "Combined index size {} should be at least as large as the largest individual index size {}",
        combined_size,
        max_individual_size
    );
}

#[test]
fn test_index_diff() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let result_path = temp_dir.path().join("result.bin");

    // Create test FASTAs
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build indexes
    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    // Diff second index from first
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("diff")
        .arg("-o")
        .arg(&result_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .assert()
        .success();

    // Check diffed file exists
    assert!(result_path.exists());

    // The diffed should be smaller than or equal to the first index
    let result_size = fs::metadata(&result_path).unwrap().len();
    let bin1_size = fs::metadata(&bin1_path).unwrap().len();

    assert!(
        result_size <= bin1_size,
        "Result index size {} should be less than or equal to the first index size {}",
        result_size,
        bin1_size
    );
}

#[test]
fn test_index_diff_three_methods() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let result_index_path = temp_dir.path().join("result_index.bin");
    let result_fastx_path = temp_dir.path().join("result_fastx.bin");
    let result_stdin_path = temp_dir.path().join("result_stdin.bin");

    // Create test FASTAs with overlapping content
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build indexes
    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    // Method 1: Index + Index diff
    let output1 = Command::cargo_bin("deacon")
        .unwrap()
        .arg("index")
        .arg("diff")
        .arg("-o")
        .arg(&result_index_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .output()
        .unwrap();
    assert!(output1.status.success());

    // Method 2: Index + FASTX file diff (with explicit k,w)
    let output2 = Command::cargo_bin("deacon")
        .unwrap()
        .arg("index")
        .arg("diff")
        .arg("-k")
        .arg("31")
        .arg("-w")
        .arg("15")
        .arg("-o")
        .arg(&result_fastx_path)
        .arg(&bin1_path)
        .arg(&fasta2_path)
        .output()
        .unwrap();
    assert!(output2.status.success());

    // Method 3: Index + FASTX stdin diff (auto-detect k,w)
    let fasta2_content = fs::read(&fasta2_path).unwrap();
    let output3 = Command::cargo_bin("deacon")
        .unwrap()
        .arg("index")
        .arg("diff")
        .arg("-o")
        .arg(&result_stdin_path)
        .arg(&bin1_path)
        .arg("-")
        .write_stdin(fasta2_content)
        .output()
        .unwrap();
    assert!(output3.status.success());

    // All three result files should exist
    assert!(result_index_path.exists());
    assert!(result_fastx_path.exists());
    assert!(result_stdin_path.exists());

    // Parse the number of remaining minimizers from stderr output
    fn extract_remaining_count(stderr: &[u8]) -> usize {
        let stderr_str = String::from_utf8_lossy(stderr);
        for line in stderr_str.lines() {
            if line.contains("remaining") {
                // Look for pattern like "Removed X minimizers, Y remaining"
                if let Some(parts) = line.split_once("remaining") {
                    let before_remaining = parts.0;
                    // Look for the last number before "remaining"
                    for word in before_remaining.split_whitespace().rev() {
                        // Try to parse the word, removing trailing comma if present
                        let clean_word = word.trim_end_matches(',');
                        if let Ok(count) = clean_word.parse::<usize>() {
                            return count;
                        }
                    }
                }
            }
        }
        panic!(
            "Could not extract remaining minimizer count from stderr: {}",
            stderr_str
        );
    }

    let remaining1 = extract_remaining_count(&output1.stderr);
    let remaining2 = extract_remaining_count(&output2.stderr);
    let remaining3 = extract_remaining_count(&output3.stderr);

    // All three methods should produce the same number of remaining minimizers
    assert_eq!(
        remaining1, remaining2,
        "Index+Index ({}) and Index+FASTX ({}) should have same remaining count",
        remaining1, remaining2
    );
    assert_eq!(
        remaining1, remaining3,
        "Index+Index ({}) and Index+FASTX stdin ({}) should have same remaining count",
        remaining1, remaining3
    );

    // Verify all result files have the same size (they should be identical)
    let size1 = fs::metadata(&result_index_path).unwrap().len();
    let size2 = fs::metadata(&result_fastx_path).unwrap().len();
    let size3 = fs::metadata(&result_stdin_path).unwrap().len();

    assert_eq!(
        size1, size2,
        "Index+Index and Index+FASTX should produce same file size"
    );
    assert_eq!(
        size1, size3,
        "Index+Index and Index+FASTX stdin should produce same file size"
    );
}

#[test]
fn test_index_diff_auto_detect_parameters() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let result_auto_path = temp_dir.path().join("result_auto.bin");
    let result_explicit_path = temp_dir.path().join("result_explicit.bin");

    // Create test FASTAs
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build index with default parameters (k=31, w=15)
    build_index(&fasta1_path, &bin1_path);

    // Method 1: Auto-detect k,w from first index
    let output_auto = Command::cargo_bin("deacon")
        .unwrap()
        .arg("index")
        .arg("diff")
        .arg("-o")
        .arg(&result_auto_path)
        .arg(&bin1_path)
        .arg(&fasta2_path)
        .output()
        .unwrap();
    assert!(output_auto.status.success());

    // Method 2: Explicitly specify k,w (should match index defaults)
    let output_explicit = Command::cargo_bin("deacon")
        .unwrap()
        .arg("index")
        .arg("diff")
        .arg("-k")
        .arg("31")
        .arg("-w")
        .arg("15")
        .arg("-o")
        .arg(&result_explicit_path)
        .arg(&bin1_path)
        .arg(&fasta2_path)
        .output()
        .unwrap();
    assert!(output_explicit.status.success());

    // Both should produce identical results
    let auto_content = fs::read(&result_auto_path).unwrap();
    let explicit_content = fs::read(&result_explicit_path).unwrap();

    assert_eq!(
        auto_content, explicit_content,
        "Auto-detected and explicit parameters should produce identical results"
    );
}

#[test]
fn test_index_dump() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");
    let dump_path = temp_dir.path().join("dump.fa");

    // All As
    let test_sequence = ">test\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    fs::write(&fasta_path, test_sequence).unwrap();

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg(&fasta_path)
        .arg("-o")
        .arg(&bin_path)
        .assert()
        .success();

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("dump")
        .arg(&bin_path)
        .arg("-o")
        .arg(&dump_path)
        .assert()
        .success();

    // Check dump contains exactly one minimizer
    let dump_content = fs::read_to_string(&dump_path).unwrap();
    let lines: Vec<&str> = dump_content.trim().lines().collect();

    assert_eq!(lines.len(), 2, "Should have one record");
    assert_eq!(lines[0], ">1", "Header should be '>1'");
    assert_eq!(lines[1], "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "Just 31 As");
}

#[test]
fn test_index_intersect() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let intersect_path = temp_dir.path().join("intersect.bin");

    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("intersect")
        .arg("-o")
        .arg(&intersect_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .assert()
        .success();

    assert!(intersect_path.exists());

    // The intersection should smaller or equal to either index size
    let intersect_size = fs::metadata(&intersect_path).unwrap().len();
    let bin1_size = fs::metadata(&bin1_path).unwrap().len();
    let bin2_size = fs::metadata(&bin2_path).unwrap().len();

    assert!(
        intersect_size <= bin1_size,
        "Intersection size {} should be <= first index size {}",
        intersect_size,
        bin1_size
    );
    assert!(
        intersect_size <= bin2_size,
        "Intersection size {} should be <= second index size {}",
        intersect_size,
        bin2_size
    );
}
