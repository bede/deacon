#![cfg(feature = "server")]
//! Tests for the Deacon server functionality
//! Basically a duplicate of the contents of `filter_tests.rs`,
//! but swapping to use the server filtering instead of local.
use assert_cmd::Command;
use std::fs;
use std::fs::File;
use std::path::Path;
use std::process::Command as StdCommand;
use tempfile::tempdir;

/// Run a given command with a server running on a random port.
/// The server will be started in a separate thread, and killed after the command completes.
/// The server will be started with the provided index path and a random port.
///
/// # Arguments:
/// * `$index_path` - The path to the index file to be used by the server.
/// * `$cmd` - A closure that takes a string slice representing the port number.
macro_rules! run_with_server {
    ($index_path: expr, $cmd: expr) => {
        let random_port = rand::random_range(3000..10000);
        let rt = tokio::runtime::Runtime::new().unwrap();
        let server_thread = rt.spawn(async move {
            deacon::server::run_server($index_path, random_port).await;
        });
        let random_port = random_port.to_string();
        let mut retries = 0;
        // Ensure the server is actually running
        while retries < 10 {
            let client = reqwest::blocking::Client::new();

            let response = client
                .get("http://0.0.0.0:".to_owned() + &random_port + "/")
                .send();
            match response {
                Ok(response) => {
                    // Check if the response indicates a match
                    if response.status().is_success() {
                        break;
                    } else {
                        // Request worked, but returned an error, so wait and retry
                        std::thread::sleep(std::time::Duration::from_millis(100));
                    }
                }
                Err(_) => {
                    // If the request fails, wait and retry
                    std::thread::sleep(std::time::Duration::from_millis(100));
                }
            }
            retries += 1;
        }
        if retries == 10 {
            panic!("Server did not start in time");
        }
        eprintln!("Started server on port {random_port}");
        // Run the provided command after starting the server
        $cmd(&random_port);

        // Kill the server thread
        server_thread.abort();
    };
}

fn create_test_fasta(path: &Path) {
    let fasta_content = ">seq1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n>seq2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
    fs::write(path, fasta_content).unwrap();
}

fn create_test_fasta_aaa(path: &Path) {
    let fasta_content = ">seq1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    fs::write(path, fasta_content).unwrap();
}

fn create_test_fastq(path: &Path) {
    let fastq_content = "@seq1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@seq2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_paired_fastq(path1: &Path, path2: &Path) {
    let fastq_content1 = "@read1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    let fastq_content2 = "@read1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

fn build_index(fasta_path: &Path, bin_path: &Path) {
    let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
        .arg("index")
        .arg("build")
        .arg(fasta_path)
        .output()
        .expect("Failed to execute command");

    fs::write(bin_path, output.stdout).expect("Failed to write index file");
    assert!(output.status.success(), "Index build command failed");
}

fn create_test_fasta_sc2(path: &Path) {
    let fasta_content =
        ">mn908947.3_0:60\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n";
    fs::write(path, fasta_content).unwrap();
}

fn create_test_fastq_sc2_fwd(path: &Path) {
    let fastq_content = "@mn908947.3_0:60_fwd\n\
        ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_fastq_sc2_rev(path: &Path) {
    let fastq_content = "@mn908947.3_0:60_rev\n\
        AGATCTACAAGAGATCGAAAGTTGGTTGGTTTGTTACCTGGGAAGGTATAAACCTTTAAT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

fn create_test_paired_fastq_sc2_fwd(path1: &Path, path2: &Path) {
    let fastq_content1 = "@mn908947.3_0:60_fwd\n\
        ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    let fastq_content2 = "@mn908947.3_60:120_fwd\n\
        GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

fn create_test_paired_fastq_sc2_rev(path1: &Path, path2: &Path) {
    let fastq_content1 = "@mn908947.3_0:60_rev\n\
        AGATCTACAAGAGATCGAAAGTTGGTTGGTTTGTTACCTGGGAAGGTATAAACCTTTAAT\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    let fastq_content2 = "@mn908947.3_60:120_rev\n\
        AGTGCACTAAGCATGCAGCCGAGTGACAGCCACACAGATTTTAAAGTTCGTTTAGAGAAC\n\
        +\n\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(path1, fastq_content1).unwrap();
    fs::write(path2, fastq_content2).unwrap();
}

#[test]
fn test_filter_to_file() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta_aaa(&fasta_path);
    create_test_fastq(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Run filtering command
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .arg("--summary")
            .arg(&summary_path)
            .assert()
            .success();
    });

    // Check output and report creation
    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    // With new default behavior: sequences without matches are filtered out (sequences too short for k=31)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.is_empty(),
        "Output file should be empty - sequences too short for minimizers"
    );
}

#[test]
fn test_filter_to_file_gzip() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.gz");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Check gzipped output file creation
    assert!(output_path.exists(), "Gzipped output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Gzipped output file is empty"
    );
}

#[test]
fn test_filter_to_file_zstd() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.zst");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Check that zstd output file was created
    assert!(output_path.exists(), "Zstd output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Zstd output file is empty"
    );
}

#[test]
fn test_filter_to_file_xz() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.xz");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Check that xz output file was created
    assert!(output_path.exists(), "XZ output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "XZ output file is empty"
    );
}

#[test]
fn test_filter_deplete_flag() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_depleted.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with deplete flag wasn't created"
    );
}

#[test]
fn test_filter_rename() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_renamed.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--rename")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with rename flag wasn't created"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.contains("@1\n") || output_content.contains("@2\n"),
        "Output does not contain renamed sequences"
    );
}

#[test]
fn test_filter_min_matches() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_min_matches.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("2")
            .arg("--rel-threshold")
            .arg("0.01")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with min_matches parameter wasn't created"
    );
}

#[test]
fn test_filter_prefix_length() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_prefix.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--prefix-length")
            .arg("6")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with prefix_length parameter wasn't created"
    );
}

#[test]
fn test_filter_paired() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Run filtering command with paired-end reads (using -a 1 so short sequences pass through)
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Check output creation
    assert!(output_path.exists(), "Output file wasn't created");

    // Validate output content (should be interleaved)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(!output_content.is_empty(), "Output file is empty");
}

#[test]
fn test_filter_paired_with_deplete() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_depleted.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with deplete flag wasn't created"
    );
}

#[test]
fn test_filter_paired_with_rename() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_renamed.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--rename")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with rename flag wasn't created"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.contains("@1\n") && output_content.contains("@2\n"),
        "Output does not contain renamed sequences"
    );
}

#[test]
fn test_filter_paired_with_min_matches() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_min_matches.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("2")
            .arg("--rel-threshold")
            .arg("0.01")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with min_matches parameter wasn't created"
    );
}

#[test]
fn test_interleaved_paired_reads_stdin() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let interleaved_fastq_path = temp_dir.path().join("interleaved_reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    // Create test files
    create_test_fasta(&fasta_path);

    let interleaved_content =
        "@read1/1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read1/2\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            .to_owned()
            + "@read2/1\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2/2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(&interleaved_fastq_path, interleaved_content).unwrap();

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Test piping interleaved file to stdin for processing
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"));
        let output = cmd
            .arg("client")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg("-") // stdin for input
            .arg("-") // stdin for input2 (signals interleaved mode)
            .arg("--output")
            .arg(&output_path)
            .stdin(File::open(&interleaved_fastq_path).unwrap())
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success(), "Command failed");
        assert!(output_path.exists(), "Output file wasn't created");
    });

    // Validate output content (should contain processed reads)
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(!output_content.is_empty(), "Output file is empty");
}

#[test]
fn test_single_read_stdin() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta(&fasta_path);

    let fastq_content = "@read1\nACGTGCATAGCTGCATGCATGCATGCATGCATGCATGCAATGCAACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2\nTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTGCAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    fs::write(&fastq_path, fastq_content).unwrap();

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Test single-end stdin
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"));
        let output = cmd
            .arg("client")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg("-") // stdin
            .arg("--output")
            .arg(&output_path)
            .stdin(File::open(&fastq_path).unwrap())
            .output()
            .expect("Failed to execute command");

        assert!(
            output.status.success(),
            "Command failed for single-read stdin"
        );
    });
    assert!(
        output_path.exists(),
        "Output file wasn't created for single-read stdin"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        !output_content.is_empty(),
        "Output file is empty for single-read stdin"
    );
    assert!(
        output_content.contains("read1"),
        "read1 not found in output"
    );
    assert!(
        output_content.contains("read2"),
        "read2 not found in output"
    );
}

#[test]
fn test_filter_filtration_fwd() {
    // Tests filtering with forward reads from SC2
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta_sc2(&fasta_path);
    create_test_fastq_sc2_fwd(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .arg("--summary")
            .arg(&summary_path)
            .arg("--abs-threshold")
            .arg("1")
            .arg("--rel-threshold")
            .arg("0.01")
            .assert()
            .success();
    });

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_filtration_rev() {
    // Tests filtering with reverse read from SC2
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let summary_path = temp_dir.path().join("summary.json");

    create_test_fasta_sc2(&fasta_path);
    create_test_fastq_sc2_rev(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .arg("--summary")
            .arg(&summary_path)
            .assert()
            .success();
    });

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_paired_filtration_fwd() {
    // Tests that both reads are filtered when a forward read matches the SC2 ref
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta_sc2(&fasta_path);
    create_test_paired_fastq_sc2_fwd(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(output_path.exists(), "Output file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[test]
fn test_filter_paired_filtration_rev() {
    // Tests that both reads are filtered when a reverse read matches the SC2 ref
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");

    create_test_fasta_sc2(&fasta_path);
    create_test_paired_fastq_sc2_rev(&fastq_path1, &fastq_path2);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(output_path.exists(), "Output file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(output_content.is_empty(), "Output file is not empty");
}

#[cfg(test)]
mod output2_tests {
    use assert_cmd::Command;
    use std::fs;
    use std::path::Path;
    use std::process::Command as StdCommand;
    use tempfile::tempdir;

    fn create_test_fasta(path: &Path) {
        let fasta_content = ">seq1\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA\n>seq2\nCGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC\n";
        fs::write(path, fasta_content).unwrap();
    }

    fn create_test_paired_fastq(path1: &Path, path2: &Path) {
        let fastq_content1 = "@read1\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2\nCGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        let fastq_content2 = "@read1\nTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTG\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n@read2\nTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        fs::write(path1, fastq_content1).unwrap();
        fs::write(path2, fastq_content2).unwrap();
    }

    fn build_index(fasta_path: &Path, bin_path: &Path) {
        let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
            .arg("index")
            .arg("build")
            .arg(fasta_path)
            .output()
            .expect("Failed to execute command");

        fs::write(bin_path, output.stdout).expect("Failed to write index file");
        assert!(output.status.success(), "Index build command failed");
    }

    #[test]
    fn test_filter_paired_with_output2() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path1 = temp_dir.path().join("reads_1.fastq");
        let fastq_path2 = temp_dir.path().join("reads_2.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path1 = temp_dir.path().join("filtered_1.fastq.gz");
        let output_path2 = temp_dir.path().join("filtered_2.fastq.gz");

        create_test_fasta(&fasta_path);
        create_test_paired_fastq(&fastq_path1, &fastq_path2);

        build_index(&fasta_path, &bin_path);
        assert!(bin_path.exists(), "Index file wasn't created");

        // Run filtering command with separate output files
        run_with_server!(bin_path.clone(), |port: &str| {
            let mut cmd = Command::cargo_bin("deacon").unwrap();
            cmd.arg("client")
                .arg("http://0.0.0.0:".to_owned() + port)
                .arg(&fastq_path1)
                .arg(&fastq_path2)
                .arg("--output")
                .arg(&output_path1)
                .arg("--output2")
                .arg(&output_path2)
                .assert()
                .success();
        });

        // Check both output files were created
        assert!(output_path1.exists(), "First output file wasn't created");
        assert!(output_path2.exists(), "Second output file wasn't created");

        // Validate output content
        assert!(
            fs::metadata(&output_path1).unwrap().len() > 0,
            "First gzipped output file is empty"
        );
        assert!(
            fs::metadata(&output_path2).unwrap().len() > 0,
            "Second gzipped output file is empty"
        );

        // Actually decompress and check if there are reads
        use flate2::read::GzDecoder;
        use std::fs::File;
        use std::io::Read;

        let file1 = File::open(&output_path1).unwrap();
        let mut gz1 = GzDecoder::new(file1);
        let mut contents1 = String::new();
        gz1.read_to_string(&mut contents1).unwrap();

        let file2 = File::open(&output_path2).unwrap();
        let mut gz2 = GzDecoder::new(file2);
        let mut contents2 = String::new();
        gz2.read_to_string(&mut contents2).unwrap();

        println!(
            "Output2 test - Output1 length: {}, Output2 length: {}",
            contents1.len(),
            contents2.len()
        );
        println!(
            "Output2 test - Output1 preview: {:?}",
            &contents1.chars().take(100).collect::<String>()
        );
        println!(
            "Output2 test - Output2 preview: {:?}",
            &contents2.chars().take(100).collect::<String>()
        );
    }

    #[test]
    fn test_filter_paired_with_output2_gzip() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path1 = temp_dir.path().join("reads_1.fastq");
        let fastq_path2 = temp_dir.path().join("reads_2.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path1 = temp_dir.path().join("filtered_1.fastq.gz");
        let output_path2 = temp_dir.path().join("filtered_2.fastq.gz");

        create_test_fasta(&fasta_path);
        create_test_paired_fastq(&fastq_path1, &fastq_path2);
        build_index(&fasta_path, &bin_path);
        run_with_server!(bin_path.clone(), |port: &str| {
            let mut cmd = Command::cargo_bin("deacon").unwrap();
            cmd.arg("client")
                .arg("http://0.0.0.0:".to_owned() + port)
                .arg(&fastq_path1)
                .arg(&fastq_path2)
                .arg("--output")
                .arg(&output_path1)
                .arg("--output2")
                .arg(&output_path2)
                .assert()
                .success();
        });

        // Check both gzipped output files were created
        assert!(
            output_path1.exists(),
            "First gzipped output file wasn't created"
        );
        assert!(
            output_path2.exists(),
            "Second gzipped output file wasn't created"
        );

        assert!(
            fs::metadata(&output_path1).unwrap().len() > 0,
            "First gzipped output file is empty"
        );
        assert!(
            fs::metadata(&output_path2).unwrap().len() > 0,
            "Second gzipped output file is empty"
        );

        // Actually decompress and check if there are reads
        use flate2::read::GzDecoder;
        use std::fs::File;
        use std::io::Read;

        let file1 = File::open(&output_path1).unwrap();
        let mut gz1 = GzDecoder::new(file1);
        let mut contents1 = String::new();
        gz1.read_to_string(&mut contents1).unwrap();

        let file2 = File::open(&output_path2).unwrap();
        let mut gz2 = GzDecoder::new(file2);
        let mut contents2 = String::new();
        gz2.read_to_string(&mut contents2).unwrap();

        println!(
            "Gzip test - Output1 length: {}, Output2 length: {}",
            contents1.len(),
            contents2.len()
        );
        println!(
            "Gzip test - Output1 preview: {:?}",
            &contents1.chars().take(100).collect::<String>()
        );
        println!(
            "Gzip test - Output2 preview: {:?}",
            &contents2.chars().take(100).collect::<String>()
        );
    }

    #[test]
    fn test_filter_single_input_with_output2_warning() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("ref.fasta");
        let fastq_path = temp_dir.path().join("reads.fastq");
        let bin_path = temp_dir.path().join("ref.bin");
        let output_path = temp_dir.path().join("filtered.fastq");
        let output_path2 = temp_dir.path().join("filtered_2.fastq");

        create_test_fasta(&fasta_path);
        let fastq_content = "@seq1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n";
        fs::write(&fastq_path, fastq_content).unwrap();
        build_index(&fasta_path, &bin_path);

        // Run filtering command with output2 but no second input (should warn)
        run_with_server!(bin_path.clone(), |port: &str| {
            let mut cmd = Command::cargo_bin("deacon").unwrap();
            cmd.arg("client")
                .arg("http://0.0.0.0:".to_owned() + port)
                .arg(&fastq_path)
                .arg("--output")
                .arg(&output_path)
                .arg("-O")
                .arg(&output_path2)
                .assert()
                .success()
                .stderr(predicates::str::contains("Warning"));
        });

        // Check only the first output file was created
        assert!(output_path.exists(), "First output file wasn't created");
        assert!(
            !output_path2.exists(),
            "Second output file shouldn't be created for single input"
        );
    }
}

#[test]
fn test_shared_minimizer_counted_once() {
    // Catch bug where the same minimizer in different paired mates is counted twice
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fasta_path1 = temp_dir.path().join("reads_1.fasta");
    let fasta_path2 = temp_dir.path().join("reads_2.fasta");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fasta");
    let summary_path = temp_dir.path().join("summary.json");

    let ref_content = ">reference\nACGTACGTACGTACGTTGCATGCATGCATGCATAAGGTTAAGGTTAAGGTTAAGGTTCCCGGGCCCGGGCCCGGGCCCGGGATATATATATATATATATGCGCGCGCGCGCGCGCGC\n";
    fs::write(&fasta_path, ref_content).unwrap();

    // Create 120bp ref
    let ref_content = ">reference\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    fs::write(&fasta_path, ref_content).unwrap();

    // Create paired reads (80bp each) where both contain the same 60bp region from the reference
    // Shared region: first 60bp of reference (ACGT repeated 15 times)
    let fasta_content1 = ">read1/1\n\
        AAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAAAAAAA\n";

    let fasta_content2 = ">read1/2\n\
        TTTTTTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTT\n";

    fs::write(&fasta_path1, fasta_content1).unwrap();
    fs::write(&fasta_path2, fasta_content2).unwrap();

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // If shared minimizers are counted once (correct): total hits = 1, pair kept (1 < 2)
    // If shared minimizers are counted twice (bug): total hits = 2+, pair filtered (2+ >= 2)
    // Using --deplete to restore original behavior for this bug test
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--deplete")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fasta_path1)
            .arg(&fasta_path2)
            .arg("--output")
            .arg(&output_path)
            .arg("--summary")
            .arg(&summary_path)
            .arg("--abs-threshold")
            .arg("2")
            .arg("--rel-threshold")
            .arg("0.01") // Critical parameter: any pair with 2+ hits gets filtered
            .assert()
            .success();
    });

    assert!(output_path.exists(), "Output file wasn't created");
    assert!(summary_path.exists(), "Summary file wasn't created");

    let output_content = fs::read_to_string(&output_path).unwrap();
    let summary_content = fs::read_to_string(&summary_path).unwrap();
    let summary: serde_json::Value = serde_json::from_str(&summary_content).unwrap();

    // The reads should be kept because shared minimizers should only count once
    assert!(
        !output_content.is_empty(),
        "Read pair should be kept in output because shared minimizers should only count once. \
         Current implementation incorrectly counts them multiple times and filters the pair."
    );

    // Additional verification using the JSON summary
    let seqs_out = summary["seqs_out"].as_u64().unwrap();
    assert_eq!(
        seqs_out, 2,
        "Expected 2 sequences in output (both reads of the pair should be kept) \
         but got {seqs_out}. This indicates shared minimizers were double-counted."
    );
}

#[test]
fn test_filter_proportional_threshold() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_proportional.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("1")
            .arg("--rel-threshold")
            .arg("0.5") // 50% proportional threshold
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with proportional threshold wasn't created"
    );
}

#[test]
fn test_filter_proportional_paired() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path1 = temp_dir.path().join("reads_1.fastq");
    let fastq_path2 = temp_dir.path().join("reads_2.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_proportional_paired.fastq");

    create_test_fasta(&fasta_path);
    create_test_paired_fastq(&fastq_path1, &fastq_path2);
    build_index(&fasta_path, &bin_path);
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("1")
            .arg("--rel-threshold")
            .arg("0.3") // 30% proportional threshold
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path1)
            .arg(&fastq_path2)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output file with proportional threshold for paired reads wasn't created"
    );
}

#[test]
fn test_filter_edge_case_proportional_values() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_edge.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path.clone());

    // Work around moving values into closures for server running
    let bin_path1 = bin_path.clone();
    let bin_path2 = bin_path.clone();

    // Test with 0.0 (should pass everything)
    run_with_server!(bin_path1, |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("1")
            .arg("--rel-threshold")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Test with 1.0 (very strict)
    let output_path_strict = temp_dir.path().join("filtered_strict.fastq");
    run_with_server!(bin_path2, |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("--abs-threshold")
            .arg("1")
            .arg("--rel-threshold")
            .arg("1.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("--output")
            .arg(&output_path_strict)
            .assert()
            .success();
    });

    assert!(
        output_path.exists(),
        "Output with 0.0 threshold wasn't created"
    );
    assert!(
        output_path_strict.exists(),
        "Output with 1.0 threshold wasn't created"
    );
}

#[test]
fn test_multiline_fasta_matching() {
    let temp_dir = tempdir().unwrap();
    let ref_path = temp_dir.path().join("ref.fasta");
    let query_path = temp_dir.path().join("query.fasta");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("output.fasta");

    let reference_fasta = ">ref\nACGTTTAAGGCCAACCACACACACACACATT\n";
    let query_fasta = ">query\nACGTTTAAGGCCAACC\nACACACACACACATT\n";

    fs::write(&ref_path, reference_fasta).unwrap();
    fs::write(&query_path, query_fasta).unwrap();

    // Build index with k=31, w=1
    let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
        .arg("index")
        .arg("build")
        .arg("-k")
        .arg("31")
        .arg("-w")
        .arg("1")
        .arg(&ref_path)
        .output()
        .expect("Failed to execute index command");

    fs::write(&bin_path, output.stdout).expect("Failed to write index file");
    assert!(output.status.success(), "Index build command failed");

    // client with -a 1
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("-a")
            .arg("1")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&query_path)
            .arg("-o")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Verify that mid record newline doesn't break match
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        !output_content.is_empty(),
        "Multiline FASTA should match indexed sequence"
    );
    assert!(
        output_content.contains(">query"),
        "Output should contain query header"
    );
    assert!(
        output_content.contains("ACGTTTAAGGCCAACCACACACACACACATT"),
        "Output should contain the full sequence"
    );
}

#[test]
fn test_newline_mapping_bug() {
    let temp_dir = tempdir().unwrap();
    let ref_path = temp_dir.path().join("reference.fa");
    let query_path = temp_dir.path().join("query.fa");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("output.fa");

    // Create reference file with sequence split across lines
    // The newlines should be stripped but if they're not, they'll be mapped to 'C'
    let ref_content = ">reference\nAAAAA\nAAAAA\nAAAAA\nAAAAA\n";
    fs::write(&ref_path, ref_content).unwrap();

    // Create query file with Cs where newlines would be
    let query_content = ">query\nAAAAACAAAAACAAAAACAAAAA\n";
    fs::write(&query_path, query_content).unwrap();

    // Build index with k=5, w=5 (k+w-1 must be odd: 5+5-1=9, odd ✓)
    let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
        .arg("index")
        .arg("build")
        .arg("-k")
        .arg("5")
        .arg("-w")
        .arg("5")
        .arg(&ref_path)
        .output()
        .expect("Failed to execute index command");

    fs::write(&bin_path, output.stdout).expect("Failed to write index file");
    assert!(output.status.success(), "Index build command failed");

    // client query against index
    run_with_server!(bin_path.clone(), |port: &str| {
        let mut cmd = Command::cargo_bin("deacon").unwrap();
        cmd.arg("client")
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&query_path)
            .arg("-o")
            .arg(&output_path)
            .assert()
            .success();
    });

    // Read filtered output
    let output_str = fs::read_to_string(&output_path).unwrap();

    // If newlines are being mapped to C, the query would match
    // The bug would cause the reference "AAAAA\nAAAAA\nAAAAA\nAAAAA" to become
    // "AAAAACAAAAACAAAAACAAAAA" after mapping newlines to C
    // So if the bug exists, the query would match and be filtered (kept with deplete=false)

    // With the bug, we'd expect a match. Without the bug, no match.
    if output_str.contains(">query") {
        panic!(
            "BUG DETECTED: Query matched due to newlines being mapped to 'C'. Output: {output_str}"
        );
    }

    println!("Test passed - no false matches from newline mapping");
}

#[test]
fn test_large_kmer_filter() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");
    let fastq_path = temp_dir.path().join("test.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);

    // Index with k=41 (u128 code path)
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg("-k")
        .arg("41")
        .arg("-w")
        .arg("15")
        .arg(&fasta_path)
        .arg("-o")
        .arg(&bin_path)
        .assert()
        .success();

    // Test filtering with our k=41 index
    run_with_server!(bin_path.clone(), |port: &str| {
        let output = Command::cargo_bin("deacon")
            .unwrap()
            .arg("client")
            .arg("http://0.0.0.0:".to_owned() + port)
            .arg(&fastq_path)
            .arg("-a")
            .arg("1")
            .arg("-r")
            .arg("0.0")
            .output()
            .unwrap();

        assert!(output.status.success(), "client command failed");

        // Should retain both seqs
        let stdout = String::from_utf8_lossy(&output.stdout);
        let num_sequences = stdout.lines().filter(|line| line.starts_with('@')).count();
        assert_eq!(num_sequences, 2, "Should retain both sequences");
    });
}
