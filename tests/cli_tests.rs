use assert_cmd::cargo;
use predicates::str;
use std::fs;
use std::process::{Child, Command as StdCommand};
use std::thread;
use std::time::Duration;
use tempfile::tempdir;

#[test]
fn test_version() {
    let mut cmd = cargo::cargo_bin_cmd!("deacon");
    cmd.arg("--version")
        .assert()
        .success()
        .stdout(str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn test_no_args() {
    let mut cmd = cargo::cargo_bin_cmd!("deacon");
    cmd.assert().failure().stderr(str::contains("Usage"));
}

#[test]
fn test_server_mode() {
    let temp_dir = tempdir().unwrap();
    let ref_fasta = temp_dir.path().join("ref.fa");
    let index_path = temp_dir.path().join("ref.idx");
    let test_fasta = temp_dir.path().join("test.fa");
    let output_path = temp_dir.path().join("out.fa");

    // Create reference and build index
    fs::write(
        &ref_fasta,
        ">ref\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n",
    )
    .unwrap();
    let index_output = StdCommand::new(cargo::cargo_bin!("deacon"))
        .arg("index")
        .arg("build")
        .arg(&ref_fasta)
        .output()
        .unwrap();
    fs::write(&index_path, index_output.stdout).unwrap();

    // Create test fasta
    fs::write(
        &test_fasta,
        ">test\nATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT\n",
    )
    .unwrap();

    // Start server
    let mut server: Child = StdCommand::new(cargo::cargo_bin!("deacon"))
        .arg("server")
        .arg("start")
        .spawn()
        .unwrap();

    thread::sleep(Duration::from_millis(500));

    // Filter via server
    StdCommand::new(cargo::cargo_bin!("deacon"))
        .arg("--use-server")
        .arg("filter")
        .arg(&index_path)
        .arg(&test_fasta)
        .arg("-o")
        .arg(&output_path)
        .arg("-a")
        .arg("1")
        .arg("-r")
        .arg("0")
        .output()
        .unwrap();

    assert!(output_path.exists());

    // Stop server
    StdCommand::new(cargo::cargo_bin!("deacon"))
        .arg("--use-server")
        .arg("server")
        .arg("stop")
        .output()
        .unwrap();

    let _ = server.kill();
    let _ = fs::remove_file("deacon_server_socket");
}
