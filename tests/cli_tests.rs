use assert_cmd::cargo;
use deacon::{
    ComplexityAlgorithm, FilterKernel, FilterParams, MinimizerSet, RapidHashSet,
    validate_unit_interval,
};
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
fn unit_interval_validation() {
    for value in [0.0, 1.0] {
        assert!(validate_unit_interval("threshold", value).is_ok());
    }
    for value in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(validate_unit_interval("threshold", value).is_err());
    }

    // Library entry points must apply the same validation
    let mut set = MinimizerSet::U64(RapidHashSet::default());
    assert!(
        set.retain_complexity(31, ComplexityAlgorithm::Kdust, f32::NAN, false)
            .is_err()
    );
    assert!(
        FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                abs_threshold: 1,
                rel_threshold: f64::NAN,
                prefix_length: 0,
            },
        )
        .is_err()
    );
}

#[test]
fn cli_rejects_invalid_thresholds() {
    // So clap does not read a negative value as a flag
    for args in [
        ["filter", "missing.idx", "--rel-threshold=1.1"],
        ["filter", "missing.idx", "--rel-threshold=-0.1"],
        ["filter", "missing.idx", "--complexity-threshold=NaN"],
        ["index", "filter", "--complexity-threshold=inf"],
    ] {
        cargo::cargo_bin_cmd!("deacon")
            .args(args)
            .assert()
            .failure()
            .stderr(str::contains("must be between 0.0 and 1.0 inclusive"));
    }
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
    let _ = server.wait();
    let _ = fs::remove_file("deacon_server_socket");
}
