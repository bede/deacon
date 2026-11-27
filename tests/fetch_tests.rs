use assert_cmd::assert::Assert;
use assert_cmd::cargo;
use predicates::prelude::PredicateBooleanExt;
use predicates::str;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::path::Path;
use std::thread;
use std::time::{Duration, Instant};
use tempfile::tempdir;

const EXACT: &str = "panhuman-1.k31w15.idx";
const FUSE: &str = "panhuman-1.k31w15.pidx";
const MANIFEST: &[u8] = br#"{
  "schema_version": 1,
  "default_index": "different name.idx",
  "indexes": {
    "panhuman-1.k31w15": {
      "k": 31,
      "w": 15,
      "default_asset": "panhuman-1.k31w15.idx",
      "assets": {
        "panhuman-1.k31w15.idx": {
          "size_bytes": 3,
          "sha256": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
          "opts": {"algorithm": "kdust", "c": 0.6}
        },
        "panhuman-1.k31w15.pidx": {
          "size_bytes": 5,
          "sha256": "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
          "opts": {"fingerprint_bits": 16, "c": 0.7}
        }
      }
    },
    "panhuman-1.k31w1": {
      "k": 31,
      "w": 1,
      "default_asset": "different name.idx",
      "assets": {
        "different name.idx": {
          "size_bytes": 3,
          "sha256": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        },
        "different name.pidx": {
          "size_bytes": 5,
          "sha256": "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        }
      }
    }
  },
  "aliases": {"panhuman-1": "panhuman-1.k31w15", "small-default": "panhuman-1.k31w1"}
}"#;

fn run(cwd: &Path, args: &[&str], manifest: &[u8], asset: Option<(&str, &[u8])>) -> Assert {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    listener.set_nonblocking(true).unwrap();
    let url = format!(
        "http://{}/catalog/test.json",
        listener.local_addr().unwrap()
    );
    thread::scope(|scope| {
        let server = scope.spawn(move || {
            for (path, body) in std::iter::once(("test.json", manifest)).chain(asset) {
                let start = Instant::now();
                let mut stream = loop {
                    match listener.accept() {
                        Ok((stream, _)) => break stream,
                        Err(err)
                            if err.kind() == std::io::ErrorKind::WouldBlock
                                && start.elapsed() < Duration::from_secs(5) =>
                        {
                            thread::sleep(Duration::from_millis(5));
                        }
                        Err(err) => panic!("No request for {path}: {err}"),
                    }
                };
                stream
                    .set_read_timeout(Some(Duration::from_secs(5)))
                    .unwrap();
                let mut reader = BufReader::new(&mut stream);
                let mut line = String::new();
                reader.read_line(&mut line).unwrap();
                assert_eq!(line.trim_end(), format!("GET /catalog/{path} HTTP/1.1"));
                loop {
                    line.clear();
                    if reader.read_line(&mut line).unwrap() == 0 || line == "\r\n" {
                        break;
                    }
                }
                write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    body.len()
                )
                .unwrap();
                stream.write_all(body).unwrap();
            }
        });
        let result = cargo::cargo_bin_cmd!("deacon")
            .env("DEACON_REPOSITORY_URL", &url)
            .current_dir(cwd)
            .arg("index")
            .args(args)
            .assert();
        server.join().unwrap();
        result
    })
}

#[test]
fn resolves_alias_group_and_filename() {
    for name in ["panhuman-1", "panhuman-1.k31w15", EXACT] {
        let dir = tempdir().unwrap();
        run(
            dir.path(),
            &["fetch", name],
            MANIFEST,
            Some((EXACT, b"abc")),
        )
        .success()
        .stderr(str::contains("Using default index").not());
        assert_eq!(fs::read(dir.path().join(EXACT)).unwrap(), b"abc");
    }

    let dir = tempdir().unwrap();
    run(
        dir.path(),
        &["fetch", FUSE],
        MANIFEST,
        Some((FUSE, b"hello")),
    )
    .success();
    assert_eq!(fs::read(dir.path().join(FUSE)).unwrap(), b"hello");

    for name in ["panhuman-1.k31w1", "different name.idx"] {
        let dir = tempdir().unwrap();
        run(
            dir.path(),
            &["fetch", name],
            MANIFEST,
            Some(("different%20name.idx", b"abc")),
        )
        .success();
        assert_eq!(
            fs::read(dir.path().join("different name.idx")).unwrap(),
            b"abc"
        );
    }
}

#[test]
fn fetch_uses_manifest_default_when_no_name_is_given() {
    let dir = tempdir().unwrap();
    run(
        dir.path(),
        &["fetch"],
        MANIFEST,
        Some(("different%20name.idx", b"abc")),
    )
    .success()
    .stderr(str::contains("Using default index different name.idx"));
    assert_eq!(
        fs::read(dir.path().join("different name.idx")).unwrap(),
        b"abc"
    );
    run(dir.path(), &["fetch", "--quiet"], MANIFEST, None)
        .success()
        .stderr(str::contains("Using default index").not());

    let mut manifest: serde_json::Value = serde_json::from_slice(MANIFEST).unwrap();
    manifest.as_object_mut().unwrap().remove("default_index");
    let body = serde_json::to_vec(&manifest).unwrap();
    run(dir.path(), &["fetch"], &body, None)
        .failure()
        .stderr(str::contains("The manifest has no default index"));
    run(dir.path(), &["list"], &body, None)
        .success()
        .stdout(str::contains("default fetch:").not());
    run(dir.path(), &["fetch", EXACT], &body, Some((EXACT, b"abc"))).success();

    for name in ["small-default", "panhuman-1.k31w1", "unknown.idx"] {
        manifest["default_index"] = name.into();
        run(
            dir.path(),
            &["fetch"],
            &serde_json::to_vec(&manifest).unwrap(),
            None,
        )
        .failure()
        .stderr(str::contains("default_index must name an asset filename"));
    }
}

#[test]
fn list_uses_only_the_manifest_and_orders_options() {
    let dir = tempdir().unwrap();
    run(dir.path(), &["list"], MANIFEST, None)
        .success()
        .stdout(str::starts_with(
            "panhuman-1.k31w1 (aliases: small-default, default fetch: different name.idx)\n",
        ));

    let mut manifest: serde_json::Value = serde_json::from_slice(MANIFEST).unwrap();
    manifest["default_index"] = EXACT.into();
    manifest["indexes"]["panhuman-1.k31w15"]["assets"][FUSE]["size_bytes"] =
        1_845_493_803_u64.into();
    run(
        dir.path(),
        &["list"],
        &serde_json::to_vec(&manifest).unwrap(),
        None,
    )
    .success()
    .stdout(str::starts_with(
        "panhuman-1.k31w15 (aliases: panhuman-1, default fetch: .idx)\n    panhuman-1.k31w15.idx",
    ))
    .stdout(str::contains("(default)").not())
    .stdout(str::contains("3 B  k=31 w=15 c=0.6 algorithm=\"kdust\"\n"))
    .stdout(str::contains(
        "1.8 GB  k=31 w=15 c=0.7 fingerprint_bits=16\npanhuman-1.k31w1 (aliases: small-default)\n",
    ));
}

#[test]
fn verification_protects_existing_outputs() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("output.idx");
    let args = ["fetch", "panhuman-1", "-o", "output.idx"];
    let forced = ["fetch", "panhuman-1", "-o", "output.idx", "--force"];

    fs::write(&path, b"abc").unwrap();
    run(dir.path(), &args, MANIFEST, None)
        .success()
        .stderr(str::contains("skipping download"));
    fs::write(&path, b"bad").unwrap();
    run(dir.path(), &args, MANIFEST, None)
        .failure()
        .stderr(str::contains("--force to replace"));
    run(dir.path(), &forced, MANIFEST, Some((EXACT, b"xyz")))
        .failure()
        .stderr(str::contains("SHA-256 mismatch"));
    run(dir.path(), &forced, MANIFEST, Some((EXACT, b"ab")))
        .failure()
        .stderr(str::contains("Size mismatch"));
    assert_eq!(fs::read(&path).unwrap(), b"bad");
    assert_eq!(dir.path().read_dir().unwrap().count(), 1);

    run(dir.path(), &forced, MANIFEST, Some((EXACT, b"abc"))).success();
    assert_eq!(fs::read(&path).unwrap(), b"abc");
}

#[test]
fn explicit_group_default_and_manifest_errors() {
    let dir = tempdir().unwrap();
    let mut manifest: serde_json::Value = serde_json::from_slice(MANIFEST).unwrap();
    manifest["default_index"] = EXACT.into();
    manifest["indexes"]["panhuman-1.k31w15"]["default_asset"] = FUSE.into();
    let body = serde_json::to_vec(&manifest).unwrap();
    run(
        dir.path(),
        &["fetch", "panhuman-1"],
        &body,
        Some((FUSE, b"hello")),
    )
    .success();
    run(dir.path(), &["fetch"], &body, Some((EXACT, b"abc"))).success();
    run(dir.path(), &["fetch", "unknown"], MANIFEST, None)
        .failure()
        .stderr(str::contains("Unknown index \"unknown\""));
    run(dir.path(), &["list"], br#"{"schema_version":2}"#, None)
        .failure()
        .stderr(str::contains("Unsupported index manifest schema version 2"));
    manifest["aliases"]["panhuman-1"] = "unknown".into();
    run(
        dir.path(),
        &["list"],
        &serde_json::to_vec(&manifest).unwrap(),
        None,
    )
    .failure()
    .stderr(str::contains("points to unknown index group"));
}

#[test]
fn multi_asset_groups_require_an_explicit_default() {
    let dir = tempdir().unwrap();
    let mut manifest: serde_json::Value = serde_json::from_slice(MANIFEST).unwrap();
    let group = manifest["indexes"]["panhuman-1.k31w15"]
        .as_object_mut()
        .unwrap();
    group.remove("default_asset");
    run(
        dir.path(),
        &["list"],
        &serde_json::to_vec(&manifest).unwrap(),
        None,
    )
    .failure()
    .stderr(str::contains(
        "needs one asset or an explicit default_asset",
    ));

    manifest["indexes"]["panhuman-1.k31w15"]["assets"]
        .as_object_mut()
        .unwrap()
        .remove(FUSE);
    run(
        dir.path(),
        &["fetch", "panhuman-1"],
        &serde_json::to_vec(&manifest).unwrap(),
        Some((EXACT, b"abc")),
    )
    .success();
}
