use deacon::{
    FilterConfig, FilterKernel, FilterParams, Index, IndexKind, MinimizerVec, Minimizers,
    filter_files, load_index, write_index,
};
use std::{fs, io::Cursor, sync::Arc};

const REFERENCE: &[u8] = b"ACGTTGCAAGGCTTAACCGGTTACGATCGATCGGATCCTAGCTAGCTTAACCGGATCGTAACGTGCTAGCATCGATGCTAGCTACGATCGATCG";

fn reference_index(k: u8, w: u8) -> Index {
    let minimizers = Minimizers::new(k, w).unwrap().compute(REFERENCE).clone();
    Index::from_minimizers(k, w, minimizers).unwrap()
}

#[test]
fn shared_workers_score_distinct_hits_for_both_widths() {
    for (k, w) in [(7, 3), (33, 3), (61, 3)] {
        let index = Arc::new(reference_index(k, w));
        let params = FilterParams {
            abs_threshold: 1,
            ..FilterParams::default()
        };
        let mut filter = FilterKernel::new(Arc::clone(&index), params).unwrap();
        let mut deplete = FilterKernel::new(
            Arc::clone(&index),
            FilterParams {
                deplete: true,
                ..params
            },
        )
        .unwrap();

        for read in [
            REFERENCE,
            b"ACG",
            b"",
            b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        ] {
            let fast = filter.classify_read(read);
            let score = filter.score_read(read);
            let diagnostics = filter.classify_read_with_diagnostics(read);
            assert_eq!(score.is_match, fast.is_match);
            assert_eq!(score.keep, fast.keep);
            assert_eq!(score, diagnostics.score);
            assert_eq!(score.hit_count, diagnostics.hit_kmers.len());
            let inverse = deplete.score_read(read);
            assert_eq!(score.is_match, inverse.is_match);
            assert_ne!(score.keep, inverse.keep);
        }

        // Duplicate mates count once.
        assert_eq!(
            filter.score_read(REFERENCE),
            filter.score_pair(REFERENCE, REFERENCE)
        );
        let diagnostics = filter.classify_pair_with_diagnostics(REFERENCE, REFERENCE);
        assert_eq!(diagnostics.score, filter.score_pair(REFERENCE, REFERENCE));
        assert_eq!(diagnostics.hit_kmers.len(), diagnostics.score.hit_count);
        let mut worker = FilterKernel::new(Arc::clone(&index), params).unwrap();
        let expected = filter.score_read(REFERENCE);
        std::thread::scope(|scope| {
            assert_eq!(
                scope
                    .spawn(move || worker.score_read(REFERENCE))
                    .join()
                    .unwrap(),
                expected
            );
        });
    }
}

#[test]
fn exact_writer_preserves_index_and_roundtrips_both_widths() {
    for (k, w) in [(7, 3), (33, 3), (61, 3)] {
        let index = reference_index(k, w);
        let count = index.len();
        let mut bytes = Vec::new();
        write_index(&index, &mut bytes).unwrap();
        assert_eq!(index.len(), count);
        let mut reader = Cursor::new([bytes.as_slice(), bytes.as_slice()].concat());
        let loaded = load_index(&mut reader).unwrap();
        assert_eq!(reader.position(), bytes.len() as u64);
        assert_eq!(load_index(&mut reader).unwrap().len(), count);
        assert_eq!(loaded.kind(), IndexKind::Exact);
        assert_eq!(loaded.len(), count);
        let mut filter = FilterKernel::new(Arc::new(loaded), FilterParams::default()).unwrap();
        assert!(filter.classify_read(REFERENCE).is_match);
        let mut second = Vec::new();
        write_index(&index, &mut second).unwrap();
        assert_eq!(bytes, second);
    }
}

#[test]
fn invalid_options_fail_before_outputs_are_touched() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fa");
    let output = directory.path().join("out.fa");
    fs::write(&input, [b">read\n".as_slice(), REFERENCE, b"\n"].concat()).unwrap();
    fs::write(&output, b"preserve this").unwrap();
    let index = Arc::new(reference_index(7, 3));
    let mut config = FilterConfig::new(input);
    config.output_path = Some(output.clone());
    for (params, expected_error) in [
        (
            FilterParams {
                abs_threshold: 0,
                ..FilterParams::default()
            },
            "abs_threshold must be at least 1",
        ),
        (
            FilterParams {
                rel_threshold: f64::NAN,
                ..FilterParams::default()
            },
            "relative threshold must be between 0.0 and 1.0 inclusive",
        ),
    ] {
        config.params = params;
        let error = filter_files(Arc::clone(&index), "reference", None, &config).unwrap_err();
        assert!(error.to_string().contains(expected_error), "{error:#}");
        assert_eq!(fs::read(&output).unwrap(), b"preserve this");
    }
    for (k, w) in [(0, 2), (2, 2), (3, 0), (31, 66)] {
        assert!(Index::from_minimizers(k, w, MinimizerVec::U64(vec![])).is_err());
        assert!(Minimizers::new(k, w).is_err());
    }
    assert!(Index::from_minimizers(33, 3, MinimizerVec::U64(vec![])).is_err());
    assert!(Index::from_minimizers(7, 3, MinimizerVec::U64(vec![1 << 14])).is_err());
}

#[test]
fn repeated_file_filtering_leaves_global_rayon_pool_to_caller() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.fa");
    let output = directory.path().join("out.fa");
    fs::write(&input, [b">read\n".as_slice(), REFERENCE, b"\n"].concat()).unwrap();
    let index = Arc::new(reference_index(7, 3));
    let mut config = FilterConfig::new(input);
    config.output_path = Some(output);
    for threads in [1, 3] {
        config.threads = threads;
        let summary = filter_files(Arc::clone(&index), "reference", None, &config).unwrap();
        assert_eq!(summary.seqs_in, 1);
        assert_eq!(summary.seqs_out, 1);
    }
    // Filtering must leave the global pool uninitialized.
    rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build_global()
        .unwrap();
    assert_eq!(rayon::current_num_threads(), 2);
}

#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn paired_non_utf8_paths_work_with_compression_and_summary() {
    use std::{ffi::OsString, os::unix::ffi::OsStringExt};
    let directory = tempfile::tempdir().unwrap();
    let path = |suffix: &[u8]| {
        directory.path().join(OsString::from_vec(
            [b"read-\xff".as_slice(), suffix].concat(),
        ))
    };
    let (input1, input2, output1, output2) = (
        path(b"1.fa"),
        path(b"2.fa"),
        path(b"1.fa.gz"),
        path(b"2.fa.gz"),
    );
    fs::write(
        &input1,
        [b">read/1\n".as_slice(), REFERENCE, b"\n"].concat(),
    )
    .unwrap();
    fs::write(
        &input2,
        [b">read/2\n".as_slice(), REFERENCE, b"\n"].concat(),
    )
    .unwrap();
    let mut config = FilterConfig::new(input1.clone());
    config.input2_path = Some(input2.clone());
    config.output_path = Some(output1.clone());
    config.output2_path = Some(output2.clone());
    config.summary_path = Some(directory.path().join("summary.json"));
    config.check_pairs = true;
    config.threads = 3;
    let summary =
        filter_files(Arc::new(reference_index(7, 3)), "reference", None, &config).unwrap();
    assert_eq!(summary.seqs_out, 2);
    for (output, input) in [(output1, input1), (output2, input2)] {
        use std::io::Read;
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(output).unwrap())
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, fs::read(input).unwrap());
    }
    let saved: serde_json::Value =
        serde_json::from_slice(&fs::read(config.summary_path.unwrap()).unwrap()).unwrap();
    assert_eq!(saved["seqs_out"], 2);
}
