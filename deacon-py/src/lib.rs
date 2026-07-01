//! Python bindings for deacon: load a minimizer index once, filter many files.
#![allow(clippy::too_many_arguments)]

use std::path::{Path, PathBuf};
use std::sync::Arc;

use ::deacon::{
    FilterRunConfig, IndexHeader, MinimizerSet, index_fetch, load_index_from_path_auto,
    run_with_index,
};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn to_pyerr(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// A loaded minimizer index, reusable across many `filter` calls.
#[pyclass(frozen)]
struct Index {
    label: String,
    k: u8,
    w: u8,
    minimizers: Arc<MinimizerSet>,
}

#[pymethods]
impl Index {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let (minimizers, header) = load_index_from_path_auto(Path::new(path)).map_err(to_pyerr)?;
        Ok(Index {
            label: path.to_string(),
            k: header.kmer_length(),
            w: header.window_size(),
            minimizers: Arc::new(minimizers),
        })
    }

    /// Download a prebuilt index, then load and return it.
    #[staticmethod]
    #[pyo3(signature = (name="panhuman-1", k=31, w=15, output=None))]
    fn fetch(name: &str, k: u8, w: u8, output: Option<String>) -> PyResult<Self> {
        let out_path = output
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from(format!("{name}.k{k}w{w}.idx")));
        index_fetch(name, k, w, Some(&out_path)).map_err(to_pyerr)?;
        Index::new(&out_path.to_string_lossy())
    }

    /// Index metadata: k, w, format and minimizer/key count.
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let format = match &*self.minimizers {
            MinimizerSet::U64(_) => "exact-u64",
            MinimizerSet::U128(_) => "exact-u128",
            MinimizerSet::Fuse(_) => "bff",
        };
        let d = PyDict::new(py);
        d.set_item("k", self.k)?;
        d.set_item("w", self.w)?;
        d.set_item("format", format)?;
        d.set_item("count", self.minimizers.len())?;
        Ok(d)
    }

    #[pyo3(signature = (
        fastq,
        fastq2=None,
        deplete=false,
        rename=false,
        rename_random=false,
        output=None,
        output2=None,
        abs_threshold=2,
        rel_threshold=0.01,
        prefix_length=0,
        output_fasta=false,
        threads=8,
        compression_level=2,
        compression_threads=0,
        debug=false,
        quiet=true,
    ))]
    fn filter(
        &self,
        py: Python<'_>,
        fastq: String,
        fastq2: Option<String>,
        deplete: bool,
        rename: bool,
        rename_random: bool,
        output: Option<String>,
        output2: Option<String>,
        abs_threshold: usize,
        rel_threshold: f64,
        prefix_length: usize,
        output_fasta: bool,
        threads: u16,
        compression_level: u8,
        compression_threads: u16,
        debug: bool,
        quiet: bool,
    ) -> PyResult<Py<PyAny>> {
        let cfg = FilterRunConfig {
            input_path: fastq,
            input2_path: fastq2,
            output_path: output.map(PathBuf::from),
            output2_path: output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary_path: None,
            deplete,
            rename,
            rename_random,
            output_fasta,
            threads,
            compression_level,
            compression_threads,
            debug,
            quiet,
            index_label: self.label.clone(),
        };

        let mins = Arc::clone(&self.minimizers);
        let (k, w) = (self.k, self.w);
        let summary = py
            .detach(|| run_with_index(&mins, &IndexHeader::new(k, w), &cfg))
            .map_err(to_pyerr)?;
        Ok(pythonize::pythonize(py, &summary)?.unbind())
    }
}

// Declarative module form so pyo3 introspection can link members (see scripts/gen-py-stubs.sh).
#[pymodule]
mod _deacon {
    use pyo3::prelude::*;

    #[pymodule_export]
    use super::Index;

    #[pymodule_init]
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        // Fail cleanly (Python exception) on CPUs lacking the compiled SIMD features.
        ensure_simd::ensure_simd();
        m.add("__version__", env!("CARGO_PKG_VERSION"))?;
        Ok(())
    }
}
