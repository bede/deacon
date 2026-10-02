//! Python bindings for deacon: load a minimizer index once, filter many files.
#![allow(clippy::too_many_arguments)]

use std::path::{Path, PathBuf};
use std::sync::Arc;

use ::deacon::{
    DEFAULT_CBQ_BLOCK_SIZE_MIB, FilterConfig, FilterParams, Index as DeaconIndex, IndexKind,
    filter_files, index_fetch, load_filter_index,
};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

// Keep the literal in the Python signature/stub while detecting a changed core default at compile time.
const _: () = assert!(DEFAULT_CBQ_BLOCK_SIZE_MIB == 16);

fn to_pyerr(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// A loaded minimizer index, reusable across many `filter` calls.
#[pyclass(frozen)]
struct Index {
    label: String,
    index: Arc<DeaconIndex>,
}

impl Index {
    fn load(path: &Path, complexity_threshold: Option<f32>) -> PyResult<Self> {
        let index = load_filter_index(path, complexity_threshold).map_err(to_pyerr)?;
        Ok(Index {
            label: path.to_string_lossy().into_owned(),
            index: Arc::new(index),
        })
    }
}

#[pymethods]
impl Index {
    #[new]
    #[pyo3(signature = (path, /, *, complexity_threshold=None))]
    fn new(path: PathBuf, complexity_threshold: Option<f32>) -> PyResult<Self> {
        Self::load(&path, complexity_threshold)
    }

    /// Download a prebuilt index, then load and return it.
    #[staticmethod]
    #[pyo3(signature = (*, name="panhuman-1", k=31, w=15, output=None, complexity_threshold=None))]
    fn fetch(
        name: &str,
        k: u8,
        w: u8,
        output: Option<PathBuf>,
        complexity_threshold: Option<f32>,
    ) -> PyResult<Self> {
        let path = index_fetch(name, k, w, output.as_deref()).map_err(to_pyerr)?;
        Index::load(&path, complexity_threshold)
    }

    /// Index metadata: k, w, format and minimizer/key count.
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let format = match self.index.kind() {
            IndexKind::Exact if self.index.kmer_length() <= 32 => "exact-u64",
            IndexKind::Exact => "exact-u128",
            IndexKind::Fuse { .. } => "bff",
        };
        let d = PyDict::new(py);
        d.set_item("k", self.index.kmer_length())?;
        d.set_item("w", self.index.window_size())?;
        d.set_item("format", format)?;
        d.set_item("count", self.index.len())?;
        Ok(d)
    }

    #[pyo3(signature = (
        input,
        /,
        *,
        input2=None,
        interleaved=false,
        check_pairs=false,
        deplete=false,
        rename=false,
        output=None,
        output2=None,
        summary=None,
        abs_threshold=2,
        rel_threshold=0.01,
        prefix_length=0,
        discard_quality=false,
        ordered=false,
        threads=8,
        compression_level=2,
        compression_threads=0,
        cbq_block_size=16,
        quiet=true,
        debug=None,
    ))]
    fn filter(
        &self,
        py: Python<'_>,
        input: PathBuf,
        input2: Option<PathBuf>,
        interleaved: bool,
        check_pairs: bool,
        deplete: bool,
        rename: bool,
        output: Option<PathBuf>,
        output2: Option<PathBuf>,
        summary: Option<PathBuf>,
        abs_threshold: usize,
        rel_threshold: f64,
        prefix_length: usize,
        discard_quality: bool,
        ordered: bool,
        threads: u16,
        compression_level: u8,
        compression_threads: u16,
        cbq_block_size: u16,
        quiet: bool,
        debug: Option<PathBuf>,
    ) -> PyResult<Py<PyDict>> {
        let cfg = FilterConfig {
            input_path: input,
            input2_path: input2,
            interleaved,
            check_pairs,
            output_path: output,
            output2_path: output2,
            params: FilterParams {
                abs_threshold,
                rel_threshold,
                prefix_length,
                deplete,
            },
            summary_path: summary,
            rename,
            discard_quality,
            ordered,
            threads,
            compression_level,
            cbq_block_size,
            compression_threads,
            debug,
            progress: !quiet,
        };

        cfg.validate()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let index = Arc::clone(&self.index);
        let summary = py
            .detach(|| filter_files(index, &self.label, None, &cfg))
            .map_err(to_pyerr)?;
        Ok(pythonize::pythonize(py, &summary)?
            .cast_into::<PyDict>()?
            .unbind())
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
