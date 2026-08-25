//! Python bindings for deacon: load a minimizer index once, filter many files.
#![allow(clippy::too_many_arguments)]

use std::path::{Path, PathBuf};
use std::sync::Arc;

use ::deacon::{
    ComplexityAlgorithm, DEFAULT_CBQ_BLOCK_SIZE_MIB, FilterRunConfig, IndexHeader, MinimizerSet,
    index_fetch, load_index_from_path_auto, run_with_index,
};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

// Keep the literal in the Python signature/stub while detecting a changed core default at compile time.
const _: () = assert!(DEFAULT_CBQ_BLOCK_SIZE_MIB == 16);

fn to_pyerr(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

fn path_to_string(path: PathBuf, argument: &str) -> PyResult<String> {
    path.into_os_string().into_string().map_err(|_| {
        PyValueError::new_err(format!("{argument} must be representable as valid UTF-8"))
    })
}

/// A loaded minimizer index, reusable across many `filter` calls.
#[pyclass(frozen)]
struct Index {
    label: String,
    k: u8,
    w: u8,
    minimizers: Arc<MinimizerSet>,
}

impl Index {
    fn load(path: &Path, complexity_threshold: Option<f32>) -> PyResult<Self> {
        let (mut minimizers, header) = load_index_from_path_auto(path).map_err(to_pyerr)?;
        // Discard low-complexity index minimizers once at load (kdust); reused across filters.
        if let Some(threshold) = complexity_threshold {
            if matches!(minimizers, MinimizerSet::Fuse(_)) {
                return Err(PyRuntimeError::new_err(
                    "complexity filtering is not supported on BFF indexes; use an exact index",
                ));
            }
            minimizers
                .retain_complexity(
                    header.kmer_length(),
                    ComplexityAlgorithm::Kdust,
                    threshold,
                    false,
                )
                .map_err(to_pyerr)?;
        }
        Ok(Index {
            label: path.to_string_lossy().into_owned(),
            k: header.kmer_length(),
            w: header.window_size(),
            minimizers: Arc::new(minimizers),
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
        let out_path = output.unwrap_or_else(|| PathBuf::from(format!("{name}.k{k}w{w}.idx")));
        index_fetch(name, k, w, Some(&out_path)).map_err(to_pyerr)?;
        Index::load(&out_path, complexity_threshold)
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
        inverse_output=None,
        inverse_output2=None,
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
        debug=false,
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
        inverse_output: Option<PathBuf>,
        inverse_output2: Option<PathBuf>,
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
        debug: bool,
    ) -> PyResult<Py<PyDict>> {
        if interleaved && input2.is_some() {
            return Err(PyValueError::new_err(
                "interleaved cannot be combined with input2 (interleaved input is a single file/stream)",
            ));
        }
        if abs_threshold == 0 {
            return Err(PyValueError::new_err("abs_threshold must be at least 1"));
        }
        if !(1..=1024).contains(&cbq_block_size) {
            return Err(PyValueError::new_err(
                "cbq_block_size must be between 1 and 1024 MiB inclusive",
            ));
        }

        let input = path_to_string(input, "input")?;
        let input2 = input2
            .map(|path| path_to_string(path, "input2"))
            .transpose()?;
        let output2 = output2
            .map(|path| path_to_string(path, "output2"))
            .transpose()?;
        let inverse_output2 = inverse_output2
            .map(|path| path_to_string(path, "inverse_output2"))
            .transpose()?;

        let cfg = FilterRunConfig {
            input_path: input,
            input2_path: input2,
            interleaved,
            check_pairs,
            output_path: output,
            output2_path: output2,
            inverse_output_path: inverse_output,
            inverse_output2_path: inverse_output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary_path: summary,
            deplete,
            rename,
            discard_quality,
            ordered,
            threads,
            compression_level,
            cbq_block_size,
            compression_threads,
            debug,
            quiet,
            index_label: self.label.clone(),
        };

        let mins = Arc::clone(&self.minimizers);
        let (k, w) = (self.k, self.w);
        let summary = py
            .detach(|| run_with_index(mins, &IndexHeader::new(k, w), &cfg))
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
