//! This module provides filtering functionality for processing data.
//! Basic switching of which `run` to use based on whether the executable is compiled with the
//! `server` feature or not.
//!
#[cfg(not(feature = "server"))]
pub use crate::local_filter::run;
#[cfg(feature = "server")]
pub use crate::remote_filter::run;
