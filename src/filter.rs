#[cfg(not(feature = "server"))]
pub use crate::local_filter::run;
#[cfg(feature = "server")]
pub use crate::remote_filter::run;
