//! Regenerate the .pyi type stubs from the compiled module. Run via scripts/gen-py-stubs.sh.

use std::path::{Path, PathBuf};

fn main() {
    let lib = PathBuf::from(
        std::env::args()
            .nth(1)
            .expect("usage: gen_stubs <path-to-compiled-cdylib>"),
    );
    let module = pyo3_introspection::introspect_cdylib(&lib, "_deacon")
        .expect("failed to introspect cdylib");

    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("python/deacon");
    for (name, content) in pyo3_introspection::module_stub_files(&module) {
        // root module comes back as __init__.pyi, ship it as the submodule stub _deacon.pyi
        let name = if name == Path::new("__init__.pyi") {
            PathBuf::from("_deacon.pyi")
        } else {
            name
        };
        let dest = out_dir.join(&name);
        std::fs::create_dir_all(dest.parent().unwrap()).unwrap();
        std::fs::write(&dest, content).unwrap();
        println!("wrote {}", dest.display());
    }
}
