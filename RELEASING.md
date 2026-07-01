# Releasing

Tagging triggers a GitHub release with CLI binaries, does crates.io upload, and builds PyPI wheels.

1. Bump `[workspace.package] version` in `Cargo.toml`; update changelog.
2. `cargo fmt --check && cargo clippy && cargo test`.
3. Regenerate Python stubs and confirm none are stale:
   `./scripts/gen-py-stubs.sh && git diff --exit-code deacon-py/python/deacon/_deacon.pyi`.
4. Merge to `main`, then `git tag 0.17.0 && git push origin 0.17.0`.
5. Approve the `pypi` environment if prompted; publish the draft GitHub release.
6. Stable only: merge the bioconda autobump PR.

**Python-only prerelease:** set an explicit `version = "0.17.0rc1"` in `deacon-py/Cargo.toml` (overriding the workspace version) and run `release-pypi.yaml` via workflow_dispatch. No tag, so CLI/crate are untouched. Installs only with `pip install --pre`. Revert to `version.workspace = true` afterwards.

**One-time setup:** PyPI trusted/pending publisher (project `deacon`, repo `bede/deacon`, workflow `release-pypi.yaml`, env `pypi`); GitHub env `pypi`; secret `CARGO_REGISTRY_TOKEN`.
