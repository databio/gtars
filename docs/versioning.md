# Versioning policy for gtars and its crates
Because `gtars` is organized as a [workspace](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html) with multiple crates, we need to ensure that we have a clear versioning policy for the different crates, bindings, and command line tools.

The versioning and tagging scheme was influenced by several other Rust projects we admire/use:

- `polars`: https://github.com/pola-rs/polars/tags
- `bigtools`: https://github.com/jackh726/bigtools/releases
- `noodles`: https://github.com/zaeleus/noodles/tags?after=noodles-fastq-0.19.0

## Versioning scheme
We follow [Semantic Versioning](https://semver.org/) for all crates, bindings, and command line tools, with two tiers:

- **Internal library crates** (e.g. `gtars-core`, `gtars-uniwig`, `gtars-refget`) version independently. Patch versions are bumped as needed and do not need to stay in sync across crates.
- **Wrapper/binding/CLI crates** (`gtars`, `gtars-cli`, `gtars-py`, `gtars-r`, `gtars-js`) sync their **minor version** to the main `gtars` crate at each release. For example, when `gtars` releases `0.7.0`, all wrapper crates with changes are bumped to `0.7.0` as well. This makes it easy to tell which release a wrapper corresponds to.

Patch versions within wrappers can still increment independently between minor releases (e.g. `0.7.0` → `0.7.1`), but the next minor release will re-sync them.

Finally, we also have a specific tagging scheme for releases on GitHub.

## Tagging scheme
- Core Rust Crates: `vX.Y.Z` (e.g., `v0.5.2`). This tag signifies a release of the core library crates to crates.io. It should trigger the cargo publish workflow.
- Python Bindings: `py-X.Y.Z` (e.g., `py-0.3.1`). This tag signifies a release of the Python package to PyPI.
- CLI: `cli-X.Y.Z` (e.g., `cli-1.1.0`). This tag signifies a release of the CLI binaries, which will be attached to a GitHub Release.
- WASM Bindings: `wasm-X.Y.Z` (e.g., `wasm-0.1.5`). This tag signifies a release of the WASM package to npm.

## Automated release workflow

Pushing a `vX.Y.Z` tag (e.g. `v0.7.0`) triggers the `release-all.yml` GitHub Actions workflow, which orchestrates the entire release process automatically:

1. **Rust crates** are published to crates.io in dependency order (leaf crates first, then dependents, then aggregators), with pauses between waves for index propagation. Already-published crates are skipped.
2. **Python bindings**, **CLI binaries**, and **WASM bindings** are built and published in parallel once the Rust crates are done.
3. Component tags (`py-X.Y.Z`, `cli-X.Y.Z`, `wasm-X.Y.Z`) are created automatically after all publishes succeed.

Each individual publish workflow (`rust-publish.yml`, `build-python-bindings.yml`, `build-binaries.yml`, `build-wasm-bindings.yml`) can still be triggered manually via `workflow_dispatch` for one-off publishes.

## An example scenario:

Say we fix a bug in `uniwig`, we will bump its version by a single patch `x.x.1` inside its `Cargo.toml`:
```diff
// gtars-uniwig/Cargo.toml
- version="0.5.0"
+ version="0.5.1"
```
We publish this to `crates.io` using `cargo publish`.

Then we will bump this version accordingly in the `gtars` wrapper crate **and** bump the crates version (since it got a new uniwig)
```diff
// gtars/Cargo.toml
- gtars-uniwig = { version="0.5.0" }
+ gtars-uniwig = { version="0.5.1" }
```
```diff
// gtars/Cargo.toml
- version="0.5.11"
+ version="0.5.12"
```
We will publish this to `crates.io` using `cargo publish`

Finally, because `uniwig` is a tool used in the command line interface, we will bump the version of `gtars-uniwig` in `gtars-cli` to the most recent version with the bug fix similarly to before. Then we will bump the version of `gtars-cli` as a whole (a single patch since its a simple bug fix in `gtars-uniwig`).
```diff
- gtars-uniwig = { version="0.5.0" }
+ gtars-uniwig = { version="0.5.1" }
```
```diff
- version="0.4.1"
+ version="0.4.2"
```
We will publish this to `crates.io` using `cargo publish`