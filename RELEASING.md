# Releasing

gtars is a monorepo of independently-versioned components. Each component is
released **on its own**, when you decide it ships — there is no longer a single
"release everything" event.

## The rule

A release is a GitHub Release whose **tag names the component and its version**:

```
<package>-v<version>
```

`<package>` is the exact workspace package name (the `name` in its `Cargo.toml`,
or the R package). Cutting that release publishes **only** that component, and
the Release itself stands as the commit-pinned, browsable marker for it.

| Tag | Publishes to | Attaches asset? |
|-----|--------------|-----------------|
| `gtars-python-v<ver>` | PyPI | no |
| `gtars-r-v<ver>` | (built in Bioconductor container) | yes — R binary tarball |
| `gtars-cli-v<ver>` | crates.io **and** CLI binaries | yes — CLI binaries |
| `gtars-wasm-v<ver>` | npm (`@databio/gtars`) | no |
| `gtars-node-v<ver>` | npm (`@databio/gtars-node`) | no |
| `gtars-core-v<ver>` | crates.io | no |
| `gtars-refget-v<ver>` | crates.io | no |
| `gtars-v<ver>` | crates.io (meta crate) | no |
| `gtars-<crate>-v<ver>` | crates.io | no |

The `-v` delimiter is load-bearing: `gtars-v0.9.1` (the meta crate) does **not**
collide with `gtars-core-v0.5.7` etc., because the character after `gtars-`
differs (`v` vs `c`). Keep the `-v`.

## How to cut a release

1. Bump the component's version (e.g. `gtars-r/DESCRIPTION` + `gtars-r/src/rust/Cargo.toml`,
   or the crate's `Cargo.toml`). Merge to `master`.
2. Create the GitHub Release (e.g. with `releasegh`), naming the tag
   `<package>-v<version>` and writing release notes.
3. Only that component's job runs; everything else is skipped.

## Bindings vs crates

- **Bindings (python, r, cli, wasm, node)** build from source via in-repo path
  dependencies, so they always compile current monorepo source. There is **no
  dependency cascade** — a binding release is fully self-contained.
- **crates.io crates** depend on each other. Because dependencies are caret
  ranges (`version = "0.5.6"` means `>=0.5.6, <0.6`), a patch/minor bump of a
  dependency is picked up automatically downstream — no action needed. Only a
  **breaking** bump (e.g. `gtars-core` 0.5 → 0.6) requires re-releasing the
  dependents, **in dependency order**: cut `gtars-core-v0.6.0`, wait for
  crates.io indexing, then `gtars-refget-v…`, then `gtars-v…`. A single-crate
  release will fail if its dependencies aren't already on crates.io.

## Manual / emergency publishing

`rust-publish.yml` still supports `workflow_dispatch` with per-target toggles
(`python`, `cli`, `wasm`, `node`, `r`) and a `crates` selector
(`all` / `none` / a single crate). Use this to publish without cutting a
Release, or to republish the entire workspace (`crates: all`) in an emergency.
