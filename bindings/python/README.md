# gtars

This is a python wrapper around the `gtars` crate. It provides an easy interface for using `gtars` in python. It is currently in early development, and as such, it does not have a lot of functionality yet, but new tools are being worked on right now.

## Installation

You can get `gtars` from PyPI:

```bash
pip install gtars
```

## Usage

Import the package, and use the tools:
```python
import gtars as gt

gt.prune_universe(...)
```
## Developer docs

To build for development:

```bash
cd bindings/python
maturin build --release
```

Then install the local wheel that was just built:

```
version=`grep '^version =' Cargo.toml | cut -d '"' -f 2`
pip install --force-reinstall target/wheels/gtars-${version}-cp312-cp312-manylinux_2_38_x86_64.whl
```
