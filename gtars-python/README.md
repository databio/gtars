# gtars

This is a Python package that wraps the `gtars` crate so you can call gtars code from Python.

Documentation for Python bindings is hosted at: https://docs.bedbase.org/gtars/

## Brief instructions

To install the development version, you'll have to build it locally. The easy way to build/install is to do this, which reads the configuration in `pyproject.toml` for easy building/installing:

```
python -m pip install -e .
```

Or, you can have more control over it in two steps using `maturin` directly.
Build Python bindings like this:

```console
cd gtars-python
python_version=$(python --version | awk '{print $2}' | cut -d '.' -f1-2 )
maturin build --interpreter $python_version  --release
```


Then install the local wheel that was just built:

```console
gtars_version=`grep '^version =' Cargo.toml | cut -d '"' -f 2`
python_version_nodot=$(python --version | awk '{print $2}' | cut -d '.' -f1-2 | tr -d '.')
ll ../target/wheels/gtars*
wheel_path=$(find ../target/wheels/gtars-${gtars_version}-cp${python_version_nodot}-cp${python_version_nodot}-manylinux*.whl)
echo $wheel_path
pip install --force-reinstall ${wheel_path}
```


## Importing into python

Once installed, you can import and use the package in Python. For example:

```
from gtars import refget
sc = refget.digest_fasta("../tests/data/fasta/base.fa")
sc2 = refget.load_fasta("../tests/data/fasta/base.fa")
```
