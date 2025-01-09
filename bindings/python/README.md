# gtars

This is a Python package that wraps the `gtars` crate so you can call gtars code from Python.

Documentation for Python bindings is hosted at: https://docs.bedbase.org/gtars/

## Brief instructions

To install the development version, you'll have to build it locally. Build Python bindings like this:

```console
cd bindings/python
maturin build --interpreter 3.11  --release
```

Then install the local wheel that was just built:

```console
gtars_version=`grep '^version =' Cargo.toml | cut -d '"' -f 2`
python_version=$(python --version | awk '{print $2}' | cut -d '.' -f1-2 | tr -d '.')
wheel_path=$(find target/wheels/gtars-${gtars_version}-cp${python_version}-cp${python_version}-*.whl)
pip install --force-reinstall ${wheel_path}
```