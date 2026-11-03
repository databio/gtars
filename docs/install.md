# 💿 Installing gtars

gtars is available in many forms to support your chosen compute environment and framework. We've extended our core rust implementation into many forms:

- **[Rust Library](#rust-library)**: The core rust implementation can be used as a library in other Rust projects, allowing for deep integration and customization.
- **[Command Line Interface (CLI)](#command-line-interface-cli)**: The core rust implementation can be installed via cargo. This is the most flexible way to use gtars, and is suitable for any environment that supports rust.
- **[Python Package](#python-package)**: A Python package is available via pip, making it easy to integrate gtars into Python-based workflows and applications.
- **[Wasm/JS Package](#wasmjs-package)**: A WebAssembly (Wasm) package is available for JavaScript environments, allowing gtars to be used in web applications and other JS-based projects.
- **[R bindings (package)](#r-package)**: R bindings are available for R users, enabling the use of gtars within R scripts and applications.

Each environment has its own installation instructions, which are detailed below.

## Rust Library
To use gtars as a Rust library, add the following to your `Cargo.toml`. We feature-gate optional dependencies to keep the base install small.

```toml
[dependencies]
gtars = { version = "0.5.0", features = ["uniwig", "tokenizers"] }
```

Then, in your Rust code, you can import and use gtars as follows:

```rust
use gtars::uniwig;
// Your code here
```

For more details, refer to the [gtars Rust documentation](https://docs.rs/gtars).

## Command Line Interface (CLI)
To install the gtars CLI, you need to have Rust and Cargo installed. You can then install gtars using the following command:
```bash
cargo install gtars-cli
```

Similarly, we feature-gate binary dependencies maximize compatibility and minimize install size. You can specify features during installation like so:

```bash
cargo install gtars-cli --features "uniwig tokenizers"
```

Additionally, you can add all features with:
```bash
cargo install gtars-cli --all-features
```
Once installed, you can verify the installation by running:

```bash
gtars --help
```

## Python Package
To install the gtars Python package, you can use pip. The package is available on PyPI and can be installed with the following command:
```bash
pip install gtars
```

You can then import and use gtars in your Python code as follows:

```python
from gtars.tokenizers import Tokenizer

tokenizer = Tokenizer.from_pretrained("databio/atacformer-base-hg38")
```

## Wasm/JS Package
To use gtars in a JavaScript environment, you can install the Wasm package via npm. First, ensure you have Node.js and npm installed. Then, you can install the gtars package with the following command:
```bash
npm install @databio/gtars
```

You can then import and use gtars in your JavaScript code as follows:

```ts
import { Overlapper } from '@databio/gtars';

const universe = [
    ['chr1', 100, 200],
    ['chr1', 150, 250],
    ['chr2', 300, 400],
]

const overlapper = new Overlapper(universe, 'ailist');
console.log(`Using backend: ${overlapper.get_backend()}`);

const query = ['chr1', 180, 220];
const overlaps = overlapper.find(query);

console.log(`Overlaps for ${query}:`, overlaps);
```

## R Package
COMING SOON.