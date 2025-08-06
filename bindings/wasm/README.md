# gtars wasm
This directory contains the WebAssembly (Wasm) bindings for the gtars crate. It enables you to use gtars functionalities in a web environment, allowing for efficient genomic data processing directly in the browser (i.e. no server required!!!)

## Building the Wasm Bindings
To build the Wasm bindings, you need to have Rust and the wasm32 target installed. We use the `wasm-pack` tool to build the project. You can install it via:

```bash
cargo install wasm-pack
```
Then, you can build the Wasm bindings with the following command:

```bash
wasm-pack build --target web --release
```

Ths will generate the Wasm package in the `pkg` directory, which can be used in your web applications:

```bash
npm install path/to/gtars/bindings/wasm/pkg
```
