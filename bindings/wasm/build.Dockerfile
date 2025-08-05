FROM rustlang/rust:nightly-alpine3.20 as wasmbuild

# install protoc and other deps
RUN apk update && apk add --no-cache curl musl-dev clang openssl-dev build-base

# https://github.com/rust-lang/rust/issues/115450#issuecomment-1717228111a
ENV OPENSSL_DIR=/usr
ENV RUSTFLAGS="--Z wasm_c_abi=spec"

# # required for building the WASM module
RUN curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

# # compile and build the web-assembly project
WORKDIR /usr/src/gtars/bindings/wasm
COPY bindings/wasm/ /usr/src/gtars/bindings/wasm/
COPY gtars/ /usr/src/gtars/gtars/

RUN wasm-pack build --target web --release

# CMD ["sleep", "infinity"]