FROM rust:alpine3.22 as wasmbuild

# install protoc and other deps
RUN apk update && apk add --no-cache build-base

# required for building the WASM module
RUN curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
RUN rustup target add wasm32-unknown-unknown
RUN cargo install -f wasm-bindgen-cli

# compile and build the web-assembly project
WORKDIR /usr/src/gtars/
COPY ./ /usr/src/gtars/

RUN wasm-pack build gtars-wasm --target web --release

FROM alpine:3.20 as output
# copy the output from the build stage
COPY --from=wasmbuild /usr/src/gtars/gtars-wasm/pkg /output

CMD ["sh", "-c", "cp -r /output/* /host/"]