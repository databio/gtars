wasm:
	wasm-pack build --target web --release gtars-wasm
	jq '.name = "@databio/gtars"' gtars-wasm/pkg/package.json > tmp.json && mv tmp.json gtars-wasm/pkg/package.json

test:
	cargo test --all --workspace -- --nocapture

fmt:
	cargo fmt --all -- --check