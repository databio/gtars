wasm:
	wasm-pack build --target web --release gtars-wasm
	# Ship the hand-written JS layer (RemoteRefgetStore: byte-range + OPFS refget
	# client) alongside the generated bindings, exposed as @databio/gtars/remote.
	cp gtars-wasm/js/remote-refget-store.js gtars-wasm/js/remote-refget-store.d.ts gtars-wasm/pkg/
	jq '.name = "@databio/gtars" | .repository = {"type": "git", "url": "https://github.com/databio/gtars"} | .files += ["remote-refget-store.js", "remote-refget-store.d.ts"] | .exports = {".": {"types": ("./" + .types), "default": ("./" + .main)}, "./remote": {"types": "./remote-refget-store.d.ts", "default": "./remote-refget-store.js"}}' gtars-wasm/pkg/package.json > tmp.json && mv tmp.json gtars-wasm/pkg/package.json

test:
	cargo test --all --workspace -- --nocapture

test-r:
	bash gtars-r/test-r.sh

fmt:
	cargo fmt --all -- --check