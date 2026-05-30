# @databio/gtars-node

Native Node.js bindings for [gtars](https://github.com/databio/gtars) (genomic tools and region sets). Provides a `RefgetStore` class for working with GA4GH refget sequence collections — loading FASTA files, retrieving sequences by digest, and comparing collections.

Built with [NAPI-RS](https://napi.rs). The main package auto-detects your platform and loads the correct native binary.

## Supported platforms

- Linux x64 (glibc)
- macOS ARM64 (Apple Silicon)

## Quick start

```javascript
const { RefgetStore } = require("@databio/gtars-node")
const store = RefgetStore.openRemote("/tmp/rgcache", "https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle/")
console.log(store.stats())                // { nSequences: 580742, nCollections: 205, ... }
console.log(store.listCollections(0, 5))  // first 5 collections
const seq = store.listCollections(0, 1)[0]
console.log(store.getCollectionMetadata(seq.digest))
```

## Installation

```bash
npm install @databio/gtars-node
```

## Usage

```javascript
const { RefgetStore } = require('@databio/gtars-node')

// Create an in-memory store and load a FASTA file
const store = RefgetStore.inMemory()
store.addFasta('/path/to/genome.fa.gz')

// Inspect what was loaded
store.stats()            // { nSequences, nCollections, storageMode, ... }
store.listSequences()    // [{ name, length, sha512T24U, md5 }, ...]
store.listCollections()  // [{ digest, nSequences, namesDigest, ... }, ...]

// Retrieve a sequence by its digest
store.getSequence('iYtREV555dUFKg2_agSJW6suquUyPpMw')

// Retrieve a substring (0-indexed, exclusive end)
store.getSubstring('iYtREV555dUFKg2_agSJW6suquUyPpMw', 2, 7)

// Get metadata for a specific sequence or collection
store.getSequenceMetadata(digest)
store.getCollectionMetadata(digest)

// Compare two sequence collections (returns a JSON string)
const result = JSON.parse(store.compare(digestA, digestB))

// Persist to disk, then reload
store.write('/path/to/store')
const reloaded = RefgetStore.openLocal('/path/to/store')

// Open a remote store (fetches rgstore.json from the URL, caches locally)
const remote = RefgetStore.openRemote('/local/cache', 'https://example.com/path/to/store')
remote.listCollections()
```

For ESM modules, use `createRequire`:

```javascript
import { createRequire } from 'node:module'
const require = createRequire(import.meta.url)
const { RefgetStore } = require('@databio/gtars-node')
```

## Opening a remote store

`openRemote` lets you load a refget store hosted on any static file server. The remote URL must serve a directory created by `store.write()`, with `rgstore.json` at the root:

```
https://example.com/my-store/
  rgstore.json        # store manifest
  sequences/...       # sequence data files
```

The first argument is a local cache directory where fetched data is stored:

```javascript
const store = RefgetStore.openRemote('/tmp/rgcache', 'https://example.com/my-store')
store.listCollections()
store.getSequence('iYtREV555dUFKg2_agSJW6suquUyPpMw')
```

To publish a store for remote access, build it locally and upload the output directory to any static host (S3, GitHub Pages, nginx, etc.):

```javascript
const store = RefgetStore.inMemory()
store.addFasta('/path/to/genome.fa.gz')
store.write('/path/to/output')   // upload this directory
```

## API reference

### Factory methods

| Method | Description |
|--------|-------------|
| `RefgetStore.inMemory()` | Create an empty in-memory store. Load data with `addFasta()`. |
| `RefgetStore.openLocal(path)` | Open a store previously saved with `write()`. |
| `RefgetStore.openRemote(cachePath, url)` | Fetch a store from a remote URL with local caching. |

### Sequence retrieval

Sequences are lazy-loaded automatically — when you call `getSequence` or `getSubstring` on a store opened from disk or remote, the sequence data is fetched on first access and cached.

| Method | Returns |
|--------|---------|
| `getSequence(digest)` | Full sequence string for the given digest |
| `getSubstring(digest, start, end)` | Substring (0-indexed, exclusive end) |
| `getSequenceByName(collectionDigest, name)` | Sequence looked up by collection digest and sequence name |

### Listing and metadata

| Method | Returns |
|--------|---------|
| `listSequences()` | `SequenceMetadata[]` — name, length, sha512T24U, md5 for each sequence |
| `listCollections(page?, pageSize?)` | `CollectionMetadata[]` — supports optional pagination |
| `getSequenceMetadata(digest)` | `SequenceMetadata \| null` |
| `getCollectionMetadata(digest)` | `CollectionMetadata \| null` |
| `stats()` | `StoreStats` — counts and storage mode |

**Note:** `listSequences()` returns all sequences at once. For large stores (500k+ sequences), use `stats()` to check the count first and prefer `getSequenceMetadata(digest)` for individual lookups.

### Streaming

`streamSequence` returns a Node `stream.Readable` that emits ASCII sequence bytes as they are decoded from the underlying store — no full-sequence buffering on the JS side. This is the recommended path for HTTP refget servers that want to pipe sequence bytes directly to the response body with bounded memory usage.

| Method | Returns |
|--------|---------|
| `streamSequence(digest, start?, end?)` | `stream.Readable` yielding ASCII bases |

Example — an Express handler for `GET /sequence/:digest`:

```javascript
const { RefgetStore } = require('@databio/gtars-node')
const express = require('express')

const store = RefgetStore.openLocal('/path/to/store')
const app = express()

app.get('/sequence/:digest', (req, res) => {
  const stream = store.streamSequence(req.params.digest)
  stream.on('error', (err) => {
    res.status(/not found/i.test(err.message) ? 404 : 500).send(err.message)
  })
  res.setHeader('Content-Type', 'text/vnd.ga4gh.refget.v2.0.0+plain')
  stream.pipe(res)
})
```

### Mutation and persistence

| Method | Description |
|--------|-------------|
| `addFasta(fastaPath)` | Load sequences from a FASTA file (gzipped or plain) |
| `write(path)` | Save the store to a directory on disk |
| `compare(digestA, digestB)` | Compare two collections; returns a JSON string |
| `ensureDecoded(digest)` | Decode a sequence into the cache |
| `clearDecodedCache()` | Free memory used by decoded sequences |

## Development

```bash
# Build the native module (requires Rust toolchain)
npm run build

# Run tests
npm test
```

## License

MIT
