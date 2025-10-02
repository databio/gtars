# gtars-js
`gtars-js` is the wasm bindings for our [`gtars` Rust library](https://github.com/databio/gtars). It provides a JavaScript interface to the core functionality of `gtars`, allowing you to use it in web applications or other JavaScript environments.

## Installation
You can install `gtars-js` via npm:

```bash
npm install @databio/gtars-js
```

## Usage
Here's a simple example of how to use the `Overlapper` class from `gtars-js`:
```ts
import { Overlapper } from '@databio/gtars-js';

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