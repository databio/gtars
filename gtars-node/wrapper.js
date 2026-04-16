/* gtars-node public entry point.
 *
 * Subclasses the auto-generated napi `RefgetStore` to:
 *   - expose `streamSequence(...)` returning a Node `stream.Readable`
 *     bridged from the underlying callback-based `streamSequenceInternal`.
 *   - coerce `number | bigint` arguments at the JS boundary so callers
 *     don't have to write `BigInt(...)` for chromosome-sized values that
 *     fit comfortably in a JS Number.
 *
 * The generated factory methods (`openLocal`, `openRemote`, `inMemory`)
 * return raw native instances; we rewrap them via prototype reassignment
 * so `instanceof RefgetStore` works for callers.
 */

const { Readable } = require('node:stream')
const native = require('./index.js')

const NativeRefgetStore = native.RefgetStore

// Coerce number|bigint -> bigint (for napi BigInt boundary). null/undefined
// pass through so callers can still omit optional bounds.
function toBigInt(v) {
  if (v == null) return null
  return typeof v === 'bigint' ? v : BigInt(v)
}

class RefgetStore extends NativeRefgetStore {
  static openLocal(path) {
    return Object.setPrototypeOf(NativeRefgetStore.openLocal(path), RefgetStore.prototype)
  }
  static openRemote(cachePath, url) {
    return Object.setPrototypeOf(
      NativeRefgetStore.openRemote(cachePath, url),
      RefgetStore.prototype,
    )
  }
  static inMemory() {
    return Object.setPrototypeOf(NativeRefgetStore.inMemory(), RefgetStore.prototype)
  }

  // Coerce number|bigint bounds for the BigInt-typed napi getSubstring.
  getSubstring(digest, start, end) {
    return super.getSubstring(digest, BigInt(start), BigInt(end))
  }

  streamSequence(digest, start, end) {
    const stream = new Readable({
      read() {
        // Push-driven; backpressure is applied on the Rust side via the
        // blocking threadsafe-function call.
      },
    })

    const onChunk = (buf) => {
      stream.push(buf)
    }
    const onError = (msg) => {
      stream.destroy(msg instanceof Error ? msg : new Error(String(msg)))
    }
    const onEnd = () => {
      stream.push(null)
    }

    try {
      super.streamSequenceInternal(
        digest,
        toBigInt(start),
        toBigInt(end),
        onChunk,
        onError,
        onEnd,
      )
    } catch (err) {
      // Synchronous setup errors (e.g. sequence not found) — propagate as
      // a stream 'error' event on the next tick so callers who've attached
      // listeners after the call still receive it. Always wrap in Error
      // so the surface matches the async onError path.
      const wrapped = err instanceof Error ? err : new Error(String(err))
      queueMicrotask(() => stream.destroy(wrapped))
    }

    return stream
  }
}

module.exports = {
  ...native,
  RefgetStore,
}
