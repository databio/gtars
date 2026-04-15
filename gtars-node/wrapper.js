/* gtars-node public entry point.
 *
 * Imports the auto-generated napi loader (`index.js`) and adds a
 * `streamSequence(digest, start?, end?)` method that returns a Node
 * `stream.Readable`, bridging the underlying `streamSequenceInternal`
 * callback-based Rust API.
 */

const { Readable } = require('node:stream')
const native = require('./index.js')

const NativeRefgetStore = native.RefgetStore

// Patch streamSequence onto the native class so instances returned from
// factory methods (openLocal/openRemote/inMemory) already have it.
NativeRefgetStore.prototype.streamSequence = function streamSequence(
  digest,
  start,
  end,
) {
  const stream = new Readable({
    read() {
      // Push-driven; nothing to do here. Backpressure is applied on the
      // Rust side via a blocking threadsafe-function call.
    },
  })

  const onChunk = (buf) => {
    stream.push(buf)
  }
  const onError = (msg) => {
    stream.destroy(new Error(msg))
  }
  const onEnd = () => {
    stream.push(null)
  }

  try {
    this.streamSequenceInternal(
      digest,
      start ?? null,
      end ?? null,
      onChunk,
      onError,
      onEnd,
    )
  } catch (err) {
    // Synchronous setup errors (e.g. sequence not found) — propagate as
    // a stream 'error' event on the next tick so callers who've attached
    // listeners after the call still receive it.
    queueMicrotask(() => stream.destroy(err))
  }

  return stream
}

module.exports = {
  RefgetStore: NativeRefgetStore,
}
