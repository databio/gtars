import { Readable } from 'node:stream'

// Re-export the auto-generated types (interfaces, plus the native class
// type) so consumers see a single surface from the entry point.
export * from './index.js'

import { RefgetStore as NativeRefgetStore } from './index.js'

export declare class RefgetStore extends NativeRefgetStore {
  static openLocal(path: string): RefgetStore
  static openRemote(cachePath: string, url: string): RefgetStore
  static inMemory(): RefgetStore

  /** Accepts plain `number` or `bigint` for ergonomic call sites. */
  getSubstring(digest: string, start: number | bigint, end: number | bigint): string

  /**
   * Stream a (sub)sequence as a Node `Readable` of ASCII bases.
   *
   * Bytes flow from the underlying Rust store via a background thread
   * with bounded memory usage (constant w.r.t. sequence length). Use
   * `.pipe(res)` to forward directly to an HTTP response body, or
   * `stream.destroy()` to abort early — the worker thread will exit
   * promptly via an internal abort flag.
   */
  streamSequence(
    digest: string,
    start?: number | bigint | null,
    end?: number | bigint | null,
  ): Readable
}
