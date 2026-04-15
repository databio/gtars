import { Readable } from 'node:stream'

export interface SequenceMetadata {
  name: string
  length: number
  sha512T24U: string
  md5: string
}

export interface CollectionMetadata {
  digest: string
  nSequences: number
  namesDigest: string
  sequencesDigest: string
  lengthsDigest: string
}

export interface StoreStats {
  nSequences: number
  nSequencesLoaded: number
  nCollections: number
  nCollectionsLoaded: number
  storageMode: string
}

export declare class RefgetStore {
  static openLocal(path: string): RefgetStore
  static openRemote(cachePath: string, url: string): RefgetStore
  static inMemory(): RefgetStore
  getSequence(digest: string): string
  getSubstring(digest: string, start: number, end: number): string
  getSequenceByName(collectionDigest: string, name: string): string
  listCollections(page?: number | null, pageSize?: number | null): Array<CollectionMetadata>
  listSequences(): Array<SequenceMetadata>
  getCollectionMetadata(digest: string): CollectionMetadata | null
  getSequenceMetadata(digest: string): SequenceMetadata | null
  stats(): StoreStats
  ensureDecoded(digest: string): void
  clearDecodedCache(): void
  addFasta(fastaPath: string): void
  write(path: string): void
  /** Compare two sequence collections by digest. Returns JSON string. */
  compare(digestA: string, digestB: string): string
  /**
   * Stream a (sub)sequence as a Node `Readable` of ASCII bases.
   *
   * Bytes flow from the underlying Rust store via a background thread,
   * with bounded memory usage (constant w.r.t. sequence length). Use
   * `.pipe(res)` to forward directly to an HTTP response body.
   */
  streamSequence(digest: string, start?: number, end?: number): Readable
}
