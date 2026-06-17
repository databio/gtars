// Type declarations for the RemoteRefgetStore byte-range + OPFS refget client.

/** Metadata for one sequence in a refgetstore index. */
export interface SeqMeta {
  /** Sequence name (e.g. "chr1"). */
  name: string;
  /** Sequence length in bases. */
  length: number;
  /** Alphabet ("dna2bit", "dna3bit", "dnaio"/"dnaiupac", or a plain alphabet). */
  alphabet: string;
  /** sha512t24u digest (the GA4GH/VRS sequence accession). */
  sha512t24u: string;
  /** MD5 digest (may be empty). */
  md5: string;
}

/** The two wasm primitives RemoteRefgetStore needs (import them from @databio/gtars). */
export interface RefgetWasm {
  encodedByteRange(start: number, end: number, alphabet: string): number[];
  decodeEncodedRange(
    bytes: Uint8Array,
    byteOffset: number,
    start: number,
    end: number,
    alphabet: string,
  ): string;
}

export interface RemoteRefgetStoreOptions {
  /** Base URL of the served refgetstore directory. */
  baseUrl: string;
  /** The wasm module (or just the two required functions). */
  wasm: RefgetWasm;
  /** OPFS subdirectory for the byte/window cache (default "refget"). */
  opfsDir?: string;
  /** Base padding applied on each side of a range fetch to coalesce reads. */
  windowPadBases?: number;
  /** Override the `.seq` path template (default "sequences/%s2/%s.seq"). */
  seqTemplate?: string;
  /** Override the collection path template (default "collections/%s.rgsi"). */
  collTemplate?: string;
  /** Optional logger. */
  log?: (text: string, isError?: boolean) => void;
}

/**
 * Byte-range + OPFS-caching reader for a served gtars refgetstore. Reads regions
 * via HTTP `Range:` requests decoded by the wasm primitives, caching fetched
 * windows in OPFS. Never downloads a whole chromosome unless `prefetchSequence`
 * is called.
 */
export class RemoteRefgetStore {
  constructor(opts: RemoteRefgetStoreOptions);

  readonly base: string;
  byName: Map<string, SeqMeta>;
  byDigest: Map<string, SeqMeta>;
  collectionDigest?: string;

  /** Fetch the store manifest and adopt its path templates. */
  prepare(): Promise<any>;
  /** Bind a collection by digest; loads its sequence index. Returns seq count. */
  openCollection(collectionDigest: string): Promise<number>;
  /** Load a flat store's top-level sequence index. Returns seq count. */
  loadSequenceIndex(indexName?: string): Promise<number>;
  /** Register a pre-parsed sequence list. Returns seq count. */
  setSequences(seqMetas: SeqMeta[]): number;
  /** Resolve a sequence name to its metadata (exact match), or null. */
  resolve(name: string): SeqMeta | null;

  /** Read bases [start, end) via a byte-range fetch + decode (OPFS-cached). */
  getSubstring(name: string, start: number, end: number): Promise<string>;
  /** Download the whole encoded sequence once (OPFS-cached). Returns bytes downloaded. */
  prefetchSequence(name: string): Promise<number>;
  /**
   * Return the whole encoded `.seq` bytes for `name` (download once, OPFS-cached).
   * The bridge to the synchronous wasm VRS compute path (feed into
   * `RefgetStore.add_encoded_sequence`).
   */
  getSequenceBytes(name: string): Promise<{ bytes: Uint8Array; cached: boolean; meta: SeqMeta }>;
}

/** Parse a sequences.rgsi index into SeqMeta records. */
export function parseRgsi(text: string): SeqMeta[];
