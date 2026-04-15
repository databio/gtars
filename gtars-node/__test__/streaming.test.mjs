import { describe, it, before } from 'node:test'
import assert from 'node:assert/strict'
import { Writable } from 'node:stream'
import { pipeline } from 'node:stream/promises'
import { RefgetStore } from '../wrapper.js'

const FULL_DIGEST = 'iYtREV555dUFKg2_agSJW6suquUyPpMw'
const FULL_SEQ = 'TTGGGGAA'

async function collect(stream) {
  const chunks = []
  for await (const chunk of stream) {
    chunks.push(chunk)
  }
  return Buffer.concat(chunks)
}

describe('RefgetStore streaming', () => {
  let store

  before(() => {
    store = RefgetStore.inMemory()
    store.addFasta('../tests/data/fasta/base.fa.gz')
  })

  it('streams the full sequence and matches getSequence', async () => {
    const stream = store.streamSequence(FULL_DIGEST)
    const buf = await collect(stream)
    assert.equal(buf.toString('ascii'), FULL_SEQ)
    assert.equal(buf.toString('ascii'), store.getSequence(FULL_DIGEST))
  })

  it('streams a substring matching getSubstring', async () => {
    const stream = store.streamSequence(FULL_DIGEST, 2, 7)
    const buf = await collect(stream)
    assert.equal(buf.toString('ascii'), store.getSubstring(FULL_DIGEST, 2, 7))
  })

  it('emits an error event for unknown digest', async () => {
    const stream = store.streamSequence('nonexistent')
    await assert.rejects(collect(stream), /not found/i)
  })

  it('emits an error event for out-of-bounds range', async () => {
    const stream = store.streamSequence(FULL_DIGEST, 100, 200)
    await assert.rejects(collect(stream), /range|length|bounds|stream/i)
  })

  it('supports piping to a slow writable (backpressure smoke test)', async () => {
    const stream = store.streamSequence(FULL_DIGEST)
    const received = []
    const slow = new Writable({
      highWaterMark: 2,
      write(chunk, _enc, cb) {
        received.push(chunk)
        setImmediate(cb)
      },
    })
    await pipeline(stream, slow)
    const joined = Buffer.concat(received).toString('ascii')
    assert.equal(joined, FULL_SEQ)
  })
})
