import { describe, it, before } from 'node:test'
import assert from 'node:assert/strict'
import { RefgetStore } from '../index.js'

describe('RefgetStore', () => {
  let store

  before(() => {
    store = RefgetStore.inMemory()
    store.addFasta('../tests/data/fasta/base.fa.gz')
  })

  it('lists sequences after import', () => {
    const seqs = store.listSequences()
    assert.ok(seqs.length > 0)
    assert.ok(seqs[0].sha512T24U)
    assert.ok(seqs[0].name)
  })

  it('lists collections', () => {
    const collections = store.listCollections()
    assert.ok(collections.length > 0)
  })

  it('retrieves a substring', () => {
    const sub = store.getSubstring('iYtREV555dUFKg2_agSJW6suquUyPpMw', 2, 7)
    assert.equal(sub, 'GGGGA')
  })

  it('reports stats', () => {
    const stats = store.stats()
    assert.ok(stats.nSequences > 0)
  })

  it('throws on invalid digest', () => {
    assert.throws(() => store.getSequence('nonexistent'), /not found/i)
  })

  it('throws on invalid path', () => {
    assert.throws(() => RefgetStore.openLocal('/nonexistent/path'))
  })

})
