# gtars Module Vision

## Conceptual Organization

The gtars modules follow a clear separation of concerns:

### 1. gtars-overlaprs
**Core overlap computation infrastructure**
- All interval overlap algorithms live here
- Single source of truth for overlap operations
- Other modules should never reimplement overlap logic

### 2. gtars-scoring
**Wraps overlaprs to produce X-by-peak matrices**
- Handles all matrix generation use cases:
  - Cell-by-peak (single-cell)
  - Sample-by-peak (bulk)
  - Pseudobulk-by-peak (aggregated)
- Focuses on matrix data structures and I/O formats

### 3. gtars-tokenizers
**Wraps overlaprs to produce tokens for ML**
- Converts genomic regions to vocabulary tokens
- Designed for transformer models and ML pipelines
- Manages token vocabularies and encoding strategies

## Key Principle

**All overlap computation flows through gtars-overlaprs**. Higher-level modules add domain-specific value:
- scoring adds matrix generation
- tokenizers adds ML tokenization
- Both delegate overlap computation to overlaprs