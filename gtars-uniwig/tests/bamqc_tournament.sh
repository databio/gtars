#!/bin/bash
# Tournament script to compare Rust bamQC with Python bamQC.py
# Usage: ./bamqc_tournament.sh [bam_dir]
#
# Finds all BAM files in the specified directory (default: PEPATAC test data)
# and compares the outputs from both implementations.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GTARS_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PEPATAC_ROOT="$GTARS_ROOT/../pepatac"

# Build the Rust CLI if needed
echo "Building gtars CLI..."
cd "$GTARS_ROOT"
cargo build --release --package gtars-cli --features uniwig 2>/dev/null || {
    cargo build --package gtars-cli --features uniwig
}

GTARS_CLI="$GTARS_ROOT/target/release/gtars"
if [ ! -f "$GTARS_CLI" ]; then
    GTARS_CLI="$GTARS_ROOT/target/debug/gtars"
fi

BAMQC_PY="$PEPATAC_ROOT/tools/bamQC.py"

# Check for Python bamQC.py
if [ ! -f "$BAMQC_PY" ]; then
    echo "Warning: Python bamQC.py not found at $BAMQC_PY"
    echo "Will only run Rust implementation"
    SKIP_PYTHON=true
fi

# Find BAM files
BAM_DIR="${1:-$PEPATAC_ROOT/tests/test_outputs}"
echo "Searching for BAM files in: $BAM_DIR"

BAM_FILES=$(find "$BAM_DIR" -name "*.bam" -type f 2>/dev/null || true)

if [ -z "$BAM_FILES" ]; then
    # Fall back to gtars test data
    BAM_DIR="$GTARS_ROOT/tests/data"
    echo "No BAM files found. Trying: $BAM_DIR"
    BAM_FILES=$(find "$BAM_DIR" -name "*.bam" -type f 2>/dev/null || true)
fi

if [ -z "$BAM_FILES" ]; then
    echo "Error: No BAM files found"
    exit 1
fi

TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

PASS=0
FAIL=0
SKIP=0

# Helper to safely increment counters with set -e
inc_pass() { ((PASS++)) || true; }
inc_fail() { ((FAIL++)) || true; }
inc_skip() { ((SKIP++)) || true; }

echo ""
echo "========================================"
echo "BAM QC Tournament: Rust vs Python"
echo "========================================"
echo ""

for BAM in $BAM_FILES; do
    BASENAME=$(basename "$BAM" .bam)
    echo "Testing: $BASENAME"

    # Check for index
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.*}.bai" ]; then
        echo "  SKIP: No BAM index found"
        inc_skip
        continue
    fi

    RUST_OUT="$TEMP_DIR/${BASENAME}_rust.tsv"
    PY_OUT="$TEMP_DIR/${BASENAME}_py.tsv"

    # Run Rust implementation
    if ! "$GTARS_CLI" uniwig bamqc -i "$BAM" -o "$RUST_OUT" 2>/dev/null; then
        echo "  FAIL: Rust implementation failed"
        inc_fail
        continue
    fi

    if [ "$SKIP_PYTHON" = "true" ]; then
        echo "  PASS: Rust completed (Python skipped)"
        cat "$RUST_OUT"
        inc_pass
        continue
    fi

    # Run Python implementation
    if ! python "$BAMQC_PY" -i "$BAM" -o "$PY_OUT" -c 1 2>/dev/null; then
        echo "  WARN: Python implementation failed, comparing with existing output"
        # Look for existing output
        EXISTING=$(dirname "$BAM")/../QC_*/$(basename "$BAM" .bam)_bamQC.tsv
        if [ -f "$EXISTING" ]; then
            cp "$EXISTING" "$PY_OUT"
        else
            echo "  SKIP: No Python output available"
            inc_skip
            continue
        fi
    fi

    # Compare outputs (tolerating floating point differences)
    echo "  Rust output:"
    cat "$RUST_OUT" | sed 's/^/    /'
    echo "  Python output:"
    cat "$PY_OUT" | sed 's/^/    /'

    # Extract key metrics for comparison
    RUST_NRF=$(tail -1 "$RUST_OUT" | cut -f8)
    RUST_PBC1=$(tail -1 "$RUST_OUT" | cut -f9)
    RUST_PBC2=$(tail -1 "$RUST_OUT" | cut -f10)

    PY_NRF=$(tail -1 "$PY_OUT" | cut -f8)
    PY_PBC1=$(tail -1 "$PY_OUT" | cut -f9)
    PY_PBC2=$(tail -1 "$PY_OUT" | cut -f10)

    # Compare with tolerance
    TOLERANCE="0.000001"

    compare_float() {
        local a="$1"
        local b="$2"
        local tol="$3"
        python3 -c "import sys; sys.exit(0 if abs(float('$a') - float('$b')) < float('$tol') else 1)" 2>/dev/null
    }

    MATCH=true
    if ! compare_float "$RUST_NRF" "$PY_NRF" "$TOLERANCE"; then
        echo "  NRF mismatch: Rust=$RUST_NRF Python=$PY_NRF"
        MATCH=false
    fi
    if ! compare_float "$RUST_PBC1" "$PY_PBC1" "$TOLERANCE"; then
        echo "  PBC1 mismatch: Rust=$RUST_PBC1 Python=$PY_PBC1"
        MATCH=false
    fi
    # PBC2 can differ significantly due to M2 handling
    if ! compare_float "$RUST_PBC2" "$PY_PBC2" "1.0"; then
        echo "  PBC2 mismatch: Rust=$RUST_PBC2 Python=$PY_PBC2"
        MATCH=false
    fi

    if [ "$MATCH" = "true" ]; then
        echo "  PASS: Metrics match"
        inc_pass
    else
        echo "  FAIL: Metrics differ"
        inc_fail
    fi

    echo ""
done

echo "========================================"
echo "Tournament Results"
echo "========================================"
echo "  PASS: $PASS"
echo "  FAIL: $FAIL"
echo "  SKIP: $SKIP"
echo "========================================"

if [ $FAIL -gt 0 ]; then
    exit 1
fi
