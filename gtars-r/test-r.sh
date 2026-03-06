#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

bulker exec databio/nsheff -- R CMD INSTALL --no-multiarch "$SCRIPT_DIR"
bulker exec databio/nsheff -- Rscript -e "testthat::test_dir('${SCRIPT_DIR}/tests')"
