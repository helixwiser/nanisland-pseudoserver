#!/bin/bash
# Run the full Aisland pipeline end-to-end.
# Usage: bash pipeline/run_all.sh [config/params.yaml]

set -euo pipefail
cd "$(dirname "$0")/.."

CONFIG="${1:-config/params.yaml}"

echo "=== Aisland Pipeline ==="
echo "Config: $CONFIG"
echo ""

echo "--- Step 1: Parse MAF files ---"
python pipeline/01_parse_maf.py "$CONFIG"
echo ""

echo "--- Step 2: Build count matrices ---"
python pipeline/02_build_matrix.py "$CONFIG"
echo ""

echo "--- Step 3: Differential analysis ---"
python pipeline/03_differential.py "$CONFIG"
echo ""

echo "--- Step 4: Fetch FASTA (optional) ---"
python pipeline/04_fetch_fasta.py "$CONFIG"
echo ""

echo "=== Pipeline complete ==="
