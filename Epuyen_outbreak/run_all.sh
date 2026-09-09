#!/bin/bash
# Master pipeline for the Epuyen Andes virus outbreak replication analysis.
# Runs steps 01-07 in sequence. Idempotent - safe to re-run from scratch.
set -euo pipefail

CODON_PY="${CODON_PY:-python3}"
HANTA_R="${HANTA_R:-Rscript}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

echo "=== $(date) : Starting pipeline ==="

echo "--- 01: near-cognate tRNA mapping ---"
"$CODON_PY" pipeline/01_trna_near_cognate.py

echo "--- 02: codon usage metrics (CAI/ENC) ---"
"$CODON_PY" pipeline/02_codon_usage_metrics.py

echo "--- 03: co-AI ---"
"$CODON_PY" pipeline/03_coai.py

echo "--- 04: folding energy (mFE) ---"
"$CODON_PY" pipeline/04_folding_energy.py

echo "--- 05: synonymous randomization + t-test ---"
"$CODON_PY" pipeline/05_randomization.py

echo "--- 06: pn/ps polymorphism + mutation classification ---"
"$CODON_PY" pipeline/06_pnps_mutations.py

echo "--- 07: R statistics + figures ---"
"$HANTA_R" pipeline/07_stats_figures.R

echo "=== $(date) : Pipeline complete ==="
