#!/usr/bin/env bash
set -euo pipefail

pairs=(
  "17780 17792"
  "17780 17793"
  "17781 17792"
  "17781 17793"
)

for pair in "${pairs[@]}"; do
  set -- $pair
  r1="$1"
  r2="$2"

  echo "============================================================"
  echo "Running QE pair ${r1}/${r2}"
  echo "============================================================"

  perl -0pi -e "s/^run1:.*/run1: ${r1}/m; s/^run2:.*/run2: ${r2}/m" \
    configs/qe_analysis_kin0.yaml configs/qe_analysis_kin1.yaml

  ./run_qe_analysis.sh --build --config configs/qe_analysis_kin0.yaml
  ./run_qe_analysis.sh --config configs/qe_analysis_kin1.yaml
done
