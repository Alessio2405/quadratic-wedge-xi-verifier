#!/usr/bin/env bash
# Downloads the certified lambda_n intervals (n <= 12016) published by Liu & Mao.
# Needed only for the numerical illustration (Section 9), not for the proof.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data
if [ ! -d data/riemann-xi-infinite-log-concavity ]; then
  git clone --depth 1 https://github.com/YXLiu97/riemann-xi-infinite-log-concavity.git data/riemann-xi-infinite-log-concavity
fi
echo "data/riemann-xi-infinite-log-concavity/supplement/downloads/riemann_xi_infinite_log_concavity/finite_certificate_current_review_lambda.json"
