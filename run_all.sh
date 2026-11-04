#!/usr/bin/env bash
# Runs every check. Usage: ./run_all.sh [--with-data]
set -euo pipefail
cd "$(dirname "$0")"
OUT=output; mkdir -p $OUT
echo "== 1. Arb certificate (REQUIRED by the proof)";      python3 verifier/certificate.py      | tee $OUT/certificate.txt
echo "== 2. Symbolic algebra checks";                       python3 verifier/sym_check.py        | tee $OUT/sym_check.txt
echo "== 3. Dilation response (Gate B) numeric check";      python3 verifier/gateB_numeric.py    | tee $OUT/gateB_numeric.txt
echo "== 4. Response scaling study";                        python3 verifier/response_scaling.py | tee $OUT/response_scaling.txt
echo "== 5. Taylor coefficients via Cauchy integrals (~1 min)"; python3 verifier/taylor_check.py | tee $OUT/taylor_check.txt
if [ "${1:-}" = "--with-data" ]; then
  LAMBDA=$(./fetch_data.sh | tail -1)
  echo "== 6. Curvature windows on real xi data";          python3 -I verifier/xi_data.py "$LAMBDA" | tee $OUT/xi_data.txt
  echo "== 7. Whitened perturbation 3-pt vs 4-pt";         (cd verifier && python3 whitened_experiment.py "../$LAMBDA") | tee $OUT/whitened_experiment.txt
fi
echo "Done. Compare $OUT/ with expected_output/."
