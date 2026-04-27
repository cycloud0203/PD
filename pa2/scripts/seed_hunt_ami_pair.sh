#!/usr/bin/env bash
# Reproducible seed sweep for ami33 and ami49 only.
# Set FP_SEED in the environment for each trial (see src/floorplanner.cpp).
#
# Usage:
#   bash scripts/seed_hunt_ami_pair.sh
# Optional:
#   ALPHA=0.5 SEED_START=0 TRIALS=20 FP_BIN=./bin/fp
#   INPUT_DIR=input_pa2  OUT_ROOT=output/seed_ami  LOG_ROOT=logs/seed_ami

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

FP_BIN="${FP_BIN:-./bin/fp}"
EVAL="${EVAL:-./evaluator/evaluator.sh}"
INPUT_DIR="${INPUT_DIR:-input_pa2}"
OUT_ROOT="${OUT_ROOT:-output/seed_ami}"
LOG_ROOT="${LOG_ROOT:-logs/seed_ami}"

ALPHA="${ALPHA:-0.5}"
SEED_START="${SEED_START:-0}"
TRIALS="${TRIALS:-20}"
CASES=(ami33 ami49)

if [[ ! -x "$FP_BIN" ]]; then
  echo "error: missing $FP_BIN (run: make)" 1>&2
  exit 1
fi
if [[ ! -x "$EVAL" ]]; then
  echo "error: missing $EVAL" 1>&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_ROOT}_${STAMP}"
LOG_DIR="${LOG_ROOT}_${STAMP}"
mkdir -p "$RUN_DIR" "$LOG_DIR"

CSV_PATH="${LOG_DIR}/results.csv"
{
  echo "case,fp_seed,alpha,fp_exit,fp_sec,rpt,fp_log,eval_exit,eval_sec,normalized_cost,public_score_0_5,eval_log"
} > "$CSV_PATH"

for CASE in "${CASES[@]}"; do
  for ((i = 0; i < TRIALS; i++)); do
    SEED=$((SEED_START + i))
    rpt_path="${RUN_DIR}/${CASE}_seed${SEED}.rpt"
    fp_log="${LOG_DIR}/fp_${CASE}_seed${SEED}.log"
    eval_log="${LOG_DIR}/eval_${CASE}_seed${SEED}.log"
    echo "=== ${CASE} FP_SEED=${SEED} ==="
    t0=$(date +%s)
    if env FP_SEED="${SEED}" \
      "$FP_BIN" "$ALPHA" "${INPUT_DIR}/${CASE}.block" "${INPUT_DIR}/${CASE}.nets" "$rpt_path" >"$fp_log" 2>&1; then
      fp_exit=0
    else
      fp_exit=$?
    fi
    t1=$(date +%s)
    fp_sec=$((t1 - t0))

    t2=$(date +%s)
    if "$EVAL" "${INPUT_DIR}/${CASE}.block" "${INPUT_DIR}/${CASE}.nets" "$rpt_path" "$ALPHA" >"$eval_log" 2>&1; then
      eval_exit=0
    else
      eval_exit=$?
    fi
    t3=$(date +%s)
    eval_sec=$((t3 - t2))

    norm="NA"
    escore="NA"
    if [[ -f "$eval_log" ]]; then
      norm_line="$(awk -F '=' '/Your normalized cost =/ {line=$0} END{print line}' "$eval_log" 2>/dev/null || true)"
      if [[ -n "${norm_line:-}" ]]; then
        after="${norm_line#*= }"
        after="${after#*= }"
        after="${after#"${after%%[![:space:]]*}"}"
        norm="${after%%[[:space:]]*}"
      fi
      if [[ "$ALPHA" == "0.5" ]]; then
        escore="$(
          awk '/^\[Eval\] SCORE/ {m=1; next} m==1 {print; exit}' "$eval_log" 2>/dev/null \
            | sed -n 's/.*= *\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' \
            | head -n 1
        )"
        if [[ -z "${escore:-}" ]]; then
          escore="$(
            grep -E '^[[:space:]]*=' "$eval_log" 2>/dev/null \
              | tail -n 1 \
              | sed -n 's/.*= *\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'
          )"
        fi
      fi
    fi
    echo "${CASE},${SEED},${ALPHA},${fp_exit},${fp_sec},${rpt_path},${fp_log},${eval_exit},${eval_sec},${norm},${escore},${eval_log}" >>"$CSV_PATH"
  done
done

echo
echo "Wrote: $CSV_PATH"
echo "Reports: $RUN_DIR"
echo
echo "Best per case (lowest normalized cost for alpha=${ALPHA}):"
awk -F, 'NR>1 && $9!="" && $9!="NA" {
  c=$1; n=$9+0
  if (!(c in best) || n < best[c]) { best[c]=n; line[c]=$0 }
}
END { for (c in line) print line[c] }' "$CSV_PATH" | sort -t, -k1,1

