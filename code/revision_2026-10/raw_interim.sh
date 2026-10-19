#!/bin/bash
# Interim or final analysis of the raw-signal derivation: participant-level results -> summary -> BMI models.
# Usage: raw_interim.sh <output-folder> [prefix|all] [R_MAX]
#   prefix (default): random sample = first R participants of analytic_random_order.tsv that are all processed; all: every participant with four series.
# Paths (adapt to the machine): the working-tree root, the folder with the raw-signal outputs and the slim PAXMIN caches.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
export SRI_PROJECT_ROOT="${SRI_PROJECT_ROOT:-$(cd "$here/../../.." && pwd)}"
export RAW_ROOT="${RAW_ROOT:-/Volumes/LaCie/NHANES_1112_1314}"
export PAXMIN_1112="${PAXMIN_1112:-$SRI_PROJECT_ROOT/01_recompute_SRI/intermediate/slim_paxmin/PAXMIN_slim_2011_12.rds}"
export PAXMIN_1314="${PAXMIN_1314:-$SRI_PROJECT_ROOT/01_recompute_SRI/intermediate/slim_paxmin/PAXMIN_slim_2013_14.rds}"
export EXCLUDED="${EXCLUDED:-$here/excluded_seqn.txt}"   # participants whose SWaN run returned no result (excluded from the prefix)
out="${1:?output folder}"; mkdir -p "$out"; export OUT_DIR="$out" IN_DIR="$out" MODE="${2:-prefix}"; [ -n "$3" ] && export R_MAX="$3"
[ -f "$RAW_ROOT/_swan/failed.log" ] && echo "note: $RAW_ROOT/_swan/failed.log lists SWaN runs without result"
Rscript "$here/raw_derivation_participants.R"          > "$out/log_participants.txt" 2>&1
Rscript "$here/raw_derivation_summary.R"               > "$out/log_summary.txt" 2>&1
SRI_COL=sri_A Rscript "$here/raw_derivation_models.R"  > "$out/log_models_A.txt" 2>&1
SRI_COL=sri_B Rscript "$here/raw_derivation_models.R"  > "$out/log_models_B.txt" 2>&1
echo "done: $out"
