#!/bin/bash
# Author: Anicet Ebou
# Tested only on Linux Ubuntu 24.04
#
# Directory layout created at the current working directory:
#   results/   all result CSVs
#   figures/   all figures, prefixed 01_, 02_, ... in order of appearance
#   rds/       checkpoint .rds files (one per step; used both for the
#              skip-if-already-run mechanism below and, for 07, as the
#              only inputs figures are generated from)
#   data/      raw/intermediate data for step 05 (genome selection table,
#              downloaded FASTAs, FastANI results)
#   logs/      one timestamped log per step, plus a run-level summary
#
# Usage:
#   ./00_orchestrator.sh                 # run everything, skip completed
#                                         # steps (checkpoint-based)
#   ./00_orchestrator.sh --force         # ignore all checkpoints, redo
#                                         # every step from scratch
#   ./00_orchestrator.sh --from 04       # start from step 04 onward
#                                         # (earlier steps' checkpoints
#                                         # must already exist)
#   ./00_orchestrator.sh --only 05       # run only step 05
#
# Environment variables consulted by individual steps:
#   RUN_FULL_GRID   (03_realistic_sketch_calib.R): set to FALSE to stop
#                   after the runtime probe rather than committing to the
#                   full grid on the first run
#   FORCE_RERUN_05  (05_real_genomes.sh): set to TRUE to force
#                   every shell sub-part (A/B/C) inside step 06 to redo,
#                   even if their individual outputs already exist
#                   (automatically set when --force is passed here)

# deliberately NOT -e: a single failed step should not
# silently abort logging/summary generation
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

RESULTS_DIR="results"
FIGURES_DIR="figures"
RDS_DIR="rds"
DATA_DIR="data"
LOGS_DIR="logs"

mkdir -p "$RESULTS_DIR" "$FIGURES_DIR" "$RDS_DIR" "$DATA_DIR" "$LOGS_DIR"

RUN_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
PIPELINE_LOG="${LOGS_DIR}/pipeline_${RUN_TIMESTAMP}.log"
TIMING_CSV="${LOGS_DIR}/timing_summary.csv"

# argument parsing
FORCE_ALL="FALSE"
FROM_STEP=""
ONLY_STEP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE_ALL="TRUE"; shift ;;
    --from)  FROM_STEP="$2"; shift 2 ;;
    --only)  ONLY_STEP="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ "$FORCE_ALL" == "TRUE" ]]; then
  export FORCE_RERUN_06="TRUE"
fi

log() {
  local msg="[orchestrator] $(date '+%Y-%m-%d %H:%M:%S') $1"
  echo "$msg" | tee -a "$PIPELINE_LOG"
}

# Initialize timing CSV with header if it doesn't exist yet
if [[ ! -f "$TIMING_CSV" ]]; then
  echo "run_timestamp,step,status,elapsed_seconds,log_file" > "$TIMING_CSV"
fi

STEP_ORDER=(01 02 03 04 05 06)
declare -A STEP_SCRIPT=(
  [01]="01_wald_approx.R"
  [02]="02_idealized_binom.R"
  [03]="03_realistic_sketch_calib.R"
  [04]="04_staircase_structure.R"
  [05]="05_real_genomes.sh"
  [06]="06_figures.R"
)
declare -A STEP_CHECKPOINT=(
  [01]="${RDS_DIR}/01_wald_approx.rds"
  [02]="${RDS_DIR}/02_idealized_binom.rds"
  [03]="${RDS_DIR}/03_realistic_sketch_calib.rds"
  [04]="${RDS_DIR}/04_staircase_structure.rds"
  [05]="${RDS_DIR}/05_real_genomes.rds"
  [06]="${FIGURES_DIR}/06_real_genome_bias_vs_simulated.pdf"
)
declare -A STEP_KIND=(
  [01]="R" [02]="R" [03]="R" [04]="R" [05]="SH" [06]="R"
)

# step selection (--from / --only)

steps_to_run=()
if [[ -n "$ONLY_STEP" ]]; then
  steps_to_run=("$ONLY_STEP")
elif [[ -n "$FROM_STEP" ]]; then
  started="FALSE"
  for s in "${STEP_ORDER[@]}"; do
    if [[ "$s" == "$FROM_STEP" ]]; then started="TRUE"; fi
    if [[ "$started" == "TRUE" ]]; then steps_to_run+=("$s"); fi
  done
else
  steps_to_run=("${STEP_ORDER[@]}")
fi

log "=== Pipeline run started ==="
log "Steps to run: ${steps_to_run[*]}"
log "Force rerun: $FORCE_ALL"
log "Pipeline log: $PIPELINE_LOG"

overall_start=$(date +%s)
any_failed="FALSE"

# run one step, with skip-if-checkpoint-exists and timing

run_step() {
  local step="$1"
  local script="${STEP_SCRIPT[$step]}"
  local checkpoint="${STEP_CHECKPOINT[$step]}"
  local kind="${STEP_KIND[$step]}"
  local step_log="${LOGS_DIR}/step${step}_${RUN_TIMESTAMP}.log"

  if [[ ! -f "$script" ]]; then
    log "ERROR: step $step script '$script' not found in $SCRIPT_DIR, aborting this step..."
    echo "${RUN_TIMESTAMP},${step},MISSING_SCRIPT,0,${step_log}" >> "$TIMING_CSV"
    any_failed="TRUE"
    return
  fi

  if [[ -f "$checkpoint" && "$FORCE_ALL" != "TRUE" ]]; then
    log "Step $step ($script): checkpoint '$checkpoint' exists, SKIPPING..."
    echo "${RUN_TIMESTAMP},${step},SKIPPED,0,${step_log}" >> "$TIMING_CSV"
    return
  fi

  log "Step $step ($script): starting..."
  local t0 t1 elapsed status
  t0=$(date +%s)

  if [[ "$kind" == "R" ]]; then
    Rscript "$script" > "$step_log" 2>&1
    status=$?
  else
    bash "$script" > "$step_log" 2>&1
    status=$?
  fi

  t1=$(date +%s)
  elapsed=$((t1 - t0))

  # Also append this step's full output into the master pipeline log
  {
    echo "----- BEGIN step $step output ($script) -----"
    cat "$step_log"
    echo "----- END step $step output -----"
  } >> "$PIPELINE_LOG"

  if [[ $status -ne 0 ]]; then
    log "Step $step ($script): FAILED after ${elapsed}s (exit code $status). See $step_log"
    echo "${RUN_TIMESTAMP},${step},FAILED,${elapsed},${step_log}" >> "$TIMING_CSV"
    any_failed="TRUE"
  elif [[ ! -f "$checkpoint" ]]; then
    log "Step $step ($script): completed in ${elapsed}s but expected checkpoint '$checkpoint' was NOT created -- check $step_log"
    echo "${RUN_TIMESTAMP},${step},NO_CHECKPOINT,${elapsed},${step_log}" >> "$TIMING_CSV"
    any_failed="TRUE"
  else
    log "Step $step ($script): completed successfully in ${elapsed}s."
    echo "${RUN_TIMESTAMP},${step},SUCCESS,${elapsed},${step_log}" >> "$TIMING_CSV"
  fi
}

# main loop

for step in "${steps_to_run[@]}"; do
  if [[ -z "${STEP_SCRIPT[$step]+x}" ]]; then
    log "ERROR: unknown step '$step' (valid: ${STEP_ORDER[*]}) -- skipping."
    any_failed="TRUE"
    continue
  fi
  run_step "$step"
done

overall_end=$(date +%s)
overall_elapsed=$((overall_end - overall_start))

log "=== Pipeline run finished in ${overall_elapsed}s ==="
log "Timing summary: $TIMING_CSV"
log "Full log: $PIPELINE_LOG"

if [[ "$any_failed" == "TRUE" ]]; then
  log "One or more steps FAILED or did not produce their expected checkpoint -- see above."
  exit 1
fi

log "All requested steps completed successfully."
exit 0
