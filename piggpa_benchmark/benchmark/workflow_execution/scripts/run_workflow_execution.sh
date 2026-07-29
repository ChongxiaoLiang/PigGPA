#!/usr/bin/env bash
# PigGPA workflow_execution benchmark orchestrator.
#
# Thin driver that runs every per-task wrapper run_wfNNN_<skill>.sh in
# numeric order, then aggregates the per-task JSONs into
# workflow_execution_results.json via lib/aggregate_results.py.
#
# Usage:
#   bash run_workflow_execution.sh                # run all 26 tasks
#   bash run_workflow_execution.sh wf001          # run only tasks matching "wf001"
#   bash run_workflow_execution.sh 'wf01[0-7]'    # regex filter
#
# Environment overrides (with defaults):
#   PIGGPA_SKILLS_DIR  (default: /public/share/likui/liangcx/bole/skills/piggpa-G)
#   TESTDATA_DIR       (default: /public/share/likui/liangcx/bole/testdata)
#   RESULTS_DIR        (default: <this script's dir>/results)
#   APP_BIN            (default: /public/share/likui/liangcx/software/miniconda3/envs/sys_tools/bin)
#   RSCRIPT_BIN        (default: /public/share/likui/liangcx/software/miniconda3/envs/R/bin/Rscript)
#   PYTHON_BIN         (default: /public/share/likui/liangcx/software/miniconda3/envs/py_analysis/bin/python)

set -uo pipefail
ulimit -v unlimited

# ------------------------------------------------------------------
# Resolve directories and binaries
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PIGGPA_SKILLS_DIR="${PIGGPA_SKILLS_DIR:-/public/share/likui/liangcx/bole/skills/piggpa-G}"
TESTDATA_DIR="${TESTDATA_DIR:-/public/share/likui/liangcx/bole/testdata}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
APP_BIN="${APP_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/sys_tools/bin}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/R/bin/Rscript}"
PYTHON_BIN="${PYTHON_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/py_analysis/bin/python}"

# Export for child per-task scripts to consume
export PIGGPA_SKILLS_DIR TESTDATA_DIR RESULTS_DIR APP_BIN RSCRIPT_BIN PYTHON_BIN

mkdir -p "$RESULTS_DIR"

# Optional regex filter (default: match all wfNNN_* scripts)
FILTER="${1:-wf.*}"

echo "=== PigGPA workflow_execution benchmark orchestrator ==="
echo "Skills dir : $PIGGPA_SKILLS_DIR"
echo "Testdata   : $TESTDATA_DIR"
echo "Results dir: $RESULTS_DIR"
echo "APP_BIN    : $APP_BIN"
echo "Rscript    : $RSCRIPT_BIN"
echo "Python     : $PYTHON_BIN"
echo "Filter     : $FILTER"
echo

# ------------------------------------------------------------------
# Iterate per-task scripts in numeric order
# ------------------------------------------------------------------
# Disable -e so one failing task doesn't abort the whole benchmark.
# Each run_wfNNN_*.sh is responsible for capturing its own exit code
# and writing its own per-task JSON regardless of success/failure.
set +e
declare -a RUN_SCRIPTS=()
while IFS= read -r script; do
    RUN_SCRIPTS+=("$script")
done < <(ls "$SCRIPT_DIR"/run_wf[0-9][0-9][0-9]_*.sh 2>/dev/null | sort)

if [[ ${#RUN_SCRIPTS[@]} -eq 0 ]]; then
    echo "ERROR: no run_wfNNN_*.sh scripts found in $SCRIPT_DIR" >&2
    exit 1
fi

started_at="$(date +%s.%N)"
n_run=0
n_passed=0
n_failed=0
for script in "${RUN_SCRIPTS[@]}"; do
    base="$(basename "$script")"
    if [[ ! "$base" =~ $FILTER ]]; then
        echo "--- skip (filter): $base"
        continue
    fi
    n_run=$((n_run + 1))
    echo ">>> [$(date +%H:%M:%S)] running $base"
    bash "$script"
    rc=$?
    if [[ $rc -eq 0 ]]; then
        n_passed=$((n_passed + 1))
        echo "<<< $base OK (rc=0)"
    else
        n_failed=$((n_failed + 1))
        echo "<<< $base FAILED (rc=$rc) — continuing"
    fi
done
ended_at="$(date +%s.%N)"
total_elapsed="$(awk "BEGIN{printf \"%.2f\", $ended_at - $started_at}")"

echo
echo "=== Orchestrator summary ==="
echo "  ran    : $n_run"
echo "  passed : $n_passed"
echo "  failed : $n_failed"
echo "  wall   : ${total_elapsed} s"
echo

# ------------------------------------------------------------------
# Aggregate per-task JSONs into a single results file
# ------------------------------------------------------------------
echo "=== Aggregating results ==="
"$PYTHON_BIN" "$SCRIPT_DIR/lib/aggregate_results.py" \
    --results-dir "$RESULTS_DIR" \
    --output "$SCRIPT_DIR/workflow_execution_results.json" \
    --testdata-dir "$TESTDATA_DIR" \
    --skills-dir "$PIGGPA_SKILLS_DIR"

echo
echo "=== Done. Aggregate: $SCRIPT_DIR/workflow_execution_results.json ==="
