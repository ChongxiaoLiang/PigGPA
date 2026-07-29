#!/usr/bin/env bash
# WF-003: plink-gwas — PLINK linear regression GWAS on reference (2238 samples)
set -uo pipefail
set +H
ulimit -v unlimited

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIGGPA_SKILLS_DIR="${PIGGPA_SKILLS_DIR:-/public/share/likui/liangcx/bole/skills/piggpa-G}"
TESTDATA_DIR="${TESTDATA_DIR:-/public/share/likui/liangcx/bole/testdata}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
APP_BIN="${APP_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/sys_tools/bin}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/R/bin/Rscript}"
PYTHON_BIN="${PYTHON_BIN:-/public/share/likui/liangcx/software/miniconda3/envs/py_analysis/bin/python}"

WF_ID="wf003"
WF_SKILL="plink-gwas"
WF_TASK="PLINK linear regression GWAS on reference + GWAStestphenotype"
WF_ENTRY="$PIGGPA_SKILLS_DIR/plink-gwas/scripts/PLINK-GWAS-linear.sh"
WF_OUTDIR="$RESULTS_DIR/$WF_ID-$WF_SKILL"
WF_LOG="$WF_OUTDIR/$WF_ID-$WF_SKILL.log"
WF_JSON="$WF_OUTDIR/$WF_ID-$WF_SKILL.json"

mkdir -p "$WF_OUTDIR"

INPUT_FILES=(
    "$TESTDATA_DIR/reference.bed"
    "$TESTDATA_DIR/reference.bim"
    "$TESTDATA_DIR/reference.fam"
    "$TESTDATA_DIR/GWAStestphenotype.fam"
)
export WF_INPUT_FILES="$(printf "%s\n" "${INPUT_FILES[@]}")"

echo ">>> [$WF_ID] $WF_SKILL starting at $(date +%H:%M:%S)"

start_t="$(date +%s.%N)"
cd "$WF_OUTDIR"
# === Skill invocation ===
# Note: uses `bash -c 'exec bash "$0" "$@"'` pattern. The previous `setsid bash ...`
# form was broken: when the calling shell is a process-group leader, setsid
# forks and returns immediately (rc=0) without waiting for the child, which
# produced false "passed" results with elapsed ~0.006s. The bash -c exec form
# runs the skill in the foreground of the current shell so $? captures the
# real exit code. If the multi-threaded C/C++ binary still receives SIGKILL
# (rc=137) from the session cgroup, that is accepted as a system-level limit.
export PATH="$APP_BIN:$PATH"
bash -c 'exec bash "$0" "$@"' "$WF_ENTRY" \
    --genotype_prefix "$TESTDATA_DIR/reference" \
    --phenotype_file "$TESTDATA_DIR/GWAStestphenotype.fam" \
    --output_prefix "plink_gwas" \
    --app-bin "$APP_BIN" \
    --rscript-bin "$RSCRIPT_BIN" \
    --threads 4 > "$WF_LOG" 2>&1 < /dev/null
WF_RC=$?
end_t="$(date +%s.%N)"
WF_ELAPSED="$(awk "BEGIN{printf \"%.3f\", $end_t - $start_t}")"

if awk "BEGIN{exit !($WF_ELAPSED < 60)}"; then
    WF_ELAPSED_HUMAN="$(printf "%.2f s" "$WF_ELAPSED")"
elif awk "BEGIN{exit !($WF_ELAPSED < 3600)}"; then
    m="$(awk "BEGIN{printf \"%d\", $WF_ELAPSED / 60}")"
    s="$(awk "BEGIN{printf \"%.0f\", $WF_ELAPSED - $m * 60}")"
    WF_ELAPSED_HUMAN="$m min $s s"
else
    h="$(awk "BEGIN{printf \"%d\", $WF_ELAPSED / 3600}")"
    rem="$(awk "BEGIN{printf \"%.0f\", $WF_ELAPSED - $h * 3600}")"
    m="$(awk "BEGIN{printf \"%d\", $rem / 60}")"
    s="$(awk "BEGIN{printf \"%.0f\", $rem - $m * 60}")"
    WF_ELAPSED_HUMAN="$h h $m min $s s"
fi

if [[ $WF_RC -eq 0 ]]; then
    WF_STATUS="passed"
else
    WF_STATUS="failed"
fi

WF_OUTPUT_FILES=""
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    base="$(basename "$f")"
    [[ "$base" == "${WF_ID}-${WF_SKILL}.log" ]] && continue
    [[ "$base" == "${WF_ID}-${WF_SKILL}.json" ]] && continue
    WF_OUTPUT_FILES+="${f}"$'\n'
done < <(find "$WF_OUTDIR" -maxdepth 2 -type f 2>/dev/null | sort)

export WF_ID WF_SKILL WF_TASK WF_ENTRY WF_RC WF_ELAPSED WF_ELAPSED_HUMAN WF_STATUS WF_LOG WF_JSON WF_OUTPUT_FILES
"$PYTHON_BIN" "$SCRIPT_DIR/lib/write_task_json.py"

echo "<<< [$WF_ID] $WF_SKILL rc=$WF_RC elapsed=${WF_ELAPSED_HUMAN} status=$WF_STATUS"
exit $WF_RC
