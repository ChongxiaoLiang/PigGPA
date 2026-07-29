#!/usr/bin/env bash
# WF-026: internal-impute — multi-chip internal imputation (chip_20 target + ref_78), chr1
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

WF_ID="wf026"
WF_SKILL="internal-impute"
WF_TASK="V1 multi-chip internal imputation: gdp_chip_20 (target) + gdp_ref_78 (ref), chr1"
WF_ENTRY="$PIGGPA_SKILLS_DIR/internal-impute/scripts/internal_imputation.py"
WF_OUTDIR="$RESULTS_DIR/$WF_ID-$WF_SKILL"
WF_LOG="$WF_OUTDIR/$WF_ID-$WF_SKILL.log"
WF_JSON="$WF_OUTDIR/$WF_ID-$WF_SKILL.json"

mkdir -p "$WF_OUTDIR"

INPUT_FILES=(
    "$TESTDATA_DIR/gdp_chip_20.bed"
    "$TESTDATA_DIR/gdp_chip_20.bim"
    "$TESTDATA_DIR/gdp_chip_20.fam"
    "$TESTDATA_DIR/gdp_ref_78.bed"
    "$TESTDATA_DIR/gdp_ref_78.bim"
    "$TESTDATA_DIR/gdp_ref_78.fam"
)
export WF_INPUT_FILES="$(printf "%s\n" "${INPUT_FILES[@]}")"

echo ">>> [$WF_ID] $WF_SKILL starting at $(date +%H:%M:%S)"

start_t="$(date +%s.%N)"
cd "$WF_OUTDIR"
# === Skill invocation (V1: chip_20 target + ref_78 secondary, single chromosome) ===
# Bug fix applied: internal_imputation.py now reads BEAGLE_MEM_MB and BEAGLE_NTHREADS
# from env (was hardcoded 32000/16, causing SIGKILL on test environment). Reduce to
# 8000MB/4 threads for small chr1 test data (170K SNP × 98 samples after merge).
export PLINK_BIN="$APP_BIN/plink" PLINK2_BIN="$APP_BIN/plink2" BCFTOOLS_BIN="$APP_BIN/bcftools" BGZIP_BIN="$APP_BIN/bgzip"
export BEAGLE_JAR="/public/share/likui/liangcx/software/beagle.jar"
export BEAGLE_MEM_MB=8000
export BEAGLE_NTHREADS=4
export PATH="$APP_BIN:$PATH"
"$PYTHON_BIN" "$WF_ENTRY" \
    -i "$TESTDATA_DIR/gdp_chip_20" \
    -o "$WF_OUTDIR" \
    --input-secondary "$TESTDATA_DIR/gdp_ref_78" \
    --chrom 1 \
    --dosage-format plink > "$WF_LOG" 2>&1
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
