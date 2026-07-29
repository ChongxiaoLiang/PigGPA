#!/usr/bin/env bash
# WF-005: genomic-selection — 6 GS models (BayesA/B/C/BRR/BL/GBLUP) on reference_2k
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

WF_ID="wf005"
WF_SKILL="genomic-selection"
WF_TASK="6 GS models (BayesA/B/C/BRR/BL/GBLUP) on reference_2k with 3-fold CV"
WF_ENTRY="$PIGGPA_SKILLS_DIR/genomic-selection/scripts/GS6models.R"
WF_OUTDIR="$RESULTS_DIR/$WF_ID-$WF_SKILL"
WF_LOG="$WF_OUTDIR/$WF_ID-$WF_SKILL.log"
WF_JSON="$WF_OUTDIR/$WF_ID-$WF_SKILL.json"

mkdir -p "$WF_OUTDIR"

INPUT_FILES=(
    "$TESTDATA_DIR/reference_2k.bed"
    "$TESTDATA_DIR/reference_2k.bim"
    "$TESTDATA_DIR/reference_2k.fam"
)
export WF_INPUT_FILES="$(printf "%s\n" "${INPUT_FILES[@]}")"

echo ">>> [$WF_ID] $WF_SKILL starting at $(date +%H:%M:%S)"

start_t="$(date +%s.%N)"
cd "$WF_OUTDIR"
# === Skill invocation ===
# Step 1: generate .raw file via makerawfile.sh (plink --recodeA)
export PATH="$APP_BIN:$PATH"
bash "$PIGGPA_SKILLS_DIR/genomic-selection/scripts/makerawfile.sh" \
    --bfile "$TESTDATA_DIR/reference_2k" \
    --out "$WF_OUTDIR/gs_raw" \
    --app-bin "$APP_BIN" > "$WF_LOG" 2>&1
MAKERAW_RC=$?
if [[ $MAKERAW_RC -ne 0 ]]; then
    echo "ERROR: makerawfile.sh failed with rc=$MAKERAW_RC" >> "$WF_LOG"
    WF_RC=$MAKERAW_RC
else
    # Step 2: run GS6models.R
    "$RSCRIPT_BIN" "$WF_ENTRY" \
        --input "$WF_OUTDIR/gs_raw.raw" \
        --outdir "$WF_OUTDIR" \
        --output_prefix "gs_6models" \
        --folds 3 \
        --seed 123 \
        --bglr_iter 500 \
        --bglr_burnin 100 >> "$WF_LOG" 2>&1
    WF_RC=$?
fi
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
