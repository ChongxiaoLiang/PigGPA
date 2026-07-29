#!/bin/bash
# =============================================================================
# Task: T20 - LD Score Regression
# Function: ld_score_regression
# Source log: HIBLUP_results/15/ldreg_result.log
# Description: Runs LD score regression (--ldreg) to estimate SNP-heritability
#              and other regression-based statistics from GWAS summary
#              statistics. The LD scores come from the T18 ld_calculation
#              step (ldscore_result.ldsc), and the summary statistics come
#              from sumstat_hiblup.txt (prepared separately).
#              NOTE: This command FAILED at runtime with
#              "Error: not a float value '1:29726' at 2th row and 8th column
#              of the file [sumstat_hiblup.txt]". The 8th column of the
#              summary-statistics file contained a "chr:pos" string instead
#              of a numeric value, which HIBLUP's parser rejected. The
#              command is reproduced here verbatim from the log for
#              reproducibility / traceability; fixing the input file format
#              is outside the scope of this wrapper.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - Summary stats:  sumstat_hiblup.txt   (local)
#   - LD scores:      /public/share/likui/hanyu/testresult/HIBLUP/14/ldscore_result.ldsc
# Outputs (in working directory):
#   - ldreg_result.*   (regression results; this run FAILED, see note above)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
SUMSTAT=sumstat_hiblup.txt
LDS=/public/share/likui/hanyu/testresult/HIBLUP/14/ldscore_result.ldsc

# --- LD score regression (NOTE: failed at runtime, see header) ---
$HIBLUP \
  --ldreg \
  --sumstat $SUMSTAT \
  --lds $LDS \
  --thread 4 \
  --out ldreg_result
