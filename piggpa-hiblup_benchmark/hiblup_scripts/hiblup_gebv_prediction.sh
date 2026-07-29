#!/bin/bash
# =============================================================================
# Task: T16 - GEBV Prediction on New Individuals
# Function: gebv_prediction
# Source log: HIBLUP_results/12/prediction_result.log
# Description: Predicts genomic estimated breeding values (GEBVs) for
#              held-out individuals (--pred) using the SNP effects computed
#              in the snp_effect_calculation step (HIBLUP_results/11/).
#              The --score file (snp_effect.snpeff) is loaded and applied to
#              the genotype matrix of the prediction-set individuals to
#              produce per-individual BVs.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:   /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Pred samples:  pred_samples.txt    (local)
#   - SNP list:      chr1_snps.txt       (local)
#   - SNP effects:   /public/share/likui/hanyu/testresult/HIBLUP/11/snp_effect.snpeff
# Outputs (in working directory):
#   - prediction_result.bv   (predicted BVs for the prediction-set individuals)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
KEEP=pred_samples.txt
SNPS=chr1_snps.txt
SCORE=/public/share/likui/hanyu/testresult/HIBLUP/11/snp_effect.snpeff

# --- Predict GEBVs for the held-out individuals ---
$HIBLUP \
  --pred \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --score $SCORE \
  --thread 4 \
  --out prediction_result
