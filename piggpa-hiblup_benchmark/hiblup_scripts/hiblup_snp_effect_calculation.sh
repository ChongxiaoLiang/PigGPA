#!/bin/bash
# =============================================================================
# Task: SNP Effect Calculation (HIBLUP function 3.18)
# Function: snp_effect_calculation
# Source logs (HIBLUP_results/11/1/):
#   - single_trait.log       (GBLUP training to obtain GEBVs)
#   - snp_effect.log         (SNP effect back-calculation from GEBVs)
#   - snp_effect_test.log    (FAILED: --gebv not specified)
# Description: Two-step pipeline that first fits a single-trait GBLUP model
#              on a 1,000-individual / chr1-SNP subset to estimate GEBVs, then
#              back-calculates per-SNP additive effects from those GEBVs
#              (--snp-effect). The resulting snp_effect.snpeff file is the
#              input for downstream GEBV prediction (T16, see
#              hiblup_gebv_prediction.sh).
#              NOTE: This directory (HIBLUP_results/11/) was originally marked
#              "?" in the task mapping. Inspection of its analysis report
#              (分析报告.md) and log files shows it is the HIBLUP 3.18
#              "SNP effect calculation" task. It is the prerequisite step
#              feeding T16 (gebv_prediction), and is
#              therefore kept as a standalone wrapper script.
#              NOTE: snp_effect_test.log attempted to run --snp-effect without
#              --gebv and FAILED with "'--gbvfile' should be specified with
#              prior calculated GEBVs". It is reproduced here commented out
#              for traceability.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:    /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: phenotype.txt        (local)
#   - Keep list:      keep_samples.txt     (local)
#   - Extract SNPs:   extract_snps.txt     (local)
# Outputs (in working directory):
#   - single_trait.vars / .beta / .rand    (variance components, fixed, GEBVs)
#   - snp_effect.snpeff                    (per-SNP additive effects)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=phenotype.txt
KEEP=keep_samples.txt
SNPS=extract_snps.txt

# --- (1) Single-trait GBLUP to obtain GEBVs (column 'GA' in .rand) ---
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --pheno-pos 2 \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --out single_trait \
  --thread 4

# --- (2) Back-calculate SNP effects from GEBVs ---
$HIBLUP \
  --snp-effect \
  --bfile $BFILE \
  --gebv single_trait.rand \
  --keep $KEEP \
  --extract $SNPS \
  --out snp_effect \
  --thread 4

# =============================================================================
# (FAILED) snp_effect_test: attempted --snp-effect without --gebv.
# Error: "'--gbvfile' should be specified with prior calculated GEBVs for
# calculating SNP effect." Reproduced here for traceability only.
# =============================================================================
# $HIBLUP --snp-effect --bfile $BFILE --keep $KEEP --extract $SNPS \
#   --out snp_effect_test --thread 4
