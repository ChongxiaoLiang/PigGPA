#!/bin/bash
# =============================================================================
# Task: T9 - Repeated Records Model (Permanent Environment + Additive Genetic)
# Function: repeated_records_model
# Source log: HIBLUP_results/8/repeated_model.log
# Description: Fits a single-trait model with repeated records. The model
#              includes a permanent-environment random effect (--rand 1, the
#              individual ID column) on top of the additive-genetic GA effect,
#              plus two discrete covariates (--dcovar 2,3). Phenotype column
#              is at position 5 in the long-format phenotype file
#              (phenotype_long.txt). This is the HIBLUP-side counterpart to
#              piggpa-G's repeated_records_model.py for direct T9 comparison.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:    /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: phenotype_long.txt   (local, long format with repeated records)
#   - Train samples:  train_samples.txt    (local)
#   - SNP list:       chr1_snps.txt        (local)
# Outputs (in working directory):
#   - repeated_model.vars    (Va, Vpe, Ve, h^2, pe^2 and SEs)
#   - repeated_model.beta    (fixed-effect estimates)
#   - repeated_model.rand    (random-effect estimates: PE + GA)
#   - repeated_model.anova   (ANOVA table)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=phenotype_long.txt
KEEP=train_samples.txt
SNPS=chr1_snps.txt

# --- Repeated-records model: GA + PE (rand column 1), covariates 2,3 ---
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $KEEP \
  --extract $SNPS \
  --pheno-pos 5 \
  --dcovar 2,3 \
  --rand 1 \
  --add \
  --thread 4 \
  --float-prec \
  --out repeated_model
