#!/bin/bash
# =============================================================================
# Task: T8 - Single-Trait Variance Component Estimation (GBLUP)
# Function: single_trait_model
# Source log: HIBLUP_results/7/single_trait_gblup.log
# Description: Fits a single-trait GBLUP model (GA random effect) on chr1 SNPs
#              to estimate additive-genetic variance, residual variance and
#              heritability (h^2). This is the HIBLUP-side counterpart to
#              piggpa-G's single_trait_model.py for direct T8 comparison.
#              Note: phenotype file is phenotype_hiblup.txt (re-formatted with
#              the HIBLUP-required header), and sample list is samples.txt.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:   /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: phenotype_hiblup.txt  (local)
#   - Sample list:   samples.txt            (local)
#   - SNP list:      chr1_snps.txt          (local)
# Outputs (in working directory):
#   - single_trait_gblup.vars   (Va, Ve, h^2 and SEs)
#   - single_trait_gblup.beta   (fixed-effect estimates)
#   - single_trait_gblup.rand   (random-effect GEBVs)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=phenotype_hiblup.txt
KEEP=samples.txt
SNPS=chr1_snps.txt

# --- Single-trait GBLUP (GA random effect, phenotype at column 2) ---
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --pheno-pos 2 \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --thread 4 \
  --out single_trait_gblup
