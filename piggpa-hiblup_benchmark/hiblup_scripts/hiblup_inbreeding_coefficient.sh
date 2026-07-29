#!/bin/bash
# =============================================================================
# Task: T5 - Inbreeding Coefficient and Relationship Coefficient
# Function: inbreeding_coefficient
# Source logs:
#   - HIBLUP_results/4/chr1_inbreeding.log                  (--ibc)
#   - HIBLUP_results/4/chr1_relationship_coefficient.log    (--rc)
# Description: Computes inbreeding coefficients (--ibc) and relationship
#              coefficients (--rc) for 1,000 individuals on chr1 SNPs.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile: /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Keep list:   /public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
#   - SNP list:    /public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt
# Outputs (in working directory):
#   - chr1_inbreeding.ibc                  (inbreeding coefficients, Fhat1/Fhat2/Fhat3)
#   - chr1_relationship_coefficient.rc     (pairwise relationship coefficients)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
KEEP=/public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
SNPS=/public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt

# --- Inbreeding coefficient (--ibc) ---
$HIBLUP \
  --ibc \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --out chr1_inbreeding

# --- Relationship coefficient (--rc) ---
$HIBLUP \
  --rc \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --out chr1_relationship_coefficient
