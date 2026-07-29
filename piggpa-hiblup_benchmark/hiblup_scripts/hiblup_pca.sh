#!/bin/bash
# =============================================================================
# Task: T6 - Principal Component Analysis (PCA)
# Function: pca
# Source log: HIBLUP_results/5/chr1_pca.log
# Description: Performs PCA on the genomic relationship matrix of 1,000
#              individuals (chr1 SNPs), computing the top 10 PCs.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile: /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Keep list:   /public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
#   - SNP list:    /public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt
# Outputs (in working directory):
#   - chr1_pca.pc   (PC scores per individual)
#   - chr1_pca.pcp  (PC proportion / variance explained)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
KEEP=/public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
SNPS=/public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt

# --- PCA with 10 principal components ---
$HIBLUP \
  --pca \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --npc 10 \
  --out chr1_pca
