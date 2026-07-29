#!/bin/bash
# =============================================================================
# Task: T18 - Linkage Disequilibrium (LD) Calculation
# Function: ld_calculation
# Source logs (HIBLUP_results/14/):
#   - ld_result.log      (pairwise LD r calculation, --ld)
#   - ldscore_result.log (per-SNP LD scores, --ldscore)
#   - convert.log        (FAILED: --trans-xrm could not find ld_result.id)
# Description: Computes pairwise LD (correlation r) between SNPs in windows
#              of 1 Mb (--ld), and per-SNP LD scores (--ldscore) which are
#              later used by the LD score regression step (T20,
#              hiblup_ld_score_regression.sh). Both commands operate on
#              1,000 individuals restricted to chr1 SNPs.
#              NOTE: convert.log attempted to convert the LD .bin to text
#              with --trans-xrm, but FAILED with "can not open the file
#              [ld_result.id]". The --ld output uses .info (not .id) as the
#              companion index, so --trans-xrm is not applicable here. The
#              conversion was instead done with a custom Python script
#              (convert_ld_bin_to_txt.py / convert_ld_bin_to_txt_v2.py) in
#              the original benchmark; that script is NOT part of the
#              HIBLUP binary and is therefore commented out here. The
#              convert.log command itself is reproduced for traceability.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:   /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Train samples: train_samples.txt    (local)
#   - SNP list:      chr1_snps.txt        (local)
# Outputs (in working directory):
#   - ld_result.bin / ld_result.info    (pairwise LD r, binary + SNP index)
#   - ldscore_result.ldsc               (per-SNP LD scores)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
KEEP=train_samples.txt
SNPS=chr1_snps.txt

# --- (1) Pairwise LD (r) in 1 Mb windows ---
$HIBLUP \
  --ld \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --thread 4 \
  --out ld_result

# --- (2) Per-SNP LD scores (input for T20 LD score regression) ---
$HIBLUP \
  --ldscore \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --thread 4 \
  --out ldscore_result

# =============================================================================
# (FAILED) convert.log: attempt to dump LD .bin to text via --trans-xrm.
# Error: "can not open the file [ld_result.id]" -- --ld output uses .info, not
# .id, so --trans-xrm is not applicable. Conversion was done with a separate
# Python script in the original benchmark; not reproduced here.
# =============================================================================
# $HIBLUP --trans-xrm --xrm ld_result --write-txt --out ld_result_txt
