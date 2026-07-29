#!/bin/bash
# =============================================================================
# Task: T10 - Multi-Trait Model (3 traits)
# Function: multi_trait_model
# Source logs (HIBLUP_results/9/):
#   - multi_trait.log    (main multi-trait analysis)
#   - geno_code1.log     (genotype coding for verification, --trans-geno)
#   - grm_verify.log     (GA matrix dump to TXT for verification, --write-txt)
# Description: Fits a 3-trait GBLUP model (--multi-trait) with GA random
#              effect. Traits are at pheno columns 2, 3, 4. Each trait shares
#              discrete covariates at columns 5 and 6 (--dcovar 5,6 5,6 5,6).
#              Two auxiliary verification commands are included:
#                (a) --trans-geno with --code-method 1: exports the genotype
#                    matrix coded as 0/1/2 (VanRaden method I) for cross-check
#                    against piggpa-G;
#                (b) --make-xrm --write-txt: dumps the GA matrix as a plain
#                    text square matrix for cross-check against piggpa-G.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:    /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: /public/share/likui/hanyu/testdata/In-silico-data/t10/2/simulated_phenotypes_multi_trait.txt
#   - Train samples:  train_samples.txt   (local)
#   - SNP list:       chr1_snps.txt       (local)
# Outputs (in working directory):
#   - multi_trait.vars            (3x3 Va, Ve covariance matrices, h^2 per trait)
#   - multi_trait.T1/T2/T3.*      (per-trait beta/rand/anova)
#   - multi_trait.covars          (covariate estimates)
#   - geno_code1.geno.A.txt       (genotype matrix 0/1/2, verification)
#   - grm_verify.GA.bin/.id/.txt  (GA matrix dump, verification)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=/public/share/likui/hanyu/testdata/In-silico-data/t10/2/simulated_phenotypes_multi_trait.txt
KEEP=train_samples.txt
SNPS=chr1_snps.txt

# --- (1) Multi-trait GBLUP: 3 traits (cols 2,3,4), covariates 5,6 per trait ---
$HIBLUP \
  --multi-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $KEEP \
  --extract $SNPS \
  --pheno-pos 2 3 4 \
  --dcovar 5,6 5,6 5,6 \
  --add \
  --thread 4 \
  --float-prec \
  --out multi_trait

# --- (2) Verification: genotype coding (VanRaden method I = 0/1/2) ---
$HIBLUP \
  --trans-geno \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --code-method 1 \
  --out geno_code1

# --- (3) Verification: dump GA matrix to TXT (square format) ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --write-txt \
  --thread 4 \
  --out grm_verify
