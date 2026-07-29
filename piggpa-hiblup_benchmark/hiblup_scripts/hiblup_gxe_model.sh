#!/bin/bash
# =============================================================================
# Task: T12 - Gene-Environment (GxE) Interaction Model
# Function: gxe_model
# Source logs (HIBLUP_results/10/):
#   - GA.log          (build GA matrix used by the GxE random effect)
#   - gxe_xrm2.log    (build GxE-augmented XRM with --rand-gxe + --write-txt)
#   - gxe_model.log   (main GxE single-trait fit, converged)
#   - grm_check.log   (GRM verification, uses relative paths from project root)
#   - gxe_he.log      (FAILED: --algorithm HE not supported for --rand-gxe)
#   - gxe_hi.log      (FAILED: --algorithm HI not supported for --rand-gxe)
# Description: Fits a single-trait GxE model. Environmental covariates are
#              columns 2 and 3 of the phenotype file (qcovar), and the
#              interaction random effect is built on the pre-computed GA
#              matrix (--rand-gxe 2:GA.GA,3:GA.GA). The phenotype is at
#              column 4. Two auxiliary commands build the GA matrix and the
#              GxE-augmented XRM (with --write-txt for verification).
#              NOTE: The gxe_he.log and gxe_hi.log attempts to use
#              --algorithm HE / --algorithm HI together with --rand-gxe
#              FAILED with "The following arguments were not expected".
#              They are reproduced here for completeness but are commented
#              out by default (HIBLUP does not support HE/HI for GxE models).
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:    /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: /public/share/likui/hanyu/testdata/In-silico-data/t12/2/simulated_phenotypes_env.txt
#   - Train samples:  train_samples.txt   (local)
#   - SNP list:       chr1_snps.txt       (local)
# Outputs (in working directory):
#   - GA.GA.bin / GA.GA.id              (GA matrix, prerequisite for --rand-gxe)
#   - gxe_xrm2.GA.bin / .GA.id / .txt   (GxE-augmented XRM, verification)
#   - gxe_model.vars / .beta / .rand    (main GxE model results)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=/public/share/likui/hanyu/testdata/In-silico-data/t12/2/simulated_phenotypes_env.txt
KEEP=train_samples.txt
SNPS=chr1_snps.txt

# --- (1) Build GA matrix (prerequisite for --rand-gxe) ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --thread 4 \
  --float-prec \
  --out GA

# --- (2) Build GxE-augmented XRM with TXT output (for verification) ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --pheno $PHENO \
  --rand-gxe 2:GA.GA,3:GA.GA \
  --add \
  --write-txt \
  --thread 4 \
  --out gxe_xrm2

# --- (3) Main GxE single-trait model (converged) ---
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $KEEP \
  --extract $SNPS \
  --pheno-pos 4 \
  --qcovar 2,3 \
  --add \
  --rand-gxe 2:GA.GA,3:GA.GA \
  --thread 4 \
  --float-prec \
  --out gxe_model

# =============================================================================
# Below commands are reproduced from log files but FAILED at runtime and are
# therefore commented out. HIBLUP v1.6.0 does NOT support --algorithm HE or
# --algorithm HI together with --rand-gxe (GxE). Kept here for traceability.
# =============================================================================
# --- (FAILED) gxe_he: --algorithm HE not supported with --rand-gxe ---
# $HIBLUP --single-trait --bfile $BFILE --pheno $PHENO --keep $KEEP \
#   --extract $SNPS --pheno-pos 4 --qcovar 2,3 --add \
#   --rand-gxe 2:GA.GA,3:GA.GA --thread 4 --algorithm HE --out gxe_he
#
# --- (FAILED) gxe_hi: --algorithm HI not supported with --rand-gxe ---
# $HIBLUP --single-trait --bfile $BFILE --pheno $PHENO --keep $KEEP \
#   --extract $SNPS --pheno-pos 4 --qcovar 2,3 --add \
#   --rand-gxe 2:GA.GA,3:GA.GA --thread 4 --algorithm HI --out gxe_hi
#
# --- (verification, relative-path variant) grm_check: GRM dump to TXT ---
# NOTE: original log used relative paths from the project root; paths adjusted
# here to match the absolute layout used by the other commands in this script.
# $HIBLUP --make-xrm --bfile $BFILE --keep $KEEP --extract $SNPS \
#   --add --write-txt --thread 4 --out grm_check
