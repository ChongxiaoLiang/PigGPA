#!/bin/bash
# =============================================================================
# Task: T7 - Five-Model BLUP Breeding Value Prediction
# Function: blup_prediction
# Source logs (HIBLUP_results/6/1/):
#   - blup_train.log,      blup_pred.log,      blup_snp_effect.log
#   - gblup_train.log,     gblup_pred.log,     gblup_snp_effect.log
#   - lm_train.log,        lm_pred.log,        lm_snp_effect.log
#   - pblup_train.log
#   - ssblup_train.log,    ssblup_pred.log,    ssblup_snp_effect.log
# Description: Reproduces the full 5-model breeding-value prediction pipeline
#              for T7. For each model (BLUP/GBLUP/LM/PBLUP/SSBLUP):
#                (a) train single-trait model to estimate variance components
#                    and GEBVs (--single-trait);
#                (b) back-calculate SNP effects from GEBVs (--snp-effect);
#                (c) predict BVs on held-out individuals (--pred).
#              PBLUP has no genotype data and therefore no SNP-effect / pred
#              steps (it only produces training GEBVs from the pedigree).
#              Differences between models:
#                * BLUP   : GA only, vc-method = EMAI
#                * GBLUP  : GA only, default AI-REML
#                * LM     : GA only, vc-method = HE (Henderson's REML-like)
#                * PBLUP  : pedigree only, no genotype
#                * SSBLUP : HA (alpha=0.05 hybrid), pedigree + genotype
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:    /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Phenotype file: /public/share/likui/hanyu/testdata/In-silico-data/phenotypes.txt
#   - Pedigree file:  pedigree_hiblup.txt  (local, prepared in HIBLUP_results/6/1/)
#   - Train samples:  train_samples.txt    (local)
#   - Pred samples:   pred_samples.txt     (local)
#   - SNP list:       chr1_snps.txt        (local)
# Outputs (per model, in working directory):
#   - {model}_train.vars / .beta / .rand     (variance components, fixed, random effects)
#   - {model}_snp_effect.snpeff              (SNP effects)
#   - {model}_pred.bv                        (predicted breeding values)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
PHENO=/public/share/likui/hanyu/testdata/In-silico-data/phenotypes.txt
PED=pedigree_hiblup.txt
TRAIN=train_samples.txt
PRED=pred_samples.txt
SNPS=chr1_snps.txt

# =============================================================================
# Model 1: BLUP (GBLUP with EMAI variance-component method)
# =============================================================================
# Train
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $TRAIN \
  --extract $SNPS \
  --pheno-pos 2 \
  --add \
  --vc-method EMAI \
  --thread 4 \
  --float-prec \
  --out blup_train

# SNP effect from GEBVs
$HIBLUP \
  --snp-effect \
  --bfile $BFILE \
  --keep $TRAIN \
  --extract $SNPS \
  --gebv blup_train.rand \
  --thread 4 \
  --out blup_snp_effect

# Predict
$HIBLUP \
  --pred \
  --bfile $BFILE \
  --keep $PRED \
  --extract $SNPS \
  --score blup_snp_effect.snpeff \
  --thread 4 \
  --out blup_pred

# =============================================================================
# Model 2: GBLUP (default AI-REML)
# =============================================================================
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $TRAIN \
  --extract $SNPS \
  --pheno-pos 2 \
  --add \
  --thread 4 \
  --float-prec \
  --out gblup_train

$HIBLUP \
  --snp-effect \
  --bfile $BFILE \
  --keep $TRAIN \
  --extract $SNPS \
  --gebv gblup_train.rand \
  --thread 4 \
  --out gblup_snp_effect

$HIBLUP \
  --pred \
  --bfile $BFILE \
  --keep $PRED \
  --extract $SNPS \
  --score gblup_snp_effect.snpeff \
  --thread 4 \
  --out gblup_pred

# =============================================================================
# Model 3: LM (GA with HE variance-component method)
# =============================================================================
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --keep $TRAIN \
  --extract $SNPS \
  --pheno-pos 2 \
  --add \
  --vc-method HE \
  --he-pred \
  --thread 4 \
  --float-prec \
  --out lm_train

$HIBLUP \
  --snp-effect \
  --bfile $BFILE \
  --keep $TRAIN \
  --extract $SNPS \
  --gebv lm_train.rand \
  --thread 4 \
  --out lm_snp_effect

$HIBLUP \
  --pred \
  --bfile $BFILE \
  --keep $PRED \
  --extract $SNPS \
  --score lm_snp_effect.snpeff \
  --thread 4 \
  --out lm_pred

# =============================================================================
# Model 4: PBLUP (pedigree-only, no genotype -> no SNP effect / pred steps)
# =============================================================================
$HIBLUP \
  --single-trait \
  --pheno $PHENO \
  --pedigree $PED \
  --keep $TRAIN \
  --pheno-pos 2 \
  --thread 4 \
  --float-prec \
  --out pblup_train

# =============================================================================
# Model 5: SSBLUP (hybrid HA, alpha=0.05: 0.95*G + 0.05*A)
# =============================================================================
$HIBLUP \
  --single-trait \
  --bfile $BFILE \
  --pheno $PHENO \
  --pedigree $PED \
  --keep $TRAIN \
  --extract $SNPS \
  --pheno-pos 2 \
  --add \
  --alpha 0.05 \
  --thread 4 \
  --float-prec \
  --out ssblup_train

$HIBLUP \
  --snp-effect \
  --bfile $BFILE \
  --keep $TRAIN \
  --extract $SNPS \
  --gebv ssblup_train.rand \
  --thread 4 \
  --out ssblup_snp_effect

$HIBLUP \
  --pred \
  --bfile $BFILE \
  --keep $PRED \
  --extract $SNPS \
  --score ssblup_snp_effect.snpeff \
  --thread 4 \
  --out ssblup_pred
