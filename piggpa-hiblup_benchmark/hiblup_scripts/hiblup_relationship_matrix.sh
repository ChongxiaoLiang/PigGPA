#!/bin/bash
# =============================================================================
# Task: T3 - Relationship Matrix Construction
# Function: relationship_matrix
# Source logs:
#   - HIBLUP_results/3/chr1_GA.log  (Genomic Additive relationship matrix)
#   - HIBLUP_results/3/chr1_GD.log  (Genomic Dominance relationship matrix)
#   - HIBLUP_results/3/chr1_HA.log  (Hybrid/pedigree+genomic Additive matrix)
#   - HIBLUP_results/3/chr1_PA.log  (Pedigree Additive relationship matrix)
# Description: Constructs four relationship matrices using HIBLUP --make-xrm:
#              GA (genomic additive), GD (genomic dominance),
#              HA (hybrid additive, pedigree + genotype, alpha=0.05),
#              PA (pedigree additive, genotype-free).
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile:   /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - Keep list:     /public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
#   - SNP list:      /public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt
#   - Pedigree file: /public/share/likui/hanyu/testdata/In-silico-data/pedigree_1000.txt
# Outputs (per command, in working directory):
#   - chr1_GA.GA.bin / chr1_GA.GA.id  (GA matrix)
#   - chr1_GD.GD.bin / chr1_GD.GD.id  (GD matrix)
#   - chr1_HA.HA.bin / chr1_HA.HA.id  (HA matrix)
#   - chr1_PA.PA.bin / chr1_PA.PA.id  (PA matrix)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
KEEP=/public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt
SNPS=/public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt
PED=/public/share/likui/hanyu/testdata/In-silico-data/pedigree_1000.txt

# --- GA: Genomic Additive relationship matrix ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --add \
  --out chr1_GA

# --- GD: Genomic Dominance relationship matrix ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --dom \
  --out chr1_GD

# --- HA: Hybrid Additive (pedigree + genomic), alpha=0.05 ---
$HIBLUP \
  --make-xrm \
  --bfile $BFILE \
  --keep $KEEP \
  --extract $SNPS \
  --pedigree $PED \
  --add \
  --alpha 0.05 \
  --out chr1_HA

# --- PA: Pedigree Additive relationship matrix (no genotype) ---
$HIBLUP \
  --make-xrm \
  --pedigree $PED \
  --add \
  --out chr1_PA
