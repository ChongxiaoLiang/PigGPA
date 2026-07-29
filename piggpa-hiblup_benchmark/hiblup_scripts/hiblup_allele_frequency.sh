#!/bin/bash
# =============================================================================
# Task: T1/T2 - Allele/Genotype Frequency and Heterozygosity/Homozygosity
# Function: allele_frequency
# Source logs:
#   - HIBLUP_results/1.2/chr1_allele_freq.log
#   - HIBLUP_results/1.2/chr1_geno_freq.log
#   - HIBLUP_results/1.2/chr1_heterozygosity.log
#   - HIBLUP_results/1.2/chr1_homozygosity.log
# Description: Calculates allele frequency, genotype frequency, heterozygosity
#              and homozygosity for chromosome 1 SNPs using HIBLUP.
# HIBLUP version: v1.6.0 (2025-09-29 Release)
# Binary: /public/share/likui/hanyu/software/bin/hiblup
# Inputs:
#   - PLINK bfile: /public/share/likui/hanyu/testdata/In-silico-data/simulated_population
#   - SNP list:    /public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt
# Outputs (per command, in working directory):
#   - chr1_allele_freq.afreq        (allele frequencies)
#   - chr1_geno_freq.gfreq          (genotype frequencies)
#   - chr1_heterozygosity.hete      (heterozygosity per individual)
#   - chr1_homozygosity.homo        (homozygosity per individual)
# =============================================================================
set -e

HIBLUP=/public/share/likui/hanyu/software/bin/hiblup
BFILE=/public/share/likui/hanyu/testdata/In-silico-data/simulated_population
SNPS=/public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt

# --- T1: Allele frequency ---
$HIBLUP \
  --allele-freq \
  --bfile $BFILE \
  --extract $SNPS \
  --out chr1_allele_freq

# --- T1: Genotype frequency ---
$HIBLUP \
  --geno-freq \
  --bfile $BFILE \
  --extract $SNPS \
  --out chr1_geno_freq

# --- T2: Heterozygosity ---
$HIBLUP \
  --hete \
  --bfile $BFILE \
  --extract $SNPS \
  --out chr1_heterozygosity

# --- T2: Homozygosity ---
$HIBLUP \
  --homo \
  --bfile $BFILE \
  --extract $SNPS \
  --out chr1_homozygosity
