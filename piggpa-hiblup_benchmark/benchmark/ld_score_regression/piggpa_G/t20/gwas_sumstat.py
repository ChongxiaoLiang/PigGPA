#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS Summary Statistics Generation Script
Generate summary statistics including SNP effects, p-values, sample size
Using samples 1001-2000, chromosome 1
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings('ignore')

try:
    from pandas_plink import read_plink
except ImportError:
    print("Installing pandas-plink...")
    os.system("pip install pandas-plink")
    from pandas_plink import read_plink

INPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data"
OUTPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t20"
BED_FILE = os.path.join(INPUT_DIR, "simulated_population")
PHENOTYPE_FILE = os.path.join(INPUT_DIR, "phenotypes.txt")

SAMPLE_START = 1000
SAMPLE_END = 2000
TARGET_CHROM = "1"

os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 12

NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377']


def calculate_maf(genotype):
    maf = np.mean(genotype, axis=0) / 2
    maf = np.minimum(maf, 1 - maf)
    return maf


def t_distribution_cdf(t, df):
    x = df / (df + t**2)
    return 0.5 * (1 + np.sign(t) * np.sqrt(1 - x * (df / (df + 4)) * (1 + x * (df / (df + 6)) * (1 + x * (df / (df + 8))))))


def perform_gwas(genotype, phenotype):
    n_snps = genotype.shape[1]
    n_samples = len(phenotype)
    
    effects = np.zeros(n_snps)
    standard_errors = np.zeros(n_snps)
    p_values = np.zeros(n_snps)
    t_stats = np.zeros(n_snps)
    
    y_centered = phenotype - np.mean(phenotype)
    
    print(f"  Performing GWAS: {n_snps} SNPs, {n_samples} samples")
    
    for i in range(n_snps):
        x = genotype[:, i]
        
        if np.var(x) == 0:
            effects[i] = 0
            standard_errors[i] = np.nan
            p_values[i] = 1
            t_stats[i] = 0
            continue
        
        x_centered = x - np.mean(x)
        
        slope = np.sum(x_centered * y_centered) / np.sum(x_centered ** 2)
        
        y_pred = np.mean(phenotype) + slope * x_centered
        residuals = phenotype - y_pred
        
        mse = np.sum(residuals ** 2) / (n_samples - 2)
        
        se = np.sqrt(mse / np.sum(x_centered ** 2))
        
        t_stat = slope / se if se > 0 else 0
        
        abs_t = abs(t_stat)
        p_value = 2 * (1 - t_distribution_cdf(abs_t, n_samples - 2))
        
        effects[i] = slope
        standard_errors[i] = se
        p_values[i] = p_value
        t_stats[i] = t_stat
        
        if (i + 1) % 1000 == 0:
            print(f"    Processing {i+1}/{n_snps} SNPs...")
    
    return effects, standard_errors, p_values, t_stats


def plot_manhattan(snp_info, p_values, output_path):
    print("  Plotting Manhattan plot...")
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    neg_log_p = -np.log10(p_values + 1e-300)
    
    positions = snp_info['pos'].values / 1e6
    
    scatter = ax.scatter(positions, neg_log_p, c=neg_log_p, cmap='RdYlBu_r', 
                         alpha=0.7, s=10)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(p-value)', fontsize=12)
    
    bonferroni_threshold = -np.log10(0.05 / len(p_values))
    ax.axhline(y=bonferroni_threshold, color='red', linestyle='--', linewidth=1,
               label=f'Bonferroni (p={0.05/len(p_values):.2e})')
    
    ax.set_xlabel('Position (Mb)', fontsize=14)
    ax.set_ylabel('-log10(p-value)', fontsize=14)
    ax.set_title(f'GWAS Manhattan Plot - Chromosome {TARGET_CHROM}', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Manhattan plot saved: {output_path}")


def plot_qq(p_values, output_path):
    print("  Plotting QQ plot...")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    valid_mask = ~np.isnan(p_values) & (p_values > 0) & (p_values <= 1)
    p_values_valid = p_values[valid_mask]
    
    if len(p_values_valid) == 0:
        ax.text(0.5, 0.5, 'No valid p-values', ha='center', va='center', fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    QQ plot saved: {output_path}")
        return
    
    p_values_sorted = np.sort(p_values_valid)
    n = len(p_values_sorted)
    
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(p_values_sorted)
    
    ax.scatter(expected, observed, alpha=0.5, s=5, c=NATURE_COLORS[0])
    
    max_val = max(np.max(expected), np.max(observed)) * 1.1
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1, label='Expected')
    
    ax.set_xlabel('Expected -log10(p)', fontsize=14)
    ax.set_ylabel('Observed -log10(p)', fontsize=14)
    ax.set_title('QQ Plot', fontsize=16)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    lambda_gc = np.median(p_values_sorted) / 0.5
    ax.text(0.05, 0.95, f'lambda = {lambda_gc:.3f}', transform=ax.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    QQ plot saved: {output_path}")


def plot_effect_distribution(effects, p_values, output_path):
    print("  Plotting effect distribution...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1 = axes[0]
    ax1.hist(effects, bins=50, color=NATURE_COLORS[1], edgecolor='black', alpha=0.7)
    ax1.axvline(x=np.mean(effects), color='red', linestyle='--', 
               label=f'Mean={np.mean(effects):.4f}')
    ax1.set_xlabel('Effect Size', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=14)
    ax1.set_title('Distribution of Effect Sizes', fontsize=16)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[1]
    ax2.scatter(np.abs(effects), -np.log10(p_values + 1e-300), 
               alpha=0.5, s=10, c=NATURE_COLORS[2])
    ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', 
               label='p=0.05')
    ax2.set_xlabel('|Effect Size|', fontsize=14)
    ax2.set_ylabel('-log10(p-value)', fontsize=14)
    ax2.set_title('Effect Size vs Significance', fontsize=16)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Effect distribution saved: {output_path}")


def main():
    print("=" * 70)
    print("GWAS Summary Statistics Generation")
    print(f"Sample range: {SAMPLE_START+1}-{SAMPLE_END}")
    print(f"Chromosome: {TARGET_CHROM}")
    print("=" * 70)
    
    print("\n[1/5] Reading PLINK files...")
    bim, fam, bed = read_plink(BED_FILE, verbose=False)
    print(f"  Total samples: {len(fam)}")
    print(f"  Total SNPs: {len(bim)}")
    
    print(f"\n[2/5] Filtering chromosome {TARGET_CHROM} SNPs...")
    chr_mask = bim['chrom'].astype(str) == TARGET_CHROM
    chr_snps = bim[chr_mask].copy()
    snp_indices = chr_snps.index.tolist()
    snp_info = chr_snps.reset_index(drop=True)
    print(f"  Chromosome {TARGET_CHROM} SNPs: {len(snp_indices)}")
    
    print(f"\n[3/5] Extracting genotype data (samples {SAMPLE_START+1}-{SAMPLE_END})...")
    genotype = bed[snp_indices, SAMPLE_START:SAMPLE_END].compute().T
    genotype = np.nan_to_num(genotype, nan=0)
    print(f"  Samples: {genotype.shape[0]}, SNPs: {genotype.shape[1]}")
    
    print(f"\n[4/5] Loading phenotype data...")
    pheno_df = pd.read_csv(PHENOTYPE_FILE, sep='\t')
    pheno_df['ID'] = pheno_df['ID'].astype(int)
    pheno_subset = pheno_df[(pheno_df['ID'] > SAMPLE_START) & (pheno_df['ID'] <= SAMPLE_END)]
    pheno_subset = pheno_subset.set_index('ID').sort_index()
    phenotype = pheno_subset['Phenotype'].values
    print(f"  Phenotypes: {len(phenotype)}")
    print(f"  Phenotype mean: {np.mean(phenotype):.4f}")
    print(f"  Phenotype std: {np.std(phenotype):.4f}")
    
    print(f"\n[5/5] Performing GWAS analysis...")
    effects, standard_errors, p_values, t_stats = perform_gwas(genotype, phenotype)
    
    maf = calculate_maf(genotype)
    
    n_samples = len(phenotype)
    
    print(f"\n  GWAS results:")
    print(f"    Significant SNPs (p<0.05): {np.sum(p_values < 0.05)}")
    print(f"    Significant SNPs (p<0.01): {np.sum(p_values < 0.01)}")
    print(f"    Significant SNPs (Bonferroni): {np.sum(p_values < 0.05/len(p_values))}")
    print(f"    Effect range: [{np.min(effects):.4f}, {np.max(effects):.4f}]")
    print(f"    Effect mean: {np.mean(effects):.6f}")
    
    print("\n[6/5] Generating visualizations...")
    plot_manhattan(snp_info, p_values, os.path.join(OUTPUT_DIR, "manhattan_plot.pdf"))
    plot_qq(p_values, os.path.join(OUTPUT_DIR, "qq_plot.pdf"))
    plot_effect_distribution(effects, p_values, os.path.join(OUTPUT_DIR, "effect_distribution.pdf"))
    
    print("\n[7/5] Saving results...")
    
    sumstat_df = pd.DataFrame({
        'SNP': snp_info['snp'].values,
        'CHR': snp_info['chrom'].values,
        'BP': snp_info['pos'].values,
        'A1': snp_info['a1'].values,
        'A2': snp_info['a0'].values,
        'MAF': maf,
        'BETA': effects,
        'SE': standard_errors,
        'T': t_stats,
        'P': p_values,
        'N': n_samples
    })
    
    output_file = os.path.join(OUTPUT_DIR, "sumstat.csv")
    sumstat_df.to_csv(output_file, index=False)
    print(f"  Summary statistics saved: {output_file}")
    
    output_summary = os.path.join(OUTPUT_DIR, "gwas_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("GWAS Summary Statistics\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sample range: {SAMPLE_START+1}-{SAMPLE_END}\n")
        f.write(f"Chromosome: {TARGET_CHROM}\n")
        f.write(f"SNPs: {len(snp_indices)}\n")
        f.write(f"Sample size: {n_samples}\n\n")
        
        f.write("GWAS Results:\n")
        f.write(f"  Significant SNPs (p<0.05): {np.sum(p_values < 0.05)}\n")
        f.write(f"  Significant SNPs (p<0.01): {np.sum(p_values < 0.01)}\n")
        f.write(f"  Significant SNPs (Bonferroni): {np.sum(p_values < 0.05/len(p_values))}\n\n")
        
        f.write("Effect Statistics:\n")
        f.write(f"  Effect range: [{np.min(effects):.4f}, {np.max(effects):.4f}]\n")
        f.write(f"  Effect mean: {np.mean(effects):.6f}\n")
        f.write(f"  Effect std: {np.std(effects):.4f}\n\n")
        
        f.write("Phenotype Statistics:\n")
        f.write(f"  Phenotype mean: {np.mean(phenotype):.4f}\n")
        f.write(f"  Phenotype std: {np.std(phenotype):.4f}\n")
    print(f"  Summary saved: {output_summary}")
    
    print("\n" + "=" * 70)
    print("GWAS Summary Statistics Generation Complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. sumstat.csv - GWAS summary statistics (SNP, effect, p-value, N)")
    print(f"  2. gwas_summary.txt - Analysis summary")
    print(f"  3. manhattan_plot.pdf - Manhattan plot")
    print(f"  4. qq_plot.pdf - QQ plot")
    print(f"  5. effect_distribution.pdf - Effect distribution")


if __name__ == "__main__":
    main()
