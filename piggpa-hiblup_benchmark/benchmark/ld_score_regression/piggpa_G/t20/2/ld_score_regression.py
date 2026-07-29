#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LD Score Regression Script - Complete Version
Estimate heritability using LD Score regression
Using pre-computed LD scores and summary statistics
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings('ignore')

LD_SCORE_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/t19/2/ld_score_results.csv"
SUMSTAT_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/t20/sumstat.csv"
OUTPUT_DIR = "/public/share/likui/liangcx/bole/bole_benchmark/piggpa-hiblup_benchmark_upload/benchmark/ld_score_regression/piggpa_G/t20/2"

os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 12

NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377']


def load_data():
    """Load LD scores and summary statistics"""
    print("  Loading LD scores...")
    ld_df = pd.read_csv(LD_SCORE_FILE)
    print(f"    LD scores: {len(ld_df)} SNPs")
    
    print("  Loading summary statistics...")
    sumstat_df = pd.read_csv(SUMSTAT_FILE)
    print(f"    Summary stats: {len(sumstat_df)} SNPs")
    
    merged = pd.merge(ld_df, sumstat_df, left_on='SNP_ID', right_on='SNP', how='inner')
    print(f"    Merged: {len(merged)} SNPs")
    
    return merged


def calculate_chi2_stats(df):
    """Calculate chi-squared statistics from summary statistics"""
    chi2_stats = np.zeros(len(df))
    
    for i, row in df.iterrows():
        beta = row['BETA']
        se = row['SE']
        
        if pd.isna(se) or se == 0:
            chi2_stats[i] = 0
        else:
            z = beta / se
            chi2_stats[i] = z ** 2
    
    return chi2_stats


def ld_score_regression(chi2_stats, ld_scores, n_samples, n_snps):
    """LD Score regression to estimate heritability"""
    valid_mask = (~np.isnan(chi2_stats)) & (~np.isnan(ld_scores)) & (chi2_stats > 0) & (ld_scores > 0)
    chi2_valid = chi2_stats[valid_mask]
    ld_valid = ld_scores[valid_mask]
    
    if len(chi2_valid) < 10:
        return None, None, None, None, None
    
    x = ld_valid - np.mean(ld_valid)
    y = chi2_valid - np.mean(chi2_valid)
    
    x_mean = np.mean(ld_valid)
    y_mean = np.mean(chi2_valid)
    
    slope = np.sum(x * y) / np.sum(x ** 2)
    intercept = y_mean - slope * x_mean
    
    residuals = chi2_valid - (intercept + slope * ld_valid)
    mse = np.sum(residuals ** 2) / (len(chi2_valid) - 2)
    
    var_slope = mse / np.sum(x ** 2)
    se_slope = np.sqrt(var_slope)
    
    var_intercept = mse * (1/len(chi2_valid) + x_mean**2 / np.sum(x**2))
    se_intercept = np.sqrt(var_intercept)
    
    h2 = slope * n_snps / n_samples
    h2_se = se_slope * n_snps / n_samples
    
    if h2_se > 0:
        z_score = h2 / h2_se
        h2_pval = 2 * (1 - 0.5 * (1 + np.sign(z_score) * (1 - np.exp(-0.717 * z_score**2 - 0.416 * z_score))))
    else:
        h2_pval = 1.0
    
    return h2, h2_se, intercept, se_intercept, h2_pval


def plot_forest(h2, h2_se, intercept, se_intercept, h2_pval, output_path):
    """Plot forest plot for heritability estimate"""
    print("  Plotting forest plot...")
    
    fig, ax = plt.subplots(figsize=(10, 4))
    
    estimates = ['Heritability (h²)', 'Intercept']
    values = [h2, intercept]
    ses = [h2_se, se_intercept]
    
    y_pos = np.arange(len(estimates))
    
    ci_lower = [h2 - 1.96 * h2_se, intercept - 1.96 * se_intercept]
    ci_upper = [h2 + 1.96 * h2_se, intercept + 1.96 * se_intercept]
    
    colors = [NATURE_COLORS[0], NATURE_COLORS[1]]
    
    for i, (est, val, se, lower, upper, color) in enumerate(zip(estimates, values, ses, ci_lower, ci_upper, colors)):
        ax.errorbar(val, i, xerr=[[val - lower], [upper - val]], 
                   fmt='o', color=color, markersize=10, capsize=5, capthick=2, linewidth=2)
    
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=1, color='gray', linestyle=':', alpha=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(estimates, fontsize=14)
    ax.set_xlabel('Estimate', fontsize=14)
    ax.set_title('LD Score Regression Estimates', fontsize=16)
    ax.grid(True, alpha=0.3, axis='x')
    
    textstr = f'h² = {h2:.4f} ± {h2_se:.4f}\nIntercept = {intercept:.4f} ± {se_intercept:.4f}\nP-value = {h2_pval:.2e}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Forest plot saved: {output_path}")


def plot_h2_vs_maf(df, h2, output_path):
    """Plot heritability vs MAF"""
    print("  Plotting h² vs MAF...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1 = axes[0]
    
    if 'MAF_x' in df.columns:
        maf = df['MAF_x'].values
    elif 'MAF_y' in df.columns:
        maf = df['MAF_y'].values
    else:
        maf = df['MAF'].values
    
    chi2 = df['Chi2'].values
    
    valid_mask = (~np.isnan(maf)) & (~np.isnan(chi2)) & (maf > 0) & (chi2 > 0)
    maf_valid = maf[valid_mask]
    chi2_valid = chi2[valid_mask]
    
    ax1.scatter(maf_valid, chi2_valid, alpha=0.3, s=5, c=NATURE_COLORS[2])
    ax1.set_xlabel('Minor Allele Frequency (MAF)', fontsize=14)
    ax1.set_ylabel('χ² Statistic', fontsize=14)
    ax1.set_title('χ² vs MAF', fontsize=16)
    ax1.grid(True, alpha=0.3)
    
    maf_bins = np.linspace(0, 0.5, 11)
    maf_bin_labels = [f'{maf_bins[i]:.1f}-{maf_bins[i+1]:.1f}' for i in range(len(maf_bins)-1)]
    maf_binned = np.digitize(maf_valid, maf_bins) - 1
    maf_binned = np.clip(maf_binned, 0, len(maf_bins)-2)
    
    mean_chi2_by_maf = [np.mean(chi2_valid[maf_binned == i]) if np.sum(maf_binned == i) > 0 else 0 
                        for i in range(len(maf_bins)-1)]
    
    ax2 = axes[1]
    x_pos = np.arange(len(maf_bin_labels))
    bars = ax2.bar(x_pos, mean_chi2_by_maf, color=NATURE_COLORS[3], edgecolor='black', alpha=0.7)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(maf_bin_labels, rotation=45, ha='right')
    ax2.set_xlabel('MAF Bin', fontsize=14)
    ax2.set_ylabel('Mean χ²', fontsize=14)
    ax2.set_title('Mean χ² by MAF Bin', fontsize=16)
    ax2.grid(True, alpha=0.3, axis='y')
    
    ax2.axhline(y=np.mean(chi2_valid), color='red', linestyle='--', 
               label=f'Overall mean = {np.mean(chi2_valid):.2f}')
    ax2.legend(fontsize=10)
    
    plt.suptitle(f'Heritability Analysis (h² = {h2:.4f})', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    h² vs MAF plot saved: {output_path}")


def plot_manhattan(df, output_path):
    """Plot Manhattan plot"""
    print("  Plotting Manhattan plot...")
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    p_values = df['P'].values
    valid_mask = ~np.isnan(p_values) & (p_values > 0) & (p_values <= 1)
    
    positions = df['BP'].values[valid_mask] / 1e6
    p_valid = p_values[valid_mask]
    
    neg_log_p = -np.log10(p_valid)
    
    scatter = ax.scatter(positions, neg_log_p, c=neg_log_p, cmap='RdYlBu_r', 
                         alpha=0.7, s=10)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(p-value)', fontsize=12)
    
    if len(p_valid) > 0:
        bonferroni_threshold = -np.log10(0.05 / len(p_valid))
        ax.axhline(y=bonferroni_threshold, color='red', linestyle='--', linewidth=1,
                   label=f'Bonferroni (p={0.05/len(p_valid):.2e})')
    
    ax.set_xlabel('Position (Mb)', fontsize=14)
    ax.set_ylabel('-log10(p-value)', fontsize=14)
    ax.set_title('Manhattan Plot from LD Score Regression', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Manhattan plot saved: {output_path}")


def plot_ld_score_regression(chi2_stats, ld_scores, h2, intercept, slope, output_path):
    """Plot LD Score regression"""
    print("  Plotting LD Score regression...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    valid_mask = (~np.isnan(chi2_stats)) & (~np.isnan(ld_scores)) & (chi2_stats > 0) & (ld_scores > 0)
    chi2_valid = chi2_stats[valid_mask]
    ld_valid = ld_scores[valid_mask]
    
    ax.scatter(ld_valid, chi2_valid, alpha=0.3, s=5, c=NATURE_COLORS[0])
    
    x_line = np.linspace(ld_valid.min(), ld_valid.max(), 100)
    y_line = intercept + slope * x_line
    ax.plot(x_line, y_line, 'r-', lw=2, label=f'Regression line')
    
    ax.axhline(y=1, color='green', linestyle='--', alpha=0.5, label='Null (χ²=1)')
    
    ax.set_xlabel('LD Score', fontsize=14)
    ax.set_ylabel('χ² Statistic', fontsize=14)
    ax.set_title(f'LD Score Regression (h² = {h2:.4f})', fontsize=16)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    textstr = f'Slope = {slope:.4f}\nIntercept = {intercept:.4f}\nh² = {h2:.4f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    LD Score regression plot saved: {output_path}")


def main():
    print("=" * 70)
    print("LD Score Regression Analysis")
    print("=" * 70)
    
    print("\n[1/5] Loading data...")
    df = load_data()
    
    print("\n[2/5] Calculating chi-squared statistics...")
    chi2_stats = calculate_chi2_stats(df)
    df['Chi2'] = chi2_stats
    
    n_samples = df['N'].iloc[0] if 'N' in df.columns else 1000
    n_snps = len(df)
    
    print(f"  Sample size: {n_samples}")
    print(f"  Number of SNPs: {n_snps}")
    
    print("\n[3/5] Performing LD Score regression...")
    ld_scores = df['LD_Score'].values
    
    h2, h2_se, intercept, se_intercept, h2_pval = ld_score_regression(
        chi2_stats, ld_scores, n_samples, n_snps)
    
    if h2 is None:
        print("  Error: Not enough valid data for regression")
        return
    
    slope = h2 * n_samples / n_snps
    
    print(f"\n  Results:")
    print(f"    Heritability (h²): {h2:.4f} ± {h2_se:.4f}")
    print(f"    Intercept: {intercept:.4f} ± {se_intercept:.4f}")
    print(f"    P-value: {h2_pval:.2e}")
    
    print("\n[4/5] Generating visualizations...")
    plot_forest(h2, h2_se, intercept, se_intercept, h2_pval, 
                os.path.join(OUTPUT_DIR, "forest_plot.pdf"))
    plot_h2_vs_maf(df, h2, os.path.join(OUTPUT_DIR, "h2_vs_maf.pdf"))
    plot_manhattan(df, os.path.join(OUTPUT_DIR, "manhattan_plot.pdf"))
    plot_ld_score_regression(chi2_stats, ld_scores, h2, intercept, slope,
                             os.path.join(OUTPUT_DIR, "ld_score_regression.pdf"))
    
    print("\n[5/5] Saving results...")
    
    results_df = pd.DataFrame({
        'Item': ['Heritability', 'Intercept'],
        'Estimate': [h2, intercept],
        'SE': [h2_se, se_intercept],
        'P_value': [h2_pval, np.nan]
    })
    
    output_results = os.path.join(OUTPUT_DIR, "ld_score_regression_results.csv")
    results_df.to_csv(output_results, index=False)
    print(f"  Results saved: {output_results}")
    
    output_summary = os.path.join(OUTPUT_DIR, "ld_score_regression_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("LD Score Regression Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sample size: {n_samples}\n")
        f.write(f"Number of SNPs: {n_snps}\n\n")
        
        f.write("Results:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Item: Heritability (h²)\n")
        f.write(f"    Estimate: {h2:.4f}\n")
        f.write(f"    SE: {h2_se:.4f}\n")
        f.write(f"    P-value: {h2_pval:.2e}\n\n")
        f.write(f"  Item: Intercept\n")
        f.write(f"    Estimate: {intercept:.4f}\n")
        f.write(f"    SE: {se_intercept:.4f}\n\n")
        
        f.write("Interpretation:\n")
        f.write("-" * 40 + "\n")
        if h2 > 0.5:
            f.write("  High heritability - trait is strongly influenced by genetics\n")
        elif h2 > 0.2:
            f.write("  Moderate heritability - trait has both genetic and environmental influences\n")
        else:
            f.write("  Low heritability - trait is primarily influenced by environment\n")
        
        if intercept > 1.5:
            f.write("  High intercept suggests potential population stratification or cryptic relatedness\n")
        elif intercept > 1.1:
            f.write("  Slightly elevated intercept - some confounding may be present\n")
        else:
            f.write("  Intercept close to 1 - minimal confounding\n")
    print(f"  Summary saved: {output_summary}")
    
    detailed_results = pd.DataFrame({
        'Item': ['h2', 'Intercept'],
        'Estimate': [h2, intercept],
        'SE': [h2_se, se_intercept],
        'P_value': [h2_pval, np.nan]
    })
    
    detailed_file = os.path.join(OUTPUT_DIR, "heritability_estimates.csv")
    detailed_results.to_csv(detailed_file, index=False)
    print(f"  Heritability estimates saved: {detailed_file}")
    
    print("\n" + "=" * 70)
    print("LD Score Regression Analysis Complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. ld_score_regression_results.csv - Regression results")
    print(f"  2. ld_score_regression_summary.txt - Analysis summary")
    print(f"  3. heritability_estimates.csv - Heritability estimates")
    print(f"  4. forest_plot.pdf - Forest plot")
    print(f"  5. h2_vs_maf.pdf - Heritability vs MAF")
    print(f"  6. manhattan_plot.pdf - Manhattan plot")
    print(f"  7. ld_score_regression.pdf - LD Score regression plot")


if __name__ == "__main__":
    main()
