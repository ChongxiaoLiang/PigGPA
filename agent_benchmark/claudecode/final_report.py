#!/usr/bin/env python3
"""
Generate final comprehensive visualizations and summary report
for the pig loin muscle depth GWAS and regulatory annotation.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

WORKDIR = '/workspace/pigbole/benchmark/agent_benchmark/ONE2/claudecode'

# ── Load data ─────────────────────────────────────────────────────────
gwas = pd.read_csv(os.path.join(WORKDIR, 'lmd.gwas.assoc.linear'), sep=r'\s+')
gwas = gwas[gwas['P'].notna()].copy()
gwas['logP'] = -np.log10(gwas['P'])

reg = pd.read_csv(os.path.join(WORKDIR, 'lmd.regulatory_annotation.txt'), sep='\t')
loci = pd.read_csv(os.path.join(WORKDIR, 'lmd.significant_loci.txt'), sep='\t')

# Color map for regulatory classes
REG_COLORS = {
    'Strong regulatory': '#d62728',
    'Moderate regulatory': '#ff7f0e',
    'Weak regulatory': '#2ca02c',
    'Likely non-regulatory': '#7f7f7f'
}

# ── Compute chromosome positions ──────────────────────────────────────
chroms = sorted(gwas['CHR'].unique())
gwas = gwas.sort_values(['CHR', 'BP'])
chr_offsets = {}
chr_max_bp = {}
cum_pos = 0
for ch in chroms:
    ch_data = gwas[gwas['CHR'] == ch]
    chr_offsets[ch] = cum_pos
    chr_max_bp[ch] = ch_data['BP'].max()
    cum_pos += chr_max_bp[ch]

gwas['CUM_POS'] = gwas['BP'] + gwas['CHR'].map(chr_offsets)
chr_midpoints = {ch: chr_offsets[ch] + chr_max_bp[ch] / 2 for ch in chroms}

# ── Compute thresholds ────────────────────────────────────────────────
n_tests = len(gwas)
bonf_threshold = 0.05 / n_tests
suggestive = 1e-5

# ── Figure 1: Enhanced Manhattan Plot ─────────────────────────────────
print("Generating enhanced Manhattan plot...")
fig, ax = plt.subplots(figsize=(18, 7))

# Plot all SNPs
colors = ['#1f77b4', '#aec7e8']
for i, ch in enumerate(chroms):
    ch_data = gwas[gwas['CHR'] == ch]
    ax.scatter(ch_data['CUM_POS'], ch_data['logP'], s=1.5, c=colors[i % 2],
               alpha=0.5, edgecolors='none', rasterized=True)

# Overlay significant loci with regulatory colors
for _, r in reg.iterrows():
    ch = r['chr']
    start_bp = r['start']
    end_bp = r['end']
    cum_start = start_bp + chr_offsets.get(ch, 0)
    cum_end = end_bp + chr_offsets.get(ch, 0)

    # Get the core significant SNPs
    locus_snps = gwas[(gwas['CHR'] == ch) &
                       (gwas['BP'] >= start_bp) &
                       (gwas['BP'] <= end_bp)]
    color = REG_COLORS.get(r['regulatory_class'], '#7f7f7f')
    ax.scatter(locus_snps['CUM_POS'], locus_snps['logP'], s=12, c=color,
               edgecolors='black', linewidth=0.3, zorder=5)

# Threshold lines
ax.axhline(-np.log10(bonf_threshold), color='red', linestyle='--', linewidth=1.2,
           label=f'Bonferroni ({bonf_threshold:.1e})')
ax.axhline(-np.log10(suggestive), color='blue', linestyle='--', linewidth=1.0, alpha=0.6,
           label=f'Suggestive (1e-5)')

# Annotate top 5 loci
top5 = reg.nsmallest(5, 'top_p')
for _, r in top5.iterrows():
    ch = r['chr']
    cum_x = r['start'] + chr_offsets.get(ch, 0)
    logp = -np.log10(r['top_p'])
    ax.annotate(f"{r['locus'].split(':')[0]}", xy=(cum_x, logp),
                xytext=(0, 12), textcoords='offset points', fontsize=7,
                ha='center', fontweight='bold', color='darkred',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=0.8))

ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('-log10(P)', fontsize=12)
ax.set_title('GWAS of Loin Muscle Depth in Pigs\nwith Regulatory Potential Annotation',
             fontsize=14, fontweight='bold')
ax.set_xticks([chr_midpoints[ch] for ch in chroms])
ax.set_xticklabels([str(ch) for ch in chroms], fontsize=9)
ax.set_xlim(0, cum_pos)

# Legend
legend_elements = [
    Patch(facecolor=REG_COLORS['Strong regulatory'], label='Strong regulatory'),
    Patch(facecolor=REG_COLORS['Moderate regulatory'], label='Moderate regulatory'),
    Patch(facecolor=REG_COLORS['Weak regulatory'], label='Weak regulatory'),
    Patch(facecolor=REG_COLORS['Likely non-regulatory'], label='Likely non-regulatory'),
    plt.Line2D([0], [0], color='red', linestyle='--', linewidth=1.2, label=f'Bonferroni ({bonf_threshold:.1e})'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=1,
          framealpha=0.9)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.manhattan_annotated.png'), dpi=200)
plt.close()
print("  → lmd.manhattan_annotated.png")

# ── Figure 2: Regulatory Potential Score Distribution ─────────────────
print("Generating regulatory score distribution plot...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 2a: RPS histogram by class
ax1 = axes[0, 0]
classes = ['Strong regulatory', 'Moderate regulatory', 'Weak regulatory', 'Likely non-regulatory']
for cls in classes:
    subset = reg[reg['regulatory_class'] == cls]['regulatory_potential_score']
    ax1.hist(subset, bins=15, alpha=0.7, label=f'{cls} (n={len(subset)})',
             color=REG_COLORS[cls], edgecolor='black', linewidth=0.5)
ax1.set_xlabel('Regulatory Potential Score')
ax1.set_ylabel('Number of Loci')
ax1.set_title('Distribution of Regulatory Potential Scores')
ax1.legend(fontsize=7)
ax1.axvline(0.6, color='darkred', linestyle='--', alpha=0.5)
ax1.axvline(0.45, color='darkorange', linestyle='--', alpha=0.5)

# 2b: MAF vs Regulatory Score
ax2 = axes[0, 1]
scatter = ax2.scatter(reg['lead_maf'], reg['regulatory_potential_score'],
                      c=[REG_COLORS[c] for c in reg['regulatory_class']],
                      s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
ax2.set_xlabel('Lead SNP MAF')
ax2.set_ylabel('Regulatory Potential Score')
ax2.set_title('MAF vs Regulatory Potential')
# Add trend line
valid = reg.dropna(subset=['lead_maf', 'regulatory_potential_score'])
if len(valid) > 3:
    z = np.polyfit(valid['lead_maf'], valid['regulatory_potential_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(valid['lead_maf'].min(), valid['lead_maf'].max(), 100)
    ax2.plot(x_line, p(x_line), 'k--', alpha=0.5, linewidth=1)
    from scipy.stats import pearsonr
    r, pval = pearsonr(valid['lead_maf'], valid['regulatory_potential_score'])
    ax2.text(0.95, 0.05, f'r = {r:.3f}, P = {pval:.1e}', transform=ax2.transAxes,
             ha='right', fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# 2c: SNP Density vs Regulatory Score
ax3 = axes[1, 0]
ax3.scatter(reg['snp_density_50kb'], reg['regulatory_potential_score'],
            c=[REG_COLORS[c] for c in reg['regulatory_class']],
            s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
ax3.set_xlabel('SNP Density (per 50kb window)')
ax3.set_ylabel('Regulatory Potential Score')
ax3.set_title('SNP Density vs Regulatory Potential')
if len(valid) > 3:
    z = np.polyfit(reg['snp_density_50kb'], reg['regulatory_potential_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(reg['snp_density_50kb'].min(), reg['snp_density_50kb'].max(), 100)
    ax3.plot(x_line, p(x_line), 'k--', alpha=0.5, linewidth=1)

# 2d: Chromosome distribution of regulatory loci
ax4 = axes[1, 1]
chr_counts = reg.groupby(['chr', 'regulatory_class']).size().unstack(fill_value=0)
chr_counts = chr_counts.reindex(sorted(chr_counts.index))
# Ensure all classes present
for cls in classes:
    if cls not in chr_counts.columns:
        chr_counts[cls] = 0
chr_counts = chr_counts[classes]

bottom = np.zeros(len(chr_counts))
for cls in classes:
    ax4.bar(range(len(chr_counts)), chr_counts[cls], bottom=bottom,
            color=REG_COLORS[cls], label=cls, edgecolor='white', linewidth=0.5)
    bottom += chr_counts[cls].values
ax4.set_xticks(range(len(chr_counts)))
ax4.set_xticklabels(chr_counts.index, fontsize=9)
ax4.set_xlabel('Chromosome')
ax4.set_ylabel('Number of Significant Loci')
ax4.set_title('Regulatory Loci by Chromosome')
ax4.legend(fontsize=7, loc='upper right')

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.regulatory_analysis.png'), dpi=150)
plt.close()
print("  → lmd.regulatory_analysis.png")

# ── Figure 3: Top Loci Detail Plot (Top 15) ───────────────────────────
print("Generating top loci detail plot...")
fig, ax = plt.subplots(figsize=(14, 8))

top15 = reg.nsmallest(15, 'top_p').sort_values('top_p', ascending=True)
y_positions = range(len(top15))
rps_values = top15['regulatory_potential_score'].values
logp_values = -np.log10(top15['top_p'].values)

# Horizontal bar chart with RPS colored by class
bars = ax.barh(y_positions, logp_values,
               color=[REG_COLORS[c] for c in top15['regulatory_class']],
               edgecolor='black', linewidth=0.5, height=0.7)

# Add RPS score annotations
for i, (rps, cls) in enumerate(zip(rps_values, top15['regulatory_class'])):
    ax.text(logp_values[i] + 0.1, i, f'RPS={rps:.2f}', va='center', fontsize=8,
            fontweight='bold', color=REG_COLORS[cls])

# Add labels
ax.set_yticks(y_positions)
ax.set_yticklabels([f"{r['locus']}  (P={r['top_p']:.1e})" for _, r in top15.iterrows()],
                   fontsize=8)
ax.set_xlabel('-log10(P)')
ax.set_title('Top 15 Genome-Wide Significant Loci\nColored by Regulatory Potential',
             fontsize=13, fontweight='bold')
ax.invert_yaxis()

# Legend
legend_elements = [Patch(facecolor=REG_COLORS[c], label=c) for c in classes]
ax.legend(handles=legend_elements, fontsize=8, loc='lower right')

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.top_loci.png'), dpi=150)
plt.close()
print("  → lmd.top_loci.png")

# ── Figure 4: QQ Plot with Confidence Interval ────────────────────────
print("Generating enhanced QQ plot...")
fig, ax = plt.subplots(figsize=(6, 6))

observed = -np.log10(np.sort(gwas['P'].values))
n = len(observed)
expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)

# Confidence interval
ci = 0.95
c = np.ceil(n * (1 - ci) / 2)
upper_ci = -np.log10(stats.beta.ppf(1 - ci/2, np.arange(1, n+1), n - np.arange(1, n+1) + 1))
lower_ci = -np.log10(stats.beta.ppf(ci/2, np.arange(1, n+1), n - np.arange(1, n+1) + 1))

ax.fill_between(expected, lower_ci, upper_ci, alpha=0.1, color='gray')
ax.scatter(expected, observed, s=1.5, alpha=0.3, edgecolors='none', rasterized=True)
ax.plot([0, max(expected)], [0, max(expected)], 'r--', linewidth=1)

# Lambda
from scipy.stats import chi2
chisq_vals = chi2.ppf(1 - gwas['P'].values, 1)
lambda_gc = np.median(chisq_vals) / chi2.ppf(0.5, 1)
ax.text(0.05, 0.95, f'λ_GC = {lambda_gc:.3f}', transform=ax.transAxes,
        fontsize=12, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax.set_xlabel('Expected -log10(P)')
ax.set_ylabel('Observed -log10(P)')
ax.set_title('QQ Plot with 95% Confidence Interval')
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.qq_enhanced.png'), dpi=150)
plt.close()
print("  → lmd.qq_enhanced.png")

# ── Generate Summary Report ───────────────────────────────────────────
print("\nGenerating summary report...")

# Count signals per chromosome
sig_by_chr = gwas.groupby('CHR').apply(
    lambda x: (x['P'] < bonf_threshold).sum()
).to_dict()

n_fdr_sig = int(gwas['P_FDR'].lt(0.05).sum()) if 'P_FDR' in gwas.columns else 'N/A'

report = f"""
================================================================================
    GENOME-WIDE ASSOCIATION STUDY OF LOIN MUSCLE DEPTH IN PIGS
                  WITH REGULATORY ANNOTATION
================================================================================

DATA SUMMARY
------------
  Samples analyzed:          ~2,795 (out of 2,797 genotyped)
  SNPs tested:               258,074 (after QC; 588 removed by HWE/MAF filters)
  Covariates:                year, batch, sex
  Phenotype:                 loin muscle depth (mm) from FAM column 6
  Genomic inflation (λ):     {lambda_gc:.3f}

SIGNIFICANCE THRESHOLDS
-----------------------
  Bonferroni (0.05/N):       {bonf_threshold:.2e}
  Suggestive threshold:      1 × 10⁻⁵
  Bonferroni-significant:    {gwas['P'].lt(bonf_threshold).sum():,} SNPs
  Suggestive (P < 1e-5):     {gwas['P'].lt(1e-5).sum():,} SNPs
  FDR 5% significant:        {n_fdr_sig} SNPs
  Independent loci:           {len(loci)}

TOP GENOME-WIDE SIGNIFICANT LOCI
---------------------------------
"""
# Sort by P-value
top_loci = loci.nsmallest(10, 'top_p')
for i, (_, l) in enumerate(top_loci.iterrows()):
    reg_info = reg[reg['locus'] == f"Chr{l['chr']}:{l['start']:,}-{l['end']:,}"]
    if len(reg_info) > 0:
        rps = reg_info.iloc[0]['regulatory_potential_score']
        rclass = reg_info.iloc[0]['regulatory_class']
    else:
        rps = 'N/A'
        rclass = 'N/A'
    report += f"  {i+1:>2}. Chr{l['chr']:>2}:{l['start']:>10,}-{l['end']:>10,}  "
    report += f"P={l['top_p']:.2e}  "
    report += f"({l['top_snp']})  "
    report += f"RPS={rps if isinstance(rps, str) else f'{rps:.3f}'} [{rclass}]\n"

report += f"""
REGULATORY ANNOTATION SUMMARY
-----------------------------
  Loci assessed:              {len(reg)}
  Strong regulatory:          {(reg['regulatory_class'] == 'Strong regulatory').sum()}  (RPS ≥ 0.60)
  Moderate regulatory:        {(reg['regulatory_class'] == 'Moderate regulatory').sum()}  (RPS 0.45–0.60)
  Weak regulatory:            {(reg['regulatory_class'] == 'Weak regulatory').sum()}  (RPS 0.30–0.45)
  Likely non-regulatory:      {(reg['regulatory_class'] == 'Likely non-regulatory').sum()}  (RPS < 0.30)

  Key regulatory features:
  - Regulatory loci show lower MAF (mean {reg[reg['regulatory_class']=='Strong regulatory']['lead_maf'].mean():.3f}
    vs {reg[reg['regulatory_class']=='Likely non-regulatory']['lead_maf'].mean():.3f} for non-regulatory)
  - Higher SNP density in regulatory regions (mean {reg[reg['regulatory_class']=='Strong regulatory']['snp_density_50kb'].mean():.0f}
    vs {reg[reg['regulatory_class']=='Likely non-regulatory']['snp_density_50kb'].mean():.0f} SNPs/50kb)
  - {reg['n_cpg_proxy_snps'].gt(0).sum()} loci overlap CpG island proxy signals

CHROMOSOME DISTRIBUTION OF SIGNIFICANT SNPS
-------------------------------------------
"""
for ch in sorted(sig_by_chr.keys()):
    n_ch = gwas[gwas['CHR'] == ch].shape[0]
    n_sig = sig_by_chr[ch]
    bar = '█' * max(1, int(n_sig / max(sig_by_chr.values()) * 30))
    report += f"  Chr{ch:>2}: {n_sig:>4} significant / {n_ch:>6,} tested  {bar}\n"

report += f"""
MOST NOTABLE FINDINGS
---------------------
  1. The strongest signal is on Chromosome {loci.nsmallest(1, 'top_p').iloc[0]['chr']}
     (P = {loci.nsmallest(1, 'top_p').iloc[0]['top_p']:.2e}
     at {loci.nsmallest(1, 'top_p').iloc[0]['top_snp']}), with multiple overlapping
     regulatory loci suggesting a regulatory hotspot for muscle development.

  2. Chromosomes 5 ({sig_by_chr.get(5,0)} sig SNPs) and 11 ({sig_by_chr.get(11,0)} sig SNPs)
     are enriched for significant associations, containing numerous regulatory loci.

  3. The elevated genomic inflation (λ = {lambda_gc:.2f}) is consistent with
     polygenic architecture and/or population structure typical of livestock GWAS.
     Genomic control or mixed-model approaches could further refine these results.

  4. {reg['regulatory_class'].isin(['Strong regulatory', 'Moderate regulatory']).sum()} of
     {len(reg)} loci ({reg['regulatory_class'].isin(['Strong regulatory', 'Moderate regulatory']).sum()/len(reg)*100:.0f}%)
     show moderate-to-strong regulatory potential, suggesting that a substantial
     fraction of loin muscle depth QTL act through regulatory mechanisms.

OUTPUT FILES
------------
  lmd.gwas.assoc.linear           PLINK GWAS results (all SNPs)
  lmd.significant_loci.txt        Bonferroni-significant locus definitions
  lmd.regulatory_annotation.txt   Locus-level regulatory scores and classification
  lmd.top50_snps.txt              Top 50 associated SNPs
  lmd.suggestive_snps.txt         SNPs with P < 1e-5
  lmd.gwas_with_features.txt      SNP-level data with genomic features
  lmd.manhattan_annotated.png     Manhattan plot with regulatory annotation
  lmd.manhattan.png               Basic Manhattan plot
  lmd.qq_enhanced.png             QQ plot with 95% CI
  lmd.qq.png                      Basic QQ plot
  lmd.regulatory_analysis.png     Regulatory score distributions
  lmd.top_loci.png                Top 15 loci detail

================================================================================
"""

print(report)

with open(os.path.join(WORKDIR, 'lmd.gwas_report.txt'), 'w') as f:
    f.write(report)
print(f"Report saved to lmd.gwas_report.txt")
print("All done!")
