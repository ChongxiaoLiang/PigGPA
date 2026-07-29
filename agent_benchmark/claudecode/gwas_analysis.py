#!/usr/bin/env python3
"""
GWAS post-processing: genomic inflation, multiple-testing correction,
locus definition, Manhattan plot, QQ plot, and regulatory annotation.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

WORKDIR = '/workspace/pigbole/benchmark/agent_benchmark/ONE2/claudecode'
GWAS_FILE = os.path.join(WORKDIR, 'lmd.gwas.assoc.linear')

# ── 1. Load GWAS results ──────────────────────────────────────────────
print("Loading GWAS results...")
gwas = pd.read_csv(GWAS_FILE, delim_whitespace=True)
print(f"  {len(gwas):,} tests loaded")

# Remove NA p-values
gwas = gwas[gwas['P'].notna()].copy()
print(f"  {len(gwas):,} tests with non-NA p-values")

# ── 2. Genomic inflation factor (lambda) ──────────────────────────────
chisq = stats.chi2.ppf(1 - gwas['P'].values, 1)
lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)
print(f"\nGenomic inflation factor (lambda_GC): {lambda_gc:.4f}")

# ── 3. Multiple testing correction ────────────────────────────────────
n_tests = len(gwas)
bonf_threshold = 0.05 / n_tests
print(f"\nBonferroni threshold (0.05/{n_tests}): {bonf_threshold:.2e}")

# FDR correction (Benjamini-Hochberg)
reject_fdr, pvals_fdr, _, _ = multipletests(gwas['P'].values, alpha=0.05, method='fdr_bh')
gwas['P_FDR'] = pvals_fdr
gwas['significant_bonf'] = gwas['P'] < bonf_threshold
gwas['significant_fdr'] = reject_fdr

n_bonf = gwas['significant_bonf'].sum()
n_fdr = gwas['significant_fdr'].sum()
print(f"  Significant SNPs (Bonferroni): {n_bonf}")
print(f"  Significant SNPs (FDR 5%):    {n_fdr}")

# ── 4. Suggestive threshold ───────────────────────────────────────────
suggestive = 1e-5
n_suggestive = (gwas['P'] < suggestive).sum()
print(f"  Suggestive SNPs (P < 1e-5):   {n_suggestive}")

# ── 5. Define loci ────────────────────────────────────────────────────
# Use a simple window-based clumping: group significant SNPs within 500 kb
print("\n--- Defining genomic loci ---")
WINDOW_BP = 500_000  # 500 kb window

sig_snps = gwas[gwas['significant_bonf']].sort_values(['CHR', 'BP']).copy()

loci = []
if len(sig_snps) > 0:
    current_locus = {
        'chr': sig_snps.iloc[0]['CHR'],
        'start': sig_snps.iloc[0]['BP'],
        'end': sig_snps.iloc[0]['BP'],
        'snps': [sig_snps.iloc[0]['SNP']],
        'top_snp': sig_snps.iloc[0]['SNP'],
        'top_p': sig_snps.iloc[0]['P'],
        'top_bp': sig_snps.iloc[0]['BP'],
    }
    for _, row in sig_snps.iloc[1:].iterrows():
        if row['CHR'] == current_locus['chr'] and row['BP'] - current_locus['end'] <= WINDOW_BP:
            current_locus['end'] = row['BP']
            current_locus['snps'].append(row['SNP'])
            if row['P'] < current_locus['top_p']:
                current_locus['top_snp'] = row['SNP']
                current_locus['top_p'] = row['P']
                current_locus['top_bp'] = row['BP']
        else:
            loci.append(current_locus)
            current_locus = {
                'chr': row['CHR'],
                'start': row['BP'],
                'end': row['BP'],
                'snps': [row['SNP']],
                'top_snp': row['SNP'],
                'top_p': row['P'],
                'top_bp': row['BP'],
            }
    loci.append(current_locus)

print(f"Identified {len(loci)} independent genome-wide significant loci:\n")
for i, locus in enumerate(loci):
    print(f"  Locus {i+1}: Chr{locus['chr']}:{locus['start']:,}-{locus['end']:,} "
          f"({len(locus['snps'])} SNPs, top: {locus['top_snp']} P={locus['top_p']:.2e})")

# Save loci to file
loci_df = pd.DataFrame(loci)
loci_df.to_csv(os.path.join(WORKDIR, 'lmd.significant_loci.txt'), sep='\t', index=False)
print(f"\nLoci saved to lmd.significant_loci.txt")

# Also save suggestive loci (P < 1e-5) for broader annotation
suggestive_snps = gwas[gwas['P'] < suggestive].copy()
suggestive_snps.to_csv(os.path.join(WORKDIR, 'lmd.suggestive_snps.txt'), sep='\t', index=False)
print(f"Suggestive SNPs saved to lmd.suggestive_snps.txt ({len(suggestive_snps)} SNPs)")

# ── 6. Manhattan Plot ─────────────────────────────────────────────────
print("\nGenerating Manhattan plot...")
gwas['logP'] = -np.log10(gwas['P'])

# Assign cumulative position for x-axis
chroms = sorted(gwas['CHR'].unique())
gwas['CHR_IDX'] = gwas['CHR'].map({c: i for i, c in enumerate(chroms)})
gwas = gwas.sort_values(['CHR_IDX', 'BP'])

chr_offsets = {}
cumulative_pos = 0
chr_midpoints = {}
for ch in chroms:
    ch_data = gwas[gwas['CHR'] == ch]
    chr_offsets[ch] = cumulative_pos
    chr_midpoints[ch] = cumulative_pos + ch_data['BP'].max() / 2
    cumulative_pos += ch_data['BP'].max()

gwas['CUM_POS'] = gwas['BP'] + gwas['CHR'].map(chr_offsets)

fig, ax = plt.subplots(figsize=(16, 6))

colors = ['#1f77b4', '#ff7f0e']
for i, ch in enumerate(chroms):
    ch_data = gwas[gwas['CHR'] == ch]
    ax.scatter(ch_data['CUM_POS'], ch_data['logP'], s=2, c=colors[i % 2],
               alpha=0.6, edgecolors='none')

# Significance thresholds
ax.axhline(-np.log10(bonf_threshold), color='red', linestyle='--', linewidth=1,
           label=f'Bonferroni ({bonf_threshold:.1e})')
ax.axhline(-np.log10(suggestive), color='blue', linestyle='--', linewidth=0.8,
           label=f'Suggestive (1e-5)')

# Highlight top loci
for locus in loci:
    top_ch = locus['chr']
    top_pos = gwas[(gwas['CHR'] == top_ch) & (gwas['BP'] == locus['top_bp'])]
    if len(top_pos) > 0:
        ax.scatter(top_pos['CUM_POS'], top_pos['logP'], s=80, c='red',
                   edgecolors='darkred', linewidth=1, zorder=5)

ax.set_xlabel('Chromosome')
ax.set_ylabel('-log10(P)')
ax.set_title(f'GWAS: Loin Muscle Depth\n(lambda_GC = {lambda_gc:.3f}, {n_bonf} Bonferroni-significant SNPs)')
ax.set_xticks([chr_midpoints[ch] for ch in chroms])
ax.set_xticklabels([str(ch) for ch in chroms], fontsize=8)
ax.set_xlim(0, cumulative_pos)
ax.legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.manhattan.png'), dpi=150)
plt.close()
print("  Manhattan plot saved to lmd.manhattan.png")

# ── 7. QQ Plot ────────────────────────────────────────────────────────
print("Generating QQ plot...")
fig, ax = plt.subplots(figsize=(6, 6))

observed = -np.log10(np.sort(gwas['P'].values))
expected = -np.log10((np.arange(1, len(observed) + 1)) / (len(observed) + 1))

ax.scatter(expected, observed, s=2, alpha=0.4, edgecolors='none')
ax.plot([0, max(expected)], [0, max(expected)], 'r--', linewidth=1)
ax.set_xlabel('Expected -log10(P)')
ax.set_ylabel('Observed -log10(P)')
ax.set_title(f'QQ Plot (lambda_GC = {lambda_gc:.3f})')
ax.set_aspect('equal')

# Add lambda annotation
ax.text(0.05, 0.95, f'λ_GC = {lambda_gc:.4f}', transform=ax.transAxes,
        fontsize=11, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'lmd.qq.png'), dpi=150)
plt.close()
print("  QQ plot saved to lmd.qq.png")

# ── 8. Summary table of top SNPs ──────────────────────────────────────
top_n = 50
top_snps = gwas.nsmallest(top_n, 'P')[
    ['CHR', 'SNP', 'BP', 'A1', 'NMISS', 'BETA', 'STAT', 'P', 'P_FDR', 'significant_bonf']
].copy()
top_snps['P'] = top_snps['P'].apply(lambda x: f'{x:.3e}')
top_snps['P_FDR'] = top_snps['P_FDR'].apply(lambda x: f'{x:.3e}')
top_snps.to_csv(os.path.join(WORKDIR, 'lmd.top50_snps.txt'), sep='\t', index=False)
print(f"\nTop {top_n} SNPs saved to lmd.top50_snps.txt")

# ── 9. Summary statistics ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("GWAS SUMMARY")
print("=" * 70)
print(f"  Samples:                    {gwas['NMISS'].max()}")
print(f"  SNPs tested:                {n_tests:,}")
print(f"  Genomic inflation (λ):      {lambda_gc:.4f}")
print(f"  Bonferroni threshold:       {bonf_threshold:.2e}")
print(f"  Bonferroni-significant:     {n_bonf}")
print(f"  FDR-significant (5%):       {n_fdr}")
print(f"  Suggestive (P < 1e-5):      {n_suggestive}")
print(f"  Independent loci:           {len(loci)}")
print()

if len(loci) > 0:
    print("GENOME-WIDE SIGNIFICANT LOCI:")
    print("-" * 70)
    for i, locus in enumerate(loci):
        print(f"  {i+1}. Chr{locus['chr']}:{locus['start']:,}-{locus['end']:,} "
              f"[{locus['start']:,}-{locus['end']:,} bp]")
        print(f"      {len(locus['snps'])} SNPs, top={locus['top_snp']} P={locus['top_p']:.2e}")

# ── 10. Per-chromosome summary ────────────────────────────────────────
print("\nPER-CHROMOSOME SIGNIFICANT SNP COUNTS:")
for ch in chroms:
    n_ch = (gwas['CHR'] == ch).sum()
    n_sig = ((gwas['CHR'] == ch) & gwas['significant_bonf']).sum()
    print(f"  Chr{ch:>2}: {n_ch:>7,} SNPs, {n_sig:>4} significant")

print("\nDone! All output files are in:", WORKDIR)
