#!/usr/bin/env python3
import os
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import eigh
from pandas_plink import read_plink
import warnings
warnings.filterwarnings('ignore')

BFILE = "/public/share/likui/hanyu/testdata/In-silico-data/simulated_population"
KEEP_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/keep_1000_samples.txt"
SNP_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/chr1_snps.txt"
BASE_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t6/1"
NPC = 10

NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377',
                 '#BBCC33', '#CC3311', '#EE7733', '#999999']

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12


def read_genotype(bfile, keep_file, snp_file):
    keep_ids = pd.read_csv(keep_file, header=None)[0].astype(str).tolist()
    snp_ids = pd.read_csv(snp_file, header=None)[0].tolist()
    bim, fam, G_dask = read_plink(bfile, verbose=False)
    snp_set = set(snp_ids)
    snp_idx = bim[bim['snp'].isin(snp_set)].index.tolist()
    sample_mask = fam['iid'].astype(str).isin(set(keep_ids))
    sample_indices = fam[sample_mask].index.tolist()
    G = G_dask[snp_idx, :][:, sample_indices].compute()
    G_out = G.T.copy()
    for j in range(G_out.shape[1]):
        valid = ~np.isnan(G_out[:, j])
        if np.sum(valid) > 0 and np.sum(~valid) > 0:
            G_out[~valid, j] = np.nanmean(G_out[valid, j])
    sample_iids = fam[sample_mask]['iid'].values.astype(int).tolist()
    return G_out, sample_iids, len(snp_idx)


def calc_grm(G):
    p = np.mean(G, axis=0) / 2
    Z = G - 2 * p
    sum_2pq = 2 * np.sum(p * (1 - p))
    GRM_vr = (Z @ Z.T) / sum_2pq
    diag_mean = np.mean(np.diag(GRM_vr))
    GRM = GRM_vr / diag_mean
    return GRM


def pca_from_grm(grm, n_components):
    eigenvalues, eigenvectors = eigh(grm)
    sorted_idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_idx]
    eigenvectors = eigenvectors[:, sorted_idx]

    eigenvalues_top = eigenvalues[:n_components]
    eigenvectors_top = eigenvectors[:, :n_components]

    total_var = np.sum(eigenvalues)
    std_dev = np.sqrt(np.maximum(eigenvalues_top, 0))
    prop_var = eigenvalues_top / total_var
    cum_var = np.cumsum(prop_var)

    pc_scores = eigenvectors_top * np.sqrt(np.maximum(eigenvalues_top, 0))

    return eigenvalues_top, eigenvectors_top, std_dev, prop_var, cum_var, pc_scores


def save_pc(filepath, sample_ids, pc_scores, n_components):
    with open(filepath, 'w') as f:
        header = "id\t" + "\t".join([f"PC{i+1}" for i in range(n_components)])
        f.write(header + "\n")
        for i, sid in enumerate(sample_ids):
            vals = "\t".join(f"{pc_scores[i, j]:.6g}" for j in range(n_components))
            f.write(f"{sid}\t{vals}\n")


def save_pcp(filepath, std_dev, prop_var, cum_var, n_components):
    with open(filepath, 'w') as f:
        header = "Components\t" + "\t".join([f"PC{i+1}" for i in range(n_components)])
        f.write(header + "\n")
        f.write("Standard deviation\t" + "\t".join(f"{std_dev[i]:.5g}" for i in range(n_components)) + "\n")
        f.write("Proportion of Variance\t" + "\t".join(f"{prop_var[i]:.6g}" for i in range(n_components)) + "\n")
        f.write("Cumulative Proportion\t" + "\t".join(f"{cum_var[i]:.7g}" for i in range(n_components)) + "\n")


def save_log(filepath, start_time, n_samples, n_snps, n_components, std_dev, prop_var, cum_var):
    elapsed = time.time() - start_time
    with open(filepath, 'w') as f:
        f.write("PCA Analysis Log\n")
        f.write("=" * 60 + "\n")
        f.write(f"Analysis started: {time.strftime('%a %b %d %H:%M:%S %Y', time.localtime(start_time))}\n")
        f.write(f"\nModel: PCA from GRM (VanRaden + Su normalization)\n")
        f.write(f"Samples: {n_samples}, SNPs: {n_snps}, PCs: {n_components}\n")
        f.write(f"\nPrincipal Components analysis ...\n")
        f.write(f"{n_components} PCs have been calculated successfully.\n")
        f.write("Overview of the importance for the top PCs:\n")
        header = f"{'':>24s}" + "".join(f"{'PC'+str(i+1):>10s}" for i in range(min(5, n_components)))
        f.write(header + "\n")
        sd_line = f"{'Standard deviation':>24s}" + "".join(f"{std_dev[i]:10.4f}" for i in range(min(5, n_components)))
        f.write(sd_line + "\n")
        pv_line = f"{'Proportion of Variance':>24s}" + "".join(f"{prop_var[i]:10.4f}" for i in range(min(5, n_components)))
        f.write(pv_line + "\n")
        cp_line = f"{'Cumulative Proportion':>24s}" + "".join(f"{cum_var[i]:10.4f}" for i in range(min(5, n_components)))
        f.write(cp_line + "\n")
        f.write(f"\nAnalysis finished: {time.strftime('%a %b %d %H:%M:%S %Y')}\n")
        mins = int(elapsed) // 60
        secs = int(elapsed) % 60
        f.write(f"Total running time: 0h{mins}m{secs}s\n")


def plot_scree_plot(std_dev, prop_var, cum_var, output_file, n_pcs):
    n_pcs = min(n_pcs, len(std_dev))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('PCA Scree Plot and Variance Explained', fontsize=18, fontweight='bold')

    ax1 = axes[0]
    ax1.bar(range(1, n_pcs + 1), std_dev[:n_pcs] ** 2, color=NATURE_COLORS[0],
            edgecolor='black', alpha=0.8)
    ax1.plot(range(1, n_pcs + 1), std_dev[:n_pcs] ** 2, 'o-', color=NATURE_COLORS[1],
             linewidth=2, markersize=6)
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Eigenvalue')
    ax1.set_title('Eigenvalue Scree Plot')
    ax1.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax1)

    ax2 = axes[1]
    ax2.bar(range(1, n_pcs + 1), prop_var[:n_pcs] * 100,
            color=NATURE_COLORS[2], edgecolor='black', alpha=0.8)
    ax2.set_xlabel('Principal Component')
    ax2.set_ylabel('Variance Explained (%)')
    ax2.set_title('Individual Variance Explained')
    ax2.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax2)

    ax3 = axes[2]
    ax3.plot(range(1, n_pcs + 1), cum_var[:n_pcs] * 100,
             'o-', color=NATURE_COLORS[3], linewidth=2, markersize=6)
    ax3.axhline(y=80, color=NATURE_COLORS[1], linestyle='--', linewidth=1.5, label='80% threshold')
    ax3.set_xlabel('Principal Component')
    ax3.set_ylabel('Cumulative Variance (%)')
    ax3.set_title('Cumulative Variance Explained')
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax3.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax3)

    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_pc_scatter(pc_scores, pc_x, pc_y, var_x, var_y, output_file):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(pc_scores[:, pc_x - 1], pc_scores[:, pc_y - 1],
               c=NATURE_COLORS[0], alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
    ax.set_xlabel(f'PC{pc_x} ({var_x * 100:.2f}%)')
    ax.set_ylabel(f'PC{pc_y} ({var_y * 100:.2f}%)')
    ax.set_title(f'PCA: PC{pc_x} vs PC{pc_y}')
    ax.grid(True, alpha=0.15, linestyle='--')
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_pc_density(pc_scores, output_file, n_pcs=4):
    from scipy.stats import norm
    n_pcs = min(n_pcs, pc_scores.shape[1])
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Principal Component Density Distribution', fontsize=18, fontweight='bold')
    for idx in range(min(n_pcs, 4)):
        ax = axes[idx // 2, idx % 2]
        data = pc_scores[:, idx]
        ax.hist(data, bins=30, color=NATURE_COLORS[idx], edgecolor='black', alpha=0.7, density=True)
        mu, std = norm.fit(data)
        x = np.linspace(data.min(), data.max(), 100)
        ax.plot(x, norm.pdf(x, mu, std), 'k-', linewidth=2,
                label=f'Normal fit\nμ={mu:.2f}, σ={std:.2f}')
        ax.set_xlabel(f'PC{idx + 1} Score')
        ax.set_ylabel('Density')
        ax.set_title(f'PC{idx + 1} Distribution')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        ax.grid(True, alpha=0.15, linestyle='--')
        sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def main():
    os.makedirs(BASE_DIR, exist_ok=True)
    start_time = time.time()

    print("=" * 70)
    print("  PCA Analysis (matching HIBLUP)")
    print("  Method: GRM eigenvalue decomposition (VanRaden + Su normalization)")
    print("=" * 70)

    print("\n[1] Reading genotype data...")
    G_geno, sample_ids, n_snps = read_genotype(BFILE, KEEP_FILE, SNP_FILE)
    n = len(sample_ids)
    print(f"  Samples: {n}, SNPs: {n_snps}")

    print("\n[2] Computing GRM (VanRaden + Su normalization)...")
    grm = calc_grm(G_geno)
    print(f"  GRM dim: {grm.shape}, diag mean: {np.mean(np.diag(grm)):.6f}")

    print("\n[3] PCA from GRM eigenvalue decomposition...")
    eigenvalues, eigenvectors, std_dev, prop_var, cum_var, pc_scores = \
        pca_from_grm(grm, NPC)
    print(f"  Top 5 PC standard deviations: {std_dev[:5]}")
    print(f"  Top 5 PC variance explained: {prop_var[:5] * 100}%")
    print(f"  Cumulative variance (top 5): {cum_var[:5] * 100}%")

    print("\n[4] Saving output files (HIBLUP format)...")

    save_pc(os.path.join(BASE_DIR, "chr1_pca.pc"), sample_ids, pc_scores, NPC)
    print(f"  Saved: chr1_pca.pc")

    save_pcp(os.path.join(BASE_DIR, "chr1_pca.pcp"), std_dev, prop_var, cum_var, NPC)
    print(f"  Saved: chr1_pca.pcp")

    save_log(os.path.join(BASE_DIR, "chr1_pca.log"), start_time, n, n_snps, NPC,
             std_dev, prop_var, cum_var)
    print(f"  Saved: chr1_pca.log")

    print("\n[5] Generating visualization plots...")

    plot_scree_plot(std_dev, prop_var, cum_var,
                    os.path.join(BASE_DIR, "PCA_scree_plot.pdf"), NPC)

    plot_pc_scatter(pc_scores, 1, 2, prop_var[0], prop_var[1],
                    os.path.join(BASE_DIR, "PCA_PC1_vs_PC2.pdf"))

    plot_pc_scatter(pc_scores, 1, 3, prop_var[0], prop_var[2],
                    os.path.join(BASE_DIR, "PCA_PC1_vs_PC3.pdf"))

    plot_pc_scatter(pc_scores, 2, 3, prop_var[1], prop_var[2],
                    os.path.join(BASE_DIR, "PCA_PC2_vs_PC3.pdf"))

    plot_pc_density(pc_scores, os.path.join(BASE_DIR, "PCA_density.pdf"), n_pcs=4)

    print(f"\n  Results saved to: {BASE_DIR}")
    print("\n" + "=" * 70)
    print("  Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
