#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA遗传结构分析脚本
支持基于基因型或关系矩阵的主成分分析
遵循科研绘图规范：DejaVu Sans字体、PDF输出、Nature/Science配色
修改版：只抽取染色体1的SNP，使用前1000个样本
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from scipy.linalg import eigh, svd
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

try:
    from pandas_plink import read_plink
except ImportError:
    print("正在安装 pandas-plink...")
    os.system("pip install pandas-plink")
    from pandas_plink import read_plink

INPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data"
OUTPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t6"
BED_FILE = os.path.join(INPUT_DIR, "simulated_population")
PHENOTYPE_FILE = os.path.join(INPUT_DIR, "phenotypes.txt")
TARGET_CHROMOSOME = "1"
TARGET_SAMPLES = 1000

os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 18

NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377', '#BBCC33', '#CC3311', '#EE7733', '#999999']


class PCAAnalyzer:
    """PCA遗传结构分析器（复刻HIBLUP）"""
    
    def __init__(self):
        self.eigenvalues = None
        self.eigenvectors = None
        self.explained_variance_ratio = None
        self.cumulative_variance = None
        self.pc_scores = None
        self.n_components = None
        self.sample_ids = None
    
    def fit_from_genotype(self, genotype_matrix, sample_ids, n_components=20, 
                          method='svd', center=True, scale=True):
        print("从基因型数据进行PCA分析...")
        
        n_snps, n_samples = genotype_matrix.shape
        self.sample_ids = list(sample_ids)
        self.n_components = min(n_components, n_samples - 1)
        
        print(f"  样本数: {n_samples}")
        print(f"  SNP数: {n_snps}")
        print(f"  目标主成分数: {self.n_components}")
        
        print("  预处理基因型数据...")
        G = genotype_matrix.astype(np.float64)
        
        snp_means = np.nanmean(G, axis=1, keepdims=True)
        G_centered = np.where(np.isnan(G), 0, G - snp_means)
        
        if scale:
            snp_vars = np.nanvar(G, axis=1)
            snp_stds = np.sqrt(snp_vars + 1e-10)
            snp_stds = np.where(snp_stds < 1e-10, 1, snp_stds)
            G_scaled = G_centered / snp_stds[:, np.newaxis]
        else:
            G_scaled = G_centered
        
        valid_snps = ~np.all(G_scaled == 0, axis=1)
        G_filtered = G_scaled[valid_snps, :]
        print(f"  有效SNP数: {G_filtered.shape[0]}")
        
        if method == 'svd':
            print("  使用SVD方法计算主成分...")
            U, S, Vt = svd(G_filtered, full_matrices=False)
            
            self.eigenvalues = (S ** 2) / (n_samples - 1)
            self.eigenvectors = Vt.T
            
            total_variance = np.sum(self.eigenvalues)
            self.explained_variance_ratio = self.eigenvalues / total_variance
            self.cumulative_variance = np.cumsum(self.explained_variance_ratio)
            
            self.pc_scores = self.eigenvectors[:, :self.n_components] * np.sqrt(self.eigenvalues[:self.n_components])
        
        else:
            print("  计算样本间协方差矩阵...")
            cov_matrix = np.dot(G_filtered.T, G_filtered) / G_filtered.shape[0]
            
            print("  特征值分解...")
            eigenvalues, eigenvectors = eigh(cov_matrix)
            
            sorted_idx = np.argsort(eigenvalues)[::-1]
            self.eigenvalues = eigenvalues[sorted_idx]
            self.eigenvectors = eigenvectors[:, sorted_idx]
            
            total_variance = np.sum(self.eigenvalues)
            self.explained_variance_ratio = self.eigenvalues / total_variance
            self.cumulative_variance = np.cumsum(self.explained_variance_ratio)
            
            self.pc_scores = self.eigenvectors[:, :self.n_components] * np.sqrt(self.eigenvalues[:self.n_components])
        
        print(f"  前10个主成分方差解释比例: {self.explained_variance_ratio[:10]}")
        
        return self
    
    def fit_from_grm(self, grm, sample_ids, n_components=20):
        print("从基因组关系矩阵进行PCA分析...")
        
        n_samples = grm.shape[0]
        self.sample_ids = list(sample_ids)
        self.n_components = min(n_components, n_samples - 1)
        
        print(f"  样本数: {n_samples}")
        print(f"  目标主成分数: {self.n_components}")
        
        print("  特征值分解...")
        eigenvalues, eigenvectors = eigh(grm)
        
        sorted_idx = np.argsort(eigenvalues)[::-1]
        self.eigenvalues = eigenvalues[sorted_idx]
        self.eigenvectors = eigenvectors[:, sorted_idx]
        
        total_variance = np.sum(np.abs(self.eigenvalues))
        self.explained_variance_ratio = np.abs(self.eigenvalues) / total_variance
        self.cumulative_variance = np.cumsum(self.explained_variance_ratio)
        
        self.pc_scores = self.eigenvectors[:, :self.n_components] * np.sqrt(np.abs(self.eigenvalues[:self.n_components]))
        
        print(f"  前10个主成分方差解释比例: {self.explained_variance_ratio[:10]}")
        
        return self
    
    def get_pc_scores(self, n_pcs=None):
        if n_pcs is None:
            n_pcs = self.n_components
        
        pc_df = pd.DataFrame(
            self.pc_scores[:, :n_pcs],
            index=self.sample_ids,
            columns=[f'PC{i+1}' for i in range(n_pcs)]
        )
        return pc_df
    
    def get_variance_explained(self, n_pcs=None):
        if n_pcs is None:
            n_pcs = self.n_components
        
        var_df = pd.DataFrame({
            'PC': [f'PC{i+1}' for i in range(n_pcs)],
            'Eigenvalue': self.eigenvalues[:n_pcs],
            'Variance_Explained': self.explained_variance_ratio[:n_pcs],
            'Cumulative_Variance': self.cumulative_variance[:n_pcs]
        })
        return var_df
    
    def transform_to_pc_space(self, pc_x=1, pc_y=2):
        return self.pc_scores[:, pc_x-1], self.pc_scores[:, pc_y-1]
    
    def calculate_pc_correlations(self, covariates):
        correlations = {}
        
        for pc_idx in range(min(self.n_components, 10)):
            pc_name = f'PC{pc_idx + 1}'
            correlations[pc_name] = {}
            
            for cov_name, cov_values in covariates.items():
                cov_arr = np.array(cov_values)
                pc_values = self.pc_scores[:, pc_idx]
                
                valid_mask = ~np.isnan(cov_arr) & ~np.isnan(pc_values)
                if np.sum(valid_mask) > 2:
                    corr, pval = pearsonr(pc_values[valid_mask], cov_arr[valid_mask])
                    correlations[pc_name][cov_name] = {'correlation': corr, 'p_value': pval}
                else:
                    correlations[pc_name][cov_name] = {'correlation': np.nan, 'p_value': np.nan}
        
        return correlations


def calculate_grm(genotype_matrix, sample_ids):
    n_snps, n_samples = genotype_matrix.shape
    
    print("计算基因组关系矩阵...")
    grm = np.zeros((n_samples, n_samples))
    snp_count = 0
    
    batch_size = 500
    total_batches = (n_snps + batch_size - 1) // batch_size
    
    for batch_idx in range(total_batches):
        start = batch_idx * batch_size
        end = min((batch_idx + 1) * batch_size, n_snps)
        
        if (batch_idx + 1) % 5 == 0 or batch_idx == total_batches - 1:
            print(f"  进度: {end}/{n_snps} ({100*end/n_snps:.1f}%)")
        
        batch = genotype_matrix[start:end, :]
        
        for i in range(batch.shape[0]):
            g = batch[i, :]
            valid_mask = ~np.isnan(g)
            g_valid = g[valid_mask]
            
            if len(g_valid) < 10:
                continue
            
            p = np.nanmean(g) / 2
            if p == 0 or p == 1 or p < 1e-10 or p > 1 - 1e-10:
                continue
            
            g_centered = np.zeros(n_samples)
            g_centered[valid_mask] = (g_valid - 2 * p) / np.sqrt(2 * p * (1 - p))
            
            grm += np.outer(g_centered, g_centered)
            snp_count += 1
    
    if snp_count > 0:
        grm = grm / snp_count
    
    print(f"  有效SNP数: {snp_count}")
    
    return grm


def plot_scree_plot(eigenvalues, explained_variance, cumulative_variance, output_file, n_pcs=20):
    n_pcs = min(n_pcs, len(eigenvalues))
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('PCA Scree Plot and Variance Explained', fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')
    
    ax1 = axes[0]
    ax1.bar(range(1, n_pcs + 1), eigenvalues[:n_pcs], color=NATURE_COLORS[0], edgecolor='black', alpha=0.8)
    ax1.plot(range(1, n_pcs + 1), eigenvalues[:n_pcs], 'o-', color=NATURE_COLORS[1], linewidth=2, markersize=6)
    ax1.set_xlabel('Principal Component', fontsize=14)
    ax1.set_ylabel('Eigenvalue', fontsize=14)
    ax1.set_title('Eigenvalue Scree Plot', fontsize=16)
    ax1.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax1)
    
    ax2 = axes[1]
    ax2.bar(range(1, n_pcs + 1), explained_variance[:n_pcs] * 100, 
            color=NATURE_COLORS[2], edgecolor='black', alpha=0.8)
    ax2.set_xlabel('Principal Component', fontsize=14)
    ax2.set_ylabel('Variance Explained (%)', fontsize=14)
    ax2.set_title('Individual Variance Explained', fontsize=16)
    ax2.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax2)
    
    ax3 = axes[2]
    ax3.plot(range(1, n_pcs + 1), cumulative_variance[:n_pcs] * 100, 
             'o-', color=NATURE_COLORS[3], linewidth=2, markersize=6)
    ax3.axhline(y=80, color=NATURE_COLORS[1], linestyle='--', linewidth=1.5, label='80% threshold')
    ax3.axhline(y=90, color=NATURE_COLORS[4], linestyle='--', linewidth=1.5, label='90% threshold')
    ax3.set_xlabel('Principal Component', fontsize=14)
    ax3.set_ylabel('Cumulative Variance (%)', fontsize=14)
    ax3.set_title('Cumulative Variance Explained', fontsize=16)
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax3.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax3)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  碎石图已保存: {output_file}")


def plot_pc_scatter(pc_scores, sample_ids, pc_x=1, pc_y=2, 
                    group_labels=None, output_file=None, 
                    var_explained_x=None, var_explained_y=None):
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    pc_x_values = pc_scores[:, pc_x - 1]
    pc_y_values = pc_scores[:, pc_y - 1]
    
    if group_labels is not None:
        unique_groups = list(set(group_labels))
        colors = [NATURE_COLORS[i % len(NATURE_COLORS)] for i in range(len(unique_groups))]
        group_colors = dict(zip(unique_groups, colors))
        
        for group in unique_groups:
            mask = np.array([g == group for g in group_labels])
            ax.scatter(pc_x_values[mask], pc_y_values[mask], 
                      c=group_colors[group], label=group, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
        
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, title='Group')
    else:
        ax.scatter(pc_x_values, pc_y_values, c=NATURE_COLORS[0], alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
    
    xlabel = f'PC{pc_x}'
    ylabel = f'PC{pc_y}'
    if var_explained_x is not None:
        xlabel += f' ({var_explained_x * 100:.2f}%)'
    if var_explained_y is not None:
        ylabel += f' ({var_explained_y * 100:.2f}%)'
    
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(f'PCA: PC{pc_x} vs PC{pc_y}', fontsize=16)
    ax.grid(True, alpha=0.15, linestyle='--')
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
        print(f"  散点图已保存: {output_file}")
    
    plt.close()


def plot_pc_density(pc_scores, output_file, n_pcs=4):
    n_pcs = min(n_pcs, pc_scores.shape[1])
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Principal Component Density Distribution', fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')
    
    for idx in range(min(n_pcs, 4)):
        ax = axes[idx // 2, idx % 2]
        
        data = pc_scores[:, idx]
        ax.hist(data, bins=30, color=NATURE_COLORS[idx], edgecolor='black', alpha=0.7, density=True)
        
        from scipy.stats import norm
        mu, std = norm.fit(data)
        x = np.linspace(data.min(), data.max(), 100)
        ax.plot(x, norm.pdf(x, mu, std), 'k-', linewidth=2, label=f'Normal fit\nμ={mu:.2f}, σ={std:.2f}')
        
        ax.set_xlabel(f'PC{idx + 1} Score', fontsize=14)
        ax.set_ylabel('Density', fontsize=14)
        ax.set_title(f'PC{idx + 1} Distribution', fontsize=16)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        ax.grid(True, alpha=0.15, linestyle='--')
        sns.despine(ax=ax)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  密度分布图已保存: {output_file}")


def main():
    print("=" * 70)
    print("PCA遗传结构分析（复刻HIBLUP）")
    print(f"只分析染色体 {TARGET_CHROMOSOME} 的SNP")
    print(f"只使用前 {TARGET_SAMPLES} 个样本")
    print("=" * 70)
    
    print("\n[1/6] 读取PLINK文件...")
    bim, fam, bed = read_plink(BED_FILE, verbose=False)
    
    chrom_mask = bim['chrom'].astype(str) == TARGET_CHROMOSOME
    chr1_snp_indices = np.where(chrom_mask)[0]
    
    n_total_samples = len(fam)
    n_samples = min(TARGET_SAMPLES, n_total_samples)
    sample_ids = list(fam['iid'].values[:n_samples])
    
    print(f"  染色体{TARGET_CHROMOSOME}的SNP数量: {len(chr1_snp_indices)}")
    print(f"  总样本数量: {n_total_samples}")
    print(f"  使用的样本数量: {n_samples}")
    
    print("\n[2/6] 加载基因型数据...")
    genotype_matrix = bed[chr1_snp_indices, :n_samples].compute()
    
    print("\n[3/6] 进行PCA分析...")
    pca = PCAAnalyzer()
    pca.fit_from_genotype(genotype_matrix, sample_ids, n_components=20, method='svd')
    
    print("\n[4/6] 保存结果...")
    
    output_pc_scores = os.path.join(OUTPUT_DIR, "PCA_scores.csv")
    pc_scores_df = pca.get_pc_scores()
    pc_scores_df.to_csv(output_pc_scores)
    print(f"  主成分得分已保存: {output_pc_scores}")
    
    output_variance = os.path.join(OUTPUT_DIR, "PCA_variance_explained.csv")
    variance_df = pca.get_variance_explained()
    variance_df.to_csv(output_variance, index=False)
    print(f"  方差解释信息已保存: {output_variance}")
    
    output_eigenvalues = os.path.join(OUTPUT_DIR, "PCA_eigenvalues.csv")
    eigen_df = pd.DataFrame({
        'PC': [f'PC{i+1}' for i in range(len(pca.eigenvalues))],
        'Eigenvalue': pca.eigenvalues,
        'Variance_Explained': pca.explained_variance_ratio,
        'Cumulative_Variance': pca.cumulative_variance
    })
    eigen_df.to_csv(output_eigenvalues, index=False)
    print(f"  特征值已保存: {output_eigenvalues}")
    
    output_summary = os.path.join(OUTPUT_DIR, "PCA_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("PCA遗传结构分析摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"分析染色体: {TARGET_CHROMOSOME}\n")
        f.write(f"总样本数: {n_samples}\n")
        f.write(f"总SNP数: {len(chr1_snp_indices)}\n")
        f.write(f"计算主成分数: {pca.n_components}\n\n")
        
        f.write("方差解释统计:\n")
        f.write(f"  PC1解释方差: {pca.explained_variance_ratio[0] * 100:.2f}%\n")
        f.write(f"  PC2解释方差: {pca.explained_variance_ratio[1] * 100:.2f}%\n")
        f.write(f"  PC1-PC2累计: {pca.cumulative_variance[1] * 100:.2f}%\n")
        f.write(f"  前5个PC累计: {pca.cumulative_variance[4] * 100:.2f}%\n")
        f.write(f"  前10个PC累计: {pca.cumulative_variance[9] * 100:.2f}%\n\n")
        
        f.write("特征值统计:\n")
        f.write(f"  最大特征值: {pca.eigenvalues[0]:.4f}\n")
        f.write(f"  最小特征值: {pca.eigenvalues[-1]:.4f}\n\n")
        
        f.write("主成分得分统计:\n")
        for i in range(min(5, pca.n_components)):
            pc_scores = pca.pc_scores[:, i]
            f.write(f"  PC{i+1}: 均值={np.mean(pc_scores):.4f}, 标准差={np.std(pc_scores):.4f}\n")
    print(f"  摘要已保存: {output_summary}")
    
    print("\n[5/6] 从GRM进行PCA分析（对比验证）...")
    grm = calculate_grm(genotype_matrix, sample_ids)
    
    pca_grm = PCAAnalyzer()
    pca_grm.fit_from_grm(grm, sample_ids, n_components=20)
    
    output_grm = os.path.join(OUTPUT_DIR, "GRM.csv")
    pd.DataFrame(grm, index=sample_ids, columns=sample_ids).to_csv(output_grm)
    print(f"  GRM已保存: {output_grm}")
    
    output_pc_scores_grm = os.path.join(OUTPUT_DIR, "PCA_scores_from_GRM.csv")
    pca_grm.get_pc_scores().to_csv(output_pc_scores_grm)
    print(f"  GRM来源的主成分得分已保存: {output_pc_scores_grm}")
    
    print("\n[6/6] 生成可视化图表...")
    
    plot_scree_plot(pca.eigenvalues, pca.explained_variance_ratio, 
                   pca.cumulative_variance,
                   os.path.join(OUTPUT_DIR, "PCA_scree_plot.pdf"), n_pcs=20)
    
    plot_pc_scatter(pca.pc_scores, sample_ids, pc_x=1, pc_y=2,
                   var_explained_x=pca.explained_variance_ratio[0],
                   var_explained_y=pca.explained_variance_ratio[1],
                   output_file=os.path.join(OUTPUT_DIR, "PCA_PC1_vs_PC2.pdf"))
    
    plot_pc_scatter(pca.pc_scores, sample_ids, pc_x=1, pc_y=3,
                   var_explained_x=pca.explained_variance_ratio[0],
                   var_explained_y=pca.explained_variance_ratio[2],
                   output_file=os.path.join(OUTPUT_DIR, "PCA_PC1_vs_PC3.pdf"))
    
    plot_pc_scatter(pca.pc_scores, sample_ids, pc_x=2, pc_y=3,
                   var_explained_x=pca.explained_variance_ratio[1],
                   var_explained_y=pca.explained_variance_ratio[2],
                   output_file=os.path.join(OUTPUT_DIR, "PCA_PC2_vs_PC3.pdf"))
    
    plot_pc_density(pca.pc_scores, os.path.join(OUTPUT_DIR, "PCA_density.pdf"), n_pcs=4)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(f'PCA Comparison: Genotype vs GRM Method (Chr{TARGET_CHROMOSOME}, n={n_samples})', 
                 fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')
    
    ax1 = axes[0]
    ax1.scatter(pca.pc_scores[:, 0], pca.pc_scores[:, 1], 
               c=NATURE_COLORS[0], alpha=0.6, s=40, label='Genotype method')
    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio[0] * 100:.2f}%)', fontsize=14)
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio[1] * 100:.2f}%)', fontsize=14)
    ax1.set_title('PCA from Genotype', fontsize=16)
    ax1.grid(True, alpha=0.15, linestyle='--')
    ax1.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax1.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    sns.despine(ax=ax1)
    
    ax2 = axes[1]
    ax2.scatter(pca_grm.pc_scores[:, 0], pca_grm.pc_scores[:, 1], 
               c=NATURE_COLORS[1], alpha=0.6, s=40, label='GRM method')
    ax2.set_xlabel(f'PC1 ({pca_grm.explained_variance_ratio[0] * 100:.2f}%)', fontsize=14)
    ax2.set_ylabel(f'PC2 ({pca_grm.explained_variance_ratio[1] * 100:.2f}%)', fontsize=14)
    ax2.set_title('PCA from GRM', fontsize=16)
    ax2.grid(True, alpha=0.15, linestyle='--')
    ax2.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax2.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    sns.despine(ax=ax2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "PCA_method_comparison.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  方法比较图已保存: {os.path.join(OUTPUT_DIR, 'PCA_method_comparison.pdf')}")
    
    print("\n" + "=" * 70)
    print("分析完成!")
    print("=" * 70)
    print(f"\n输出文件列表:")
    print(f"  1. {output_pc_scores}")
    print(f"  2. {output_variance}")
    print(f"  3. {output_eigenvalues}")
    print(f"  4. {output_summary}")
    print(f"  5. {output_grm}")
    print(f"  6. {output_pc_scores_grm}")
    print(f"  7. {os.path.join(OUTPUT_DIR, 'PCA_scree_plot.pdf')}")
    print(f"  8. {os.path.join(OUTPUT_DIR, 'PCA_PC1_vs_PC2.pdf')}")
    print(f"  9. {os.path.join(OUTPUT_DIR, 'PCA_PC1_vs_PC3.pdf')}")
    print(f"  10. {os.path.join(OUTPUT_DIR, 'PCA_PC2_vs_PC3.pdf')}")
    print(f"  11. {os.path.join(OUTPUT_DIR, 'PCA_density.pdf')}")
    print(f"  12. {os.path.join(OUTPUT_DIR, 'PCA_method_comparison.pdf')}")


if __name__ == "__main__":
    main()
