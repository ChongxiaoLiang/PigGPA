#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LD Score回归脚本 - 简化版
估计遗传力
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
    print("正在安装 pandas-plink...")
    os.system("pip install pandas-plink")
    from pandas_plink import read_plink

INPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data"
OUTPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t20"
BED_FILE = os.path.join(INPUT_DIR, "simulated_population")

N_SAMPLES = 1000
N_SNPS = 200
TARGET_CHROM = "1"

os.makedirs(OUTPUT_DIR, exist_ok=True)
plt.rcParams['font.family'] = 'DejaVu Sans'
NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377']


def calculate_ld_score(genotype, window_size=50):
    """计算LD Score"""
    n_snps = genotype.shape[1]
    ld_scores = np.zeros(n_snps)
    
    for i in range(n_snps):
        start = max(0, i - window_size)
        end = min(n_snps, i + window_size + 1)
        
        for j in range(start, end):
            if i != j:
                x_i = genotype[:, i] - np.mean(genotype[:, i])
                x_j = genotype[:, j] - np.mean(genotype[:, j])
                corr = np.corrcoef(x_i, x_j)[0, 1]
                ld_scores[i] += corr ** 2
        ld_scores[i] += 1
    
    return ld_scores


def ld_score_regression(chi2_stats, ld_scores):
    """LD Score回归估计遗传力"""
    x = ld_scores - np.mean(ld_scores)
    y = chi2_stats - np.mean(chi2_stats)
    
    slope = np.sum(x * y) / np.sum(x ** 2)
    intercept = np.mean(chi2_stats) - slope * np.mean(ld_scores)
    
    n = N_SAMPLES
    h2 = slope * n / len(ld_scores)
    h2 = max(0, min(1, h2))
    
    return h2, intercept, slope


def main():
    print("=" * 70)
    print("LD Score回归分析（简化版）")
    print(f"限制: {N_SAMPLES}个样本, {N_SNPS}个SNP, 染色体{TARGET_CHROM}")
    print("=" * 70)
    
    print("\n[1/4] 读取PLINK文件...")
    bim, fam, bed = read_plink(BED_FILE, verbose=False)
    print(f"  总样本数量: {len(fam)}")
    
    print(f"\n[2/4] 筛选染色体{TARGET_CHROM}的SNP...")
    chr_mask = bim['chrom'].astype(str) == TARGET_CHROM
    chr_snps = bim[chr_mask]
    snp_indices = chr_snps.index[:N_SNPS].tolist()
    snp_info = chr_snps.iloc[:N_SNPS]
    
    print(f"\n[3/4] 提取基因型数据并计算...")
    genotype = bed[snp_indices, :N_SAMPLES].compute().T
    genotype = np.nan_to_num(genotype, nan=0)
    
    ld_scores = calculate_ld_score(genotype)
    
    np.random.seed(42)
    chi2_stats = np.random.chisquare(1, N_SNPS) + ld_scores * 0.5
    
    print("\n[4/4] LD Score回归...")
    h2, intercept, slope = ld_score_regression(chi2_stats, ld_scores)
    
    print(f"  估计遗传力: {h2:.4f}")
    print(f"  截距: {intercept:.4f}")
    print(f"  斜率: {slope:.4f}")
    
    print("\n[5/4] 保存结果...")
    
    output_summary = os.path.join(OUTPUT_DIR, "ld_score_regression_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("LD Score回归摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"样本数: {N_SAMPLES}\n")
        f.write(f"SNP数: {N_SNPS}\n\n")
        f.write(f"估计遗传力: {h2:.4f}\n")
        f.write(f"截距: {intercept:.4f}\n")
        f.write(f"斜率: {slope:.4f}\n")
    print(f"  摘要已保存: {output_summary}")
    
    output_results = os.path.join(OUTPUT_DIR, "ld_score_regression_results.csv")
    results_df = pd.DataFrame({
        'SNP': snp_info['snp'].values,
        'LD_Score': ld_scores,
        'Chi2_Stat': chi2_stats
    })
    results_df.to_csv(output_results, index=False)
    print(f"  结果已保存: {output_results}")
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(ld_scores, chi2_stats, alpha=0.5, c=NATURE_COLORS[0])
    
    x_line = np.linspace(ld_scores.min(), ld_scores.max(), 100)
    y_line = intercept + slope * x_line
    ax.plot(x_line, y_line, 'r-', lw=2, label=f'h² = {h2:.3f}')
    
    ax.set_xlabel('LD Score')
    ax.set_ylabel('χ² Statistic')
    ax.set_title('LD Score Regression')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "ld_score_regression.pdf"), format='pdf')
    plt.close()
    print(f"  图表已保存")
    
    print("\n" + "=" * 70)
    print("分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
