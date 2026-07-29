#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算PLINK格式数据的等位基因频率和基因型频率并进行可视化
遵循科研绘图规范：Arial字体、PDF输出、Nature/Science配色
修改版：只抽取染色体1的SNP，使用前1000个样本
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

try:
    from pandas_plink import read_plink
except ImportError:
    print("正在安装 pandas-plink...")
    os.system("pip install pandas-plink")
    from pandas_plink import read_plink

INPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data"
OUTPUT_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t1"
BED_FILE = os.path.join(INPUT_DIR, "simulated_population")
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

def calculate_allele_frequency(genotypes):
    genotypes_valid = genotypes[~np.isnan(genotypes)]
    if len(genotypes_valid) == 0:
        return np.nan, np.nan, 0
    
    total_alleles = len(genotypes_valid) * 2
    minor_allele_count = np.sum(genotypes_valid)
    major_allele_count = total_alleles - minor_allele_count
    
    maf = minor_allele_count / total_alleles
    major_af = major_allele_count / total_alleles
    call_rate = len(genotypes_valid) / len(genotypes)
    
    return maf, major_af, call_rate

def calculate_genotype_frequency(genotypes):
    genotypes_valid = genotypes[~np.isnan(genotypes)]
    total = len(genotypes_valid)
    
    if total == 0:
        return {'AA': np.nan, 'AB': np.nan, 'BB': np.nan, 'missing': 1.0}
    
    aa_count = np.sum(genotypes_valid == 0)
    ab_count = np.sum(genotypes_valid == 1)
    bb_count = np.sum(genotypes_valid == 2)
    missing = np.sum(np.isnan(genotypes))
    
    return {
        'AA': aa_count / total,
        'AB': ab_count / total,
        'BB': bb_count / total,
        'missing': missing / len(genotypes)
    }

def sort_chromosomes(chrom_list):
    def chrom_sort_key(x):
        x_str = str(x)
        if x_str.isdigit():
            return (0, int(x_str))
        else:
            return (1, x_str)
    
    return sorted(chrom_list, key=chrom_sort_key)

def main():
    print("=" * 60)
    print("开始读取PLINK文件...")
    print(f"只分析染色体 {TARGET_CHROMOSOME} 的SNP")
    print(f"只使用前 {TARGET_SAMPLES} 个样本")
    print("=" * 60)
    
    bim, fam, bed = read_plink(BED_FILE, verbose=False)
    
    chrom_mask = bim['chrom'].astype(str) == TARGET_CHROMOSOME
    chr1_snp_indices = np.where(chrom_mask)[0]
    n_chr1_snps = len(chr1_snp_indices)
    
    n_total_samples = len(fam)
    n_samples = min(TARGET_SAMPLES, n_total_samples)
    
    print(f"染色体{TARGET_CHROMOSOME}的SNP数量: {n_chr1_snps}")
    print(f"总样本数量: {n_total_samples}")
    print(f"使用的样本数量: {n_samples}")
    
    print("\n正在计算等位基因频率和基因型频率...")
    
    results = []
    batch_size = 500
    total_batches = (n_chr1_snps + batch_size - 1) // batch_size
    
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_chr1_snps)
        
        if (batch_idx + 1) % 5 == 0 or batch_idx == total_batches - 1:
            print(f"处理进度: {end_idx}/{n_chr1_snps} ({100*end_idx/n_chr1_snps:.1f}%)")
        
        snp_indices = chr1_snp_indices[start_idx:end_idx]
        bed_batch = bed[snp_indices, :n_samples].compute()
        
        for i, snp_idx in enumerate(snp_indices):
            genotypes = bed_batch[i, :]
            
            maf, major_af, call_rate = calculate_allele_frequency(genotypes)
            gt_freq = calculate_genotype_frequency(genotypes)
            
            results.append({
                'CHR': bim['chrom'].values[snp_idx],
                'SNP': bim['snp'].values[snp_idx],
                'POS': bim['pos'].values[snp_idx],
                'A1': bim['a1'].values[snp_idx],
                'A2': bim['a0'].values[snp_idx],
                'MAF': maf,
                'Major_AF': major_af,
                'Call_Rate': call_rate,
                'AA_Freq': gt_freq['AA'],
                'AB_Freq': gt_freq['AB'],
                'BB_Freq': gt_freq['BB'],
                'Missing_Rate': gt_freq['missing']
            })
    
    df_results = pd.DataFrame(results)
    
    print("\n" + "=" * 60)
    print("保存结果...")
    print("=" * 60)
    
    output_csv = os.path.join(OUTPUT_DIR, "allele_genotype_frequency.csv")
    df_results.to_csv(output_csv, index=False)
    print(f"结果已保存至: {output_csv}")
    
    output_summary = os.path.join(OUTPUT_DIR, "frequency_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("等位基因频率和基因型频率分析摘要\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"分析染色体: {TARGET_CHROMOSOME}\n")
        f.write(f"总样本数: {n_samples}\n")
        f.write(f"总SNP数: {n_chr1_snps}\n\n")
        f.write("等位基因频率统计:\n")
        f.write(f"  MAF均值: {df_results['MAF'].mean():.4f}\n")
        f.write(f"  MAF中位数: {df_results['MAF'].median():.4f}\n")
        f.write(f"  MAF标准差: {df_results['MAF'].std():.4f}\n")
        f.write(f"  MAF范围: [{df_results['MAF'].min():.4f}, {df_results['MAF'].max():.4f}]\n\n")
        f.write("基因型频率统计:\n")
        f.write(f"  AA频率均值: {df_results['AA_Freq'].mean():.4f}\n")
        f.write(f"  AB频率均值: {df_results['AB_Freq'].mean():.4f}\n")
        f.write(f"  BB频率均值: {df_results['BB_Freq'].mean():.4f}\n\n")
        f.write("数据质量:\n")
        f.write(f"  平均缺失率: {df_results['Missing_Rate'].mean():.4f}\n")
        f.write(f"  平均检出率: {df_results['Call_Rate'].mean():.4f}\n")
        f.write(f"  MAF < 0.01 的SNP数: {len(df_results[df_results['MAF'] < 0.01])}\n")
        f.write(f"  MAF < 0.05 的SNP数: {len(df_results[df_results['MAF'] < 0.05])}\n")
    print(f"摘要已保存至: {output_summary}")
    
    print("\n" + "=" * 60)
    print("生成可视化图表...")
    print("=" * 60)
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f'Allele and Genotype Frequency Analysis (Chr{TARGET_CHROMOSOME}, n={n_samples})', 
                 fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')
    
    ax1 = axes[0, 0]
    maf_values = df_results['MAF'].dropna()
    ax1.hist(maf_values, bins=50, color=NATURE_COLORS[0], edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Minor Allele Frequency (MAF)', fontsize=14)
    ax1.set_ylabel('Number of SNPs', fontsize=14)
    ax1.set_title('MAF Distribution', fontsize=16)
    ax1.axvline(x=0.01, color=NATURE_COLORS[1], linestyle='--', linewidth=2, label='MAF=0.01')
    ax1.axvline(x=0.05, color=NATURE_COLORS[2], linestyle='--', linewidth=2, label='MAF=0.05')
    ax1.grid(True, alpha=0.15, linestyle='--')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    sns.despine(ax=ax1)
    
    ax2 = axes[0, 1]
    bp = ax2.boxplot([df_results['AA_Freq'].dropna(), 
                     df_results['AB_Freq'].dropna(), 
                     df_results['BB_Freq'].dropna()],
                    tick_labels=['AA (A1A1)', 'AB (A1A2)', 'BB (A2A2)'],
                    patch_artist=True,
                    boxprops=dict(facecolor=NATURE_COLORS[0], alpha=0.7),
                    medianprops=dict(color=NATURE_COLORS[1], linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5))
    ax2.set_ylabel('Frequency', fontsize=14)
    ax2.set_title('Genotype Frequency Distribution', fontsize=16)
    ax2.grid(True, alpha=0.15, linestyle='--', axis='y')
    
    ax3 = axes[0, 2]
    geno_means = [df_results['AA_Freq'].mean(), 
                  df_results['AB_Freq'].mean(), 
                  df_results['BB_Freq'].mean()]
    geno_labels = ['AA', 'AB', 'BB']
    bars = ax3.bar(geno_labels, geno_means, color=[NATURE_COLORS[0], NATURE_COLORS[1], NATURE_COLORS[2]], 
                   edgecolor='black', alpha=0.8)
    ax3.set_ylabel('Mean Frequency', fontsize=14)
    ax3.set_title('Average Genotype Frequency', fontsize=16)
    ax3.grid(True, alpha=0.15, linestyle='--', axis='y')
    for bar, val in zip(bars, geno_means):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{val:.3f}', ha='center', va='bottom', fontsize=12)
    sns.despine(ax=ax3)
    
    ax4 = axes[1, 0]
    call_rate_values = df_results['Call_Rate'].dropna()
    ax4.hist(call_rate_values, bins=50, color=NATURE_COLORS[3], edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Call Rate', fontsize=14)
    ax4.set_ylabel('Number of SNPs', fontsize=14)
    ax4.set_title('SNP Call Rate Distribution', fontsize=16)
    ax4.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax4)
    
    ax5 = axes[1, 1]
    ax5.scatter(df_results['POS'], df_results['MAF'], alpha=0.3, s=5, c=NATURE_COLORS[4])
    ax5.set_xlabel('Position', fontsize=14)
    ax5.set_ylabel('MAF', fontsize=14)
    ax5.set_title(f'MAF along Chromosome {TARGET_CHROMOSOME}', fontsize=16)
    ax5.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax5)
    
    ax6 = axes[1, 2]
    maf_categories = ['MAF < 0.01', '0.01-0.05', '0.05-0.1', 'MAF >= 0.1']
    maf_counts = [
        len(df_results[df_results['MAF'] < 0.01]),
        len(df_results[(df_results['MAF'] >= 0.01) & (df_results['MAF'] < 0.05)]),
        len(df_results[(df_results['MAF'] >= 0.05) & (df_results['MAF'] < 0.1)]),
        len(df_results[df_results['MAF'] >= 0.1])
    ]
    colors_pie = [NATURE_COLORS[1], NATURE_COLORS[0], NATURE_COLORS[3], NATURE_COLORS[4]]
    wedges, texts, autotexts = ax6.pie(maf_counts, autopct='%1.1f%%',
                                        colors=colors_pie, startangle=90, 
                                        pctdistance=0.75,
                                        textprops={'fontsize': 11})
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_fontweight('bold')
    
    ax6.legend(wedges, maf_categories, 
               loc='center left', 
               bbox_to_anchor=(1, 0.5), 
               frameon=False,
               fontsize=11)
    ax6.set_title('MAF Category Distribution', fontsize=16)
    
    plt.tight_layout()
    output_fig = os.path.join(OUTPUT_DIR, "frequency_visualization.pdf")
    plt.savefig(output_fig, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"可视化图表已保存至: {output_fig}")
    
    fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
    fig2.suptitle(f'Detailed Genotype Frequency Analysis (Chr{TARGET_CHROMOSOME})', 
                  fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')
    
    ax_geno1 = axes2[0]
    sample_size = min(1000, len(df_results))
    sample_df = df_results.sample(n=sample_size, random_state=42)
    x = range(sample_size)
    ax_geno1.bar(x, sample_df['AA_Freq'], label='AA', color=NATURE_COLORS[0], alpha=0.8)
    ax_geno1.bar(x, sample_df['AB_Freq'], bottom=sample_df['AA_Freq'], label='AB', color=NATURE_COLORS[1], alpha=0.8)
    ax_geno1.bar(x, sample_df['BB_Freq'], bottom=sample_df['AA_Freq']+sample_df['AB_Freq'], 
                 label='BB', color=NATURE_COLORS[2], alpha=0.8)
    ax_geno1.set_xlabel('SNP (Random Sample)', fontsize=14)
    ax_geno1.set_ylabel('Genotype Frequency', fontsize=14)
    ax_geno1.set_title(f'Genotype Frequency Stack Plot (n={sample_size})', fontsize=16)
    ax_geno1.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax_geno1.grid(True, alpha=0.15, linestyle='--', axis='y')
    sns.despine(ax=ax_geno1)
    
    ax_geno2 = axes2[1]
    ax_geno2.scatter(df_results['MAF'], df_results['AB_Freq'], alpha=0.4, s=2, c=NATURE_COLORS[0])
    ax_geno2.set_xlabel('MAF', fontsize=14)
    ax_geno2.set_ylabel('Heterozygous Frequency (AB)', fontsize=14)
    ax_geno2.set_title('MAF vs Heterozygous Frequency', fontsize=16)
    
    x_theory = np.linspace(0, 0.5, 100)
    y_theory = 2 * x_theory * (1 - x_theory)
    ax_geno2.plot(x_theory, y_theory, color=NATURE_COLORS[1], linestyle='-', linewidth=2.5, 
                   label='Hardy-Weinberg: 2p(1-p)')
    ax_geno2.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax_geno2.grid(True, alpha=0.15, linestyle='--')
    sns.despine(ax=ax_geno2)
    
    plt.tight_layout()
    output_fig2 = os.path.join(OUTPUT_DIR, "genotype_frequency_detail.pdf")
    plt.savefig(output_fig2, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"基因型频率详细图已保存至: {output_fig2}")
    
    print("\n" + "=" * 60)
    print("分析完成!")
    print("=" * 60)
    print(f"\n输出文件:")
    print(f"  1. {output_csv}")
    print(f"  2. {output_summary}")
    print(f"  3. {output_fig}")
    print(f"  4. {output_fig2}")

if __name__ == "__main__":
    main()
