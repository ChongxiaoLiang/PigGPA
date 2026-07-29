#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GBLUP模型脚本 - 与HIBLUP输出格式一致
使用EMAI估计方差组分，GA关系矩阵

注意：GBLUP (Genomic BLUP) 与 BLUP 在当前实现中等价。
- GBLUP 使用基因组关系矩阵 G (VanRaden方法)，基于SNP标记计算个体间关系。
- 传统 BLUP 在动物育种中通常指使用系谱关系矩阵 A 的 PBLUP，
  但在本项目中 BLUP_model.py 同样使用 G 矩阵，因此两者实现一致。
- 两者的核心区别在于关系矩阵来源：GBLUP 基于基因组信息，PBLUP 基于系谱信息。

用法:
  python GBLUP_model.py \
    --bfile /path/to/plink_prefix \
    --pheno /path/to/phenotypes.txt \
    --out /path/to/output_dir \
    [--chrom 1,2,3] \
    [--pred-id-file /path/to/pred_ids.txt] \
    [--pheno-col Phenotype] \
    [--id-col ID]
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

try:
    from pandas_plink import read_plink
except ImportError:
    os.system("pip install pandas_plink")
    from pandas_plink import read_plink

PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)


def parse_args():
    parser = argparse.ArgumentParser(
        description='GBLUP模型分析 (与HIBLUP格式一致) - 使用EMAI估计方差组分',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--bfile', required=True, help='PLINK二进制文件前缀 (不含.bed/.bim/.fam)')
    parser.add_argument('--pheno', required=True, help='表型文件路径 (需包含ID列和表型列)')
    parser.add_argument('--out', required=True, help='输出目录')
    parser.add_argument('--chrom', default=None, help='指定染色体，逗号分隔 (默认: 全部染色体)')
    parser.add_argument('--pred-id-file', default=None, help='预测样本ID文件，每行一个ID (默认: 不进行预测)')
    parser.add_argument('--pheno-col', default='Phenotype', help='表型列名 (默认: Phenotype)')
    parser.add_argument('--id-col', default='ID', help='样本ID列名 (默认: ID)')
    parser.add_argument('--font-size', type=int, default=12, help='Font size for figures (default: 12)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for figure output (default: 300)')
    return parser.parse_args()


def calc_grm(genotype):
    p = np.nanmean(genotype, axis=0) / 2
    Z = genotype - 2 * p
    Z = np.nan_to_num(Z, nan=0.0)
    denom = 2 * p * (1 - p)
    denom[denom < 1e-10] = 1e-10
    G = (Z @ Z.T) / np.sum(denom)
    return G, p


def emai_reml(y, K, max_iter=50, tol=1e-4, em_rounds=5):
    n = len(y)
    y_c = y - np.mean(y)
    var_g = np.var(y_c) * 0.5
    var_e = np.var(y_c) * 0.5
    converged = False
    n_iter = 0
    for iteration in range(max_iter):
        n_iter = iteration + 1
        V = K * var_g + np.eye(n) * var_e
        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            V_inv = np.linalg.pinv(V)
        ones = np.ones((n, 1))
        Vinv_ones = V_inv @ ones
        XtVinvX = ones.T @ Vinv_ones
        if XtVinvX > 0:
            P = V_inv - Vinv_ones @ Vinv_ones.T / XtVinvX
        else:
            P = V_inv
        Py = P @ y_c
        PK = P @ K

        if iteration < em_rounds:
            u_hat = var_g * K @ Py
            var_g_new = (u_hat @ Py + var_g * (1 - var_g * np.trace(PK))) / n
            e_hat = y_c - u_hat
            var_e_new = (e_hat @ P @ y_c * var_e + var_e * (1 - var_e * np.trace(P))) / n
        else:
            score_g = -0.5 * np.trace(PK) + 0.5 * (y_c @ PK @ Py)
            score_e = -0.5 * np.trace(P) + 0.5 * (y_c @ P @ Py)
            ai_gg = 0.5 * np.trace(PK @ PK)
            ai_ge = 0.5 * np.trace(PK @ P)
            ai_ee = 0.5 * np.trace(P @ P)
            AI = np.array([[ai_gg, ai_ge], [ai_ge, ai_ee]])
            score = np.array([score_g, score_e])
            try:
                delta = np.linalg.solve(AI, score)
            except np.linalg.LinAlgError:
                delta = np.linalg.lstsq(AI, score, rcond=None)[0]
            step = 1.0
            var_g_new = var_g
            var_e_new = var_e
            for _ in range(10):
                candidate_g = var_g + step * delta[0]
                candidate_e = var_e + step * delta[1]
                if candidate_g > 0 and candidate_e > 0:
                    var_g_new = candidate_g
                    var_e_new = candidate_e
                    break
                step *= 0.5
            else:
                var_g_new = max(1e-10, var_g + 0.001 * score_g)
                var_e_new = max(1e-10, var_e + 0.001 * score_e)

        var_g_new = max(1e-10, var_g_new)
        var_e_new = max(1e-10, var_e_new)
        if abs(var_g_new - var_g) < tol and abs(var_e_new - var_e) < tol:
            converged = True
            var_g = var_g_new
            var_e = var_e_new
            break
        var_g = var_g_new
        var_e = var_e_new
    return var_g, var_e, n_iter, converged


def main():
    args = parse_args()
    plt.rcParams.update({
        'font.size': args.font_size,
        'font.family': 'DejaVu Sans',
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
    })
    os.makedirs(args.out, exist_ok=True)

    print("=" * 70)
    print("GBLUP模型分析 (与HIBLUP格式一致)")
    print("=" * 70)
    print(f"  PLINK文件: {args.bfile}")
    print(f"  表型文件: {args.pheno}")
    print(f"  输出目录: {args.out}")
    print(f"  染色体: {args.chrom if args.chrom else '全部'}")

    print("\n[1] 读取数据...")
    bim, fam, bed = read_plink(args.bfile, verbose=False)

    if args.chrom:
        chrom_list = [c.strip() for c in args.chrom.split(',')]
        snp_mask = bim['chrom'].astype(str).isin(chrom_list)
    else:
        snp_mask = np.ones(len(bim), dtype=bool)
    filtered_snps = bim[snp_mask]
    snp_indices = filtered_snps.index.tolist()
    snp_ids = filtered_snps['snp'].values
    snp_a1 = filtered_snps['a0'].values
    snp_a2 = filtered_snps['a1'].values
    print(f"  SNP数量: {len(snp_indices)}")

    pheno_df = pd.read_csv(args.pheno, sep='\t')
    pheno_df[args.id_col] = pheno_df[args.id_col].astype(str)
    pheno_df = pheno_df.dropna(subset=[args.pheno_col])
    pheno_df = pheno_df.set_index(args.id_col).sort_index()

    fam_ids = fam['iid'].astype(str).values
    train_id_to_fam_idx = {}
    for idx, fid in enumerate(fam_ids):
        if fid in pheno_df.index:
            train_id_to_fam_idx[fid] = idx
    matched_train_ids = sorted(train_id_to_fam_idx.keys(), key=lambda x: train_id_to_fam_idx[x])
    train_fam_indices = [train_id_to_fam_idx[tid] for tid in matched_train_ids]
    pheno_sub = pheno_df.loc[matched_train_ids]
    y = pheno_sub[args.pheno_col].values.astype(float)
    train_ids = matched_train_ids
    print(f"  训练样本: {len(y)}")

    genotype = bed[snp_indices, :][:, train_fam_indices].compute().T
    genotype = np.nan_to_num(genotype, nan=0)
    print(f"  训练基因型: {genotype.shape}")

    print("\n[2] 构建GRM...")
    G, p_train = calc_grm(genotype)
    print(f"  GRM维度: {G.shape}")

    print("\n[3] EMAI-REML估计方差组分...")
    var_g, var_e, n_iter, converged = emai_reml(y, G)
    total_var = var_g + var_e
    h2 = var_g / total_var
    print(f"  V(GA)={var_g:.6f}, V(e)={var_e:.6f}, h2={h2:.6f}")
    print(f"  迭代次数: {n_iter}, 收敛: {converged}")

    print("\n[4] 求解混合模型方程...")
    n = len(y)
    mu = np.mean(y)
    lambda_ratio = var_g / var_e if var_e > 0 else 1.0
    C = G * lambda_ratio + np.eye(n)
    rhs = lambda_ratio * G @ (y - mu)
    try:
        ga_values = np.linalg.solve(C, rhs)
    except np.linalg.LinAlgError:
        ga_values = np.linalg.lstsq(C, rhs, rcond=None)[0]
    residuals = y - mu - ga_values

    print("\n[5] 计算SNP效应值...")
    Z = genotype - 2 * p_train
    sum_2pq = 2 * np.sum(p_train * (1 - p_train))
    snp_effects = (Z.T @ ga_values) / sum_2pq
    freq_a1 = np.mean(genotype, axis=0) / 2

    pred_bv = None
    pred_ids = None
    if args.pred_id_file:
        print("\n[6] 预测新数据...")
        pred_id_list = []
        with open(args.pred_id_file, 'r') as f:
            for line in f:
                pid = line.strip()
                if pid:
                    pred_id_list.append(pid)

        pred_id_to_fam_idx = {}
        for idx, fid in enumerate(fam_ids):
            if fid in pred_id_list:
                pred_id_to_fam_idx[fid] = idx
        matched_pred_ids = sorted(pred_id_to_fam_idx.keys(), key=lambda x: pred_id_to_fam_idx[x])
        pred_fam_indices = [pred_id_to_fam_idx[pid] for pid in matched_pred_ids]

        genotype_new = bed[snp_indices, :][:, pred_fam_indices].compute().T
        genotype_new = np.nan_to_num(genotype_new, nan=0)
        Z_new = genotype_new - 2 * p_train
        pred_bv = Z_new @ snp_effects
        pred_ids = matched_pred_ids
        print(f"  预测样本: {len(pred_ids)}, 均值: {np.mean(pred_bv):.4f}")

    print("\n[7] 保存结果 (HIBLUP格式)...")
    vars_file = os.path.join(args.out, "gblup_train.vars")
    with open(vars_file, 'w') as f:
        f.write(f"# EMAI-REML: iterations={n_iter}, converged={converged}\n")
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        f.write(f"GA\t{var_g:.6f}\tNA\t{h2:.6f}\tNA\tNA\n")
        f.write(f"e\t{var_e:.6f}\tNA\t{1-h2:.6f}\tNA\tNA\n")

    beta_file = os.path.join(args.out, "gblup_train.beta")
    with open(beta_file, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        f.write(f"mu\t{mu:.8f}\tNA\n")

    rand_file = os.path.join(args.out, "gblup_train.rand")
    with open(rand_file, 'w') as f:
        f.write("ID\tGA\tresiduals\n")
        for i, sid in enumerate(train_ids):
            f.write(f"{sid}\t{ga_values[i]:.6f}\t{residuals[i]:.6f}\n")

    snpeff_file = os.path.join(args.out, "gblup_snp_effect.snpeff")
    with open(snpeff_file, 'w') as f:
        f.write("id\ta1\ta2\tfreq_a1\tadd_a1\n")
        for i in range(len(snp_indices)):
            f.write(f"{snp_ids[i]}\t{snp_a1[i]}\t{snp_a2[i]}\t{freq_a1[i]:.4f}\t{snp_effects[i]:.8f}\n")

    if pred_bv is not None:
        pred_file = os.path.join(args.out, "gblup_pred.bv")
        with open(pred_file, 'w') as f:
            f.write("id\tadd_a1\n")
            for i, pid in enumerate(pred_ids):
                f.write(f"{pid}\t{pred_bv[i]:.8f}\n")

    print("\n[8] 生成可视化图...")
    train_pred = mu + ga_values
    train_corr = np.corrcoef(y, train_pred)[0, 1] if np.std(train_pred) > 0 else 0

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    ax1.scatter(y, train_pred, alpha=0.5, c=PALETTE_A[0])
    ax1.plot([y.min(), y.max()], [y.min(), y.max()], color=PALETTE_A[4], linestyle='--')
    ax1.set_xlabel('Observed Phenotype')
    ax1.set_ylabel('Predicted Phenotype')
    ax1.set_title(f'Training Set (r={train_corr:.3f})')

    ax2 = axes[1]
    if pred_bv is not None:
        ax2.hist(pred_bv, bins=30, edgecolor='#333333', alpha=0.7, color=PALETTE_A[0])
        ax2.axvline(x=np.mean(pred_bv), color=PALETTE_A[4], linestyle='--', label=f'Mean={np.mean(pred_bv):.2f}')
        ax2.set_xlabel('Predicted Breeding Value')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Predicted BV Distribution')
        leg = ax2.legend()
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')
    else:
        ax2.hist(ga_values, bins=30, edgecolor='#333333', alpha=0.7, color=PALETTE_A[0])
        ax2.axvline(x=np.mean(ga_values), color=PALETTE_A[4], linestyle='--', label=f'Mean={np.mean(ga_values):.2f}')
        ax2.set_xlabel('Estimated Breeding Value (GA)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Breeding Value Distribution')
        leg = ax2.legend()
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')

    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out, f"GBLUP_analysis.{fmt}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  图表已保存: GBLUP_analysis.pdf/png")

    print(f"  所有结果已保存到: {args.out}")
    print("\n" + "=" * 70)
    print("GBLUP模型分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
