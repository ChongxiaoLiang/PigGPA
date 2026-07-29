#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SSBLUP模型脚本 - 与HIBLUP输出格式一致
使用EMAI估计方差组分，HA混合关系矩阵

用法:
  python SSBLUP_model.py \
    --bfile /path/to/plink_prefix \
    --pheno /path/to/phenotypes.txt \
    --pedigree /path/to/pedigree.txt \
    --out /path/to/output_dir \
    [--chrom 1,2,3] \
    [--pred-id-file /path/to/pred_ids.txt] \
    [--pheno-col Phenotype] \
    [--id-col ID] \
    [--g-weight 0.95] \
    [--a-weight 0.05]
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
        description='SSBLUP模型分析 (与HIBLUP格式一致) - 使用EMAI估计方差组分，HA混合关系矩阵',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--bfile', required=True, help='PLINK二进制文件前缀 (不含.bed/.bim/.fam)')
    parser.add_argument('--pheno', required=True, help='表型文件路径 (需包含ID列和表型列)')
    parser.add_argument('--pedigree', required=True, help='系谱文件路径 (需包含ID, Sire, Dam列)')
    parser.add_argument('--out', required=True, help='输出目录')
    parser.add_argument('--chrom', default=None, help='指定染色体，逗号分隔 (默认: 全部染色体)')
    parser.add_argument('--pred-id-file', default=None, help='预测样本ID文件，每行一个ID (默认: 不进行预测)')
    parser.add_argument('--pheno-col', default='Phenotype', help='表型列名 (默认: Phenotype)')
    parser.add_argument('--id-col', default='ID', help='样本ID列名 (默认: ID)')
    parser.add_argument('--g-weight', type=float, default=0.95, help='G矩阵权重 (默认: 0.95)')
    parser.add_argument('--a-weight', type=float, default=0.05, help='A矩阵权重 (默认: 0.05)')
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


def _topological_sort(ped_df, all_ids):
    id_set = set(all_ids)
    visited = set()
    result = []
    def visit(ind):
        if ind in visited:
            return
        visited.add(ind)
        row = ped_df[ped_df['id'] == ind]
        if len(row) > 0:
            sire = str(row.iloc[0].get('sire', '0'))
            dam = str(row.iloc[0].get('dam', '0'))
            if sire != '0' and sire in id_set:
                visit(sire)
            if dam != '0' and dam in id_set:
                visit(dam)
        result.append(ind)
    for ind in all_ids:
        visit(ind)
    return result


def build_pedigree_matrix(pedigree_file, sample_ids, id_col='ID'):
    ped_df = pd.read_csv(pedigree_file, sep='\t')
    ped_df[id_col] = ped_df[id_col].astype(str)
    ped_df['Sire'] = ped_df['Sire'].astype(str)
    ped_df['Dam'] = ped_df['Dam'].astype(str)
    ped_df = ped_df.rename(columns={id_col: 'id', 'Sire': 'sire', 'Dam': 'dam'})

    all_ids_set = set(sample_ids)
    for _, row in ped_df.iterrows():
        all_ids_set.add(str(row['id']))

    all_ids = sorted(all_ids_set)
    n_all = len(all_ids)
    id_to_idx = {ind: i for i, ind in enumerate(all_ids)}

    A_full = np.eye(n_all)
    sorted_all = _topological_sort(ped_df, all_ids)

    for ind in sorted_all:
        i = id_to_idx[ind]
        row = ped_df[ped_df['id'] == ind]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        sire = str(row.get('sire', '0'))
        dam = str(row.get('dam', '0'))
        sire_idx = id_to_idx.get(sire, None)
        dam_idx = id_to_idx.get(dam, None)
        if sire_idx is not None and dam_idx is not None:
            A_full[i, i] = 1 + 0.5 * A_full[sire_idx, dam_idx]
        for j_id in sorted_all:
            j = id_to_idx[j_id]
            if j >= i:
                break
            val = 0.0
            if sire_idx is not None:
                val += A_full[sire_idx, j]
            if dam_idx is not None:
                val += A_full[dam_idx, j]
            A_full[i, j] = 0.5 * val
            A_full[j, i] = A_full[i, j]

    train_indices = [id_to_idx[sid] for sid in sample_ids if sid in id_to_idx]
    A = A_full[np.ix_(train_indices, train_indices)]

    off_diag = A.copy()
    np.fill_diagonal(off_diag, 0)
    diag_mean = np.mean(np.diag(A))
    off_diag_max = np.max(np.abs(off_diag))
    print(f"  A矩阵: 维度={A.shape}, 对角线均值={diag_mean:.4f}, 非对角线最大绝对值={off_diag_max:.4f}")

    return A


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
    print("SSBLUP模型分析 (与HIBLUP格式一致)")
    print("=" * 70)
    print(f"  PLINK文件: {args.bfile}")
    print(f"  表型文件: {args.pheno}")
    print(f"  系谱文件: {args.pedigree}")
    print(f"  输出目录: {args.out}")
    print(f"  染色体: {args.chrom if args.chrom else '全部'}")
    print(f"  H矩阵权重: G={args.g_weight}, A={args.a_weight}")

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

    print("\n[2] 构建GRM和系谱矩阵...")
    G, p_train = calc_grm(genotype)
    A = build_pedigree_matrix(args.pedigree, train_ids, args.id_col)

    print("\n[3] 构建HA混合关系矩阵...")
    g_scale = np.mean(np.diag(A)) / np.mean(np.diag(G)) if np.mean(np.diag(G)) > 0 else 1.0
    G_scaled = G * g_scale
    print(f"  G矩阵缩放因子: {g_scale:.4f}")
    H = args.g_weight * G_scaled + args.a_weight * A
    print(f"  HA矩阵维度: {H.shape}, 对角线均值: {np.mean(np.diag(H)):.4f}")

    print("\n[4] EMAI-REML估计方差组分...")
    var_g, var_e, n_iter, converged = emai_reml(y, H)
    total_var = var_g + var_e
    h2 = var_g / total_var
    print(f"  V(HA)={var_g:.6f}, V(e)={var_e:.6f}, h2={h2:.6f}")
    print(f"  迭代次数: {n_iter}, 收敛: {converged}")

    print("\n[5] 求解混合模型方程...")
    n = len(y)
    mu = np.mean(y)
    lambda_ratio = var_g / var_e if var_e > 0 else 1.0
    C = H * lambda_ratio + np.eye(n)
    rhs = lambda_ratio * H @ (y - mu)
    try:
        ha_values = np.linalg.solve(C, rhs)
    except np.linalg.LinAlgError:
        ha_values = np.linalg.lstsq(C, rhs, rcond=None)[0]
    residuals = y - mu - ha_values

    print("\n[6] 计算SNP效应值...")
    Z = genotype - 2 * p_train
    sum_2pq = 2 * np.sum(p_train * (1 - p_train))
    snp_effects = (Z.T @ ha_values) / sum_2pq
    freq_a1 = np.mean(genotype, axis=0) / 2

    pred_bv = None
    pred_ids = None
    if args.pred_id_file:
        print("\n[7] 预测新数据...")
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

    print("\n[8] 保存结果 (HIBLUP格式)...")
    vars_file = os.path.join(args.out, "ssblup_train.vars")
    with open(vars_file, 'w') as f:
        f.write("# EMAI-REML: iterations={}, converged={}\n".format(n_iter, converged))
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        f.write(f"HA\t{var_g:.6f}\tNA\t{h2:.6f}\tNA\tNA\n")
        f.write(f"e\t{var_e:.6f}\tNA\t{1-h2:.6f}\tNA\tNA\n")

    beta_file = os.path.join(args.out, "ssblup_train.beta")
    with open(beta_file, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        f.write(f"mu\t{mu:.8f}\tNA\n")

    rand_file = os.path.join(args.out, "ssblup_train.rand")
    with open(rand_file, 'w') as f:
        f.write("ID\tHA\tresiduals\n")
        for i, sid in enumerate(train_ids):
            f.write(f"{sid}\t{ha_values[i]:.6f}\t{residuals[i]:.6f}\n")

    snpeff_file = os.path.join(args.out, "ssblup_snp_effect.snpeff")
    with open(snpeff_file, 'w') as f:
        f.write("id\ta1\ta2\tfreq_a1\tadd_a1\n")
        for i in range(len(snp_indices)):
            f.write(f"{snp_ids[i]}\t{snp_a1[i]}\t{snp_a2[i]}\t{freq_a1[i]:.4f}\t{snp_effects[i]:.8f}\n")

    if pred_bv is not None:
        pred_file = os.path.join(args.out, "ssblup_pred.bv")
        with open(pred_file, 'w') as f:
            f.write("id\tadd_a1\n")
            for i, pid in enumerate(pred_ids):
                f.write(f"{pid}\t{pred_bv[i]:.8f}\n")

    print("\n[9] 生成可视化图...")
    train_pred = mu + ha_values
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
        ax2.hist(ha_values, bins=30, edgecolor='#333333', alpha=0.7, color=PALETTE_A[0])
        ax2.axvline(x=np.mean(ha_values), color=PALETTE_A[4], linestyle='--', label=f'Mean={np.mean(ha_values):.2f}')
        ax2.set_xlabel('Estimated Breeding Value (HA)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Breeding Value Distribution')
        leg = ax2.legend()
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')

    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out, f"SSBLUP_analysis.{fmt}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  图表已保存: SSBLUP_analysis.pdf/png")

    print(f"  所有结果已保存到: {args.out}")
    print("\n" + "=" * 70)
    print("SSBLUP模型分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
