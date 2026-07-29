#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PBLUP模型脚本 - 与HIBLUP输出格式一致
使用EMAI估计方差组分，PA系谱关系矩阵

用法:
  python PBLUP_model.py \
    --pheno /path/to/phenotypes.txt \
    --pedigree /path/to/pedigree.txt \
    --out /path/to/output_dir \
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

PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)


def parse_args():
    parser = argparse.ArgumentParser(
        description='PBLUP模型分析 (与HIBLUP格式一致) - 基于系谱关系矩阵',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--pheno', required=True, help='表型文件路径 (需包含ID列和表型列)')
    parser.add_argument('--pedigree', required=True, help='系谱文件路径 (需包含ID, Sire, Dam列)')
    parser.add_argument('--out', required=True, help='输出目录')
    parser.add_argument('--pheno-col', default='Phenotype', help='表型列名 (默认: Phenotype)')
    parser.add_argument('--id-col', default='ID', help='样本ID列名 (默认: ID)')
    parser.add_argument('--font-size', type=int, default=12, help='Font size for figures (default: 12)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for figure output (default: 300)')
    return parser.parse_args()


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
    if off_diag_max < 0.01:
        print("  警告: A矩阵接近单位矩阵，系谱信息量不足，遗传力估计可能不可靠")

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
    print("PBLUP模型分析 (与HIBLUP格式一致)")
    print("=" * 70)
    print(f"  表型文件: {args.pheno}")
    print(f"  系谱文件: {args.pedigree}")
    print(f"  输出目录: {args.out}")

    print("\n[1] 读取数据...")
    pheno_df = pd.read_csv(args.pheno, sep='\t')
    pheno_df[args.id_col] = pheno_df[args.id_col].astype(str)
    pheno_df = pheno_df.dropna(subset=[args.pheno_col])
    pheno_df = pheno_df.set_index(args.id_col).sort_index()
    y = pheno_df[args.pheno_col].values.astype(float)
    train_ids = pheno_df.index.tolist()
    print(f"  训练样本: {len(y)}")

    print("\n[2] 构建系谱关系矩阵...")
    A = build_pedigree_matrix(args.pedigree, train_ids, args.id_col)

    print("\n[3] EMAI-REML估计方差组分...")
    var_g, var_e, n_iter, converged = emai_reml(y, A)
    print(f"  迭代次数: {n_iter}, 收敛: {converged}")
    total_var = var_g + var_e
    h2 = var_g / total_var
    print(f"  V(PA)={var_g:.6f}, V(e)={var_e:.6f}, h2={h2:.6f}")

    print("\n[4] 求解混合模型方程...")
    n = len(y)
    mu = np.mean(y)
    lambda_ratio = var_g / var_e if var_e > 0 else 1.0
    C = A * lambda_ratio + np.eye(n)
    rhs = lambda_ratio * A @ (y - mu)
    try:
        pa_values = np.linalg.solve(C, rhs)
    except np.linalg.LinAlgError:
        pa_values = np.linalg.lstsq(C, rhs, rcond=None)[0]
    residuals = y - mu - pa_values

    print("\n[5] 保存结果 (HIBLUP格式)...")
    vars_file = os.path.join(args.out, "pblup_train.vars")
    with open(vars_file, 'w') as f:
        f.write("# EMAI-REML: iterations={}, converged={}\n".format(n_iter, converged))
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        f.write(f"PA\t{var_g:.6f}\tNA\t{h2:.6f}\tNA\tNA\n")
        f.write(f"e\t{var_e:.6f}\tNA\t{1-h2:.6f}\tNA\tNA\n")

    beta_file = os.path.join(args.out, "pblup_train.beta")
    with open(beta_file, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        f.write(f"mu\t{mu:.8f}\tNA\n")

    rand_file = os.path.join(args.out, "pblup_train.rand")
    with open(rand_file, 'w') as f:
        f.write("ID\tPA\tresiduals\n")
        for i, sid in enumerate(train_ids):
            f.write(f"{sid}\t{pa_values[i]:.6f}\t{residuals[i]:.6f}\n")

    # Bilingual summary files
    output_summary_en = os.path.join(args.out, "PBLUP_summary.txt")
    with open(output_summary_en, 'w') as f:
        f.write("PBLUP Model Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Method: BLUP (Pedigree-based)\n")
        f.write(f"Samples: {n}\n\n")
        f.write(f"Variance components:\n")
        f.write(f"  V(PA): {var_g:.6f}\n")
        f.write(f"  V(e):  {var_e:.6f}\n")
        f.write(f"  h2:    {h2:.6f}\n\n")
        f.write(f"Mean (mu): {mu:.8f}\n")
    print(f"  Summary (EN) saved: {output_summary_en}")

    output_summary_zh = os.path.join(args.out, "PBLUP_summary-zh.txt")
    with open(output_summary_zh, 'w') as f:
        f.write("PBLUP模型分析摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"方法: BLUP (基于系谱)\n")
        f.write(f"样本数: {n}\n\n")
        f.write(f"方差组分:\n")
        f.write(f"  V(PA): {var_g:.6f}\n")
        f.write(f"  V(e):  {var_e:.6f}\n")
        f.write(f"  h2:    {h2:.6f}\n\n")
        f.write(f"均值 (mu): {mu:.8f}\n")
    print(f"  摘要 (中文) 已保存: {output_summary_zh}")

    print("\n[6] 生成可视化图...")
    train_pred = mu + pa_values
    train_corr = np.corrcoef(y, train_pred)[0, 1] if np.std(train_pred) > 0 else 0

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    ax1.scatter(y, train_pred, alpha=0.5, c=PALETTE_A[0])
    ax1.plot([y.min(), y.max()], [y.min(), y.max()], color=PALETTE_A[4], linestyle='--')
    ax1.set_xlabel('Observed Phenotype')
    ax1.set_ylabel('Predicted Phenotype')
    ax1.set_title(f'Training Set (r={train_corr:.3f})')

    ax2 = axes[1]
    ax2.hist(pa_values, bins=30, edgecolor='#333333', alpha=0.7, color=PALETTE_A[0])
    ax2.axvline(x=np.mean(pa_values), color=PALETTE_A[4], linestyle='--', label=f'Mean={np.mean(pa_values):.2f}')
    ax2.set_xlabel('Estimated Breeding Value (PA)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Breeding Value Distribution')
    leg = ax2.legend()
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    leg.get_frame().set_edgecolor('none')

    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out, f"PBLUP_analysis.{fmt}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  图表已保存: PBLUP_analysis.pdf/png")

    print(f"  所有结果已保存到: {args.out}")
    print("\n" + "=" * 70)
    print("PBLUP模型分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
