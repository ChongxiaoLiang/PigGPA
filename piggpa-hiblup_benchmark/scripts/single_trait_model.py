#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
单性状模型脚本 - 复刻HIBLUP --single-trait
输出格式与HIBLUP完全一致，便于对比验证

模型: Phenotype = mu + HA + e  (SSGBLUP)
  - HA: 混合关系矩阵效应 (H = w*G + (1-w)*A)
  - 使用AI-REML估计方差组分
  - 使用真实表型数据

用法:
  python single_trait_model.py \
    --bfile /path/to/simulated_population \
    --pheno /path/to/phenotypes.txt \
    --pedigree /path/to/pedigree.txt \
    --out /path/to/output_dir \
    [--g-weight 0.95]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 14

PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Single Trait Model (Replicating HIBLUP --single-trait)',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--bfile', required=True,
                        help='PLINK bed/bim/fam prefix path')
    parser.add_argument('--pheno', required=True,
                        help='Phenotype file path (tab-separated: ID, Phenotype, BreedingValue)')
    parser.add_argument('--pedigree', required=True,
                        help='Pedigree file path (tab-separated: ID, Sire, Dam, ...)')
    parser.add_argument('--out', required=True,
                        help='Output directory')
    parser.add_argument('--g-weight', type=float, default=0.95,
                        help='Weight for G matrix in HA construction (default: 0.95)')
    parser.add_argument('--em-max-iter', type=int, default=200,
                        help='Max iterations for EM-REML (default: 200)')
    parser.add_argument('--font-size', type=int, default=12, help='Font size for figures (default: 12)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for figure output (default: 300)')
    return parser.parse_args()


def read_plink_all(bfile):
    print(f"Reading PLINK files (all chromosomes, all samples)...")
    sys.stdout.flush()

    bim = pd.read_csv(bfile + ".bim", sep='\t', header=None,
                      names=['chrom', 'snp', 'cm', 'pos', 'a1', 'a2'])
    fam = pd.read_csv(bfile + ".fam", sep=r'\s+', header=None,
                      names=['fid', 'iid', 'pid', 'mid', 'sex', 'pheno'])

    n_snps = len(bim)
    snp_ids = bim['snp'].values
    sample_ids = fam['iid'].values
    n_samples = len(fam)

    print(f"  Total samples: {n_samples}, Total SNPs: {n_snps}")

    try:
        from pandas_plink import read_plink
        _, _, bed = read_plink(bfile, verbose=False)
        M = bed[:, :].compute().T.astype(float)
    except ImportError:
        import struct
        with open(bfile + ".bed", 'rb') as f:
            magic = f.read(3)
            n_total_samples = len(fam)
            n_total_snps = len(bim)

            M = np.zeros((n_samples, n_snps))

            for j in range(n_total_snps):
                f.seek(3 + j * ((n_total_samples + 3) // 4))
                bytes_per_snp = (n_total_samples + 3) // 4
                raw = f.read(bytes_per_snp)

                for i in range(n_samples):
                    byte_idx = i // 4
                    bit_idx = (i % 4) * 2
                    genotype = (raw[byte_idx] >> bit_idx) & 3

                    if genotype == 0:
                        M[i, j] = 2
                    elif genotype == 1:
                        M[i, j] = np.nan
                    elif genotype == 2:
                        M[i, j] = 1
                    elif genotype == 3:
                        M[i, j] = 0
    except Exception as e:
        print(f"  Warning: Error reading PLINK bed file: {e}")
        import struct
        with open(bfile + ".bed", 'rb') as f:
            magic = f.read(3)
            n_total_samples = len(fam)
            n_total_snps = len(bim)

            M = np.zeros((n_samples, n_snps))

            for j in range(n_total_snps):
                f.seek(3 + j * ((n_total_samples + 3) // 4))
                bytes_per_snp = (n_total_samples + 3) // 4
                raw = f.read(bytes_per_snp)

                for i in range(n_samples):
                    byte_idx = i // 4
                    bit_idx = (i % 4) * 2
                    genotype = (raw[byte_idx] >> bit_idx) & 3

                    if genotype == 0:
                        M[i, j] = 2
                    elif genotype == 1:
                        M[i, j] = np.nan
                    elif genotype == 2:
                        M[i, j] = 1
                    elif genotype == 3:
                        M[i, j] = 0

    p = np.nanmean(M, axis=0) / 2
    M = np.where(np.isnan(M), 2 * p, M)

    print(f"  Genotype matrix shape: {M.shape}")
    sys.stdout.flush()
    return M, sample_ids, snp_ids


def compute_grm_vanraden(M):
    print("Computing GRM (VanRaden method)...")
    sys.stdout.flush()

    n, m = M.shape

    p = np.nansum(M, axis=0) / (2 * n)
    p = np.clip(p, 1e-10, 1 - 1e-10)

    Z = M - 2 * p

    denom = 2 * np.sum(p * (1 - p))
    G = (Z @ Z.T) / denom

    diag_mean = np.mean(np.diag(G))
    G = G / diag_mean

    print(f"  GRM diagonal mean: {np.mean(np.diag(G)):.6f}")
    print(f"  GRM off-diagonal mean: {np.mean(G - np.diag(np.diag(G))):.6f}")
    sys.stdout.flush()
    return G


def compute_pedigree_A(pedigree_file, sample_ids):
    print("Computing pedigree relationship matrix A...")
    sys.stdout.flush()

    ped_df = pd.read_csv(pedigree_file, sep='\t')
    if 'ID' in ped_df.columns:
        ped_df = ped_df.rename(columns={'ID': 'Individual', 'Sire': 'Sire', 'Dam': 'Dam'})

    all_ids = ped_df['Individual'].values
    id_to_idx = {int(id): i for i, id in enumerate(all_ids)}
    n = len(all_ids)

    def _topological_sort(ped_df, all_ids):
        id_set = set(all_ids)
        visited = set()
        result = []
        def visit(ind):
            if ind in visited:
                return
            visited.add(ind)
            row = ped_df[ped_df['Individual'] == ind]
            if len(row) > 0:
                sire = int(row.iloc[0]['Sire']) if pd.notna(row.iloc[0]['Sire']) else 0
                dam = int(row.iloc[0]['Dam']) if pd.notna(row.iloc[0]['Dam']) else 0
                if sire != 0 and sire in id_set:
                    visit(sire)
                if dam != 0 and dam in id_set:
                    visit(dam)
            result.append(ind)
        for ind in all_ids:
            visit(ind)
        return result

    sorted_ids = _topological_sort(ped_df, all_ids)

    A = np.zeros((n, n))

    for ind in sorted_ids:
        i = id_to_idx[ind]
        row = ped_df[ped_df['Individual'] == ind]
        sire = int(row.iloc[0]['Sire']) if len(row) > 0 and pd.notna(row.iloc[0]['Sire']) else 0
        dam = int(row.iloc[0]['Dam']) if len(row) > 0 and pd.notna(row.iloc[0]['Dam']) else 0

        sire_idx = id_to_idx.get(sire) if sire != 0 else None
        dam_idx = id_to_idx.get(dam) if dam != 0 else None

        if sire_idx is not None and dam_idx is not None:
            A[i, i] = 1 + 0.5 * A[sire_idx, dam_idx]
        else:
            A[i, i] = 1.0

        for j_id in sorted_ids:
            j = id_to_idx[j_id]
            if j == i:
                break
            val = 0.0
            if sire_idx is not None:
                val += A[sire_idx, j]
            if dam_idx is not None:
                val += A[dam_idx, j]
            A[i, j] = 0.5 * val
            A[j, i] = A[i, j]

    sample_ids_int = [int(s) for s in sample_ids]
    sample_indices = [id_to_idx[s] for s in sample_ids_int if s in id_to_idx]

    A_sub = A[np.ix_(sample_indices, sample_indices)]

    print(f"  A matrix diagonal mean: {np.mean(np.diag(A_sub)):.6f}")
    print(f"  A matrix off-diagonal mean: {np.mean(A_sub - np.diag(np.diag(A_sub))):.6f}")
    sys.stdout.flush()
    return A_sub


def construct_HA_matrix(G, A, w=0.95):
    print(f"Constructing HA matrix (w={w})...")
    sys.stdout.flush()

    diag_G = np.mean(np.diag(G))
    offdiag_G = np.mean(G - np.diag(np.diag(G)))
    diag_A = np.mean(np.diag(A))
    offdiag_A = np.mean(A - np.diag(np.diag(A)))

    print(f"  Mean diagonal and Off-diagonal of A11: {diag_A:.0f} {offdiag_A:.6f}")
    print(f"  Mean diagonal and Off-diagonal of G11: {diag_G:.6f} {offdiag_G:.8f}")

    alpha = 0.999
    n = G.shape[0]
    G_adj = alpha * G + (1 - alpha) * np.eye(n)
    print(f"  Adjust G11: G11* = {alpha} * G11 + {1 - alpha} * I")

    HA = w * G_adj + (1 - w) * A
    print(f"  Weight of A11 and G11: Gw = {w} * G11 + {1 - w} * A11")
    print(f"  HA matrix construction accomplished.")
    sys.stdout.flush()
    return HA


def ai_reml(y, K, max_iter=20, tol=1e-8):
    print("Variance components estimation using: AI(20)")
    sys.stdout.flush()

    n = len(y)
    X = np.ones((n, 1))

    var_k = np.var(y) * 0.5
    var_e = np.var(y) * 0.5

    prev_loglik = -np.inf
    convergence_history = []

    print(f"The matrix V has a dimension of {n} x {n}.")
    print("Running ...")
    print(f"{'Alg.':<8} {'Iter.':<8} {'LogL.':<12} {'V(HA)':<12} {'V(e)':<12}")
    sys.stdout.flush()

    for iteration in range(max_iter):
        V = var_k * K + var_e * np.eye(n)

        try:
            V_inv = np.linalg.inv(V)
        except:
            V_inv = np.linalg.pinv(V)

        V_inv_X = V_inv @ X
        XtV_inv_X = X.T @ V_inv_X

        try:
            XtV_inv_X_inv = np.linalg.inv(XtV_inv_X)
        except:
            XtV_inv_X_inv = np.linalg.pinv(XtV_inv_X)

        P = V_inv - V_inv_X @ XtV_inv_X_inv @ V_inv_X.T
        Py = P @ y

        sign_v, logdet_v = np.linalg.slogdet(V)
        log_det_V = logdet_v if sign_v > 0 else np.log(1e-300)
        sign_x, logdet_x = np.linalg.slogdet(XtV_inv_X)
        log_det_X = logdet_x if sign_x > 0 else np.log(1e-300)
        loglik = -0.5 * (log_det_V + log_det_X + y.T @ Py)

        print(f"[AI]     {iteration+1:<8} {loglik:<12.2f} {var_k:<12.5f} {var_e:<12.5f}")
        sys.stdout.flush()

        convergence_history.append({
            'iteration': iteration + 1,
            'log_likelihood': loglik,
            'var_k': var_k,
            'var_e': var_e
        })

        if abs(loglik - prev_loglik) < tol and iteration > 0:
            print("[Converged?] Yes!")
            sys.stdout.flush()
            break
        prev_loglik = loglik

        PK = P @ K
        score_k = -0.5 * (np.trace(PK) - y.T @ PK @ Py)
        score_e = -0.5 * (np.trace(P) - y.T @ P @ Py)

        info_kk = 0.5 * np.trace(PK @ PK)
        info_ee = 0.5 * np.trace(P @ P)
        info_ke = 0.5 * np.trace(PK @ P)

        info_matrix = np.array([[info_kk, info_ke], [info_ke, info_ee]])
        score_vec = np.array([score_k, score_e])

        try:
            info_inv = np.linalg.inv(info_matrix)
        except:
            info_inv = np.linalg.pinv(info_matrix)

        updates = info_inv @ score_vec

        step = 1.0
        var_k_new = var_k + step * updates[0]
        var_e_new = var_e + step * updates[1]

        if var_k_new < 0:
            var_k_new = 1e-10
        if var_e_new < 0:
            var_e_new = 1e-10

        var_k = var_k_new
        var_e = var_e_new

    total_var = var_k + var_e
    h2 = var_k / total_var if total_var > 0 else 0

    h2_se = 0.0
    try:
        V = var_k * K + var_e * np.eye(n)
        V_inv = np.linalg.inv(V)
        V_inv_X = V_inv @ X
        XtV_inv_X = X.T @ V_inv_X
        XtV_inv_X_inv = np.linalg.inv(XtV_inv_X)
        P = V_inv - V_inv_X @ XtV_inv_X_inv @ V_inv_X.T

        PK = P @ K
        info_kk = 0.5 * np.trace(PK @ PK)
        info_ee = 0.5 * np.trace(P @ P)
        info_ke = 0.5 * np.trace(PK @ P)

        info_matrix = np.array([[info_kk, info_ke], [info_ke, info_ee]])
        info_inv = np.linalg.inv(info_matrix)

        var_k_se = np.sqrt(max(0, info_inv[0, 0]))
        var_e_se = np.sqrt(max(0, info_inv[1, 1]))

        dh2_dvk = var_e / (total_var ** 2)
        dh2_dve = -var_k / (total_var ** 2)
        grad = np.array([dh2_dvk, dh2_dve])
        h2_se = np.sqrt(max(0, grad @ info_inv @ grad))
    except:
        var_k_se = 0
        var_e_se = 0
        h2_se = 0

    from scipy.stats import chi2
    h2_pval = 1 - chi2.cdf(h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0
    e_h2 = var_e / total_var if total_var > 0 else 0
    e_h2_pval = 1 - chi2.cdf(e_h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0

    beta = XtV_inv_X_inv @ (X.T @ V_inv @ y)
    beta_se = np.sqrt(np.diag(XtV_inv_X_inv))

    V = var_k * K + var_e * np.eye(n)
    V_inv = np.linalg.inv(V)
    V_inv_X = V_inv @ X
    XtV_inv_X = X.T @ V_inv_X
    XtV_inv_X_inv = np.linalg.inv(XtV_inv_X)
    P = V_inv - V_inv_X @ XtV_inv_X_inv @ V_inv_X.T

    ha_effects = var_k * K @ V_inv @ (y - X @ beta)
    residuals = y - X @ beta - ha_effects

    return {
        'var_k': var_k, 'var_e': var_e,
        'var_k_se': var_k_se, 'var_e_se': var_e_se,
        'h2': h2, 'h2_se': h2_se,
        'h2_pval': h2_pval, 'e_h2_pval': e_h2_pval,
        'beta': beta, 'beta_se': beta_se,
        'ha_effects': ha_effects, 'residuals': residuals,
        'loglik': loglik, 'converged': True, 'n_iter': iteration + 1,
        'convergence_history': convergence_history
    }


def em_reml(y, K, max_iter=200, tol=1e-6):
    """EM-REML方差组分估计

    参数:
        y: 表型向量
        K: 亲缘关系矩阵 (HA)
        max_iter: 最大迭代次数 (默认200)
        tol: 收敛阈值
    """
    print(f"Variance components estimation using: EM({max_iter})")
    sys.stdout.flush()

    n = len(y)
    X = np.ones((n, 1))

    var_k = np.var(y) * 0.5
    var_e = np.var(y) * 0.5

    prev_loglik = -np.inf
    convergence_history = []
    n_iter = 0  # 初始化迭代计数器，确保非收敛时也有正确值
    loglik = -np.inf
    converged = False

    print(f"The matrix V has a dimension of {n} x {n}.")
    print("Running ...")
    print(f"{'Alg.':<8} {'Iter.':<8} {'LogL.':<12} {'V(HA)':<12} {'V(e)':<12}")
    sys.stdout.flush()

    for iteration in range(max_iter):
        n_iter = iteration + 1  # 每次迭代更新计数器

        V = var_k * K + var_e * np.eye(n)

        try:
            V_inv = np.linalg.inv(V)
        except:
            V_inv = np.linalg.pinv(V)

        V_inv_X = V_inv @ X
        XtV_inv_X = X.T @ V_inv_X

        try:
            XtV_inv_X_inv = np.linalg.inv(XtV_inv_X)
        except:
            XtV_inv_X_inv = np.linalg.pinv(XtV_inv_X)

        P = V_inv - V_inv_X @ XtV_inv_X_inv @ V_inv_X.T
        Py = P @ y

        sign_v, logdet_v = np.linalg.slogdet(V)
        log_det_V = logdet_v if sign_v > 0 else np.log(1e-300)
        sign_x, logdet_x = np.linalg.slogdet(XtV_inv_X)
        log_det_X = logdet_x if sign_x > 0 else np.log(1e-300)
        loglik = -0.5 * (log_det_V + log_det_X + y.T @ Py)

        # 每次迭代的logL写入日志
        print(f"[EM]     {iteration+1:<8} {loglik:<12.2f} {var_k:<12.5f} {var_e:<12.5f}")
        sys.stdout.flush()

        convergence_history.append({
            'iteration': iteration + 1,
            'log_likelihood': loglik,
            'var_k': var_k,
            'var_e': var_e
        })

        if abs(loglik - prev_loglik) < tol and iteration > 0:
            print("[Converged?] Yes!")
            sys.stdout.flush()
            converged = True
            break
        prev_loglik = loglik

        # EM更新步骤
        PK = P @ K
        var_k_new = (var_k**2 * y.T @ PK @ Py +
                     var_k * np.trace(K @ V_inv) -
                     var_k**2 * np.trace(PK)) / n
        var_e_new = (var_e**2 * y.T @ P @ Py +
                     var_e * np.trace(V_inv) -
                     var_e**2 * np.trace(P)) / n

        var_k = max(1e-10, float(var_k_new))
        var_e = max(1e-10, float(var_e_new))

    if not converged:
        print(f"[Converged?] No (reached max_iter={max_iter})")
        sys.stdout.flush()

    # 确保n_iter被正确返回（无论是否收敛）
    total_var = var_k + var_e
    h2 = var_k / total_var if total_var > 0 else 0

    # 计算标准误（使用最后一次迭代的信息矩阵）
    h2_se = 0.0
    var_k_se = 0.0
    var_e_se = 0.0
    try:
        PK = P @ K
        info_kk = 0.5 * np.trace(PK @ PK)
        info_ee = 0.5 * np.trace(P @ P)
        info_ke = 0.5 * np.trace(PK @ P)

        info_matrix = np.array([[info_kk, info_ke], [info_ke, info_ee]])
        info_inv = np.linalg.inv(info_matrix)

        var_k_se = np.sqrt(max(0, info_inv[0, 0]))
        var_e_se = np.sqrt(max(0, info_inv[1, 1]))

        dh2_dvk = var_e / (total_var ** 2)
        dh2_dve = -var_k / (total_var ** 2)
        grad = np.array([dh2_dvk, dh2_dve])
        h2_se = np.sqrt(max(0, grad @ info_inv @ grad))
    except:
        pass

    from scipy.stats import chi2
    h2_pval = 1 - chi2.cdf(h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0
    e_h2 = var_e / total_var if total_var > 0 else 0
    e_h2_pval = 1 - chi2.cdf(e_h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0

    beta = XtV_inv_X_inv @ (X.T @ V_inv @ y)
    beta_se = np.sqrt(np.diag(XtV_inv_X_inv))

    ha_effects = var_k * K @ V_inv @ (y - X @ beta)
    residuals = y - X @ beta - ha_effects

    return {
        'var_k': var_k, 'var_e': var_e,
        'var_k_se': var_k_se, 'var_e_se': var_e_se,
        'h2': h2, 'h2_se': h2_se,
        'h2_pval': h2_pval, 'e_h2_pval': e_h2_pval,
        'beta': beta, 'beta_se': beta_se,
        'ha_effects': ha_effects, 'residuals': residuals,
        'loglik': loglik, 'converged': converged, 'n_iter': n_iter,
        'convergence_history': convergence_history
    }


def he_regression(y, K):
    """Haseman-Elston回归方差组分估计

    参数:
        y: 表型向量
        K: 亲缘关系矩阵 (HA)
    """
    print("Variance components estimation using: HE")
    sys.stdout.flush()

    n = len(y)
    X = np.ones((n, 1))

    print(f"The matrix V has a dimension of {n} x {n}.")
    print("Running ...")
    sys.stdout.flush()

    y_centered = y - np.mean(y)
    yy = np.outer(y_centered, y_centered)
    yy_vec = yy[np.triu_indices(n, k=1)]

    K_vec = K[np.triu_indices(n, k=1)]
    I_vec = np.ones(len(yy_vec))
    K_matrix = np.column_stack([K_vec, I_vec])

    # 初次尝试：无正则化的最小二乘
    try:
        coeffs = np.linalg.lstsq(K_matrix, yy_vec, rcond=None)[0]
    except:
        coeffs = np.zeros(2)

    var_k = max(1e-10, float(coeffs[0]))
    var_e = max(1e-10, float(coeffs[1]))

    # 当残差方差<1e-6时，使用ridge正则化重新估计，避免退化解
    if var_e < 1e-6:
        print(f"  [HE] 检测到退化解 (V(e)={var_e:.2e}), 应用ridge正则化 (lambda=1e-3)...")
        sys.stdout.flush()

        lambda_ridge = 1e-3
        # ridge回归: (K^T K + lambda*I)^-1 K^T y
        KtK = K_matrix.T @ K_matrix
        KtK_ridge = KtK + lambda_ridge * np.eye(KtK.shape[0])
        Kty = K_matrix.T @ yy_vec

        try:
            coeffs_ridge = np.linalg.solve(KtK_ridge, Kty)
        except:
            coeffs_ridge = np.linalg.lstsq(KtK_ridge, Kty, rcond=None)[0]

        var_k = max(1e-10, float(coeffs_ridge[0]))
        var_e = max(1e-10, float(coeffs_ridge[1]))

        print(f"  [HE] Ridge正则化后: V(HA)={var_k:.6f}, V(e)={var_e:.6f}")
        sys.stdout.flush()

    total_var = var_k + var_e
    h2 = var_k / total_var if total_var > 0 else 0

    # 确保h2 < 1.0（防止退化）
    if h2 >= 1.0:
        h2 = 0.999
        var_k = h2 * total_var
        var_e = (1 - h2) * total_var

    # 估计标准误
    h2_se = 0.0
    var_k_se = 0.0
    var_e_se = 0.0
    try:
        V = var_k * K + var_e * np.eye(n)
        V_inv = np.linalg.inv(V)
        V_inv_X = V_inv @ X
        XtV_inv_X = X.T @ V_inv_X
        XtV_inv_X_inv = np.linalg.inv(XtV_inv_X)
        P = V_inv - V_inv_X @ XtV_inv_X_inv @ V_inv_X.T

        PK = P @ K
        info_kk = 0.5 * np.trace(PK @ PK)
        info_ee = 0.5 * np.trace(P @ P)
        info_ke = 0.5 * np.trace(PK @ P)

        info_matrix = np.array([[info_kk, info_ke], [info_ke, info_ee]])
        info_inv = np.linalg.inv(info_matrix)

        var_k_se = np.sqrt(max(0, info_inv[0, 0]))
        var_e_se = np.sqrt(max(0, info_inv[1, 1]))

        dh2_dvk = var_e / (total_var ** 2)
        dh2_dve = -var_k / (total_var ** 2)
        grad = np.array([dh2_dvk, dh2_dve])
        h2_se = np.sqrt(max(0, grad @ info_inv @ grad))
    except:
        pass

    from scipy.stats import chi2
    h2_pval = 1 - chi2.cdf(h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0
    e_h2 = var_e / total_var if total_var > 0 else 0
    e_h2_pval = 1 - chi2.cdf(e_h2 / h2_se**2 if h2_se > 0 else 0, 1) if h2_se > 0 else 1.0

    beta = XtV_inv_X_inv @ (X.T @ V_inv @ y)
    beta_se = np.sqrt(np.diag(XtV_inv_X_inv))

    ha_effects = var_k * K @ V_inv @ (y - X @ beta)
    residuals = y - X @ beta - ha_effects

    print(f"[HE]     h2={h2:.4f}, V(HA)={var_k:.6f}, V(e)={var_e:.6f}")
    sys.stdout.flush()

    return {
        'var_k': var_k, 'var_e': var_e,
        'var_k_se': var_k_se, 'var_e_se': var_e_se,
        'h2': h2, 'h2_se': h2_se,
        'h2_pval': h2_pval, 'e_h2_pval': e_h2_pval,
        'beta': beta, 'beta_se': beta_se,
        'ha_effects': ha_effects, 'residuals': residuals,
        'loglik': None, 'converged': True, 'n_iter': 1,
        'convergence_history': []
    }


def plot_variance_components(result, output_file, args):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Variance Component Estimation (AI-REML)', fontsize=20, fontweight='bold', fontfamily='DejaVu Sans')

    ax1 = axes[0, 0]
    components = ['HA', 'Residual']
    variances = [result['var_k'], result['var_e']]
    var_ses = [result['var_k_se'], result['var_e_se']]
    x = np.arange(len(components))
    bars = ax1.bar(x, variances, width=0.5, color=[PALETTE_A[0], PALETTE_A[1]],
                   edgecolor='#333333', alpha=0.8, yerr=var_ses, capsize=5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(components)
    ax1.set_ylabel('Variance', fontsize=14)
    ax1.set_title('Estimated Variance Components', fontsize=16)
    for bar, val in zip(bars, variances):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.4f}', ha='center', va='bottom', fontsize=12)
    sns.despine(ax=ax1)

    ax2 = axes[0, 1]
    h2 = result['h2']
    h2_se = result['h2_se']
    e_h2 = result['var_e'] / (result['var_k'] + result['var_e'])
    h2_labels = ['HA (h2)', 'Residual (1-h2)']
    h2_values = [h2, e_h2]
    h2_errors = [h2_se, h2_se]
    x = np.arange(len(h2_labels))
    bars = ax2.bar(x, h2_values, width=0.5, color=[PALETTE_A[2], PALETTE_A[3]],
                   edgecolor='#333333', alpha=0.8, yerr=h2_errors, capsize=5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(h2_labels)
    ax2.set_ylabel('Heritability', fontsize=14)
    ax2.set_title('Heritability Estimation', fontsize=16)
    for bar, val in zip(bars, h2_values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.4f}', ha='center', va='bottom', fontsize=12)
    sns.despine(ax=ax2)

    ax3 = axes[1, 0]
    history = result.get('convergence_history', [])
    if history:
        iters = [h['iteration'] for h in history]
        logliks = [h['log_likelihood'] for h in history]
        ax3.plot(iters, logliks, 'o-', color=PALETTE_A[0], linewidth=2, markersize=6)
        ax3.set_xlabel('Iteration', fontsize=14)
        ax3.set_ylabel('Log-Likelihood', fontsize=14)
        ax3.set_title('Convergence History', fontsize=16)
    sns.despine(ax=ax3)

    ax4 = axes[1, 1]
    if history:
        iters = [h['iteration'] for h in history]
        var_ks = [h['var_k'] for h in history]
        var_es = [h['var_e'] for h in history]
        ax4.plot(iters, var_ks, 'o-', color=PALETTE_A[0], linewidth=2, markersize=6, label='V(HA)')
        ax4.plot(iters, var_es, 's-', color=PALETTE_A[1], linewidth=2, markersize=6, label='V(e)')
        ax4.set_xlabel('Iteration', fontsize=14)
        ax4.set_ylabel('Variance', fontsize=14)
        ax4.set_title('Variance Components Trace', fontsize=16)
        leg = ax4.legend(frameon=False)
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')
    sns.despine(ax=ax4)

    plt.tight_layout()
    out_dir = os.path.dirname(output_file)
    basename = os.path.splitext(os.path.basename(output_file))[0]
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(out_dir, f"{basename}.{fmt}"), dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  Variance components plot saved: {os.path.join(out_dir, basename)}.pdf/png")


def plot_breeding_values(result, sample_ids, true_bv, output_file, args):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Breeding Values Distribution', fontsize=20, fontweight='bold', fontfamily='DejaVu Sans')

    ha_effects = result['ha_effects']
    residuals = result['residuals']

    ax1 = axes[0, 0]
    ax1.hist(ha_effects, bins=30, color=PALETTE_A[0], edgecolor='#333333', alpha=0.7)
    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax1.axvline(x=np.mean(ha_effects), color=PALETTE_A[1], linestyle='--', linewidth=2,
               label=f'Mean={np.mean(ha_effects):.3f}')
    ax1.set_xlabel('Estimated HA Effect', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=14)
    ax1.set_title('Estimated HA Effects', fontsize=16)
    leg1 = ax1.legend(frameon=False)
    leg1.get_frame().set_linewidth(0)
    leg1.get_frame().set_facecolor('none')
    leg1.get_frame().set_edgecolor('none')
    sns.despine(ax=ax1)

    ax2 = axes[0, 1]
    if true_bv is not None:
        ax2.hist(true_bv, bins=30, color=PALETTE_A[2], edgecolor='#333333', alpha=0.7)
        ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1)
        ax2.axvline(x=np.mean(true_bv), color=PALETTE_A[1], linestyle='--', linewidth=2,
                   label=f'Mean={np.mean(true_bv):.3f}')
        ax2.set_xlabel('True Breeding Value', fontsize=14)
        ax2.set_ylabel('Frequency', fontsize=14)
        ax2.set_title('True Breeding Values', fontsize=16)
        leg2 = ax2.legend(frameon=False)
        leg2.get_frame().set_linewidth(0)
        leg2.get_frame().set_facecolor('none')
        leg2.get_frame().set_edgecolor('none')
    else:
        sorted_bv = np.sort(ha_effects)[::-1]
        top_n = min(20, len(sorted_bv))
        ax2.barh(range(top_n), sorted_bv[:top_n], color=PALETTE_A[2], edgecolor='#333333', alpha=0.7)
        ax2.set_xlabel('Estimated Breeding Value', fontsize=14)
        ax2.set_ylabel('Rank', fontsize=14)
        ax2.set_title('Top 20 Breeding Values', fontsize=16)
        ax2.invert_yaxis()
    sns.despine(ax=ax2)

    ax3 = axes[1, 0]
    if true_bv is not None:
        ax3.scatter(true_bv, ha_effects, alpha=0.5, s=20, c=PALETTE_A[0],
                   edgecolors='#333333', linewidth=0.3)
        lims_min = min(true_bv.min(), ha_effects.min())
        lims_max = max(true_bv.max(), ha_effects.max())
        ax3.plot([lims_min, lims_max], [lims_min, lims_max], color='#333333', linestyle='--', linewidth=1.5, label='y=x')
        corr = np.corrcoef(true_bv, ha_effects)[0, 1]
        ax3.set_xlabel('True Breeding Value', fontsize=14)
        ax3.set_ylabel('Estimated HA Effect', fontsize=14)
        ax3.set_title(f'True vs Estimated (r={corr:.4f})', fontsize=16)
        leg3 = ax3.legend(frameon=False)
        leg3.get_frame().set_linewidth(0)
        leg3.get_frame().set_facecolor('none')
        leg3.get_frame().set_edgecolor('none')
    else:
        ax3.scatter(ha_effects, residuals, alpha=0.5, s=20, c=PALETTE_A[0],
                   edgecolors='#333333', linewidth=0.3)
        ax3.axhline(y=0, color='gray', linestyle='--', linewidth=1)
        ax3.set_xlabel('Estimated Breeding Value', fontsize=14)
        ax3.set_ylabel('Residual', fontsize=14)
        ax3.set_title('Breeding Value vs Residual', fontsize=16)
    sns.despine(ax=ax3)

    ax4 = axes[1, 1]
    ax4.hist(residuals, bins=30, color=PALETTE_A[4], edgecolor='#333333', alpha=0.7)
    ax4.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax4.set_xlabel('Residual', fontsize=14)
    ax4.set_ylabel('Frequency', fontsize=14)
    ax4.set_title('Residuals Distribution', fontsize=16)
    sns.despine(ax=ax4)

    plt.tight_layout()
    out_dir = os.path.dirname(output_file)
    basename = os.path.splitext(os.path.basename(output_file))[0]
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(out_dir, f"{basename}.{fmt}"), dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  Breeding values plot saved: {os.path.join(out_dir, basename)}.pdf/png")


def plot_variance_pie(result, output_file, args):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Variance Component Proportions', fontsize=18, fontweight='bold', fontfamily='DejaVu Sans')

    labels = ['HA', 'Residual']
    values = [result['var_k'], result['var_e']]
    total = sum(values)
    percentages = [v / total * 100 for v in values]
    colors = [PALETTE_A[0], PALETTE_A[1]]

    ax1 = axes[0]
    wedges, texts, autotexts = ax1.pie(values, labels=labels, autopct='%1.1f%%',
                                        colors=colors, startangle=90, textprops={'fontsize': 12})
    ax1.set_title('Variance Components', fontsize=16)

    ax2 = axes[1]
    bars = ax2.barh(labels, percentages, color=colors, edgecolor='#333333', alpha=0.8)
    ax2.set_xlabel('Percentage (%)', fontsize=14)
    ax2.set_title('Variance Proportions', fontsize=16)

    for bar, pct in zip(bars, percentages):
        ax2.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{pct:.1f}%', va='center', fontsize=12)

    sns.despine(ax=ax2)

    plt.tight_layout()
    out_dir = os.path.dirname(output_file)
    basename = os.path.splitext(os.path.basename(output_file))[0]
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(out_dir, f"{basename}.{fmt}"), dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  Variance pie plot saved: {os.path.join(out_dir, basename)}.pdf/png")


def save_vars(result, output_file):
    with open(output_file, 'w') as f:
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        f.write(f"HA\t{result['var_k']:.6f}\t{result['var_k_se']:.6f}\t"
                f"{result['h2']:.6e}\t{result['h2_se']:.6f}\t{result['h2_pval']:.6e}\n")
        e_h2 = result['var_e'] / (result['var_k'] + result['var_e'])
        f.write(f"e\t{result['var_e']:.6f}\t{result['var_e_se']:.6f}\t"
                f"{e_h2:.6f}\t{result['h2_se']:.6f}\t{result['e_h2_pval']:.6e}\n")
    print(f"Results of estimated variance components have been saved in the file [{output_file}].")


def save_beta(result, output_file):
    with open(output_file, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        f.write(f"mu\t{result['beta'][0]:.8g}\t{result['beta_se'][0]:.6g}\n")
    print(f"Coefficients of all covariates and fixed effects are saved in file [{output_file}].")


def save_rand(sample_ids, result, output_file):
    with open(output_file, 'w') as f:
        f.write("ID\tHA\tresiduals\n")
        for i, sid in enumerate(sample_ids):
            f.write(f"{sid}\t{result['ha_effects'][i]:.8g}\t{result['residuals'][i]:.6f}\n")
    print(f"Random effects of all individuals are saved in file [{output_file}].")


def main():
    args = parse_args()
    plt.rcParams.update({
        'font.size': args.font_size,
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
    })
    bfile = args.bfile
    pheno_file = args.pheno
    pedigree_file = args.pedigree
    output_dir = args.out
    g_weight = args.g_weight

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("Single Trait Model (Replicating HIBLUP --single-trait)")
    print(f"Model: Phenotype = mu + HA + e")
    print(f"BFILE: {bfile}")
    print(f"Phenotype: {pheno_file}")
    print(f"Pedigree: {pedigree_file}")
    print(f"Output: {output_dir}")
    print(f"G weight: {g_weight}")
    print("=" * 70)
    sys.stdout.flush()

    print("\n[1/5] Reading PLINK genotype data...")
    sys.stdout.flush()
    M, sample_ids, snp_ids = read_plink_all(bfile)

    print(f"\n[2/5] Computing GRM (VanRaden)...")
    sys.stdout.flush()
    G = compute_grm_vanraden(M)

    print(f"\n[3/5] Computing pedigree A matrix...")
    sys.stdout.flush()
    A = compute_pedigree_A(pedigree_file, sample_ids)

    print(f"\n[4/5] Constructing HA matrix and reading phenotype...")
    sys.stdout.flush()
    HA = construct_HA_matrix(G, A, w=g_weight)

    pheno_df = pd.read_csv(pheno_file, sep='\t')
    pheno_df.columns = [c.strip() for c in pheno_df.columns]

    if len(pheno_df.columns) >= 3 and 'BreedingValue' in pheno_df.columns:
        true_bv_col = 'BreedingValue'
    elif len(pheno_df.columns) >= 3:
        true_bv_col = pheno_df.columns[2]
    else:
        true_bv_col = None

    pheno_df.rename(columns={pheno_df.columns[0]: 'ID', pheno_df.columns[1]: 'Phenotype'}, inplace=True)

    sample_ids_int = [int(s) for s in sample_ids]
    pheno_df['ID'] = pheno_df['ID'].astype(int)

    pheno_subset = pheno_df[pheno_df['ID'].isin(sample_ids_int)].copy()
    pheno_subset = pheno_subset.set_index('ID').loc[sample_ids_int].reset_index()

    y = pheno_subset['Phenotype'].values
    true_bv = pheno_subset[true_bv_col].values if true_bv_col and true_bv_col in pheno_subset.columns else None

    print(f"  {len(y)} effective records without missing will be used for analysis.")
    print(f"  Phenotype mean: {np.mean(y):.6f}, std: {np.std(y):.6f}")
    sys.stdout.flush()

    print(f"\n[5/5] Fitting single trait model (AI-REML)...")
    sys.stdout.flush()
    result = ai_reml(y, HA, max_iter=20)

    print(f"\nSummary of estimated variances and heritability:")
    print(f"       Var Var_SE       h2  h2_SE h2_Pr(Chisq)")
    print(f" HA {result['var_k']:.4f} {result['var_k_se']:.4f} {result['h2']:.2e} {result['h2_se']:.4f} {result['h2_pval']:.6e}")
    e_h2 = result['var_e'] / (result['var_k'] + result['var_e'])
    print(f"  e {result['var_e']:.4f} {result['var_e_se']:.4f} {e_h2:.4f} {result['h2_se']:.4f} {result['e_h2_pval']:.2e}")

    print(f"\nEstimate fixed effects and random effects ...")
    print(f"Summary of estimated coefficients:")
    print(f" Levels Estimation     SE")
    print(f"     mu {result['beta'][0]:.8g} {result['beta_se'][0]:.4f}")
    sys.stdout.flush()

    # EM-REML 方差组分估计
    print(f"\n--- EM-REML (max_iter={args.em_max_iter}) ---")
    sys.stdout.flush()
    result_em = em_reml(y, HA, max_iter=args.em_max_iter)

    # HE 回归方差组分估计
    print(f"\n--- HE Regression ---")
    sys.stdout.flush()
    result_he = he_regression(y, HA)

    vars_file = os.path.join(output_dir, "single_trait.vars")
    beta_file = os.path.join(output_dir, "single_trait.beta")
    rand_file = os.path.join(output_dir, "single_trait.rand")

    save_vars(result, vars_file)
    save_beta(result, beta_file)
    save_rand(sample_ids, result, rand_file)

    print(f"\n[6/6] Generating visualization plots...")
    sys.stdout.flush()
    plot_variance_components(result, os.path.join(output_dir, "variance_components"), args)
    plot_breeding_values(result, sample_ids, true_bv, os.path.join(output_dir, "breeding_values"), args)
    plot_variance_pie(result, os.path.join(output_dir, "variance_pie"), args)

    # Bilingual summary files (AI-REML 为主结果)
    en_summary_path = os.path.join(output_dir, "single_trait_summary.txt")
    zh_summary_path = os.path.join(output_dir, "single_trait_summary-zh.txt")
    e_h2 = result['var_e'] / (result['var_k'] + result['var_e'])
    with open(en_summary_path, 'w') as f:
        f.write("Single Trait Model Summary\n")
        f.write(f"Model: Phenotype = mu + HA + e\n")
        f.write(f"V(HA): {result['var_k']:.6f} (SE={result['var_k_se']:.6f})\n")
        f.write(f"V(e): {result['var_e']:.6f} (SE={result['var_e_se']:.6f})\n")
        f.write(f"Heritability (h2): {result['h2']:.6f} (SE={result['h2_se']:.6f})\n")
        f.write(f"mu: {result['beta'][0]:.8g} (SE={result['beta_se'][0]:.6g})\n")
        f.write(f"Converged: {result['converged']} in {result['n_iter']} iterations\n")
        f.write(f"Output: variance_components.pdf/png, breeding_values.pdf/png, variance_pie.pdf/png\n")
    with open(zh_summary_path, 'w') as f:
        f.write("单性状模型分析摘要\n")
        f.write("模型: 表型 = mu + HA + e\n")
        f.write(f"V(HA): {result['var_k']:.6f} (SE={result['var_k_se']:.6f})\n")
        f.write(f"V(e): {result['var_e']:.6f} (SE={result['var_e_se']:.6f})\n")
        f.write(f"遗传力 (h2): {result['h2']:.6f} (SE={result['h2_se']:.6f})\n")
        f.write(f"mu: {result['beta'][0]:.8g} (SE={result['beta_se'][0]:.6g})\n")
        f.write(f"收敛: {result['converged']}, 迭代{result['n_iter']}次\n")
        f.write("输出: variance_components.pdf/png, breeding_values.pdf/png, variance_pie.pdf/png\n")

    # 三方法对比摘要 (包含AI-REML, EM-REML, HE的方差组分、遗传力、迭代次数)
    allmethods_summary_path = os.path.join(output_dir, "variance_component_summary.txt")
    all_results = {'AI_REML': result, 'EM_REML': result_em, 'HE': result_he}
    with open(allmethods_summary_path, 'w') as f:
        f.write("单性状模型方差组分分析摘要 (三方法对比)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"样本数: {len(y)}\n")
        f.write(f"G weight: {g_weight}\n\n")
        for method, res in all_results.items():
            f.write(f"{method}方法:\n")
            f.write(f"  收敛: {res.get('converged', 'N/A')}\n")
            # 确保EM-REML迭代次数被正确记录到摘要文件
            f.write(f"  迭代次数: {res.get('n_iter', 'N/A')}\n")
            if res.get('loglik') is not None:
                f.write(f"  对数似然: {res['loglik']:.4f}\n")
            f.write("  方差组分:\n")
            f.write(f"    V(HA): {res['var_k']:.6f} (SE={res['var_k_se']:.6f})\n")
            f.write(f"    V(e): {res['var_e']:.6f} (SE={res['var_e_se']:.6f})\n")
            f.write("  遗传力:\n")
            f.write(f"    h2: {res['h2']:.6f} (SE={res['h2_se']:.6f})\n")
            e_h2_val = res['var_e'] / (res['var_k'] + res['var_e']) if (res['var_k'] + res['var_e']) > 0 else 0
            f.write(f"    e_h2: {e_h2_val:.6f}\n")
            f.write("\n")
    print(f"  三方法对比摘要已保存: {allmethods_summary_path}")

    print(f"\nAnalysis finished!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. {vars_file}")
    print(f"  2. {beta_file}")
    print(f"  3. {rand_file}")
    print(f"  4. {os.path.join(output_dir, 'variance_components')}.pdf/png")
    print(f"  5. {os.path.join(output_dir, 'breeding_values')}.pdf/png")
    print(f"  6. {os.path.join(output_dir, 'variance_pie')}.pdf/png")
    print(f"  7. {en_summary_path}")
    print(f"  8. {zh_summary_path}")
    print(f"  9. {allmethods_summary_path}")


if __name__ == "__main__":
    main()
