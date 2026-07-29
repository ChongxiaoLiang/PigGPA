#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
重复记录模型脚本 - 与HIBLUP/8输出文件构成一致
模型: weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e
输出文件:
  - repeated_model.vars  方差组分
  - repeated_model.beta  固定效应
  - repeated_model.anova 方差分析表
  - repeated_model.rand  随机效应
  - repeated_model.log   分析日志
  - variance_components_bar.pdf  方差组分图
  - breeding_value_distribution.pdf  育种值分布图
  - fixed_effects_forest.pdf  固定效应森林图

用法:
  python repeated_records_model.py \
    --bfile /path/to/simulated_population \
    --pheno /path/to/phenotype_long.txt \
    --snp-file /path/to/chr1_snps.txt \
    --sample-file /path/to/train_samples.txt \
    --out /path/to/output_dir
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import cholesky, cho_solve
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap
PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)

import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 10
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Repeated Records Model (Replicating HIBLUP --repeat-trait)',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--bfile', required=True,
                        help='PLINK bed/bim/fam prefix path')
    parser.add_argument('--pheno', required=True,
                        help='Phenotype file path (long format, tab-separated)')
    parser.add_argument('--snp-file', required=True,
                        help='SNP list file path (one SNP ID per line)')
    parser.add_argument('--sample-file', required=True,
                        help='Sample list file path (one sample ID per line)')
    parser.add_argument('--out', required=True,
                        help='Output directory')
    parser.add_argument('--font-size', type=int, default=12, help='Font size for figures (default: 12)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for figure output (default: 300)')
    parser.add_argument('--pheno-col', default='weight', help='Phenotype column name in phenotype file (default: weight)')
    return parser.parse_args()


def read_plink(bfile, sample_ids, snp_ids):
    bim = pd.read_csv(f"{bfile}.bim", sep='\t', header=None,
                      names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
    fam = pd.read_csv(f"{bfile}.fam", sep='\t', header=None,
                      names=['fid', 'iid', 'pid', 'mid', 'sex', 'pheno'])

    snp_set = set(snp_ids)
    snp_mask = bim['snp'].isin(snp_set)
    bim_filtered = bim[snp_mask]
    snp_indices = bim_filtered.index.tolist()
    chr1_bim = bim_filtered.reset_index(drop=True)
    n_snps = len(snp_indices)
    print(f"  SNP count: {n_snps}")

    sample_iids = [int(x) for x in sample_ids]
    sample_mask = fam['iid'].isin(sample_iids)
    sample_indices = fam[sample_mask].index.tolist()
    fam_subset = fam[sample_mask].reset_index(drop=True)
    n_samples = len(sample_indices)
    n_samples_total = len(fam)
    print(f"  Sample count: {n_samples}")

    n_bytes_per_snp = (n_samples_total + 3) // 4
    genotypes = np.zeros((n_snps, n_samples), dtype=np.int8)

    with open(f"{bfile}.bed", 'rb') as f:
        magic = f.read(3)
        if magic != b'\x6c\x1b\x01':
            raise ValueError("Invalid PLINK bed file")
        for i, snp_idx in enumerate(snp_indices):
            f.seek(3 + snp_idx * n_bytes_per_snp)
            bytes_data = f.read(n_bytes_per_snp)
            all_geno = np.zeros(n_samples_total, dtype=np.int8)
            for j, byte in enumerate(bytes_data):
                for k in range(4):
                    idx = j * 4 + k
                    if idx < n_samples_total:
                        all_geno[idx] = (byte >> k) & 3
            genotypes[i, :] = all_geno[sample_indices]

    geno_codes = np.array([0, 1, -1, 2], dtype=np.float64)
    G = geno_codes[genotypes]

    G_valid = np.where(G == -1, np.nan, G)
    n_valid = np.sum(~np.isnan(G_valid), axis=1)
    p = np.nansum(G, axis=1) / (2 * n_valid)
    p = np.where(np.isnan(p), 0, p)

    for i in range(G.shape[0]):
        valid = ~np.isnan(G_valid[i, :])
        if np.sum(valid) > 0 and np.sum(~valid) > 0:
            G[i, ~valid] = 2 * p[i]

    G = G.T
    return G, chr1_bim, fam_subset, p


def calc_grm(G, p):
    Z = G - 2 * p
    sum_2pq = 2 * np.sum(p * (1 - p))
    GRM = (Z @ Z.T) / sum_2pq
    diag_mean = np.mean(np.diag(GRM))
    GRM = GRM / diag_mean
    return GRM, Z


def build_design_matrices(pheno_long, unique_ids, pheno_col='weight'):
    n_obs = len(pheno_long)
    n_ind = len(unique_ids)
    id_to_idx = {id_: i for i, id_ in enumerate(unique_ids)}

    X = np.zeros((n_obs, 5))
    X[:, 0] = 1
    X[:, 1] = (pheno_long['sex'] == 'M').astype(float).values

    season_dummies = pd.get_dummies(pheno_long['season'], prefix='season')
    for i, col in enumerate(['season_Spring', 'season_Summer', 'season_Winter']):
        if col in season_dummies.columns:
            X[:, 2 + i] = season_dummies[col].values

    Z = np.zeros((n_obs, n_ind))
    row_idx = np.arange(n_obs)
    col_idx = np.array([id_to_idx[id_] for id_ in pheno_long['ID'].values])
    Z[row_idx, col_idx] = 1

    y = pheno_long[pheno_col].values.astype(float)
    return y, X, Z


def ai_reml(y, X, Z, G, max_iter=20, tol=1e-8):
    n = len(y)
    q = Z.shape[1]
    I_q = np.eye(q)

    ZtZ = Z.T @ Z
    ZtX = Z.T @ X
    Zty = Z.T @ y

    var_y = np.var(y)
    var_pe = var_y * 0.5
    var_ga = var_y * 0.3
    var_e = var_y * 0.2

    log_lines = []
    log_lines.append(f"Total {n} records will be predicted.")
    log_lines.append(f"Total 3 variance components need to be estimated.")
    log_lines.append(f"Variance components estimation using: AI({max_iter})")
    log_lines.append(f"The matrix V has a dimension of {n} x {n}.")
    log_lines.append("Running ...")
    log_lines.append("Alg.\tIter.\tLogL.\tV(ID)\tV(GA)\tV(e)")

    print(f"  Initial: V(ID)={var_pe:.5f}, V(GA)={var_ga:.5f}, V(e)={var_e:.5f}")

    for iteration in range(max_iter):
        Sigma = var_pe * I_q + var_ga * G
        Sigma_inv = np.linalg.inv(Sigma)
        M = Sigma_inv + ZtZ / var_e
        W = np.linalg.inv(M)
        WZtZ = W @ ZtZ
        WZtX = W @ ZtX

        Vinv_y = (y - Z @ (W @ Zty) / var_e) / var_e
        Vinv_X = (X - Z @ WZtX / var_e) / var_e

        XtVX = X.T @ Vinv_X
        XtVy = X.T @ Vinv_y
        try:
            XtVX_inv = np.linalg.inv(XtVX)
        except:
            XtVX_inv = np.linalg.pinv(XtVX)
        beta = XtVX_inv @ XtVy

        r = y - X @ beta
        Ztr = Zty - ZtX @ beta
        Py = (r - Z @ (W @ Ztr) / var_e) / var_e

        try:
            L = cholesky(Sigma, lower=True)
            logdet_Sigma = 2 * np.sum(np.log(np.diag(L)))
        except:
            sign, logdet_Sigma = np.linalg.slogdet(Sigma)

        logdet_W = -np.linalg.slogdet(M)[1]
        logdet_V = logdet_Sigma - logdet_W + n * np.log(var_e)
        logL = -0.5 * (n * np.log(2 * np.pi) + logdet_V + r @ Py)

        ZtPy = (Ztr - WZtZ @ Ztr / var_e) / var_e

        K_pe_Py = Z @ ZtPy
        K_ga_Py = Z @ (G @ ZtPy)

        Vinv_K_pe_Py = Z @ (ZtPy - WZtZ @ ZtPy / var_e) / var_e
        PK_pe_Py = Vinv_K_pe_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_K_pe_Py))

        Vinv_K_ga_Py = Z @ (G @ ZtPy - WZtZ @ G @ ZtPy / var_e) / var_e
        PK_ga_Py = Vinv_K_ga_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_K_ga_Py))

        Vinv_Py = (Py - Z @ (W @ ZtPy) / var_e) / var_e
        PK_e_Py = Vinv_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_Py))

        ai_pe_pe = 0.5 * PK_pe_Py @ K_pe_Py
        ai_pe_ga = 0.5 * PK_pe_Py @ K_ga_Py
        ai_pe_e = 0.5 * PK_pe_Py @ Py
        ai_ga_ga = 0.5 * PK_ga_Py @ K_ga_Py
        ai_ga_e = 0.5 * PK_ga_Py @ Py
        ai_e_e = 0.5 * PK_e_Py @ Py

        ai_mat = np.array([
            [ai_pe_pe, ai_pe_ga, ai_pe_e],
            [ai_pe_ga, ai_ga_ga, ai_ga_e],
            [ai_pe_e, ai_ga_e, ai_e_e]
        ])

        ZtZ2 = ZtZ @ ZtZ
        WZtZ2 = W @ ZtZ2
        tr_Vinv_K_pe = (np.trace(ZtZ) - np.trace(WZtZ2) / var_e) / var_e
        GZtZ = G @ ZtZ
        WGZtZ2 = W @ (G @ ZtZ2)
        tr_Vinv_K_ga = (np.trace(GZtZ) - np.trace(WGZtZ2) / var_e) / var_e
        tr_Vinv = (n - np.trace(WZtZ) / var_e) / var_e

        ZtA = (ZtX - ZtZ @ WZtX / var_e) / var_e
        tr_XtVX_inv_K_pe = np.trace(XtVX_inv @ (ZtA.T @ ZtA))
        tr_XtVX_inv_K_ga = np.trace(XtVX_inv @ (ZtA.T @ G @ ZtA))
        tr_XtVX_inv_K_e = np.trace(XtVX_inv @ (Vinv_X.T @ Vinv_X))

        tr_PK_pe = tr_Vinv_K_pe - tr_XtVX_inv_K_pe
        tr_PK_ga = tr_Vinv_K_ga - tr_XtVX_inv_K_ga
        tr_PK_e = tr_Vinv - tr_XtVX_inv_K_e

        Py_K_pe_Py = Py @ K_pe_Py
        Py_K_ga_Py = Py @ K_ga_Py
        Py_Py = Py @ Py

        score_pe = 0.5 * (Py_K_pe_Py - tr_PK_pe)
        score_ga = 0.5 * (Py_K_ga_Py - tr_PK_ga)
        score_e = 0.5 * (Py_Py - tr_PK_e)

        score_vec = np.array([score_pe, score_ga, score_e])

        try:
            delta = np.linalg.solve(ai_mat, score_vec)
        except:
            delta = np.linalg.lstsq(ai_mat, score_vec, rcond=None)[0]

        new_var_pe = var_pe + delta[0]
        new_var_ga = var_ga + delta[1]
        new_var_e = var_e + delta[2]

        step = 1.0
        halvings = 0
        while (new_var_pe < 0.01 or new_var_ga < 0.01 or new_var_e < 0.01 or
               new_var_pe > 200 or new_var_ga > 200 or new_var_e > 200) and halvings < 10:
            step *= 0.5
            new_var_pe = var_pe + step * delta[0]
            new_var_ga = var_ga + step * delta[1]
            new_var_e = var_e + step * delta[2]
            halvings += 1

        log_line = f"[AI]\t{iteration+1}\t{logL:.2f}\t{new_var_pe:.5f}\t{new_var_ga:.5f}\t{new_var_e:.5f}"
        log_lines.append(log_line)
        print(f"  [AI] {iteration+1}\t{logL:.2f}\t{new_var_pe:.5f}\t{new_var_ga:.5f}\t{new_var_e:.5f}")

        abs_change = max(abs(new_var_pe - var_pe), abs(new_var_ga - var_ga), abs(new_var_e - var_e))
        if abs_change < tol and iteration >= 5:
            print(f"  Converged!")
            var_pe, var_ga, var_e = new_var_pe, new_var_ga, new_var_e
            break

        var_pe, var_ga, var_e = new_var_pe, new_var_ga, new_var_e

    converged = abs_change < tol and iteration >= 5
    if not converged:
        log_lines.append("[Converged?] No! More iterations are required.")
    else:
        log_lines.append("[Converged?] Yes!")

    print(f"  Final: V(ID)={var_pe:.6f}, V(GA)={var_ga:.6f}, V(e)={var_e:.6f}")

    Sigma = var_pe * I_q + var_ga * G
    Sigma_inv = np.linalg.inv(Sigma)
    M = Sigma_inv + ZtZ / var_e
    W = np.linalg.inv(M)
    WZtZ = W @ ZtZ
    WZtX = W @ ZtX

    Vinv_y = (y - Z @ (W @ Zty) / var_e) / var_e
    Vinv_X = (X - Z @ WZtX / var_e) / var_e

    XtVX = X.T @ Vinv_X
    try:
        XtVX_inv = np.linalg.inv(XtVX)
    except:
        XtVX_inv = np.linalg.pinv(XtVX)
    beta = XtVX_inv @ (X.T @ Vinv_y)
    beta_se = np.sqrt(np.diag(XtVX_inv))

    r = y - X @ beta
    Ztr = Zty - ZtX @ beta
    ZtPy = (Ztr - WZtZ @ Ztr / var_e) / var_e

    pe_hat = var_pe * ZtPy
    ga_hat = var_ga * G @ ZtPy

    fitted = X @ beta + Z @ pe_hat + Z @ ga_hat
    residuals = y - fitted

    Py = (r - Z @ (W @ Ztr) / var_e) / var_e

    K_pe_Py = Z @ ZtPy
    K_ga_Py = Z @ (G @ ZtPy)

    Vinv_K_pe_Py = Z @ (ZtPy - WZtZ @ ZtPy / var_e) / var_e
    PK_pe_Py = Vinv_K_pe_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_K_pe_Py))

    Vinv_K_ga_Py = Z @ (G @ ZtPy - WZtZ @ G @ ZtPy / var_e) / var_e
    PK_ga_Py = Vinv_K_ga_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_K_ga_Py))

    Vinv_Py = (Py - Z @ (W @ ZtPy) / var_e) / var_e
    PK_e_Py = Vinv_Py - Vinv_X @ (XtVX_inv @ (X.T @ Vinv_Py))

    ai_pe_pe = 0.5 * PK_pe_Py @ K_pe_Py
    ai_pe_ga = 0.5 * PK_pe_Py @ K_ga_Py
    ai_pe_e = 0.5 * PK_pe_Py @ Py
    ai_ga_ga = 0.5 * PK_ga_Py @ K_ga_Py
    ai_ga_e = 0.5 * PK_ga_Py @ Py
    ai_e_e = 0.5 * PK_e_Py @ Py

    ai_mat_final = np.array([
        [ai_pe_pe, ai_pe_ga, ai_pe_e],
        [ai_pe_ga, ai_ga_ga, ai_ga_e],
        [ai_pe_e, ai_ga_e, ai_e_e]
    ])

    try:
        cov_mat = np.linalg.inv(ai_mat_final)
    except:
        cov_mat = np.linalg.pinv(ai_mat_final)

    se_var_pe = np.sqrt(max(0, cov_mat[0, 0]))
    se_var_ga = np.sqrt(max(0, cov_mat[1, 1]))
    se_var_e = np.sqrt(max(0, cov_mat[2, 2]))

    total_var = var_pe + var_ga + var_e
    h2_pe = var_pe / total_var
    h2_ga = var_ga / total_var
    h2_e = var_e / total_var

    Vp2 = total_var ** 2
    grad_pe = np.array([(var_ga + var_e) / Vp2, -var_pe / Vp2, -var_pe / Vp2])
    grad_ga = np.array([-var_ga / Vp2, (var_pe + var_e) / Vp2, -var_ga / Vp2])
    grad_e = np.array([-var_e / Vp2, -var_e / Vp2, (var_pe + var_ga) / Vp2])

    se_h2_pe = np.sqrt(max(0, grad_pe @ cov_mat @ grad_pe))
    se_h2_ga = np.sqrt(max(0, grad_ga @ cov_mat @ grad_ga))
    se_h2_e = np.sqrt(max(0, grad_e @ cov_mat @ grad_e))

    wald_pe = (var_pe / se_var_pe) ** 2 if se_var_pe > 0 else 0
    wald_ga = (var_ga / se_var_ga) ** 2 if se_var_ga > 0 else 0
    wald_e = (var_e / se_var_e) ** 2 if se_var_e > 0 else 0
    p_value_pe = 0.5 * stats.chi2.sf(wald_pe, df=1)
    p_value_ga = 0.5 * stats.chi2.sf(wald_ga, df=1)
    p_value_e = 0.5 * stats.chi2.sf(wald_e, df=1)

    return (var_pe, var_ga, var_e,
            se_var_pe, se_var_ga, se_var_e,
            h2_pe, h2_ga, h2_e,
            se_h2_pe, se_h2_ga, se_h2_e,
            p_value_pe, p_value_ga, p_value_e,
            beta, beta_se, pe_hat, ga_hat, residuals,
            log_lines)


def compute_anova(y, X, var_e):
    n = len(y)
    p_x = X.shape[1]

    X_mu = X[:, :1]
    beta_mu = np.linalg.lstsq(X_mu, y, rcond=None)[0]
    rss_mu = np.sum((y - X_mu @ beta_mu) ** 2)

    X_mu_sex = X[:, :2]
    beta_mu_sex = np.linalg.lstsq(X_mu_sex, y, rcond=None)[0]
    rss_mu_sex = np.sum((y - X_mu_sex @ beta_mu_sex) ** 2)

    beta_full = np.linalg.lstsq(X, y, rcond=None)[0]
    rss_full = np.sum((y - X @ beta_full) ** 2)

    ss_sex = rss_mu - rss_mu_sex
    ss_season = rss_mu_sex - rss_full

    df_sex = 1
    df_season = 3
    df_e = n - p_x

    ms_sex = ss_sex / df_sex
    ms_season = ss_season / df_season

    ms_e = var_e
    ss_e = var_e * df_e

    f_sex = ms_sex / ms_e
    f_season = ms_season / ms_e

    p_sex = stats.f.sf(f_sex, df_sex, df_e)
    p_season = stats.f.sf(f_season, df_season, df_e)

    return {
        'sex': {'df': df_sex, 'ss': ss_sex, 'ms': ms_sex, 'f': f_sex, 'p': p_sex},
        'season': {'df': df_season, 'ss': ss_season, 'ms': ms_season, 'f': f_season, 'p': p_season},
        'e': {'df': df_e, 'ss': ss_e, 'ms': ms_e}
    }


def save_vars(filepath, var_pe, var_ga, var_e, se_var_pe, se_var_ga, se_var_e,
              h2_pe, h2_ga, h2_e, se_h2_pe, se_h2_ga, se_h2_e,
              p_value_pe, p_value_ga, p_value_e):
    with open(filepath, 'w') as f:
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        f.write(f"ID\t{var_pe:.6g}\t{se_var_pe:.6g}\t{h2_pe:.6g}\t{se_h2_pe:.6g}\t{p_value_pe:.6g}\n")
        f.write(f"GA\t{var_ga:.6g}\t{se_var_ga:.6g}\t{h2_ga:.6g}\t{se_h2_ga:.6g}\t{p_value_ga:.6g}\n")
        f.write(f"e\t{var_e:.6g}\t{se_var_e:.6g}\t{h2_e:.6g}\t{se_h2_e:.6g}\t{p_value_e:.6g}\n")


def save_beta(filepath, beta, beta_se):
    level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
    with open(filepath, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        for name, est, se in zip(level_names, beta, beta_se):
            f.write(f"{name}\t{est:.6g}\t{se:.6g}\n")


def save_rand(filepath, pheno_long, unique_ids, pe_hat, ga_hat, residuals):
    id_to_idx = {id_: i for i, id_ in enumerate(unique_ids)}
    with open(filepath, 'w') as f:
        f.write("ID\tID\tGA\tresiduals\n")
        for idx, row in pheno_long.iterrows():
            i = id_to_idx[row['ID']]
            f.write(f"{row['ID']}\t{pe_hat[i]:.6g}\t{ga_hat[i]:.6g}\t{residuals[idx]:.6g}\n")


def save_anova(filepath, anova_results):
    with open(filepath, 'w') as f:
        f.write("Factors\tDf\tSumSq\tMeanSq\tF\tPr(>F)\n")
        r = anova_results['sex']
        p_str = "0" if r['p'] < 1e-15 else f"{r['p']:.6g}"
        f.write(f"sex\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t{r['f']:.6g}\t{p_str}\n")
        r = anova_results['season']
        p_str = "0" if r['p'] < 1e-15 else f"{r['p']:.6g}"
        f.write(f"season\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t{r['f']:.6g}\t{p_str}\n")
        r = anova_results['e']
        f.write(f"e\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t.\t.\n")


def save_log(filepath, log_lines, var_pe, var_ga, var_e,
             se_var_pe, se_var_ga, se_var_e,
             h2_pe, h2_ga, h2_e,
             se_h2_pe, se_h2_ga, se_h2_e,
             p_value_pe, p_value_ga, p_value_e,
             beta, beta_se, n_snps, n_individuals):
    with open(filepath, 'w') as f:
        f.write("#=============================================================#\n")
        f.write("#                    REPEATED RECORDS MODEL                    #\n")
        f.write("#=============================================================#\n\n")
        f.write("Model:\n")
        f.write("  weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e\n\n")
        f.write("Method: AI-REML (Woodbury identity)\n")
        f.write(f"GRM: VanRaden method, {n_snps} markers, {n_individuals} individuals\n\n")
        for line in log_lines:
            f.write(line + "\n")
        f.write(f"\nSummary of estimated variances and heritability:\n")
        f.write(f"        Var Var_SE     h2  h2_SE h2_Pr(Chisq)\n")
        f.write(f" ID {var_pe:.4f} {se_var_pe:.4f} {h2_pe:.4f} {se_h2_pe:.4f} {p_value_pe:.6g}\n")
        f.write(f" GA {var_ga:.4f} {se_var_ga:.4f} {h2_ga:.4f} {se_h2_ga:.4f} {p_value_ga:.6g}\n")
        f.write(f"  e {var_e:.4f} {se_var_e:.4f} {h2_e:.4f} {se_h2_e:.4f} {p_value_e:.6g}\n")
        f.write(f"\nEstimate fixed effects and random effects ...\n")
        f.write(f"Summary of estimated coefficients:\n")
        f.write(f"        Levels Estimation     SE\n")
        level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
        for name, est, se in zip(level_names, beta, beta_se):
            f.write(f"  {name:>15s} {est:10.4f} {se:.4f}\n")
        f.write(f"\nAnalysis finished.\n")


def create_visualizations(var_pe, var_ga, var_e, se_var_pe, se_var_ga, se_var_e,
                          h2_pe, h2_ga, h2_e,
                          beta, beta_se, ga_hat, output_dir, args):
    fs = args.font_size

    # ---- Figure 1: Variance Components Bar ----
    fig, ax = plt.subplots(figsize=(10, 6))
    variances = [var_pe, var_ga, var_e]
    se_values = [se_var_pe, se_var_ga, se_var_e]
    labels = ['Permanent\nEnvironment\n(ID)', 'Additive\nGenetic\n(GA)', 'Residual\n(e)']
    colors = [PALETTE_A[0], PALETTE_A[2], PALETTE_A[7]]

    bars = ax.bar(labels, variances, color=colors, edgecolor='#333333', linewidth=1.5)
    for i, (bar, se) in enumerate(zip(bars, se_values)):
        height = bar.get_height()
        ax.errorbar(bar.get_x() + bar.get_width() / 2, height, yerr=se * 1.96,
                    fmt='none', color='black', capsize=5, capthick=2)
        ax.text(bar.get_x() + bar.get_width() / 2, height + se * 2,
                f'{variances[i]:.2f}\n+/-{se * 1.96:.2f}',
                ha='center', va='bottom', fontsize=fs + 2, fontweight='bold')

    ax.set_ylabel('Variance Component', fontsize=fs + 4, fontweight='bold')
    ax.tick_params(axis='both', labelsize=fs + 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    total_var = var_pe + var_ga + var_e
    h2 = var_ga / total_var
    repeatability = (var_ga + var_pe) / total_var
    # Annotation outside right side, no background, no frame, vertically centered
    ax.text(1.02, 0.5, f'h2 = {h2:.3f}\nRepeatability = {repeatability:.3f}',
            transform=ax.transAxes, fontsize=fs + 2, verticalalignment='center',
            horizontalalignment='left')

    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        plt.savefig(os.path.join(output_dir, f"variance_components_bar.{fmt}"),
                    dpi=300, bbox_inches='tight')
    plt.close()

    # ---- Figure 2: Breeding Value Distribution ----
    fig, ax = plt.subplots(figsize=(10, 6))
    n_hist, bins_hist, patches = ax.hist(ga_hat, bins=30, color=PALETTE_A[0], edgecolor='#333333',
                                         alpha=0.7, density=True)
    mu, sigma = np.mean(ga_hat), np.std(ga_hat)
    x = np.linspace(min(ga_hat), max(ga_hat), 100)
    ax.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', linewidth=2,
            label=f'Normal fit\nmu={mu:.2f}, sigma={sigma:.2f}')

    stat, p_value = stats.shapiro(ga_hat[:min(5000, len(ga_hat))])

    ax.set_xlabel('Breeding Value (GA)', fontsize=fs + 4, fontweight='bold')
    ax.set_ylabel('Density', fontsize=fs + 4, fontweight='bold')
    ax.tick_params(axis='both', labelsize=fs + 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend outside right side, no background, no frame, vertically centered, smaller font
    leg = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.55),
                    borderaxespad=0, frameon=False, fontsize=fs)
    # Shapiro annotation aligned with legend text, below legend, smaller font
    ax.text(1.08, 0.38, f'Shapiro-Wilk p={p_value:.4f}',
            transform=ax.transAxes, fontsize=fs, verticalalignment='center',
            horizontalalignment='left')

    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        plt.savefig(os.path.join(output_dir, f"breeding_value_distribution.{fmt}"),
                    dpi=300, bbox_inches='tight')
    plt.close()

    # ---- Figure 3: Fixed Effects Forest Plot ----
    fig, ax = plt.subplots(figsize=(12, 6))
    factor_names = ['mu (Intercept)', 'sex_M (Male vs Female)', 'season_Spring',
                    'season_Summer', 'season_Winter']
    n_factors = len(beta)
    y_pos = np.arange(n_factors)
    colors = [PALETTE_A[0], PALETTE_A[2], PALETTE_A[7], PALETTE_A[4], PALETTE_A[5]]

    for i, (name, est, se) in enumerate(zip(factor_names, beta, beta_se)):
        ci_lower = est - 1.96 * se
        ci_upper = est + 1.96 * se
        ax.errorbar(est, i, xerr=[[est - ci_lower], [ci_upper - est]],
                    fmt='o', color=colors[i], markersize=10, capsize=5, capthick=2,
                    ecolor=colors[i], elinewidth=2)
        p_val = 2 * (1 - stats.norm.cdf(abs(est) / se)) if se > 0 else 1.0
        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
        ax.text(ci_upper + 0.5, i, f'{est:.2f} +/- {se:.2f}{sig}',
                va='center', fontsize=fs + 2, fontweight='bold')

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(factor_names, fontsize=fs + 2, fontweight='bold')
    ax.set_xlabel('Effect Size (Coefficient)', fontsize=fs + 4, fontweight='bold')
    ax.tick_params(axis='x', labelsize=fs + 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Significance annotation below x-axis label, centered
    sig_text = '* p<0.05, ** p<0.01, *** p<0.001'
    fig.text(0.5, 0.01, sig_text, ha='center', fontsize=fs + 1, fontstyle='italic')

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    for fmt in ['pdf', 'png']:
        plt.savefig(os.path.join(output_dir, f"fixed_effects_forest.{fmt}"),
                    dpi=300, bbox_inches='tight')
    plt.close()

    print("  Visualizations saved.")


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
    bfile = args.bfile
    pheno_file = args.pheno
    snp_file = args.snp_file
    sample_file = args.sample_file
    output_dir = args.out

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("  Repeated Records Model Analysis")
    print("  Model: weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e")
    print("  Method: AI-REML (Woodbury identity)")
    print(f"  BFILE: {bfile}")
    print(f"  Phenotype: {pheno_file}")
    print(f"  SNP file: {snp_file}")
    print(f"  Sample file: {sample_file}")
    print(f"  Output: {output_dir}")
    print("=" * 70)

    print("\n[1] Reading SNP and sample lists...")
    snp_ids = pd.read_csv(snp_file, header=None)[0].tolist()
    sample_ids = pd.read_csv(sample_file, header=None)[0].astype(str).tolist()
    print(f"  SNPs: {len(snp_ids)}, Samples: {len(sample_ids)}")

    print("\n[2] Reading genotype data...")
    G_geno, bim, fam, p_geno = read_plink(bfile, sample_ids, snp_ids)
    print(f"  Genotype matrix: {G_geno.shape}")

    print("\n[3] Computing GRM (GA)...")
    G, Z_geno = calc_grm(G_geno, p_geno)
    print(f"  GA dim: {G.shape}, diagonal mean: {np.mean(np.diag(G)):.6f}")

    print("\n[4] Reading phenotype data (long format)...")
    pheno_long = pd.read_csv(pheno_file, sep='\t')
    pheno_long['ID'] = pheno_long['ID'].astype(int)
    train_ids = [int(x) for x in sample_ids]
    pheno_long = pheno_long[pheno_long['ID'].isin(train_ids)].reset_index(drop=True)
    unique_ids = sorted(pheno_long['ID'].unique())
    print(f"  Observation records: {len(pheno_long)}, Individuals: {len(unique_ids)}")

    print("\n[5] Building design matrices...")
    y, X, Z = build_design_matrices(pheno_long, unique_ids, args.pheno_col)
    print(f"  y dim: {y.shape}, X dim: {X.shape}, Z dim: {Z.shape}")

    print("\n[6] AI-REML variance component estimation...")
    print("  Alg.\tIter.\tLogL.\tV(ID)\tV(GA)\tV(e)")
    (var_pe, var_ga, var_e,
     se_var_pe, se_var_ga, se_var_e,
     h2_pe, h2_ga, h2_e,
     se_h2_pe, se_h2_ga, se_h2_e,
     p_value_pe, p_value_ga, p_value_e,
     beta, beta_se, pe_hat, ga_hat, residuals,
     log_lines) = ai_reml(y, X, Z, G)

    total_var = var_pe + var_ga + var_e
    h2 = var_ga / total_var
    repeatability = (var_ga + var_pe) / total_var

    print(f"\n  Variance components: V(ID)={var_pe:.6f}, V(GA)={var_ga:.6f}, V(e)={var_e:.6f}")
    print(f"  Heritability h2={h2:.6f}, Repeatability={repeatability:.6f}")

    print("\n[7] Estimating fixed effects and random effects...")
    level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
    for name, est, se in zip(level_names, beta, beta_se):
        print(f"  {name}: {est:.6g} +/- {se:.6g}")

    print("\n[8] ANOVA...")
    anova_results = compute_anova(y, X, var_e)
    print(f"  sex: F={anova_results['sex']['f']:.2f}, p={anova_results['sex']['p']:.2e}")
    print(f"  season: F={anova_results['season']['f']:.2f}, p={anova_results['season']['p']:.2e}")

    print("\n[9] Saving output files...")
    save_vars(os.path.join(output_dir, "repeated_model.vars"),
              var_pe, var_ga, var_e, se_var_pe, se_var_ga, se_var_e,
              h2_pe, h2_ga, h2_e, se_h2_pe, se_h2_ga, se_h2_e,
              p_value_pe, p_value_ga, p_value_e)

    save_beta(os.path.join(output_dir, "repeated_model.beta"), beta, beta_se)

    save_rand(os.path.join(output_dir, "repeated_model.rand"),
              pheno_long, unique_ids, pe_hat, ga_hat, residuals)

    save_anova(os.path.join(output_dir, "repeated_model.anova"), anova_results)

    save_log(os.path.join(output_dir, "repeated_model.log"),
             log_lines, var_pe, var_ga, var_e,
             se_var_pe, se_var_ga, se_var_e,
             h2_pe, h2_ga, h2_e,
             se_h2_pe, se_h2_ga, se_h2_e,
             p_value_pe, p_value_ga, p_value_e,
             beta, beta_se, len(snp_ids), len(unique_ids))

    # Bilingual summary files
    output_summary_en = os.path.join(output_dir, "repeated_records_summary.txt")
    with open(output_summary_en, 'w') as f:
        f.write("Repeated Records Model Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Model: weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e\n")
        f.write(f"Method: AI-REML (Woodbury identity)\n")
        f.write(f"SNPs: {len(snp_ids)}, Individuals: {len(unique_ids)}\n\n")
        f.write(f"Variance components:\n")
        f.write(f"  V(ID):  {var_pe:.6f} (SE={se_var_pe:.6f})\n")
        f.write(f"  V(GA):  {var_ga:.6f} (SE={se_var_ga:.6f})\n")
        f.write(f"  V(e):   {var_e:.6f} (SE={se_var_e:.6f})\n\n")
        f.write(f"Heritability: h2 = {h2:.6f}\n")
        f.write(f"Repeatability:    = {repeatability:.6f}\n\n")
        f.write("Fixed effects:\n")
        level_names_list = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
        for name, est, se in zip(level_names_list, beta, beta_se):
            f.write(f"  {name}: {est:.6g} +/- {se:.6g}\n")
    print(f"  Summary (EN) saved: {output_summary_en}")

    output_summary_zh = os.path.join(output_dir, "repeated_records_summary-zh.txt")
    with open(output_summary_zh, 'w') as f:
        f.write("重复记录模型分析摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"模型: weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e\n")
        f.write(f"方法: AI-REML (Woodbury恒等式)\n")
        f.write(f"SNP数: {len(snp_ids)}, 个体数: {len(unique_ids)}\n\n")
        f.write(f"方差组分:\n")
        f.write(f"  V(ID):  {var_pe:.6f} (SE={se_var_pe:.6f})\n")
        f.write(f"  V(GA):  {var_ga:.6f} (SE={se_var_ga:.6f})\n")
        f.write(f"  V(e):   {var_e:.6f} (SE={se_var_e:.6f})\n\n")
        f.write(f"遗传力: h2 = {h2:.6f}\n")
        f.write(f"重复力:    = {repeatability:.6f}\n\n")
        f.write("固定效应:\n")
        level_names_list = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
        for name, est, se in zip(level_names_list, beta, beta_se):
            f.write(f"  {name}: {est:.6g} +/- {se:.6g}\n")
    print(f"  摘要 (中文) 已保存: {output_summary_zh}")

    print("\n[10] Creating visualizations...")
    create_visualizations(var_pe, var_ga, var_e, se_var_pe, se_var_ga, se_var_e,
                          h2_pe, h2_ga, h2_e,
                          beta, beta_se, ga_hat, output_dir, args)

    print(f"\n  Results saved to: {output_dir}")
    print("\n" + "=" * 70)
    print("  Analysis complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. repeated_model.vars    - Variance component estimates")
    print(f"  2. repeated_model.beta    - Fixed effect estimates")
    print(f"  3. repeated_model.anova   - ANOVA table")
    print(f"  4. repeated_model.rand    - Random effect estimates")
    print(f"  5. repeated_model.log     - Analysis log")
    print(f"  6. variance_components_bar.pdf/png    - Variance components plot")
    print(f"  7. breeding_value_distribution.pdf/png - Breeding value distribution")
    print(f"  8. fixed_effects_forest.pdf/png       - Fixed effects forest plot")


if __name__ == "__main__":
    main()
