#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import eigh
from scipy.optimize import minimize
from pandas_plink import read_plink
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BFILE = "/public/share/likui/hanyu/testdata/In-silico-data/simulated_population"
PHENO_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/t10/2/simulated_phenotypes_multi_trait.txt"
SNP_FILE = "/public/share/likui/hanyu/testresult/HIBLUP/9/chr1_snps.txt"
SAMPLE_FILE = "/public/share/likui/hanyu/testresult/HIBLUP/9/train_samples.txt"
BASE_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t10/3"
T = 3
TRAIT_NAMES = ['T1', 'T2', 'T3']


def param_idx(i, j, t):
    return i * (2 * t - i - 1) // 2 + j


def read_genotype(bfile, snp_file, sample_file):
    snp_ids = pd.read_csv(snp_file, header=None)[0].tolist()
    sample_ids = pd.read_csv(sample_file, header=None)[0].astype(str).tolist()
    bim, fam, G_dask = read_plink(bfile, verbose=False)
    snp_idx = bim[bim['snp'].isin(set(snp_ids))].index.tolist()
    sample_iids = [int(x) for x in sample_ids]
    sample_mask = fam['iid'].astype(int).isin(sample_iids)
    sample_indices = fam[sample_mask].index.tolist()
    G = G_dask[snp_idx, :][:, sample_indices].compute()
    G_out = G.T.copy()
    for j in range(G_out.shape[1]):
        valid = ~np.isnan(G_out[:, j])
        if np.sum(valid) > 0 and np.sum(~valid) > 0:
            G_out[~valid, j] = np.nanmean(G_out[valid, j])
    return G_out, bim[bim['snp'].isin(set(snp_ids))].reset_index(drop=True), \
           fam[sample_mask].reset_index(drop=True)


def calc_grm(G):
    p = np.mean(G, axis=0) / 2
    Z = G - 2 * p
    sum_2pq = 2 * np.sum(p * (1 - p))
    GRM_vr = (Z @ Z.T) / sum_2pq
    diag_mean = np.mean(np.diag(GRM_vr))
    GRM = GRM_vr / diag_mean
    return GRM, p


def build_X_single(pheno):
    n = len(pheno)
    X = np.zeros((n, 5))
    X[:, 0] = 1
    X[:, 1] = (pheno['sex'] == 'M').astype(float).values
    season_dummies = pd.get_dummies(pheno['season'], prefix='season')
    for i, col in enumerate(['season_Spring', 'season_Summer', 'season_Winter']):
        if col in season_dummies.columns:
            X[:, 2 + i] = season_dummies[col].values
    return X


def chol_params_to_matrices(params, t):
    n_chol = t * (t + 1) // 2
    L_G = np.zeros((t, t))
    idx = 0
    for i in range(t):
        for j in range(i + 1):
            if i == j:
                L_G[i, j] = np.exp(params[idx])
            else:
                L_G[i, j] = params[idx]
            idx += 1
    G0 = L_G @ L_G.T

    L_R = np.zeros((t, t))
    idx = 0
    for i in range(t):
        for j in range(i + 1):
            if i == j:
                L_R[i, j] = np.exp(params[n_chol + idx])
            else:
                L_R[i, j] = params[n_chol + idx]
            idx += 1
    R0 = L_R @ L_R.T

    return G0, R0


def matrices_to_chol_params(G0, R0, t):
    L_G = np.linalg.cholesky(G0)
    params_G = []
    for i in range(t):
        for j in range(i + 1):
            if i == j:
                params_G.append(np.log(L_G[i, j]))
            else:
                params_G.append(L_G[i, j])

    L_R = np.linalg.cholesky(R0)
    params_R = []
    for i in range(t):
        for j in range(i + 1):
            if i == j:
                params_R.append(np.log(L_R[i, j]))
            else:
                params_R.append(L_R[i, j])

    return np.array(params_G + params_R)


def reml_logL_from_matrices(G0, R0, eigenvalues, y_tilde, X_tilde, n, t):
    nt = n * t
    V_blocks = eigenvalues[:, None, None] * G0[None, :, :] + R0[None, :, :]
    try:
        V_inv_blocks = np.linalg.inv(V_blocks)
        signs, logdets = np.linalg.slogdet(V_blocks)
    except np.linalg.LinAlgError:
        return -1e10, None, None, None, None

    logdet_V = np.sum(logdets)

    V_tilde_inv = np.zeros((nt, nt))
    for k in range(n):
        V_tilde_inv[k * t:(k + 1) * t, k * t:(k + 1) * t] = V_inv_blocks[k]

    XtVX = X_tilde.T @ V_tilde_inv @ X_tilde
    try:
        XtVX_inv = np.linalg.inv(XtVX)
    except np.linalg.LinAlgError:
        return -1e10, None, None, None, None

    beta = XtVX_inv @ (X_tilde.T @ V_tilde_inv @ y_tilde)
    r_tilde = y_tilde - X_tilde @ beta
    P_tilde_y = V_tilde_inv @ r_tilde

    sign_xvx, logdet_XtVX = np.linalg.slogdet(XtVX)
    logL = -0.5 * (logdet_V + logdet_XtVX + r_tilde @ P_tilde_y)

    return logL, V_inv_blocks, V_tilde_inv, XtVX_inv, beta


def reml_neg_logL(params, eigenvalues, y_tilde, X_tilde, n, t, p_total):
    G0, R0 = chol_params_to_matrices(params, t)
    logL, _, _, _, _ = reml_logL_from_matrices(
        G0, R0, eigenvalues, y_tilde, X_tilde, n, t)
    if logL == -1e10:
        return 1e10
    return -logL


def make_E(i, j, t):
    E = np.zeros((t, t))
    E[i, j] = 1.0
    if i != j:
        E[j, i] = 1.0
    return E


def compute_ai_matrix_and_score(G0, R0, eigenvalues, y_tilde, X_tilde,
                                 n, t, E_list, n_g0):
    nt = n * t
    n_params = len(E_list) * 2

    logL, V_inv_blocks, V_tilde_inv, XtVX_inv, beta = \
        reml_logL_from_matrices(G0, R0, eigenvalues, y_tilde, X_tilde, n, t)

    r_tilde = y_tilde - X_tilde @ beta
    P_tilde_y = V_tilde_inv @ r_tilde

    K_Py = np.zeros((n_params, nt))
    for a in range(n_params):
        for k in range(n):
            if a < n_g0:
                E_a = E_list[a]
                K_Py[a, k * t:(k + 1) * t] = eigenvalues[k] * (E_a @ P_tilde_y[k * t:(k + 1) * t])
            else:
                E_a = E_list[a - n_g0]
                K_Py[a, k * t:(k + 1) * t] = E_a @ P_tilde_y[k * t:(k + 1) * t]

    PK_Py = np.zeros((n_params, nt))
    for a in range(n_params):
        Vinv_KPy = V_tilde_inv @ K_Py[a]
        PK_Py[a] = Vinv_KPy - V_tilde_inv @ X_tilde @ (XtVX_inv @ (X_tilde.T @ Vinv_KPy))

    ai_mat = np.zeros((n_params, n_params))
    for a in range(n_params):
        for b in range(a, n_params):
            ai_mat[a, b] = 0.5 * PK_Py[a] @ K_Py[b]
            ai_mat[b, a] = ai_mat[a, b]

    score_vec = np.zeros(n_params)
    V_inv_X = V_tilde_inv @ X_tilde
    for a in range(n_params):
        tr_Vinv_Ka = 0.0
        for k in range(n):
            if a < n_g0:
                E_a = E_list[a]
                tr_Vinv_Ka += eigenvalues[k] * np.sum(V_inv_blocks[k] * E_a)
            else:
                E_a = E_list[a - n_g0]
                tr_Vinv_Ka += np.sum(V_inv_blocks[k] * E_a)

        Ka_VinvX = np.zeros_like(V_inv_X)
        for k in range(n):
            if a < n_g0:
                E_a = E_list[a]
                Ka_VinvX[k * t:(k + 1) * t, :] = eigenvalues[k] * (E_a @ V_inv_X[k * t:(k + 1) * t, :])
            else:
                E_a = E_list[a - n_g0]
                Ka_VinvX[k * t:(k + 1) * t, :] = E_a @ V_inv_X[k * t:(k + 1) * t, :]

        Vinv_Ka_VinvX = V_tilde_inv @ Ka_VinvX
        tr_XtVX_inv = np.trace(XtVX_inv @ (X_tilde.T @ Vinv_Ka_VinvX))
        tr_PKa = tr_Vinv_Ka - tr_XtVX_inv
        score_vec[a] = 0.5 * (P_tilde_y @ K_Py[a] - tr_PKa)

    return ai_mat, score_vec, logL, P_tilde_y, V_tilde_inv, XtVX_inv, beta


def compute_anova(R0_ii, n, p_each, y_i, X_i):
    beta_full = np.linalg.lstsq(X_i, y_i, rcond=None)[0]
    rss_full = np.sum((y_i - X_i @ beta_full) ** 2)
    X_mu_season = X_i[:, [0, 2, 3, 4]]
    rss_mu_season = np.sum((y_i - X_mu_season @ np.linalg.lstsq(X_mu_season, y_i, rcond=None)[0]) ** 2)
    X_mu_sex = X_i[:, :2]
    rss_mu_sex = np.sum((y_i - X_mu_sex @ np.linalg.lstsq(X_mu_sex, y_i, rcond=None)[0]) ** 2)
    ss_sex = rss_mu_season - rss_full
    ss_season = rss_mu_sex - rss_full
    df_sex = 1
    df_season = 3
    df_e = n - p_each
    ms_sex = ss_sex / df_sex
    ms_season = ss_season / df_season
    ms_e = R0_ii
    ss_e = R0_ii * df_e
    f_sex = ms_sex / ms_e
    f_season = ms_season / ms_e
    p_sex = stats.f.sf(f_sex, df_sex, df_e)
    p_season = stats.f.sf(f_season, df_season, df_e)
    return {
        'sex': {'df': df_sex, 'ss': ss_sex, 'ms': ms_sex, 'f': f_sex, 'p': p_sex},
        'season': {'df': df_season, 'ss': ss_season, 'ms': ms_season, 'f': f_season, 'p': p_season},
        'e': {'df': df_e, 'ss': ss_e, 'ms': ms_e}
    }


def format_sci(x):
    if x == 0:
        return "0"
    if abs(x) >= 0.001 and abs(x) < 1e6:
        return f"{x:.6g}"
    return f"{x:.5e}"


def save_beta(filepath, beta, beta_se, trait_idx, p_each):
    level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
    start = trait_idx * p_each
    end = (trait_idx + 1) * p_each
    with open(filepath, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        for name, est, se in zip(level_names, beta[start:end], beta_se[start:end]):
            f.write(f"{name}\t{est:.6g}\t{se:.6g}\n")


def save_rand(filepath, ids, ga_vals, resid_vals):
    with open(filepath, 'w') as f:
        f.write("ID\tGA\tresiduals\n")
        for id_, ga, res in zip(ids, ga_vals, resid_vals):
            f.write(f"{id_}\t{ga:.6g}\t{res:.6g}\n")


def save_anova(filepath, anova_results):
    with open(filepath, 'w') as f:
        f.write("Factors\tDf\tSumSq\tMeanSq\tF\tPr(>F)\n")
        r = anova_results['sex']
        p_str = "0" if r['p'] < 1e-15 else format_sci(r['p'])
        f.write(f"sex\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t{r['f']:.6g}\t{p_str}\n")
        r = anova_results['season']
        p_str = "0" if r['p'] < 1e-15 else format_sci(r['p'])
        f.write(f"season\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t{r['f']:.6g}\t{p_str}\n")
        r = anova_results['e']
        f.write(f"e\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t.\t.\n")


NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377',
                 '#BBCC33', '#CC3311', '#EE7733', '#999999']


def plot_variance_components(G0, R0, trait_names, output_dir):
    t = len(trait_names)
    var_ga = np.diag(G0)
    var_e = np.diag(R0)
    x = np.arange(t)
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    bars1 = ax.bar(x - width / 2, var_ga, width, label='V(GA)', color=NATURE_COLORS[0], edgecolor='black')
    bars2 = ax.bar(x + width / 2, var_e, width, label='V(e)', color=NATURE_COLORS[1], edgecolor='black')
    ax.set_xlabel('Trait')
    ax.set_ylabel('Variance')
    ax.set_title('Variance Components by Trait')
    ax.set_xticks(x)
    ax.set_xticklabels(trait_names)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "variance_components.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_heritability(h2_list, se_h2_list, trait_names, output_dir):
    t = len(trait_names)
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(t)
    ax.bar(x, h2_list, color=NATURE_COLORS[2], edgecolor='black', alpha=0.8)
    ax.errorbar(x, h2_list, yerr=se_h2_list, fmt='none', ecolor='black', capsize=5, linewidth=1.5)
    ax.set_xlabel('Trait')
    ax.set_ylabel('Heritability (h²)')
    ax.set_title('Heritability Estimates with Standard Errors')
    ax.set_xticks(x)
    ax.set_xticklabels(trait_names)
    ax.set_ylim(0, min(1.0, max(h2_list) + max(se_h2_list) * 2))
    ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "heritability.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_correlation_heatmap(matrix, trait_names, title, output_file):
    fig, ax = plt.subplots(figsize=(6, 5))
    mask = np.zeros_like(matrix, dtype=bool)
    np.fill_diagonal(mask, True)
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(matrix, mask=mask, annot=True, fmt='.3f', cmap=cmap,
                vmin=-1, vmax=1, center=0, square=True, linewidths=1,
                xticklabels=trait_names, yticklabels=trait_names, ax=ax,
                cbar_kws={"shrink": 0.8})
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_beta_forest(beta, beta_se, trait_names, p_each, output_dir):
    level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
    t = len(trait_names)
    n_effects = len(level_names)

    fig, axes = plt.subplots(1, t, figsize=(5 * t, 6), sharey=True)
    if t == 1:
        axes = [axes]

    for i in range(t):
        ax = axes[i]
        start = i * p_each
        end = (i + 1) * p_each
        b = beta[start:end]
        se = beta_se[start:end]
        y_pos = np.arange(n_effects)

        ax.errorbar(b, y_pos, xerr=1.96 * se, fmt='o', color=NATURE_COLORS[i],
                    ecolor=NATURE_COLORS[i], capsize=4, linewidth=1.5, markersize=6)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(level_names)
        ax.set_xlabel('Estimate')
        ax.set_title(trait_names[i])
        sns.despine(ax=ax)

    fig.suptitle('Fixed Effects Estimates (95% CI)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "beta_forest.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_random_effects_distribution(u_hat_list, e_hat_list, trait_names, output_dir):
    t = len(trait_names)
    fig, axes = plt.subplots(2, t, figsize=(5 * t, 8))

    for i in range(t):
        ax_ga = axes[0, i] if t > 1 else axes[0]
        ax_e = axes[1, i] if t > 1 else axes[1]

        ax_ga.hist(u_hat_list[i], bins=30, color=NATURE_COLORS[i], edgecolor='black', alpha=0.7, density=True)
        from scipy.stats import norm
        mu, std = norm.fit(u_hat_list[i])
        x = np.linspace(u_hat_list[i].min(), u_hat_list[i].max(), 100)
        ax_ga.plot(x, norm.pdf(x, mu, std), 'k-', linewidth=2)
        ax_ga.set_title(f'{trait_names[i]} GA')
        ax_ga.set_xlabel('GA effect')
        ax_ga.set_ylabel('Density')
        sns.despine(ax=ax_ga)

        ax_e.hist(e_hat_list[i], bins=30, color=NATURE_COLORS[i + 3] if i + 3 < len(NATURE_COLORS) else NATURE_COLORS[-1],
                  edgecolor='black', alpha=0.7, density=True)
        mu, std = norm.fit(e_hat_list[i])
        x = np.linspace(e_hat_list[i].min(), e_hat_list[i].max(), 100)
        ax_e.plot(x, norm.pdf(x, mu, std), 'k-', linewidth=2)
        ax_e.set_title(f'{trait_names[i]} Residuals')
        ax_e.set_xlabel('Residual')
        ax_e.set_ylabel('Density')
        sns.despine(ax=ax_e)

    fig.suptitle('Random Effects Distribution', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "random_effects_distribution.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_anova_bar(anova_results_list, trait_names, output_dir):
    t = len(trait_names)
    fig, axes = plt.subplots(1, t, figsize=(5 * t, 5))
    if t == 1:
        axes = [axes]

    for i in range(t):
        ax = axes[i]
        r = anova_results_list[i]
        factors = ['sex', 'season', 'e']
        ss = [r['sex']['ss'], r['season']['ss'], r['e']['ss']]
        colors = [NATURE_COLORS[0], NATURE_COLORS[2], NATURE_COLORS[-1]]
        ax.bar(factors, ss, color=colors, edgecolor='black', alpha=0.8)
        ax.set_ylabel('Sum of Squares')
        ax.set_title(trait_names[i])
        sns.despine(ax=ax)

    fig.suptitle('ANOVA Sum of Squares', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "anova_sumSq.pdf"), format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


def main():
    os.makedirs(BASE_DIR, exist_ok=True)

    print("=" * 70)
    print("  Multi-trait Model Analysis")
    print("  Model: T_i = 1 + sex(F) + season(F) + GA(R[G]) + e")
    print("  Method: REML (scipy.optimize + Cholesky parameterization)")
    print("=" * 70)

    print("\n[1] Reading SNP and sample lists...")
    snp_ids = pd.read_csv(SNP_FILE, header=None)[0].tolist()
    sample_ids = pd.read_csv(SAMPLE_FILE, header=None)[0].astype(str).tolist()
    print(f"  SNPs: {len(snp_ids)}, Samples: {len(sample_ids)}")

    print("\n[2] Reading genotype data (pandas_plink)...")
    G_geno, bim, fam = read_genotype(BFILE, SNP_FILE, SAMPLE_FILE)
    print(f"  Genotype matrix: {G_geno.shape}")

    print("\n[3] Computing GRM (VanRaden + Su normalization)...")
    GA, p_train = calc_grm(G_geno)
    print(f"  GA dim: {GA.shape}, diag mean: {np.mean(np.diag(GA)):.6f}")

    print("\n[4] Reading phenotype data...")
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno['ID'] = pheno['ID'].astype(int)
    train_ids = [int(x) for x in sample_ids]
    pheno = pheno[pheno['ID'].isin(train_ids)].copy()
    fam_iids = fam['iid'].values.astype(int)
    pheno = pheno.set_index('ID').loc[fam_iids].reset_index()
    print(f"  Individuals: {len(pheno)}")

    print("\n[5] Building design matrices...")
    y_list = []
    X_list = []
    for trait in TRAIT_NAMES:
        y_list.append(pheno[trait].values.astype(float))
        X_list.append(build_X_single(pheno))
    n = len(pheno)
    p_each = X_list[0].shape[1]
    p_total = p_each * T
    print(f"  Each trait: y({n},), X({n},{p_each})")

    print("\n[6] Eigenvalue decomposition and data transformation...")
    eigenvalues, U = eigh(GA)

    Uty_list = [U.T @ y_list[i] for i in range(T)]
    UtX_list = [U.T @ X_list[i] for i in range(T)]

    y_tilde = np.vstack(Uty_list).T.flatten()

    X_tilde = np.zeros((n * T, p_total))
    for i in range(T):
        rows = np.arange(n) * T + i
        X_tilde[rows, i * p_each:(i + 1) * p_each] = UtX_list[i]

    print(f"  y_tilde: {y_tilde.shape}, X_tilde: {X_tilde.shape}")

    print("\n[7] Verifying LogL at HIBLUP's parameters...")
    G0_hiblup = np.array([
        [29.4388, 42.8831, 4.8201],
        [42.8831, 63.7777, 7.4902],
        [4.8201, 7.4902, 0.9579]
    ])
    R0_hiblup = np.array([
        [85.6507, 73.9569, 16.5609],
        [73.9569, 104.218, 15.6874],
        [16.5609, 15.6874, 5.5582]
    ])
    logL_hiblup, _, _, _, _ = reml_logL_from_matrices(
        G0_hiblup, R0_hiblup, eigenvalues, y_tilde, X_tilde, n, T)
    print(f"  LogL at HIBLUP params: {logL_hiblup:.2f}")
    print(f"  HIBLUP reported LogL: ~-6264.39 (iter 4) to -6266.47 (iter 20)")

    print("\n[8] REML optimization (scipy.optimize.minimize + Cholesky)...")
    var_y = [np.var(y_list[i]) for i in range(T)]
    G0_init = np.zeros((T, T))
    R0_init = np.zeros((T, T))
    for i in range(T):
        G0_init[i, i] = var_y[i] * 0.3
        R0_init[i, i] = var_y[i] * 0.7
    for i in range(T):
        for j in range(i + 1, T):
            r_pheno = np.corrcoef(y_list[i], y_list[j])[0, 1]
            G0_init[i, j] = G0_init[j, i] = r_pheno * np.sqrt(G0_init[i, i] * G0_init[j, j])
            R0_init[i, j] = R0_init[j, i] = r_pheno * np.sqrt(R0_init[i, i] * R0_init[j, j]) * 0.5

    params_init = matrices_to_chol_params(G0_init, R0_init, T)
    neg_logL_init = reml_neg_logL(params_init, eigenvalues, y_tilde, X_tilde, n, T, p_total)
    print(f"  Initial LogL: {-neg_logL_init:.2f}")

    iter_count = [0]
    def callback_func(xk):
        iter_count[0] += 1
        if iter_count[0] % 10 == 0 or iter_count[0] <= 5:
            neg_ll = reml_neg_logL(xk, eigenvalues, y_tilde, X_tilde, n, T, p_total)
            G0_c, R0_c = chol_params_to_matrices(xk, T)
            line = f"  [L-BFGS-B] {iter_count[0]}\t{-neg_ll:.2f}\t"
            for i in range(T):
                line += f"{G0_c[i,i]:.5f}\t{R0_c[i,i]:.5f}\t"
            print(line)

    result = minimize(
        reml_neg_logL,
        params_init,
        args=(eigenvalues, y_tilde, X_tilde, n, T, p_total),
        method='L-BFGS-B',
        callback=callback_func,
        options={'maxiter': 1000, 'ftol': 1e-12, 'gtol': 1e-8, 'disp': False}
    )

    G0_opt, R0_opt = chol_params_to_matrices(result.x, T)
    logL_opt = -result.fun

    print(f"\n  L-BFGS-B result: success={result.success}, msg={result.message}")
    print(f"  LogL: {logL_opt:.2f}")
    for i in range(T):
        total = G0_opt[i, i] + R0_opt[i, i]
        h2 = G0_opt[i, i] / total if total > 0 else 0
        print(f"  T{i+1}: V(GA)={G0_opt[i,i]:.4f}, V(e)={R0_opt[i,i]:.4f}, h2={h2:.4f}")

    if logL_hiblup > logL_opt:
        print(f"\n  WARNING: HIBLUP params give higher LogL ({logL_hiblup:.2f}) than optimizer ({logL_opt:.2f})")
        print("  Trying Nelder-Mead from HIBLUP params...")
        params_hiblup = matrices_to_chol_params(G0_hiblup, R0_hiblup, T)
        result2 = minimize(
            reml_neg_logL,
            params_hiblup,
            args=(eigenvalues, y_tilde, X_tilde, n, T, p_total),
            method='Nelder-Mead',
            options={'maxiter': 50000, 'xatol': 1e-10, 'fatol': 1e-10, 'disp': False}
        )
        G0_nm, R0_nm = chol_params_to_matrices(result2.x, T)
        logL_nm = -result2.fun
        print(f"  Nelder-Mead LogL: {logL_nm:.2f}")
        if logL_nm > logL_opt:
            print("  Nelder-Mead found better solution, using it.")
            G0_opt, R0_opt = G0_nm, R0_nm
            logL_opt = logL_nm

    print("\n[9] Computing AI matrix for standard errors...")
    n_g0 = T * (T + 1) // 2
    E_list = []
    for i in range(T):
        for j in range(i, T):
            E_list.append(make_E(i, j, T))

    ai_mat, score_vec, logL_final, P_tilde_y, V_tilde_inv, XtVX_inv, beta = \
        compute_ai_matrix_and_score(G0_opt, R0_opt, eigenvalues, y_tilde,
                                     X_tilde, n, T, E_list, n_g0)

    try:
        cov_mat = np.linalg.inv(ai_mat)
    except np.linalg.LinAlgError:
        cov_mat = np.linalg.pinv(ai_mat)

    se_theta = np.sqrt(np.maximum(0, np.diag(cov_mat)))

    print(f"  Final LogL: {logL_final:.2f}")
    print(f"  Score norm: {np.linalg.norm(score_vec):.6f}")

    print("\n  Parameter SE verification:")
    for i in range(T):
        ga_idx = param_idx(i, i, T)
        r_idx = n_g0 + param_idx(i, i, T)
        print(f"  T{i+1}_GA: idx={ga_idx}, SE={se_theta[ga_idx]:.4f}")
        print(f"  T{i+1}_e:  idx={r_idx}, SE={se_theta[r_idx]:.4f}")

    print("\n[10] Computing heritability SE...")
    se_h2_list = []
    for i in range(T):
        total = G0_opt[i, i] + R0_opt[i, i]
        h2 = G0_opt[i, i] / total if total > 0 else 0
        grad_h2 = np.zeros(len(se_theta))
        a_idx = param_idx(i, i, T)
        b_idx = n_g0 + param_idx(i, i, T)
        if total > 0:
            grad_h2[a_idx] = R0_opt[i, i] / total ** 2
            grad_h2[b_idx] = -G0_opt[i, i] / total ** 2
        se_h2 = np.sqrt(max(0, grad_h2 @ cov_mat @ grad_h2))
        se_h2_list.append(se_h2)

    print(f"\n  Variance component estimates:")
    for i in range(T):
        total = G0_opt[i, i] + R0_opt[i, i]
        h2 = G0_opt[i, i] / total if total > 0 else 0
        print(f"  T{i+1}: V(GA)={G0_opt[i,i]:.4f}, V(e)={R0_opt[i,i]:.4f}, h2={h2:.4f} +/- {se_h2_list[i]:.4f}")

    print(f"\n  Genetic correlations:")
    for i in range(T):
        for j in range(i + 1, T):
            r_g = G0_opt[i, j] / np.sqrt(G0_opt[i, i] * G0_opt[j, j]) \
                if G0_opt[i, i] > 0 and G0_opt[j, j] > 0 else 0
            print(f"  T{i+1}-T{j+1}: r_g={r_g:.4f}")

    print("\n[11] Estimating fixed effects and random effects...")
    beta_se = np.sqrt(np.abs(np.diag(XtVX_inv)))

    nt = n * T
    u_tilde = np.zeros(nt)
    for k in range(n):
        u_tilde[k * T:(k + 1) * T] = eigenvalues[k] * G0_opt @ P_tilde_y[k * T:(k + 1) * T]

    u_hat_list = []
    for i in range(T):
        u_hat_i = u_tilde[i::T]
        u_hat_list.append(U @ u_hat_i)

    e_hat_list = []
    for i in range(T):
        e_hat_i = y_list[i] - X_list[i] @ beta[i * p_each:(i + 1) * p_each] - u_hat_list[i]
        e_hat_list.append(e_hat_i)

    for i in range(T):
        level_names = ['mu', 'sex_M', 'season_Spring', 'season_Summer', 'season_Winter']
        start = i * p_each
        end = (i + 1) * p_each
        print(f"  T{i+1}:")
        for name, est, se in zip(level_names, beta[start:end], beta_se[start:end]):
            print(f"    {name}: {est:.6g} +/- {se:.6g}")

    print("\n[12] ANOVA...")
    anova_results_list = []
    for i in range(T):
        anova_r = compute_anova(
            R0_opt[i, i], n, p_each,
            y_list[i], X_list[i])
        anova_results_list.append(anova_r)
        print(f"  T{i+1}: sex SS={anova_r['sex']['ss']:.2f}, season SS={anova_r['season']['ss']:.2f}")

    print("\n[13] Saving output files...")

    with open(os.path.join(BASE_DIR, "multi_trait.vars"), 'w') as f:
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        for i in range(T):
            var_ga = G0_opt[i, i]
            ga_idx = param_idx(i, i, T)
            se_ga = se_theta[ga_idx]
            total = G0_opt[i, i] + R0_opt[i, i]
            h2 = var_ga / total if total > 0 else 0
            wald_h2 = (h2 / se_h2_list[i]) ** 2 if se_h2_list[i] > 0 else 0
            p_val = stats.chi2.sf(wald_h2, df=1)
            f.write(f"T{i+1}_GA\t{var_ga:.6g}\t{se_ga:.6g}\t{h2:.6g}\t{se_h2_list[i]:.6g}\t{format_sci(p_val)}\n")

            r_idx = n_g0 + param_idx(i, i, T)
            se_e = se_theta[r_idx]
            h2_e = R0_opt[i, i] / total if total > 0 else 0
            wald_h2_e = (h2_e / se_h2_list[i]) ** 2 if se_h2_list[i] > 0 else 0
            p_val_e = stats.chi2.sf(wald_h2_e, df=1)
            f.write(f"T{i+1}_e\t{R0_opt[i,i]:.6g}\t{se_e:.6g}\t{h2_e:.6g}\t{se_h2_list[i]:.6g}\t{format_sci(p_val_e)}\n")

    with open(os.path.join(BASE_DIR, "multi_trait.covars"), 'w') as f:
        f.write("Item\tCOVar\tCOVar_SE\tr\tr_SE\tr_Pr(Chisq)\n")
        for i in range(T):
            for j in range(i + 1, T):
                cov_idx = param_idx(i, j, T)
                cov_ga = G0_opt[i, j]
                se_cov_ga = se_theta[cov_idx]
                r_g = cov_ga / np.sqrt(G0_opt[i, i] * G0_opt[j, j]) \
                    if G0_opt[i, i] > 0 and G0_opt[j, j] > 0 else 0
                a_idx = param_idx(i, i, T)
                b_idx = param_idx(j, j, T)
                grad_r = np.zeros(len(se_theta))
                if G0_opt[i, i] > 0 and G0_opt[j, j] > 0:
                    denom = np.sqrt(G0_opt[i, i] * G0_opt[j, j])
                    grad_r[cov_idx] = 1.0 / denom
                    grad_r[a_idx] = -cov_ga / (2 * G0_opt[i, i] * denom)
                    grad_r[b_idx] = -cov_ga / (2 * G0_opt[j, j] * denom)
                se_r = np.sqrt(max(0, grad_r @ cov_mat @ grad_r))
                wald_r = (r_g / se_r) ** 2 if se_r > 0 else 0
                p_val = stats.chi2.sf(wald_r, df=1)
                f.write(f"T{i+1}:T{j+1}_GA\t{cov_ga:.6g}\t{se_cov_ga:.6g}\t{r_g:.6g}\t{se_r:.6g}\t{format_sci(p_val)}\n")

        for i in range(T):
            for j in range(i + 1, T):
                cov_idx = n_g0 + param_idx(i, j, T)
                cov_e = R0_opt[i, j]
                se_cov_e = se_theta[cov_idx]
                r_e = cov_e / np.sqrt(R0_opt[i, i] * R0_opt[j, j]) \
                    if R0_opt[i, i] > 0 and R0_opt[j, j] > 0 else 0
                a_idx_r = n_g0 + param_idx(i, i, T)
                b_idx_r = n_g0 + param_idx(j, j, T)
                grad_r = np.zeros(len(se_theta))
                if R0_opt[i, i] > 0 and R0_opt[j, j] > 0:
                    denom = np.sqrt(R0_opt[i, i] * R0_opt[j, j])
                    grad_r[cov_idx] = 1.0 / denom
                    grad_r[a_idx_r] = -cov_e / (2 * R0_opt[i, i] * denom)
                    grad_r[b_idx_r] = -cov_e / (2 * R0_opt[j, j] * denom)
                se_r = np.sqrt(max(0, grad_r @ cov_mat @ grad_r))
                wald_r = (r_e / se_r) ** 2 if se_r > 0 else 0
                p_val = stats.chi2.sf(wald_r, df=1)
                f.write(f"T{i+1}:T{j+1}_e\t{cov_e:.6g}\t{se_cov_e:.6g}\t{r_e:.6g}\t{se_r:.6g}\t{format_sci(p_val)}\n")

    ids = pheno['ID'].values
    for i in range(T):
        save_beta(os.path.join(BASE_DIR, f"multi_trait.T{i+1}.beta"),
                  beta, beta_se, i, p_each)
        save_rand(os.path.join(BASE_DIR, f"multi_trait.T{i+1}.rand"),
                  ids, u_hat_list[i], e_hat_list[i])
        save_anova(os.path.join(BASE_DIR, f"multi_trait.T{i+1}.anova"),
                   anova_results_list[i])

    with open(os.path.join(BASE_DIR, "multi_trait.log"), 'w') as f:
        f.write("Multi-trait Model Analysis Log\n")
        f.write("=" * 50 + "\n")
        f.write(f"Model: T_i = 1 + sex(F) + season(F) + GA(R[G]) + e\n")
        f.write(f"Method: REML (scipy.optimize L-BFGS-B + Cholesky parameterization)\n")
        f.write(f"GRM: VanRaden + Su normalization\n")
        f.write(f"SNPs: {len(snp_ids)}, Samples: {n}, Traits: {T}\n")
        f.write(f"\nFinal LogL: {logL_final:.6f}\n")
        f.write(f"\nVariance components:\n")
        for i in range(T):
            total = G0_opt[i, i] + R0_opt[i, i]
            h2 = G0_opt[i, i] / total if total > 0 else 0
            f.write(f"  T{i+1}_GA: {G0_opt[i,i]:.6f}\n")
            f.write(f"  T{i+1}_e:  {R0_opt[i,i]:.6f}\n")
            f.write(f"  T{i+1}_h2: {h2:.6f}\n")
        f.write(f"\nCovariance components (genetic):\n")
        for i in range(T):
            for j in range(i + 1, T):
                f.write(f"  T{i+1}:T{j+1}_GA: {G0_opt[i,j]:.6f}\n")
        f.write(f"\nCovariance components (residual):\n")
        for i in range(T):
            for j in range(i + 1, T):
                f.write(f"  T{i+1}:T{j+1}_e: {R0_opt[i,j]:.6f}\n")

    print(f"\n  Results saved to: {BASE_DIR}")

    print("\n[14] Generating visualization plots...")
    h2_list = []
    for i in range(T):
        total = G0_opt[i, i] + R0_opt[i, i]
        h2_list.append(G0_opt[i, i] / total if total > 0 else 0)

    plot_variance_components(G0_opt, R0_opt, TRAIT_NAMES, BASE_DIR)
    print("  Saved: variance_components.pdf")

    plot_heritability(h2_list, se_h2_list, TRAIT_NAMES, BASE_DIR)
    print("  Saved: heritability.pdf")

    r_g_matrix = np.zeros((T, T))
    for i in range(T):
        r_g_matrix[i, i] = 1.0
        for j in range(i + 1, T):
            r_g_matrix[i, j] = G0_opt[i, j] / np.sqrt(G0_opt[i, i] * G0_opt[j, j]) \
                if G0_opt[i, i] > 0 and G0_opt[j, j] > 0 else 0
            r_g_matrix[j, i] = r_g_matrix[i, j]
    plot_correlation_heatmap(r_g_matrix, TRAIT_NAMES, 'Genetic Correlation',
                             os.path.join(BASE_DIR, "genetic_correlation.pdf"))
    print("  Saved: genetic_correlation.pdf")

    r_e_matrix = np.zeros((T, T))
    for i in range(T):
        r_e_matrix[i, i] = 1.0
        for j in range(i + 1, T):
            r_e_matrix[i, j] = R0_opt[i, j] / np.sqrt(R0_opt[i, i] * R0_opt[j, j]) \
                if R0_opt[i, i] > 0 and R0_opt[j, j] > 0 else 0
            r_e_matrix[j, i] = r_e_matrix[i, j]
    plot_correlation_heatmap(r_e_matrix, TRAIT_NAMES, 'Residual Correlation',
                             os.path.join(BASE_DIR, "residual_correlation.pdf"))
    print("  Saved: residual_correlation.pdf")

    plot_beta_forest(beta, beta_se, TRAIT_NAMES, p_each, BASE_DIR)
    print("  Saved: beta_forest.pdf")

    plot_random_effects_distribution(u_hat_list, e_hat_list, TRAIT_NAMES, BASE_DIR)
    print("  Saved: random_effects_distribution.pdf")

    plot_anova_bar(anova_results_list, TRAIT_NAMES, BASE_DIR)
    print("  Saved: anova_sumSq.pdf")

    print("\n" + "=" * 70)
    print("  Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
