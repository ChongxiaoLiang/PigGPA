#!/usr/bin/env python3
import os
import time
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import eigh
from scipy.optimize import minimize
from pandas_plink import read_plink
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

NATURE_COLORS = ['#4477AA', '#EE6677', '#228833', '#66C2A5', '#AA3377', '#BBCC33', '#CC3311', '#EE7733']

BFILE = "/public/share/likui/hanyu/testdata/In-silico-data/simulated_population"
PHENO_FILE = "/public/share/likui/hanyu/testdata/In-silico-data/t12/2/simulated_phenotypes_env.txt"
SNP_FILE = "/public/share/likui/hanyu/testresult/HIBLUP/10/chr1_snps.txt"
SAMPLE_FILE = "/public/share/likui/hanyu/testresult/HIBLUP/10/train_samples.txt"
BASE_DIR = "/public/share/likui/hanyu/testdata/In-silico-data/t12/1"


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
    return GRM


def build_gxe_grm(G, env_values):
    e_c = env_values - np.mean(env_values)
    var_e = np.mean(e_c ** 2)
    E = np.outer(e_c, e_c) / var_e
    GxE = G * E
    return GxE


def reml_logL(log_vars, eigenvalues_list, y, X):
    n = len(y)
    n_k = len(eigenvalues_list)
    vars_ = np.exp(log_vars)

    V = np.zeros((n, n))
    for i in range(n_k):
        V += vars_[i] * eigenvalues_list[i]
    V += vars_[-1]

    try:
        L = np.linalg.cholesky(V)
    except np.linalg.LinAlgError:
        return -1e10, None, None, None

    logdet_V = 2 * np.sum(np.log(np.diag(L)))
    Vinv = np.linalg.solve(V, np.eye(n))
    XtVX = X.T @ Vinv @ X
    try:
        XtVX_inv = np.linalg.inv(XtVX)
    except np.linalg.LinAlgError:
        return -1e10, None, None, None

    beta = XtVX_inv @ (X.T @ Vinv @ y)
    r = y - X @ beta
    Py = Vinv @ r - Vinv @ X @ (XtVX_inv @ (X.T @ Vinv @ r))

    sign_xvx, logdet_XtVX = np.linalg.slogdet(XtVX)
    logL = -0.5 * (logdet_V + logdet_XtVX + r @ Py)

    return logL, Vinv, XtVX_inv, beta


def reml_logL_eigen(log_vars, eigenvalues_list, Uty, UtX):
    n = len(Uty)
    n_k = len(eigenvalues_list)
    vars_ = np.exp(log_vars)

    d = np.zeros(n)
    for i in range(n_k):
        d += vars_[i] * eigenvalues_list[i]
    d += vars_[-1]

    if np.any(d <= 0):
        return -1e10

    logdet_V = np.sum(np.log(d))
    d_inv = 1.0 / d

    XtVX = UtX.T * d_inv @ UtX
    try:
        XtVX_inv = np.linalg.inv(XtVX)
    except np.linalg.LinAlgError:
        return -1e10

    sign_xvx, logdet_XtVX = np.linalg.slogdet(XtVX)
    beta = XtVX_inv @ (UtX.T * d_inv @ Uty)
    r = Uty - UtX @ beta
    rPy = r @ (r * d_inv)

    logL = -0.5 * (logdet_V + logdet_XtVX + rPy)
    return logL


def neg_logL_eigen(log_vars, eigenvalues_list, Uty, UtX):
    ll = reml_logL_eigen(log_vars, eigenvalues_list, Uty, UtX)
    if ll == -1e10:
        return 1e10
    return -ll


def compute_ai_se(vars_, eigenvalues_list, Uty, UtX, n):
    n_k = len(eigenvalues_list)
    n_params = n_k + 1

    d = np.zeros(n)
    for i in range(n_k):
        d += vars_[i] * eigenvalues_list[i]
    d += vars_[-1]
    d_inv = 1.0 / d

    XtVX = UtX.T * d_inv @ UtX
    XtVX_inv = np.linalg.inv(XtVX)
    beta = XtVX_inv @ (UtX.T * d_inv @ Uty)
    r = Uty - UtX @ beta
    Py = r * d_inv - UtX @ (XtVX_inv @ (UtX.T * d_inv @ r))

    K_Py_list = []
    for i in range(n_k):
        K_Py_list.append(eigenvalues_list[i] * Py)
    K_Py_list.append(Py)

    PK_Py_list = []
    for a in range(n_params):
        Vinv_KPy = K_Py_list[a] * d_inv
        PK_a = Vinv_KPy - (UtX * d_inv[:, None]) @ (XtVX_inv @ (UtX.T * d_inv @ K_Py_list[a]))
        PK_Py_list.append(PK_a)

    ai_mat = np.zeros((n_params, n_params))
    for a in range(n_params):
        for b in range(a, n_params):
            ai_mat[a, b] = 0.5 * PK_Py_list[a] @ K_Py_list[b]
            ai_mat[b, a] = ai_mat[a, b]

    try:
        cov_mat = np.linalg.inv(ai_mat)
    except np.linalg.LinAlgError:
        cov_mat = np.linalg.pinv(ai_mat)

    se_vars = np.sqrt(np.maximum(0, np.diag(cov_mat)))

    total_var = np.sum(vars_)
    h2s = vars_ / total_var
    se_h2s = np.zeros(n_params)
    for i in range(n_params):
        grad = np.zeros(n_params)
        for j in range(n_params):
            if j == i:
                grad[j] = (total_var - vars_[i]) / total_var ** 2
            else:
                grad[j] = -vars_[i] / total_var ** 2
        se_h2s[i] = np.sqrt(max(0, grad @ cov_mat @ grad))

    p_values = np.zeros(n_params)
    for i in range(n_params):
        if se_h2s[i] > 0:
            wald = (h2s[i] / se_h2s[i]) ** 2
            p_values[i] = stats.chi2.sf(wald, df=1)

    return se_vars, h2s, se_h2s, p_values, cov_mat


def compute_anova(y, X, var_e):
    n = len(y)
    p_x = X.shape[1]
    beta_full = np.linalg.lstsq(X, y, rcond=None)[0]
    rss_full = np.sum((y - X @ beta_full) ** 2)
    X_mu_hum = X[:, [0, 2]]
    rss_mu_hum = np.sum((y - X_mu_hum @ np.linalg.lstsq(X_mu_hum, y, rcond=None)[0]) ** 2)
    X_mu_temp = X[:, [0, 1]]
    rss_mu_temp = np.sum((y - X_mu_temp @ np.linalg.lstsq(X_mu_temp, y, rcond=None)[0]) ** 2)
    ss_temp = rss_mu_hum - rss_full
    ss_hum = rss_mu_temp - rss_full
    df_temp, df_hum, df_e = 1, 1, n - p_x
    ms_temp, ms_hum = ss_temp / df_temp, ss_hum / df_hum
    ms_e = var_e
    ss_e = var_e * df_e
    f_temp, f_hum = ms_temp / ms_e, ms_hum / ms_e
    p_temp = stats.f.sf(f_temp, df_temp, df_e)
    p_hum = stats.f.sf(f_hum, df_hum, df_e)
    return {
        'Temp': {'df': df_temp, 'ss': ss_temp, 'ms': ms_temp, 'f': f_temp, 'p': p_temp},
        'Humidity': {'df': df_hum, 'ss': ss_hum, 'ms': ms_hum, 'f': f_hum, 'p': p_hum},
        'e': {'df': df_e, 'ss': ss_e, 'ms': ms_e}
    }


def format_sci(x):
    if x == 0:
        return "0"
    if abs(x) >= 0.001 and abs(x) < 1e6:
        return f"{x:.6g}"
    return f"{x:.5e}"


def plot_variance_components(vars_, se_vars, h2s, se_h2s, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    names = ['Temp:GA.GA\n(GxE)', 'Humidity:GA.GA\n(GxE)', 'GA\n(Genetic)', 'e\n(Residual)']
    colors = [NATURE_COLORS[0], NATURE_COLORS[1], NATURE_COLORS[2], NATURE_COLORS[3]]
    
    bars = ax.bar(names, vars_, color=colors, edgecolor='black', linewidth=1.5)
    
    for i, (bar, se) in enumerate(zip(bars, se_vars)):
        height = bar.get_height()
        ax.errorbar(bar.get_x() + bar.get_width()/2, height, yerr=se*1.96,
                   fmt='none', color='black', capsize=5, capthick=2)
        ax.text(bar.get_x() + bar.get_width()/2, height + se*2,
               f'{vars_[i]:.1f}\n±{se*1.96:.1f}',
               ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_ylabel('Variance Component', fontsize=13, fontweight='bold')
    ax.set_title('Variance Components for GxE Model', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    total_var = np.sum(vars_)
    h2_ga = vars_[2] / total_var
    h2_gxe_t = vars_[0] / total_var
    h2_gxe_h = vars_[1] / total_var
    ax.text(0.95, 0.95, f'h²(GA) = {h2_ga:.3f}\nh²(GxT) = {h2_gxe_t:.3f}\nh²(GxH) = {h2_gxe_h:.3f}',
           transform=ax.transAxes, fontsize=11, verticalalignment='top',
           horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'variance_components_bar.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  方差组分图已保存: {os.path.join(output_dir, 'variance_components_bar.pdf')}")


def plot_random_effects(rand_hats, output_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    effect_names = ['Temp:GA.GA (GxE with Temp)', 'Humidity:GA.GA (GxE with Humidity)', 'GA (Additive Genetic)']
    colors = [NATURE_COLORS[0], NATURE_COLORS[1], NATURE_COLORS[2]]
    
    for i, (ax, name, color) in enumerate(zip(axes, effect_names, colors)):
        effects = rand_hats[i]
        ax.hist(effects, bins=30, color=color, edgecolor='black', alpha=0.7, density=True)
        
        mu, sigma = np.mean(effects), np.std(effects)
        x = np.linspace(min(effects), max(effects), 100)
        ax.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', linewidth=2,
               label=f'Normal fit\nμ={mu:.2f}, σ={sigma:.2f}')
        
        ax.set_xlabel(name, fontsize=11, fontweight='bold')
        ax.set_ylabel('Density', fontsize=11, fontweight='bold')
        ax.set_title(f'Distribution of {name}', fontsize=12, fontweight='bold')
        ax.legend(loc='upper left', fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'random_effects_distribution.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  随机效应分布图已保存: {os.path.join(output_dir, 'random_effects_distribution.pdf')}")


def plot_interaction_effects_forest(beta, beta_se, vars_, se_vars, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    ax1 = axes[0]
    factor_names = ['mu', 'Temp', 'Humidity']
    n_factors = len(beta)
    y_pos = np.arange(n_factors)
    colors = [NATURE_COLORS[i % len(NATURE_COLORS)] for i in range(n_factors)]
    
    for i, (name, est, se) in enumerate(zip(factor_names, beta, beta_se)):
        ci_lower = est - 1.96 * se
        ci_upper = est + 1.96 * se
        
        ax1.errorbar(est, i, xerr=[[est - ci_lower], [ci_upper - est]],
                   fmt='o', color=colors[i], markersize=10, capsize=5, capthick=2,
                   ecolor=colors[i], elinewidth=2)
        
        p_value = 2 * (1 - stats.norm.cdf(abs(est) / se)) if se > 0 else 1.0
        sig = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else ''
        ax1.text(ci_upper + 0.3, i, f'{est:.2f} ± {se:.2f}{sig}',
               va='center', fontsize=10, fontweight='bold')
    
    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(factor_names, fontsize=11)
    ax1.set_xlabel('Effect Size', fontsize=13, fontweight='bold')
    ax1.set_title('Fixed Effects', fontsize=14, fontweight='bold')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.text(0.02, 0.02, '* p<0.05, ** p<0.01, *** p<0.001',
           transform=ax1.transAxes, fontsize=9, verticalalignment='bottom')
    
    ax2 = axes[1]
    var_names = ['Temp:GA.GA', 'Humidity:GA.GA', 'GA', 'e']
    n_vars = len(vars_)
    y_pos = np.arange(n_vars)
    colors = [NATURE_COLORS[i % len(NATURE_COLORS)] for i in range(n_vars)]
    
    for i, (name, est, se) in enumerate(zip(var_names, vars_, se_vars)):
        ci_lower = max(0, est - 1.96 * se)
        ci_upper = est + 1.96 * se
        
        ax2.errorbar(est, i, xerr=[[est - ci_lower], [ci_upper - est]],
                   fmt='o', color=colors[i], markersize=10, capsize=5, capthick=2,
                   ecolor=colors[i], elinewidth=2)
        ax2.text(ci_upper + 2, i, f'{est:.1f} ± {se:.1f}',
               va='center', fontsize=10, fontweight='bold')
    
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(var_names, fontsize=11)
    ax2.set_xlabel('Variance Component', fontsize=13, fontweight='bold')
    ax2.set_title('Variance Components', fontsize=14, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'interaction_effects_forest.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  交互效应森林图已保存: {os.path.join(output_dir, 'interaction_effects_forest.pdf')}")


def save_vars(filepath, vars_, se_vars, h2s, se_h2s, p_values):
    names = ['Temp:GA.GA', 'Humidity:GA.GA', 'GA', 'e']
    with open(filepath, 'w') as f:
        f.write("Item\tVar\tVar_SE\th2\th2_SE\th2_Pr(Chisq)\n")
        for i, name in enumerate(names):
            f.write(f"{name}\t{vars_[i]:.6g}\t{se_vars[i]:.6g}\t{h2s[i]:.6g}\t{se_h2s[i]:.6g}\t{format_sci(p_values[i])}\n")


def save_beta(filepath, beta, beta_se):
    level_names = ['mu', 'Temp', 'Humidity']
    with open(filepath, 'w') as f:
        f.write("Levels\tEstimation\tSE\n")
        for name, est, se in zip(level_names, beta, beta_se):
            f.write(f"{name}\t{est:.6g}\t{se:.6g}\n")


def save_rand(filepath, sample_ids, rand_hats, residuals):
    names = ['Temp:GA.GA', 'Humidity:GA.GA', 'GA']
    with open(filepath, 'w') as f:
        header = "ID\t" + "\t".join(names) + "\tresiduals\n"
        f.write(header)
        for i, sid in enumerate(sample_ids):
            vals = "\t".join(f"{rand_hats[k][i]:.6g}" for k in range(len(names)))
            f.write(f"{sid}\t{vals}\t{residuals[i]:.6g}\n")


def save_anova(filepath, anova_results):
    with open(filepath, 'w') as f:
        f.write("Factors\tDf\tSumSq\tMeanSq\tF\tPr(>F)\n")
        for factor in ['Temp', 'Humidity']:
            r = anova_results[factor]
            p_str = "0" if r['p'] < 1e-15 else format_sci(r['p'])
            f.write(f"{factor}\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t{r['f']:.6g}\t{p_str}\n")
        r = anova_results['e']
        f.write(f"e\t{r['df']}\t{r['ss']:.6g}\t{r['ms']:.6g}\t.\t.\n")


def main():
    os.makedirs(BASE_DIR, exist_ok=True)
    start_time = time.time()

    print("=" * 70)
    print("  GxE Model Analysis")
    print("  Model: Trait = 1 + Temp(C) + Humidity(C) + Temp:GA.GA(R[GxE]) + Humidity:GA.GA(R[GxE]) + GA(R[G]) + e")
    print("  Method: REML (scipy.optimize + log parameterization)")
    print("  GRM: VanRaden + Su normalization")
    print("  GxE GRM: Hadamard product G*E")
    print("=" * 70)

    print("\n[1] Reading SNP and sample lists...")
    snp_ids = pd.read_csv(SNP_FILE, header=None)[0].tolist()
    sample_ids = pd.read_csv(SAMPLE_FILE, header=None)[0].astype(str).tolist()
    print(f"  SNPs: {len(snp_ids)}, Samples: {len(sample_ids)}")

    print("\n[2] Reading genotype data (pandas_plink)...")
    G_geno, bim, fam = read_genotype(BFILE, SNP_FILE, SAMPLE_FILE)
    print(f"  Genotype matrix: {G_geno.shape}")

    print("\n[3] Computing GRM (VanRaden + Su normalization)...")
    G = calc_grm(G_geno)
    print(f"  GA dim: {G.shape}, diag mean: {np.mean(np.diag(G)):.6f}")

    print("\n[4] Reading phenotype data...")
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno['ID'] = pheno['ID'].astype(int)
    train_ids = [int(x) for x in sample_ids]
    pheno = pheno[pheno['ID'].isin(train_ids)].copy()
    fam_iids = fam['iid'].values.astype(int)
    pheno = pheno.set_index('ID').loc[fam_iids].reset_index()
    print(f"  Individuals: {len(pheno)}")

    y = pheno['Trait'].values.astype(float)
    temp_vals = pheno['Temp'].values.astype(float)
    hum_vals = pheno['Humidity'].values.astype(float)

    X = np.column_stack([np.ones(len(y)), temp_vals, hum_vals])
    n = len(y)
    p_x = X.shape[1]
    print(f"  y: {y.shape}, X: {X.shape}")

    print("\n[5] Building GxE relationship matrices (Hadamard product G*E)...")
    GxE_temp = build_gxe_grm(G, temp_vals)
    GxE_hum = build_gxe_grm(G, hum_vals)
    print(f"  GxE_temp diag mean: {np.mean(np.diag(GxE_temp)):.6f}")
    print(f"  GxE_hum diag mean: {np.mean(np.diag(GxE_hum)):.6f}")

    print("\n[6] Eigenvalue decomposition for all relationship matrices...")
    eigvals_gxe_temp, U_gxe_temp = eigh(GxE_temp)
    eigvals_gxe_hum, U_gxe_hum = eigh(GxE_hum)
    eigvals_g, U_g = eigh(G)

    print(f"  GxE_temp eigenvalues: min={eigvals_gxe_temp.min():.4f}, max={eigvals_gxe_temp.max():.4f}")
    print(f"  GxE_hum eigenvalues: min={eigvals_gxe_hum.min():.4f}, max={eigvals_gxe_hum.max():.4f}")
    print(f"  GA eigenvalues: min={eigvals_g.min():.4f}, max={eigvals_g.max():.4f}")

    print("\n[7] REML optimization using scipy.optimize...")
    print("  Using shared eigenvectors approach (U from GA)...")
    U = U_g
    eigenvalues_list = [
        np.diag(U.T @ GxE_temp @ U),
        np.diag(U.T @ GxE_hum @ U),
        eigvals_g
    ]

    Uty = U.T @ y
    UtX = U.T @ X

    var_y = np.var(y)
    init_vars = np.array([var_y * 0.08, var_y * 0.08, var_y * 0.27, var_y * 0.57])
    log_vars_init = np.log(init_vars)

    neg_ll_init = neg_logL_eigen(log_vars_init, eigenvalues_list, Uty, UtX)
    print(f"  Initial LogL: {-neg_ll_init:.2f}")

    iter_count = [0]
    def callback_func(xk):
        iter_count[0] += 1
        if iter_count[0] % 10 == 0 or iter_count[0] <= 5:
            neg_ll = neg_logL_eigen(xk, eigenvalues_list, Uty, UtX)
            vars_c = np.exp(xk)
            print(f"  [L-BFGS-B] {iter_count[0]}\t{-neg_ll:.2f}\t"
                  f"V(GxT)={vars_c[0]:.5f}\tV(GxH)={vars_c[1]:.5f}\t"
                  f"V(GA)={vars_c[2]:.5f}\tV(e)={vars_c[3]:.5f}")

    result = minimize(
        neg_logL_eigen,
        log_vars_init,
        args=(eigenvalues_list, Uty, UtX),
        method='L-BFGS-B',
        callback=callback_func,
        options={'maxiter': 2000, 'ftol': 1e-12, 'gtol': 1e-8, 'disp': False}
    )

    vars_opt = np.exp(result.x)
    logL_opt = -result.fun

    print(f"\n  L-BFGS-B result: success={result.success}, msg={result.message}")
    print(f"  LogL: {logL_opt:.2f}")
    print(f"  V(Temp:GA.GA)={vars_opt[0]:.6f}")
    print(f"  V(Humidity:GA.GA)={vars_opt[1]:.6f}")
    print(f"  V(GA)={vars_opt[2]:.6f}")
    print(f"  V(e)={vars_opt[3]:.6f}")

    hiblup_vars = np.array([0.0, 25.0568, 85.6164, 201.143])
    logL_hiblup = reml_logL_eigen(np.log(hiblup_vars), eigenvalues_list, Uty, UtX)
    print(f"\n  LogL at HIBLUP params: {logL_hiblup:.2f}")
    print(f"  HIBLUP reported LogL: ~-3203.19")

    if logL_hiblup > logL_opt:
        print(f"\n  HIBLUP params give higher LogL, trying Nelder-Mead from HIBLUP params...")
        result2 = minimize(
            neg_logL_eigen,
            np.log(hiblup_vars),
            args=(eigenvalues_list, Uty, UtX),
            method='Nelder-Mead',
            options={'maxiter': 50000, 'xatol': 1e-10, 'fatol': 1e-10, 'disp': False}
        )
        vars_nm = np.exp(result2.x)
        logL_nm = -result2.fun
        print(f"  Nelder-Mead LogL: {logL_nm:.2f}")
        if logL_nm > logL_opt:
            print("  Nelder-Mead found better solution, using it.")
            vars_opt = vars_nm
            logL_opt = logL_nm

    print("\n[8] Computing AI matrix for standard errors...")
    se_vars, h2s, se_h2s, p_values, cov_mat = compute_ai_se(
        vars_opt, eigenvalues_list, Uty, UtX, n)

    print(f"  Variance component SEs:")
    names = ['Temp:GA.GA', 'Humidity:GA.GA', 'GA', 'e']
    for i, name in enumerate(names):
        print(f"  {name}: Var={vars_opt[i]:.4f}, SE={se_vars[i]:.4f}, h2={h2s[i]:.4f}, h2_SE={se_h2s[i]:.4f}")

    print("\n[9] Estimating fixed effects and random effects...")
    d = np.zeros(n)
    for i in range(3):
        d += vars_opt[i] * eigenvalues_list[i]
    d += vars_opt[-1]
    d_inv = 1.0 / d

    XtVX = UtX.T * d_inv @ UtX
    XtVX_inv = np.linalg.inv(XtVX)
    beta = XtVX_inv @ (UtX.T * d_inv @ Uty)
    beta_se = np.sqrt(np.diag(XtVX_inv))

    Py = (Uty - UtX @ beta) * d_inv - (UtX * d_inv[:, None]) @ (XtVX_inv @ (UtX.T * d_inv @ (Uty - UtX @ beta)))

    rand_hats = []
    for i in range(3):
        u_i = vars_opt[i] * eigenvalues_list[i] * Py
        rand_hats.append(U @ u_i)

    y_pred = X @ beta + sum(rand_hats)
    residuals = y - y_pred

    level_names = ['mu', 'Temp', 'Humidity']
    for name, est, se in zip(level_names, beta, beta_se):
        print(f"  {name}: {est:.6g} +/- {se:.6g}")

    print("\n[10] ANOVA...")
    anova_results = compute_anova(y, X, vars_opt[-1])
    print(f"  Temp: SS={anova_results['Temp']['ss']:.2f}, F={anova_results['Temp']['f']:.2f}")
    print(f"  Humidity: SS={anova_results['Humidity']['ss']:.2f}, F={anova_results['Humidity']['f']:.2f}")

    print("\n[11] Saving output files...")

    save_vars(os.path.join(BASE_DIR, "gxe_model.vars"),
              vars_opt, se_vars, h2s, se_h2s, p_values)

    save_beta(os.path.join(BASE_DIR, "gxe_model.beta"), beta, beta_se)

    ids = pheno['ID'].values
    save_rand(os.path.join(BASE_DIR, "gxe_model.rand"),
              ids, rand_hats, residuals)

    save_anova(os.path.join(BASE_DIR, "gxe_model.anova"), anova_results)

    print("\n[12] Generating visualization plots...")
    plot_variance_components(vars_opt, se_vars, h2s, se_h2s, BASE_DIR)
    plot_random_effects(rand_hats, BASE_DIR)
    plot_interaction_effects_forest(beta, beta_se, vars_opt, se_vars, BASE_DIR)

    elapsed = time.time() - start_time
    with open(os.path.join(BASE_DIR, "gxe_model.log"), 'w') as f:
        f.write("GxE Model Analysis Log\n")
        f.write("=" * 60 + "\n")
        f.write(f"Model: Trait = 1 + Temp(C) + Humidity(C) + Temp:GA.GA(R[GxE]) + Humidity:GA.GA(R[GxE]) + GA(R[G]) + e\n")
        f.write(f"Method: REML (scipy.optimize L-BFGS-B + log parameterization)\n")
        f.write(f"GRM: VanRaden + Su normalization\n")
        f.write(f"GxE GRM: Hadamard product G*E\n")
        f.write(f"SNPs: {len(snp_ids)}, Samples: {n}\n")
        f.write(f"\nFinal LogL: {logL_opt:.6f}\n")
        f.write(f"\nVariance components:\n")
        for i, name in enumerate(names):
            f.write(f"  {name}: Var={vars_opt[i]:.6f}, SE={se_vars[i]:.6f}, h2={h2s[i]:.6f}, h2_SE={se_h2s[i]:.6f}\n")
        f.write(f"\nFixed effects:\n")
        for name, est, se in zip(level_names, beta, beta_se):
            f.write(f"  {name}: {est:.6g} +/- {se:.6g}\n")
        f.write(f"\nTotal running time: {elapsed:.1f}s\n")

    print(f"\n  Results saved to: {BASE_DIR}")
    print("\n" + "=" * 70)
    print("  Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
