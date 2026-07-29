#!/usr/bin/env python3
"""
Generate comprehensive correlation overview figure for piggpa-G vs HIBLUP benchmark.

Single 3x5 multi-panel PDF (fig_correlation_overview.pdf) showing:
  - 11 scatter plots (piggpa-G vs HIBLUP final results per module)
  - 4 bar charts (T3 matrix r, T8/T9 variance components, T20 LD Score regression)
  - 0 annotation panels

Layout (row, col):
  (0,0) Allele Frequency   (0,1) Inbreeding Coef.  (0,2) PCA PC1   (0,3) BLUP        (0,4) SNP Effect
  (1,0) GEBV               (1,1) Multi-Trait GA     (1,2) GxE GA    (1,3) LD r        (1,4) T3 bars
  (2,0) T8 bars            (2,1) T9 bars            (2,2) T20 bars  (2,3) Homozygosity  (2,4) Heterozygosity

Usage:
    conda run -n py_analysis python figures/scripts/plot_cns_figures.py

All paths are absolute so the script can be run from any directory.
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Global CNS (Nature/Cell/Science) style configuration
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42       # embed TrueType (editable text in PDF)
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10           # >= 8pt for all text
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.titlesize"] = 18
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["axes.linewidth"] = 1.0     # 1.0-1.5pt axis lines
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["savefig.dpi"] = 300

# Nature-inspired colour palette
COLOR_BLUE = "#4C72B0"    # piggpa-G
COLOR_GREEN = "#55A868"   # HIBLUP
COLOR_RED = "#C44E52"
COLOR_ORANGE = "#DD8452"
COLOR_GRAY = "#8C8C8C"
COLOR_PURPLE = "#8172B3"

# ---------------------------------------------------------------------------
# Paths (absolute)
# ---------------------------------------------------------------------------
BASE_DIR = "/public/share/likui/liangcx/bole/bole_benchmark/piggpa-hiblup_benchmark_upload"
BENCH_DIR = os.path.join(BASE_DIR, "benchmark")
FIG_DIR = os.path.join(BASE_DIR, "figures")

# ---------------------------------------------------------------------------
# Hardcoded verified data for bar charts
# ---------------------------------------------------------------------------

# T3 - Relationship matrix comparison
# Default HRM formula is now 0.95*G_adj + 0.05*A (HIBLUP HA equivalent)
T3_DATA = {
    "pairs": ["PRM vs PA", "GRM vs GA", "HRM vs HA"],
    "r": [1.0000, 0.9990, 0.9991],
}

# T8 - Single-trait variance components (AI-REML)
T8_DATA = {
    "components": ["V(G)", "V(e)"],
    "piggpa_G": [0.0, 1.018681],
    "hiblup": [0.0, 0.991223],
}

# T9 - Repeated records model variance components
T9_DATA = {
    "components": ["V(ID)", "V(GA)", "V(e)"],
    "piggpa_G": [85.5030, 19.5169, 14.1266],
    "hiblup": [79.7363, 24.8079, 14.1266],
}


# ===========================================================================
# Data loading functions (one per module)
# ===========================================================================
def load_allele_frequency():
    """Load allele frequency data and align by SNP ID with allele correction.

    piggpa-G MAF column is the frequency of A1.
    HIBLUP freq_a1 is the frequency of a1.
    If piggpa-G A1 != HIBLUP a1, compare piggpa-G MAF vs (1 - HIBLUP freq_a1).
    """
    pg_path = os.path.join(BENCH_DIR, "allele_frequency/piggpa_G/t1/allele_genotype_frequency.csv")
    hb_path = os.path.join(BENCH_DIR, "allele_frequency/hiblup/1.2/chr1_allele_freq.afreq")
    pg = pd.read_csv(pg_path)
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    # Merge on SNP
    merged = pg.merge(hb, on="SNP", how="inner")
    # Allele alignment: if piggpa-G A1 == HIBLUP a1 -> direct; else flip HIBLUP
    same_allele = merged["A1"] == merged["a1"]
    hb_freq = np.where(same_allele, merged["freq_a1"], 1.0 - merged["freq_a1"])
    pg_freq = merged["MAF"].values
    return hb_freq, pg_freq


def load_inbreeding():
    """Load inbreeding coefficients. Match by ID (both start at 1)."""
    pg_path = os.path.join(BENCH_DIR, "inbreeding_coefficient/piggpa_G/t5/inbreeding_coefficients.csv")
    hb_path = os.path.join(BENCH_DIR, "inbreeding_coefficient/hiblup/4/chr1_inbreeding.ibc")
    pg = pd.read_csv(pg_path)
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    pg = pg.rename(columns={"Sample_ID": "id"})
    merged = pg.merge(hb, on="id", how="inner")
    return merged["ibc"].values, merged["Inbreeding_GRM"].values


def load_pca():
    """Load PCA scores. Match by ID. Returns (hb_pc1, pg_pc1, hb_pc2, pg_pc2)."""
    pg_path = os.path.join(BENCH_DIR, "pca/piggpa_G/t6/1/chr1_pca.pc")
    hb_path = os.path.join(BENCH_DIR, "pca/hiblup/5/chr1_pca.pc")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="id", how="inner", suffixes=("_pg", "_hb"))
    return (merged["PC1_hb"].values, merged["PC1_pg"].values,
            merged["PC2_hb"].values, merged["PC2_pg"].values)


def load_blup():
    """Load BLUP breeding values. Match by id (both start at 5001)."""
    pg_path = os.path.join(BENCH_DIR, "blup_prediction/piggpa_G/GBLUP/gblup_pred.bv")
    hb_path = os.path.join(BENCH_DIR, "blup_prediction/hiblup/gblup_pred.bv")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="id", how="inner", suffixes=("_pg", "_hb"))
    return merged["add_a1_hb"].values, merged["add_a1_pg"].values


def load_snp_effect():
    """Load SNP effects. Match by SNP id. Both files have add_a1 column."""
    pg_path = os.path.join(BENCH_DIR, "snp_effect/piggpa_G/t11/snp_effect.snpeff")
    hb_path = os.path.join(BENCH_DIR, "snp_effect/hiblup/11/snp_effect.snpeff")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="id", how="inner", suffixes=("_pg", "_hb"))
    return merged["add_a1_hb"].values, merged["add_a1_pg"].values


def load_gebv():
    """Load GEBV predictions. Both tools now use same SNP effects and prediction samples (1001-2000)."""
    pg_path = os.path.join(BENCH_DIR, "gebv_prediction/piggpa_G/t16/prediction_result.bv")
    hb_path = os.path.join(BENCH_DIR, "gebv_prediction/hiblup/12/prediction_result.bv")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="id", how="inner", suffixes=("_pg", "_hb"))
    return merged["add_a1_hb"].values, merged["add_a1_pg"].values


def load_multi_trait():
    """Load multi-trait model GA. Match by ID (both start at 1001)."""
    pg_path = os.path.join(BENCH_DIR, "multi_trait_model/piggpa_G/t10/3/multi_trait.T1.rand")
    hb_path = os.path.join(BENCH_DIR, "multi_trait_model/hiblup/9/multi_trait.T1.rand")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="ID", how="inner", suffixes=("_pg", "_hb"))
    return merged["GA_hb"].values, merged["GA_pg"].values


def load_gxe():
    """Load GxE model GA. Match by ID (both start at 1001)."""
    pg_path = os.path.join(BENCH_DIR, "gxe_model/piggpa_G/t12/1/gxe_model.rand")
    hb_path = os.path.join(BENCH_DIR, "gxe_model/hiblup/10/gxe_model.rand")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    merged = pg.merge(hb, on="ID", how="inner", suffixes=("_pg", "_hb"))
    return merged["GA_hb"].values, merged["GA_pg"].values


def load_ld():
    """Load LD r values. Match by SNP pair (order-independent).
    Both tools now output Pearson r (not r^2). Filter self-pairs.
    """
    pg_path = os.path.join(BENCH_DIR, "ld_calculation/piggpa_G/t18/ld_result_all.txt")
    hb_path = os.path.join(BENCH_DIR, "ld_calculation/hiblup/14/ld_result_all.txt")
    pg = pd.read_csv(pg_path, sep=r"\s+", engine="python")
    hb = pd.read_csv(hb_path, sep=r"\s+", engine="python")
    # Filter self-pairs from both
    pg = pg[pg["SNP_i"] != pg["SNP_j"]].copy()
    hb = hb[hb["SNP_i"] != hb["SNP_j"]].copy()
    # Create order-independent pair keys
    pg["key"] = pg.apply(lambda r: tuple(sorted([r["SNP_i"], r["SNP_j"]])), axis=1)
    hb["key"] = hb.apply(lambda r: tuple(sorted([r["SNP_i"], r["SNP_j"]])), axis=1)
    merged = pg.merge(hb[["key", "LD_r"]], on="key", how="inner", suffixes=("_pg", "_hb"))
    return merged["LD_r_hb"].values, merged["LD_r_pg"].values


def load_ldreg_data():
    """Load LD Score regression comparison data (T20).

    Returns a dict with slope and h2 values for both tools:
      {"slope": {"piggpa_G": float, "hiblup": float},
       "h2":    {"piggpa_G": float, "hiblup": float}}
    """
    path = os.path.join(BENCH_DIR,
                        "ld_score_regression/hiblup/15/ldreg_comparison_results.csv")
    df = pd.read_csv(path)
    # Normalize the Parameter column for lookup
    df["Parameter_lower"] = df["Parameter"].str.lower().str.strip()

    out = {"slope": {"piggpa_G": np.nan, "hiblup": np.nan},
           "h2": {"piggpa_G": np.nan, "hiblup": np.nan}}

    # Derived slope row
    slope_row = df[df["Parameter_lower"].str.contains("derived slope", na=False)]
    if len(slope_row) > 0:
        row = slope_row.iloc[0]
        out["slope"]["piggpa_G"] = float(row["piggpa_G"])
        out["slope"]["hiblup"] = float(row["HIBLUP"])

    # h2 (heritability) row
    h2_row = df[df["Parameter_lower"].str.startswith("h2 (heritability", na=False)]
    if len(h2_row) > 0:
        row = h2_row.iloc[0]
        out["h2"]["piggpa_G"] = float(row["piggpa_G"])
        out["h2"]["hiblup"] = float(row["HIBLUP"])

    return out


def load_homozygosity_data():
    """Load homozygosity comparison data (T2). Match by Sample_ID.
    piggpa-G Homozygosity_Rate vs HIBLUP (a1a1 + a2a2).
    """
    path = os.path.join(BENCH_DIR, "homozygosity_heterozygosity/homo_hete_merged.csv")
    df = pd.read_csv(path)
    return df["Homo_hiblup"].values, df["Homo_piggpa"].values


def load_heterozygosity_data():
    """Load heterozygosity comparison data (T2). Match by Sample_ID.
    piggpa-G Heterozygosity_Rate vs HIBLUP a1a2.
    """
    path = os.path.join(BENCH_DIR, "homozygosity_heterozygosity/homo_hete_merged.csv")
    df = pd.read_csv(path)
    return df["Hete_hiblup"].values, df["Hete_piggpa"].values


# ===========================================================================
# Plotting helpers
# ===========================================================================
def _safe_corr(x, y):
    """Compute Pearson r and Spearman rho, returning (r, rho) or (nan, nan) on failure."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return np.nan, np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r, _ = stats.pearsonr(x, y)
        rho, _ = stats.spearmanr(x, y)
    return float(r), float(rho)


def plot_scatter(ax, x, y, title, xlabel="HIBLUP", ylabel="piggpa-G", note=None, diagonal=True):
    """Draw a CNS-style scatter panel with optional y=x reference line and correlation stats."""
    ax.set_title(title, fontweight="bold", fontsize=18)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)

    # Near-constant detection
    if n == 0 or np.std(y) < 1e-10 or np.std(x) < 1e-10:
        ax.text(0.5, 0.5, "Near-zero (h^2~0)\nvalues near zero" if n > 0 else "Data unavailable",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=14, color=COLOR_GRAY, style="italic")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        return {"r": np.nan, "rho": np.nan, "n": n, "note": "near-zero"}

    r, rho = _safe_corr(x, y)

    # Scatter
    ax.scatter(x, y, s=12, color=COLOR_BLUE, alpha=0.4, edgecolors="none", zorder=3)

    # Axis limits
    x_lo, x_hi = float(np.min(x)), float(np.max(x))
    y_lo, y_hi = float(np.min(y)), float(np.max(y))
    x_pad = (x_hi - x_lo) * 0.05 if x_hi > x_lo else 1.0
    y_pad = (y_hi - y_lo) * 0.05 if y_hi > y_lo else 1.0

    if diagonal:
        lo = min(x_lo, y_lo)
        hi = max(x_hi, y_hi)
        pad = (hi - lo) * 0.05 if hi > lo else 1.0
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                linestyle="--", color=COLOR_GRAY, linewidth=1.0, alpha=0.7, zorder=2)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
    else:
        ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

    # Stats annotation (top-left)
    txt = f"r = {r:.4f}\n" + chr(961) + f" = {rho:.4f}\nn = {n}"
    if note:
        txt += f"\n{note}"
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  edgecolor=COLOR_GRAY, alpha=0.85))

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return {"r": r, "rho": rho, "n": n, "note": note}


def plot_data_unavailable(ax, title):
    """Show 'Data unavailable' in a panel."""
    ax.set_title(title, fontweight="bold", fontsize=18)
    ax.text(0.5, 0.5, "Data unavailable", transform=ax.transAxes,
            ha="center", va="center", fontsize=16, color=COLOR_RED, style="italic")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_t3_bars(ax):
    """T3 relationship matrix correlation bar chart."""
    ax.set_title("T3 Relationship Matrix r", fontweight="bold", fontsize=18)
    pairs = T3_DATA["pairs"]
    r_vals = T3_DATA["r"]
    # Color by consistency level
    colors = []
    for rv in r_vals:
        if rv >= 1.0:
            colors.append(COLOR_GREEN)
        elif rv >= 0.99:
            colors.append(COLOR_BLUE)
        elif rv >= 0.95:
            colors.append(COLOR_ORANGE)
        else:
            colors.append(COLOR_RED)
    x = np.arange(len(pairs))
    bars = ax.bar(x, r_vals, color=colors, edgecolor="black", linewidth=0.8, width=0.62)
    # Threshold lines
    ax.axhline(0.99, color=COLOR_GRAY, linestyle="--", linewidth=1.0, alpha=0.7)
    ax.axhline(0.95, color=COLOR_RED, linestyle="--", linewidth=1.0, alpha=0.7)
    ax.text(len(pairs) - 0.45, 0.9905, "r=0.99", fontsize=12, color=COLOR_GRAY,
            va="bottom", ha="right")
    ax.text(len(pairs) - 0.45, 0.9505, "r=0.95", fontsize=12, color=COLOR_RED,
            va="bottom", ha="right")
    # Annotate r values above bars
    for bar, rv in zip(bars, r_vals):
        ax.text(bar.get_x() + bar.get_width() / 2.0, rv + 0.0008,
                f"{rv:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=20, ha="right", fontsize=16)
    ax.set_ylabel("Pearson r", fontsize=16)
    ax.set_ylim(0.94, 1.005)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_grouped_bars(ax, data, title, decimals=4):
    """Grouped bar chart for variance components (T8/T9)."""
    ax.set_title(title, fontweight="bold", fontsize=18)
    components = data["components"]
    pg_vals = data["piggpa_G"]
    hb_vals = data["hiblup"]
    x = np.arange(len(components))
    w = 0.35
    bars1 = ax.bar(x - w / 2, pg_vals, w, label="piggpa-G", color=COLOR_BLUE,
                   edgecolor="black", linewidth=0.8)
    bars2 = ax.bar(x + w / 2, hb_vals, w, label="HIBLUP", color=COLOR_GREEN,
                   edgecolor="black", linewidth=0.8)
    # Annotate values above bars
    ymax = max(max(pg_vals), max(hb_vals))
    offset = ymax * 0.01
    fmt = f"{{:.{decimals}f}}"
    for bar, val in zip(bars1, pg_vals):
        ax.text(bar.get_x() + bar.get_width() / 2.0, val + offset,
                fmt.format(val), ha="center", va="bottom", fontsize=10, color=COLOR_BLUE,
                fontweight="bold")
    for bar, val in zip(bars2, hb_vals):
        ax.text(bar.get_x() + bar.get_width() / 2.0, val + offset,
                fmt.format(val), ha="center", va="bottom", fontsize=10, color=COLOR_GREEN,
                fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(components, fontsize=16)
    ax.set_ylabel("Variance estimate", fontsize=16)
    ax.set_ylim(0, ymax * 1.25)
    ax.legend(loc="upper right", fontsize=12, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ldreg_bar(ax):
    """T20 LD Score Regression grouped bar chart.

    Two metric groups (Slope, h^2) with two bars per group (piggpa-G, HIBLUP).
    Slopes are directly comparable; h^2 >1 for both tools due to strong
    genetic signal in in-silico data.
    """
    ax.set_title("T20 LD Score Regression", fontweight="bold", fontsize=18)
    data = load_ldreg_data()

    groups = ["Slope", "h" + chr(178)]
    pg_vals = [data["slope"]["piggpa_G"], data["h2"]["piggpa_G"]]
    hb_vals = [data["slope"]["hiblup"], data["h2"]["hiblup"]]

    x = np.arange(len(groups))
    w = 0.35
    bars1 = ax.bar(x - w / 2, pg_vals, w, label="piggpa-G", color=COLOR_BLUE,
                   edgecolor="black", linewidth=0.8)
    bars2 = ax.bar(x + w / 2, hb_vals, w, label="HIBLUP", color=COLOR_GREEN,
                   edgecolor="black", linewidth=0.8)

    # Annotate values above bars
    ymax = max(max([v for v in pg_vals if np.isfinite(v)]),
               max([v for v in hb_vals if np.isfinite(v)]))
    offset = ymax * 0.015
    for bar, val in zip(bars1, pg_vals):
        if np.isfinite(val):
            ax.text(bar.get_x() + bar.get_width() / 2.0, val + offset,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=10,
                    color=COLOR_BLUE, fontweight="bold")
    for bar, val in zip(bars2, hb_vals):
        if np.isfinite(val):
            ax.text(bar.get_x() + bar.get_width() / 2.0, val + offset,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=10,
                    color=COLOR_GREEN, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=16)
    ax.set_ylabel("Estimate", fontsize=16)
    ax.set_ylim(0, ymax * 1.25)
    ax.legend(loc="upper left", fontsize=12, frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_summary_panel(ax, results):
    """Render a summary table of all modules and their correlation metrics."""
    ax.axis("off")
    ax.set_title("Summary: Pearson r / Spearman " + chr(961) + " across modules",
                 fontweight="bold", fontsize=11, loc="left")

    # Build table rows
    col_labels = ["Module", "Type", "n", "Pearson r", "Spearman " + chr(961), "Note"]
    rows = []
    for r in results:
        if r.get("r") is None:
            # bar chart row
            rows.append([r["name"], r.get("type", "-"), "-", "-", "-", r.get("note", "")])
        else:
            rv = r["r"]
            rho = r["rho"]
            rv_str = f"{rv:.4f}" if not np.isnan(rv) else "N/A"
            rho_str = f"{rho:.4f}" if not np.isnan(rho) else "N/A"
            n_str = str(r["n"]) if r["n"] is not None else "-"
            rows.append([r["name"], r.get("type", "scatter"), n_str, rv_str, rho_str,
                         r.get("note", "")])

    table = ax.table(cellText=rows, colLabels=col_labels, loc="center",
                     cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.5)

    # Style header
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor("#4C72B0")
        cell.set_text_props(color="white", fontweight="bold")

    # Style body: highlight r values
    n_rows = len(rows)
    for i in range(1, n_rows + 1):
        for j in range(len(col_labels)):
            cell = table[i, j]
            cell.set_edgecolor("#D0D0D0")
        # zebra striping
        if i % 2 == 0:
            for j in range(len(col_labels)):
                table[i, j].set_facecolor("#F5F7FA")
        # colour the r cell by magnitude
        r_cell = table[i, 3]
        try:
            rv = float(r_cell.get_text().get_text())
            if rv >= 0.99:
                r_cell.set_facecolor("#D5E8D4")
            elif rv >= 0.95:
                r_cell.set_facecolor("#FFE6CC")
            elif not np.isnan(rv):
                r_cell.set_facecolor("#F8CECC")
        except (ValueError, TypeError):
            pass

    # Legend explanation below the table
    legend_txt = (
        "Legend:  " + chr(9632) + " piggpa-G (#4C72B0, blue)   "
        + chr(9632) + " HIBLUP (#55A868, green)   "
        + "dashed gray line = y = x reference\n"
        "r = Pearson correlation; " + chr(961) + " = Spearman rank correlation; "
        "n = number of matched observations.\n"
        "Bar charts use hardcoded verified values (T3/T8/T9). Scatter panels use "
        "raw benchmark outputs with ID/allele alignment."
    )
    ax.text(0.5, 0.02, legend_txt, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=7.5, color="#333333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#FAFAFA",
                      edgecolor=COLOR_GRAY, alpha=0.9))


# ===========================================================================
# Main figure assembly
# ===========================================================================
def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    print("Generating comprehensive correlation overview figure...", flush=True)
    print(f"  Output directory: {FIG_DIR}", flush=True)

    fig = plt.figure(figsize=(24, 15))
    # 3x5 grid: 11 scatter + 4 bar cells
    gs = GridSpec(3, 5, figure=fig, hspace=0.50, wspace=0.40,
                  left=0.05, right=0.97, top=0.93, bottom=0.05)

    results = []  # collect summary metrics

    # ------------------------------------------------------------------
    # Row 0
    # ------------------------------------------------------------------
    # (0,0) Allele Frequency
    ax = fig.add_subplot(gs[0, 0])
    try:
        x, y = load_allele_frequency()
        res = plot_scatter(ax, x, y, "Allele Frequency")
        res["name"] = "Allele Frequency (T1/T2)"
        results.append(res)
        print(f"  [Allele Frequency] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "Allele Frequency")
        results.append({"name": "Allele Frequency (T1/T2)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [Allele Frequency] ERROR: {e}", file=sys.stderr, flush=True)

    # (0,1) Inbreeding Coefficient
    ax = fig.add_subplot(gs[0, 1])
    try:
        x, y = load_inbreeding()
        res = plot_scatter(ax, x, y, "Inbreeding Coef.", diagonal=False)
        res["name"] = "Inbreeding Coef. (T5)"
        results.append(res)
        print(f"  [Inbreeding] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "Inbreeding Coef.")
        results.append({"name": "Inbreeding Coef. (T5)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [Inbreeding] ERROR: {e}", file=sys.stderr, flush=True)

    # (0,2) PCA PC1 (PC2 removed; load_pca still returns 4 values, ignore PC2)
    ax_pc1 = fig.add_subplot(gs[0, 2])
    try:
        hb1, pg1, _hb2, _pg2 = load_pca()
        # PCA sign-flip handling: if r < 0, flip piggpa-G signs and recompute
        r1, rho1 = _safe_corr(pg1, hb1)
        note1 = None
        if not np.isnan(r1) and r1 < 0:
            pg1 = -pg1
            note1 = "sign flipped"
            r1, rho1 = _safe_corr(pg1, hb1)
        res1 = plot_scatter(ax_pc1, hb1, pg1, "PCA PC1", note=note1, diagonal=False)
        res1["r"], res1["rho"] = r1, rho1
        res1["name"] = "PCA PC1 (T6)"
        results.append(res1)
        print(f"  [PCA PC1] r={r1:.4f}, rho={rho1:.4f}, n={res1['n']}, note={note1}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax_pc1, "PCA PC1")
        results.append({"name": "PCA PC1 (T6)", "r": np.nan, "rho": np.nan, "n": None,
                        "note": f"load error: {e}"})
        print(f"  [PCA] ERROR: {e}", file=sys.stderr, flush=True)

    # (0,3) BLUP Breeding Value
    ax = fig.add_subplot(gs[0, 3])
    try:
        x, y = load_blup()
        res = plot_scatter(ax, x, y, "BLUP Breeding Value", diagonal=False)
        res["name"] = "BLUP Breeding Value (T7)"
        results.append(res)
        print(f"  [BLUP] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "BLUP Breeding Value")
        results.append({"name": "BLUP Breeding Value (T7)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [BLUP] ERROR: {e}", file=sys.stderr, flush=True)

    # (0,4) SNP Effect
    ax = fig.add_subplot(gs[0, 4])
    try:
        x, y = load_snp_effect()
        res = plot_scatter(ax, x, y, "SNP Effect (T11)")
        res["name"] = "SNP Effect (T11)"
        results.append(res)
        print(f"  [SNP Effect] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "SNP Effect (T11)")
        results.append({"name": "SNP Effect (T11)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [SNP Effect] ERROR: {e}", file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    # Row 1
    # ------------------------------------------------------------------
    # (1,0) GEBV Prediction
    ax = fig.add_subplot(gs[1, 0])
    try:
        x, y = load_gebv()
        # GEBV sign-flip handling: allele coding convention difference (A1 vs A2 count)
        # x=HIBLUP (A2 count), y=piggpa-G (A1 count); flip y if r < 0
        r_gebv, rho_gebv = _safe_corr(x, y)
        note_gebv = None
        if not np.isnan(r_gebv) and r_gebv < 0:
            y = -y
            note_gebv = "sign flipped (allele coding)"
            r_gebv, rho_gebv = _safe_corr(x, y)
        res = plot_scatter(ax, x, y, "GEBV Prediction (T16)", note=note_gebv, diagonal=False)
        res["r"], res["rho"] = r_gebv, rho_gebv
        res["name"] = "GEBV Prediction (T16)"
        results.append(res)
        print(f"  [GEBV] r={r_gebv:.4f}, rho={rho_gebv:.4f}, n={res['n']}, note={note_gebv}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "GEBV Prediction (T16)")
        results.append({"name": "GEBV Prediction (T16)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [GEBV] ERROR: {e}", file=sys.stderr, flush=True)

    # (1,1) Multi-Trait GA
    ax = fig.add_subplot(gs[1, 1])
    try:
        x, y = load_multi_trait()
        note = "scale discrepancy" if (np.std(y) < 1e-6 or np.std(x) < 1e-6) else None
        res = plot_scatter(ax, x, y, "Multi-Trait GA (T10)", note=note, diagonal=False)
        res["name"] = "Multi-Trait GA (T10)"
        results.append(res)
        print(f"  [Multi-Trait] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "Multi-Trait GA (T10)")
        results.append({"name": "Multi-Trait GA (T10)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [Multi-Trait] ERROR: {e}", file=sys.stderr, flush=True)

    # (1,2) GxE Model GA
    ax = fig.add_subplot(gs[1, 2])
    try:
        x, y = load_gxe()
        res = plot_scatter(ax, x, y, "GxE Model GA (T12)", diagonal=False)
        res["name"] = "GxE Model GA (T12)"
        results.append(res)
        print(f"  [GxE] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "GxE Model GA (T12)")
        results.append({"name": "GxE Model GA (T12)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [GxE] ERROR: {e}", file=sys.stderr, flush=True)

    # (1,3) LD r
    ax = fig.add_subplot(gs[1, 3])
    try:
        x, y = load_ld()
        res = plot_scatter(ax, x, y, "LD r (T18)", diagonal=False)
        res["name"] = "LD r (T18)"
        results.append(res)
        print(f"  [LD] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "LD r (T18)")
        results.append({"name": "LD r (T18)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [LD] ERROR: {e}", file=sys.stderr, flush=True)

    # (1,4) T3 Relationship Matrix r (bar chart)
    ax = fig.add_subplot(gs[1, 4])
    plot_t3_bars(ax)
    results.append({"name": "T3 Relationship Matrix r", "r": None, "rho": None,
                    "n": None, "type": "bar", "note": "4 matrix pairs; all r>=0.99"})
    print("  [T3 bars] rendered (hardcoded verified values)", flush=True)

    # ------------------------------------------------------------------
    # Row 2
    # ------------------------------------------------------------------
    # (2,0) T8 Variance Components (grouped bar chart)
    ax = fig.add_subplot(gs[2, 0])
    plot_grouped_bars(ax, T8_DATA, "T8 Variance Components")
    results.append({"name": "T8 Variance Components", "r": None, "rho": None,
                    "n": None, "type": "bar", "note": "V(G)=0 both; V(e) diff~0.027"})
    print("  [T8 bars] rendered (hardcoded verified values)", flush=True)

    # (2,1) T9 Variance Components (grouped bar chart)
    ax = fig.add_subplot(gs[2, 1])
    plot_grouped_bars(ax, T9_DATA, "T9 Variance Components", decimals=3)
    results.append({"name": "T9 Variance Components", "r": None, "rho": None,
                    "n": None, "type": "bar", "note": "V(ID)/V(GA)/V(e) comparison"})
    print("  [T9 bars] rendered (hardcoded verified values)", flush=True)

    # (2,2) LD Score Regression bar chart (T20)
    ax = fig.add_subplot(gs[2, 2])
    try:
        plot_ldreg_bar(ax)
        results.append({"name": "T20 LD Score Regression", "r": None, "rho": None,
                        "n": None, "type": "bar", "note": "Slope and h2 comparison"})
        print("  [T20 LD Score Reg] rendered (slope comparison)", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "T20 LD Score Regression")
        print(f"  [T20 LD Score Reg] ERROR: {e}", file=sys.stderr, flush=True)

    # (2,3) Homozygosity (T2)
    ax = fig.add_subplot(gs[2, 3])
    try:
        x, y = load_homozygosity_data()
        res = plot_scatter(ax, x, y, "Homozygosity (T2)", diagonal=False)
        res["name"] = "Homozygosity (T2)"
        results.append(res)
        print(f"  [Homozygosity] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "Homozygosity (T2)")
        results.append({"name": "Homozygosity (T2)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [Homozygosity] ERROR: {e}", file=sys.stderr, flush=True)

    # (2,4) Heterozygosity (T2)
    ax = fig.add_subplot(gs[2, 4])
    try:
        x, y = load_heterozygosity_data()
        res = plot_scatter(ax, x, y, "Heterozygosity (T2)", diagonal=False)
        res["name"] = "Heterozygosity (T2)"
        results.append(res)
        print(f"  [Heterozygosity] r={res['r']:.4f}, rho={res['rho']:.4f}, n={res['n']}", flush=True)
    except Exception as e:
        plot_data_unavailable(ax, "Heterozygosity (T2)")
        results.append({"name": "Heterozygosity (T2)", "r": np.nan, "rho": np.nan,
                        "n": None, "note": f"load error: {e}"})
        print(f"  [Heterozygosity] ERROR: {e}", file=sys.stderr, flush=True)

    # ------------------------------------------------------------------
    # Super title & save
    # ------------------------------------------------------------------
    fig.suptitle("piggpa-G vs HIBLUP: Comprehensive Correlation Overview",
                 fontsize=20, fontweight="bold", y=0.98)

    out_path = os.path.join(FIG_DIR, "fig_correlation_overview.pdf")
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    size_kb = os.path.getsize(out_path) / 1024.0
    print(f"\nFigure generated successfully: {out_path}", flush=True)
    print(f"  File size: {size_kb:.1f} KB", flush=True)

    # Print summary table of all scatter modules
    print("\n===== Correlation Summary (scatter modules) =====", flush=True)
    print(f"{'Module':<30} {'n':>6} {'Pearson r':>12} {'Spearman rho':>14}  Note", flush=True)
    print("-" * 90, flush=True)
    for r in results:
        if r.get("r") is None:
            continue
        rv = r["r"]
        rho = r["rho"]
        rv_str = f"{rv:.4f}" if not np.isnan(rv) else "N/A"
        rho_str = f"{rho:.4f}" if not np.isnan(rho) else "N/A"
        n_str = str(r["n"]) if r["n"] is not None else "-"
        note = r.get("note", "") or ""
        print(f"{r['name']:<30} {n_str:>6} {rv_str:>12} {rho_str:>14}  {note}", flush=True)
    print("=" * 90, flush=True)

    return out_path


if __name__ == "__main__":
    main()
