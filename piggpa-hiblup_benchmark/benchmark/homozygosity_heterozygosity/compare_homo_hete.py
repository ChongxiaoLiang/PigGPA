#!/usr/bin/env python3
"""
Step: Compare piggpa-G vs HIBLUP homozygosity/heterozygosity rates.

This script benchmarks the concordance between piggpa-G and HIBLUP (v1.2) on
homozygosity and heterozygosity rate calculations using in-silico simulated
data (chr1, 1000 samples subset of the 10000-individual HIBLUP reference set).

Logic:
  1. Load piggpa-G CSV (Sample_ID, Homozygosity_Rate, Heterozygosity_Rate).
  2. Load HIBLUP homo file (id, a1a1, a2a2), filter to ID 1-1000,
     Homo_hiblup = a1a1 + a2a2.
  3. Load HIBLUP hete file (id, a1a2), filter to ID 1-1000,
     Hete_hiblup = a1a2.
  4. Merge all on Sample_ID.
  5. For BOTH homozygosity and heterozygosity compute:
     Pearson r, Spearman rho, MSE, Max Diff.
  6. Print results to stdout (flush=True) and save:
     - homo_hete_comparison_results.csv (summary)
     - homo_hete_merged.csv (per-individual)

Usage:
    conda run -n py_analysis python compare_homo_hete.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# ---------------------------------------------------------------------------
# Absolute paths
# ---------------------------------------------------------------------------
BASE_DIR = "/public/share/likui/liangcx/bole/bole_benchmark/piggpa-hiblup_benchmark_upload/benchmark/homozygosity_heterozygosity"

PIGGPA_CSV = os.path.join(BASE_DIR, "piggpa_G", "t2", "sample_homozygosity_heterozygosity.csv")
HIBLUP_HOMO = os.path.join(BASE_DIR, "hiblup", "1.2", "chr1_homozygosity.homo")
HIBLUP_HETE = os.path.join(BASE_DIR, "hiblup", "1.2", "chr1_heterozygosity.hete")

OUT_SUMMARY = os.path.join(BASE_DIR, "homo_hete_comparison_results.csv")
OUT_MERGED = os.path.join(BASE_DIR, "homo_hete_merged.csv")


def compute_metrics(x, y):
    """Return Pearson r, Spearman rho, MSE, Max Diff between two arrays."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    mse = float(np.mean((x - y) ** 2))
    max_diff = float(np.max(np.abs(x - y)))
    return float(r), float(rho), mse, max_diff


def main():
    print("=" * 70, flush=True)
    print("piggpa-G vs HIBLUP homozygosity/heterozygosity benchmark", flush=True)
    print("=" * 70, flush=True)

    # -----------------------------------------------------------------------
    # 1. Load piggpa-G
    # -----------------------------------------------------------------------
    print("\n[1] Loading piggpa-G CSV ...", flush=True)
    df_piggpa = pd.read_csv(PIGGPA_CSV)
    print(f"    Loaded {len(df_piggpa)} rows, columns: {list(df_piggpa.columns)}", flush=True)
    print(f"    Sample_ID range: {df_piggpa['Sample_ID'].min()} - {df_piggpa['Sample_ID'].max()}", flush=True)
    df_piggpa = df_piggpa[["Sample_ID", "Homozygosity_Rate", "Heterozygosity_Rate"]].copy()
    df_piggpa = df_piggpa.rename(columns={
        "Homozygosity_Rate": "Homo_piggpa",
        "Heterozygosity_Rate": "Hete_piggpa",
    })

    # -----------------------------------------------------------------------
    # 2. Load HIBLUP homozygosity
    # -----------------------------------------------------------------------
    print("\n[2] Loading HIBLUP homozygosity file ...", flush=True)
    df_homo = pd.read_csv(HIBLUP_HOMO, sep=r"\s+", engine="python")
    print(f"    Loaded {len(df_homo)} rows, columns: {list(df_homo.columns)}", flush=True)
    df_homo = df_homo[df_homo["id"].between(1, 1000)].copy()
    df_homo["Homo_hiblup"] = df_homo["a1a1"] + df_homo["a2a2"]
    df_homo = df_homo[["id", "Homo_hiblup"]].rename(columns={"id": "Sample_ID"})
    print(f"    After filter (ID 1-1000): {len(df_homo)} rows", flush=True)

    # -----------------------------------------------------------------------
    # 3. Load HIBLUP heterozygosity
    # -----------------------------------------------------------------------
    print("\n[3] Loading HIBLUP heterozygosity file ...", flush=True)
    df_hete = pd.read_csv(HIBLUP_HETE, sep=r"\s+", engine="python")
    print(f"    Loaded {len(df_hete)} rows, columns: {list(df_hete.columns)}", flush=True)
    df_hete = df_hete[df_hete["id"].between(1, 1000)].copy()
    df_hete["Hete_hiblup"] = df_hete["a1a2"]
    df_hete = df_hete[["id", "Hete_hiblup"]].rename(columns={"id": "Sample_ID"})
    print(f"    After filter (ID 1-1000): {len(df_hete)} rows", flush=True)

    # -----------------------------------------------------------------------
    # 4. Merge
    # -----------------------------------------------------------------------
    print("\n[4] Merging datasets on Sample_ID ...", flush=True)
    df_merged = df_piggpa.merge(df_homo, on="Sample_ID", how="inner")
    df_merged = df_merged.merge(df_hete, on="Sample_ID", how="inner")
    df_merged = df_merged[["Sample_ID", "Homo_piggpa", "Homo_hiblup", "Hete_piggpa", "Hete_hiblup"]]
    print(f"    Merged rows: {len(df_merged)}", flush=True)
    print(f"    Columns: {list(df_merged.columns)}", flush=True)

    # -----------------------------------------------------------------------
    # 5. Sample-level sanity check (Sample 1)
    # -----------------------------------------------------------------------
    print("\n[5] Sample-level sanity check (Sample_ID=1) ...", flush=True)
    s1 = df_merged[df_merged["Sample_ID"] == 1].iloc[0]
    print(f"    Homo_piggpa   = {s1['Homo_piggpa']:.10f}", flush=True)
    print(f"    Homo_hiblup   = {s1['Homo_hiblup']:.10f}", flush=True)
    print(f"    Hete_piggpa   = {s1['Hete_piggpa']:.10f}", flush=True)
    print(f"    Hete_hiblup   = {s1['Hete_hiblup']:.10f}", flush=True)

    # -----------------------------------------------------------------------
    # 6. Compute comparison metrics
    # -----------------------------------------------------------------------
    print("\n[6] Computing comparison metrics ...", flush=True)
    r_homo, rho_homo, mse_homo, maxdiff_homo = compute_metrics(
        df_merged["Homo_piggpa"], df_merged["Homo_hiblup"]
    )
    r_hete, rho_hete, mse_hete, maxdiff_hete = compute_metrics(
        df_merged["Hete_piggpa"], df_merged["Hete_hiblup"]
    )

    summary_rows = [
        {
            "Module": "Homozygosity",
            "Pearson_r": r_homo,
            "Spearman_rho": rho_homo,
            "MSE": mse_homo,
            "Max_Diff": maxdiff_homo,
            "n": len(df_merged),
        },
        {
            "Module": "Heterozygosity",
            "Pearson_r": r_hete,
            "Spearman_rho": rho_hete,
            "MSE": mse_hete,
            "Max_Diff": maxdiff_hete,
            "n": len(df_merged),
        },
    ]
    df_summary = pd.DataFrame(summary_rows)

    # -----------------------------------------------------------------------
    # 7. Print results
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70, flush=True)
    print("Comparison Results", flush=True)
    print("=" * 70, flush=True)
    print(df_summary.to_string(index=False), flush=True)
    print("\n" + "-" * 70, flush=True)
    print("Homozygosity:", flush=True)
    print(f"  Pearson r    = {r_homo:.10f}", flush=True)
    print(f"  Spearman rho = {rho_homo:.10f}", flush=True)
    print(f"  MSE          = {mse_homo:.10e}", flush=True)
    print(f"  Max Diff     = {maxdiff_homo:.10e}", flush=True)
    print("\nHeterozygosity:", flush=True)
    print(f"  Pearson r    = {r_hete:.10f}", flush=True)
    print(f"  Spearman rho = {rho_hete:.10f}", flush=True)
    print(f"  MSE          = {mse_hete:.10e}", flush=True)
    print(f"  Max Diff     = {maxdiff_hete:.10e}", flush=True)

    # Pass/fail check
    print("\n" + "-" * 70, flush=True)
    pass_homo = r_homo >= 0.99
    pass_hete = r_hete >= 0.99
    print(f"Homozygosity Pearson r >= 0.99 : {'PASS' if pass_homo else 'FAIL'}", flush=True)
    print(f"Heterozygosity Pearson r >= 0.99: {'PASS' if pass_hete else 'FAIL'}", flush=True)

    # -----------------------------------------------------------------------
    # 8. Save outputs
    # -----------------------------------------------------------------------
    print("\n[7] Saving output files ...", flush=True)
    df_summary.to_csv(OUT_SUMMARY, index=False)
    df_merged.to_csv(OUT_MERGED, index=False)
    print(f"    Saved summary : {OUT_SUMMARY}", flush=True)
    print(f"    Saved merged  : {OUT_MERGED}", flush=True)
    print(f"    Merged rows (excl header): {len(df_merged)}", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
