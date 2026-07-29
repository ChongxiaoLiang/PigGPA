#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step T18: 连锁不平衡(LD)计算

计算 1Mb 窗口内所有配对的 LD r（Pearson 相关系数），与 HIBLUP 对齐。

核心目标：
  1. 使用全部 5556 个 chr1 SNP
  2. 使用 train_samples.txt 中的 1000 个样本（1001-2000）
  3. 计算 1Mb 窗口内所有配对（含自配对）的 LD r（Pearson 相关系数）
  4. 输出格式与 HIBLUP 完全一致：Window\tSNP_i\tSNP_j\tLD_r\tPos_i\tPos_j\tDistance

HIBLUP 窗口规则：
  window = floor(pos / 1_000_000) + 1
  每个窗口内输出所有 i <= j 的配对（上三角含对角线）
  自配对 LD_r = 1.0, Distance = 0
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from pandas_plink import read_plink

# ============================================================
# 配置
# ============================================================
BENCH = "/public/share/likui/liangcx/bole/bole_benchmark/piggpa-hiblup_benchmark_upload"
UNIFIED_DIR = os.path.join(BENCH, "unified_testdata")

BED_FILE = os.path.join(UNIFIED_DIR, "simulated_population")  # 使用 symlink
CHR1_SNPS_FILE = os.path.join(UNIFIED_DIR, "chr1_snps.txt")
TRAIN_SAMPLES_FILE = os.path.join(UNIFIED_DIR, "train_samples.txt")

OUTPUT_DIR = os.path.join(BENCH, "benchmark", "ld_calculation", "piggpa_G", "t18")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "ld_result_all.txt")

WINDOW_SIZE = 1_000_000  # 1Mb 窗口


def log(msg):
    """带时间戳的打印"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log("=" * 70)
    log("LD 计算 - HIBLUP 对齐版 (t18)")
    log("=" * 70)

    # ----------------------------------------------------------
    # [1/5] 读取 PLINK 文件
    # ----------------------------------------------------------
    log("[1/5] 读取 PLINK 文件...")
    t0 = time.time()
    bim, fam, bed = read_plink(BED_FILE, verbose=False)
    log(f"  总样本数: {len(fam)}, 总 SNP 数: {len(bim)} (耗时 {time.time()-t0:.1f}s)")

    # ----------------------------------------------------------
    # [2/5] 筛选 chr1 SNP（按 chr1_snps.txt）
    # ----------------------------------------------------------
    log("[2/5] 筛选 chr1 SNP...")
    chr1_snp_ids = pd.read_csv(CHR1_SNPS_FILE, header=None)[0].astype(str).tolist()
    chr1_snp_set = set(chr1_snp_ids)
    log(f"  chr1_snps.txt 中 SNP 数: {len(chr1_snp_ids)}")

    bim_snp = bim['snp'].astype(str)
    snp_mask = bim_snp.isin(chr1_snp_set)
    # 保留原始 bim 行索引（对应 bed 行号），再按位置排序
    snp_info = bim[snp_mask].copy()
    snp_info['snp'] = snp_info['snp'].astype(str)
    snp_info = snp_info.reset_index().rename(columns={'index': 'bim_idx'})
    snp_info = snp_info.sort_values('pos').reset_index(drop=True)
    snp_info['pos'] = snp_info['pos'].astype(int)
    # 排序后的 bed 行索引列表（保证基因型列顺序与 snp_info 一致）
    snp_indices = snp_info['bim_idx'].tolist()
    log(f"  在 bim 中匹配到的 SNP 数: {len(snp_indices)}")
    log(f"  排序后 SNP 位置范围: {snp_info['pos'].min()} - {snp_info['pos'].max()}")

    # ----------------------------------------------------------
    # [3/5] 筛选样本（按 train_samples.txt）
    # ----------------------------------------------------------
    log("[3/5] 筛选样本...")
    train_ids = pd.read_csv(TRAIN_SAMPLES_FILE, header=None)[0].astype(str).tolist()
    train_id_set = set(train_ids)
    log(f"  train_samples.txt 中样本数: {len(train_ids)}")

    fam_iid = fam['iid'].astype(str)
    sample_mask = fam_iid.isin(train_id_set)
    sample_indices = fam.index[sample_mask].tolist()
    log(f"  在 fam 中匹配到的样本数: {len(sample_indices)}")

    # ----------------------------------------------------------
    # [4/5] 提取基因型并计算窗口内 LD
    # ----------------------------------------------------------
    log("[4/5] 提取基因型数据并计算 LD...")
    t0 = time.time()
    # bed[snp_indices, :][:, sample_indices] -> shape (n_snps, n_samples)
    # 转置为 (n_samples, n_snps)
    genotype = bed[snp_indices, :][:, sample_indices].compute().T
    genotype = np.nan_to_num(genotype, nan=0)
    n_samples, n_snps = genotype.shape
    log(f"  基因型矩阵: {n_samples} 样本 x {n_snps} SNP (提取耗时 {time.time()-t0:.1f}s)")

    # 计算每个 SNP 的窗口号
    positions = snp_info['pos'].values.astype(np.int64)
    snp_names = snp_info['snp'].values.astype(str)
    windows = (positions // WINDOW_SIZE) + 1

    # 标准化基因型（每列减均值除标准差）
    log("  标准化基因型矩阵...")
    geno_centered = genotype - genotype.mean(axis=0, keepdims=True)
    geno_std = genotype.std(axis=0, ddof=0)  # 总体标准差
    # 避免除零：方差为0的SNP（单态型）std设为1，标准化后为全0，相关系数为0
    geno_std_safe = np.where(geno_std > 0, geno_std, 1.0)
    geno_stdized = geno_centered / geno_std_safe  # shape (n_samples, n_snps)

    # 逐窗口计算 LD
    log("  逐窗口计算 LD 配对...")
    t0 = time.time()
    results = []  # 每个元素: (window, snp_i, snp_j, ld_r, pos_i, pos_j, distance)
    unique_windows = np.unique(windows)
    n_windows = len(unique_windows)
    log(f"  窗口数: {n_windows}")

    for wi, w in enumerate(unique_windows, 1):
        w_mask = windows == w
        w_idx = np.where(w_mask)[0]
        n_w = len(w_idx)
        if n_w == 0:
            continue

        # 提取该窗口的标准化基因型 (n_samples, n_w)
        G = geno_stdized[:, w_idx]  # (n_samples, n_w)
        # 相关系数矩阵 = G^T @ G / (n_samples - 1) ... 实际上标准化后用 /n 即为 Pearson r
        # Pearson r = sum((x-xbar)(y-ybar)) / (sqrt(sum(x-xbar)^2) * sqrt(sum(y-ybar)^2))
        # 标准化后每列 sum((x-xbar)^2) = n_samples (因为除以了 std=sqrt(sum/n))
        # 所以 r = (G^T @ G) / n_samples
        corr = (G.T @ G) / n_samples  # (n_w, n_w)

        w_positions = positions[w_idx]
        w_names = snp_names[w_idx]

        # 提取上三角（含对角线）
        for ii in range(n_w):
            for jj in range(ii, n_w):
                pos_i = int(w_positions[ii])
                pos_j = int(w_positions[jj])
                if ii == jj:
                    ld_r = 1.0
                else:
                    ld_r = float(corr[ii, jj])
                # 限制到 [-1, 1]（浮点误差可能导致略超）
                if ld_r > 1.0:
                    ld_r = 1.0
                elif ld_r < -1.0:
                    ld_r = -1.0
                distance = pos_j - pos_i
                results.append((int(w), w_names[ii], w_names[jj],
                                ld_r, pos_i, pos_j, distance))

        if wi % 20 == 0 or wi == n_windows:
            log(f"    窗口 {w} ({wi}/{n_windows}), SNP数={n_w}, 累计配对={len(results)}")

    log(f"  LD 计算完成，总配对数: {len(results)} (耗时 {time.time()-t0:.1f}s)")

    # ----------------------------------------------------------
    # [5/5] 输出结果
    # ----------------------------------------------------------
    log("[5/5] 输出结果...")
    t0 = time.time()
    out_df = pd.DataFrame(results, columns=[
        'Window', 'SNP_i', 'SNP_j', 'LD_r', 'Pos_i', 'Pos_j', 'Distance'
    ])
    out_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    log(f"  已保存: {OUTPUT_FILE}")
    log(f"  总行数: {len(out_df)} (含数据行，不含表头)")
    log(f"  自配对数: {(out_df['SNP_i'] == out_df['SNP_j']).sum()}")
    log(f"  LD_r 范围: [{out_df['LD_r'].min():.6f}, {out_df['LD_r'].max():.6f}]")
    log(f"  写入耗时: {time.time()-t0:.1f}s")

    log("=" * 70)
    log("分析完成!")
    log("=" * 70)


if __name__ == "__main__":
    main()
