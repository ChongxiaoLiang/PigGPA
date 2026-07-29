#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
关系矩阵构建：PRM, GRM, HRM, 环境随机效应矩阵
灵活版：通过命令行参数指定输入输出
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
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

try:
    from pandas_plink import read_plink
except ImportError:
    print("错误: 未安装 pandas-plink，请先运行: pip install pandas-plink")
    sys.exit(1)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 14

# Unified color palettes
PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)


class Pedigree:
    def __init__(self, sample_ids, sire_ids=None, dam_ids=None):
        self.sample_ids = list(sample_ids)
        self.n = len(self.sample_ids)
        self.sample_idx = {sid: i for i, sid in enumerate(self.sample_ids)}
        self.sire_ids = sire_ids if sire_ids else [None] * self.n
        self.dam_ids = dam_ids if dam_ids else [None] * self.n

    def _topological_sort(self):
        id_set = set(self.sample_ids)
        visited = set()
        result = []
        def visit(ind):
            if ind in visited:
                return
            visited.add(ind)
            idx = self.sample_idx.get(ind)
            if idx is not None:
                sire = self.sire_ids[idx]
                dam = self.dam_ids[idx]
                if sire and sire in id_set:
                    visit(sire)
                if dam and dam in id_set:
                    visit(dam)
            result.append(ind)
        for ind in self.sample_ids:
            visit(ind)
        return result

    def calculate_prm(self):
        sorted_ids = self._topological_sort()
        n = self.n
        prm = np.zeros((n, n))
        for ind in sorted_ids:
            i = self.sample_idx[ind]
            sire = self.sire_ids[i]
            dam = self.dam_ids[i]
            sire_idx = self.sample_idx.get(sire) if sire else None
            dam_idx = self.sample_idx.get(dam) if dam else None
            if sire_idx is not None and dam_idx is not None:
                prm[i, i] = 1 + 0.5 * prm[sire_idx, dam_idx]
            else:
                prm[i, i] = 1.0
            for j_id in sorted_ids:
                j = self.sample_idx[j_id]
                if j == i:
                    break
                val = 0.0
                if sire_idx is not None:
                    val += prm[sire_idx, j]
                if dam_idx is not None:
                    val += prm[dam_idx, j]
                prm[i, j] = 0.5 * val
                prm[j, i] = prm[i, j]
        return prm


class RelationshipMatrixBuilder:
    def __init__(self, genotype_matrix, sample_ids, bim=None, maf_threshold=0.0):
        self.genotype_matrix = genotype_matrix
        self.sample_ids = sample_ids
        self.n_samples = len(sample_ids)
        self.n_snps = genotype_matrix.shape[0]
        self.bim = bim
        self.maf_threshold = maf_threshold

    def _preprocess_genotype(self):
        """GRM 预处理：计算等位基因频率、过滤 SNP、填充缺失值、中心化"""
        M = self.genotype_matrix.T
        p = np.nanmean(M, axis=0) / 2
        valid_snp = (p >= self.maf_threshold) & (p <= 1 - self.maf_threshold)
        M = M[:, valid_snp]
        p = p[valid_snp]
        p = np.where(np.isnan(p), 0, p)
        M = np.where(np.isnan(M), 2 * p, M)
        Z = M - 2 * p
        snp_count = int(np.sum(valid_snp))
        return M, Z, p, valid_snp, snp_count

    def calculate_grm_vanraden(self):
        print("  计算GRM (VanRaden方法)...")
        M, Z, p, valid_snp, snp_count = self._preprocess_genotype()
        denom = 2 * np.sum(p * (1 - p))
        grm = (Z @ Z.T) / denom if denom > 0 else np.zeros((self.n_samples, self.n_samples))
        print(f"  有效SNP数: {snp_count}")
        return grm, snp_count

    def calculate_grm_yang(self):
        print("  计算GRM (Yang方法)...")
        M, Z, p, valid_snp, snp_count = self._preprocess_genotype()
        n = self.n_samples
        n_valid = n - np.isnan(self.genotype_matrix.T[:, valid_snp]).sum(axis=0)
        n_valid = np.where(n_valid < 2, 2, n_valid)
        denom_per_snp = 2 * p * (1 - p) * n_valid / (n_valid - 1)
        denom_per_snp = np.where(denom_per_snp == 0, 1, denom_per_snp)
        Z_scaled = Z / np.sqrt(denom_per_snp)
        grm = (Z_scaled @ Z_scaled.T) / snp_count if snp_count > 0 else np.zeros((n, n))
        print(f"  有效SNP数: {snp_count}")
        return grm, snp_count

    def calculate_hrm(self, grm, prm=None):
        """构建混合关系矩阵 (HRM)
        默认公式: HRM = 0.95×G_adj + 0.05×A (与HIBLUP HA公式等价)
        其中 G_adj = 0.999×G + 0.001×I, 使对角线均值≈1.0
        """
        if prm is None:
            prm = np.eye(self.n_samples)
        identity_mat = np.eye(self.n_samples)
        G_adj = 0.999 * grm + 0.001 * identity_mat
        return 0.95 * G_adj + 0.05 * prm

    def calculate_hrm_ssblup(self, grm, prm, lambda_g=0.5):
        return lambda_g * grm + (1 - lambda_g) * prm


class EnvironmentalMatrix:
    def __init__(self, n_samples, sample_ids):
        self.n_samples = n_samples
        self.sample_ids = sample_ids

    def create_identity_matrix(self):
        return np.eye(self.n_samples)

    def create_block_matrix(self, group_labels):
        group_labels = np.asarray(group_labels)
        return (group_labels[:, None] == group_labels[None, :]).astype(float)

    def create_spatial_matrix(self, coordinates, decay_rate=0.1):
        coords = np.array(coordinates)
        dist_matrix = cdist(coords, coords, metric='euclidean')
        return np.exp(-decay_rate * dist_matrix)

    def create_temporal_matrix(self, time_points, decay_rate=0.1):
        time_arr = np.array(time_points).reshape(-1, 1)
        time_diff = np.abs(time_arr - time_arr.T)
        return np.exp(-decay_rate * time_diff)


def calculate_matrix_statistics(matrix):
    diag = np.diag(matrix)
    upper_tri = matrix[np.triu_indices(len(matrix), k=1)]
    return {
        'mean': np.mean(matrix), 'std': np.std(matrix),
        'min': np.min(matrix), 'max': np.max(matrix),
        'diag_mean': np.mean(diag), 'diag_std': np.std(diag),
        'off_diag_mean': np.mean(upper_tri), 'off_diag_std': np.std(upper_tri)
    }


def plot_matrix_heatmap(matrix, ax, title, vmin=None, vmax=None, cmap=WARM_COOL_CMAP):
    if vmin is None:
        vmin = np.min(matrix)
    if vmax is None:
        vmax = np.max(matrix)
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Sample Index')
    ax.set_title(title)
    ax.tick_params(direction='out')
    return im


def plot_distribution(data, ax, title, xlabel, color, bins=50):
    ax.hist(data, bins=bins, color=color, edgecolor='none', alpha=0.85)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out')
    sns.despine(ax=ax)


def main():
    parser = argparse.ArgumentParser(description='关系矩阵构建')
    parser.add_argument('-i', '--input', required=True, help='PLINK文件前缀(如 /path/to/data)')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('--pedigree', default=None, help='系谱文件路径(可选)')
    parser.add_argument('--chrom', default=None, help='指定染色体(默认全部)')
    parser.add_argument('--maf', type=float, default=0.0, help='MAF过滤阈值(默认0.0，不过滤)')
    parser.add_argument('--lambda_g', type=float, default=0.7, help='SSBLUP中GRM的权重lambda_g(默认0.7)')
    parser.add_argument('--spatial_coords', default=None, help='空间坐标文件路径(CSV, n_samples行2列, 可选)')
    parser.add_argument('--time_points', default=None, help='时间点文件路径(CSV, n_samples行1列, 可选)')
    parser.add_argument('--env_block', default=None, help='区组标签文件路径(CSV, n_samples行1列, 可选)')
    parser.add_argument('--font_size', type=int, default=12, help='字体大小(默认12)')
    parser.add_argument('--dpi', type=int, default=300, help='图片DPI(默认300)')
    args = parser.parse_args()
    plt.rcParams.update({
        'font.size': args.font_size,
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
    })

    bed_file = args.input
    output_dir = args.output
    pedigree_file = args.pedigree
    target_chrom = args.chrom

    os.makedirs(output_dir, exist_ok=True)

    # 输入文件校验
    for ext in ['.bed', '.bim', '.fam']:
        fpath = bed_file + ext
        if not os.path.exists(fpath):
            print(f"错误: PLINK文件不存在: {fpath}")
            sys.exit(1)

    print("=" * 70)
    print("关系矩阵构建: PRM, GRM, HRM, 环境随机效应矩阵")
    print(f"输入: {bed_file}")
    print(f"输出: {output_dir}")
    if target_chrom:
        print(f"指定染色体: {target_chrom}")
    else:
        print("分析全部染色体")
    print(f"MAF阈值: {args.maf}")
    print(f"HRM公式: 0.95×G_adj + 0.05×A (与HIBLUP HA等价, G_adj=0.999×G+0.001×I)")
    print(f"SSBLUP参数: lambda_g={args.lambda_g}")
    print("=" * 70)

    print("\n[1/5] 读取PLINK文件...")
    bim, fam, bed = read_plink(bed_file, verbose=False)

    if target_chrom:
        chrom_mask = bim['chrom'].astype(str) == target_chrom
        snp_indices = np.where(chrom_mask)[0]
    else:
        snp_indices = np.arange(len(bim))

    n_snps = len(snp_indices)
    n_samples = len(fam)
    sample_ids = fam['iid'].values

    print(f"  SNP数量: {n_snps}")
    print(f"  样本数量: {n_samples}")

    print("\n[2/5] 加载基因型数据...")
    genotype_matrix = bed[snp_indices, :].compute()

    print("\n[3/5] 构建各类关系矩阵...")

    print("\n  === 构建家系关系矩阵 (PRM) ===")
    has_pedigree = False
    family_ids = None
    if pedigree_file and os.path.exists(pedigree_file):
        try:
            pedigree_df = pd.read_csv(pedigree_file, sep='\t')
            pedigree_df = pedigree_df[pedigree_df['ID'].isin(sample_ids)]
            # 使用 set_index 优化查找性能
            ped_indexed = pedigree_df.set_index('ID')
            sire_ids = []
            dam_ids = []
            family_ids = np.ones(n_samples, dtype=int)
            for i, sid in enumerate(sample_ids):
                if sid in ped_indexed.index:
                    row = ped_indexed.loc[sid]
                    sire = row['Sire'] if isinstance(row, pd.Series) else row['Sire'].values[0]
                    dam = row['Dam'] if isinstance(row, pd.Series) else row['Dam'].values[0]
                    sire_ids.append(str(sire) if sire != 0 and str(sire) in sample_ids else None)
                    dam_ids.append(str(dam) if dam != 0 and str(dam) in sample_ids else None)
                    # 使用确定性映射构建 family_ids
                    if sire != 0 and str(sire) in sample_ids:
                        try:
                            family_ids[i] = int(sire) % n_samples + 1
                        except (ValueError, TypeError):
                            family_ids[i] = abs(hash(str(sire))) % n_samples + 1
                else:
                    sire_ids.append(None)
                    dam_ids.append(None)
            pedigree = Pedigree(sample_ids, sire_ids, dam_ids)
            prm = pedigree.calculate_prm()
            print(f"  PRM构建完成: {prm.shape}")
            has_pedigree = True
        except Exception as e:
            print(f"  无法加载家系数据: {e}")

    if not has_pedigree:
        print("  警告: 未提供系谱文件，使用模拟家系数据！结果仅供演示，不代表真实遗传关系。")
        np.random.seed(42)
        n_families = max(1, n_samples // 10)
        family_ids = np.random.randint(1, n_families + 1, n_samples)
        sire_ids = []
        dam_ids = []
        for i in range(n_samples):
            if family_ids[i] <= n_families // 2:
                sire_ids.append(None)
                dam_ids.append(None)
            else:
                parent_family = family_ids[i] % (n_families // 2) + 1
                sire_ids.append(f"SIRE_{parent_family}")
                dam_ids.append(f"DAM_{parent_family}")
        pedigree = Pedigree(sample_ids, sire_ids, dam_ids)
        prm = pedigree.calculate_prm()
        print(f"  PRM构建完成: {prm.shape}")

    print("\n  === 构建基因组关系矩阵 (GRM) ===")
    builder = RelationshipMatrixBuilder(genotype_matrix, sample_ids, bim, maf_threshold=args.maf)
    grm_vanraden, snp_count_vr = builder.calculate_grm_vanraden()
    print(f"  VanRaden GRM构建完成: 有效SNP数 {snp_count_vr}")
    grm_yang, snp_count_yang = builder.calculate_grm_yang()
    print(f"  Yang GRM构建完成: 有效SNP数 {snp_count_yang}")

    print("\n  === 构建混合关系矩阵 (HRM) ===")
    hrm = builder.calculate_hrm(grm_vanraden, prm)
    hrm_diag_mean = float(np.mean(np.diag(hrm)))
    print(f"  HRM构建完成: HA = 0.95*G_adj + 0.05*A (G_adj = 0.999*G + 0.001*I)")
    print(f"  HRM对角线均值: {hrm_diag_mean:.4f} (期望≈1.0)")
    hrm_ssblup = builder.calculate_hrm_ssblup(grm_vanraden, prm, lambda_g=args.lambda_g)
    print(f"  HRM (SSBLUP) 构建完成 (lambda_g={args.lambda_g})")

    print("\n  === 构建环境随机效应矩阵 ===")
    env_builder = EnvironmentalMatrix(n_samples, sample_ids)
    env_identity = env_builder.create_identity_matrix()
    print(f"  单位矩阵构建完成")

    # 区组矩阵：优先使用用户提供的区组文件，其次使用系谱家系，最后随机
    if args.env_block and os.path.exists(args.env_block):
        block_labels = pd.read_csv(args.env_block, header=None).values.flatten()
        env_block = env_builder.create_block_matrix(block_labels)
        print(f"  区组矩阵构建完成 (使用用户提供的区组标签)")
    elif has_pedigree and family_ids is not None:
        env_block = env_builder.create_block_matrix(family_ids)
        print(f"  区组矩阵构建完成 (使用系谱家系分组)")
    else:
        block_labels = np.random.randint(1, max(2, n_samples // 10) + 1, n_samples)
        env_block = env_builder.create_block_matrix(block_labels)
        print(f"  区组矩阵构建完成 (警告: 使用随机分组，仅供演示)")

    # 空间矩阵：优先使用用户提供的坐标文件
    if args.spatial_coords and os.path.exists(args.spatial_coords):
        coords = pd.read_csv(args.spatial_coords, header=None).values
        env_spatial = env_builder.create_spatial_matrix(coords, decay_rate=0.05)
        print(f"  空间矩阵构建完成 (使用用户提供的坐标)")
    else:
        print("  警告: 未提供空间坐标，使用随机坐标生成空间矩阵，结果仅供演示！")
        coords = np.random.rand(n_samples, 2) * 100
        env_spatial = env_builder.create_spatial_matrix(coords, decay_rate=0.05)
        print(f"  空间矩阵构建完成 (随机坐标)")

    # 时间矩阵：优先使用用户提供的时间点文件
    if args.time_points and os.path.exists(args.time_points):
        time_points = pd.read_csv(args.time_points, header=None).values.flatten()
        env_temporal = env_builder.create_temporal_matrix(time_points, decay_rate=0.2)
        print(f"  时间矩阵构建完成 (使用用户提供的时间点)")
    else:
        print("  警告: 未提供时间点，使用随机时间生成时间矩阵，结果仅供演示！")
        time_points = np.random.randint(1, 11, n_samples)
        env_temporal = env_builder.create_temporal_matrix(time_points, decay_rate=0.2)
        print(f"  时间矩阵构建完成 (随机时间)")

    print("\n[4/5] 保存结果...")
    pd.DataFrame(prm, index=sample_ids, columns=sample_ids).to_csv(os.path.join(output_dir, "PRM.csv"))
    print(f"  PRM已保存")
    pd.DataFrame(grm_vanraden, index=sample_ids, columns=sample_ids).to_csv(os.path.join(output_dir, "GRM_VanRaden.csv"))
    print(f"  VanRaden GRM已保存")
    pd.DataFrame(grm_yang, index=sample_ids, columns=sample_ids).to_csv(os.path.join(output_dir, "GRM_Yang.csv"))
    print(f"  Yang GRM已保存")
    pd.DataFrame(hrm, index=sample_ids, columns=sample_ids).to_csv(os.path.join(output_dir, "HRM.csv"))
    print(f"  HRM已保存 (公式: 0.95*G_adj + 0.05*A, 对角线均值≈{np.mean(np.diag(hrm)):.4f})")
    pd.DataFrame(hrm_ssblup, index=sample_ids, columns=sample_ids).to_csv(os.path.join(output_dir, "HRM_SSBLUP.csv"))
    print(f"  HRM SSBLUP已保存")
    np.savez(os.path.join(output_dir, "Environmental_Matrices.npz"),
             identity=env_identity, block=env_block, spatial=env_spatial, temporal=env_temporal)
    print(f"  环境矩阵已保存")

    stats_prm = calculate_matrix_statistics(prm)
    stats_grm_vr = calculate_matrix_statistics(grm_vanraden)
    stats_grm_yang = calculate_matrix_statistics(grm_yang)
    stats_hrm = calculate_matrix_statistics(hrm)
    stats_hrm_ssblup = calculate_matrix_statistics(hrm_ssblup)

    output_summary = os.path.join(output_dir, "relationship_matrix_summary.txt")
    with open(output_summary, 'w') as f:
        f.write("Relationship Matrix Construction Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Chromosome: {target_chrom if target_chrom else 'All'}\n")
        f.write(f"Total samples: {n_samples}\n")
        f.write(f"Total SNPs: {n_snps}\n")
        f.write(f"MAF threshold: {args.maf}\n")
        f.write(f"Valid SNPs (VanRaden): {snp_count_vr}\n")
        f.write(f"Valid SNPs (Yang): {snp_count_yang}\n")
        f.write(f"HRM formula: 0.95*G_adj + 0.05*A (HIBLUP HA equivalent, G_adj=0.999*G+0.001*I)\n")
        f.write(f"SSBLUP parameters: lambda_g={args.lambda_g}\n")
        f.write("\n")
        matrix_stats_list = [("PRM", stats_prm), ("GRM-VanRaden", stats_grm_vr),
                             ("GRM-Yang", stats_grm_yang), ("HRM", stats_hrm),
                             ("HRM-SSBLUP", stats_hrm_ssblup)]
        for name, stats in matrix_stats_list:
            f.write(f"{name} Statistics:\n")
            f.write(f"  Mean: {stats['mean']:.4f}\n")
            f.write(f"  Std: {stats['std']:.4f}\n")
            f.write(f"  Range: [{stats['min']:.4f}, {stats['max']:.4f}]\n")
            f.write(f"  Diagonal mean: {stats['diag_mean']:.4f}\n")
            f.write(f"  Off-diagonal mean: {stats['off_diag_mean']:.4f}\n\n")
    print(f"  Summary saved: {output_summary}")

    output_summary_zh = os.path.join(output_dir, "relationship_matrix_summary-zh.txt")
    with open(output_summary_zh, 'w') as f:
        f.write("关系矩阵构建摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"分析染色体: {target_chrom if target_chrom else '全部'}\n")
        f.write(f"总样本数: {n_samples}\n")
        f.write(f"总SNP数: {n_snps}\n")
        f.write(f"MAF阈值: {args.maf}\n")
        f.write(f"有效SNP数 (VanRaden): {snp_count_vr}\n")
        f.write(f"有效SNP数 (Yang): {snp_count_yang}\n")
        f.write(f"HRM公式: 0.95*G_adj + 0.05*A (与HIBLUP HA等价, G_adj=0.999*G+0.001*I)\n")
        f.write(f"SSBLUP参数: lambda_g={args.lambda_g}\n")
        f.write("\n")
        matrix_stats_list_zh = [("PRM", stats_prm), ("GRM-VanRaden", stats_grm_vr),
                                ("GRM-Yang", stats_grm_yang), ("HRM", stats_hrm),
                                ("HRM-SSBLUP", stats_hrm_ssblup)]
        for name, stats in matrix_stats_list_zh:
            f.write(f"{name} 统计:\n")
            f.write(f"  均值: {stats['mean']:.4f}\n")
            f.write(f"  标准差: {stats['std']:.4f}\n")
            f.write(f"  范围: [{stats['min']:.4f}, {stats['max']:.4f}]\n")
            f.write(f"  对角线均值: {stats['diag_mean']:.4f}\n")
            f.write(f"  非对角线均值: {stats['off_diag_mean']:.4f}\n\n")
    print(f"  中文摘要已保存: {output_summary_zh}")

    print("\n[5/5] 生成可视化图表...")
    fig1, axes1 = plt.subplots(2, 3, figsize=(18, 12))
    fig1.suptitle(f'Relationship Matrix Overview (n={n_samples})', fontweight='bold', fontfamily='DejaVu Sans')
    im1 = plot_matrix_heatmap(prm, axes1[0, 0], 'Pedigree Relationship Matrix (PRM)', cmap=WARM_COOL_CMAP)
    plt.colorbar(im1, ax=axes1[0, 0], shrink=0.8)
    im2 = plot_matrix_heatmap(grm_vanraden, axes1[0, 1], 'GRM-VanRaden', vmin=-0.5, vmax=1.5, cmap=WARM_COOL_CMAP)
    plt.colorbar(im2, ax=axes1[0, 1], shrink=0.8)
    im3 = plot_matrix_heatmap(grm_yang, axes1[0, 2], 'GRM-Yang', vmin=-0.5, vmax=1.5, cmap=WARM_COOL_CMAP)
    plt.colorbar(im3, ax=axes1[0, 2], shrink=0.8)
    im4 = plot_matrix_heatmap(hrm, axes1[1, 0], 'HRM (0.95*G_adj+0.05*A)', cmap=WARM_COOL_CMAP)
    plt.colorbar(im4, ax=axes1[1, 0], shrink=0.8)
    im5 = plot_matrix_heatmap(hrm_ssblup, axes1[1, 1], 'HRM (SSBLUP)', cmap=WARM_COOL_CMAP)
    plt.colorbar(im5, ax=axes1[1, 1], shrink=0.8)
    kinship_prm = prm[np.triu_indices(n_samples, k=1)]
    kinship_grm = grm_vanraden[np.triu_indices(n_samples, k=1)]
    axes1[1, 2].scatter(kinship_prm, kinship_grm, alpha=0.3, s=5, c=PALETTE_A[0])
    axes1[1, 2].plot([-0.5, 1.5], [-0.5, 1.5], color='#333333', linestyle='--', linewidth=1.5)
    axes1[1, 2].set_xlabel('PRM Kinship')
    axes1[1, 2].set_ylabel('GRM Kinship')
    axes1[1, 2].set_title('PRM vs GRM Comparison')
    axes1[1, 2].spines['top'].set_visible(False)
    axes1[1, 2].spines['right'].set_visible(False)
    axes1[1, 2].tick_params(direction='out')
    sns.despine(ax=axes1[1, 2])
    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig1.savefig(os.path.join(output_dir, f"relationship_matrix_overview.{fmt}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  总览图已保存")

    fig3, axes3 = plt.subplots(2, 2, figsize=(14, 12))
    fig3.suptitle('Matrix Distribution Analysis', fontweight='bold', fontfamily='DejaVu Sans')
    plot_distribution(np.diag(prm), axes3[0, 0], 'PRM Diagonal Distribution', 'Self-Relationship', PALETTE_A[0])
    plot_distribution(np.diag(grm_vanraden), axes3[0, 1], 'GRM Diagonal Distribution', 'Self-Relationship', PALETTE_A[1])
    plot_distribution(prm[np.triu_indices(n_samples, k=1)], axes3[1, 0], 'PRM Off-Diagonal Distribution', 'Kinship', PALETTE_A[2])
    plot_distribution(grm_vanraden[np.triu_indices(n_samples, k=1)], axes3[1, 1], 'GRM Off-Diagonal Distribution', 'Kinship', PALETTE_A[3])
    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig3.savefig(os.path.join(output_dir, f"matrix_distribution.{fmt}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  分布分析图已保存")

    print("\n" + "=" * 70)
    print("分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
