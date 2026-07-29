#!/usr/bin/env python3
"""
Step: Compare piggpa-G relationship matrices against HIBLUP's matrices
用法: python compare_matrices.py

比较3对矩阵:
1. piggpa-G PRM vs HIBLUP PA (期望 r=1.0)
2. piggpa-G GRM_VanRaden vs HIBLUP GA (期望 r≈0.99)
3. piggpa-G HRM (默认 0.95*G_adj+0.05*A, 等价HIBLUP HA) vs HIBLUP HA (期望 r>0.999, 对角线均值≈1.0)
"""

import os
import sys
import numpy as np
import pandas as pd

# 路径配置
PIGGPA_DIR = "/public/share/likui/liangcx/bole/HIBLUP_benchmark/new_benchmark/t3"
HIBLUP_DIR = "/public/share/likui/hanyu/testresult/HIBLUP/3/1"
OUTPUT_CSV = "/public/share/likui/liangcx/bole/HIBLUP_benchmark/new_benchmark/t3/matrix_comparison_results.csv"


def read_piggpa_csv(csv_file):
    """读取piggpa-G生成的CSV矩阵(带行列名)"""
    print(f"  Loading piggpa-G matrix: {csv_file}", flush=True)
    df = pd.read_csv(csv_file, index_col=0)
    # 列名可能被读取为字符串, 转为字符串以匹配
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    print(f"    shape={df.shape}, diag_mean={np.mean(np.diag(df.values)):.6f}", flush=True)
    return df


def read_hiblup_txt(txt_file, id_file):
    """读取HIBLUP文本矩阵(无表头, 带独立ID文件)"""
    print(f"  Loading HIBLUP matrix: {txt_file}", flush=True)
    ids = pd.read_csv(id_file, header=None)[0].astype(str).tolist()
    matrix = pd.read_csv(txt_file, sep=r'\s+', header=None).values
    df = pd.DataFrame(matrix, index=ids, columns=ids)
    print(f"    shape={df.shape}, diag_mean={np.mean(np.diag(df.values)):.6f}", flush=True)
    return df


def align_matrices(df1, df2, name1, name2):
    """按共同样本ID对齐两个矩阵"""
    common_ids = [i for i in df1.index if i in df2.index]
    print(f"  Common samples: {len(common_ids)} (piggpa={len(df1)}, hiblup={len(df2)})", flush=True)
    if len(common_ids) == 0:
        raise ValueError(f"No common samples between {name1} and {name2}")
    m1 = df1.loc[common_ids, common_ids].values
    m2 = df2.loc[common_ids, common_ids].values
    return m1, m2


def compute_metrics(m1, m2):
    """计算Pearson r, MSE, Max Diff, 对角线均值"""
    # 排除NaN
    mask = ~(np.isnan(m1) | np.isnan(m2))
    m1_flat = m1[mask]
    m2_flat = m2[mask]

    # Pearson相关系数
    if m1_flat.std() == 0 or m2_flat.std() == 0:
        r = 1.0 if np.allclose(m1_flat, m2_flat) else 0.0
    else:
        r = float(np.corrcoef(m1_flat, m2_flat)[0, 1])

    # MSE
    mse = float(np.mean((m1_flat - m2_flat) ** 2))
    # Max Diff
    max_diff = float(np.max(np.abs(m1_flat - m2_flat)))
    # 对角线均值
    diag1 = np.diag(m1)
    diag2 = np.diag(m2)
    diag_mean1 = float(np.mean(diag1))
    diag_mean2 = float(np.mean(diag2))

    return {
        'r': r,
        'MSE': mse,
        'Max_Diff': max_diff,
        'piggpa_G_diag_mean': diag_mean1,
        'HIBLUP_diag_mean': diag_mean2,
    }


def classify_level(r, mse):
    """根据相关系数和MSE判定一致性等级"""
    # 使用容差判定IDENTICAL, 避免浮点精度问题
    if np.isclose(r, 1.0, atol=1e-12) and mse < 1e-20:
        return 'IDENTICAL'
    elif r >= 0.99 and mse < 1e-4:
        return 'GOOD'
    elif r >= 0.95:
        return 'ACCEPTABLE'
    else:
        return 'DIVERGENT'


def main():
    print("=" * 70, flush=True)
    print("piggpa-G vs HIBLUP 关系矩阵比对", flush=True)
    print("=" * 70, flush=True)

    # 读取HIBLUP矩阵
    print("\n[1/3] 读取HIBLUP矩阵...", flush=True)
    hiblup_pa = read_hiblup_txt(os.path.join(HIBLUP_DIR, "PA_txt.txt"),
                                os.path.join(HIBLUP_DIR, "PA_txt.id.txt"))
    hiblup_ga = read_hiblup_txt(os.path.join(HIBLUP_DIR, "GA_txt.txt"),
                                os.path.join(HIBLUP_DIR, "GA_txt.id.txt"))
    hiblup_ha = read_hiblup_txt(os.path.join(HIBLUP_DIR, "HA_txt.txt"),
                                os.path.join(HIBLUP_DIR, "HA_txt.id.txt"))

    # 读取piggpa-G矩阵
    print("\n[2/3] 读取piggpa-G矩阵...", flush=True)
    piggpa_prm = read_piggpa_csv(os.path.join(PIGGPA_DIR, "PRM.csv"))
    piggpa_grm = read_piggpa_csv(os.path.join(PIGGPA_DIR, "GRM_VanRaden.csv"))
    piggpa_hrm = read_piggpa_csv(os.path.join(PIGGPA_DIR, "HRM.csv"))

    # 3对比对
    print("\n[3/3] 执行3对比对...", flush=True)
    pairs = [
        ('PRM_vs_PA', piggpa_prm, hiblup_pa, 'PRM', 'PA'),
        ('GRM_VanRaden_vs_GA', piggpa_grm, hiblup_ga, 'GRM_VanRaden', 'GA'),
        ('HRM_vs_HA', piggpa_hrm, hiblup_ha, 'HRM', 'HA'),
    ]

    results = []
    for pair_name, df_p, df_h, name_p, name_h in pairs:
        print(f"\n--- {pair_name}: piggpa-G {name_p} vs HIBLUP {name_h} ---", flush=True)
        m_p, m_h = align_matrices(df_p, df_h, name_p, name_h)
        metrics = compute_metrics(m_p, m_h)
        level = classify_level(metrics['r'], metrics['MSE'])
        metrics['Matrix_Pair'] = pair_name
        metrics['Level'] = level
        results.append(metrics)
        print(f"  r={metrics['r']:.6f}, MSE={metrics['MSE']:.6e}, "
              f"Max_Diff={metrics['Max_Diff']:.6f}, "
              f"piggpa_diag={metrics['piggpa_G_diag_mean']:.6f}, "
              f"hiblup_diag={metrics['HIBLUP_diag_mean']:.6f}, "
              f"Level={level}", flush=True)

    # 写入CSV
    df_results = pd.DataFrame(results)[
        ['Matrix_Pair', 'Level', 'r', 'MSE', 'Max_Diff',
         'piggpa_G_diag_mean', 'HIBLUP_diag_mean']
    ]
    df_results.to_csv(OUTPUT_CSV, index=False)
    print(f"\n结果已保存至: {OUTPUT_CSV}", flush=True)
    print("\n" + "=" * 70, flush=True)
    print("比对结果汇总:", flush=True)
    print("=" * 70, flush=True)
    print(df_results.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
