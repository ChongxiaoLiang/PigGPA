#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HIBLUP LD矩阵二进制文件转换脚本（窗口格式）
将ld_result.bin转换为文本格式
"""

import numpy as np
import pandas as pd
import os

def read_ld_binary_windowed(bin_file, info_file, output_prefix):
    """
    读取HIBLUP生成的按窗口存储的LD矩阵二进制文件并转换为文本格式
    
    参数:
        bin_file: 二进制文件路径
        info_file: SNP信息文件路径
        output_prefix: 输出文件前缀
    """
    
    print("=" * 70)
    print("HIBLUP LD矩阵转换工具（窗口格式）")
    print("=" * 70)
    
    print("\n[第一步] 读取SNP信息...")
    snp_info = pd.read_csv(info_file, sep='\t')
    n_snp = len(snp_info)
    print(f"  总SNP数量: {n_snp}")
    
    print("\n[第二步] 分析窗口结构...")
    windows = snp_info.groupby('WindIndex')
    n_windows = len(windows)
    print(f"  窗口数量: {n_windows}")
    
    window_info = []
    for wind_idx, wind_data in windows:
        n_snps_in_window = len(wind_data)
        n_ld_values = n_snps_in_window * (n_snps_in_window + 1) // 2
        window_info.append({
            'window_index': wind_idx,
            'n_snps': n_snps_in_window,
            'n_ld_values': n_ld_values,
            'snp_indices': wind_data.index.tolist()
        })
        print(f"  窗口 {wind_idx}: {n_snps_in_window} SNPs, {n_ld_values} LD值")
    
    total_ld_values = sum([w['n_ld_values'] for w in window_info])
    print(f"\n  总LD值数量: {total_ld_values}")
    
    print("\n[第三步] 读取二进制矩阵...")
    file_size = os.path.getsize(bin_file)
    print(f"  文件大小: {file_size} bytes")
    print(f"  预期大小: {total_ld_values * 4} bytes")
    
    with open(bin_file, 'rb') as f:
        data = np.frombuffer(f.read(), dtype=np.float32)
    
    print(f"  实际读取: {len(data)} 个浮点数")
    
    if len(data) != total_ld_values:
        print(f"  警告: 数据大小不匹配!")
        print(f"  尝试按实际数据解析...")
    
    print("\n[第四步] 解析LD矩阵...")
    ld_results = []
    offset = 0
    
    for w_info in window_info:
        wind_idx = w_info['window_index']
        n_snps = w_info['n_snps']
        snp_indices = w_info['snp_indices']
        n_ld = w_info['n_ld_values']
        
        ld_values = data[offset:offset+n_ld]
        offset += n_ld
        
        ld_matrix = np.zeros((n_snps, n_snps), dtype=np.float32)
        idx = 0
        for i in range(n_snps):
            for j in range(i, n_snps):
                ld_matrix[i, j] = ld_values[idx]
                ld_matrix[j, i] = ld_matrix[i, j]
                idx += 1
        
        for i in range(n_snps):
            for j in range(i, n_snps):
                snp_i_idx = snp_indices[i]
                snp_j_idx = snp_indices[j]
                snp_i = snp_info.iloc[snp_i_idx]['SNP']
                snp_j = snp_info.iloc[snp_j_idx]['SNP']
                ld_value = ld_matrix[i, j]
                
                ld_results.append({
                    'Window': wind_idx,
                    'SNP_i': snp_i,
                    'SNP_j': snp_j,
                    'LD_r': ld_value,
                    'Pos_i': snp_info.iloc[snp_i_idx]['Pos'],
                    'Pos_j': snp_info.iloc[snp_j_idx]['Pos'],
                    'Distance': abs(snp_info.iloc[snp_i_idx]['Pos'] - snp_info.iloc[snp_j_idx]['Pos'])
                })
    
    ld_df = pd.DataFrame(ld_results)
    print(f"  解析完成，共 {len(ld_df)} 个LD值对")
    
    print("\n[第五步] 保存结果...")
    
    output_all = f"{output_prefix}_all.txt"
    print(f"  保存所有LD值到: {output_all}")
    ld_df.to_csv(output_all, sep='\t', index=False)
    
    output_summary = f"{output_prefix}_summary.txt"
    print(f"  保存统计摘要到: {output_summary}")
    with open(output_summary, 'w') as f:
        f.write("LD矩阵统计摘要\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"总SNP数: {n_snp}\n")
        f.write(f"窗口数: {n_windows}\n")
        f.write(f"LD值对数: {len(ld_df)}\n\n")
        
        ld_values = ld_df['LD_r'].values
        f.write("LD值分布:\n")
        f.write(f"  最小值: {np.min(ld_values):.6f}\n")
        f.write(f"  最大值: {np.max(ld_values):.6f}\n")
        f.write(f"  均值: {np.mean(ld_values):.6f}\n")
        f.write(f"  中位数: {np.median(ld_values):.6f}\n")
        f.write(f"  标准差: {np.std(ld_values):.6f}\n\n")
        
        f.write("LD值分位数:\n")
        for q in [0.25, 0.5, 0.75, 0.9, 0.95, 0.99]:
            f.write(f"  {int(q*100)}%分位数: {np.percentile(ld_values, q*100):.6f}\n")
        
        f.write("\n高LD SNP对统计:\n")
        for threshold in [0.2, 0.5, 0.8, 0.9]:
            count = np.sum(ld_values >= threshold)
            f.write(f"  LD >= {threshold}: {count} ({count/len(ld_values)*100:.2f}%)\n")
        
        f.write("\n距离分布:\n")
        distances = ld_df['Distance'].values
        f.write(f"  最小距离: {np.min(distances):.0f} bp\n")
        f.write(f"  最大距离: {np.max(distances):.0f} bp\n")
        f.write(f"  平均距离: {np.mean(distances):.0f} bp\n")
    
    output_high_ld = f"{output_prefix}_high_ld.txt"
    print(f"  保存高LD值对到: {output_high_ld}")
    high_ld_df = ld_df[ld_df['LD_r'] >= 0.8].sort_values('LD_r', ascending=False)
    high_ld_df.to_csv(output_high_ld, sep='\t', index=False)
    
    print("\n[第六步] 生成窗口统计...")
    window_stats = []
    for w_info in window_info:
        wind_idx = w_info['window_index']
        wind_ld = ld_df[ld_df['Window'] == wind_idx]
        
        window_stats.append({
            'Window': wind_idx,
            'N_SNPs': w_info['n_snps'],
            'N_LD_pairs': len(wind_ld),
            'Mean_LD': wind_ld['LD_r'].mean(),
            'Max_LD': wind_ld['LD_r'].max(),
            'Median_LD': wind_ld['LD_r'].median(),
            'High_LD_count': len(wind_ld[wind_ld['LD_r'] >= 0.8])
        })
    
    window_stats_df = pd.DataFrame(window_stats)
    output_window_stats = f"{output_prefix}_window_stats.txt"
    print(f"  保存窗口统计到: {output_window_stats}")
    window_stats_df.to_csv(output_window_stats, sep='\t', index=False)
    
    print("\n" + "=" * 70)
    print("转换完成!")
    print("=" * 70)
    print(f"\n输出文件:")
    print(f"  1. 所有LD值: {output_all}")
    print(f"  2. 统计摘要: {output_summary}")
    print(f"  3. 高LD值对: {output_high_ld}")
    print(f"  4. 窗口统计: {output_window_stats}")
    print(f"\n数据预览:")
    print(ld_df.head(20).to_string(index=False))


if __name__ == "__main__":
    bin_file = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result.bin"
    info_file = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result.info"
    output_prefix = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result"
    
    read_ld_binary_windowed(bin_file, info_file, output_prefix)
