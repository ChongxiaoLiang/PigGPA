#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HIBLUP LD矩阵二进制文件转换脚本
将ld_result.bin转换为文本格式
"""

import numpy as np
import pandas as pd
import os

def read_ld_binary(bin_file, info_file, output_file):
    """
    读取HIBLUP生成的LD矩阵二进制文件并转换为文本格式
    
    参数:
        bin_file: 二进制文件路径
        info_file: SNP信息文件路径
        output_file: 输出文本文件路径
    """
    
    print("=" * 70)
    print("HIBLUP LD矩阵转换工具")
    print("=" * 70)
    
    print("\n[第一步] 读取SNP信息...")
    snp_info = pd.read_csv(info_file, sep='\t')
    n_snp = len(snp_info)
    print(f"  SNP数量: {n_snp}")
    print(f"  前5个SNP:")
    print(snp_info.head())
    
    print("\n[第二步] 读取二进制矩阵...")
    file_size = os.path.getsize(bin_file)
    print(f"  文件大小: {file_size} bytes")
    
    expected_size_full = n_snp * n_snp * 4
    expected_size_triangle = (n_snp * (n_snp + 1) // 2) * 4
    
    print(f"  预期大小(完整矩阵): {expected_size_full} bytes")
    print(f"  预期大小(三角矩阵): {expected_size_triangle} bytes")
    
    with open(bin_file, 'rb') as f:
        data = f.read()
    
    print(f"  实际读取: {len(data)} bytes")
    
    if len(data) == expected_size_full:
        print("  检测到完整矩阵格式")
        ld_matrix = np.frombuffer(data, dtype=np.float32).reshape(n_snp, n_snp)
    elif len(data) == expected_size_triangle:
        print("  检测到三角矩阵格式")
        ld_matrix = np.zeros((n_snp, n_snp), dtype=np.float32)
        idx = 0
        for i in range(n_snp):
            for j in range(i, n_snp):
                ld_matrix[i, j] = np.frombuffer(data[idx*4:(idx+1)*4], dtype=np.float32)[0]
                ld_matrix[j, i] = ld_matrix[i, j]
                idx += 1
    else:
        print(f"  警告: 文件大小不匹配，尝试作为完整矩阵读取...")
        try:
            ld_matrix = np.frombuffer(data, dtype=np.float32).reshape(n_snp, n_snp)
        except:
            print("  错误: 无法解析二进制文件")
            return
    
    print(f"  矩阵维度: {ld_matrix.shape}")
    print(f"  对角线均值: {np.mean(np.diag(ld_matrix)):.6f}")
    print(f"  非对角线均值: {np.mean(ld_matrix[~np.eye(n_snp, dtype=bool)]):.6f}")
    print(f"  矩阵范围: [{np.min(ld_matrix):.6f}, {np.max(ld_matrix):.6f}]")
    
    print("\n[第三步] 保存为文本格式...")
    
    output_triangle = output_file.replace('.txt', '_triangle.txt')
    print(f"  保存三角矩阵格式到: {output_triangle}")
    with open(output_triangle, 'w') as f:
        f.write("SNP_i\tSNP_j\tLD_r\n")
        for i in range(n_snp):
            for j in range(i, n_snp):
                snp_i = snp_info.iloc[i]['SNP']
                snp_j = snp_info.iloc[j]['SNP']
                ld_value = ld_matrix[i, j]
                f.write(f"{snp_i}\t{snp_j}\t{ld_value:.6f}\n")
    
    print(f"  三角矩阵已保存")
    
    print("\n[第四步] 生成统计摘要...")
    ld_values = ld_matrix[~np.eye(n_snp, dtype=bool)]
    
    stats_file = output_file.replace('.txt', '_stats.txt')
    with open(stats_file, 'w') as f:
        f.write("LD矩阵统计摘要\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"SNP总数: {n_snp}\n")
        f.write(f"LD值对数: {len(ld_values)}\n\n")
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
        high_ld = np.sum(ld_values > 0.8)
        f.write(f"  LD > 0.8 的SNP对数: {high_ld}\n")
        f.write(f"  占比: {high_ld/len(ld_values)*100:.2f}%\n")
    
    print(f"  统计摘要已保存: {stats_file}")
    
    print("\n[第五步] 保存前100个SNP的LD矩阵示例...")
    sample_size = min(100, n_snp)
    sample_file = output_file.replace('.txt', '_sample100.txt')
    
    with open(sample_file, 'w') as f:
        header = ['SNP'] + snp_info.iloc[:sample_size]['SNP'].tolist()
        f.write('\t'.join(header) + '\n')
        
        for i in range(sample_size):
            row = [snp_info.iloc[i]['SNP']]
            for j in range(sample_size):
                row.append(f"{ld_matrix[i, j]:.6f}")
            f.write('\t'.join(row) + '\n')
    
    print(f"  示例矩阵已保存: {sample_file}")
    
    print("\n" + "=" * 70)
    print("转换完成!")
    print("=" * 70)
    print(f"\n输出文件:")
    print(f"  1. 三角矩阵: {output_triangle}")
    print(f"  2. 统计摘要: {stats_file}")
    print(f"  3. 示例矩阵: {sample_file}")


if __name__ == "__main__":
    bin_file = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result.bin"
    info_file = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result.info"
    output_file = "/public/share/likui/hanyu/testresult/HIBLUP/14/ld_result.txt"
    
    read_ld_binary(bin_file, info_file, output_file)
