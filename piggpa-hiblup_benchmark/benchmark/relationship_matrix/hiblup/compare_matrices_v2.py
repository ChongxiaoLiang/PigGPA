#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
关系矩阵比对脚本
比较Python和HIBLUP构建的三种关系矩阵的一致性
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

def read_hiblup_txt_matrix(txt_file, id_file):
    """读取HIBLUP转换的文本矩阵"""
    ids = pd.read_csv(id_file, header=None)[0].tolist()
    matrix = pd.read_csv(txt_file, sep=r'\s+', header=None).values
    return matrix, ids

def read_python_csv_matrix(csv_file):
    """读取Python生成的CSV矩阵（带行列名）"""
    df = pd.read_csv(csv_file, index_col=0)
    matrix = df.values
    return matrix

def compare_matrices(matrix1, matrix2, name1, name2):
    """比较两个矩阵的一致性"""
    
    if matrix1.shape != matrix2.shape:
        return {
            'status': 'ERROR',
            'message': f'矩阵维度不匹配: {matrix1.shape} vs {matrix2.shape}',
            'correlation': np.nan,
            'mse': np.nan,
            'max_diff': np.nan,
            'mean_diff': np.nan,
            'identical': False
        }
    
    mask = ~(np.isnan(matrix1) | np.isnan(matrix2))
    
    m1_flat = matrix1[mask]
    m2_flat = matrix2[mask]
    
    corr = np.corrcoef(m1_flat, m2_flat)[0, 1]
    mse = np.mean((m1_flat - m2_flat) ** 2)
    max_diff = np.max(np.abs(m1_flat - m2_flat))
    mean_diff = np.mean(np.abs(m1_flat - m2_flat))
    rmse = np.sqrt(mse)
    
    diag1 = np.diag(matrix1)
    diag2 = np.diag(matrix2)
    diag_corr = np.corrcoef(diag1, diag2)[0, 1] if np.std(diag1) > 0 and np.std(diag2) > 0 else 1.0
    diag_mse = np.mean((diag1 - diag2) ** 2)
    
    off_diag_mask = ~np.eye(matrix1.shape[0], dtype=bool)
    off_diag1 = matrix1[off_diag_mask]
    off_diag2 = matrix2[off_diag_mask]
    off_diag_corr = np.corrcoef(off_diag1, off_diag2)[0, 1] if np.std(off_diag1) > 0 and np.std(off_diag2) > 0 else 1.0
    off_diag_mse = np.mean((off_diag1 - off_diag2) ** 2)
    
    identical = np.allclose(matrix1, matrix2, rtol=1e-5, atol=1e-8)
    
    if identical:
        status = 'IDENTICAL'
    elif corr > 0.9999:
        status = 'EXCELLENT'
    elif corr > 0.999:
        status = 'VERY_GOOD'
    elif corr > 0.99:
        status = 'GOOD'
    elif corr > 0.9:
        status = 'ACCEPTABLE'
    else:
        status = 'POOR'
    
    return {
        'status': status,
        'correlation': corr,
        'mse': mse,
        'rmse': rmse,
        'max_diff': max_diff,
        'mean_diff': mean_diff,
        'diag_correlation': diag_corr,
        'diag_mse': diag_mse,
        'off_diag_correlation': off_diag_corr,
        'off_diag_mse': off_diag_mse,
        'matrix1_mean': np.mean(matrix1),
        'matrix1_std': np.std(matrix1),
        'matrix1_diag_mean': np.mean(diag1),
        'matrix1_off_diag_mean': np.mean(off_diag1),
        'matrix2_mean': np.mean(matrix2),
        'matrix2_std': np.std(matrix2),
        'matrix2_diag_mean': np.mean(diag2),
        'matrix2_off_diag_mean': np.mean(off_diag2),
        'identical': identical
    }

def create_comparison_plot(matrix1, matrix2, name1, name2, output_file):
    """创建矩阵比对可视化图"""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    ax1 = axes[0, 0]
    im1 = ax1.imshow(matrix1, cmap='RdBu_r', aspect='auto')
    ax1.set_title(f'{name1}')
    plt.colorbar(im1, ax=ax1)
    
    ax2 = axes[0, 1]
    im2 = ax2.imshow(matrix2, cmap='RdBu_r', aspect='auto')
    ax2.set_title(f'{name2}')
    plt.colorbar(im2, ax=ax2)
    
    ax3 = axes[0, 2]
    diff = matrix1 - matrix2
    im3 = ax3.imshow(diff, cmap='RdBu_r', aspect='auto')
    ax3.set_title('Difference (Python - HIBLUP)')
    plt.colorbar(im3, ax=ax3)
    
    ax4 = axes[1, 0]
    mask = ~(np.isnan(matrix1) | np.isnan(matrix2))
    m1_flat = matrix1[mask]
    m2_flat = matrix2[mask]
    ax4.scatter(m1_flat, m2_flat, alpha=0.1, s=1)
    min_val = min(m1_flat.min(), m2_flat.min())
    max_val = max(m1_flat.max(), m2_flat.max())
    ax4.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
    corr = np.corrcoef(m1_flat, m2_flat)[0, 1]
    ax4.set_xlabel(name1)
    ax4.set_ylabel(name2)
    ax4.set_title(f'Correlation: {corr:.6f}')
    ax4.grid(True, alpha=0.3)
    
    ax5 = axes[1, 1]
    ax5.hist(diff[mask], bins=50, edgecolor='black', alpha=0.7)
    ax5.axvline(x=0, color='r', linestyle='--', lw=2)
    ax5.set_xlabel('Difference')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Difference Distribution')
    ax5.grid(True, alpha=0.3)
    
    ax6 = axes[1, 2]
    diag_diff = np.diag(matrix1) - np.diag(matrix2)
    off_diag_diff = matrix1[~np.eye(matrix1.shape[0], dtype=bool)] - matrix2[~np.eye(matrix2.shape[0], dtype=bool)]
    ax6.hist(diag_diff, bins=30, alpha=0.7, label='Diagonal', edgecolor='black')
    ax6.hist(off_diag_diff, bins=30, alpha=0.7, label='Off-diagonal', edgecolor='black')
    ax6.axvline(x=0, color='r', linestyle='--', lw=2)
    ax6.set_xlabel('Difference')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Diagonal vs Off-diagonal Differences')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', dpi=300)
    plt.close()

def main():
    print("=" * 80)
    print("关系矩阵一致性比对分析")
    print("=" * 80)
    
    python_dir = "/public/share/likui/hanyu/testdata/In-silico-data/t3"
    hiblup_dir = "/public/share/likui/hanyu/testresult/HIBLUP/3/1"
    output_dir = "/public/share/likui/hanyu/testresult/HIBLUP/16"
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\nPython目录: {python_dir}")
    print(f"HIBLUP目录: {hiblup_dir}")
    print(f"输出目录: {output_dir}")
    
    matrix_pairs = [
        ('PRM', 'PA', '系谱关系矩阵', 'PRM.csv', 'PA_txt.txt', 'PA_txt.id.txt'),
        ('GRM_VanRaden', 'GA', '基因组关系矩阵(VanRaden)', 'GRM_VanRaden.csv', 'GA_txt.txt', 'GA_txt.id.txt'),
        ('HRM', 'HA', '混合关系矩阵', 'HRM.csv', 'HA_txt.txt', 'HA_txt.id.txt'),
    ]
    
    results = []
    
    print("\n" + "=" * 80)
    print("开始比对...")
    print("=" * 80)
    
    for py_name, hiblup_name, cn_name, csv_file, txt_file, id_file in matrix_pairs:
        print(f"\n{'='*60}")
        print(f"比对: {cn_name}")
        print(f"Python矩阵: {py_name}")
        print(f"HIBLUP矩阵: {hiblup_name}")
        print(f"{'='*60}")
        
        csv_path = os.path.join(python_dir, csv_file)
        txt_path = os.path.join(hiblup_dir, txt_file)
        id_path = os.path.join(hiblup_dir, id_file)
        
        if not os.path.exists(csv_path):
            print(f"  警告: Python文件不存在 - {csv_path}")
            continue
        
        if not os.path.exists(txt_path):
            print(f"  警告: HIBLUP文件不存在 - {txt_path}")
            continue
        
        print(f"  读取Python CSV文件: {csv_file}")
        matrix_python = read_python_csv_matrix(csv_path)
        print(f"    维度: {matrix_python.shape}")
        print(f"    均值: {np.mean(matrix_python):.6f}")
        print(f"    标准差: {np.std(matrix_python):.6f}")
        print(f"    对角线均值: {np.mean(np.diag(matrix_python)):.6f}")
        print(f"    非对角线均值: {np.mean(matrix_python[~np.eye(matrix_python.shape[0], dtype=bool)]):.6f}")
        
        print(f"  读取HIBLUP文本文件: {txt_file}")
        matrix_hiblup, ids = read_hiblup_txt_matrix(txt_path, id_path)
        print(f"    维度: {matrix_hiblup.shape}")
        print(f"    样本数: {len(ids)}")
        print(f"    均值: {np.mean(matrix_hiblup):.6f}")
        print(f"    标准差: {np.std(matrix_hiblup):.6f}")
        print(f"    对角线均值: {np.mean(np.diag(matrix_hiblup)):.6f}")
        print(f"    非对角线均值: {np.mean(matrix_hiblup[~np.eye(matrix_hiblup.shape[0], dtype=bool)]):.6f}")
        
        print(f"  开始比对...")
        result = compare_matrices(matrix_python, matrix_hiblup, py_name, hiblup_name)
        
        print(f"\n  比对结果:")
        print(f"    状态: {result['status']}")
        print(f"    完全一致: {result['identical']}")
        print(f"    整体相关系数: {result['correlation']:.10f}")
        print(f"    MSE: {result['mse']:.10e}")
        print(f"    RMSE: {result['rmse']:.10e}")
        print(f"    最大差异: {result['max_diff']:.10e}")
        print(f"    平均差异: {result['mean_diff']:.10e}")
        print(f"    对角线相关系数: {result['diag_correlation']:.10f}")
        print(f"    对角线MSE: {result['diag_mse']:.10e}")
        print(f"    非对角线相关系数: {result['off_diag_correlation']:.10f}")
        print(f"    非对角线MSE: {result['off_diag_mse']:.10e}")
        
        plot_file = os.path.join(output_dir, f"{py_name}_vs_{hiblup_name}_comparison.pdf")
        print(f"  生成可视化图表: {plot_file}")
        create_comparison_plot(matrix_python, matrix_hiblup, py_name, hiblup_name, plot_file)
        
        results.append({
            'python_name': py_name,
            'hiblup_name': hiblup_name,
            'chinese_name': cn_name,
            **result
        })
    
    print("\n" + "=" * 80)
    print("比对结果汇总")
    print("=" * 80)
    
    summary_df = pd.DataFrame(results)
    print("\n" + summary_df[['python_name', 'hiblup_name', 'chinese_name', 'status', 'identical', 
                           'correlation', 'mse', 'max_diff']].to_string(index=False))
    
    output_file = os.path.join(output_dir, "matrix_comparison_report.txt")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("关系矩阵一致性比对报告\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("一、数据来源\n")
        f.write("-" * 40 + "\n")
        f.write(f"Python目录: {python_dir}\n")
        f.write(f"HIBLUP目录: {hiblup_dir}\n\n")
        
        f.write("二、比对结果汇总\n")
        f.write("-" * 40 + "\n")
        f.write(summary_df[['python_name', 'hiblup_name', 'chinese_name', 'status', 'identical', 
                           'correlation', 'mse', 'max_diff']].to_string(index=False))
        f.write("\n\n")
        
        f.write("三、详细分析\n")
        f.write("-" * 40 + "\n")
        
        for result in results:
            f.write(f"\n{result['chinese_name']} ({result['python_name']} vs {result['hiblup_name']}):\n")
            f.write(f"  状态: {result['status']}\n")
            f.write(f"  完全一致: {result['identical']}\n")
            f.write(f"  整体相关系数: {result['correlation']:.10f}\n")
            f.write(f"  MSE: {result['mse']:.10e}\n")
            f.write(f"  RMSE: {result['rmse']:.10e}\n")
            f.write(f"  最大差异: {result['max_diff']:.10e}\n")
            f.write(f"  平均差异: {result['mean_diff']:.10e}\n")
            f.write(f"  对角线相关系数: {result['diag_correlation']:.10f}\n")
            f.write(f"  对角线MSE: {result['diag_mse']:.10e}\n")
            f.write(f"  非对角线相关系数: {result['off_diag_correlation']:.10f}\n")
            f.write(f"  非对角线MSE: {result['off_diag_mse']:.10e}\n")
            f.write(f"\n  Python矩阵统计:\n")
            f.write(f"    均值: {result['matrix1_mean']:.6f}\n")
            f.write(f"    标准差: {result['matrix1_std']:.6f}\n")
            f.write(f"    对角线均值: {result['matrix1_diag_mean']:.6f}\n")
            f.write(f"    非对角线均值: {result['matrix1_off_diag_mean']:.6f}\n")
            f.write(f"\n  HIBLUP矩阵统计:\n")
            f.write(f"    均值: {result['matrix2_mean']:.6f}\n")
            f.write(f"    标准差: {result['matrix2_std']:.6f}\n")
            f.write(f"    对角线均值: {result['matrix2_diag_mean']:.6f}\n")
            f.write(f"    非对角线均值: {result['matrix2_off_diag_mean']:.6f}\n")
        
        f.write("\n四、结论\n")
        f.write("-" * 40 + "\n")
        
        all_identical = all(r['identical'] for r in results)
        all_status = [r['status'] for r in results]
        correlations = [r['correlation'] for r in results]
        
        if all_identical:
            f.write("所有关系矩阵完全一致!\n")
            f.write("Python和HIBLUP构建的关系矩阵在数值精度范围内完全相同。\n")
        elif all(s in ['IDENTICAL', 'EXCELLENT'] for s in all_status):
            f.write("所有关系矩阵高度一致!\n")
            f.write("Python和HIBLUP构建的关系矩阵具有极高的相关性。\n")
        elif all(s in ['IDENTICAL', 'EXCELLENT', 'VERY_GOOD'] for s in all_status):
            f.write("所有关系矩阵一致性很好!\n")
            f.write("Python和HIBLUP构建的关系矩阵具有很高的相关性。\n")
        else:
            f.write("部分关系矩阵存在差异，需要检查!\n")
        
        f.write(f"\n统计摘要:\n")
        f.write(f"  平均相关系数: {np.mean(correlations):.10f}\n")
        f.write(f"  最低相关系数: {np.min(correlations):.10f}\n")
        f.write(f"  最高相关系数: {np.max(correlations):.10f}\n")
        f.write(f"  完全一致的矩阵数: {sum(r['identical'] for r in results)}/{len(results)}\n")
        
        f.write("\n五、方法说明\n")
        f.write("-" * 40 + "\n")
        f.write("比对方法:\n")
        f.write("  1. 相关系数: 衡量两个矩阵的整体线性相关性\n")
        f.write("  2. MSE (均方误差): 衡量两个矩阵的平均偏差\n")
        f.write("  3. RMSE (均方根误差): MSE的平方根，与原始数据同单位\n")
        f.write("  4. 最大差异: 两个矩阵元素的最大绝对差异\n")
        f.write("  5. 平均差异: 两个矩阵元素的平均绝对差异\n")
        f.write("\n判断标准:\n")
        f.write("  - IDENTICAL: 完全一致 (数值精度范围内)\n")
        f.write("  - EXCELLENT: 相关系数 > 0.9999\n")
        f.write("  - VERY_GOOD: 相关系数 > 0.999\n")
        f.write("  - GOOD: 相关系数 > 0.99\n")
        f.write("  - ACCEPTABLE: 相关系数 > 0.9\n")
        f.write("  - POOR: 相关系数 <= 0.9\n")
    
    print(f"\n报告已保存: {output_file}")
    
    csv_output = os.path.join(output_dir, "matrix_comparison_results.csv")
    summary_df.to_csv(csv_output, index=False)
    print(f"结果CSV已保存: {csv_output}")
    
    print("\n" + "=" * 80)
    print("比对完成!")
    print("=" * 80)
    
    if all(r['identical'] for r in results):
        print("\n结论: 所有关系矩阵完全一致!")
    elif all(r['correlation'] > 0.999 for r in results):
        print("\n结论: 所有关系矩阵高度一致!")
    else:
        print("\n结论: 部分关系矩阵存在差异，请查看详细报告。")


if __name__ == "__main__":
    main()
