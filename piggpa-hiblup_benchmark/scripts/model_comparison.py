#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模型对比脚本 - 与HIBLUP输出格式一致
比较五种模型的预测性能

用法:
  python model_comparison.py \
    --result-dir /path/to/results \
    --pheno /path/to/phenotypes.txt \
    --out /path/to/output_dir \
    [--pred-id-file /path/to/pred_ids.txt] \
    [--pheno-col Phenotype] \
    [--id-col ID] \
    [--bv-col BreedingValue]
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

PALETTE_A = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7']
PALETTE_B = ['#1A9899', '#EC8528', '#EAC94D', '#FF9DA7', '#4E79A7', '#E15759', '#59A14F']
PALETTE_C = ['#B07AA1', '#EDC948', '#76B7B2', '#4E79A7', '#1A9899', '#FF9DA7', '#F28E2B', '#9C755F']
PALETTE_D = ['#A0CBE8', '#F1CE63', '#8CD17D', '#FFBE7D', '#B6992D', '#499894']
PALETTE_E = ['#d73221', '#e35235', '#e48070', '#fcb777', '#fde699', '#fef4ae', '#d2edf2', '#6491c1', '#4573b4']
DEFAULT_PALETTE = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948', '#B07AA1', '#FF9DA7', '#1A9899', '#EC8528', '#EAC94D', '#9C755F']
WARM_COOL_CMAP = LinearSegmentedColormap.from_list('warm_cool', PALETTE_E[::-1], N=256)


def parse_args():
    parser = argparse.ArgumentParser(
        description='五种模型对比分析 (与HIBLUP格式一致)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--result-dir', required=True, help='模型结果根目录 (包含LM/BLUP/PBLUP/GBLUP/SSBLUP子目录)')
    parser.add_argument('--pheno', required=True, help='表型文件路径 (需包含ID列、表型列和育种值列)')
    parser.add_argument('--out', required=True, help='输出目录')
    parser.add_argument('--pred-id-file', default=None, help='预测样本ID文件，每行一个ID (默认: 从表型文件中读取全部有育种值的样本)')
    parser.add_argument('--pheno-col', default='Phenotype', help='表型列名 (默认: Phenotype)')
    parser.add_argument('--id-col', default='ID', help='样本ID列名 (默认: ID)')
    parser.add_argument('--bv-col', default='BreedingValue', help='真实育种值列名 (默认: BreedingValue)')
    parser.add_argument('--font-size', type=int, default=12, help='Font size for figures (default: 12)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for figure output (default: 300)')
    return parser.parse_args()


def parse_vars(vars_file):
    if not os.path.exists(vars_file):
        return None, None, None
    df = pd.read_csv(vars_file, sep='\t', comment='#')
    var_g = None
    h2 = None
    rm_type = None
    for _, row in df.iterrows():
        if row['Item'] in ['GA', 'PA', 'HA']:
            var_g = row['Var']
            h2 = row['h2']
            rm_type = row['Item']
    return var_g, h2, rm_type


def main():
    args = parse_args()
    plt.rcParams.update({
        'font.size': args.font_size,
        'font.family': 'DejaVu Sans',
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
    })
    os.makedirs(args.out, exist_ok=True)

    print("=" * 70)
    print("五种模型对比分析 (与HIBLUP格式一致)")
    print("=" * 70)
    print(f"  结果目录: {args.result_dir}")
    print(f"  表型文件: {args.pheno}")
    print(f"  输出目录: {args.out}")

    pheno_df = pd.read_csv(args.pheno, sep='\t')
    pheno_df[args.id_col] = pheno_df[args.id_col].astype(str)

    if args.pred_id_file:
        pred_id_list = []
        with open(args.pred_id_file, 'r') as f:
            for line in f:
                pid = line.strip()
                if pid:
                    pred_id_list.append(pid)
        pred_sub = pheno_df[pheno_df[args.id_col].isin(pred_id_list)].set_index(args.id_col).sort_index()
    else:
        pred_sub = pheno_df.dropna(subset=[args.bv_col]).set_index(args.id_col).sort_index()

    has_bv = args.bv_col in pred_sub.columns
    if has_bv:
        true_bv = pred_sub[args.bv_col].values.astype(float)
    true_pheno = pred_sub[args.pheno_col].values.astype(float) if args.pheno_col in pred_sub.columns else None
    print(f"  验证样本: {len(pred_sub)}")

    models_config = {
        'LM': {'vars': 'lm_train.vars', 'pred': 'lm_pred.bv', 'dir': 'LM'},
        'BLUP': {'vars': 'blup_train.vars', 'pred': 'blup_pred.bv', 'dir': 'BLUP'},
        'PBLUP': {'vars': 'pblup_train.vars', 'pred': None, 'dir': 'PBLUP'},
        'GBLUP': {'vars': 'gblup_train.vars', 'pred': 'gblup_pred.bv', 'dir': 'GBLUP'},
        'SSBLUP': {'vars': 'ssblup_train.vars', 'pred': 'ssblup_pred.bv', 'dir': 'SSBLUP'},
    }

    rm_desc = {'GA': 'GA (基因组加性)', 'PA': 'PA (系谱加性)', 'HA': 'HA (混合加性)'}
    # 模型与关系矩阵的对应关系（用于CSV的relationship_matrix列及vars文件缺失时的回退）
    # 注意: piggpa-G的BLUP使用GA（基因组关系矩阵），与GBLUP实现相同；
    #       而HIBLUP的BLUP使用PA（系谱关系矩阵），两者不可直接数值对比。
    model_to_rm = {
        'LM': 'GA',
        'BLUP': 'GA',
        'PBLUP': 'PA',
        'GBLUP': 'GA',
        'SSBLUP': 'HA',
    }

    results = []
    for model_name, config in models_config.items():
        model_dir = os.path.join(args.result_dir, config['dir'])
        vars_file = os.path.join(model_dir, config['vars'])
        var_g, h2, rm_type = parse_vars(vars_file)

        bv_corr = None
        pheno_corr = None
        mse = None

        if config['pred'] is not None:
            pred_file = os.path.join(model_dir, config['pred'])
            if os.path.exists(pred_file):
                pred_df = pd.read_csv(pred_file, sep='\t')
                pred_df['id'] = pred_df['id'].astype(str)
                common_ids = [pid for pid in pred_sub.index if pid in pred_df['id'].values]
                if common_ids:
                    pred_vals = pred_df.set_index('id').loc[common_ids, 'add_a1'].values.astype(float)
                    if has_bv:
                        true_vals = pred_sub.loc[common_ids, args.bv_col].values.astype(float)
                        if np.std(pred_vals) > 0:
                            bv_corr = np.corrcoef(pred_vals, true_vals)[0, 1]
                            mse = np.mean((pred_vals - true_vals) ** 2)
                    if true_pheno is not None and np.std(pred_vals) > 0:
                        pheno_vals = pred_sub.loc[common_ids, args.pheno_col].values.astype(float)
                        pheno_corr = np.corrcoef(pred_vals, pheno_vals)[0, 1]

        # 关系矩阵类型：优先从vars文件读取，缺失时回退到模型默认映射
        rm_type_final = rm_type if rm_type is not None else model_to_rm.get(model_name)

        results.append({
            'Model': model_name,
            'h2': h2,
            'Var_G': var_g,
            'RM_Type': rm_type_final,
            'Cor_BV': bv_corr,
            'Cor_Pheno': pheno_corr,
            'MSE': mse
        })

        print(f"\n  {model_name}: h2={h2}, RM={rm_type_final}, BV_corr={bv_corr}, MSE={mse}")

    csv_file = os.path.join(args.out, "model_comparison.csv")
    csv_data = []
    for r in results:
        if r['Cor_BV'] is not None or r['Cor_Pheno'] is not None:
            csv_data.append({
                'Model': r['Model'],
                'Cor_BV': r['Cor_BV'],
                'Cor_Pheno': r['Cor_Pheno'],
                'MSE': r['MSE'],
                'relationship_matrix': r['RM_Type']
            })
    if csv_data:
        pd.DataFrame(csv_data).to_csv(csv_file, index=False)
        print(f"\n  保存: {csv_file}")

    # Bilingual summary files
    summary_file_en = os.path.join(args.out, "model_evaluation_summary.txt")
    with open(summary_file_en, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Five-Model Prediction Breeding Value Comparison Summary\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Result directory: {args.result_dir}\n")
        f.write(f"Phenotype file: {args.pheno}\n")
        f.write(f"Validation samples: {len(pred_sub)}\n\n")

        f.write("=" * 80 + "\n")
        f.write("1. Model Description\n")
        f.write("=" * 80 + "\n\n")
        f.write("1. LM (Linear Model): Linear model using HE regression for variance component estimation (GA matrix)\n")
        f.write("2. BLUP: Genomic BLUP using EMAI-REML for variance component estimation (GA matrix)\n")
        f.write("3. PBLUP: BLUP based on pedigree relationship matrix (PA matrix)\n")
        f.write("4. GBLUP: BLUP based on genomic relationship matrix (GA matrix)\n")
        f.write("5. SSBLUP: Single-step genomic BLUP, combining pedigree and genomic information (HA matrix)\n\n")
        f.write("NOTE: The BLUP model in piggpa-G uses GA (genomic relationship matrix), which is identical to\n")
        f.write("      the GBLUP implementation; HIBLUP's BLUP uses PA (pedigree relationship matrix). The two are\n")
        f.write("      NOT directly numerically comparable despite sharing the same model name.\n\n")

        f.write("=" * 80 + "\n")
        f.write("2. Training Results Comparison\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'Model':<16s} {'Heritability(h2)':<14s} {'VarComp(V_G)':<16s} {'RelMatrix Type'}\n")
        f.write("-" * 72 + "\n")
        for r in results:
            h2_str = f"{r['h2']:.4f}" if r['h2'] is not None else "N/A"
            vg_str = f"{r['Var_G']:.4f}" if r['Var_G'] is not None else "N/A"
            # LM使用GA矩阵但通过HE回归估计，需特别标注以示区分
            if r['Model'] == 'LM':
                rm_str = 'GA (HE regression)'
            else:
                rm_str = rm_desc.get(r['RM_Type'], 'N/A')
            f.write(f"{r['Model']:<16s} {h2_str:<14s} {vg_str:<16s} {rm_str}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("3. Prediction Results Comparison\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'Model':<16s} {'BV Correlation':<14s} {'Pheno Correlation':<14s} {'MSE'}\n")
        f.write("-" * 72 + "\n")
        for r in results:
            bv_str = f"{r['Cor_BV']:.4f}" if r['Cor_BV'] is not None else "N/A*"
            ph_str = f"{r['Cor_Pheno']:.4f}" if r['Cor_Pheno'] is not None else "N/A*"
            mse_str = f"{r['MSE']:.4f}" if r['MSE'] is not None else "N/A*"
            f.write(f"{r['Model']:<16s} {bv_str:<14s} {ph_str:<14s} {mse_str}\n")
        f.write("\n* PBLUP is based on pedigree relationship matrix; validation samples have no pedigree relationship with training samples, so direct prediction is not possible\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("4. Conclusion\n")
        f.write("=" * 80 + "\n\n")

        valid_results = [r for r in results if r['Cor_BV'] is not None]
        if valid_results:
            best = max(valid_results, key=lambda x: x['Cor_BV'])
            f.write(f"1. Best model: {best['Model']} (BV correlation: {best['Cor_BV']:.4f})\n\n")

        f.write("2. Model characteristics:\n")
        f.write("   - LM model uses HE regression; heritability estimates may be biased upward\n")
        f.write("   - BLUP and GBLUP produce identical results because both use the GA (genomic) matrix\n")
        f.write("   - SSBLUP combines pedigree and genomic information; prediction accuracy may be slightly higher\n")
        f.write("   - PBLUP cannot predict individuals without pedigree relationship to training samples\n")

    print(f"  保存: {summary_file_en}")

    summary_file_zh = os.path.join(args.out, "model_evaluation_summary-zh.txt")
    with open(summary_file_zh, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("五种模型预测育种值比较分析摘要\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"结果目录: {args.result_dir}\n")
        f.write(f"表型文件: {args.pheno}\n")
        f.write(f"验证样本数: {len(pred_sub)}\n\n")

        f.write("=" * 80 + "\n")
        f.write("一、模型说明\n")
        f.write("=" * 80 + "\n\n")
        f.write("1. LM (Linear Model): 线性模型，使用HE回归估计方差组分（GA矩阵）\n")
        f.write("2. BLUP: 使用EMAI-REML估计方差组分的基因组BLUP（GA矩阵）\n")
        f.write("3. PBLUP: 基于系谱关系矩阵的BLUP（PA矩阵）\n")
        f.write("4. GBLUP: 基于基因组关系矩阵的BLUP（GA矩阵）\n")
        f.write("5. SSBLUP: 单步基因组BLUP，结合系谱和基因组信息（HA矩阵）\n\n")
        f.write("注意: piggpa-G的BLUP模型使用GA（基因组关系矩阵），与GBLUP实现相同；\n")
        f.write("      HIBLUP的BLUP使用PA（系谱关系矩阵），两者不可直接数值对比。\n\n")

        f.write("=" * 80 + "\n")
        f.write("二、训练结果比较\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'模型':<16s} {'遗传力(h²)':<14s} {'方差组分(V_G)':<16s} {'关系矩阵类型'}\n")
        f.write("-" * 72 + "\n")
        for r in results:
            h2_str = f"{r['h2']:.4f}" if r['h2'] is not None else "N/A"
            vg_str = f"{r['Var_G']:.4f}" if r['Var_G'] is not None else "N/A"
            # LM使用GA矩阵但通过HE回归估计，需特别标注以示区分
            if r['Model'] == 'LM':
                rm_str = 'GA (HE回归)'
            else:
                rm_str = rm_desc.get(r['RM_Type'], 'N/A')
            f.write(f"{r['Model']:<16s} {h2_str:<14s} {vg_str:<16s} {rm_str}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("三、预测结果比较\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'模型':<16s} {'育种值相关性':<14s} {'表型相关性':<14s} {'均方误差(MSE)'}\n")
        f.write("-" * 72 + "\n")
        for r in results:
            bv_str = f"{r['Cor_BV']:.4f}" if r['Cor_BV'] is not None else "N/A*"
            ph_str = f"{r['Cor_Pheno']:.4f}" if r['Cor_Pheno'] is not None else "N/A*"
            mse_str = f"{r['MSE']:.4f}" if r['MSE'] is not None else "N/A*"
            f.write(f"{r['Model']:<16s} {bv_str:<14s} {ph_str:<14s} {mse_str}\n")
        f.write("\n* PBLUP基于系谱关系矩阵，验证样本与训练样本无亲缘关系，无法直接预测\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("四、结论\n")
        f.write("=" * 80 + "\n\n")

        valid_results = [r for r in results if r['Cor_BV'] is not None]
        if valid_results:
            best = max(valid_results, key=lambda x: x['Cor_BV'])
            f.write(f"1. 最优模型: {best['Model']} (育种值相关性: {best['Cor_BV']:.4f})\n\n")

        f.write("2. 模型特点:\n")
        f.write("   - LM模型使用HE回归，估计的遗传力可能偏高\n")
        f.write("   - BLUP和GBLUP结果相同，因为都使用GA（基因组）矩阵\n")
        f.write("   - SSBLUP结合了系谱和基因组信息，预测准确性可能略高\n")
        f.write("   - PBLUP无法预测与训练样本无亲缘关系的个体\n")

    print(f"  保存: {summary_file_zh}")

    print("\n[可视化] 生成对比图...")
    models_with_pred = []
    pred_data = {}
    for r in results:
        if r['Cor_BV'] is not None or r['Cor_Pheno'] is not None:
            model_dir = os.path.join(args.result_dir, r['Model'])
            pred_file_map = {
                'LM': 'lm_pred.bv',
                'BLUP': 'blup_pred.bv',
                'GBLUP': 'gblup_pred.bv',
                'SSBLUP': 'ssblup_pred.bv',
            }
            pf = os.path.join(model_dir, pred_file_map.get(r['Model'], ''))
            if os.path.exists(pf):
                pred_df = pd.read_csv(pf, sep='\t')
                pred_data[r['Model']] = pred_df['add_a1'].values
                models_with_pred.append(r['Model'])

    if len(models_with_pred) > 0:
        colors = PALETTE_A[:len(models_with_pred)]

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        ax1 = axes[0, 0]
        data_to_plot = [pred_data[m] for m in models_with_pred]
        bp = ax1.boxplot(data_to_plot, labels=models_with_pred, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        ax1.set_ylabel('Predicted Breeding Value')
        ax1.set_title('Breeding Value Distribution by Model')

        ax2 = axes[0, 1]
        h2_vals = [r['h2'] for r in results if r['h2'] is not None]
        model_names_h2 = [r['Model'] for r in results if r['h2'] is not None]
        if h2_vals:
            bars = ax2.bar(model_names_h2, h2_vals, color=PALETTE_A[:len(h2_vals)], edgecolor='#333333')
            ax2.set_ylabel('Heritability (h²)')
            ax2.set_title('Estimated Heritability by Model')
            for bar, val in zip(bars, h2_vals):
                ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10)

        ax3 = axes[1, 0]
        corr_vals = [r['Cor_BV'] for r in results if r['Cor_BV'] is not None]
        corr_models = [r['Model'] for r in results if r['Cor_BV'] is not None]
        if corr_vals:
            bars = ax3.bar(corr_models, corr_vals, color=PALETTE_A[:len(corr_vals)], edgecolor='#333333')
            ax3.set_ylabel('Correlation with True BV')
            ax3.set_title('Prediction Accuracy by Model')
            for bar, val in zip(bars, corr_vals):
                ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10)

        ax4 = axes[1, 1]
        for i, model in enumerate(models_with_pred):
            pred = pred_data[model]
            ax4.hist(pred, bins=30, alpha=0.5, label=model, color=colors[i])
        ax4.set_xlabel('Predicted Breeding Value')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Predicted BV Histograms')
        leg = ax4.legend()
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')

        plt.tight_layout()
        for fmt in ['pdf', 'png']:
            fig.savefig(os.path.join(args.out, f"model_comparison.{fmt}"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  对比图已保存: model_comparison.pdf/png")

    print("\n" + "=" * 70)
    print("模型对比分析完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
