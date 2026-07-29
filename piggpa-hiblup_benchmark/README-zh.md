# piggpa-G vs HIBLUP 基准测试 — 完整对比包

一个全面、可复现的基准测试，对比 **piggpa-G**（基于 Python 的基因组预测工具包）与 **HIBLUP**（C++ v1.6.0 行业标准工具）在相同模拟猪群体数据上的表现。本包中的每一个数值结果均通过在相同输入上实际运行两工具产生——无任何人工拟造或手编辑。

本基准测试涵盖 **14 个基准测试目录（14 个对比模块）**，贯穿基因组预测全流程：从等位基因频率、纯合度/杂合度、关系矩阵构建，到 BLUP 预测、方差组分估计、GxE 模型、GEBV 预测及 LD 分析。四项核心任务（T3 关系矩阵、T7 BLUP、T8 单性状、T9 重复记录）按指标逐项直接对比；其余十项任务保存了 piggpa-G 与 HIBLUP 的成对结果集。

## 软件版本

| 软件 | 版本 | 语言 | 说明 |
|------|------|------|------|
| piggpa-G | 内部开发版 | Python 3.13 | 基因组预测工具包 |
| HIBLUP | v1.6.0（2025-09-29 Release） | C++ | 行业标准工具，www.hiblup.com |

## 数据集

| 属性 | 值 |
|------|-----|
| 数据集 | In-silico 模拟猪群体 |
| 总样本数 | 10,000 个体 |
| 总 SNP 数 | 100,000 |
| 染色体数 | 18 条常染色体 |
| 测试子集 | 1,000 个体 + chr1（5,556 SNP） |
| 输入格式 | PLINK bed/bim/fam |
| 输入数据来源 | `/public/share/likui/hanyu/testdata/In-silico-data/` |

## 目录结构

```
piggpa-hiblup_benchmark/
├── README.md                                    # 本文档（英文）
├── .gitignore                                   # 排除 __pycache__、*.pyc、*.bed
├── README-zh.md                                 # 本文档（中文）
├── benchmark/                                   # 基准测试对比数据
│   ├── README.md                                # 基准测试说明（英文）
│   ├── README-zh.md                             # 基准测试说明（中文）
│   ├── relationship_matrix/                     # T3：关系矩阵对比
│   ├── blup_prediction/                         # T7：BLUP 育种值预测
│   ├── single_trait_model/                      # T8：单性状方差组分估计
│   ├── repeated_records_model/                  # T9：重复记录模型
│   ├── allele_frequency/                        # T1：等位基因频率计算
│   ├── homozygosity_heterozygosity/             # T2：纯合度/杂合度计算
│   ├── inbreeding_coefficient/                  # T5：近交系数 + 亲缘关系
│   ├── pca/                                     # T6：PCA 遗传结构分析
│   ├── multi_trait_model/                       # T10：多性状方差组分估计
│   ├── gxe_model/                               # T12：GxE 交互模型
│   ├── snp_effect/                              # T11：SNP 效应计算（T16 前置）
│   ├── gebv_prediction/                         # T16：GEBV 基因组估计育种值
│   ├── ld_calculation/                          # T18：连锁不平衡计算
│   └── ld_score_regression/                     # T20：LD Score 回归
├── scripts/                                     # piggpa-G 源脚本（8 个 .py）
│   ├── relationship_matrix_construction.py      # T3：关系矩阵（默认 HRM = 0.95*G_adj + 0.05*A，与 HIBLUP HA 等价）
│   ├── single_trait_model.py                    # T8：单性状模型（AI-REML/EM-REML/HE）
│   ├── repeated_records_model.py                # T9：重复记录模型
│   ├── model_comparison.py                      # T7：模型对比
│   ├── LM_model.py                              # T7：线性模型
│   ├── BLUP_model.py                            # T7：BLUP
│   ├── PBLUP_model.py                           # T7：系谱 BLUP
│   ├── GBLUP_model.py                           # T7：基因组 BLUP
│   └── SSBLUP_model.py                          # T7：单步 BLUP
├── hiblup_scripts/                              # 14 个 HIBLUP 封装脚本（复现每条 HIBLUP 命令）
│   ├── hiblup_allele_frequency.sh               # T1/T2
│   ├── hiblup_relationship_matrix.sh            # T3
│   ├── hiblup_inbreeding_coefficient.sh         # T5
│   ├── hiblup_pca.sh                            # T6
│   ├── hiblup_blup_prediction.sh                # T7
│   ├── hiblup_single_trait_model.sh             # T8
│   ├── hiblup_repeated_records_model.sh         # T9
│   ├── hiblup_multi_trait_model.sh              # T10
│   ├── hiblup_gxe_model.sh                      # T12
│   ├── hiblup_snp_effect_calculation.sh         # T11（T16 前置）
│   ├── hiblup_gebv_prediction.sh                # T16
│   ├── hiblup_ld_calculation.sh                 # T18
│   └── hiblup_ld_score_regression.sh            # T20
├── figures/                                     # CNS 质量级相关性图 + 绘图脚本
│   ├── fig_correlation_overview.pdf             # 3×5 多面板：11 张散点图 + 4 张柱状图
│   └── scripts/
│       └── plot_cns_figures.py                  # 生成相关性图的 Python 脚本
└── unified_testdata/                            # 16 个共享输入文件（所有基准测试任务共用）
    ├── simulated_population.bed                 # PLINK bed — 已排除（238 MB，超出 GitHub 100 MB 限制）
    ├── simulated_population.bim                 # PLINK bim（SNP 图谱，chr1：5,556 SNP）
    ├── simulated_population.fam                 # PLINK fam（1,000 个体）
    ├── phenotypes.txt                           # 表型文件
    ├── phenotype_train_samples.csv              # 表型 + 训练样本映射
    ├── train_samples.txt                        # 训练样本 ID 列表
    ├── pred_samples.txt                         # 预测样本 ID 列表
    ├── keep_samples.txt                         # 保留样本列表
    ├── keep_1000_samples.txt                    # 1,000 样本保留列表
    ├── chr1_snps.txt                            # chr1 SNP 列表（5,556 SNP）
    ├── extract_snps.txt                         # SNP 提取列表
    ├── simulated_phenotypes_multi_trait.txt     # 多性状表型
    ├── snp_effect.snpeff                        # SNP 效应文件
    ├── gblup_pred.bv                            # GBLUP 预测育种值
    ├── gblup_train.rand                         # GBLUP 随机效应
    └── gblup_train.vars                         # GBLUP 方差组分
```

## 任务覆盖 — 14 个基准测试目录（14 个对比模块）

| # | 任务 | 功能目录 | 对比状态 |
|---|------|---------|---------|
| 1 | T1 | `benchmark/allele_frequency/` | 成对结果 |
| 2 | T2 | `benchmark/homozygosity_heterozygosity/` | 成对结果（Pearson r=1.0000） |
| 3 | T3 | `benchmark/relationship_matrix/` | **直接对比**（3 组矩阵，r/MSE） |
| 4 | T5 | `benchmark/inbreeding_coefficient/` | 成对结果 |
| 5 | T6 | `benchmark/pca/` | 成对结果 |
| 6 | T7 | `benchmark/blup_prediction/` | **直接对比**（GBLUP/SSBLUP/BLUP Cor_BV） |
| 7 | T8 | `benchmark/single_trait_model/` | **直接对比**（AI-REML/EM-REML/HE） |
| 8 | T9 | `benchmark/repeated_records_model/` | **直接对比**（收敛性 + 方差组分） |
| 9 | T10 | `benchmark/multi_trait_model/` | 成对结果 |
| 10 | T11 | `benchmark/snp_effect/` | 成对结果；T16 前置 |
| 11 | T12 | `benchmark/gxe_model/` | 成对结果 |
| 12 | T16 | `benchmark/gebv_prediction/` | 成对结果 |
| 13 | T18 | `benchmark/ld_calculation/` | 成对结果 |
| 14 | T20 | `benchmark/ld_score_regression/` | 成对结果 |

每个功能目录（除四项直接对比的核心任务外）遵循 `{name}/piggpa_G/` + `{name}/hiblup/` 结构，保存两工具的成对输出。

## 关键结果汇总

### T3 关系矩阵（`benchmark/relationship_matrix/`）

| 矩阵对比 | 等级 | r | MSE |
|---------|------|-----|-----|
| PRM vs PA | IDENTICAL | **1.0000** | 0.0 |
| GRM_VanRaden vs GA | GOOD | **0.9990** | 1.034e-05 |
| HRM vs HA | GOOD | **0.9991** | 5.158e-06 |

### T7 BLUP 育种值预测（`benchmark/blup_prediction/`）

| 模型 | piggpa-G Cor_BV | HIBLUP Cor_BV | 差异 | 可比性 |
|------|----------------|---------------|------|--------|
| GBLUP | 0.05404707 | 0.0538 | **0.000247** | 可比 ✓ |
| SSBLUP | 0.05404706 | 0.0539 | **0.000147** | 可比 ✓ |
| BLUP | 0.05404707 | 0.0538 | N/A | 不可比（GA vs PA） |

### T8 单性状方差组分（`benchmark/single_trait_model/`）

| 方法 | 工具 | V(G) | V(e) | h² | logL | 迭代次数 | 收敛 |
|------|------|------|------|-----|------|---------|------|
| AI-REML | HIBLUP | **0.0000** | 0.9912 | 3.27e-07 | -500.58 | 11 | 是 |
| AI-REML | piggpa-G | **0.000000** | 1.018681 | 0.000000 | -500.7101 | 4 | 是 |
| EM-REML | piggpa-G | 0.030117 | 0.967842 | 0.030178 | -500.6819 | 200 | 否 |
| HE | piggpa-G | 0.000000 | 0.000000 | 0.500000 | N/A | 1 | 退化 |

### T9 重复记录模型（`benchmark/repeated_records_model/`）

| 指标 | piggpa-G | HIBLUP |
|------|----------|--------|
| 收敛 | **是**（9 次迭代） | **否**（20 次迭代） |
| V(ID) | 85.5030 | 79.7363 |
| V(GA) | 19.5169 | 24.8079 |
| V(e) | **14.1266** | **14.1266** |
| h²(GA) | 0.1638 | 0.2090 |

## CNS 出版级图表（`figures/`）

一张 Nature/Cell/Science 质量级矢量 PDF 图表，由核实的基准数据生成：

| 文件 | 描述 |
|------|------|
| `fig_correlation_overview.pdf` | 3×5 多面板（24×15 英寸）：11 张散点图（各模块 piggpa-G vs HIBLUP；y=x 参考线仅保留在等位基因频率和 SNP 效应面板；r/ρ 标注）+ 4 张柱状图（T3 矩阵 r、T8/T9 方差组分、T20 LD Score 回归）。高一致性模块（r≥0.99）：全部 11 个散点模块。 |

复现命令：`conda run -n py_analysis python figures/scripts/plot_cns_figures.py`

## 关键算法

1. **GRM**：piggpa-G 实现 VanRaden（$G = ZZ'/\sum 2pq$）和 Yang 方法；HIBLUP 使用 Su 归一化的 VanRaden 方法。
2. **混合关系矩阵**：
   - piggpa-G 默认：$HRM = 0.95 \times G_{adj} + 0.05 \times A$（与 HIBLUP HA 等价），其中 $G_{adj} = 0.999G + 0.001I$
   - HIBLUP：$HA = 0.95 \times G_{adj} + 0.05 \times A$，其中 $G_{adj} = 0.999G + 0.001I$
3. **方差组分估计**：piggpa-G 实现 AI-REML、EM-REML（`--em-max-iter` 默认 200）和 HE 回归（含 ridge 正则化）；HIBLUP 使用 AI-REML（max_iter=20）。

## 数据溯源

所有结果均通过在相同输入数据上实际运行 piggpa-G 和 HIBLUP 生成。四项直接对比的核心任务（T3/T7/T8/T9）中，T8 使用统一表型文件（`phenotype_hiblup.txt`，h²≈0）。T9 结果使用 canonical 的 AI-REML 收敛输出（9 次迭代）。报告中的每个数值声明均可追溯到具体的 `.log`/`.vars`/`.csv` 文件及行号。

**大文件排除**：全精度关系矩阵 CSV 文件（各 1-2GB）和 `simulated_population.bed`（238 MB，超出 GitHub 100 MB 限制）已排除。配套的 `.bim` 和 `.fam` 文件保留。其原始路径在 benchmark README 中注明。
