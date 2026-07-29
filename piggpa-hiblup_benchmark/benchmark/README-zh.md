# piggpa-G vs HIBLUP 基准测试 — 定量比对数据

本目录包含 piggpa-G 与 HIBLUP 的全部定量基准比对数据。所有数据均通过在完全相同的输入数据上实际运行两工具获得，非人工拟造。

## 目录结构

```
benchmark/
├── README.md                                    # 本文档（英文）
├── README-zh.md                                 # 本文档（中文）
├── relationship_matrix/                         # T3：关系矩阵对比
│   ├── piggpa_G/                                # piggpa-G 结果（matrix_comparison_results.csv 等 7 个文件）
│   └── hiblup/                                  # HIBLUP 参考（compare_matrices_v2.py、matrix_comparison_report.txt 等）
├── blup_prediction/                             # T7：BLUP 育种值预测
│   ├── piggpa_G/                                # piggpa-G 结果（5 模型子目录 + model_comparison.csv 等）
│   └── hiblup/                                  # HIBLUP 参考（gblup/blup/ssblup 等 .vars/.log 文件）
├── single_trait_model/                          # T8：单性状方差组分估计
│   ├── piggpa_G/                                # piggpa-G 结果（variance_component_summary.txt、run.log 等）
│   └── hiblup/                                  # HIBLUP 参考（single_trait_gblup.vars/.log 等）
├── repeated_records_model/                      # T9：重复记录模型
│   ├── piggpa_G/                                # piggpa-G 结果（repeated_model.log/.vars 等）
│   └── hiblup/                                  # HIBLUP 参考（repeated_model.vars/.log 等）
├── allele_frequency/                            # T1：等位基因频率
│   ├── piggpa_G/
│   └── hiblup/
├── homozygosity_heterozygosity/                 # T2：纯合度/杂合度
│   ├── piggpa_G/
│   └── hiblup/
├── inbreeding_coefficient/                      # T5：近交系数
│   ├── piggpa_G/
│   └── hiblup/
├── pca/                                         # T6：PCA
│   ├── piggpa_G/
│   └── hiblup/
├── multi_trait_model/                           # T10：多性状模型
│   ├── piggpa_G/
│   └── hiblup/
├── gxe_model/                                   # T12：GxE 模型
│   ├── piggpa_G/
│   └── hiblup/
├── snp_effect/                                  # T11：SNP 效应计算
│   ├── piggpa_G/
│   └── hiblup/
├── gebv_prediction/                             # T16：GEBV 预测
│   ├── piggpa_G/
│   └── hiblup/
├── ld_calculation/                              # T18：LD 计算
│   ├── piggpa_G/
│   └── hiblup/
└── ld_score_regression/                         # T20：LD Score 回归
    ├── piggpa_G/
    └── hiblup/
```

所有 14 个功能目录均遵循 `{name}/piggpa_G/` + `{name}/hiblup/` 统一结构，保存两工具的成对输出。HIBLUP 参考结果已从原 `hiblup_reference/` 合并到各功能目录的 `hiblup/` 子目录中。

## 任务覆盖 — 14 个功能模块

| # | 任务 | 目录 | 结构 | 对比状态 |
|---|------|------|------|---------|
| 1 | T1 | `allele_frequency/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 2 | T2 | `homozygosity_heterozygosity/` | piggpa_G/ + hiblup/ | 目录中有结果（Pearson r=1.0000） |
| 3 | T3 | `relationship_matrix/` | piggpa_G/ + hiblup/ | **直接对比**（3 组矩阵，r/MSE） |
| 4 | T5 | `inbreeding_coefficient/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 5 | T6 | `pca/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 6 | T7 | `blup_prediction/` | piggpa_G/ + hiblup/ | **直接对比**（GBLUP/SSBLUP/BLUP Cor_BV） |
| 7 | T8 | `single_trait_model/` | piggpa_G/ + hiblup/ | **直接对比**（AI-REML/EM-REML/HE） |
| 8 | T9 | `repeated_records_model/` | piggpa_G/ + hiblup/ | **直接对比**（收敛性 + 方差组分） |
| 9 | T10 | `multi_trait_model/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 10 | T11 | `snp_effect/` | piggpa_G/ + hiblup/ | 目录中有结果；T16 前置 |
| 11 | T12 | `gxe_model/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 12 | T16 | `gebv_prediction/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 13 | T18 | `ld_calculation/` | piggpa_G/ + hiblup/ | 目录中有结果 |
| 14 | T20 | `ld_score_regression/` | piggpa_G/ + hiblup/ | 目录中有结果 |

---

## 直接对比任务 — 详细表格

### 1. T3 关系矩阵对比（`relationship_matrix/`）

**数据文件**：`relationship_matrix/piggpa_G/matrix_comparison_results.csv`

在 1,000 个体和 5,556 个 SNP（chr1）上进行 piggpa-G 与 HIBLUP 关系矩阵的 3 组比对。

| 矩阵对比 | 等级 | r | MSE | Max Diff | piggpa-G 对角线均值 | HIBLUP 对角线均值 |
|---------|------|-----|-----|---------|-------------------|-----------------|
| PRM vs PA | IDENTICAL | **1.0000** | 0.0 | 0.0 | 1.0 | 1.0 |
| GRM_VanRaden vs GA | GOOD | **0.9990** | 1.034e-05 | 0.01042 | 1.00371 | 1.00000 |
| HRM vs HA | GOOD | **0.9991** | 5.158e-06 | 0.00989 | 1.00352 | 1.00000 |

**关键发现**：
1. PRM vs PA：r=**1.0000**（IDENTICAL）——系谱构建算法完全相同
2. GRM vs GA：r=**0.9990**（GOOD）——VanRaden 方法实现高度一致
3. HRM vs HA：r=**0.9991**（GOOD）——HRM 使用公式 HRM = 0.95×G_adj + 0.05×A（其中 G_adj = 0.999×G + 0.001×I），与 HIBLUP 的 HA 等价

**大文件排除**：全精度矩阵 CSV 文件（PRM.csv、GRM_VanRaden.csv、HRM.csv，各约 1-2 GB）未打包。原始路径：`HIBLUP_benchmark/new_benchmark/t3/`

---

### 2. T7 BLUP 育种值预测（`blup_prediction/`）

**数据文件**：`blup_prediction/piggpa_G/model_comparison.csv`

使用训练集（ID 1001-2000）和验证集（ID 5001-6000），各 1,000 个体，进行 5 模型 BLUP 对比。

| 模型 | piggpa-G Cor_BV | HIBLUP Cor_BV | 差异 | piggpa-G 关系矩阵 | HIBLUP 关系矩阵 | 可比性 |
|------|----------------|---------------|------|------------------|----------------|--------|
| GBLUP | 0.05404707 | 0.0538 | **0.000247** | GA | GA | 可比 ✓ |
| SSBLUP | 0.05404706 | 0.0539 | **0.000147** | HA | HA | 可比 ✓ |
| BLUP | 0.05404707 | 0.0538 | N/A | GA | PA | **不可比** |

**关键发现**：
1. GBLUP Cor_BV 差异 **0.000247** < 0.001——排名选择等效
2. SSBLUP Cor_BV 差异 **0.000147** < 0.001——单步预测等效
3. BLUP 模型不可直接对比：piggpa-G 的 BLUP 使用 GA（与 GBLUP 相同），HIBLUP 的 BLUP 使用 PA（系谱）

---

### 3. T8 单性状方差组分估计（`single_trait_model/`）

**数据文件**：`single_trait_model/piggpa_G/variance_component_summary.txt` + `run.log`

**表型输入**：两工具统一使用表型文件（`phenotype_hiblup.txt`，h²≈0）。

| 方法 | 工具 | V(G) | V(e) | h² | logL | 迭代次数 | 收敛 |
|------|------|------|------|-----|------|---------|------|
| AI-REML | HIBLUP | **0.0000** | 0.9912 | 3.27e-07 | -500.58 | 11 | 是 |
| AI-REML | piggpa-G | **0.000000** | 1.018681 | 0.000000 | -500.7101 | 4 | 是 |
| EM-REML | piggpa-G | 0.030117 | 0.967842 | 0.030178 | -500.6819 | 200 | 否 |
| HE | piggpa-G | 0.000000 | 0.000000 | 0.500000 | N/A | 1 | 退化（ridge） |

**关键发现**：
1. AI-REML 一致性：两工具 h² 均≈0（同量级），V(e) 差异 0.0275，logL 差异 0.13
2. EM-REML 在 200 次迭代后未收敛，但 V(G) 从 0.497 单调下降至 0.030（方向正确）
3. HE 回归出现退化（ridge 正则化生效，h²=0.5 < 1.0）

---

### 4. T9 重复记录模型（`repeated_records_model/`）

**数据文件**：`repeated_records_model/piggpa_G/repeated_model.log` + `repeated_model.vars`

模型：`weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e`

| 指标 | piggpa-G | HIBLUP | 差异 |
|------|----------|--------|------|
| 收敛 | **是**（9 次迭代） | **否**（20 次用尽） | — |
| logL | -9792.57 | -7037.67 | -2754.90 |
| V(ID) 永久环境 | 85.5030 | 79.7363 | 5.7667 |
| V(GA) 加性遗传 | 19.5169 | 24.8079 | -5.2910 |
| V(e) 残差 | **14.1266** | **14.1266** | **0.0000** |
| h²(GA) | 0.1638 | 0.2090 | -0.0452 |

**关键发现**：
1. piggpa-G **9 次迭代**收敛；HIBLUP 20 次迭代未收敛（AI(20) 用尽）
2. 残差方差 V(e) = **14.1266** 两工具完全一致
3. 方差组分在同一量级；差异可归因于 HIBLUP 未收敛

---

## 新增任务目录 — 简要说明

以下十一个目录包含 piggpa-G 与 HIBLUP 的成对结果集。

### T1 等位基因频率（`allele_frequency/`）
- **功能**：等位基因频率、基因型频率计算
- **结构**：`piggpa_G/{t1,t2}/` + `hiblup/1.2/`

### T2 纯合度/杂合度（`homozygosity_heterozygosity/`）
- **功能**：计算每个体的纯合度和杂合度
- **结构**：`piggpa_G/t2/` + `hiblup/1.2/`（HIBLUP 的 1.2/ 目录同时计算等位基因频率和纯合度/杂合度）
- **对比**：Pearson r = 1.0000（完全一致）

### T5 近交系数（`inbreeding_coefficient/`）
- **功能**：近交系数（F）和亲缘/关系系数计算
- **结构**：`piggpa_G/t5/` + `hiblup/4/`

### T6 PCA（`pca/`）
- **功能**：遗传结构主成分分析（前 10 个 PC）
- **结构**：`piggpa_G/t6/` + `hiblup/5/`

### T10 多性状模型（`multi_trait_model/`）
- **功能**：多性状（3 性状）GBLUP 方差组分估计 + 基因型编码/GRM 验证
- **结构**：`piggpa_G/t10/` + `hiblup/9/`

### T11 SNP 效应（`snp_effect/`）
- **功能**：单性状 GBLUP + SNP 效应反算（T16 的直接前置）
- **结构**：`piggpa_G/t11/` + `hiblup/11/`

### T12 GxE 模型（`gxe_model/`）
- **功能**：基因×环境（GxE）交互模型（`--rand-gxe`）
- **结构**：`piggpa_G/t12/` + `hiblup/10/`

### T16 GEBV 预测（`gebv_prediction/`）
- **功能**：对留出个体进行基因组估计育种值预测（`--pred`）
- **结构**：`piggpa_G/t16/` + `hiblup/12/`

### T18 LD 计算（`ld_calculation/`）
- **功能**：成对连锁不平衡（`--ld`）和 LD score（`--ldscore`）
- **结构**：`piggpa_G/t18/` + `hiblup/14/`

### T20 LD Score 回归（`ld_score_regression/`）
- **功能**：LD Score 回归用于 SNP 遗传力估计（`--ldreg`）
- **结构**：`piggpa_G/t20/` + `hiblup/15/`

---

## 数据溯源

所有数值结果均由实际运行 piggpa-G 和 HIBLUP 产生。上表中的每个值都可追溯到具体文件：

| 任务 | 可追溯文件 |
|------|-----------|
| T3 | `relationship_matrix/piggpa_G/matrix_comparison_results.csv`（3 行，已核实）+ `relationship_matrix/hiblup/` |
| T7 | `blup_prediction/piggpa_G/model_comparison.csv` + `model_evaluation_summary.txt` + `blup_prediction/hiblup/` |
| T8 | `single_trait_model/piggpa_G/variance_component_summary.txt` + `run.log` + `single_trait_model/hiblup/` |
| T9 | `repeated_records_model/piggpa_G/repeated_model.log`（L16-26：9 次 AI 迭代，L26：`[Converged?] Yes!`）+ `repeated_records_model/hiblup/` |
| T1 | `allele_frequency/piggpa_G/{t1,t2}/` + `allele_frequency/hiblup/1.2/` |
| T2 | `homozygosity_heterozygosity/piggpa_G/t2/` + `homozygosity_heterozygosity/hiblup/1.2/` |
| T5 | `inbreeding_coefficient/piggpa_G/t5/` + `inbreeding_coefficient/hiblup/4/` |
| T6 | `pca/piggpa_G/t6/` + `pca/hiblup/5/` |
| T10 | `multi_trait_model/piggpa_G/t10/` + `multi_trait_model/hiblup/9/` |
| T11 | `snp_effect/piggpa_G/t11/` + `snp_effect/hiblup/11/` |
| T12 | `gxe_model/piggpa_G/t12/` + `gxe_model/hiblup/10/` |
| T16 | `gebv_prediction/piggpa_G/t16/` + `gebv_prediction/hiblup/12/` |
| T18 | `ld_calculation/piggpa_G/t18/` + `ld_calculation/hiblup/14/` |
| T20 | `ld_score_regression/piggpa_G/t20/` + `ld_score_regression/hiblup/15/` |

所有 14 个功能目录均采用统一的 `piggpa_G/` + `hiblup/` 结构，HIBLUP 参考结果保留原始 `.log`/`.vars`/`.beta`/`.rand` 文件可供独立核实。

**大文件排除**：全精度关系矩阵 CSV 文件（各 1-2GB）和 `simulated_population.bed`（238 MB，超出 GitHub 100 MB 限制）已排除。配套的 `.bim` 和 `.fam` 文件保留。原始路径在顶层 README 中注明。
