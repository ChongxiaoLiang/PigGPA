# PigGPA 工作流执行基准测试

本基准测试对位于 `/public/share/likui/liangcx/bole/skills/piggpa-G/` 的全部 26 个 PigGPA-G 技能在小型代表性测试数据上进行了完整的端到端测试。每个技能被封装在一个独立的任务脚本中，用于测量运行耗时、捕获 stdout/stderr 到日志文件，并写入任务级 JSON 记录。最终的编排脚本按数字顺序运行全部 26 个封装脚本，并通过 `lib/aggregate_results.py` 将各任务 JSON 汇总为单一结果文件。在 26 个任务中，**全部 26 个通过**，工作流执行率为 **100.0%**。全部 26 个技能均能端到端成功运行。

## 目的

本基准测试的目的在于验证每一个 PigGPA-G 技能都能够在与生产环境一致的输入数据形态上端到端运行，并捕获耗时与退出码以保证可复现性。本基准测试属于冒烟测试：仅确认每个技能能接受文档化的参数、定位输入、产生预期的输出文件形态并正常退出。它既不是生产规模下的性能基准测试，也不是与已发表结果对照的正确性基准测试。

## 测试环境

| 项目 | 取值 |
|------|-------|
| 主机名 | ln03 |
| 基准测试时间戳 | 2026-07-17T03:41:23+08:00 |
| PigGPA 技能根目录 | `/public/share/likui/liangcx/bole/skills/piggpa-G` |
| 测试数据根目录 | `/public/share/likui/liangcx/bole/testdata`（清单见 `testdata_list.md`） |
| 结果根目录 | `/public/share/likui/liangcx/bole/skills/workflow_execution/results/` |
| 汇总 JSON | `/public/share/likui/liangcx/bole/skills/workflow_execution/workflow_execution_results.json` |
| 被测系统 | PigGPA (bole/skills/piggpa-G) |
| 基准测试归属 | PigGPA |

工具路径：

| 工具 | 路径 |
|------|------|
| PLINK / GCTA / ADMIXTURE / bcftools | `/public/share/likui/liangcx/software/miniconda3/envs/sys_tools/bin` |
| Rscript | `/public/share/likui/liangcx/software/miniconda3/envs/R/bin/Rscript` |
| Python | `/public/share/likui/liangcx/software/miniconda3/envs/py_analysis/bin/python` (numpy 2.4.6, pandas 2.3.3, scipy 1.16.3) |

## 测试数据输入

测试数据被有意控制在较小规模，以便快速完成端到端冒烟测试。生产规模数据可能呈现不同的耗时特征，并可能触及小型测试数据无法覆盖的代码路径。

| 测试数据文件 | 描述 | 消费者 |
|----------------|-------------|-------------|
| `selectedfortest.vcf` | 用于 PLINK 导入的 VCF | wf001 |
| `testps.bed/bim/fam` | 45 个样本，约 50K SNPs | wf002, wf006, wf008, wf009 |
| `reference.bed/bim/fam` | 2237 个样本，约 258K SNPs | wf003, wf004, wf007, wf010, wf011, wf012, wf013, wf014, wf015, wf017 |
| `reference_50k.bed/bim/fam` | reference 的 50K-SNP 子集 | wf005 |
| `reference_2k.bed/bim/fam` | reference 的 2K-SNP 子集 | wf005 |
| `simulated_population_chr1_part.bed/bim/fam` | 1000 个样本，5556 SNPs (chr1) | wf016, wf018, wf019, wf020, wf021, wf022 |
| `gdp_98samples.bed/bim/fam` | 98 个样本，约 1.75M SNPs (WGS) | wf023, wf024, wf025 |
| `gdp_chip_20.bed` | 20 样本目标芯片 | wf026 |
| `gdp_ref_78.bed` | 78 样本参考芯片 | wf026 |
| `GWAStestphenotype.fam` | 2237 行 (FID IID pheno) | wf003, wf004, wf007, wf017 |
| `snp_effect.snpeff` | TSV (id, add_a1) | wf022 |
| `pca_groups.txt` | 群体分组注释 | wf009 |
| `pheno_snp_effect.csv` | SNP 效应任务所用表型 | wf018, wf019, wf020, wf021 |
| `chr1_snps.txt` | chr1 SNP 列表 | wf011 |
| `train_samples.txt`, `pred_samples.txt`, `all_samples.txt` | 样本子集列表 | wf005, wf022 |
| `simulated_phenotypes_env_1001_2000.txt` | GxE 表型 | wf019 |
| `simulated_phenotypes_multi_trait_1001_2000.txt` | 多性状表型 | wf020 |
| `simulated_phenotypes_repeated_long.txt` | 纵向表型 | wf021 |
| `pedigree.txt` | SSGBLUP 所用系谱 | wf016 |
| `gdp_pop.txt` | 广东猪群体标签 | wf023, wf025 |

## 任务封装脚本

每个技能由一个名为 `run_wfNNN_<skill>.sh` 的任务封装脚本驱动。封装模板的特征如下：

- 在 shebang 之后立即执行 `ulimit -v unlimited`（按工作区规范解除虚拟内存上限，防止大内存任务被 OOM Kill）。
- 从环境变量解析 `PIGGPA_SKILLS_DIR`、`TESTDATA_DIR`、`RESULTS_DIR`、`APP_BIN`、`RSCRIPT_BIN`、`PYTHON_BIN`，并提供指向上述"测试环境"中所列路径的合理默认值。
- 创建任务级输出目录 `results/wfNNN-<skill>/`（例如 `results/wf001-genotype-import/`）。
- 记录开始时间戳，以文档化参数调用底层技能脚本，记录结束时间戳。
- 计算 `WF_ELAPSED`（秒，3 位小数）与 `WF_ELAPSED_HUMAN`（根据时长自动格式化为秒、"X min Y s" 或 "X h Y min Z s"）。
- 调用 `lib/write_task_json.py` 生成 `results/wfNNN-<skill>/wfNNN-<skill>.json`，包含全部指标（elapsed_seconds、status、rc、input_files、output_files、error_message）。
- 返回底层技能的退出码，使编排脚本与汇总报告反映真实的技能结果。

根据技能的实现语言，使用三种调用模式：

1. **PLINK/GCTA/ADMIXTURE shell 技能**：调用技能的 `.sh` 脚本，传入 `--genotype_prefix`、`--phenotype_file`、`--output_prefix`、`--threads`、`--app-bin` 等参数。
2. **Python 技能**：`"$PYTHON_BIN" "$WF_ENTRY" --args...`
3. **R 技能**：`"$RSCRIPT_BIN" "$WF_ENTRY" --args...`

全部 26 个封装脚本使用上述相同的调用模式。

## 运行方式

按数字顺序运行全部 26 个任务：

```bash
bash run_workflow_execution.sh
```

通过传入匹配目标 wf ID 的正则表达式运行子集（编排脚本按数字顺序遍历 `run_wfNNN_*.sh`，并按所给模式过滤）：

```bash
bash run_workflow_execution.sh wf00[1-9]
```

跳过任务执行，仅基于已有的任务级 JSON 重新生成汇总 JSON：

```bash
python lib/aggregate_results.py
```

## 结果汇总

### 汇总统计

| 指标 | 取值 |
|--------|-------|
| 任务总数 | 26 |
| 通过 | 26 |
| 失败 | 0 |
| 总耗时（秒） | 2714.82 |
| 总耗时（可读） | 45 min 15 s |
| 工作流执行率 | 100.0% |
| 失败任务 ID | （无） |

### 任务级结果

| wf_id | skill | task | 耗时 (s) | 状态 | rc |
|-------|-------|------|-------------|--------|----|
| wf001 | genotype-import | Import selectedfortest.vcf to PLINK binary format | 0.111 | passed | 0 |
| wf002 | genotype-qc | PLINK QC on testps (45 samples, MAF/missing/HWE filters) | 0.105 | passed | 0 |
| wf003 | plink-gwas | PLINK linear regression GWAS on reference + GWAStestphenotype | 12.104 | passed | 0 |
| wf004 | gcta-gwas | GCTA MLMA GWAS on reference + GWAStestphenotype (GRM-based) | 1140.557 | passed | 0 |
| wf005 | genomic-selection | 6 GS models (BayesA/B/C/BRR/BL/GBLUP) on reference_2k with 3-fold CV | 46.378 | passed | 0 |
| wf006 | admixture-analysis | ADMIXTURE CV on testps for K=2..3 (small K range for speed) | 22.512 | passed | 0 |
| wf007 | heritability-analysis | GCTA GRM + REML heritability on reference + GWAStestphenotype (Weight trait) | 83.891 | passed | 0 |
| wf008 | ld-pruning | PLINK LD pruning on testps (indep-pairwise 50 5 0.2) | 0.192 | passed | 0 |
| wf009 | pca | PCA on testps (45 samples) with population groups annotation | 9.200 | passed | 0 |
| wf010 | ld | LD decay analysis on reference (all chromosomes) via PopLDdecay (max-dist 300kb) | 157.571 | passed | 0 |
| wf011 | ld-score | LD Score calculation on reference chr1 (max-dist 100kb for feasible runtime) | 5.505 | passed | 0 |
| wf012 | allele-genotype-frequency | Allele & genotype frequency on reference chr1 (batch processing) | 11.252 | passed | 0 |
| wf013 | homozygosity-heterozygosity | SNP & sample homo/heterozygosity on reference chr1 | 14.717 | passed | 0 |
| wf014 | inbreeding-relationship | Inbreeding (GRM/homozygosity/excess) & kinship coefficients on reference chr1 | 22.655 | passed | 0 |
| wf015 | relationship-matrix | PRM/GRM (VanRaden+Yang)/HRM construction on reference chr1 | 85.707 | passed | 0 |
| wf016 | single-trait-model | Single-trait SSGBLUP (AI-REML) on simulated_population_chr1_part + pedigree | 54.853 | passed | 0 |
| wf017 | linear-mixed-model | GBLUP (EMAI-REML) on reference chr1 + GWAStestphenotype (Weight trait) | 9.677 | passed | 0 |
| wf018 | snp-effect | RR-BLUP SNP effect estimation on simulated chr1 data | 5.977 | passed | 0 |
| wf019 | gxe-model | GxE interaction model on simulated chr1 data | 5.262 | passed | 0 |
| wf020 | multi-trait-model | Multi-trait variance component analysis on simulated chr1 data | 59.515 | passed | 0 |
| wf021 | repeated-records-model | Repeated records model on simulated longitudinal chr1 data | 10.917 | passed | 0 |
| wf022 | gebv-gprs-prediction | GEBV prediction from SNP effects on simulated chr1 data | 8.780 | passed | 0 |
| wf023 | qc | QC diagnostic on Guangdong-pig 98 WGS samples | 13.889 | passed | 0 |
| wf024 | roh | ROH detection + F_ROH inbreeding on Guangdong-pig 98 WGS samples | 2.917 | passed | 0 |
| wf025 | nj-tree | NJ phylogenetic tree from IBS distance on Guangdong-pig 98 WGS samples | 4.078 | passed | 0 |
| wf026 | internal-impute | V1 multi-chip internal imputation: gdp_chip_20 (target) + gdp_ref_78 (ref), chr1 | 926.497 | passed | 0 |

## 输出物

本基准测试产出以下交付物：

- 26 个任务级日志文件，位于 `results/wfNNN-<skill>/wfNNN-<skill>.log`（完整 stdout + stderr）。
- 26 个任务级 JSON 文件，位于 `results/wfNNN-<skill>/wfNNN-<skill>.json`（elapsed_seconds、status、rc、input_files、output_files、error_message）。
- 26 个任务级输出目录，位于 `results/wfNNN-<skill>/<skill outputs>`（每个技能实际产出的文件）。
- 1 个汇总 JSON，位于 `workflow_execution_results.json`（最终交付物，由 `lib/aggregate_results.py` 生成）。
- 2 个 README 文件：本英文 README 与简体中文镜像 `README-zh.md`。

## 图

本基准测试为 20 个含可视化输出的 wf 文件夹生成了 20 张 Nature 风格多面板组图（含 A/B/C 角标）。每个图目录 `figures/wfNNN-<skill>/` 包含：
- `results/` — 从 `results/wfNNN-<skill>/` 复制的 PDF/PNG 文件（扁平化）
- `combined_figure.pdf` — 合成的多面板 PDF（含 A/B/C 角标）
- `figure_caption.txt` — 发表级图注（英文）

| 图号 | wf 文件夹 | 技能 | 面板数 | 图注标题 |
|------|-----------|------|--------|----------|
| S1 | wf004-gcta-gwas | gcta-gwas | 2 | GCTA MLMA GWAS 结果（Manhattan + Q-Q） |
| S2 | wf005-genomic-selection | genomic-selection | 1 | 6 模型基因组选择比较 |
| S3 | wf009-pca | pca | 5 | 群体结构 PCA（碎石图 + PC1-2-3 + 密度图） |
| S4 | wf010-ld | ld | 1 | LD 衰减曲线（PopLDdecay，全染色体，300 kb） |
| S5 | wf011-ld-score | ld-score | 2 | LD Score 分布 + LD score vs MAF |
| S6 | wf012-allele-genotype-frequency | allele-genotype-frequency | 2 | 等位基因 & 基因型频率（总览 + 详情） |
| S7 | wf013-homozygosity-heterozygosity | homozygosity-heterozygosity | 3 | 纯合/杂合（染色体 + 样本 + 质量） |
| S8 | wf014-inbreeding-relationship | inbreeding-relationship | 2 | 近交系数 F 分布 + 方法比较 |
| S9 | wf015-relationship-matrix | relationship-matrix | 2 | GRM 分布 + 热图总览 |
| S10 | wf016-single-trait-model | single-trait-model | 3 | SSGBLUP（EBV + 方差组分 + 饼图） |
| S11 | wf017-linear-mixed-model | linear-mixed-model | 1 | GBLUP（EMAI-REML）分析 |
| S12 | wf018-snp-effect | snp-effect | 1 | RR-BLUP SNP 效应估计 |
| S13 | wf019-gxe-model | gxe-model | 3 | GxE（森林图 + 随机效应 + 方差条形图） |
| S14 | wf020-multi-trait-model | multi-trait-model | 7 | 多性状（ANOVA + beta 森林图 + r_g + h² + 随机 + r_e + VC） |
| S15 | wf021-repeated-records-model | repeated-records-model | 3 | 重复记录模型（EBV + 固定效应 + VC） |
| S16 | wf022-gebv-gprs-prediction | gebv-gprs-prediction | 1 | GEBV/GPRS 预测 |
| S17 | wf023-qc | qc | 4 | QC 诊断（call rate + 杂合 + 缺失 + PCA） |
| S18 | wf024-roh | roh | 2 | ROH（F_ROH 条形图 + 长度分布） |
| S19 | wf025-nj-tree | nj-tree | 1 | NJ 系统发育树（IBS 距离） |
| S20 | wf026-internal-impute | internal-impute | 3 | 插补 QC（dr² + info score + 直方图） |

合成脚本：`lib/compose_figures.py`（pypdfium2 + pikepdf + reportlab；通过 XObject 放置实现矢量 PDF 嵌入——源 PDF 以矢量形式保留，不栅格化；A/B/C 面板角标通过 reportlab overlay 经 pikepdf 叠加，无背景框）。
图注来源：`lib/figure_captions.py`（20 份发表级英文图注）。
清单文件：`lib/figures_manifest.json`（源图文件完整清单）。

## 局限性与注意事项

- 测试数据被有意控制在较小规模，以使全部 26 个技能的扫描在数分钟内完成。耗时数据不代表生产规模；在全队列 WGS 或大型 SNP 芯片上的生产运行将显著更慢，并可能触及小型测试数据无法覆盖的代码路径。
- 本基准测试不修改 PigGPA 底层技能脚本。全部 26 个技能以原始实现进行测试。
