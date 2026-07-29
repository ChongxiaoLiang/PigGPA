# PigGPA 基准测试

**PigGPA**（Genomic breeding；Prediction；Annotation，猪基因组研究 AI Agent 系统）的可复现证据与基准测试数据。目录布局参照 BOLE 基准测试（ChongxiaoLiang/BOLE），但此处所有脚本与 Agent 调用均执行 **PigGPA 自有的子技能**与 Agent 循环。当前基准测试包含三个子集：全部 26 个 PigGPA-G 技能的端到端工作流执行、100 条自然语言查询的意图解析、5 个对抗场景的错误处理。

> **重要归属说明**：原始 BOLE 基准测试（https://github.com/ChongxiaoLiang/BOLE）是独立项目。BOLE 的 `bole_benchmark/` 目录中的数据是 *BOLE 的* 结果，并非 PigGPA 的结果。本 `piggpa_benchmark/` 目录包含的是 PigGPA 自身的评测数据，由运行 PigGPA 自有子技能脚本与 `piggpa chat` Agent 循环产生。两者不可混为一谈。

## 基准测试子集

| 子集 | 描述 | 运行模式 |
|--------|-------------|----------|
| `benchmark/workflow_execution/` | 在代表性测试数据上对全部 26 个 PigGPA-G 子技能进行端到端冒烟测试（详见 [workflow_execution/README-zh.md](benchmark/workflow_execution/README-zh.md)） | 真实运行，26/26 通过（100.0%） |
| `benchmark/intent_parsing/` | 100 条自然语言查询 × 8 个意图类别 × 3 个难度级别，经 `piggpa chat -q` 派发 | 真实运行，100/100 PASS（100.0%） |
| `benchmark/error_handling/` | 5 个错误场景（缺失输入、非法参数、意图歧义、非法路径、参数不兼容），经 `piggpa chat -q` 派发 | 真实运行，5/5 PASS（100%） |

## 受测 PigGPA 子技能

`workflow_execution` 子集端到端测试全部 26 个 PigGPA-G 子技能。完整的任务级清单（技能、脚本路径、耗时、状态、退出码）见 [benchmark/workflow_execution/README-zh.md](benchmark/workflow_execution/README-zh.md)。26 个技能覆盖：

- **导入与质控**：genotype-import、genotype-qc、qc
- **GWAS**：plink-gwas、gcta-gwas
- **群体遗传学**：admixture-analysis、pca、ld、ld-score、ld-pruning、allele-genotype-frequency、homozygosity-heterozygosity、inbreeding-relationship、relationship-matrix、roh、nj-tree
- **育种值与选择**：genomic-selection、single-trait-model、linear-mixed-model、multi-trait-model、repeated-records-model、gxe-model、snp-effect、gebv-gprs-prediction、heritability-analysis
- **插补**：internal-impute

## 与 BOLE 的差异

PigGPA 脚本与 BOLE 脚本的差异：

1. 环境变量重命名：`ALBA_BIN_DIR/ALBA_APP_BIN/ALBA_OUTDIR/ALBA_INDIR` → `APP_BIN`
2. 新增 `--app-bin` 参数
3. 新增相对路径告警
4. `plink-gwas` 子技能（含 `PLINK-GWAS-linear.sh`）为 PigGPA 独有（BOLE 无对应实现）
5. `relationship-matrix` Python 脚本为 PigGPA 独有

## 依赖

最新一次运行使用的工具版本与绝对路径见 [benchmark/workflow_execution/README-zh.md](benchmark/workflow_execution/README-zh.md#测试环境)。最低要求：

- **PLINK**（1.9+）、**GCTA**（1.94+）、**ADMIXTURE**、**bcftools**、**PopLDdecay**
- **R**（4.0+）含 `data.table`、`qqman`、`ggplot2`、`BGLR` 等
- **Python**（3.8+）含 `pandas`、`numpy`、`scipy`、`matplotlib`、`pypdfium2`、`pikepdf`、`reportlab`

## 可复现性

`benchmark/*/` 中的所有结果均通过实际调用 PigGPA 子技能脚本或 `piggpa chat` Agent 循环产生。无任何结果 JSON 为手写，均由各子集的运行脚本从真实 stdout/stderr 提取。

- `workflow_execution`：`benchmark/workflow_execution/scripts/` 中的任务级封装脚本可重建 `results/wfNNN-<skill>/`（日志 + 输出 + 任务级 JSON）；`lib/aggregate_results.py` 可基于任务级 JSON 重新生成汇总 JSON。
- `intent_parsing`：任务级日志位于 `benchmark/intent_parsing/logs/IP-XXX.log`；汇总位于 `benchmark/intent_parsing/intent_parsing_results.json`。
- `error_handling`：场景级日志位于 `benchmark/error_handling/logs/ERR-00X.log`；汇总位于 `benchmark/error_handling/error_handling_results.json`。

## 真实运行结果

### workflow_execution（26/26 通过，100.0%，2026-07-17）

在小型代表性测试数据上对全部 26 个 PigGPA-G 子技能进行端到端冒烟测试。每个技能由一个任务级封装脚本驱动，测量耗时、捕获 stdout/stderr 到日志文件，并写入任务级 JSON 记录。完整的任务级表格、20 张合成多面板图与方法学见 [benchmark/workflow_execution/README-zh.md](benchmark/workflow_execution/README-zh.md)。

| 指标 | 取值 |
|--------|-------|
| 任务总数 | 26 |
| 通过 | 26 |
| 失败 | 0 |
| 总耗时 | 45 min 15 s（2714.82 s） |
| 工作流执行率 | 100.0% |
| 失败任务 ID | （无） |

要点：

- `wf004 gcta-gwas`（GCTA MLMA，基于 GRM）为最慢的分析任务，耗时 1140.6 s（约 19 min）。
- `wf026 internal-impute`（V1 多芯片内部插补，chr1）为整体最慢任务，耗时 926.5 s（约 15.5 min）。
- `wf005 genomic-selection` 以 3 折交叉验证运行全部 6 个 GS 模型（BayesA/B/C/BRR/BL/GBLUP），耗时 46.4 s。
- 本基准测试为 20 个含可视化输出的 wf 文件夹生成了 20 张 Nature 风格多面板组图（含 A/B/C 角标）。

### intent_parsing（100/100 PASS，100.0%，2026-07-01/02）

通过 `piggpa chat -q "<query>" --max-turns 1 --yolo -m deepseek-v4-pro` 对 100 条自然语言查询进行真实运行，覆盖 8 个意图类别 × 3 个难度级别。运行脚本从 `piggpa chat` 日志（通过 `📚 skill <name>` 行与会话 JSON 工具消息）提取被派发的子技能名，并与 `expected_flow` 比对判定通过/失败。

- **总耗时**：7,693.2 s（约 128 min 13 s，平均 77 s/查询）
- **5 条查询经 `--max-attempts 5` 重试后通过；2 条查询经改写后通过**
- **按类别**：GWAS 28/28，QC 13/13，Heritability 10/10，Population_Structure 10/10，Genomic_Selection 10/10，LD 10/10，Import 9/9，Annotation 10/10
- **按难度**：easy 46/46，medium 35/35，hard 19/19
- **日志**：`benchmark/intent_parsing/logs/IP-XXX.log`（每查询一份，100/100 齐全）
- **结果 JSON**：`benchmark/intent_parsing/intent_parsing_results.json`
- **关键映射**（expected_flow → 实际子技能）：plink-gwas-linear → plink-gwas；geno-qc → genotype-qc；heritability-gcta-pipeline → heritability-analysis；admixture → admixture-analysis；gs-6models → genomic-selection；geno-import → genotype-import；annotation-track → pig-annotation-track

### error_handling（5/5 PASS，100%，2026-07-17）

通过 `piggpa chat -q "<query>" --max-turns 2 --yolo -m deepseek-v4-pro` 进行真实运行。5 个错误场景全部通过（ERR-002 启用重试机制 `max_attempts=5`，在第 1 次尝试即通过）：

| 场景 | 类型 | 耗时 | 状态 | PigGPA 检测到 |
|----------|------|---------|--------|-----------------|
| ERR-001 | missing_input | 85.3 s | PASS | 缺失 `--bfile`，提示用户提供文件路径 |
| ERR-002 | invalid_param | 93.9 s | PASS | MAF=0.8 超出 [0, 0.5] 范围，检测并修正（第 1/5 次尝试通过） |
| ERR-003 | ambiguous_intent | 97.0 s | PASS | 意图歧义，列出可选分析类型 |
| ERR-004 | invalid_path | 139.0 s | PASS | 路径不存在，建议核验 |
| ERR-005 | incompatible_params | 93.0 s | PASS | 性状不匹配，请求确认 |

- **总耗时**：约 508 s（约 8.5 min）
- **日志**：`benchmark/error_handling/logs/ERR-001.log` ~ `ERR-005.log`
- **结果 JSON**：`benchmark/error_handling/error_handling_results.json`（由脚本从真实 stdout/stderr 提取）

> **说明**：所有基准测试结果 JSON 均由脚本从真实 `piggpa chat -q` stdout/stderr 或任务级封装脚本的退出码生成，无任何手写 JSON。原始日志保留在 `logs/` 子目录中，可供独立核验。
