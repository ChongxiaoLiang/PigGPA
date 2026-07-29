# Agent 基准测试：PigGPA 与通用编程 Agent 在猪基因组学任务上的对比

将 **PigGPA**（猪基因组领域 Agent）与两个通用编程 Agent（**Claude Code**、**OpenClaw**）在真实的 GWAS + 调控活性预测任务（G→P→A 流程）上进行对比的基准测试。

---

## 1. 任务

> 给定 2000+ 头猪的基因型数据和眼肌深度（loin muscle depth）表型记录，
> 鉴定与眼肌深度显著关联的基因组位点，并进一步预测这些位点所在区域是否具有调控活性。

| 项目 | 详情 |
|------|------|
| 表型 | 眼肌深度（定量性状，`.fam` 第 6 列） |
| 协变量 | `year`、`batch`、`sex` |
| 基因型 | 258,662 SNPs，PLINK binary 格式 |
| 样本 | QC 后约 2,795 头猪 |
| 参考基因组 | *Sscrofa11.1* |

评分量表（总分 0–10）：**GWAS 方法学 (0–3)** + **调控活性预测 (0–2)** + **可视化产出 (0–3)** + **流程完整性 (0–2)**，详见 [`task/rubric.md`](task/rubric.md)。

---

## 2. 基准测试设置

| Agent | 沙箱 | 记忆 / Skills | LLM |
|-------|------|--------------|-----|
| PigGPA | 内置 bwrap | 项目自带 Skills（G→P→A），无外部记忆 | DeepSeek-V4-Pro |
| Claude Code | bwrap（RO 数据 + RW 自身目录 + 工具） | 空白（无 Skills、无记忆） | DeepSeek-V4-Pro |
| OpenClaw | bwrap（RO 数据 + RW 自身目录 + 工具） | 空白（仅 API key，session-memory 已禁用） | DeepSeek-V4-Pro |

三个 Agent 均连接**同一个** DeepSeek-V4-Pro 后端，因此差异体现的是 Agent 编排与领域知识，而非基座模型能力。

---

## 3. 结果

| Agent | 评分 | GWAS(0-3) | 调控预测(0-2) | 可视化(0-3) | 流程(0-2) | λ | 显著SNP(P<5e-8) |
|-------|------|-----------|--------------|------------|----------|---|-----------------|
| **PigGPA** | **10/10** | 3 | 2 | 3 | 2 | **0.9886** | 0 |
| Claude Code | 4/10 | 1 | 1 | 1 | 1 | 2.958 | 715 |
| OpenClaw | 2/10 | 1 | 1 | 0 | 0 | 2.953 | 716 |

### 核心方法论对比

- **PigGPA —— GCTA MLMA（混合线性模型 + GRM）。** 通过遗传关系矩阵（GRM）显式建模样本间遗传相关性。λ ≈ 0.99（理想），0 个全基因组显著 SNP，9 个提示性 SNP（P<1e-5）。随后执行完整 G→P→A 流程，使用基于 *Sscrofa11.1* 训练的 **pig-mutbert** 深度学习模型，从 SNP 侧翼 DNA 序列预测 7 种调控元件活性（ATAC、CTCF、enhancer、promoter、H3K27ac、H3K27me3、H3K4me1）。
- **Claude Code / OpenClaw —— PLINK `--linear`（线性回归，无混合模型）。** 两者均纳入协变量，但无法捕获隐含的群体分层，导致 λ ≈ 2.95，约 715 个"显著"SNP —— **大部分为假阳性**。两者均无法访问猪基因组 ML 模型，只能回退到启发式 / 注释性的"调控预测"（手工设计的调控潜力评分 RPS，或 Ensembl 基因重叠注释）。

---

## 4. 仓库结构

```
agent_benchmark/
├── README.md                      本文件英文版
├── README-zh.md                   本文件（中文）
├── .gitignore
├── task/
│   ├── prompt.txt                 精确 prompt（三个 Agent 完全相同）
│   └── rubric.md                  评分量表（10 分制，4 个维度）
├── scripts/
│   └── run_benchmark.sh           隔离运行脚本（bwrap + 空白状态）[已脱敏]
├── piggpa/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── outputs/                   PigGPA 的全部产出文件
├── claudecode/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── (产出文件)
├── openclaw/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── (产出文件)
└── logs/
    ├── piggpa_raw_3325629.log     PigGPA 原始 stdout
    ├── cc_raw_3324586.json        Claude Code 原始 JSON 结果
    └── oc_raw_3342818.log         OpenClaw 原始 stdout
```

---

## 5. 如何验证

### 检查核心指标（λ，基因组膨胀因子）

```bash
# PigGPA — GCTA MLMA 输出（λ ≈ 0.99）
python3 -c "
import pandas as pd, numpy as np
df = pd.read_csv('piggpa/outputs/lmd_gcta.mlma', sep='\t')
p = df['p'].dropna()
chi = -2*np.log(p)
print('PigGPA  λ =', round(np.median(chi)/0.4549, 4), ' n_sig_5e8 =', (p<5e-8).sum())
"

# Claude Code / OpenClaw — PLINK .assoc.linear
python3 -c "
import pandas as pd, numpy as np
df = pd.read_csv('claudecode/lmd.gwas.assoc.linear', sep=r'\s+')
p = pd.to_numeric(df['P'], errors='coerce').dropna()
chi = -2*np.log(p)
print('Claude  λ =', round(np.median(chi)/0.4549, 4), ' n_sig_5e8 =', (p<5e-8).sum())
"
```

### 阅读各 Agent 实际做了什么

各 Agent 目录下的 `session_log.txt` 是完整 CLI 记录（Agent 推理、工具调用、执行的命令）。

---

## 6. 结论

1. **PigGPA 以 10/10 满分胜出**，在 GWAS 方法学正确性（GCTA MLMA，λ≈0.99）、调控预测深度（pig-mutbert 真实 ML 预测）、可视化完整性和端到端流程方面均显著优于通用 Agent。
2. **通用 Agent 的核心短板**在于领域方法论：两者均使用 PLINK `--linear` 而非混合模型，导致 λ≈2.95 严重膨胀和约 715 个假阳性；调控预测方面均无法访问猪基因组 ML 模型，只能采用启发式或注释性替代方案。
3. **PigGPA 的领域优势**体现在完整 G→P→A 工作流执行能力（使用 pig-mutbert ML 模型）和领域方法保真度（GCTA MLMA 正确校正群体结构）。
