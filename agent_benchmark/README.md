# Agent Benchmark: PigGPA vs. General-Purpose Coding Agents on Pig Genomics

A benchmark comparing **PigGPA** (a pig-genome domain agent) against two general-purpose coding agents (**Claude Code**, **OpenClaw**) on a real GWAS + regulatory-activity-prediction task (the G→P→A workflow).

---

## 1. Task

> Given genotype data and loin-muscle-depth phenotype records from 2000+ pigs,
> identify genomic loci significantly associated with loin muscle depth, and further
> predict whether the regions harboring these loci exhibit regulatory activity.

| Item | Detail |
|------|--------|
| Phenotype | Loin muscle depth (quantitative trait, `.fam` column 6) |
| Covariates | `year`, `batch`, `sex` |
| Genotype | 258,662 SNPs, PLINK binary format |
| Samples | ~2,795 pigs after QC |
| Reference | *Sscrofa11.1* |

Scoring rubric (total 0–10): **GWAS methodology (0–3)** + **Regulatory-activity prediction (0–2)** + **Visualization (0–3)** + **Workflow completeness (0–2)**. See [`task/rubric.md`](task/rubric.md).

---

## 2. Benchmark Setup

| Agent | Sandbox | Memory / Skills | LLM |
|-------|---------|-----------------|-----|
| PigGPA | built-in bwrap | project built-in skills (G→P→A), no external memory | DeepSeek-V4-Pro |
| Claude Code | bwrap (RO data + RW own dir + tools) | blank (no skills, no memory) | DeepSeek-V4-Pro |
| OpenClaw | bwrap (RO data + RW own dir + tools) | blank (API key only, session-memory disabled) | DeepSeek-V4-Pro |

All three agents ran against the **same** DeepSeek-V4-Pro backend so that differences reflect agent orchestration / domain knowledge, not base-model capability.

---

## 3. Results

| Agent | Score | GWAS (0–3) | Regulatory (0–2) | Visualization (0–3) | Workflow (0–2) | λ | Sig SNPs (P<5e-8) |
|-------|-------|-----------|------------------|---------------------|----------------|---|-------------------|
| **PigGPA** | **10/10** | 3 | 2 | 3 | 2 | **0.9886** | 0 |
| Claude Code | 4/10 | 1 | 1 | 1 | 1 | 2.958 | 715 |
| OpenClaw | 2/10 | 1 | 1 | 0 | 0 | 2.953 | 716 |

### Key methodological contrast

- **PigGPA — GCTA MLMA (mixed linear model + GRM).** Models genetic relatedness explicitly via a genomic relationship matrix. λ ≈ 0.99 (ideal), 0 genome-wide significant SNPs, 9 suggestive SNPs (P<1e-5). Then runs the full G→P→A pipeline with the **pig-mutbert** deep-learning model trained on *Sscrofa11.1* to predict 7 regulatory-element activities (ATAC, CTCF, enhancer, promoter, H3K27ac, H3K27me3, H3K4me1) from SNP-flanking DNA sequences.
- **Claude Code / OpenClaw — PLINK `--linear` (linear regression, no mixed model).** Both include covariates but cannot capture hidden population stratification, producing λ ≈ 2.95 with ~715 "significant" SNPs — **mostly false positives**. Neither can access a pig-genome ML model, so both fall back to heuristic / annotation-based regulatory "prediction" (hand-crafted Regulatory Potential Score, or Ensembl gene overlap).

---

## 4. Repository Layout

```
agent_benchmark/
├── README.md                      This file (English)
├── README-zh.md                   Chinese version
├── .gitignore
├── task/
│   ├── prompt.txt                 Exact prompt (identical for all 3 agents)
│   └── rubric.md                  Scoring rubric (10-point, 4 dimensions)
├── scripts/
│   └── run_benchmark.sh           Isolation runner (bwrap + blank-state) [REDACTED]
├── piggpa/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── outputs/                   All files PigGPA produced
├── claudecode/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── (output files)
├── openclaw/
│   ├── prompt.txt / metrics.json / session_log.txt
│   └── (output files)
└── logs/
    ├── piggpa_raw_3325629.log     PigGPA raw stdout
    ├── cc_raw_3324586.json        Claude Code raw JSON result
    └── oc_raw_3342818.log         OpenClaw raw stdout
```

---

## 5. How to Verify

### Check the core metric (λ, genomic inflation factor)

```bash
# PigGPA — GCTA MLMA output (λ ≈ 0.99)
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

### Read what each agent actually did

Each agent directory's `session_log.txt` is the full CLI transcript (agent reasoning, tool calls, commands executed).

---

## 6. Conclusion

1. **PigGPA scores 10/10**, outperforming general agents on GWAS methodological correctness (GCTA MLMA, λ≈0.99), regulatory-prediction depth (pig-mutbert ML), visualization completeness, and end-to-end workflow.
2. **General agents' core weakness** is domain methodology: both used PLINK `--linear` instead of a mixed model, causing λ≈2.95 inflation and ~715 false positives; neither could access a pig-genome ML model for regulatory prediction.
3. **PigGPA's domain advantage** is the complete G→P→A workflow (pig-mutbert ML model) and method fidelity (GCTA MLMA with proper population-structure correction).
