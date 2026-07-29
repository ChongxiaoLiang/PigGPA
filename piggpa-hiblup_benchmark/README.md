# piggpa-G vs HIBLUP Benchmark — Complete Comparison Package

A comprehensive, reproducible benchmark comparing **piggpa-G** (a Python-based genomic prediction toolkit) against **HIBLUP** (the C++ v1.6.0 industry-standard tool) on identical simulated pig population data. Every numerical result in this package was produced by actually running both tools on the same inputs — nothing is fabricated or hand-edited.

The benchmark spans **14 benchmark directories (14 compared modules)** covering the full genomic-prediction pipeline, from allele frequency, homozygosity/heterozygosity, and relationship-matrix construction through BLUP prediction, variance-component estimation, GxE modelling, GEBV prediction, and LD analyses. Four core tasks (T3 relationship matrix, T7 BLUP, T8 single-trait, T9 repeated-records) are directly compared metric-by-metric; the remaining ten tasks carry paired piggpa-G and HIBLUP result sets.

## Software Versions

| Software | Version | Language | Notes |
|----------|---------|----------|-------|
| piggpa-G | Internal development version | Python 3.13 | Genomic prediction toolkit |
| HIBLUP | v1.6.0 (2025-09-29 Release) | C++ | Industry-standard tool, www.hiblup.com |

## Dataset

| Property | Value |
|----------|-------|
| Dataset | In-silico simulated pig population |
| Total individuals | 10,000 |
| Total SNPs | 100,000 |
| Chromosomes | 18 autosomes |
| Test subset | 1,000 individuals + chr1 (5,556 SNPs) |
| Input format | PLINK bed/bim/fam |
| Input data source | `/public/share/likui/hanyu/testdata/In-silico-data/` |

## Directory Structure

```
piggpa-hiblup_benchmark/
├── README.md                                    # This document (English)
├── .gitignore                                   # Excludes __pycache__, *.pyc, *.bed
├── README-zh.md                                 # This document (Chinese)
├── benchmark/                                   # Benchmark comparison data
│   ├── README.md                                # Benchmark description (English)
│   ├── README-zh.md                             # Benchmark description (Chinese)
│   ├── relationship_matrix/                     # T3: Relationship matrix comparison
│   ├── blup_prediction/                         # T7: BLUP breeding value prediction
│   ├── single_trait_model/                      # T8: Single-trait variance component estimation
│   ├── repeated_records_model/                  # T9: Repeated records model
│   ├── allele_frequency/                        # T1: Allele frequency calculation
│   ├── homozygosity_heterozygosity/             # T2: Homozygosity and heterozygosity calculation
│   ├── inbreeding_coefficient/                  # T5: Inbreeding coefficient + kinship
│   ├── pca/                                     # T6: PCA genetic structure analysis
│   ├── multi_trait_model/                       # T10: Multi-trait variance component estimation
│   ├── gxe_model/                               # T12: GxE interaction model
│   ├── snp_effect/                              # T11: SNP effect calculation (prerequisite for T16)
│   ├── gebv_prediction/                         # T16: GEBV genomic estimated breeding values
│   ├── ld_calculation/                          # T18: Linkage disequilibrium calculation
│   └── ld_score_regression/                     # T20: LD Score regression
├── scripts/                                     # piggpa-G source scripts (8 .py)
│   ├── relationship_matrix_construction.py      # T3: Relationship matrix (default HRM = 0.95*G_adj + 0.05*A, HIBLUP HA equivalent)
│   ├── single_trait_model.py                    # T8: Single-trait model (AI-REML/EM-REML/HE)
│   ├── repeated_records_model.py                # T9: Repeated records model
│   ├── model_comparison.py                      # T7: Model comparison
│   ├── LM_model.py                              # T7: Linear Model
│   ├── BLUP_model.py                            # T7: BLUP
│   ├── PBLUP_model.py                           # T7: Pedigree BLUP
│   ├── GBLUP_model.py                           # T7: Genomic BLUP
│   └── SSBLUP_model.py                          # T7: Single-step BLUP
├── hiblup_scripts/                              # 14 HIBLUP wrapper scripts (reproduce every HIBLUP command)
│   ├── hiblup_allele_frequency.sh               # T1/T2
│   ├── hiblup_relationship_matrix.sh            # T3
│   ├── hiblup_inbreeding_coefficient.sh         # T5
│   ├── hiblup_pca.sh                            # T6
│   ├── hiblup_blup_prediction.sh                # T7
│   ├── hiblup_single_trait_model.sh             # T8
│   ├── hiblup_repeated_records_model.sh         # T9
│   ├── hiblup_multi_trait_model.sh              # T10
│   ├── hiblup_gxe_model.sh                      # T12
│   ├── hiblup_snp_effect_calculation.sh         # T11 (prerequisite for T16)
│   ├── hiblup_gebv_prediction.sh                # T16
│   ├── hiblup_ld_calculation.sh                 # T18
│   └── hiblup_ld_score_regression.sh            # T20
├── figures/                                     # CNS-quality correlation figure + plot script
│   ├── fig_correlation_overview.pdf             # 3×5 multi-panel: 11 scatter + 4 bar charts
│   └── scripts/
│       └── plot_cns_figures.py                  # Python script generating the correlation figure
└── unified_testdata/                            # 16 shared input files for all benchmark tasks
    ├── simulated_population.bed                 # PLINK bed — EXCLUDED (238 MB, exceeds GitHub 100 MB limit)
    ├── simulated_population.bim                 # PLINK bim (SNP map, chr1: 5,556 SNPs)
    ├── simulated_population.fam                 # PLINK fam (1,000 individuals)
    ├── phenotypes.txt                           # Phenotype file
    ├── phenotype_train_samples.csv              # Phenotype + training sample mapping
    ├── train_samples.txt                        # Training sample ID list
    ├── pred_samples.txt                         # Prediction sample ID list
    ├── keep_samples.txt                         # Samples to keep
    ├── keep_1000_samples.txt                    # 1,000-sample keep list
    ├── chr1_snps.txt                            # chr1 SNP list (5,556 SNPs)
    ├── extract_snps.txt                         # SNP extraction list
    ├── simulated_phenotypes_multi_trait.txt     # Multi-trait phenotypes
    ├── snp_effect.snpeff                        # SNP effect file
    ├── gblup_pred.bv                            # GBLUP predicted breeding values
    ├── gblup_train.rand                         # GBLUP random effects
    └── gblup_train.vars                         # GBLUP variance components
```

## Task Coverage — 14 Benchmark Directories (14 Compared Modules)

| # | Task | Functional Directory | Comparison Status |
|---|------|---------------------|-------------------|
| 1 | T1 | `benchmark/allele_frequency/` | Paired results |
| 2 | T2 | `benchmark/homozygosity_heterozygosity/` | Paired results (Pearson r=1.0000) |
| 3 | T3 | `benchmark/relationship_matrix/` | **Directly compared** (3 matrix pairs, r/MSE) |
| 4 | T5 | `benchmark/inbreeding_coefficient/` | Paired results |
| 5 | T6 | `benchmark/pca/` | Paired results |
| 6 | T7 | `benchmark/blup_prediction/` | **Directly compared** (GBLUP/SSBLUP/BLUP Cor_BV) |
| 7 | T8 | `benchmark/single_trait_model/` | **Directly compared** (AI-REML/EM-REML/HE) |
| 8 | T9 | `benchmark/repeated_records_model/` | **Directly compared** (convergence + variance components) |
| 9 | T10 | `benchmark/multi_trait_model/` | Paired results |
| 10 | T11 | `benchmark/snp_effect/` | Paired results; prerequisite for T16 |
| 11 | T12 | `benchmark/gxe_model/` | Paired results |
| 12 | T16 | `benchmark/gebv_prediction/` | Paired results |
| 13 | T18 | `benchmark/ld_calculation/` | Paired results |
| 14 | T20 | `benchmark/ld_score_regression/` | Paired results |

Each functional directory (except the four directly-compared core tasks) follows the structure `{name}/piggpa_G/` + `{name}/hiblup/`, holding the paired outputs from both tools.

## Key Findings Summary

### T3 Relationship Matrix (`benchmark/relationship_matrix/`)

| Matrix Pair | Level | r | MSE |
|-------------|-------|-----|-----|
| PRM vs PA | IDENTICAL | **1.0000** | 0.0 |
| GRM_VanRaden vs GA | GOOD | **0.9990** | 1.034e-05 |
| HRM vs HA | GOOD | **0.9991** | 5.158e-06 |

### T7 BLUP Breeding Value Prediction (`benchmark/blup_prediction/`)

| Model | piggpa-G Cor_BV | HIBLUP Cor_BV | Difference | Comparable |
|-------|----------------|---------------|------------|------------|
| GBLUP | 0.05404707 | 0.0538 | **0.000247** | Yes ✓ |
| SSBLUP | 0.05404706 | 0.0539 | **0.000147** | Yes ✓ |
| BLUP | 0.05404707 | 0.0538 | N/A | No (GA vs PA) |

### T8 Single-Trait Variance Components (`benchmark/single_trait_model/`)

| Method | Tool | V(G) | V(e) | h² | logL | Iters | Converged |
|--------|------|------|------|-----|------|-------|-----------|
| AI-REML | HIBLUP | **0.0000** | 0.9912 | 3.27e-07 | -500.58 | 11 | Yes |
| AI-REML | piggpa-G | **0.000000** | 1.018681 | 0.000000 | -500.7101 | 4 | Yes |
| EM-REML | piggpa-G | 0.030117 | 0.967842 | 0.030178 | -500.6819 | 200 | No |
| HE | piggpa-G | 0.000000 | 0.000000 | 0.500000 | N/A | 1 | Degenerate |

### T9 Repeated Records Model (`benchmark/repeated_records_model/`)

| Metric | piggpa-G | HIBLUP |
|--------|----------|--------|
| Converged | **Yes** (9 iters) | **No** (20 iters) |
| V(ID) | 85.5030 | 79.7363 |
| V(GA) | 19.5169 | 24.8079 |
| V(e) | **14.1266** | **14.1266** |
| h²(GA) | 0.1638 | 0.2090 |

## CNS Publication Figures (`figures/`)

One Nature/Cell/Science-quality vector PDF figure generated from the verified benchmark data:

| File | Description |
|------|-------------|
| `fig_correlation_overview.pdf` | 3×5 multi-panel (24×15 in): 11 scatter plots (piggpa-G vs HIBLUP per module; y=x reference line on Allele Frequency and SNP Effect panels only; r/ρ annotations) + 4 bar charts (T3 matrix r, T8/T9 variance components, T20 LD Score regression). High-consistency modules (r≥0.99): all 11 scatter modules. |

Reproduction: `conda run -n py_analysis python figures/scripts/plot_cns_figures.py`

## Key Algorithms

1. **GRM**: piggpa-G implements VanRaden ($G = ZZ'/\sum 2pq$) and Yang methods; HIBLUP uses Su-normalized VanRaden.
2. **Hybrid Relationship Matrix**:
   - piggpa-G default: $HRM = 0.95 \times G_{adj} + 0.05 \times A$ (equivalent to HIBLUP HA), where $G_{adj} = 0.999G + 0.001I$
   - HIBLUP: $HA = 0.95 \times G_{adj} + 0.05 \times A$, where $G_{adj} = 0.999G + 0.001I$
3. **Variance Component Estimation**: piggpa-G implements AI-REML, EM-REML (`--em-max-iter` default 200), and HE regression (with ridge regularization); HIBLUP uses AI-REML (max_iter=20).

## Data Provenance

All results were generated by actually running piggpa-G and HIBLUP on identical input data. The four directly-compared core tasks (T3/T7/T8/T9) use a unified phenotype file (`phenotype_hiblup.txt`, h²≈0) for T8. T9 results use the canonical AI-REML converged output (9 iterations). Every numerical claim can be traced back to specific `.log`/`.vars`/`.csv` files with line numbers.

**Large file exclusion**: Full-precision relationship matrix CSV files (1-2GB each) and `simulated_population.bed` (238 MB, exceeds GitHub 100 MB limit) are excluded. The accompanying `.bim` and `.fam` files are retained. Original paths are documented in the benchmark README.
