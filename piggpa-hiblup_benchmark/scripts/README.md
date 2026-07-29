# piggpa-G Source Scripts

This directory contains the source scripts of piggpa-G used in the piggpa-G vs HIBLUP benchmark. All scripts produced the benchmark results in `../benchmark/`.

The full benchmark covers **14 benchmark directories** (14 compared modules; see `../benchmark/`), but the scripts in this directory only correspond to the **4 directly-compared modules** where piggpa-G and HIBLUP outputs were placed side-by-side (T3 relationship matrix, T7 BLUP prediction, T8 single-trait model, T9 repeated records). Wrapper scripts for the other 10 HIBLUP-only functions are in `../hiblup_scripts/`, and the piggpa-G source scripts for those functions live outside this upload package (see "Other piggpa-G Scripts" below).

## Directory Structure

```
scripts/
├── README.md                              # This document (English)
├── README-zh.md                           # This document (Chinese)
├── relationship_matrix_construction.py    # Relationship matrix construction (PRM/GRM/HRM, default HRM=0.95*G_adj+0.05*A)
├── single_trait_model.py                  # Single-trait variance component estimation (AI-REML/EM-REML/HE)
├── repeated_records_model.py              # Repeated records model (permanent environment + additive genetic)
├── model_comparison.py                    # 5-model BLUP comparison driver (LM/BLUP/PBLUP/GBLUP/SSBLUP)
├── LM_model.py                            # Linear Model (LM) — no relationship matrix
├── GBLUP_model.py                         # Genomic BLUP — uses GA
├── PBLUP_model.py                         # Pedigree BLUP — uses PA
└── SSBLUP_model.py                        # Single-step BLUP — uses HA (0.95G_adj + 0.05A)
```

**Note**: `BLUP_model.py` is not present because, in piggpa-G's implementation, the "BLUP" entry in the 5-model comparison reuses the GBLUP code path with GA (genomic relationship matrix). This is reflected in the `relationship_matrix` column of `blup_prediction/model_comparison.csv`.

---

## Script-to-Task Mapping

| Script | Benchmark Task | Output Directory (functional module) |
|--------|---------------|------------------|
| `relationship_matrix_construction.py` | T3 | `../benchmark/relationship_matrix/` |
| `model_comparison.py` + 4 model scripts (`LM_model.py`, `PBLUP_model.py`, `GBLUP_model.py`, `SSBLUP_model.py`) | T7 | `../benchmark/blup_prediction/` |
| `single_trait_model.py` | T8 | `../benchmark/single_trait_model/` |
| `repeated_records_model.py` | T9 | `../benchmark/repeated_records_model/` |

The four T7 sub-model scripts (`LM_model.py`, `PBLUP_model.py`, `GBLUP_model.py`, `SSBLUP_model.py`) are invoked by `model_comparison.py` and all write their outputs under `../benchmark/blup_prediction/`.

---

## Key Algorithms

### relationship_matrix_construction.py
- **PRM** (Pedigree Relationship Matrix): Henderson recursive algorithm
- **GRM** (Genomic Relationship Matrix): VanRaden method ($G = ZZ'/\sum 2pq$) and Yang method
- **HRM** (Hybrid Relationship Matrix): default formula $HRM = 0.95 \cdot G_{adj} + 0.05 \cdot A$ (equivalent to HIBLUP's HA), where $G_{adj} = 0.999 \cdot G + 0.001 \cdot I$
- The HIBLUP-compatible formula is the default HRM output. `HRM.csv` directly uses the HIBLUP-compatible formula.

### single_trait_model.py
- **AI-REML**: Average Information REML (default, max_iter=20)
- **EM-REML**: Expectation-Maximization REML (`--em-max-iter`, default 200)
- **HE regression**: Henderson's regression with ridge regularization
- Internally constructs HA using the same 0.95/0.05 formula as HIBLUP

### repeated_records_model.py
- Model: `weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e`
- AI-REML with permanent environmental effect (ID) and additive genetic effect (GA)
- Converges in 9 iterations (vs HIBLUP's 20 iterations non-converged)

### model_comparison.py + model scripts
- Drives 5 BLUP models: LM, BLUP (=GBLUP with GA), PBLUP, GBLUP, SSBLUP
- Outputs `model_comparison.csv` with `relationship_matrix` column for explicit comparability annotation
- Training set: ID 1001-2000; Validation set: ID 5001-6000

---

## Dependencies

- **Python**: 3.13+
- **NumPy**, **pandas**, **SciPy**: numerical computation
- **matplotlib**: figure generation
- **HIBLUP** (v1.6.0): reference tool, not required to run piggpa-G scripts

## Usage

Each script is self-contained and can be run independently. For exact invocation commands, refer to the run logs in `../benchmark/{relationship_matrix,blup_prediction,single_trait_model,repeated_records_model}/`.

## Other piggpa-G Scripts

The benchmark covers 14 benchmark directories (14 compared modules), but only the 4 directly-compared modules (T3, T7, T8, T9) ship their piggpa-G Python source in this `scripts/` directory. The remaining 10 functions (allele_frequency, homozygosity_heterozygosity, inbreeding, PCA, multi_trait, GxE, GEBV, LD, LD_score_regression, snp_effect) were executed HIBLUP-side only; their piggpa-G equivalents are available at:

```
/public/share/likui/hanyu/灵活版纯脚本/function1-function14/
```

The HIBLUP wrapper shell scripts for all 14 functions (including the 10 HIBLUP-only ones) are in `../hiblup_scripts/`.

## Data Provenance

These scripts produced all piggpa-G results in `../benchmark/` for the four directly-compared tasks (T3, T7, T8, T9).
