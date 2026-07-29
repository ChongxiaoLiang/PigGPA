# PigGPA Workflow Execution Benchmark

This benchmark performs a comprehensive end-to-end test of all 26 PigGPA-G skills located at `/public/share/likui/liangcx/bole/skills/piggpa-G/` against small representative test data. Each skill was wrapped in a per-task script that measures elapsed time, captures stdout/stderr to a log file, and writes a per-task JSON record. A final orchestrator ran all 26 wrappers in numeric order and aggregated the per-task JSONs into a single results file. Of the 26 tasks, **all 26 passed**, yielding a workflow execution rate of **100.0%**. All 26 skills run end-to-end successfully.

## Purpose

This benchmark exists to verify that every PigGPA-G skill runs end-to-end on real input shapes that mirror production usage, with timing and exit codes captured for reproducibility. The benchmark is a smoke test: it confirms that each skill accepts the documented arguments, locates its inputs, produces the expected output file shapes, and exits cleanly. It is not a performance benchmark at production scale, nor a correctness benchmark against published results.

## Test Environment

| Item | Value |
|------|-------|
| Hostname | ln03 |
| Benchmark timestamp | 2026-07-17T03:41:23+08:00 |
| PigGPA skills root | `/public/share/likui/liangcx/bole/skills/piggpa-G` |
| Test data root | `/public/share/likui/liangcx/bole/testdata` (see `testdata_list.md` for inventory) |
| Results root | `/public/share/likui/liangcx/bole/skills/workflow_execution/results/` |
| Aggregate JSON | `/public/share/likui/liangcx/bole/skills/workflow_execution/workflow_execution_results.json` |
| System under test | PigGPA (bole/skills/piggpa-G) |
| Benchmark owner | PigGPA |

Tool paths:

| Tool | Path |
|------|------|
| PLINK / GCTA / ADMIXTURE / bcftools | `/public/share/likui/liangcx/software/miniconda3/envs/sys_tools/bin` |
| Rscript | `/public/share/likui/liangcx/software/miniconda3/envs/R/bin/Rscript` |
| Python | `/public/share/likui/liangcx/software/miniconda3/envs/py_analysis/bin/python` (numpy 2.4.6, pandas 2.3.3, scipy 1.16.3) |

## Test Data Inputs

The test data is intentionally small to enable a fast end-to-end smoke test. Production-scale data may produce different timing characteristics and may exercise code paths that the small test data does not reach.

| Test data file | Description | Consumed by |
|----------------|-------------|-------------|
| `selectedfortest.vcf` | VCF for PLINK import | wf001 |
| `testps.bed/bim/fam` | 45 samples, ~50K SNPs | wf002, wf006, wf008, wf009 |
| `reference.bed/bim/fam` | 2237 samples, ~258K SNPs | wf003, wf004, wf007, wf010, wf011, wf012, wf013, wf014, wf015, wf017 |
| `reference_50k.bed/bim/fam` | 50K-SNP subset of reference | wf005 |
| `reference_2k.bed/bim/fam` | 2K-SNP subset of reference | wf005 |
| `simulated_population_chr1_part.bed/bim/fam` | 1000 samples, 5556 SNPs (chr1) | wf016, wf018, wf019, wf020, wf021, wf022 |
| `gdp_98samples.bed/bim/fam` | 98 samples, ~1.75M SNPs (WGS) | wf023, wf024, wf025 |
| `gdp_chip_20.bed` | 20-sample target chip | wf026 |
| `gdp_ref_78.bed` | 78-sample reference chip | wf026 |
| `GWAStestphenotype.fam` | 2237 rows (FID IID pheno) | wf003, wf004, wf007, wf017 |
| `snp_effect.snpeff` | TSV (id, add_a1) | wf022 |
| `pca_groups.txt` | Population group annotations | wf009 |
| `pheno_snp_effect.csv` | Phenotypes for SNP-effect tasks | wf018, wf019, wf020, wf021 |
| `chr1_snps.txt` | Chr1 SNP list | wf011 |
| `train_samples.txt`, `pred_samples.txt`, `all_samples.txt` | Sample subset lists | wf005, wf022 |
| `simulated_phenotypes_env_1001_2000.txt` | GxE phenotypes | wf019 |
| `simulated_phenotypes_multi_trait_1001_2000.txt` | Multi-trait phenotypes | wf020 |
| `simulated_phenotypes_repeated_long.txt` | Longitudinal phenotypes | wf021 |
| `pedigree.txt` | Pedigree for SSGBLUP | wf016 |
| `gdp_pop.txt` | Population labels for Guangdong-pig | wf023, wf025 |

## Per-Task Wrappers

Each skill is exercised by a per-task wrapper script named `run_wfNNN_<skill>.sh`. The wrapper template:

- Begins with `ulimit -v unlimited` immediately after the shebang (per workspace rule, to lift the virtual-memory cap and prevent OOM kills on large-memory tasks).
- Resolves `PIGGPA_SKILLS_DIR`, `TESTDATA_DIR`, `RESULTS_DIR`, `APP_BIN`, `RSCRIPT_BIN`, `PYTHON_BIN` from the environment, with sensible defaults pointing at the paths listed in Test Environment above.
- Creates the per-task output directory `results/wfNNN-<skill>/` (e.g. `results/wf001-genotype-import/`).
- Records a start timestamp, invokes the underlying skill script with the documented arguments, and records an end timestamp.
- Computes `WF_ELAPSED` (seconds, 3 decimal places) and `WF_ELAPSED_HUMAN` (auto-formatted as seconds, "X min Y s", or "X h Y min Z s" depending on magnitude).
- Calls `lib/write_task_json.py` to emit `results/wfNNN-<skill>/wfNNN-<skill>.json` containing all metrics (elapsed_seconds, status, rc, input_files, output_files, error_message).
- Returns the underlying skill's exit code, so the orchestrator and the aggregate report reflect the real skill outcome.

Three invocation patterns are used depending on the skill's implementation language:

1. **PLINK/GCTA/ADMIXTURE shell skills**: invoke the skill's `.sh` script with `--genotype_prefix`, `--phenotype_file`, `--output_prefix`, `--threads`, `--app-bin`, etc.
2. **Python skills**: `"$PYTHON_BIN" "$WF_ENTRY" --args...`
3. **R skills**: `"$RSCRIPT_BIN" "$WF_ENTRY" --args...`

All 26 wrapper scripts use the same invocation patterns above.

## How to Run

Run all 26 tasks in numeric order:

```bash
bash run_workflow_execution.sh
```

Run a subset by passing a regex that matches the wf IDs of interest (the orchestrator iterates `run_wfNNN_*.sh` in numeric order and filters by the supplied pattern):

```bash
bash run_workflow_execution.sh wf00[1-9]
```

Skip re-running tasks and just regenerate the aggregate JSON from the existing per-task JSONs:

```bash
python lib/aggregate_results.py
```

## Results Summary

### Aggregate Statistics

| Metric | Value |
|--------|-------|
| Total tasks | 26 |
| Passed | 26 |
| Failed | 0 |
| Total elapsed (seconds) | 2714.82 |
| Total elapsed (human) | 45 min 15 s |
| Workflow execution rate | 100.0% |
| Failed task IDs | (none) |

### Per-Task Results

| wf_id | skill | task | elapsed (s) | status | rc |
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

## Outputs

The benchmark produces the following deliverables:

- 26 per-task log files at `results/wfNNN-<skill>/wfNNN-<skill>.log` (full stdout + stderr).
- 26 per-task JSON files at `results/wfNNN-<skill>/wfNNN-<skill>.json` (elapsed_seconds, status, rc, input_files, output_files, error_message).
- 26 per-task output directories at `results/wfNNN-<skill>/<skill outputs>` (the actual files produced by each skill).
- 1 aggregate JSON at `workflow_execution_results.json` (the final deliverable, produced by `lib/aggregate_results.py`).
- 2 README files: this English README and a Simplified Chinese mirror at `README-zh.md`.

## Figures

The benchmark produces 20 composed multi-panel figures (Nature-style with A/B/C panel labels) for the 20 wf folders that contain visualization outputs. Each figure directory at `figures/wfNNN-<skill>/` contains:
- `results/` — copied PDF/PNG files from `results/wfNNN-<skill>/` (flattened)
- `combined_figure.pdf` — composed multi-panel PDF with A/B/C panel labels
- `figure_caption.txt` — publication-ready figure caption

| Fig. | wf folder | Skill | Panels | Caption title |
|------|-----------|-------|--------|---------------|
| S1 | wf004-gcta-gwas | gcta-gwas | 2 | GCTA MLMA GWAS results (Manhattan + Q-Q) |
| S2 | wf005-genomic-selection | genomic-selection | 1 | 6-model genomic selection comparison |
| S3 | wf009-pca | pca | 5 | PCA of population structure (scree + PC1-2-3 + density) |
| S4 | wf010-ld | ld | 1 | LD decay curve (PopLDdecay, all chromosomes, 300 kb) |
| S5 | wf011-ld-score | ld-score | 2 | LD Score distribution + LD score vs MAF |
| S6 | wf012-allele-genotype-frequency | allele-genotype-frequency | 2 | Allele & genotype frequency (overview + detail) |
| S7 | wf013-homozygosity-heterozygosity | homozygosity-heterozygosity | 3 | Homozygosity/heterozygosity (chr + sample + quality) |
| S8 | wf014-inbreeding-relationship | inbreeding-relationship | 2 | Inbreeding F distribution + method comparison |
| S9 | wf015-relationship-matrix | relationship-matrix | 2 | GRM distribution + heatmap overview |
| S10 | wf016-single-trait-model | single-trait-model | 3 | SSGBLUP (EBV + variance components + pie) |
| S11 | wf017-linear-mixed-model | linear-mixed-model | 1 | GBLUP (EMAI-REML) analysis |
| S12 | wf018-snp-effect | snp-effect | 1 | RR-BLUP SNP effect estimation |
| S13 | wf019-gxe-model | gxe-model | 3 | GxE (forest + random effects + variance bar) |
| S14 | wf020-multi-trait-model | multi-trait-model | 7 | Multi-trait (ANOVA + beta forest + r_g + h² + random + r_e + VC) |
| S15 | wf021-repeated-records-model | repeated-records-model | 3 | Repeated records (EBV + fixed effects + VC) |
| S16 | wf022-gebv-gprs-prediction | gebv-gprs-prediction | 1 | GEBV/GPRS prediction |
| S17 | wf023-qc | qc | 4 | QC diagnostic (call rate + het + missing + PCA) |
| S18 | wf024-roh | roh | 2 | ROH (F_ROH barplot + length distribution) |
| S19 | wf025-nj-tree | nj-tree | 1 | NJ phylogenetic tree (IBS distance) |
| S20 | wf026-internal-impute | internal-impute | 3 | Imputation QC (dr² + info score + histogram) |

Composition script: `lib/compose_figures.py` (pypdfium2 + pikepdf + reportlab; vector PDF embedding via XObject placement — source PDFs are preserved as vectors without rasterization; A/B/C panel labels are added as a reportlab overlay merged via pikepdf, no background box).
Caption source: `lib/figure_captions.py` (20 publication-ready English captions).
Manifest: `lib/figures_manifest.json` (full inventory of source figure files).

## Limitations & Caveats

- The test data is intentionally small so that the full 26-skill sweep completes in minutes. Timing figures are not representative of production scale; production runs on full-cohort WGS or large SNP arrays will be substantially slower and may exercise code paths that the small test data does not reach.
- The benchmark does NOT modify the underlying PigGPA skill scripts. All 26 skills are tested as-is with their original implementation.
