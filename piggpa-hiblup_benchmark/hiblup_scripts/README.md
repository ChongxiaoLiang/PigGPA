# HIBLUP Wrapper Scripts

Reproducible shell wrappers for every HIBLUP command run in the piggpa-G vs
HIBLUP benchmark. Each wrapper reproduces the exact HIBLUP invocation(s)
recorded in the original `.log` files under `HIBLUP_results/`, so the
benchmark can be re-run end-to-end without consulting the logs.

## About HIBLUP

HIBLUP is a C++ binary (NOT a Python script). The version used in this
benchmark is:

| Property | Value |
|----------|-------|
| Software | HIBLUP |
| Version  | v1.6.0 (2025-09-29 Release) |
| Platform | Linux x86_64 |
| Binary   | `/public/share/likui/hanyu/software/bin/hiblup` |
| Citation | Yin et al. (HIBLUP, www.hiblup.com) |

All commands were extracted verbatim from the `Commands:` section at the top
of each HIBLUP `.log` file. Multi-line commands were re-joined with `\`
continuations; absolute paths to the binary and shared inputs were preserved
exactly as recorded. Local working-directory inputs (e.g. `train_samples.txt`,
`chr1_snps.txt`, `phenotype_hiblup.txt`, `pedigree_hiblup.txt`) are
referenced by relative path because they were prepared in each task's
`HIBLUP_results/{n}/` working directory.

## Wrapper Scripts

| # | Wrapper script | Task | HIBLUP_results dir | Function | Status |
|---|----------------|------|--------------------|----------|--------|
| 1 | `hiblup_allele_frequency.sh`        | T1/T2 | `1.2` | T1: allele freq, genotype freq; T2: heterozygosity, homozygosity (single HIBLUP run) | OK |
| 2 | `hiblup_relationship_matrix.sh`     | T3    | `3`   | GA / GD / HA / PA relationship matrices (--make-xrm) | OK |
| 3 | `hiblup_inbreeding_coefficient.sh`  | T5    | `4`   | Inbreeding (--ibc) and relationship (--rc) coefficients | OK |
| 4 | `hiblup_pca.sh`                     | T6    | `5`   | PCA (top 10 PCs) | OK |
| 5 | `hiblup_blup_prediction.sh`         | T7    | `6/1` | 5-model BLUP: BLUP / GBLUP / LM / PBLUP / SSBLUP (train + SNP-effect + pred) | OK |
| 6 | `hiblup_single_trait_model.sh`      | T8    | `7`   | Single-trait GBLUP variance component estimation | OK |
| 7 | `hiblup_repeated_records_model.sh`  | T9    | `8`   | Repeated-records model (GA + permanent environment) | OK |
| 8 | `hiblup_multi_trait_model.sh`       | T10   | `9`   | 3-trait multi-trait GBLUP + genotype-coding / GRM verification | OK |
| 9 | `hiblup_gxe_model.sh`               | T12   | `10`  | GxE interaction model (--rand-gxe) | OK (gxe_he/gxe_hi failed, commented out) |
| 10 | `hiblup_snp_effect_calculation.sh` | 3.18  | `11/1` | Single-trait GBLUP + SNP-effect back-calculation (prereq for T16) | OK (snp_effect_test failed, commented out) |
| 11 | `hiblup_gebv_prediction.sh`        | T16   | `12`  | GEBV prediction on held-out individuals (--pred) | OK |
| 12 | `hiblup_ld_calculation.sh`         | T18   | `14`  | Pairwise LD (--ld) and LD scores (--ldscore) | OK (convert.log failed, commented out) |
| 13 | `hiblup_ld_score_regression.sh`    | T20   | `15`  | LD score regression (--ldreg) | FAILED at runtime (see notes) |

**13 wrapper scripts total**, covering 13 HIBLUP_results directories
(`1.2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15`) and 14 functional
modules. T1 (allele frequency) and T2 (homozygosity/heterozygosity) are
both covered by `hiblup_allele_frequency.sh` because HIBLUP's `1.2/`
directory computes both in a single run, so no separate wrapper is
needed for T2.

## Notes on directory 11

`HIBLUP_results/11/` was originally listed as "unknown" in the task mapping.
Analysis shows it is the HIBLUP 3.18
"SNP effect calculation" task: a 2-step pipeline that first fits a
single-trait GBLUP model to obtain GEBVs, then back-calculates per-SNP
additive effects (`snp_effect.snpeff`). The resulting SNP-effect file is the
direct input for T16 (`hiblup_gebv_prediction.sh`), so it is kept as a standalone wrapper.

## Commands that failed at runtime

A few log files record commands that HIBLUP rejected. They are reproduced
verbatim inside the relevant wrapper scripts but commented out, with a note
explaining the failure. This keeps the wrappers a faithful record of what
was attempted while ensuring `set -e` does not abort a re-run on a known
broken command.

| Wrapper | Failed command | Reason |
|---------|----------------|--------|
| `hiblup_gxe_model.sh`              | `gxe_he.log` / `gxe_hi.log` | `--algorithm HE` / `--algorithm HI` not supported together with `--rand-gxe` |
| `hiblup_snp_effect_calculation.sh` | `snp_effect_test.log`       | `--snp-effect` called without `--gebv` |
| `hiblup_ld_calculation.sh`         | `convert.log`               | `--trans-xrm` expects `.id` index, but `--ld` produces `.info` (used a Python script instead) |
| `hiblup_ld_score_regression.sh`    | `ldreg_result.log`          | 8th column of `sumstat_hiblup.txt` contained a `chr:pos` string, not a float (input format issue, not a HIBLUP defect) |

## How to run

```bash
# From a working directory that contains the local inputs
# (train_samples.txt, chr1_snps.txt, phenotype_hiblup.txt, etc.)
# OR from the corresponding HIBLUP_results/{n}/ directory.

bash hiblup_relationship_matrix.sh
bash hiblup_blup_prediction.sh
# ... etc.
```

Each wrapper defines `HIBLUP=/public/share/likui/hanyu/software/bin/hiblup`
and shared absolute input paths as variables at the top, so paths can be
re-pointed by editing one place if the data layout changes.

## Relationship to other directories in the upload package

| Directory | Contents |
|-----------|----------|
| `../benchmark/`         | Benchmark result data (matrices, BVs, summaries) for direct piggpa-G vs HIBLUP comparison |
| `../scripts/`           | piggpa-G Python source scripts (the other side of the benchmark) |
| `../benchmark/{task}/hiblup/` | HIBLUP reference outputs merged into each functional directory's `hiblup/` subdirectory (all 14 tasks) |
| `hiblup_scripts/` (this dir) | Shell wrappers that reproduce every HIBLUP command from the original `.log` files |
