# HIBLUP 包装脚本

piggpa-G 与 HIBLUP 基准测试中所有 HIBLUP 命令的可复现 shell 包装脚本。
每个包装脚本都严格复现了 `HIBLUP_results/` 下原始 `.log` 文件中记录的
HIBLUP 调用命令，无需查阅日志即可端到端重跑基准测试。

## 关于 HIBLUP

HIBLUP 是 C++ 二进制程序（不是 Python 脚本）。本基准测试使用的版本如下：

| 属性 | 值 |
|------|-----|
| 软件 | HIBLUP |
| 版本 | v1.6.0（2025-09-29 发布） |
| 平台 | Linux x86_64 |
| 二进制路径 | `/public/share/likui/hanyu/software/bin/hiblup` |
| 引用 | Yin 等（HIBLUP，www.hiblup.com） |

所有命令均原样提取自每个 HIBLUP `.log` 文件顶部的 `Commands:` 段落。
跨行命令以 `\` 续行符重新拼接；二进制路径和共享输入文件的绝对路径
严格保留原始记录。本地工作目录输入文件（如 `train_samples.txt`、
`chr1_snps.txt`、`phenotype_hiblup.txt`、`pedigree_hiblup.txt`）使用
相对路径引用，因为它们是在每个任务的 `HIBLUP_results/{n}/` 工作目录
中准备的。

## 包装脚本清单

| 序号 | 包装脚本 | 任务 | HIBLUP_results 目录 | 功能 | 状态 |
|------|----------|------|---------------------|------|------|
| 1 | `hiblup_allele_frequency.sh`        | T1/T2 | `1.2` | T1：等位基因频率、基因型频率；T2：杂合度、纯合度（单次 HIBLUP 运行） | 正常 |
| 2 | `hiblup_relationship_matrix.sh`     | T3    | `3`   | GA / GD / HA / PA 关系矩阵（--make-xrm） | 正常 |
| 3 | `hiblup_inbreeding_coefficient.sh`  | T5    | `4`   | 近交系数（--ibc）与关系系数（--rc） | 正常 |
| 4 | `hiblup_pca.sh`                     | T6    | `5`   | 主成分分析（前 10 个 PC） | 正常 |
| 5 | `hiblup_blup_prediction.sh`         | T7    | `6/1` | 五模型 BLUP：BLUP / GBLUP / LM / PBLUP / SSBLUP（训练 + SNP效应 + 预测） | 正常 |
| 6 | `hiblup_single_trait_model.sh`      | T8    | `7`   | 单性状 GBLUP 方差组分估计 | 正常 |
| 7 | `hiblup_repeated_records_model.sh`  | T9    | `8`   | 重复记录模型（GA + 永久环境效应） | 正常 |
| 8 | `hiblup_multi_trait_model.sh`       | T10   | `9`   | 三性状多性状 GBLUP + 基因型编码 / GRM 验证 | 正常 |
| 9 | `hiblup_gxe_model.sh`               | T12   | `10`  | 基因-环境互作模型（--rand-gxe） | 正常（gxe_he/gxe_hi 失败，已注释） |
| 10 | `hiblup_snp_effect_calculation.sh` | 3.18  | `11/1` | 单性状 GBLUP + SNP 效应反算（T16 的前置步骤） | 正常（snp_effect_test 失败，已注释） |
| 11 | `hiblup_gebv_prediction.sh`        | T16   | `12`  | 对未参与训练个体预测 GEBV（--pred） | 正常 |
| 12 | `hiblup_ld_calculation.sh`         | T18   | `14`  | 成对 LD（--ld）与 LD 评分（--ldscore） | 正常（convert.log 失败，已注释） |
| 13 | `hiblup_ld_score_regression.sh`    | T20   | `15`  | LD 评分回归（--ldreg） | 运行时失败（见说明） |

**共计 13 个包装脚本**，覆盖 13 个 HIBLUP_results 目录
（`1.2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15`）与 14 个功能模块。
其中 T1（等位基因频率）与 T2（纯合度/杂合度）均由
`hiblup_allele_frequency.sh` 覆盖，因为 HIBLUP 的 `1.2/` 目录在一次
运行中同时计算这两类结果，因此 T2 无需单独的包装脚本。

## 关于目录 11 的说明

`HIBLUP_results/11/` 在任务映射表中最初标注为"未知"。经分析确认，它是 HIBLUP 3.18 "SNP 效应
计算"任务：两步流程，先拟合单性状 GBLUP 模型获得 GEBV，再由 GEBV
反算每个 SNP 的加性效应（`snp_effect.snpeff`）。该 SNP 效应文件是
T16（`hiblup_gebv_prediction.sh`）的直接输入，因此保留为独立包装脚本。

## 运行时失败的命令

少数日志文件记录了被 HIBLUP 拒绝的命令。这些命令在对应包装脚本中
原样保留但已注释掉，并附注失败原因。这样既忠实记录了尝试过的全部
命令，又确保 `set -e` 不会在重跑时因已知失败命令而中断。

| 包装脚本 | 失败命令 | 原因 |
|----------|----------|------|
| `hiblup_gxe_model.sh`              | `gxe_he.log` / `gxe_hi.log` | `--algorithm HE` / `--algorithm HI` 不能与 `--rand-gxe` 同时使用 |
| `hiblup_snp_effect_calculation.sh` | `snp_effect_test.log`       | 调用 `--snp-effect` 时未指定 `--gebv` |
| `hiblup_ld_calculation.sh`         | `convert.log`               | `--trans-xrm` 期望 `.id` 索引文件，但 `--ld` 输出的是 `.info`（实际改用 Python 脚本转换） |
| `hiblup_ld_score_regression.sh`    | `ldreg_result.log`          | `sumstat_hiblup.txt` 第 8 列含 `chr:pos` 字符串而非浮点数（输入格式问题，非 HIBLUP 缺陷） |

## 运行方式

```bash
# 需在工作目录中包含本地输入文件
#（train_samples.txt、chr1_snps.txt、phenotype_hiblup.txt 等），
# 或直接在对应的 HIBLUP_results/{n}/ 目录下运行。

bash hiblup_relationship_matrix.sh
bash hiblup_blup_prediction.sh
# ... 以此类推
```

每个包装脚本顶部都以变量形式定义了
`HIBLUP=/public/share/likui/hanyu/software/bin/hiblup` 和共享输入文件的
绝对路径，因此当数据布局变化时只需修改一处即可重新指向。

## 与上传包中其他目录的关系

| 目录 | 内容 |
|------|------|
| `../benchmark/`         | 基准测试结果数据（矩阵、育种值、汇总），用于 piggpa-G 与 HIBLUP 直接对比 |
| `../scripts/`           | piggpa-G Python 源脚本（基准测试的另一侧） |
| `../benchmark/{task}/hiblup/` | HIBLUP 参考结果已合并到各功能目录的 `hiblup/` 子目录（全部 14 个任务） |
| `hiblup_scripts/`（本目录） | 复现原始 `.log` 文件中每一条 HIBLUP 命令的 shell 包装脚本 |
