# piggpa-G 源代码脚本

本目录包含 piggpa-G vs HIBLUP benchmark 中使用的 piggpa-G 源代码脚本。所有脚本用于生成 `../benchmark/` 中的 benchmark 结果。

完整的 benchmark 覆盖 **14 个 benchmark 目录**（14 个对比模块；见 `../benchmark/`），但本目录中的脚本仅对应其中 **4 个直接对比模块**，即 piggpa-G 与 HIBLUP 输出逐项对照的模块（T3 关系矩阵、T7 BLUP 预测、T8 单性状模型、T9 重复记录模型）。其余 10 个仅 HIBLUP 侧执行的功能的包装脚本位于 `../hiblup_scripts/`，其 piggpa-G 对应源脚本不在本上传包内（见下方"其他 piggpa-G 脚本"）。

## 目录结构

```
scripts/
├── README.md                              # 本文档（英文）
├── README-zh.md                           # 本文档（中文）
├── relationship_matrix_construction.py    # 关系矩阵构建（PRM/GRM/HRM，默认 HRM=0.95*G_adj+0.05*A）
├── single_trait_model.py                  # 单性状方差组分估计（AI-REML/EM-REML/HE）
├── repeated_records_model.py              # 重复记录模型（永久环境 + 加性遗传）
├── model_comparison.py                    # 5 模型 BLUP 对比驱动（LM/BLUP/PBLUP/GBLUP/SSBLUP）
├── LM_model.py                            # 线性模型（LM）——不使用关系矩阵
├── GBLUP_model.py                         # 基因组 BLUP——使用 GA
├── PBLUP_model.py                         # 系谱 BLUP——使用 PA
└── SSBLUP_model.py                        # 单步 BLUP——使用 HA（0.95G_adj + 0.05A）
```

**说明**：`BLUP_model.py` 不存在，因为在 piggpa-G 实现中，5 模型对比里的"BLUP"复用 GBLUP 代码路径并使用 GA（基因组关系矩阵）。此事实体现在 `blup_prediction/model_comparison.csv` 的 `relationship_matrix` 列中。

---

## 脚本与任务映射

| 脚本 | Benchmark 任务 | 输出目录（功能模块） |
|------|---------------|---------|
| `relationship_matrix_construction.py` | T3 | `../benchmark/relationship_matrix/` |
| `model_comparison.py` + 4 个模型脚本（`LM_model.py`、`PBLUP_model.py`、`GBLUP_model.py`、`SSBLUP_model.py`） | T7 | `../benchmark/blup_prediction/` |
| `single_trait_model.py` | T8 | `../benchmark/single_trait_model/` |
| `repeated_records_model.py` | T9 | `../benchmark/repeated_records_model/` |

T7 的四个子模型脚本（`LM_model.py`、`PBLUP_model.py`、`GBLUP_model.py`、`SSBLUP_model.py`）由 `model_comparison.py` 调用，输出均写入 `../benchmark/blup_prediction/`。

---

## 关键算法

### relationship_matrix_construction.py
- **PRM**（系谱关系矩阵）：Henderson 递归算法
- **GRM**（基因组关系矩阵）：VanRaden 方法（$G = ZZ'/\sum 2pq$）与 Yang 方法
- **HRM**（混合关系矩阵）：默认公式 $HRM = 0.95 \cdot G_{adj} + 0.05 \cdot A$（等价于 HIBLUP 的 HA），其中 $G_{adj} = 0.999 \cdot G + 0.001 \cdot I$
- HIBLUP 兼容公式为默认 HRM 输出。`HRM.csv` 直接使用 HIBLUP 兼容公式。

### single_trait_model.py
- **AI-REML**：平均信息 REML（默认，max_iter=20）
- **EM-REML**：期望最大化 REML（`--em-max-iter`，默认 200）
- **HE 回归**：Henderson 回归 + ridge 正则化
- 内部构造 HA 时使用与 HIBLUP 相同的 0.95/0.05 公式

### repeated_records_model.py
- 模型：`weight = 1 + sex(F) + season(F) + ID(R[E]) + GA(R[G]) + e`
- AI-REML，含永久环境效应（ID）和加性遗传效应（GA）
- 9 次迭代收敛（对比 HIBLUP 20 次未收敛）

### model_comparison.py + 模型脚本
- 驱动 5 个 BLUP 模型：LM、BLUP（=GBLUP+GA）、PBLUP、GBLUP、SSBLUP
- 输出 `model_comparison.csv`，含 `relationship_matrix` 列以明确标注可比性
- 训练集：ID 1001-2000；验证集：ID 5001-6000

---

## 依赖

- **Python**：3.13+
- **NumPy**、**pandas**、**SciPy**：数值计算
- **matplotlib**：图表生成
- **HIBLUP**（v1.6.0）：参考工具，运行 piggpa-G 脚本时非必需

## 用法

每个脚本均为独立可运行。完整调用命令可参考 `../benchmark/{relationship_matrix,blup_prediction,single_trait_model,repeated_records_model}/` 中的运行日志。

## 其他 piggpa-G 脚本

benchmark 覆盖 14 个 benchmark 目录（14 个对比模块），但仅 4 个直接对比模块（T3、T7、T8、T9）在本 `scripts/` 目录中提供 piggpa-G Python 源代码。其余 10 个功能（allele_frequency、homozygosity_heterozygosity（纯合度/杂合度）、inbreeding、PCA、multi_trait、GxE、GEBV、LD、LD_score_regression、snp_effect）仅在 HIBLUP 侧执行；其 piggpa-G 对应脚本可从以下路径获取：

```
/public/share/likui/hanyu/灵活版纯脚本/function1-function14/
```

全部 14 个功能（含 10 个仅 HIBLUP 功能）的 HIBLUP 包装 shell 脚本位于 `../hiblup_scripts/`。

## 数据溯源

这些脚本生成了 `../benchmark/` 中四个直接对比任务（T3、T7、T8、T9）的全部 piggpa-G 结果。
