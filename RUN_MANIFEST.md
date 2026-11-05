# 运行清单：论文中每个数字 / 每张图对应的脚本

> 命令按"脚本 → 输入 → 输出"给出。带 ✅ 的为本次实际执行并核对过的命令；
> 其余为脚本用途映射，命令行参数见各脚本文件头的说明或 `--help`。

## 0. 自检

| 用途 | 命令 | 结果 |
|---|---|---|
| 引擎 20 项自检 | `python engine/mr_engine.py --selftest` | ✅ 20/20 通过 |

## 1. 主分析（对应 Table 1、§3.1–3.3）

| 输出 | 脚本 | 输入 | 命令 |
|---|---|---|---|
| `d1-IBS-FM-univariate.json` | `engine/mr_engine.py` | IBS→FM 对齐后的 harmonised CSV | ✅ `python engine/mr_engine.py --input <harmonised.csv> --output <out.json> --ivw-model auto --n-boot 1000 --seed 20260922` |
| `d2-FM-IBS-univariate.json` | 同上 | FM→IBS harmonised CSV | ✅ 同上（换输入/输出路径） |
| `d2-finnrep-FM-IBS.json` | 同上 | FinnGen 复制集 | 同上 |
| MR-PRESSO 剔除后估计 | `analysis/presso_clean_ci.py` | 上述 JSON | 以 **官方 R 包 MRPRESSO** 为准（10,000 次置换），本脚本仅整理输出 |

## 2. 多变量 MR（对应 Table 2、§3.4）

| 输出 | 脚本 | 说明 |
|---|---|---|
| `d1-mv4 / d1-mv3 / d2-mv4 / d2-mv3` | `analysis/run_mvmr.py` | 四暴露（核心+BMI+抑郁+失眠）与三暴露（核心+BMI+抑郁）两个模型 |
| 暴露自身工具变量模型 | `analysis/mvmr_exposure_own_iv.py` | §3.4 "restricted to the exposure's own variants" |
| A/B/C 三模型拆解 | `analysis/mvmr_decompose.py` | 讨论中"净效应归零"的拆解依据 |

## 3. 共线性与敏感性（对应 §3.3、Supplementary S3/S4/S6）

| 输出 | 脚本 |
|---|---|
| 区域后验概率 PP.H0–H4 | `coloc/coloc_abf.py` + `coloc/extract_and_run_coloc.py` |
| 与 R/coloc 交叉核对 | `coloc/coloc_crosscheck.R` |
| Steiger 方向性 | `sensitivity/run_steiger_fill.py` |
| 敏感性汇总 | `sensitivity/build_sensitivity_summary.py` |
| 导出 R 核对数据 | `sensitivity/export_for_r_crosscheck.py` |

## 4. 图表

| 输出 | 脚本 | 命令 |
|---|---|---|
| Figure 1（双面板 A/B） | `figures/gen_figure1_combined.R` | `Rscript gen_figure1_combined.R` |
| Figure 2–4、Table 1–2 | `figures/gen_figures_tables.R` | `Rscript gen_figures_tables.R` |
| Supplementary Figure S1–S6 | `figures/gen_supp_figures.R` | `Rscript gen_supp_figures.R` |

> 作图脚本依赖同目录的 `_common.R`（样式、路径与读数），须与本目录保持一致。
