# 冻结说明

冻结时间：2026-11-05 21:07；2026-11-06 改路径解析；2026-11-08 统一换行为 LF 后重算 md5

对应稿件：`IBS-FM-bidirectional-MR_CGH-manuscript_v3_20261105.docx`（正文 CGH 口径 3,868 词）

每个文件均为**复制**，原件未改动。下表 md5 是在**本仓库 clone 下来后直接 `md5sum` 就能对上**的，
因为全仓库统一为 LF 换行（见 `.gitattributes`）；Windows 上 clone 也不会被换成 CRLF。
「原位置」一律为工程内相对路径，不含盘符与任何本地目录名。

> **改动记录**
>
> - **2026-11-06**：7 个文件里写死的本机绝对路径改为「仓库内 `engine/` + 环境变量 `MR_PROJECT_ROOT`」；
>   只动路径解析，算法、参数、随机种子（20260922）、数值逻辑一行未改；`--selftest` 仍 **20/20**。
> - **2026-11-08**：全部文本文件统一为 **LF**（此前 8 个文件在工作区是 CRLF，导致 clone 后
>   3 个文件的 md5 与本表不符）。**仅换行符，内容零改动**，已用 git 与逐文件字节比对确认。
> - 改动前的原始 md5 见**首次提交**（`git log --reverse` 第一条）。

| 分组 | 文件名 | 原位置（工程内相对路径） | 字节 | md5 |
|---|---|---|---|---|
| engine | mr_engine.py | scripts/mr_engine.py（技能库内；副本已在本仓库） | 30174 | `54f4e99721d48b1ba22cfc17300cea87` |
| coloc | coloc_abf.py | 06-敏感性/coloc_abf.py | 11307 | `150ea89aa01f285bc1b2164dbf03959f` |
| coloc | extract_and_run_coloc.py | 06-敏感性/extract_and_run_coloc.py | 10882 | `0c6ddf97faf8ec369e25ca6ec17eb948` |
| coloc | coloc_crosscheck.R | 06-敏感性/coloc-r-crosscheck/coloc_crosscheck.R | 2241 | `52f9bd4644ba0be146b9f359b9e7b922` |
| analysis | gen_results.py | 05-分析/gen_results.py | 11540 | `f899a95f2f3ad24557a21fb9a61b15e8` |
| analysis | presso_clean_ci.py | 05-分析/presso_clean_ci.py | 4843 | `a6fc86f9c8ce2cde443571a2466c87ad` |
| analysis | run_mvmr.py | 05-分析/run_mvmr.py | 4281 | `91f915c80da56b536213de384ee40f66` |
| analysis | batch3_analysis.py | 11-返修/batch3_analysis.py | 7336 | `7986ec69bee9c433b261e98492d8dae3` |
| analysis | mvmr_exposure_own_iv.py | 11-返修/mvmr_exposure_own_iv.py | 4742 | `c158cb7fe4ab9fa01f13579d82f89e07` |
| analysis | mvmr_decompose.py | 11-返修/mvmr_decompose.py | 5721 | `bc774988693eae7d599afa115d3a3872` |
| sensitivity | build_sensitivity_summary.py | 06-敏感性/build_sensitivity_summary.py | 9043 | `9e16ad5da54021f9e8c8a2adb01df0ab` |
| sensitivity | export_for_r_crosscheck.py | 06-敏感性/export_for_r_crosscheck.py | 4032 | `3ee1cffa585a6f975508f084b2f99750` |
| sensitivity | run_steiger_fill.py | 06-敏感性/run_steiger_fill.py | 4939 | `21cd9efd030e92d918eb7037f56d4de1` |
| figures | _common.R | 07-图表/R/_common.R | 7968 | `727f24b126b14c276a6823005bb9f831` |
| figures | gen_figures_tables.R | 07-图表/R/gen_figures_tables.R | 24639 | `273e6c4d7ea5cb111aa131d8782cd179` |
| figures | gen_supp_figures.R | 07-图表/R/gen_supp_figures.R | 14813 | `41d396e553e8978d0277a4024fb3ac88` |
| figures | gen_figure1_combined.R | 07-图表/R/gen_figure1_combined.R | 4091 | `291df234f1d52ee6d47d40db661d2bf1` |
