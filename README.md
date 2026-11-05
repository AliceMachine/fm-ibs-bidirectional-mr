# 复现代码：纤维肌痛遗传易感性与肠易激综合征风险的双向孟德尔随机化

本仓库为论文《Genetic liability to fibromyalgia and risk of irritable bowel syndrome:
a bidirectional Mendelian randomisation study》的分析代码，随文提交以供复现。

## 数据来源

全部为公开的 GWAS 汇总统计，不含任何个体数据，因此无需伦理审批。

| 角色 | 性状 | 数据集 | 样本量 |
|---|---|---|---|
| 结局/暴露 | 肠易激综合征（IBS） | IEU OpenGWAS `ebi-a-GCST90016564` | 486,601（53,400 例） |
| 暴露/结局 | 纤维肌痛（FM） | GWAS Catalog `GCST90838603` | 2,313,287（49,000 例） |
| 复制 | IBS（芬兰） | IEU OpenGWAS `finn-b-K11_IBS` | — |
| 协变量 | BMI / 重度抑郁 / 失眠 | `ieu-b-40` / `ieu-b-102` / `ukb-b-3957` | — |

下载后请按 `RUN_MANIFEST.md` 的输入表放置，或直接用 `--input` 指定路径。

## 环境

- Python 3.13.12，依赖见 `requirements.txt`（numpy 2.5.4、scipy 1.18.1）
- R 4.5.1，作图需 `dplyr`、`ggplot2`、`gt`、`jsonlite`、`patchwork`、`scales`
- **随机种子 20260922**（Bootstrap 与置换检验一律使用该种子）

## 快速自检

```bash
pip install -r requirements.txt
python engine/mr_engine.py --selftest      # 预期 20/20 通过
```

自检覆盖 IVW（固定/随机）、MR-Egger、加权中位数/众数、Wald、Cochran's Q 与 I²、
留一法、MR-PRESSO、Steiger、多变量 IVW，以及"缺列必须报错"等边界情形。

## 目录

```
engine/        MR 引擎（自包含，仅依赖 numpy/scipy）
coloc/         共线性（coloc.abf）实现与 R 交叉核对
analysis/      主分析、MR-PRESSO 清理、MVMR 及其拆解
sensitivity/   敏感性汇总、Steiger、导出 R 核对数据
figures/       正文与补充材料的作图脚本（R）
```

## 说明

- 引擎是自建实现，正确性由 `--selftest` 的 20 项断言保证，并与 R 的 TwoSampleMR /
  coloc 做过交叉核对；**MR-PRESSO 的结果一律以官方 R 包为准**，自建实现只用于早期探索。
- 脚本默认只读数据，不写原始文件。
- 中文路径在 Windows 上需先设 `Sys.setlocale("LC_CTYPE","English_United States.65001")`，
  作图脚本开头已处理。

## 许可

代码以 MIT 许可发布；数据版权归各原始数据集所有，使用时请引用其对应文献。
