# -*- coding: utf-8 -*-
"""
06 阶段补做：方向二（FM->IBS）与芬兰复制（FM->IBS-FinnGen）的 Steiger 方向性检验。

背景：05 阶段引擎调用时暴露侧（FM）样本量列（n_exposure）为空，
引擎因"无样本量无法估计解释方差比例"跳过了 Steiger。本脚本在 06 阶段补做：
  - FM 样本量 = 2,313,287（台账 #3 实测，md5 校验通过的全量文件）
  - IBS 主结局样本量 = 486,601（Eijsbouts 2021，台账 #1）
  - IBS 芬兰复制结局样本量 = 187,028（finn-b-K11_IBS，台账 #2）

注意（诚实声明）：
  - FM 病例数未从原文按欧洲祖先亚组拆分回填，这里用总样本量做 Steiger，
    属于二分类性状下的近似做法（TwoSampleMR 默认亦用总 N），报告中声明局限。
  - Steiger 结果是"解释方差方向"的支持性证据，不是因果性的证明（技能红线）。

预注册 §6.2 要求："Steiger filtering 做（按预设方向过滤），同时报告未过滤结果"。
本脚本对每条分析同时输出：
  - per_snp：逐 SNP 的方向判断
  - filtered IVW：剔除方向不正确 SNP 后的 IVW（若全部正确则与未过滤一致）
  - unfiltered IVW：即 05 主分析结果，不在本脚本重算，只引用
"""
import csv
import json
import os
import sys
from pathlib import Path

# ---- 路径：不依赖任何本机目录（见 README「运行环境」） ----------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "engine"))
import mr_engine as me  # noqa: E402

DATA_ROOT = Path(os.environ.get("MR_PROJECT_ROOT", REPO_ROOT))
BASE = DATA_ROOT / "06-敏感性"
BASE.mkdir(parents=True, exist_ok=True)

HARMONISE = Path(r"D:\2027\workbuddy working\021\MR演练\04-数据\harmonise")

N_FM = 2_313_287  # FM 暴露侧样本量（实测全量文件）
N_IBS = 486_601   # IBS 结局侧（Eijsbouts 2021）
N_IBS_FINN = 187_028  # IBS 结局侧（FinnGen）

TASKS = [
    {
        "tag": "d2-FM-IBS",
        "csv": HARMONISE / "FM-IBS" / "harmonised.csv",
        "n_exposure": N_FM,
        "n_outcome": N_IBS,
    },
    {
        "tag": "d2-finnrep-FM-IBS",
        "csv": HARMONISE / "FM-IBS-finnrep" / "harmonised.csv",
        "n_exposure": N_FM,
        "n_outcome": N_IBS_FINN,
    },
]


def ivw_wls(bx, by, sy, model="fixed"):
    """与引擎一致的 IVW（过原点加权最小二乘），用于 Steiger filtering 后重算。"""
    return me.ivw(bx, by, sy, model=model)


def main():
    results = {}
    for task in TASKS:
        rows = list(csv.DictReader(open(task["csv"], encoding="utf-8")))
        bx = [float(r["beta_exposure"]) for r in rows]
        sx = [float(r["se_exposure"]) for r in rows]
        by = [float(r["beta_outcome"]) for r in rows]
        sy = [float(r["se_outcome"]) for r in rows]
        snps = [r["SNP"] for r in rows]

        st = me.steiger_direction(bx, sx, by, sy, task["n_exposure"], task["n_outcome"])
        # 给 per_snp 补上 rsID，方便报告直读
        for i, item in enumerate(st.get("per_snp", [])):
            idx = item.get("SNP_index", i)
            item["rsid"] = snps[idx]

        wrong_idx = [p["SNP_index"] for p in st.get("per_snp", []) if not p["direction_correct"]]

        filtered = None
        if wrong_idx:
            keep = [i for i in range(len(bx)) if i not in wrong_idx]
            # 直接复用引擎 IVW（文献标准式 Σ(bx·by/sy²)/Σ(bx²/sy²)），
            # 保证与 05 主分析口径一致；手搓公式曾漏 bx² 因子导致估计坍缩到 0，已废弃。
            import numpy as _np
            r = me.ivw(_np.array(bx)[keep], _np.array(by)[keep], _np.array(sy)[keep], model="fixed")
            filtered = {
                "n_removed": len(wrong_idx),
                "removed_rsids": [snps[i] for i in wrong_idx],
                "ivw_beta_filtered": r["beta"],
                "ivw_se_filtered": r["se"],
                "ivw_p_filtered": r["p"],
                "model": "fixed（与 05 主分析同引擎）",
            }

        results[task["tag"]] = {
            "n_exposure": task["n_exposure"],
            "n_outcome": task["n_outcome"],
            "n_snp": len(bx),
            "steiger": st,
            "steiger_filtering": {
                "applied": bool(wrong_idx),
                "filtered_ivw": filtered,
                "unfiltered_note": "未过滤结果 = 05 主分析 IVW（全部 IV 保留）",
            },
        }

    out = BASE / "steiger-d2-completed.json"
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    for tag, r in results.items():
        st = r["steiger"]
        print(tag, "->", st.get("prop_correct"), f"({st.get('n_correct')}/{st.get('n_evaluated')})",
              "wrong:", [p["rsid"] for p in st["per_snp"] if not p["direction_correct"]])
    print("saved:", out)


if __name__ == "__main__":
    main()
