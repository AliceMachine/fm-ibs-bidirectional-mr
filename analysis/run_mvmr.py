# -*- coding: utf-8 -*-
"""MVMR 运行器：把 04 阶段的对齐矩阵喂给 mr_engine.mv_ivw，产出 4 条多变量分析的 JSON。

模型清单（与 OSF 预注册 v1.1 D4 修订一致）：
  d1-mv4 : IBS + BMI + DEP + INS -> FM   （方向一主 MVMR）
  d1-mv3 : IBS + BMI + DEP       -> FM   （D4 敏感性对照 (i)：去失眠）
  d2-mv4 : FM  + BMI + DEP + INS -> IBS  （方向二主 MVMR）
  d2-mv3 : FM  + BMI + DEP       -> IBS  （D4 敏感性对照 (i)：去失眠）

固定/随机效应规则与单变量相同（Q P<0.05 -> 乘性随机），两模型都输出并标 primary。
"""
import csv, json, os, sys, pathlib
import numpy as np

# ---- 路径：不依赖任何本机目录（见 README「运行环境」） ----------------------
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA_ROOT = pathlib.Path(os.environ.get("MR_PROJECT_ROOT", REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "engine"))
import mr_engine as eng

BASE = DATA_ROOT / "04-数据"   # 04-数据
OUT = DATA_ROOT / "05-分析"
OUT.mkdir(exist_ok=True)

MODELS = [
    # (tag, matrix_dir, exposures, outcome)
    ("d1-mv4", "MVMR-d1", ["IBS", "BMI", "DEP", "INS"], "FM"),
    ("d1-mv3", "MVMR-d1", ["IBS", "BMI", "DEP"],        "FM"),
    ("d2-mv4", "MVMR-d2", ["FM", "BMI", "DEP", "INS"],  "IBS"),
    ("d2-mv3", "MVMR-d2", ["FM", "BMI", "DEP"],         "IBS"),
]

summary = []
for tag, mdir, exps, outc in MODELS:
    rows = list(csv.DictReader(open(BASE / "harmonise" / mdir / "mvmr-matrix.csv", encoding="utf-8")))
    # 完整案例过滤：本模型用到的 暴露+结局 的 beta/se 必须齐全且 se>0
    used = exps + [outc]
    def ok(r):
        for e in used:
            try:
                if not r[f"se_{e}"] or float(r[f"se_{e}"]) <= 0:
                    return False
                if r[f"beta_{e}"] == "" :
                    return False
            except (KeyError, ValueError):
                return False
        return True
    dropped = [r["rsid"] for r in rows if not ok(r)]
    rows = [r for r in rows if ok(r)]
    rsids = [r["rsid"] for r in rows]
    BX = np.array([[float(r[f"beta_{e}"]) for e in exps] for r in rows])
    SX = np.array([[float(r[f"se_{e}"]) for e in exps] for r in rows])
    by = np.array([float(r[f"beta_{outc}"]) for r in rows])
    sy = np.array([float(r[f"se_{outc}"]) for r in rows])
    n, k = BX.shape

    fixed = eng.mv_ivw(BX, SX, by, sy, model="fixed")
    random = eng.mv_ivw(BX, SX, by, sy, model="random")
    primary = "mv_random" if fixed["q_pval"] < 0.05 else "mv_fixed"

    res = {
        "tag": tag, "matrix": mdir, "exposures": exps, "outcome": outc,
        "n_snps": n, "n_exposures": k,
        "dropped_incomplete": dropped, "n_dropped": len(dropped),
        "estimates": {
            e: {"beta": fixed["beta"][i], "se": fixed["se"][i],
                "ci_lower": fixed["beta"][i] - 1.959963984540054 * fixed["se"][i],
                "ci_upper": fixed["beta"][i] + 1.959963984540054 * fixed["se"][i]}
            for i, e in enumerate(exps)
        },
        "estimates_random": {
            e: {"beta": random["beta"][i], "se": random["se"][i],
                "ci_lower": random["beta"][i] - 1.959963984540054 * random["se"][i],
                "ci_upper": random["beta"][i] + 1.959963984540054 * random["se"][i]}
            for i, e in enumerate(exps)
        },
        "q_fixed": fixed["q"], "q_pval": fixed["q_pval"], "q_dof": fixed["q_dof"],
        "primary": primary,
        "primary_rule": "pre-registered: multiplicative random effects when Cochran Q p<0.05, else fixed",
        "f_conditional": {e: fixed["f_conditional"][str(i)] if str(i) in fixed["f_conditional"]
                          else fixed["f_conditional"][i] for i, e in enumerate(exps)},
        "snps": rsids,
    }
    out = OUT / f"{tag}-mvmr.json"
    out.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    pv = res["estimates_random" if primary == "mv_random" else "estimates"]
    print(f"{tag}: n={n} k={k} primary={primary} Qp={fixed['q_pval']:.3g}")
    for e in exps:
        v = pv[e]
        print(f"   {e:4s}: b={v['beta']:.5g} [{v['ci_lower']:.4g},{v['ci_upper']:.4g}] F_cond={res['f_conditional'][e]:.4g}")
    summary.append(tag)

print("DONE:", summary)
