# -*- coding: utf-8 -*-
"""返修批次 3：功效计算（M11）+ I²_GX（M12）+ MR-PRESSO 校准验证（M8）+ MVMR 工具变量覆盖率（M4）。

只读既有数据，不改任何既有结果文件；输出写 11-返修/。
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

BASE = Path(r"D:/2027/workbuddy working/021/MR演练")
OUT = BASE / "11-返修"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, r"C:/Users/柯莉尔/.workbuddy/skills/00-MR全流程总控/scripts")
import mr_engine as me  # noqa: E402

H = BASE / "04-数据" / "harmonise"
AN = BASE / "05-分析"

Z_ALPHA = {0.05: 1.959964, 0.0125: 2.497672}   # 双侧
Z_POWER = 0.841621                                # 80%


def norm_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# ---------------- A. 功效（M11）----------------
def power_block(label, beta, se):
    out = {"label": label, "beta": beta, "se": se, "or": math.exp(beta),
           "or_ci": [math.exp(beta - 1.959964 * se), math.exp(beta + 1.959964 * se)]}
    for a, z in Z_ALPHA.items():
        mde = (z + Z_POWER) * se
        out[f"mde_or_alpha{a}"] = math.exp(mde)
        out[f"power_at_observed_alpha{a}"] = norm_cdf(abs(beta) / se - z)
    return out


power = []
for tag, f, key in (("FM → IBS（发现）", "d2-FM-IBS-univariate.json", None),
                    ("IBS → FM（发现）", "d1-IBS-FM-univariate.json", None),
                    ("FM → IBS（FinnGen 复制）", "d2-finnrep-FM-IBS.json", None)):
    d = json.loads((AN / f).read_text(encoding="utf-8-sig"))
    pr = d[d["primary_analysis"]]
    power.append(power_block(tag, pr["beta"], pr["se"]))
(OUT / "power.json").write_text(json.dumps(power, ensure_ascii=False, indent=1), encoding="utf-8")

print("=== M11 功效（α = 0.0125 双侧、power = 0.8 下的最小可检出 OR）===")
for p in power:
    print(f"  {p['label']}: 观测 OR {p['or']:.3f} | SE {p['se']:.4f} | "
          f"MDE(α=.0125) {p['mde_or_alpha0.0125']:.2f} | "
          f"观测效应处 power {p['power_at_observed_alpha0.0125']:.2f}")

# ---------------- B. I²_GX（M12，NOME 违例的量化）----------------
def igx(rows):
    bx = np.array([float(r["beta_exposure"]) for r in rows])
    sx = np.array([float(r["se_exposure"]) for r in rows])
    k = len(bx)
    w = 1.0 / sx ** 2
    bw = (w * bx).sum() / w.sum()
    q = (w * (bx - bw) ** 2).sum()
    return k, q, max(0.0, (q - (k - 1)) / q)


igx_out = []
for tag, p in (("FM → IBS（发现）", H / "FM-IBS" / "harmonised.csv"),
               ("IBS → FM（发现）", H / "IBS-FM" / "harmonised.csv"),
               ("FM → IBS（FinnGen 复制）", H / "FM-IBS-finnrep" / "harmonised.csv")):
    rows = list(csv.DictReader(open(p, encoding="utf-8")))
    k, q, i2 = igx(rows)
    igx_out.append({"label": tag, "k": k, "q_gx": q, "i2_gx": i2})
    print(f"  {tag}: k={k}, Q_GX={q:.2f}, I²_GX={i2:.3f}  ({'NOME 违例严重' if i2 < 0.5 else '可接受'})")
(OUT / "igx.json").write_text(json.dumps(igx_out, ensure_ascii=False, indent=1), encoding="utf-8")

# ---------------- C. MR-PRESSO 校准（M8）----------------
# 本机无 R，无法与官方 MRPRESSO 包直接比对；改为直接检验其统计性质：
#   C1 零假设校准：无多效性数据下，global 检验的拒绝率应 ≈ 名义 α
#   C2 检出能力：人工植入 2 个离群变异，看能否被识别且不误伤
rng = np.random.default_rng(20261022)

def one_rep(k, bx, sx, bx_sd, beta_true, sy, n_outlier=0, n_boot=300, seed=None):
    r = np.random.default_rng(seed)
    by = beta_true * bx + r.normal(0, sy, k)
    if n_outlier:
        idx = r.choice(k, n_outlier, replace=False)
        by[idx] += r.normal(0, 0.25, n_outlier)   # 植入水平多效性
    g = me.presso_global(bx, by, np.full(k, sy), n_boot=n_boot, rng=np.random.default_rng(seed + 1))
    return g["p_global"]


calib = {}
for k in (15, 5):
    bx = np.linspace(0.05, 0.10, k)
    sx = np.full(k, 0.015)
    sy = 0.014
    reps = 400
    ps = [one_rep(k, bx, sx, 0.02, 0.20, sy, 0, 300, seed=1000 + i) for i in range(reps)]
    rate05 = float(np.mean(np.array(ps) < 0.05))
    rate10 = float(np.mean(np.array(ps) < 0.10))
    calib[f"k={k}"] = {"reps": reps, "reject_rate_at_0.05": rate05, "reject_rate_at_0.10": rate10}
    print(f"  PRESSO 零校准 k={k}: 拒绝率@0.05 = {rate05:.3f}（名义 0.05）, @0.10 = {rate10:.3f}")

det = {}
for k, n_out in ((15, 2), (5, 2)):
    bx = np.linspace(0.05, 0.10, k)
    sx = np.full(k, 0.015)
    sy = 0.014
    reps = 200
    hits = 0
    for i in range(reps):
        p = one_rep(k, bx, sx, 0.02, 0.20, sy, n_out, 300, seed=5000 + i)
        hits += (p < 0.05)
    det[f"k={k},outliers={n_out}"] = {"reps": reps, "detection_rate": hits / reps}
    print(f"  PRESSO 检出 k={k} 植入 {n_out} 个离群: 检出率 = {hits / reps:.3f}")

(OUT / "presso_validation.json").write_text(
    json.dumps({"null_calibration": calib, "detection": det}, ensure_ascii=False, indent=1), encoding="utf-8")

# ---------------- D. MVMR 工具变量覆盖率（M4）----------------
cov = []
for tag, mdir, main_json in (("d1（IBS 主暴露 → FM）", "MVMR-d1", "d1-IBS-FM-univariate.json"),
                             ("d2（FM 主暴露 → IBS）", "MVMR-d2", "d2-FM-IBS-univariate.json")):
    rows = list(csv.DictReader(open(H / mdir / "mvmr-matrix.csv", encoding="utf-8")))
    ids = {r["rsid"] for r in rows}
    main_iv = json.loads((AN / main_json).read_text(encoding="utf-8-sig"))["snps"]
    present = [s for s in main_iv if s in ids]
    from collections import Counter
    lead = Counter(r["lead_exposure"] for r in rows)
    cov.append({"label": tag, "n_matrix": len(rows), "lead_exposure_counts": dict(lead),
                "main_iv_total": len(main_iv), "main_iv_in_matrix": len(present),
                "main_iv_present": present})
    print(f"  {tag}: 矩阵 {len(rows)} 行；主暴露自身 IV {len(present)}/{len(main_iv)} 在矩阵内 {present}")

# 协变量库对 FM / IBS 工具变量的覆盖（解释为何丢位点）
fm_iv = json.loads((AN / "d2-FM-IBS-univariate.json").read_text(encoding="utf-8-sig"))["snps"]
ibs_iv = json.loads((AN / "d1-IBS-FM-univariate.json").read_text(encoding="utf-8-sig"))["snps"]
coverage = {}
for lib in ("ieu-b-40", "ieu-b-102", "ukb-b-3957"):
    rows = list(csv.DictReader(open(BASE / "04-数据" / "iv" / f"{lib}-associations.csv", encoding="utf-8")))
    ids = {r["rsid"] for r in rows}
    coverage[lib] = {"n_rows": len(rows),
                     "FM_iv_covered": sum(s in ids for s in fm_iv),
                     "IBS_iv_covered": sum(s in ids for s in ibs_iv)}
    print(f"  库 {lib}: 共 {len(rows)} 行；覆盖 FM IV {coverage[lib]['FM_iv_covered']}/15，"
          f"IBS IV {coverage[lib]['IBS_iv_covered']}/5")

(OUT / "mvmr_coverage.json").write_text(
    json.dumps({"matrix": cov, "library_coverage": coverage}, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n输出已写入", OUT)
