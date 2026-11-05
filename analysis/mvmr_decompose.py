# -*- coding: utf-8 -*-
"""定位差异来源：把"工具集"与"是否含 BMI"两个因素拆开。

三种模型，同一批数据、同一引擎：
  A 四暴露（FM, DEP, INS, BMI）→ IBS，全集工具（复现原分析，应约等于 0.138）
  B 三暴露（FM, DEP, INS）→ IBS，全集工具（只去掉 BMI）
  C 三暴露（FM, DEP, INS）→ IBS，仅 FM 自身 IV（只换工具集）

另做方向二（IBS → FM）的 B、C 对照。
输出写 11-返修/，不覆盖既有结果。
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

BASE = Path(r"D:/2027/workbuddy working/021/MR演练")
OUT = BASE / "11-返修"
IV = BASE / "04-数据" / "iv"
sys.path.insert(0, r"C:/Users/柯莉尔/.workbuddy/skills/00-MR全流程总控/scripts")
import mr_engine as me  # noqa: E402

PALIN = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
FILES = {"FM": "GCST90838603-local-associations.csv",
         "IBS": "ebi-a-GCST90016564-associations.csv",
         "DEP": "ieu-b-102-associations.csv",
         "INS": "ukb-b-3957-associations.csv",
         "BMI": "ieu-b-40-associations.csv"}


def load(name):
    d = {}
    for r in csv.DictReader(open(IV / FILES[name], encoding="utf-8")):
        d[r["rsid"]] = {"ea": (r.get("ea") or r.get("effect_allele", "")).upper(),
                        "nea": (r.get("nea") or r.get("other_allele") or r.get("non_effect_allele", "")).upper(),
                        "beta": float(r["beta"]), "se": float(r["se"]),
                        "eaf": float(r["eaf"]) if r.get("eaf") not in (None, "", "NA") else None}
    return d


D = {k: load(k) for k in FILES}


def align_to(ref, row):
    """把 row 的效应对齐到参考等位基因对。回文且缺 EAF 时丢弃。"""
    rea, rnea = ref["ea"], ref["nea"]
    if (rea, rnea) in PALIN:
        return None
    if row["ea"] == rea:
        return row["beta"], row["se"]
    if row["ea"] == rnea:
        return -row["beta"], row["se"]
    return None


def build(exposures, outcome, snp_pool=None):
    ref_name = exposures[0]
    ids = snp_pool if snp_pool is not None else set(D[ref_name])
    rows = []
    for s in ids:
        if any(s not in D[x] for x in list(exposures) + [outcome]):
            continue
        ref = D[ref_name][s]
        vec, ok = [], True
        for x in exposures:
            a = align_to(ref, D[x][s]) if x != ref_name else (ref["beta"], ref["se"])
            if a is None:
                ok = False
                break
            vec.append(a)
        if not ok:
            continue
        o = align_to(ref, D[outcome][s])
        if o is None:
            continue
        vals = [v[0] for v in vec] + [v[1] for v in vec] + [o[0], o[1]]
        if any((not np.isfinite(x)) or x <= 0 for x in [v[1] for v in vec] + [o[1]]):
            continue
        if any(not np.isfinite(x) for x in [v[0] for v in vec] + [o[0]]):
            continue
        rows.append((s, vec, o))
    BX = np.array([[v[0] for v in vec] for _, vec, _ in rows])
    SX = np.array([[v[1] for v in vec] for _, vec, _ in rows])
    by = np.array([o[0] for _, _, o in rows])
    sy = np.array([o[1] for _, _, o in rows])
    return rows, BX, SX, by, sy


def run(label, exposures, outcome, snp_pool=None):
    rows, BX, SX, by, sy = build(exposures, outcome, snp_pool)
    if len(rows) <= len(exposures):
        print(f"{label}: n={len(rows)} ≤ 暴露数 {len(exposures)}，不可识别")
        return {"label": label, "n": len(rows), "identified": False}
    res = me.mv_ivw(BX, SX, by, sy, model="fixed")
    i = 0
    b, s = res["beta"][i], res["se"][i]
    p = math.erfc(abs(b) / s / math.sqrt(2))
    out = {"label": label, "n": len(rows), "identified": True,
           "exposures": exposures, "outcome": outcome,
           "core_beta": b, "core_se": s, "core_ci": [b - 1.96 * s, b + 1.96 * s], "core_p": p,
           "core_f_cond": res["f_conditional"][0],
           "q": res["q"], "q_p": res["q_pval"],
           "all": {exposures[j]: {"beta": res["beta"][j], "se": res["se"][j],
                                  "f_cond": res["f_conditional"][j]} for j in range(len(exposures))}}
    print(f"{label}: n={len(rows)}, {exposures[0]} beta={b:+.3f} "
          f"(95% CI {b-1.96*s:+.3f} to {b+1.96*s:+.3f}), P={p:.3g}, 条件F={res['f_conditional'][0]:.1f}")
    return out


print("=== 方向二：FM → IBS ===")
pool = set((BASE / "04-数据" / "harmonise" / "MVMR-d2" / "mvmr-iv-list.txt").read_text(encoding="utf-8").split())
A = run("A 四暴露（含 BMI）全集工具", ["FM", "DEP", "INS", "BMI"], "IBS", pool)
B = run("B 三暴露（去 BMI）全集工具", ["FM", "DEP", "INS"], "IBS", pool)
FM_IVS = set(json.loads((BASE / "05-分析" / "d2-FM-IBS-univariate.json").read_text(encoding="utf-8-sig"))["snps"])
IBS_IVS = set(json.loads((BASE / "05-分析" / "d1-IBS-FM-univariate.json").read_text(encoding="utf-8-sig"))["snps"])
C = run("C 三暴露 仅 FM 自身 IV", ["FM", "DEP", "INS"], "IBS", FM_IVS)

print("\n=== 方向一：IBS → FM ===")
pool1 = set((BASE / "04-数据" / "harmonise" / "MVMR-d1" / "mvmr-iv-list.txt").read_text(encoding="utf-8").split())
A1 = run("A 四暴露（含 BMI）全集工具", ["IBS", "DEP", "INS", "BMI"], "FM", pool1)
B1 = run("B 三暴露（去 BMI）全集工具", ["IBS", "DEP", "INS"], "FM", pool1)
C1 = run("C 三暴露 仅 IBS 自身 IV", ["IBS", "DEP", "INS"], "FM", IBS_IVS)

(OUT / "mvmr_decomposition.json").write_text(
    json.dumps([A, B, C, A1, B1, C1], ensure_ascii=False, indent=1), encoding="utf-8")
print("\n已写入", OUT / "mvmr_decomposition.json")
