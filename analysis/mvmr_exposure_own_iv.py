# -*- coding: utf-8 -*-
"""M4 补救分析：剔除覆盖最差的 BMI 后，用"暴露自身工具变量"重估 MVMR 净效应。

背景：四暴露 MVMR 矩阵里，主暴露自身 IV 只剩 1/15（FM 方向）与 2/5（IBS 方向），
原因是 BMI 库（ieu-b-40，Yengo 2018，HapMap2 位点集）只覆盖 2/15 FM IV 与 3/5 IBS IV，
四库 complete-case 过滤把暴露自家工具变量几乎全数剔除。

本脚本做两件事：
  1. 数值证明：把 BMI 纳入模型时，主暴露自身 IV 剩几个（≥ 暴露数才可识别）；
  2. 若无 BMI，用三暴露（暴露 + 抑郁 + 失眠），仅以暴露自身 IV 为工具，重估净效应。
     结果只作诊断输出，不覆盖任何既有结果文件。
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

# ---- 路径：不依赖任何本机目录（见 README「运行环境」） ----------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
BASE = Path(os.environ.get("MR_PROJECT_ROOT", REPO_ROOT))
OUT = BASE / "11-返修"
H = BASE / "04-数据" / "harmonise"
IV = BASE / "04-数据" / "iv"
sys.path.insert(0, str(REPO_ROOT / "engine"))
import mr_engine as me  # noqa: E402

PALIN = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def load_cov(lib):
    d = {}
    for r in csv.DictReader(open(IV / f"{lib}-associations.csv", encoding="utf-8")):
        d[r["rsid"]] = r
    return d


def align(row, ea_ref, nea_ref, eaf_ref):
    """把协变量的效应对齐到参考等位基因。返回 (beta, se) 或 None。"""
    ea, nea = row["ea"].upper(), row["nea"].upper()
    b, s = float(row["beta"]), float(row["se"])
    if (ea_ref.upper(), nea_ref.upper()) in PALIN:
        return None                              # 回文，方向不可判 → 丢
    if ea == ea_ref.upper():
        return b, s
    if ea == nea_ref.upper():
        return -b, s
    return None                                   # 等位基因对不上 → 丢


DEP = load_cov("ieu-b-102")
INS = load_cov("ukb-b-3957")
BMI = load_cov("ieu-b-40")

report = {}

for tag, hpath, covs, exp_name in (
        ("d2_prime_FM_DEP_INS", H / "FM-IBS" / "harmonised.csv", [("DEP", DEP), ("INS", INS)], "FM"),
        ("d1_prime_IBS_DEP_INS", H / "IBS-FM" / "harmonised.csv", [("DEP", DEP), ("INS", INS)], "IBS")):
    rows = list(csv.DictReader(open(hpath, encoding="utf-8")))
    keep = []
    for r in rows:
        ea, nea = r["effect_allele"], r["other_allele"]
        vals = {}
        ok = True
        for name, lib in covs:
            if r["SNP"] not in lib:
                ok = False
                break
            a = align(lib[r["SNP"]], ea, nea, r["eaf_exposure"])
            if a is None:
                ok = False
                break
            vals[name] = a
        if ok:
            keep.append((r, vals))
    print(f"\n{tag}: 输入 {len(rows)} 个 {exp_name} IV，对齐后保留 {len(keep)} 个")

    # 同时报告：若把 BMI 也纳入，还剩几个
    with_bmi = 0
    for r, _ in keep:
        if r["SNP"] in BMI and align(BMI[r["SNP"]], r["effect_allele"], r["other_allele"], r["eaf_exposure"]):
            with_bmi += 1
    print(f"   若再纳入 BMI：仅剩 {with_bmi} 个（四暴露需要 ≥4 个才可识别）")

    if len(keep) < 4:
        report[tag] = {"n": len(keep), "note": "工具变量数不足 4，三暴露模型亦不可识别", "with_bmi": with_bmi}
        continue

    BX = np.array([[float(r["beta_exposure"])] + [vals[n][0] for n, _ in covs] for r, vals in keep])
    SX = np.array([[float(r["se_exposure"])] + [vals[n][1] for n, _ in covs] for r, vals in keep])
    by = np.array([float(r["beta_outcome"]) for r, _ in keep])
    sy = np.array([float(r["se_outcome"]) for r, _ in keep])

    res = me.mv_ivw(BX, SX, by, sy, model="fixed")
    names = [exp_name] + [n for n, _ in covs]
    print("   净效应（三暴露模型，仅用暴露自身 IV）：")
    est = {}
    for i, nm in enumerate(names):
        b, s = res["beta"][i], res["se"][i]
        z = abs(b) / s
        pv = math.erfc(z / math.sqrt(2))
        fc = res["f_conditional"][i]
        print(f"     {nm:4} beta = {b:+.3f} (95% CI {b-1.96*s:+.3f} to {b+1.96*s:+.3f}), "
              f"P = {pv:.3g}, conditional F = {fc:.1f}")
        est[nm] = {"beta": b, "se": s, "ci": [b - 1.96 * s, b + 1.96 * s],
                   "p": pv, "f_conditional": fc}
    report[tag] = {"n": len(keep), "with_bmi": with_bmi,
                   "estimates": est, "q": res["q"], "q_dof": res["q_dof"], "q_p": res["q_pval"]}

(OUT / "mvmr_exposure_own_iv.json").write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n已写入", OUT / "mvmr_exposure_own_iv.json")
