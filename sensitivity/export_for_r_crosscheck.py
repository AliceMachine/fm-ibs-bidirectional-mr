# -*- coding: utf-8 -*-
"""
导出 3 个代表区域的对齐数据供 R coloc 官方包交叉验证。

区域选择：rs2587363（d2 共定位+PRESSO离群+Steiger反向，三重印证位点）、
rs10156602（d1 PRESSO 离群、PP.H4=0.711 边缘）、rs426444（d1 主 IV）。

输出（coloc-r-crosscheck/）：
  <rs>_FM.csv / <rs>_IBS.csv   —— 对齐后行（rsid, beta, se），与 Python 引擎入参一致
  python_reference.json        —— Python 引擎的 PP（对比基准）
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import extract_and_run_coloc as ec  # noqa: E402  复用同一提取/对齐实现

OUT = HERE / "coloc-r-crosscheck"
OUT.mkdir(exist_ok=True)

TARGETS = ["rs2587363", "rs10156602", "rs426444"]
D1_IVS = ["rs426444", "rs10156602", "rs12549729", "rs9513519", "rs7106434"]
WINDOW = 500_000


def main():
    d2_ivs = ec.load_d2_ivs()
    d2_rsids = {r[0] for r in d2_ivs}

    # 第一遍：IBS 文件定位 d1 IV 位置（与主流程一致）
    ec.set_windows([])
    hits, _, cm = ec.scan_ibs(find_rsids=set(D1_IVS) | d2_rsids)
    iv_pos = {rs: ("d2", c, p) for rs, c, p in d2_ivs}
    for rs in D1_IVS:
        if rs in hits:
            iv_pos[rs] = ("d1", hits[rs][0], hits[rs][1])

    missing = [rs for rs in TARGETS if rs not in iv_pos]
    if missing:
        raise RuntimeError(f"目标区域定位失败: {missing}")

    wins = []
    for rs in TARGETS:
        _, c, p = iv_pos[rs]
        wins.append((str(c), max(1, p - WINDOW), p + WINDOW))

    # 两文件各扫一遍（只收 3 个窗口；scan_ibs 依赖全局窗口索引，须先 set_windows）
    ec.set_windows(wins)
    _, ibs_rows, cm = ec.scan_ibs(windows=wins, colmap=cm)
    fm_rows, fm_cm = ec.scan_fm(wins)
    print(f"IBS 区域行 {len(ibs_rows)}，FM 区域行 {len(fm_rows)}")

    ref = {}
    for rs in TARGETS:
        tag, c, p = iv_pos[rs]
        lo, hi = p - WINDOW, p + WINDOW
        fm_sub = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [], "label": "FM"}
        ibs_sub = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [], "label": "IBS"}
        for r in fm_rows:
            if r[fm_cm["chr"]] == str(c) and lo <= int(r[fm_cm["pos"]]) <= hi:
                fm_sub["rsid"].append(r[fm_cm["rsid"]]); fm_sub["ea"].append(r[fm_cm["ea"]])
                fm_sub["oa"].append(r[fm_cm["oa"]]); fm_sub["beta"].append(r[fm_cm["beta"]])
                fm_sub["se"].append(r[fm_cm["se"]]); fm_sub["eaf"].append(r[fm_cm["eaf"]])
        for r in ibs_rows:
            if r[cm["chr"]] == str(c) and lo <= int(r[cm["pos"]]) <= hi:
                ibs_sub["rsid"].append(r[cm["rsid"]]); ibs_sub["ea"].append(r[cm["ea"]])
                ibs_sub["oa"].append(r[cm["oa"]]); ibs_sub["beta"].append(r[cm["beta"]])
                ibs_sub["se"].append(r[cm["se"]])
                ibs_sub["eaf"].append(r[cm["eaf"]] if cm["eaf"] is not None else None)

        a_fm, a_ibs, dropped = ec.align_datasets(fm_sub, ibs_sub)
        n = len(a_fm["rsid"])
        if n < 50:
            raise RuntimeError(f"{rs} 对齐位点不足: {n}")

        # 导出对齐行（R 用 varbeta = se^2）
        for label, a in (("FM", a_fm), ("IBS", a_ibs)):
            with open(OUT / f"{rs}_{label}.csv", "w", encoding="utf-8", newline="") as f:
                f.write("rsid,beta,se\n")
                for i in range(n):
                    f.write(f"{a['rsid'][i]},{a['beta'][i]},{a['se'][i]}\n")

        r = ec.coloc_abf(a_fm, a_ibs)
        ref[rs] = {"tag": tag, "n_snp": r["n_snp"], "PP": r["PP"],
                   "lead_trait1": r["lead_trait1"], "lead_trait2": r["lead_trait2"],
                   "n_dropped_alleles": len(dropped)}
        print(rs, tag, "n =", n, "Python PP.H4 =", r["PP"]["H4"])

    (OUT / "python_reference.json").write_text(
        json.dumps(ref, ensure_ascii=False, indent=2), encoding="utf-8")
    print("saved:", OUT / "python_reference.json")


if __name__ == "__main__":
    main()
