# -*- coding: utf-8 -*-
"""B1：PRESSO 剔离群后 IVW 固定效应估计补 95%CI 与 P（决策表 A5/B1）。
方法：从 04-数据/iv/ 重构 harmonised 数据（ alleles+eaf 对齐），
先复算全集 IVW 与 05-分析 JSON 比对验证重构忠实（容差 1e-6），
再对剔离群后子集跑同一引擎 ivw(model='fixed')，输出 beta/se/CI/P。
"""
import csv, json, math, os, pathlib, sys

import numpy as np

# ---- 路径：不依赖任何本机目录（见 README「运行环境」） ----------------------
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PROJ = pathlib.Path(os.environ.get("MR_PROJECT_ROOT", REPO_ROOT))
BASE = PROJ / "05-分析"          # 原分析结果 JSON 所在目录
IV = PROJ / "04-数据" / "iv"
sys.path.insert(0, str(REPO_ROOT / "engine"))
import mr_engine as eng  # noqa: E402

def load_assoc(path, beta_col="beta", se_col="se", rsid_col="rsid",
               ea_col="ea", nea_col=None, eaf_col="eaf"):
    rows = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows[r[rsid_col]] = {
                "beta": float(r[beta_col]), "se": float(r[se_col]),
                "ea": r[ea_col].strip().upper(),
                "oa": (r[nea_col].strip().upper() if nea_col else None),
                "eaf": (float(r[eaf_col]) if r.get(eaf_col) not in (None, "", "NA") else None),
            }
    return rows

def build_arrays(exp_rows, out_rows, snps):
    bx, sx, by, sy, used = [], [], [], [], []
    for rs in snps:
        e = exp_rows.get(rs); o = out_rows.get(rs)
        if not e or not o:
            print(f"  !! {rs} 缺失，跳过"); continue
        bxo, byo = e["beta"], o["beta"]
        # 等位基因对齐：先按效应等位基因直接匹配，再按 eaf 校验/纠偏
        if o["ea"] == e["ea"]:
            pass
        elif o["oa"] and o["oa"] == e["ea"]:
            byo = -byo
        else:
            # 无法按等位基因直接对齐（eaf 只有一侧等位基因频率），用 eaf 纠偏
            if e["eaf"] is not None and o["eaf"] is not None and abs(o["eaf"] - e["eaf"]) > abs((1 - o["eaf"]) - e["eaf"]):
                byo = -byo
        bx.append(bxo); sx.append(e["se"]); by.append(byo); sy.append(o["se"]); used.append(rs)
    return bx, sx, by, sy, used

def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(b))

for tag, exp_file, exp_nea, out_file, out_nea, json_file in [
    ("d1-IBS-FM", "iv-ebi-a-GCST90016564.csv", "nea", "GCST90838603-local-associations.csv", "oa", "d1-IBS-FM-univariate.json"),
    ("d2-FM-IBS", "GCST90838603-local-associations.csv", "oa", "ebi-a-GCST90016564-associations.csv", "nea", "d2-FM-IBS-univariate.json"),
]:
    print("==", tag)
    ref = json.load(open(BASE / json_file, encoding="utf-8"))
    snps = ref["snps"]
    exp_rows = load_assoc(IV / exp_file, nea_col=exp_nea)
    # outcome 文件列名不同：d1 结局=FM（oa 列）；d2 结局=IBS（nea 列，且 eaf 列同名）
    out_rows = load_assoc(IV / out_file, nea_col=out_nea)
    bx, sx, by, sy, used = build_arrays(exp_rows, out_rows, snps)
    # 验证：全集 IVW 与 JSON 一致？
    bx, sx, by, sy = np.array(bx), np.array(sx), np.array(by), np.array(sy)
    fx = eng.ivw(bx, by, sy, model="fixed")
    rd = eng.ivw(bx, by, sy, model="random")
    ok_f = close(fx["beta"], ref["ivw_fixed"]["beta"]) and close(fx["se"], ref["ivw_fixed"]["se"])
    ok_r = close(rd["beta"], ref["ivw_random"]["beta"]) and close(rd["se"], ref["ivw_random"]["se"])
    print(f"  重构 n={len(used)}；全集 IVW固定 β={fx['beta']:.6f} vs JSON {ref['ivw_fixed']['beta']:.6f} {'OK' if ok_f else 'MISMATCH'}")
    print(f"  全集 IVW随机 SE={rd['se']:.6f} vs JSON {ref['ivw_random']['se']:.6f} {'OK' if ok_r else 'MISMATCH'}")
    if not (ok_f and ok_r):
        print("  重构不忠实，停止该方向"); continue
    # 剔离群后（PRESSO 离群索引来自 05 JSON）
    out_idx = ref["mr_presso"]["outliers"]["outlier_indices"]
    mask = np.array([i for i in range(len(bx)) if i not in out_idx])
    cl = eng.ivw(bx[mask], by[mask], sy[mask], model="fixed")
    print(f"  剔离群 {len(out_idx)} 个 → n={len(mask)}：β={cl['beta']:.6f} (vs 05 JSON distortion beta_outliers_removed {ref['mr_presso']['distortion']['beta_outliers_removed']:.6f})")
    print(f"  SE={cl['se']:.6f}  95%CI=[{cl['ci_lower']:.6f}, {cl['ci_upper']:.6f}]  P={cl['p']:.6g}")
    kept = [used[i] for i in mask]
    print(f"  保留 SNP: {kept}")
    json.dump({"tag": tag, "n_kept": len(mask), "kept_snps": kept,
               "removed_snps": [used[i] for i in out_idx],
               "ivw_fixed_clean": cl},
              open(BASE / f"presso-clean-{tag}.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
print("完成：结果写入 05-分析/presso-clean-*.json")
