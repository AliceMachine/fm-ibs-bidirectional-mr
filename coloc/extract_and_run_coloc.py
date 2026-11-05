# -*- coding: utf-8 -*-
"""
extract_and_run_coloc.py —— 共定位数据抽取与运行（06 阶段）。

做什么：
  1. 校验 IBS 全量统计量 md5（04-数据/ibs_full/GCST_IBS_FINAL.gz）
  2. 从 IBS 文件第一遍扫描定位 d1/d2 全部 IV 的 GRCh38 位置
  3. 按 IV 位置开 ±500kb 窗口（预注册共定位常规窗口），第二遍扫描抽取区域统计量
  4. FM 端（本地全量 GCST90838603.h.tsv.gz，GRCh38）同样抽取
  5. 逐 IV 跑 coloc.abf（调用同目录 coloc_abf.py），输出 JSON

诚实声明（会写进报告）：
  - FM 病例数未按欧洲祖先亚组回填，但 coloc.abf 在已给 se 时不使用病例占比，
    因此不影响后验概率计算（coloc 公式性质，非近似）。
  - 区域限制在 ±500kb；因果变异若在该窗口之外，coloc 会倾向给出 H1/H2/H3。

运行前提：IBS 下载完成且 MD5_OK。
"""
import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

from coloc_abf import align_datasets, coloc_abf

BASE = Path(r"D:\2027\workbuddy working\021\MR演练")
IBS_GZ = BASE / "04-数据" / "ibs_full" / "GCST_IBS_FINAL.gz"
IBS_MD5_EXPECT = "004f576988d426a3e151000b6123478d"
FM_GZ = BASE / "04-数据" / "fm_download" / "GCST90838603.h.tsv.gz"
OUT_DIR = BASE / "06-敏感性" / "coloc"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WINDOW = 500_000  # ±500kb

# d1 分析实际使用的 5 个 IV（IBS 暴露侧；来自 05 d1-IBS-FM-univariate.json）
D1_IVS = ["rs426444", "rs10156602", "rs12549729", "rs9513519", "rs7106434"]
# d2 的 15 个 IV（FM 暴露侧）——位置直接读 fm_iv_table.csv（GRCh38）


def md5_of(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_d2_ivs():
    rows = list(csv.DictReader(open(BASE / "04-数据" / "fm_download" / "fm_iv_table.csv")))
    return [(r["rsid"], r["chr"], int(r["pos"])) for r in rows]


def map_ibs_columns(header):
    """Catalog pre-GWAS-SSF 表头 → 统一字段名。返回 {统一名: 列号}，缺关键列则报错。"""
    names = [h.strip().lower() for h in header.rstrip("\n").split("\t")]
    def find(cands):
        for c in cands:
            if c in names:
                return names.index(c)
        return None
    colmap = {
        "rsid": find(["variant_id", "rsid", "rs_id", "snp", "h_rs_id"]),
        "chr": find(["chromosome", "chr", "h_chromosome"]),
        "pos": find(["base_pair_location", "pos", "position", "h_position"]),
        "ea": find(["effect_allele", "ea", "hm_effect_allele"]),
        "oa": find(["other_allele", "oa", "hm_other_allele"]),
        "beta": find(["beta", "hm_beta", "effect"]),
        "se": find(["standard_error", "se", "hm_standard_error"]),
        "eaf": find(["effect_allele_frequency", "eaf", "hm_effect_allele_frequency",
                     "effect_allele_af"]),
    }
    missing = [k for k in ("rsid", "chr", "pos", "ea", "oa", "beta", "se") if colmap[k] is None]
    if missing:
        raise RuntimeError(f"IBS 文件缺少关键列: {missing}; 实际表头: {names}")
    return colmap


def scan_ibs(windows=None, find_rsids=None, colmap=None):
    """单遍扫描 IBS 全量文件。

    windows: [(chr, start, end), ...] —— 收集窗口内行
    find_rsids: set(rsid) —— 返回 {rsid: (chr, pos)}（首个遇到的位置）
    """
    hits = {}
    rows = []
    wset = {(c, s, e) for c, s, e in windows} if windows else None
    with gzip.open(IBS_GZ, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        cm = colmap or map_ibs_columns(header)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if find_rsids:
                rs = p[cm["rsid"]]
                if rs in find_rsids and rs not in hits:
                    hits[rs] = (p[cm["chr"]], int(p[cm["pos"]]))
            if wset:
                c = p[cm["chr"]]
                try:
                    pos = int(p[cm["pos"]])
                except ValueError:
                    continue
                for (ws, we) in by_chr.get(c, ()):
                    if ws <= pos <= we:
                        rows.append(p)
                        break
    return hits, rows, cm


by_chr = {}


def set_windows(windows):
    """把窗口列表按染色体索引，加速扫描匹配。"""
    global by_chr
    by_chr = {}
    for c, s, e in windows:
        by_chr.setdefault(str(c), []).append((s, e))


def scan_fm(windows):
    """单遍扫描 FM 全量文件（读表头自适应列序，防硬编码错位）。"""
    set_windows(windows)
    rows = []
    with gzip.open(FM_GZ, "rt", encoding="utf-8", errors="replace") as f:
        names = [h.strip().lower() for h in f.readline().rstrip("\n").split("\t")]
        cm = {
            "rsid": names.index("rsid"),
            "chr": names.index("chromosome"),
            "pos": names.index("base_pair_location"),
            "ea": names.index("effect_allele"),
            "oa": names.index("other_allele"),
            "beta": names.index("beta"),
            "se": names.index("standard_error"),
            "eaf": names.index("effect_allele_frequency"),
        }
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(names):
                continue
            try:
                pos = int(p[cm["pos"]])
            except ValueError:
                continue
            for (ws, we) in by_chr.get(p[cm["chr"]], ()):
                if ws <= pos <= we:
                    rows.append(p)
                    break
    return rows, cm


def rows_to_dataset(rows, colidx, label):
    d = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [], "label": label}
    for p in rows:
        d["rsid"].append(p[colidx["rsid"]])
        d["ea"].append(p[colidx["ea"]])
        d["oa"].append(p[colidx["oa"]])
        d["beta"].append(p[colidx["beta"]])
        d["se"].append(p[colidx["se"]])
        d["eaf"].append(p[colidx["eaf"]] if colidx["eaf"] is not None else None)
    return d


def main():
    # 1. md5 校验
    md5 = md5_of(IBS_GZ)
    print("IBS md5:", md5)
    if md5 != IBS_MD5_EXPECT:
        print("MD5_FAIL —— 下载不完整，终止")
        sys.exit(1)
    print("MD5_OK")

    d2_ivs = load_d2_ivs()
    d1_rsids = set(D1_IVS)
    d2_rsids = {r[0] for r in d2_ivs}

    # 2. 第一遍扫描 IBS：定位 d1 IV 的 GRCh38 位置（d2 IV 位置用 fm_iv_table 的 GRCh38）
    set_windows([])  # 不收集窗口
    hits, _, cm = scan_ibs(find_rsids=d1_rsids | d2_rsids)
    print("IBS 第一遍：找到", len(hits), "/", len(d1_rsids | d2_rsids), "个 IV")
    missing = (d1_rsids | d2_rsids) - set(hits)
    if missing:
        print("警告：IBS 文件中未找到（将如实记录）:", sorted(missing))

    # d2 IV 位置：优先用本地 fm_iv_table（GRCh38 实测），d1 用 IBS 文件定位
    iv_pos = {rs: ("d2", c, p) for rs, c, p in d2_ivs}
    for rs in D1_IVS:
        if rs in hits:
            iv_pos[rs] = ("d1", hits[rs][0], hits[rs][1])

    # 3. 构建窗口（方向分别保留自己的 IV 窗口；合并重叠窗口以减少扫描量）
    all_windows = []
    for rs, (tag, c, p) in iv_pos.items():
        all_windows.append((str(c), max(1, p - WINDOW), p + WINDOW))
    set_windows(all_windows)

    # 4. 两文件各一遍扫描收集窗口行
    print("扫描 IBS 区域行 ...")
    _, ibs_rows, cm = scan_ibs(windows=all_windows, colmap=cm)
    print("IBS 区域行:", len(ibs_rows))
    print("扫描 FM 区域行 ...")
    fm_rows, fm_cm = scan_fm(all_windows)
    print("FM 区域行:", len(fm_rows))

    # IBS 列索引（row 是 list）→ FM 列索引（表头自适应）
    fm_colidx = fm_cm
    ibs_colidx = {k: cm[k] for k in cm}

    fm_all = rows_to_dataset(fm_rows, fm_colidx, "FM")
    ibs_all = rows_to_dataset(ibs_rows, ibs_colidx, "IBS")

    # 5. 逐 IV 共定位
    results = {}
    for rs, (tag, c, p) in sorted(iv_pos.items()):
        lo, hi = p - WINDOW, p + WINDOW
        # 该 IV 窗口内的行（用位置过滤已收集行）
        def subset(ds):
            out = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [],
                   "label": ds["label"]}
            # 位置与 rsid 对齐：重新从原始行过滤
            return out
        # 为省事：把整条染色体的窗口行都给 align（交集+同窗），再按窗口过滤
        fm_sub = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [], "label": "FM"}
        ibs_sub = {"rsid": [], "ea": [], "oa": [], "beta": [], "se": [], "eaf": [], "label": "IBS"}
        # FM 行列序按表头自适应（fm_cm）
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
        a1, a2, dropped = align_datasets(fm_sub, ibs_sub)
        if len(a1["rsid"]) < 50:
            results[rs] = {"tag": tag, "status": "insufficient_overlap",
                           "n_fm": len(fm_sub["rsid"]), "n_ibs": len(ibs_sub["rsid"]),
                           "n_aligned": len(a1["rsid"]), "dropped": dropped}
            print(rs, "insufficient overlap", len(a1["rsid"]))
            continue
        r = coloc_abf(a1, a2)
        r.update({"tag": tag, "iv": rs, "chr": c, "pos": p,
                  "window": [lo, hi], "n_dropped_alleles": len(dropped),
                  "dropped_examples": dropped[:10],
                  "iv_direction_correct": rs not in ("rs62100765", "rs2587363")})
        results[rs] = r
        print(rs, tag, "n=", r["n_snp"], "PP.H4=", r["PP"]["H4"],
              "coloc:", r["PP"]["H4"] > 0.8)

    out = OUT_DIR / "coloc-results.json"
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    n_ok = sum(1 for v in results.values() if v.get("PP", {}).get("H4", 0) > 0.8)
    print(f"完成：{sum(1 for v in results.values() if 'PP' in v)} 个 IV 出结果，"
          f"其中 PP.H4>0.8 的 {n_ok} 个")
    print("saved:", out)


if __name__ == "__main__":
    main()
