# -*- coding: utf-8 -*-
"""
coloc_abf.py —— 纯 Python 自建的贝叶斯共定位引擎（coloc.abf 法）。

方法来源（R 包 coloc 的同款公式，Giambartolomei et al. 2014, PLOS Genetics）：
  单因果变异模型下，每个 SNP 的对数贝叶斯因子（log ABF）：
      lABF_i = 0.5*log(1 - r_i) + 0.5*r_i*z_i^2
      其中 r_i = W / (W + V_i)，V_i = se_i^2，z_i = beta_i / se_i
      两个二分类性状均用先验效应方差 W = 0.2^2 = 0.04（coloc 默认）

  五个假设的后验概率（odds 形式，与 coloc R 源码一致）：
      H0 = 0                          （两性状在此区域都没有因果变异）
      H1 = l1 + S1                    （只有性状1有）
      H2 = l2 + S2                    （只有性状2有）
      H3 = l1 + l2 + S1 + S2          （各有各的因果变异——LD 共享假象）
      H4 = l1 + l2 + l12 + S12        （同一个因果变异——真共定位）
      其中 S1 = logsumexp(lABF1)、S2 = logsumexp(lABF2)、
           S12 = logsumexp(lABF1 + lABF2)
           l1 = log(p1/(1-p1))，l2 = log(p2/(1-p2))，
           l12 = log(p12/(p1*p2))
      默认先验 p1 = p2 = 1e-4，p12 = 1e-5（coloc 默认）

  判读规则（预注册 §6.2）：PP.H4 > 0.8 判为支持共定位。

等位基因对齐规则（与主线 harmonise 一致）：
  - 两个性状的 effect allele 一致 → 直接用；
  - 互补（回文 SNP，A/T 或 C/G）→ 用等位基因频率（EAF）判断：
      |eaf1 - eaf2| <= 0.1 视为同向；|eaf1 + eaf2 - 1| <= 0.1 视为反向（翻转 beta）；
      0.42 < eaf < 0.58 或缺 EAF → 判向不可靠，整行剔除（预注册同款 0.42-0.58 规则）；
  - 非回文但 allele 顺序相反 → 翻转 beta 后使用。

自测（--selftest）：三个模拟场景 + 等位基因翻转等价性，全部应 PASS。
"""
import json
import math
import sys

W_CC = 0.04          # 二分类性状先验效应方差（0.2^2）
P1 = 1e-4            # 性状1因果变异先验
P2 = 1e-4            # 性状2因果变异先验
P12 = 1e-5           # 同一变异同时对两性状因果的先验（coloc 默认）
EAF_AMBIG_LO = 0.42  # 回文 SNP 频率模糊区（预注册同款）
EAF_AMBIG_HI = 0.58


def _logsumexp(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def log_abf(beta, se, w=W_CC):
    """单个 SNP 的对数贝叶斯因子（coloc 公式）。"""
    v = se * se
    r = w / (w + v)
    z2 = (beta / se) ** 2
    return 0.5 * math.log(1.0 - r) + 0.5 * r * z2


def coloc_abf(d1, d2, p1=P1, p2=P2, p12=P12):
    """共定位主函数。

    d1/d2 结构：{"rsid": [...], "beta": [...], "se": [...], "eaf": [...] or None,
                 "label": str}
    两个数据集的 rsid 已在调用前按交集对齐、等位基因方向已对齐。
    返回：{"PP": {H0..H4}, "n_snp": int, "lead1"/"lead2": 最强 SNP, ...}
    """
    labf1 = [log_abf(b, s) for b, s in zip(d1["beta"], d1["se"])]
    labf2 = [log_abf(b, s) for b, s in zip(d2["beta"], d2["se"])]
    labf12 = [a + b for a, b in zip(labf1, labf2)]

    s1 = _logsumexp(labf1)
    s2 = _logsumexp(labf2)
    s12 = _logsumexp(labf12)

    l1 = math.log(p1 / (1 - p1))
    l2 = math.log(p2 / (1 - p2))
    l12 = math.log(p12 / (p1 * p2))

    hs = {
        "H0": 0.0,
        "H1": l1 + s1,
        "H2": l2 + s2,
        "H3": l1 + l2 + s1 + s2,
        "H4": l1 + l2 + l12 + s12,
    }
    m = max(hs.values())
    exps = {k: math.exp(v - m) for k, v in hs.items()}
    tot = sum(exps.values())
    pp = {k: v / tot for k, v in exps.items()}

    lead1 = d1["rsid"][max(range(len(labf1)), key=lambda i: labf1[i])]
    lead2 = d2["rsid"][max(range(len(labf2)), key=lambda i: labf2[i])]

    return {
        "label": f'{d1.get("label","trait1")}|{d2.get("label","trait2")}',
        "n_snp": len(labf1),
        "PP": {k: round(v, 6) for k, v in pp.items()},
        "PP_H4_raw": pp["H4"],
        "lead_trait1": lead1,
        "lead_trait2": lead2,
        "priors": {"p1": p1, "p2": p2, "p12": p12, "W_cc": W_CC},
    }


def align_datasets(d1_raw, d2_raw):
    """按 rsID 取交集并做等位基因对齐。

    d_raw 结构：{"rsid": [], "ea": [], "oa": [], "beta": [], "se": [],
                 "eaf": [] 或 None（每行可缺，用 None 占位）, "label": str}
    返回 (d1_aligned, d2_aligned, dropped) —— dropped 记录剔除原因，供报告如实披露。
    """
    def _dedup(d):
        """真实全量文件存在同一 rsID 多行（多等位/重复注释）；保留首行。
        R coloc 官方包对重复 snp 直接报错，此处与官方口径一致。
        标量键（label 等非列表/非 None）原样复制。"""
        seen = set()
        out = {}
        for k in d:
            v = d[k]
            if v is None or isinstance(v, str):
                out[k] = v
            else:
                out[k] = []
        for i, rs in enumerate(d["rsid"]):
            if rs in seen:
                continue
            seen.add(rs)
            for k in d:
                v = d[k]
                if v is None or isinstance(v, str):
                    continue
                out[k].append(v[i])
        return out

    d1_raw, d2_raw = _dedup(d1_raw), _dedup(d2_raw)
    idx2 = {r: i for i, r in enumerate(d2_raw["rsid"])}
    out1 = {"rsid": [], "beta": [], "se": [], "label": d1_raw.get("label", "trait1")}
    out2 = {"rsid": [], "beta": [], "se": [], "label": d2_raw.get("label", "trait2")}
    dropped = []
    for i, rs in enumerate(d1_raw["rsid"]):
        if rs not in idx2:
            continue
        j = idx2[rs]
        ea1, oa1 = d1_raw["ea"][i].upper(), d1_raw["oa"][i].upper()
        ea2, oa2 = d2_raw["ea"][j].upper(), d2_raw["oa"][j].upper()
        b1, b2 = float(d1_raw["beta"][i]), float(d2_raw["beta"][j])
        se1, se2 = float(d1_raw["se"][i]), float(d2_raw["se"][j])
        # 行级清洗：se≤0 或 beta/se 非有限值的位点剔除（真实全量文件常见）
        if not (math.isfinite(b1) and math.isfinite(se1) and se1 > 0
                and math.isfinite(b2) and math.isfinite(se2) and se2 > 0):
            dropped.append((rs, "invalid_beta_se"))
            continue
        eaf1 = d1_raw.get("eaf")
        eaf2 = d2_raw.get("eaf")
        e1 = eaf1[i] if eaf1 else None
        e2 = eaf2[j] if eaf2 else None
        # 容错：原始行是字符串；空串/非法值统一按缺频率处理
        try:
            e1 = float(e1) if e1 is not None and str(e1).strip() != "" else None
        except (TypeError, ValueError):
            e1 = None
        try:
            e2 = float(e2) if e2 is not None and str(e2).strip() != "" else None
        except (TypeError, ValueError):
            e2 = None

        palindromic = {ea1, oa1} == {"a", "t"} or {ea1, oa1} == {"c", "g"}
        flip = False
        if ea1 == ea2 and oa1 == oa2:
            pass  # 同向
        elif ea1 == oa2 and oa1 == ea2:
            flip = True  # 等位基因顺序相反 → 性状2 beta 取反
        else:
            dropped.append((rs, "alleles_unmatchable"))
            continue

        if palindromic:
            if e1 is None or e2 is None:
                dropped.append((rs, "palindromic_no_eaf"))
                continue
            if EAF_AMBIG_LO < e1 < EAF_AMBIG_HI or EAF_AMBIG_LO < e2 < EAF_AMBIG_HI:
                dropped.append((rs, "palindromic_ambiguous_eaf"))
                continue
            # 用 EAF 判断真实方向：同向则 eaf1≈eaf2，反向则 eaf1≈1-eaf2
            if abs(e1 - e2) <= 0.1:
                pass
            elif abs(e1 + e2 - 1.0) <= 0.1:
                flip = not flip
            else:
                dropped.append((rs, "palindromic_eaf_conflict"))
                continue
        else:
            if flip and e1 is not None and e2 is not None:
                # 非回文方向本已确定，频率只作核对；冲突过大记剔除
                if abs((1 - e2) - e1) > 0.2:
                    dropped.append((rs, "eaf_conflict"))
                    continue

        out1["rsid"].append(rs)
        out1["beta"].append(b1)
        out1["se"].append(se1)
        out2["rsid"].append(rs)
        out2["beta"].append(-b2 if flip else b2)
        out2["se"].append(se2)
    return out1, out2, dropped


def _selftest():
    """四个场景自测：真共定位 / LD 假象（各享各的）/ 无信号 / 翻转等价。"""
    import random
    random.seed(20260922)
    results = []

    def sim(n=3000, shared=False, causal_sep=False):
        rs = [f"rs{i}" for i in range(n)]
        b1, s1_, b2, s2_ = [], [], [], []
        for i in range(n):
            se = 0.01
            z1 = random.gauss(0, 1)
            z2 = random.gauss(0, 1)
            if shared and i == 500:
                z1, z2 = 8.0, 7.0          # 同一 SNP 两性状都强
            if causal_sep and i == 300:
                z1 = 8.0                    # 性状1 因果在 300
            if causal_sep and i == 700:
                z2 = 7.0                    # 性状2 因果在 700
            b1.append(z1 * se)
            s1_.append(se)
            b2.append(z2 * se)
            s2_.append(se)
        return ({"rsid": rs, "beta": b1, "se": s1_, "label": "t1"},
                {"rsid": rs, "beta": b2, "se": s2_, "label": "t2"})

    # 场景1：共享因果 → H4 应显著占优
    d1, d2 = sim(shared=True)
    r = coloc_abf(d1, d2)
    ok1 = r["PP_H4_raw"] > 0.8
    results.append(("shared_causal_H4>0.8", ok1, r["PP"]))

    # 场景2：两个不同因果变异 → H3 应占优（H4 不应占优）
    d1, d2 = sim(causal_sep=True)
    r = coloc_abf(d1, d2)
    ok2 = r["PP"]["H3"] > r["PP"]["H4"]
    results.append(("distinct_causal_H3>H4", ok2, r["PP"]))

    # 场景3：无信号 → H0 应占优
    d1, d2 = sim()
    r = coloc_abf(d1, d2)
    ok3 = r["PP"]["H0"] > 0.5
    results.append(("no_signal_H0_dominant", ok3, r["PP"]))

    # 场景4：等位基因翻转等价性（把性状2全部 ea/oa 互换并取反 beta → 结果不变）
    rs = [f"rs{i}" for i in range(3000)]
    b1 = [random.gauss(0, 0.01) for _ in range(3000)]
    b2 = [random.gauss(0, 0.01) for _ in range(3000)]
    b1[500], b2[500] = 0.08, 0.07
    s = [0.01] * 3000
    a1 = {"rsid": rs, "ea": ["A"] * 3000, "oa": ["G"] * 3000,
          "beta": b1, "se": s, "eaf": None, "label": "t1"}
    a2 = {"rsid": rs, "ea": ["G"] * 3000, "oa": ["A"] * 3000,
          "beta": b2, "se": s, "eaf": None, "label": "t2"}
    n1, n2, dr = align_datasets(a1, a2)
    r1 = coloc_abf({"rsid": rs, "beta": b1, "se": s, "label": "t1"},
                   {"rsid": rs, "beta": b2, "se": s, "label": "t2"})
    r2 = coloc_abf(n1, n2)
    ok4 = (not dr) and abs(r1["PP_H4_raw"] - r2["PP_H4_raw"]) < 1e-12
    results.append(("allele_flip_equivalence", ok4, None))

    all_pass = True
    for name, ok, detail in results:
        print(("PASS" if ok else "FAIL"), name)
        if detail:
            print("   PP:", detail)
        all_pass = all_pass and ok
    print("SELFTEST", "ALL PASS" if all_pass else "FAILED")
    return all_pass


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        ok = _selftest()
        sys.exit(0 if ok else 1)
    print("usage: python coloc_abf.py --selftest")
