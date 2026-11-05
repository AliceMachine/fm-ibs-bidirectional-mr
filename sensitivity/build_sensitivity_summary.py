# -*- coding: utf-8 -*-
"""
build_sensitivity_summary.py —— 06 阶段汇总器。

把 05 阶段 7 份 JSON + 06 补做的 Steiger（steiger-d2-completed.json）
+ 共定位（coloc/coloc-results.json，若已生成）自动汇总为 06-敏感性/sensitivity.json。

原则：所有数字从 JSON 机器读出，不做任何手打（沿用 gen_results.py 的教训）。
"""
import json
from pathlib import Path

BASE = Path(r"D:\2027\workbuddy working\021\MR演练")
A05 = BASE / "05-分析"
A06 = BASE / "06-敏感性"

ROUND = 6


def rd(x, n=ROUND):
    if x is None:
        return None
    try:
        return round(float(x), n)
    except (TypeError, ValueError):
        return x


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def univariate_block(path, label, direction, exposure, outcome, steiger_fill=None):
    d = load(path)
    snps = d["snps"]
    pr = d["mr_presso"]
    out_idx = pr["outliers"]["outlier_indices"]
    loo = d["leave_one_out"]
    loo_named = []
    for e in loo:
        i = e["index"]
        loo_named.append({
            "excluded_rsid": snps[i],
            "ivw_beta_without": rd(e["beta"]),
            "ivw_se_without": rd(e["se"]),
        })
    beta_lo, beta_hi = None, None
    if loo_named:
        betas = [x["ivw_beta_without"] for x in loo_named]
        beta_lo, beta_hi = min(betas), max(betas)

    st = d["steiger"]
    if "skipped" in st and steiger_fill:
        st_final = steiger_fill
        st_source = "06-补做（run_steiger_fill.py，补齐 FM 样本量后重算）"
    else:
        st_final = st
        st_source = "05-主分析（05-分析/run_mvmr.py 调用引擎自带）"

    model_name = d["primary_analysis"]  # 如 "ivw_random"
    est = d[model_name]
    dist = pr.get("distortion") or {}
    presso_out = {
        "global_p": rd(pr["global"]["p_global"]),
        "n_boot": pr["global"]["n_boot"],
        "outlier_rsids": [snps[i] for i in out_idx],
        "n_outliers": pr["outliers"]["n_outliers"],
        "beta_outlier_corrected": rd(dist.get("beta_outliers_removed")),
        "p_distortion": rd(dist.get("p_distortion")),
        "implementation": pr.get("implementation_note",
                                 "Simplified re-implementation (see Methods)"),
    }
    eg = d["mr_egger"]
    return {
        "label": label,
        "direction": direction,
        "exposure": exposure,
        "outcome": outcome,
        "n_snp": d["n_snp"],
        "snps": snps,
        "f_statistic_range": rd(d["f_statistic"]["min"]) if isinstance(d.get("f_statistic"), dict)
        else rd(d.get("f_statistic")),
        "primary": {
            "model": model_name,
            "model_rule": d["primary_model_rule"],
            "beta": rd(est["beta"]),
            "se": rd(est["se"]),
            "or": rd(math.exp(est["beta"]), 4),
            "ci_lower_or": rd(math.exp(est["ci_lower"]), 4),
            "ci_upper_or": rd(math.exp(est["ci_upper"]), 4),
            "p": rd(est["p"]),
        },
        "heterogeneity": {
            "cochran_q": rd(d["heterogeneity"]["cochran_q"]),
            "df": d["heterogeneity"]["df"],
            "p": rd(d["heterogeneity"]["p"]),
            "i_squared_percent": rd(d["heterogeneity"]["i_squared_percent"], 2),
        },
        "egger": {
            "intercept_beta": rd(eg["intercept"]["beta"]),
            "intercept_p": rd(eg["intercept"]["p"]),
            "slope_p": rd(eg["slope"]["p"]),
            "wording_rule": "截距 P>0.05 只能表述为'未检出方向性多效性证据（该检验在此功效有限）'",
        },
        "weighted_median_beta": rd(d["weighted_median"]["beta"]) if d.get("weighted_median") else None,
        "weighted_mode_beta": rd(d["weighted_mode"]["beta"]) if d.get("weighted_mode") else None,
        "mr_presso": presso_out,
        "leave_one_out": {
            "per_snp": loo_named,
            "beta_range_without_any_single_snp": [rd(beta_lo), rd(beta_hi)],
        },
        "steiger": {
            "source": st_source,
            "prop_correct": rd(st_final.get("prop_correct"), 4),
            "n_correct": st_final.get("n_correct"),
            "n_evaluated": st_final.get("n_evaluated"),
            "wrong_direction_rsids": [p["rsid"] for p in st_final.get("per_snp", [])
                                      if not p["direction_correct"]],
            "filtering_applied": st_final.get("filtering_applied") if
            isinstance(st_final, dict) and "filtering_applied" in st_final
            else (bool(st_final.get("per_snp") and
                       any(not p["direction_correct"] for p in st_final["per_snp"]))),
            "filtered_ivw": st_final.get("filtered_ivw") if isinstance(st_final, dict) else None,
            "note": "解释方差方向的支持性证据，不是因果方向的证明（预注册 §6.2 / 技能红线）",
        },
    }


import math  # noqa: E402  (放底部仅供上方 OR 计算使用——Python 顶层顺序可接受)


def mvmr_block(path):
    d = load(path)
    model_name = d.get("primary")
    est = d.get(model_name) if isinstance(model_name, str) else None
    out = {
        "tag": d["tag"],
        "exposures": d["exposures"],
        "outcome": d["outcome"],
        "n_snps": d["n_snps"],
        "n_dropped_incomplete": d.get("n_dropped"),
        "q_pval": rd(d.get("q_pval")),
        "conditional_F": {k: rd(v, 1) for k, v in d.get("f_conditional", {}).items()},
        "primary_model": d.get("primary_rule") or model_name,
    }
    if isinstance(est, dict):
        for k, v in est.items():
            out["primary_" + k] = rd(v) if isinstance(v, float) else v
    # 随机效应对照
    if d.get("estimates_random"):
        out["ivw_random"] = {k: rd(v) for k, v in d["estimates_random"].items()
                             if isinstance(v, (int, float))}
    return out


def main():
    sf = load(A06 / "steiger-d2-completed.json")

    blocks = {
        "d1_IBS_FM": univariate_block(A05 / "d1-IBS-FM-univariate.json",
                                      "方向一：IBS→FM", "IBS→FM", "IBS", "FM"),
        "d2_FM_IBS": univariate_block(A05 / "d2-FM-IBS-univariate.json",
                                      "方向二：FM→IBS", "FM→IBS", "FM", "IBS",
                                      steiger_fill=sf["d2-FM-IBS"]["steiger"]
                                      | {"filtering_applied": bool(sf["d2-FM-IBS"]["steiger_filtering"]["filtered_ivw"]),
                                         "filtered_ivw": sf["d2-FM-IBS"]["steiger_filtering"]["filtered_ivw"]}),
        "d2_finnrep": univariate_block(A05 / "d2-finnrep-FM-IBS.json",
                                       "方向二外部复制：FM→IBS（FinnGen）", "FM→IBS",
                                       "FM", "IBS(FinnGen)",
                                       steiger_fill=sf["d2-finnrep-FM-IBS"]["steiger"]
                                       | {"filtering_applied": bool(sf["d2-finnrep-FM-IBS"]["steiger_filtering"]["filtered_ivw"]),
                                          "filtered_ivw": sf["d2-finnrep-FM-IBS"]["steiger_filtering"]["filtered_ivw"]}),
    }

    mvmrs = {tag: mvmr_block(A05 / f"{tag}.json")
             for tag in ["d1-mv4-mvmr", "d1-mv3-mvmr", "d2-mv4-mvmr", "d2-mv3-mvmr"]}

    coloc_path = A06 / "coloc" / "coloc-results.json"
    coloc = None
    if coloc_path.exists():
        raw = load(coloc_path)
        coloc = {
            "method": "coloc.abf（纯 Python 自建实现，公式与 R coloc 一致；自测 4/4）",
            "priors": {"p1": 1e-4, "p2": 1e-4, "p12": 1e-5, "W_cc": 0.04},
            "window": "IV ±500kb",
            "threshold_PP_H4": 0.8,
            "per_iv": raw,
            "n_iv_with_PP_H4_gt_0.8": sum(
                1 for v in raw.values() if v.get("PP", {}).get("H4", 0) > 0.8),
        }

    summary = {
        "engine": "mr_engine.py（纯 Python 自建）+ run_steiger_fill.py + coloc_abf.py",
        "note_presso": "MR-PRESSO 为简化重实现（残差全局检验+置换零分布），非官方 R 包，论文方法节须声明",
        "univariate": blocks,
        "mvmr": mvmrs,
        "colocalisation": coloc,
        "not_done": [
            {"analysis": "中介 MR",
             "reason": "预注册 §6.2 明确不做（超出本题范围，防范围蔓延）"},
            {"analysis": "条件性外部失眠库复制（Lane 2019）",
             "reason": "双源核实不可得：完整 meta 被 23andMe 条款锁定；公开 UKB 子集与主失眠库 ukb-b-3957 同队列（不构成独立复制）"},
        ],
    }

    out = A06 / "sensitivity.json"
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print("saved:", out)
    for k, b in blocks.items():
        print(k, "-> primary OR", b["primary"]["or"],
              "PRESSO outliers:", b["mr_presso"]["outlier_rsids"],
              "Steiger:", b["steiger"]["n_correct"], "/", b["steiger"]["n_evaluated"])


if __name__ == "__main__":
    main()
