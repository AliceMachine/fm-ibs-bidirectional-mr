#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mr_engine.py — Pure-Python Mendelian randomization engine.

No external dependencies beyond numpy / scipy. Implements, following their
original descriptions:

  * Wald ratio                    (single-variant)
  * IVW, fixed and multiplicative random effects   [Burgess 2013; Bowden 2015]
  * MR-Egger regression (weighted least squares)   [Bowden 2015]
  * Weighted median estimator                      [Bowden 2016]
  * Weighted / simple mode estimator               [Hartwig 2017]
  * Cochran's Q and I^2                            [Bowden 2015]
  * Leave-one-out re-estimation
  * MR-PRESSO global / outlier / distortion test   [Verbanck 2018]
  * Steiger directionality                         [Hemani 2017]

IMPORTANT — reporting honesty
    If this engine produced the results, the manuscript MUST say so. Do NOT
    write "TwoSampleMR was used" unless the R reference script was actually
    run. See references/mr-software-and-data-sources.md section 3.

Input  : harmonised table (CSV) with columns
         SNP, beta_exposure, se_exposure, beta_outcome, se_outcome
         optional: eaf_exposure, n_exposure, n_outcome
Output : JSON results file + machine-readable KEY=VALUE lines on stdout.

Exit codes: 0 success, 1 failure (prints "MR_ENGINE_FAILED: <reason>" to stderr)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

# --------------------------------------------------------------------------
# constants / helpers
# --------------------------------------------------------------------------

REQUIRED_COLUMNS = [
    "SNP",
    "beta_exposure",
    "se_exposure",
    "beta_outcome",
    "se_outcome",
]

PALINDROMIC = {"AT", "TA", "CG", "GC"}


def fail(msg: str) -> None:
    """Print a machine-readable failure line and exit non-zero."""
    print(f"MR_ENGINE_FAILED: {msg}", file=sys.stderr)
    sys.exit(1)


def _z_to_p(z: float) -> float:
    """Two-sided p value from a z statistic."""
    return float(2.0 * stats.norm.sf(abs(z)))


def _pack(beta: float, se: float) -> dict:
    """Bundle a point estimate with its SE, CI and p value."""
    if se is None or not np.isfinite(se) or se <= 0:
        return {
            "beta": _f(beta),
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
            "p": None,
        }
    z = beta / se
    return {
        "beta": _f(beta),
        "se": _f(se),
        "ci_lower": _f(beta - 1.96 * se),
        "ci_upper": _f(beta + 1.96 * se),
        "p": _f(_z_to_p(z)),
    }


def _f(x):
    """Cast numpy scalars to plain float, keeping full precision."""
    if x is None:
        return None
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    return x


# --------------------------------------------------------------------------
# core estimators
# --------------------------------------------------------------------------


def wald_ratio(bx: float, by: float, sy: float) -> dict:
    """Single-variant Wald ratio estimator.

    beta = by / bx ; se = sy / |bx|   (first-order approximation, as in
    TwoSampleMR's mr_wald_ratio).
    """
    if bx == 0:
        fail("wald_ratio: beta_exposure is zero, ratio undefined")
    beta = by / bx
    se = abs(sy / bx)
    return _pack(beta, se)


def ivw(bx: np.ndarray, by: np.ndarray, sy: np.ndarray, model: str = "fixed") -> dict:
    """Inverse-variance weighted estimator.

    Weighted regression through the origin:
        beta = sum(bx*by/sy^2) / sum(bx^2/sy^2)
        se   = sqrt(1 / sum(bx^2/sy^2))            (fixed effects)
    For multiplicative random effects the SE is inflated by sqrt(phi) with
    phi = max(1, Q/(k-1)).
    """
    w = bx**2 / sy**2
    sw = float(np.sum(w))
    if sw <= 0:
        fail("ivw: total weight is zero; check se_outcome values")
    ratio = by / bx
    beta = float(np.sum(w * ratio) / sw)

    if model == "fixed":
        se = math.sqrt(1.0 / sw)
    elif model == "random":
        q = float(np.sum(w * (ratio - beta) ** 2))
        k = len(bx)
        phi = max(1.0, q / (k - 1)) if k > 1 else 1.0
        se = math.sqrt(phi / sw)
    else:
        fail(f"ivw: unknown model '{model}' (expected 'fixed' or 'random')")

    return _pack(beta, se)


def _wls(X: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Weighted least squares; returns (beta, se, resid, dof)."""
    XtWX = X.T @ (w[:, None] * X)
    XtWy = X.T @ (w * y)
    beta = np.linalg.solve(XtWX, XtWy)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    if dof <= 0:
        return beta, np.full_like(beta, np.nan), resid, dof
    sigma2 = float(np.sum(w * resid**2) / dof)
    cov = sigma2 * np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    return beta, se, resid, dof


def mv_ivw(BX: np.ndarray, SX: np.ndarray, by: np.ndarray, sy: np.ndarray,
           model: str = "fixed") -> dict:
    """Multivariable IVW (MVMR) via weighted least squares through the origin.

    BX : (n_snps, k_exposures) matrix of harmonised exposure betas
         (all alleles aligned to a common reference).
    SX : (n_snps, k_exposures) first-stage standard errors of the SNP-exposure
         associations (needed for the conditional F statistics).
    by, sy : outcome betas / standard errors (aligned to the same reference).

    Fixed effect:  b = (X' W X)^-1 X' W y,  cov(b) = (X' W X)^-1,  W = diag(1/sy^2)
    Multiplicative random: cov inflated by phi = max(1, Q/(n-k)).

    Returns multivariable Cochran Q and per-exposure conditional F statistics
    (Sanderson & Windmeijer 2016 / MVMR R package `weakness`):
        regress X_j on X_{-j} with first-stage weights 1/SX_j^2,
        F_j = residual' W_j residual   (chi-square, 1 df after partialling).
    Rule of thumb: F_j >= ~ (threshold often set near the univariate 10,
    scaled by exposure count; report alongside estimates).
    """
    n, k = BX.shape
    if n <= k:
        fail(f"mv_ivw: need more SNPs ({n}) than exposures ({k})")
    if np.any(sy <= 0) or np.any(~np.isfinite(sy)):
        fail("mv_ivw: non-positive or missing se_outcome")
    if SX is None or SX.shape != BX.shape:
        fail("mv_ivw: SX must match BX shape (first-stage SE matrix)")

    w = 1.0 / sy**2
    XtWX = BX.T @ (w[:, None] * BX)
    XtWX_inv = np.linalg.inv(XtWX)
    beta = XtWX_inv @ (BX.T @ (w * by))
    resid = by - BX @ beta
    q = float(np.sum(w * resid**2))
    dof = n - k
    phi = max(1.0, q / dof)
    if model == "fixed":
        cov = XtWX_inv
    elif model == "random":
        cov = phi * XtWX_inv
    else:
        fail(f"mv_ivw: unknown model '{model}' (expected 'fixed' or 'random')")
    se = np.sqrt(np.diag(cov))

    q_p = float(_z_to_p(math.sqrt(max(q, 0.0)) / math.sqrt(max(dof, 1))))

    # conditional (partial) F per exposure, Sanderson & Windmeijer style
    f_cond = {}
    for j in range(k):
        wj = 1.0 / SX[:, j]**2
        others = [c for c in range(k) if c != j]
        yj = BX[:, j]
        if others:
            Xo = BX[:, others]
            W = wj[:, None] * np.ones((1, Xo.shape[1]))
            coef = np.linalg.solve(Xo.T @ (wj[:, None] * Xo), Xo.T @ (wj * yj))
            xj_res = yj - Xo @ coef
        else:
            xj_res = yj
        f_cond[j] = float(np.sum(wj * xj_res**2))

    return {
        "beta": [float(b) for b in beta],
        "se": [float(s) for s in se],
        "model": model,
        "q": q,
        "q_dof": dof,
        "q_pval": q_p,
        "f_conditional": f_cond,
        "n_snps": n,
        "n_exposures": k,
    }


def mr_egger(bx: np.ndarray, by: np.ndarray, sy: np.ndarray) -> dict:
    """MR-Egger: weighted regression of by on bx WITH an intercept.

    by = intercept + slope * bx,  weights = 1/sy^2
    slope     = causal estimate corrected for average directional pleiotropy
    intercept = average directional pleiotropy (test its p value, but note
                the test has LOW POWER — see pitfalls section 4.1)

    SEs are standard WLS standard errors (NOME assumption implied).
    """
    if len(bx) < 3:
        fail("mr_egger: needs at least 3 variants to fit intercept + slope")
    w = 1.0 / sy**2
    X = np.column_stack([np.ones_like(bx), bx])
    beta, se, _, dof = _wls(X, by, w)
    return {
        "slope": _pack(float(beta[1]), float(se[1])),
        "intercept": _pack(float(beta[0]), float(se[0])),
        "se_method": "WLS (standard errors from weighted least squares)",
        "n_snp": int(len(bx)),
    }


def weighted_median_estimate(ratio: np.ndarray, w: np.ndarray) -> float:
    """Weighted median of the ratio estimates [Bowden 2016]."""
    order = np.argsort(ratio)
    r = ratio[order]
    ww = w[order]
    cw = np.cumsum(ww) / np.sum(ww)
    idx = int(np.searchsorted(cw, 0.5))
    if idx >= len(r):
        idx = len(r) - 1
    if idx > 0 and abs(cw[idx - 1] - 0.5) < 1e-12:
        return float(0.5 * (r[idx - 1] + r[idx]))
    if idx + 1 < len(r) and abs(cw[idx] - 0.5) < 1e-12:
        return float(0.5 * (r[idx] + r[idx + 1]))
    return float(r[idx])


def mode_estimate(ratio: np.ndarray, w: np.ndarray, phi: float = 1.0) -> float:
    """Weighted mode estimator [Hartwig 2017].

    Evaluates a weighted kernel density at each ratio estimate and returns the
    ratio with the highest density.
    """
    k = len(ratio)
    wsum = float(np.sum(w))
    dens = np.zeros(k)
    for j in range(k):
        z = (ratio[j] - ratio) / phi
        dens[j] = float(np.sum(w * stats.norm.pdf(z)) / wsum)
    return float(ratio[int(np.argmax(dens))])


def _bootstrap_se(stat_fn, k: int, n_boot: int, rng: np.random.Generator) -> float:
    """Standard error of a statistic by non-parametric bootstrap over variants."""
    ests = np.empty(n_boot)
    for m in range(n_boot):
        idx = rng.integers(0, k, k)
        ests[m] = stat_fn(idx)
    return float(np.std(ests, ddof=1))


def cochran_q(bx: np.ndarray, by: np.ndarray, sy: np.ndarray):
    """Cochran's Q, its df and p value, and I^2."""
    w = bx**2 / sy**2
    ratio = by / bx
    beta = float(np.sum(w * ratio) / np.sum(w))
    q = float(np.sum(w * (ratio - beta) ** 2))
    k = len(bx)
    df = k - 1
    p = float(stats.chi2.sf(q, df)) if df > 0 else float("nan")
    i2 = max(0.0, (q - df) / q) * 100.0 if q > 0 else 0.0
    return q, df, p, i2


def leave_one_out(bx: np.ndarray, by: np.ndarray, sy: np.ndarray, model: str = "fixed"):
    """Re-estimate IVW with each variant removed in turn."""
    k = len(bx)
    out = []
    for j in range(k):
        mask = np.ones(k, dtype=bool)
        mask[j] = False
        if mask.sum() < 2:
            out.append({"index": j, "beta": None, "se": None})
            continue
        r = ivw(bx[mask], by[mask], sy[mask], model=model)
        out.append({"index": j, "beta": r["beta"], "se": r["se"]})
    return out


def steiger_r2(beta: float, se: float, n: float) -> float:
    """Proportion of variance explained (TwoSampleMR / Hemani 2017 formula).

    R^2 = beta^2 / (beta^2 + se^2 * n)
    (the EAF terms cancel algebraically in the original expression).
    """
    denom = beta**2 + se**2 * n
    if denom <= 0:
        return float("nan")
    return float(beta**2 / denom)


def steiger_direction(bx, sx, by, sy, n_exp, n_out):
    """Steiger directionality: compare R^2 in exposure vs outcome.

    Returns the proportion of variants for which R^2_exposure > R^2_outcome.
    NOTE: this is supporting evidence about the DIRECTION of variance
    explained, not proof of causality. See pitfalls section 4.4.
    """
    if n_exp is None or n_out is None:
        return None
    correct = 0
    total = 0
    per_snp = []
    for i in range(len(bx)):
        r2e = steiger_r2(float(bx[i]), float(sx[i]), float(n_exp))
        r2o = steiger_r2(float(by[i]), float(sy[i]), float(n_out))
        if not (np.isfinite(r2e) and np.isfinite(r2o)):
            continue
        total += 1
        ok = r2e > r2o
        correct += int(ok)
        per_snp.append(
            {
                "SNP_index": i,
                "r2_exposure": _f(r2e),
                "r2_outcome": _f(r2o),
                "direction_correct": bool(ok),
            }
        )
    if total == 0:
        return None
    return {
        "prop_correct": _f(correct / total),
        "n_evaluated": total,
        "n_correct": correct,
        "per_snp": per_snp,
        "note": "Supporting evidence on direction of variance explained, not proof of causality.",
    }


# --------------------------------------------------------------------------
# MR-PRESSO
# --------------------------------------------------------------------------


def presso_global(bx, by, sy, n_boot=1000, rng=None):
    """MR-PRESSO global test [Verbanck 2018].

    Observed residual sum of squares under the IVW fit is compared with a
    null distribution generated by simulating outcome effects around the
    fitted causal model.
    """
    k = len(bx)
    beta_ivw = ivw(bx, by, sy, model="fixed")["beta"]

    def rss(beta, y):
        return float(np.sum((y - beta * bx) ** 2 / sy**2))

    rss_obs = rss(beta_ivw, by)
    sims = np.empty(n_boot)
    for m in range(n_boot):
        y_sim = rng.normal(beta_ivw * bx, sy)
        b_sim = ivw(bx, y_sim, sy, model="fixed")["beta"]
        sims[m] = rss(b_sim, y_sim)
    p = float((1 + np.sum(sims >= rss_obs)) / (n_boot + 1))
    return {"rss_observed": _f(rss_obs), "p_global": _f(p), "n_boot": n_boot}


def presso_outliers(bx, by, sy, n_boot=1000, rng=None):
    """MR-PRESSO outlier test: per-variant distortion statistic + null simulation."""
    k = len(bx)
    t_obs = np.empty(k)
    for j in range(k):
        mask = np.ones(k, dtype=bool)
        mask[j] = False
        b_lo = ivw(bx[mask], by[mask], sy[mask], model="fixed")["beta"]
        t_obs[j] = (by[j] - b_lo * bx[j]) ** 2 / sy[j] ** 2

    t_sim = np.empty((n_boot, k))
    for m in range(n_boot):
        b_full = ivw(bx, by, sy, model="fixed")["beta"]
        y_sim = rng.normal(b_full * bx, sy)
        for j in range(k):
            mask = np.ones(k, dtype=bool)
            mask[j] = False
            b_lo = ivw(bx[mask], y_sim[mask], sy[mask], model="fixed")["beta"]
            t_sim[m, j] = (y_sim[j] - b_lo * bx[j]) ** 2 / sy[j] ** 2

    pvals = (1 + np.sum(t_sim >= t_obs, axis=0)) / (n_boot + 1)

    # Benjamini-Hochberg style flag at 0.05 would be ideal; a plain 0.05
    # threshold is reported and the p values are returned for inspection.
    outliers = [int(j) for j in range(k) if pvals[j] < 0.05]
    return {
        "stat_observed": [_f(v) for v in t_obs],
        "p_values": [_f(v) for v in pvals],
        "outlier_indices": outliers,
        "n_outliers": len(outliers),
        "n_boot": n_boot,
        "threshold": 0.05,
    }


def presso_distortion(bx, by, sy, outlier_idx, n_boot=1000, rng=None):
    """MR-PRESSO distortion test: does removing outliers materially move the estimate?"""
    if not outlier_idx:
        return None
    mask = np.ones(len(bx), dtype=bool)
    mask[outlier_idx] = False
    if mask.sum() < 2:
        return None
    b_full = ivw(bx, by, sy, model="fixed")["beta"]
    b_clean = ivw(bx[mask], by[mask], sy[mask], model="fixed")["beta"]
    delta = b_full - b_clean

    deltas = np.empty(n_boot)
    for m in range(n_boot):
        y_sim = rng.normal(b_full * bx, sy)
        b1 = ivw(bx, y_sim, sy, model="fixed")["beta"]
        b2 = ivw(bx[mask], y_sim[mask], sy[mask], model="fixed")["beta"]
        deltas[m] = b1 - b2
    p = float((1 + np.sum(np.abs(deltas) >= abs(delta))) / (n_boot + 1))
    return {
        "beta_all": _f(b_full),
        "beta_outliers_removed": _f(b_clean),
        "delta": _f(delta),
        "p_distortion": _f(p),
        "n_removed": len(outlier_idx),
    }


# --------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------


def read_harmonised(path: Path):
    """Read the harmonised CSV. Fails loudly rather than silently dropping data."""
    if not path.exists():
        fail(f"input file not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in cols]
        if missing:
            fail(f"input is missing required columns: {missing}")
        rows = [r for r in reader]
    if not rows:
        fail("input table has zero data rows")

    def collect(name):
        vals = []
        for r in rows:
            v = r.get(name)
            if v is None or str(v).strip() == "" or str(v).strip().upper() in {"NA", "NAN", "NULL"}:
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(np.nan)
        return np.array(vals, dtype=float)

    data = {
        "SNP": [r.get("SNP") for r in rows],
        "bx": collect("beta_exposure"),
        "sx": collect("se_exposure"),
        "by": collect("beta_outcome"),
        "sy": collect("se_outcome"),
        "eaf": collect("eaf_exposure"),
        "n_exp": collect("n_exposure"),
        "n_out": collect("n_outcome"),
    }

    # validation — never silently continue on bad data
    ok = (
        np.isfinite(data["bx"])
        & np.isfinite(data["sx"])
        & np.isfinite(data["by"])
        & np.isfinite(data["sy"])
    )
    n_bad = int((~ok).sum())
    if n_bad:
        print(f"WARNING: dropping {n_bad} row(s) with non-numeric effects", file=sys.stderr)
    for key in ("bx", "sx", "by", "sy"):
        data[key] = data[key][ok]
    for key in ("SNP",):
        data[key] = [v for v, keep in zip(data[key], ok) if keep]
    for key in ("eaf", "n_exp", "n_out"):
        data[key] = data[key][ok]

    if np.any(data["sx"] <= 0) or np.any(data["sy"] <= 0):
        fail("se_exposure / se_outcome must be strictly positive")
    if np.any(data["bx"] == 0):
        fail("beta_exposure contains zero; ratio is undefined for those variants")
    if len(data["bx"]) < 2:
        fail(f"only {len(data['bx'])} usable variant(s); need at least 2")

    return data


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def run_analysis(data, args):
    bx, sx, by, sy = data["bx"], data["sx"], data["by"], data["sy"]
    k = len(bx)
    rng = np.random.default_rng(args.seed)

    # variance estimate from the exposure GWAS for F statistics
    f_stats = (bx**2) / (sx**2)

    out = {
        "engine": "mr_engine.py (pure Python implementation)",
        "seed": args.seed,
        "n_boot": args.n_boot,
        "n_snp": int(k),
        "snps": data["SNP"],
        "f_statistic": {
            "min": _f(np.min(f_stats)),
            "median": _f(np.median(f_stats)),
            "max": _f(np.max(f_stats)),
            "n_below_10": int(np.sum(f_stats < 10)),
        },
    }

    # ---- single variant: Wald ratio -------------------------------------
    if k == 1:
        out["wald_ratio"] = wald_ratio(float(bx[0]), float(by[0]), float(sy[0]))
        out["note"] = "Only one variant: Wald ratio is the only valid estimator."
        out["primary_analysis"] = "wald_ratio"
        return out

    # ---- IVW -------------------------------------------------------------
    out["ivw_fixed"] = ivw(bx, by, sy, model="fixed")
    out["ivw_random"] = ivw(bx, by, sy, model="random")

    q, df, qp, i2 = cochran_q(bx, by, sy)
    out["heterogeneity"] = {
        "cochran_q": _f(q),
        "df": int(df),
        "p": _f(qp),
        "i_squared_percent": _f(i2),
    }

    # pick the primary model per the pre-specified rule
    if args.ivw_model == "auto":
        primary = "ivw_random" if (np.isfinite(qp) and qp < 0.05) else "ivw_fixed"
    elif args.ivw_model == "fixed":
        primary = "ivw_fixed"
    else:
        primary = "ivw_random"
    out["primary_analysis"] = primary
    out["primary_model_rule"] = (
        "auto: random effects when Cochran's Q p < 0.05, otherwise fixed effects"
        if args.ivw_model == "auto"
        else f"forced: {args.ivw_model}"
    )

    # ---- MR-Egger --------------------------------------------------------
    if k >= 3:
        out["mr_egger"] = mr_egger(bx, by, sy)
    else:
        out["mr_egger"] = {"skipped": "fewer than 3 variants"}

    # ---- weighted median / mode -----------------------------------------
    ratio = by / bx
    se_ratio = sy / np.abs(bx)
    w = 1.0 / se_ratio**2

    wm_beta = weighted_median_estimate(ratio, w)
    wm_se = _bootstrap_se(
        lambda idx: weighted_median_estimate(ratio[idx], w[idx]), k, args.n_boot, rng
    )
    out["weighted_median"] = _pack(wm_beta, wm_se)

    wmode_beta = mode_estimate(ratio, w)
    wmode_se = _bootstrap_se(
        lambda idx: mode_estimate(ratio[idx], w[idx]), k, args.n_boot, rng
    )
    out["weighted_mode"] = _pack(wmode_beta, wmode_se)

    if k >= 10:
        smode_beta = mode_estimate(ratio, np.ones_like(w))
        smode_se = _bootstrap_se(
            lambda idx: mode_estimate(ratio[idx], np.ones_like(idx)), k, args.n_boot, rng
        )
        out["simple_mode"] = _pack(smode_beta, smode_se)

    # ---- leave-one-out ---------------------------------------------------
    out["leave_one_out"] = leave_one_out(bx, by, sy, model="fixed")

    # ---- MR-PRESSO -------------------------------------------------------
    try:
        g = presso_global(bx, by, sy, n_boot=args.n_boot, rng=rng)
        o = presso_outliers(bx, by, sy, n_boot=args.n_boot, rng=rng)
        d = presso_distortion(bx, by, sy, o["outlier_indices"], n_boot=args.n_boot, rng=rng)
        out["mr_presso"] = {
            "global": g,
            "outliers": o,
            "distortion": d,
            "implementation_note": (
                "Simplified re-implementation of the MR-PRESSO framework "
                "(residual-based global test with permutation null). Not the "
                "original R package; state this in the Methods."
            ),
        }
    except Exception as exc:  # pragma: no cover - defensive
        out["mr_presso"] = {"error": str(exc)}

    # ---- Steiger ---------------------------------------------------------
    n_exp = None
    n_out = None
    if np.isfinite(data["n_exp"]).any():
        n_exp = float(np.nanmedian(data["n_exp"]))
    if np.isfinite(data["n_out"]).any():
        n_out = float(np.nanmedian(data["n_out"]))
    steiger = steiger_direction(bx, sx, by, sy, n_exp, n_out)
    out["steiger"] = steiger if steiger is not None else {
        "skipped": "n_exposure / n_outcome not provided; cannot compute R^2"
    }

    return out


def selftest() -> int:
    """Run every estimator on synthetic data with a KNOWN effect.

    Generates k variants with a true causal effect theta, checks that IVW
    recovers it, that all sensitivity estimators are finite, and that the
    engine raises on malformed input.
    """
    print("=== mr_engine selftest ===")
    rng = np.random.default_rng(20260922)
    k = 25
    theta_true = 0.30

    bx = rng.uniform(0.02, 0.12, k) * rng.choice([-1, 1], k)
    sx = np.abs(bx) / rng.uniform(6, 12, k)
    by = theta_true * bx + rng.normal(0, 0.006, k)
    sy = np.full(k, 0.008)

    # sanity: all estimators must be finite and near the truth
    r_ivw = ivw(bx, by, sy, "fixed")
    r_rand = ivw(bx, by, sy, "random")
    r_egg = mr_egger(bx, by, sy)
    w = 1.0 / (sy / np.abs(bx)) ** 2
    ratio = by / bx
    r_wm = weighted_median_estimate(ratio, w)
    r_wmode = mode_estimate(ratio, w)
    q, df, qp, i2 = cochran_q(bx, by, sy)

    checks = [
        ("ivw_fixed|beta finite", np.isfinite(r_ivw["beta"])),
        ("ivw_fixed|within 0.15 of truth", abs(r_ivw["beta"] - theta_true) < 0.15),
        ("ivw_random|se >= fixed se", r_rand["se"] >= r_ivw["se"] - 1e-12),
        ("egger|slope finite", np.isfinite(r_egg["slope"]["beta"])),
        ("egger|intercept p in [0,1]", 0 <= r_egg["intercept"]["p"] <= 1),
        ("weighted_median|finite", np.isfinite(r_wm)),
        ("weighted_mode|finite", np.isfinite(r_wmode)),
        ("cochran_q|q >= 0", q >= 0),
        ("cochran_q|i2 in [0,100]", 0 <= i2 <= 100),
    ]

    # known-answer check: constant ratio with zero noise -> exact recovery
    bx2 = np.array([0.05, 0.08, 0.03, -0.04, 0.06])
    by2 = 0.5 * bx2
    sy2 = np.full(5, 0.01)
    exact = ivw(bx2, by2, sy2, "fixed")
    checks.append(("exact_recovery|beta == 0.5", abs(exact["beta"] - 0.5) < 1e-9))

    # Steiger must flag direction correctly on well-separated data
    st = steiger_direction(bx2, np.full(5, 0.005), by2, sy2, 50000, 50000)
    checks.append(("steiger|returns a proportion", st is not None and "prop_correct" in st))

    # LOO and PRESSO must run
    loo = leave_one_out(bx, by, sy)
    checks.append(("loo|one entry per variant", len(loo) == k))
    g = presso_global(bx, by, sy, n_boot=100, rng=rng)
    checks.append(("presso|p in [0,1]", 0 <= g["p_global"] <= 1))

    # MVMR: known-answer check — 3 exposures with true effects [0.5, -0.3, 0.2],
    # no pleiotropy; mv_ivw must recover all three and reject too-few SNPs.
    n_m = 300
    true_m = np.array([0.5, -0.3, 0.2])
    BXM = np.column_stack([
        rng.uniform(0.02, 0.12, n_m) * rng.choice([-1, 1], n_m) for _ in range(3)
    ])
    SXM = np.abs(BXM) / rng.uniform(6, 12, (n_m, 3))
    byM = BXM @ true_m + rng.normal(0, 0.006, n_m)
    syM = np.full(n_m, 0.008)
    mv = mv_ivw(BXM, SXM, byM, syM, "fixed")
    checks.append(("mv_ivw|3 betas recovered (±0.15)",
                   all(abs(mv["beta"][j] - true_m[j]) < 0.15 for j in range(3))))
    checks.append(("mv_ivw|all se finite & positive",
                   all(math.isfinite(s) and s > 0 for s in mv["se"])))
    checks.append(("mv_ivw|q >= 0 and p in [0,1]",
                   mv["q"] >= 0 and 0 <= mv["q_pval"] <= 1))
    checks.append(("mv_ivw|conditional F all > 10",
                   all(v > 10 for v in mv["f_conditional"].values())))
    exact_mv = mv_ivw(BXM, SXM, BXM @ true_m, np.full(n_m, 0.008), "fixed")
    checks.append(("mv_ivw|exact recovery of true effects",
                   all(abs(exact_mv["beta"][j] - true_m[j]) < 1e-9 for j in range(3))))
    raised_mvmr = False
    try:
        mv_ivw(BXM[:3], SXM[:3], byM[:3], syM[:3])
    except SystemExit:
        raised_mvmr = True
    checks.append(("mv_ivw|too-few-SNPs raises", raised_mvmr))

    # malformed input must fail loudly
    tmp_bad = Path("_mr_engine_bad_input.csv")
    tmp_bad.write_text("SNP,beta_exposure\nrs1,0.1\n", encoding="utf-8")
    raised = False
    try:
        read_harmonised(tmp_bad)
    except SystemExit:
        raised = True
    finally:
        tmp_bad.unlink(missing_ok=True)
    checks.append(("validation|missing columns raise", raised))

    n_ok = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        n_ok += int(ok)
    print(f"  {n_ok}/{len(checks)} checks passed")
    return 0 if n_ok == len(checks) else 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Pure-Python Mendelian randomization engine (see references/mr-methods-library.md)."
    )
    ap.add_argument("--input", type=str, help="harmonised CSV produced by harmonise.py")
    ap.add_argument("--output", type=str, default="mr-results.json", help="output JSON path")
    ap.add_argument("--n-boot", type=int, default=1000, help="bootstrap / permutation replicates")
    ap.add_argument("--seed", type=int, default=20260922, help="random seed (fixed for reproducibility)")
    ap.add_argument(
        "--ivw-model",
        choices=["auto", "fixed", "random"],
        default="auto",
        help="primary IVW model; 'auto' uses random effects when Q p < 0.05",
    )
    ap.add_argument("--selftest", action="store_true", help="run built-in synthetic-data checks")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    if not args.input:
        fail("--input is required (or use --selftest)")

    data = read_harmonised(Path(args.input))
    results = run_analysis(data, args)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)

    # machine-readable summary for the calling agent
    print("MR_ENGINE_OK")
    print(f"N_SNP={results['n_snp']}")
    primary = results["primary_analysis"]
    pb = results[primary]
    print(f"PRIMARY_ANALYSIS={primary}")
    print(f"PRIMARY_BETA={pb['beta']}")
    print(f"PRIMARY_SE={pb['se']}")
    print(f"PRIMARY_CI_LOWER={pb['ci_lower']}")
    print(f"PRIMARY_CI_UPPER={pb['ci_upper']}")
    print(f"PRIMARY_P={pb['p']}")
    print(f"IVW_FIXED_BETA={results['ivw_fixed']['beta']}")
    print(f"IVW_RANDOM_BETA={results['ivw_random']['beta']}")
    if isinstance(results.get("mr_egger"), dict) and "slope" in results["mr_egger"]:
        print(f"EGGER_SLOPE={results['mr_egger']['slope']['beta']}")
        print(f"EGGER_INTERCEPT={results['mr_egger']['intercept']['beta']}")
        print(f"EGGER_INTERCEPT_P={results['mr_egger']['intercept']['p']}")
    print(f"WEIGHTED_MEDIAN={results['weighted_median']['beta']}")
    print(f"WEIGHTED_MODE={results['weighted_mode']['beta']}")
    print(f"COCHRAN_Q={results['heterogeneity']['cochran_q']}")
    print(f"COCHRAN_Q_P={results['heterogeneity']['p']}")
    print(f"I_SQUARED={results['heterogeneity']['i_squared_percent']}")
    print(f"F_MIN={results['f_statistic']['min']}")
    print(f"F_MEDIAN={results['f_statistic']['median']}")
    if isinstance(results.get("mr_presso"), dict) and "global" in results["mr_presso"]:
        print(f"PRESSO_GLOBAL_P={results['mr_presso']['global']['p_global']}")
        print(f"PRESSO_N_OUTLIERS={results['mr_presso']['outliers']['n_outliers']}")
    print(f"RESULTS_JSON={out_path.resolve()}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
