#!/usr/bin/env python3
"""Incremental value of plasma MTBR-tau243 over p-tau217 for T status.

Added 2026-11-09 in response to simulated peer review (Revision Roadmap REV-2).

For each T-status population (tau PET n=68 and A+ tau PET n=52 primary;
combined PET/CSF n=114 and A+ PET/CSF n=67 secondary) this script fits
logistic regressions on log10 concentrations:

    M1: T ~ p-tau217
    M2: T ~ MTBR-tau243
    M3: T ~ p-tau217 + MTBR-tau243

and reports
  * likelihood-ratio tests M1 vs M3 (does MTBR-tau243 add to p-tau217?) and
    M2 vs M3 (does p-tau217 add to MTBR-tau243?);
  * apparent AUCs, and the paired patient-level bootstrap 95% CI of
    AUC(M3) - AUC(M1) with models refitted in every resample;
  * repeated stratified 5-fold cross-validated AUCs (50 repeats), so that the
    two-predictor model is not credited with in-sample optimism.

Models are unpenalized maximum likelihood (statsmodels Logit). If a fit does
not converge or shows (quasi-)complete separation, the row is flagged and
the LRT is not reported. Read-only with respect to the raw data.
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(PROJECT_DIR / "biomarker_result.xlsx")))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_VALIDATION_BOOTSTRAP", "5000"))
CV_REPEATS, CV_FOLDS = 50, 5

PTAU = "p-Tau217(≤0.758)"
MTBR = "MTBR-Tau243"


def clean_status(series: pd.Series) -> pd.Series:
    return (series.astype("string").str.strip()
            .str.replace("−", "-", regex=False).str.replace("–", "-", regex=False))


def populations(data: pd.DataFrame) -> list[tuple[str, str, pd.Series]]:
    a, t = clean_status(data["A"]), clean_status(data["T"])
    pet = pd.to_numeric(data["taupet"], errors="coerce").eq(1)
    has_t = t.isin(["T+", "T-"])
    return [
        ("T status, tau PET (n=68; primary)", "primary", has_t & pet),
        ("T status within A+, tau PET (n=52; primary)", "primary", has_t & pet & a.eq("A+")),
        ("T status, tau PET/CSF (n=114; secondary)", "secondary", has_t),
        ("T status within A+, tau PET/CSF (n=67; secondary)", "secondary", has_t & a.eq("A+")),
    ]


def fit_logit(y: np.ndarray, x: np.ndarray):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = sm.Logit(y, sm.add_constant(x)).fit(disp=0, maxiter=200)
            ok = bool(res.mle_retvals.get("converged", False))
            p = res.predict()
            if np.any(p < 1e-8) and np.any(p > 1 - 1e-8):
                ok = False  # (quasi-)separation
            return res, ok
        except Exception:  # perfect separation raises in statsmodels
            return None, False


def lin_pred(x: np.ndarray, y: np.ndarray, x_eval: np.ndarray) -> np.ndarray:
    """Near-unpenalized logistic fit used inside bootstrap/CV (robust to separation)."""
    m = LogisticRegression(C=1e6, max_iter=5000).fit(x, y)
    return m.decision_function(x_eval)


def analyse(name: str, tier: str, y: np.ndarray, X: pd.DataFrame, pop_index: int) -> dict:
    xp = X[["ptau"]].to_numpy()
    xm = X[["mtbr"]].to_numpy()
    xb = X[["ptau", "mtbr"]].to_numpy()
    r1, ok1 = fit_logit(y, xp)
    r2, ok2 = fit_logit(y, xm)
    r3, ok3 = fit_logit(y, xb)
    row = {"analysis": name, "tier": tier, "n": int(len(y)), "T+": int(y.sum()), "T-": int(len(y) - y.sum()),
           "fits_ok": bool(ok1 and ok2 and ok3)}
    if row["fits_ok"]:
        lr_add_mtbr = 2 * (r3.llf - r1.llf)
        lr_add_ptau = 2 * (r3.llf - r2.llf)
        row.update({
            "LRT_add_MTBR_to_ptau217_chi2": lr_add_mtbr, "LRT_add_MTBR_to_ptau217_p": chi2.sf(lr_add_mtbr, 1),
            "LRT_add_ptau217_to_MTBR_chi2": lr_add_ptau, "LRT_add_ptau217_to_MTBR_p": chi2.sf(lr_add_ptau, 1),
            "OR_per_log10_ptau217_in_M3": float(np.exp(r3.params[1])),
            "OR_per_log10_MTBR_in_M3": float(np.exp(r3.params[2])),
        })
    auc1 = roc_auc_score(y, lin_pred(xp, y, xp))
    auc2 = roc_auc_score(y, lin_pred(xm, y, xm))
    auc3 = roc_auc_score(y, lin_pred(xb, y, xb))
    row.update({"AUC_ptau217": auc1, "AUC_MTBR": auc2, "AUC_combined": auc3,
                "dAUC_combined_minus_ptau217": auc3 - auc1, "dAUC_combined_minus_MTBR": auc3 - auc2})

    rng = np.random.default_rng(RANDOM_SEED + 4_000_000 + pop_index)
    d1, d2 = [], []
    n = len(y)
    for _ in range(BOOTSTRAP_REPS):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.min() == yb.max():
            continue
        a1 = roc_auc_score(yb, lin_pred(xp[idx], yb, xp[idx]))
        a2 = roc_auc_score(yb, lin_pred(xm[idx], yb, xm[idx]))
        a3 = roc_auc_score(yb, lin_pred(xb[idx], yb, xb[idx]))
        d1.append(a3 - a1)
        d2.append(a3 - a2)
    row["dAUC_combined_minus_ptau217_CI"] = list(np.percentile(d1, [2.5, 97.5]))
    row["dAUC_combined_minus_MTBR_CI"] = list(np.percentile(d2, [2.5, 97.5]))
    row["bootstrap_reps_used"] = len(d1)

    cv = RepeatedStratifiedKFold(n_splits=CV_FOLDS, n_repeats=CV_REPEATS, random_state=RANDOM_SEED + pop_index)
    folds = {"ptau217": [], "MTBR": [], "combined": []}
    pooled = {k: np.zeros((CV_REPEATS, n)) for k in folds}
    for i, (tr, te) in enumerate(cv.split(xb, y)):
        rep = i // CV_FOLDS
        for key, xx in (("ptau217", xp), ("MTBR", xm), ("combined", xb)):
            pooled[key][rep, te] = lin_pred(xx[tr], y[tr], xx[te])
    for key in folds:
        aucs = [roc_auc_score(y, pooled[key][r]) for r in range(CV_REPEATS)]
        row[f"CV_AUC_{key}_mean"] = float(np.mean(aucs))
        row[f"CV_AUC_{key}_range"] = [float(np.min(aucs)), float(np.max(aucs))]
    return row


def main() -> None:
    data = pd.read_excel(INPUT, sheet_name="Sheet1")
    data[PTAU] = pd.to_numeric(data[PTAU], errors="coerce")
    data[MTBR] = pd.to_numeric(data[MTBR], errors="coerce")
    t = clean_status(data["T"])
    rows = []
    for i, (name, tier, mask) in enumerate(populations(data)):
        sub = data.loc[mask & data[PTAU].gt(0) & data[MTBR].gt(0)]
        y = clean_status(sub["T"]).eq("T+").astype(int).to_numpy()
        X = pd.DataFrame({"ptau": np.log10(sub[PTAU].to_numpy()), "mtbr": np.log10(sub[MTBR].to_numpy())})
        rows.append(analyse(name, tier, y, X, i))
    expected = [(68, 34), (52, 34), (114, 34), (67, 34)]
    got = [(r["n"], r["T+"]) for r in rows]
    if got != expected:
        raise AssertionError(f"Unexpected population sizes {got}, expected {expected}")
    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT_DIR / "incremental_value_ptau217_MTBR.csv", index=False)
    (OUTPUT_DIR / "incremental_value_ptau217_MTBR.json").write_text(
        json.dumps({"seed": RANDOM_SEED, "bootstrap_reps": BOOTSTRAP_REPS, "cv": f"{CV_REPEATS}x{CV_FOLDS}",
                    "rows": rows}, indent=2, default=float), encoding="utf-8")

    def f(x, d=3):
        return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"

    lines = ["# Incremental value of MTBR-tau243 over p-tau217 (T status)", "",
             f"Logistic models on log10 concentrations; paired bootstrap {BOOTSTRAP_REPS} resamples with refitting; "
             f"{CV_REPEATS}x{CV_FOLDS}-fold repeated CV. Seed {RANDOM_SEED}.", "",
             "| Population | n (T+/T−) | AUC p-tau217 | AUC MTBR | AUC combined | ΔAUC combined − p-tau217 (95% CI) | LRT add MTBR, p | LRT add p-tau217, p | CV AUC p-tau217 / MTBR / combined |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lo, hi = r["dAUC_combined_minus_ptau217_CI"]
        lrt_m = f"χ²={f(r.get('LRT_add_MTBR_to_ptau217_chi2'), 2)}, p={f(r.get('LRT_add_MTBR_to_ptau217_p'))}" if r["fits_ok"] else "not estimable (separation)"
        lrt_p = f"χ²={f(r.get('LRT_add_ptau217_to_MTBR_chi2'), 2)}, p={f(r.get('LRT_add_ptau217_to_MTBR_p'))}" if r["fits_ok"] else "not estimable (separation)"
        lines.append(f"| {r['analysis']} | {r['n']} ({r['T+']}/{r['T-']}) | {f(r['AUC_ptau217'])} | {f(r['AUC_MTBR'])} | "
                     f"{f(r['AUC_combined'])} | {f(r['dAUC_combined_minus_ptau217'])} ({f(lo)} to {f(hi)}) | {lrt_m} | {lrt_p} | "
                     f"{f(r['CV_AUC_ptau217_mean'])} / {f(r['CV_AUC_MTBR_mean'])} / {f(r['CV_AUC_combined_mean'])} |")
    (OUTPUT_DIR / "incremental_value_ptau217_MTBR.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
