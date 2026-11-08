#!/usr/bin/env python3
"""Audit locked numerical results against source outputs and manuscript text.

This script is deliberately read-only with respect to the analysis data.  It
checks the results that carry the main narrative, then writes a machine-readable
CSV and a concise Markdown report.  A non-zero exit status indicates version
drift that should be resolved before manuscript submission.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "outputs"
MANUSCRIPT = ROOT / "manuscript" / "ars" / "phase4_drafting" / "manuscript.md"
STORY = ROOT / "manuscript" / "00_story_outline_and_analysis_lock.md"
TOL = 5e-7


def row(df: pd.DataFrame, **selectors):
    mask = pd.Series(True, index=df.index)
    for column, value in selectors.items():
        mask &= df[column].eq(value)
    selected = df.loc[mask]
    if len(selected) != 1:
        raise AssertionError(f"Expected one row for {selectors}, found {len(selected)}")
    return selected.iloc[0]


checks: list[dict] = []


def numeric_check(section: str, item: str, actual, expected, source: str, tolerance=TOL):
    actual = float(actual)
    expected = float(expected)
    passed = math.isclose(actual, expected, rel_tol=0, abs_tol=tolerance)
    checks.append({
        "section": section,
        "item": item,
        "actual": f"{actual:.12g}",
        "expected": f"{expected:.12g}",
        "source": source,
        "status": "PASS" if passed else "FAIL",
    })


def text_check(section: str, item: str, pattern: str, documents: list[Path]):
    locations = []
    for document in documents:
        text = document.read_text(encoding="utf-8")
        if re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE):
            locations.append(str(document.relative_to(ROOT)))
    checks.append({
        "section": section,
        "item": item,
        "actual": "; ".join(locations) if locations else "not found",
        "expected": f"pattern: {pattern}",
        "source": "manuscript Markdown",
        "status": "PASS" if locations else "FAIL",
    })


# Diagnostic discrimination -------------------------------------------------
performance_path = ROOT / "outputs" / "biomarker_performance_results.csv"
performance = pd.read_csv(performance_path)
diagnostic_expectations = [
    ("A status", "MTBR-Tau243/Aβ42", 135, 0.9023266023, 0.845485, 0.950158),
    ("A status", "MTBR-Tau243", 135, 0.8985952590, 0.838271, 0.950287),
    ("A status, non-control cohort", "MTBR-Tau243/Aβ42", 101, 0.8937664618, 0.828203, 0.947193),
    ("T status, Tau PET subset", "MTBR-Tau243", 68, 0.9602076125, 0.908844, 0.995641),
    ("T status, Tau PET subset", "p-Tau217", 68, 0.9480968858, 0.887719, 0.989475),
    ("T status", "MTBR-Tau243", 114, 0.9430147059, 0.886790, 0.984192),
    ("T status", "p-Tau217", 114, 0.9341911765, 0.879983, 0.976615),
    ("T status, A+ subset", "MTBR-Tau243", 67, 0.8894830659, 0.803111, 0.959905),
    ("T status, A+ subset", "p-Tau217", 67, 0.8716577540, 0.775635, 0.946856),
]
for analysis, marker, n, auc, lower, upper in diagnostic_expectations:
    r = row(performance, Analysis=analysis, Biomarker=marker)
    label = f"{analysis}: {marker}"
    numeric_check("Diagnostic performance", f"{label}, N", r["N"], n, performance_path.name)
    numeric_check("Diagnostic performance", f"{label}, AUC", r["AUC"], auc, performance_path.name)
    numeric_check("Diagnostic performance", f"{label}, CI lower", r["AUC 95% CI lower"], lower, performance_path.name, 5e-6)
    numeric_check("Diagnostic performance", f"{label}, CI upper", r["AUC 95% CI upper"], upper, performance_path.name, 5e-6)

cv_path = OUTPUT / "diagnostic_cv_performance.csv"
cv = pd.read_csv(cv_path)
for analysis, marker, expected in [
    ("T status, Tau PET subset", "MTBR-Tau243", 0.9604353741),
    ("T status, Tau PET subset", "p-Tau217", 0.9491972789),
    ("T status", "MTBR-Tau243", 0.9437678571),
    ("T status", "p-Tau217", 0.9352678571),
    ("T status, A+ subset", "MTBR-Tau243", 0.8897142857),
]:
    r = row(cv, Analysis=analysis, Biomarker=marker)
    numeric_check("Cross-validation", f"{analysis}: {marker}, CV AUC", r["CV Auc mean"], expected, cv_path.name, 5e-6)

# A+ tau PET-only analysis (n=52, primary A+-conditioned T analysis) --------
# Now computed in the same pass as the other 5 analyses (see
# Fig3_ROC.py), so it is read from the same shared CSV
# with the "T status, A+ subset, Tau PET only" analysis key, not a separate file.
apet = performance.loc[performance["Analysis"].eq("T status, A+ subset, Tau PET only")]
for marker, auc, lower, upper in [
    ("MTBR-Tau243", 0.9411764706, 0.870258, 0.993468),
    ("p-Tau217", 0.9199346405, 0.833333, 0.984225),
]:
    r = row(apet, Biomarker=marker)
    numeric_check("A+ tau PET-only (n=52)", f"{marker}, N", r["N"], 52, performance_path.name)
    numeric_check("A+ tau PET-only (n=52)", f"{marker}, AUC", r["AUC"], auc, performance_path.name)
    numeric_check("A+ tau PET-only (n=52)", f"{marker}, CI lower", r["AUC 95% CI lower"], lower, performance_path.name, 5e-6)
    numeric_check("A+ tau PET-only (n=52)", f"{marker}, CI upper", r["AUC 95% CI upper"], upper, performance_path.name, 5e-6)

# Distributional heterogeneity ----------------------------------------------
gmm_path = OUTPUT / "bimodality_validation_summary.csv"
gmm = pd.read_csv(gmm_path)
r = row(gmm, Marker="MTBR-Tau243")
for column, expected in [
    ("N", 33), ("Observed delta BIC", 21.3784657462),
    ("Observed low n", 12), ("Observed high n", 21),
    ("Observed decision boundary", 2.2652404247),
    ("LOO minimum delta BIC", 19.5004502626),
    ("LOO maximum delta BIC", 24.8974148597),
]:
    numeric_check("Full A+T− GMM", column, r[column], expected, gmm_path.name)

reference_path = OUTPUT / "GMM_reference_stratified_summary.csv"
reference = pd.read_csv(reference_path)
for stratum, n, delta, p, q, loo_min, loo_max in [
    ("Tau PET-negative", 18, 16.1606407489, 0.0019980020, 0.0039960040, 14.0133507870, 19.9568415825),
    ("CSF-defined T-negative", 15, 5.5251427292, 0.0869130869, 0.0869130869, 3.5128263709, 8.4643973906),
]:
    r = row(reference, **{"Reference stratum": stratum})
    for column, expected in [
        ("N", n), ("Observed delta BIC", delta),
        ("Parametric-bootstrap LRT p", p), ("Bootstrap LRT BH-FDR q", q),
        ("LOO minimum delta BIC", loo_min), ("LOO maximum delta BIC", loo_max),
    ]:
        numeric_check("Reference-stratified GMM", f"{stratum}: {column}", r[column], expected, reference_path.name)

# Molecular profiles ---------------------------------------------------------
plasma_path = OUTPUT / "A+T-Subgroup_biomarkers.csv"
plasma = pd.read_csv(plasma_path)
for marker, effect, q in [
    ("p-Tau217", 1.0, 1.5862389855e-5),
    ("p-Tau231", 0.9285714286, 3.9036046736e-5),
    ("Aβ1-42", -0.6666666667, 0.0035572019),
]:
    r = row(plasma, Variable=marker)
    numeric_check("Plasma profile", f"{marker}: rank-biserial r", r["Rank-biserial r (positive means higher in High subgroup)"], effect, plasma_path.name)
    numeric_check("Plasma profile", f"{marker}: FDR q", r["FDR q (6 biomarkers)"], q, plasma_path.name)

csf_path = OUTPUT / "A+T-_CSF_Lumipulse_comparison.csv"
csf = pd.read_csv(csf_path)
for marker, effect, q in [
    ("CSF p-tau181/Aβ42", 0.875, 0.0216450216),
    ("CSF t-tau/Aβ42", 0.7708333333, 0.03496503497),
    ("CSF Aβ42/Aβ40", -0.625, 0.0715950716),
]:
    r = row(csf, Variable=marker)
    numeric_check("Exploratory CSF profile", f"{marker}: rank-biserial r", r["Rank-biserial r (positive means higher in High group)"], effect, csf_path.name)
    numeric_check("Exploratory CSF profile", f"{marker}: FDR q", r["FDR q (evaluable CSF comparisons)"], q, csf_path.name)

# Manuscript rendering checks: key values must be present in the locked prose.
docs = [MANUSCRIPT, STORY]
for item, pattern in [
    ("Primary tau PET AUCs", r"0\.960\s*\(95% CI,?\s*0\.909[–-]0\.996\).*0\.948\s*\(0\.888[–-]0\.989\)"),
    ("Operational PET/CSF AUCs", r"0\.943\s*[\(\[]95% CI,?\s*0\.887[–-]0\.984[\)\]].*0\.934\s*[\(\[](?:95% CI,?\s*)?0\.880[–-]0\.977[\)\]]"),
    ("Full GMM", r"ΔBIC\s*=\s*21\.38"),
    ("PET-stratified GMM", r"ΔBIC\s*=\s*16\.16.*p\s*=\s*0\.002.*q\s*=\s*0\.004"),
    ("CSF-stratified GMM", r"ΔBIC\s*=\s*5\.53.*p\s*=\s*q\s*=\s*0\.087"),
    ("Plasma gradient", r"p[-‑]tau217.*q\s*<\s*0\.001.*r\s*=\s*1\.00"),
    ("Exploratory CSF gradient", r"p[-‑]tau181/Aβ42.*q\s*=\s*0\.022.*t[-‑]tau/Aβ42.*q\s*=\s*0\.035"),
]:
    text_check("Manuscript text", item, pattern, docs)


audit = pd.DataFrame(checks)
csv_path = OUTPUT / "manuscript_numeric_integrity_audit.csv"
md_path = OUTPUT / "manuscript_numeric_integrity_audit.md"
audit.to_csv(csv_path, index=False)

failed = audit[audit["status"] == "FAIL"]
summary = [
    "# Manuscript numerical-integrity audit",
    "",
    f"- Checks run: **{len(audit)}**",
    f"- Passed: **{(audit['status'] == 'PASS').sum()}**",
    f"- Failed: **{len(failed)}**",
    f"- Overall status: **{'PASS' if failed.empty else 'FAIL'}**",
    "",
    "The audit covers the locked diagnostic-performance results, the A+ tau PET-only sensitivity analysis, full and reference-stratified GMM results, and the plasma/CSF molecular-profile results. It does not validate material details that remain author-confirmation items.",
    "",
]
if not failed.empty:
    summary.extend(["## Failed checks", "", failed.to_markdown(index=False), ""])
summary.extend(["## Complete audit", "", audit.to_markdown(index=False), ""])
md_path.write_text("\n".join(summary), encoding="utf-8")

print(f"Audit: {len(audit)} checks; {(audit['status'] == 'PASS').sum()} passed; {len(failed)} failed")
print(f"CSV: {csv_path}")
print(f"Markdown: {md_path}")
if not failed.empty:
    print(failed[["section", "item", "actual", "expected"]].to_string(index=False))
    raise SystemExit(1)

