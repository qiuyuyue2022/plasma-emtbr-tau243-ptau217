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
MANUSCRIPT = ROOT / "manuscript" / "full_manuscript_v1.md"
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
# Fig3_FigS2_TableS4_ROC.py), so it is read from the same shared CSV
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
    numeric_check("Pooled A+T− GMM (sensitivity)", column, r[column], expected, gmm_path.name)

# 2026-11-09 revision (Stage 4, REV-15): primary A+T- analysis = tau PET-defined T- (n=18),
# computed by 02a/02c/03 with BIOMARKER_APTM_REFERENCE=tauPET into outputs/Aplus_Tminus_tauPET_only/.
APTM = OUTPUT / "Aplus_Tminus_tauPET_only"
pet_gmm = pd.read_csv(APTM / "bimodality_validation_summary.csv")
r = row(pet_gmm, Marker="MTBR-Tau243")
for column, expected in [
    ("N", 18), ("Observed delta BIC", 16.1606407489),
    ("Parametric-bootstrap LRT p under 1-component null", 0.0029970030),
    ("Observed low n", 9), ("Observed high n", 9),
    ("Observed decision boundary", 1.9680356933),
    ("LOO minimum delta BIC", 14.0133507870), ("LOO maximum delta BIC", 19.9568415825),
]:
    numeric_check("Tau PET-defined A+T− GMM (primary)", column, r[column], expected, "Aplus_Tminus_tauPET_only/bimodality_validation_summary.csv")
r = row(pet_gmm, Marker="p-Tau217")
numeric_check("Tau PET-defined A+T− GMM (primary)", "p-tau217 delta BIC", r["Observed delta BIC"], -1.8950166066, "Aplus_Tminus_tauPET_only/bimodality_validation_summary.csv")
numeric_check("Tau PET-defined A+T− GMM (primary)", "p-tau217 bootstrap p", r["Parametric-bootstrap LRT p under 1-component null"], 0.2787212787, "Aplus_Tminus_tauPET_only/bimodality_validation_summary.csv")

uni = pd.read_csv(APTM / "unimodality_test_results.csv")
for marker, dip_p, silv_p in [("eMTBR-tau243", 0.00199, 0.0019996001), ("p-tau217", 0.16966, 0.0367926415)]:
    r = row(uni, Marker=marker)
    numeric_check("Unimodality tests (n=18)", f"{marker}: dip p", r["Hartigan dip bootstrap p"], dip_p, "Aplus_Tminus_tauPET_only/unimodality_test_results.csv", 5e-6)
    numeric_check("Unimodality tests (n=18)", f"{marker}: Silverman p", r["Silverman bootstrap p"], silv_p, "Aplus_Tminus_tauPET_only/unimodality_test_results.csv", 5e-6)

reference_path = OUTPUT / "GMM_reference_stratified_summary.csv"
reference = pd.read_csv(reference_path)
r = row(reference, **{"Reference stratum": "CSF-defined T-negative"})
for column, expected in [("N", 15), ("Observed delta BIC", 5.5251427292), ("Parametric-bootstrap LRT p", 0.0869130869)]:
    numeric_check("CSF-defined A+T− GMM (descriptive)", column, r[column], expected, reference_path.name)

plasma_path = APTM / "A+T-Subgroup_biomarkers.csv"
plasma = pd.read_csv(plasma_path)
for marker, effect, q in [
    ("p-Tau217", 1.0, 0.0012368844), ("p-Tau231", 1.0, 0.0012368844),
    ("Aβ1-42", -0.5555555556, 0.0955346627), ("GFAP", 0.5308641975, 0.0955346627),
]:
    r = row(plasma, Variable=marker)
    numeric_check("Plasma profile (n=18)", f"{marker}: rank-biserial r", r["Rank-biserial r (positive means higher in High subgroup)"], effect, "Aplus_Tminus_tauPET_only/A+T-Subgroup_biomarkers.csv")
    numeric_check("Plasma profile (n=18)", f"{marker}: FDR q", r["FDR q (6 biomarkers)"], q, "Aplus_Tminus_tauPET_only/A+T-Subgroup_biomarkers.csv")

inc_path = OUTPUT / "incremental_value_ptau217_MTBR.csv"
inc = pd.read_csv(inc_path)
for analysis, p_add_mtbr, p_add_ptau, cv_both in [
    ("T status, tau PET (n=68; primary)", 0.497903, 0.126054, 0.929135),
    ("T status within A+, tau PET (n=52; primary)", 0.729118, 0.170016, 0.888137),
    ("T status, tau PET/CSF (n=114; secondary)", 0.295687, 0.080336, 0.924971),
    ("T status within A+, tau PET/CSF (n=67; secondary)", 0.674911, 0.068093, 0.857433),
]:
    r = row(inc, analysis=analysis)
    numeric_check("Incremental value", f"{analysis}: LRT p adding MTBR", r["LRT_add_MTBR_to_ptau217_p"], p_add_mtbr, inc_path.name, 5e-6)
    numeric_check("Incremental value", f"{analysis}: LRT p adding p-tau217", r["LRT_add_ptau217_to_MTBR_p"], p_add_ptau, inc_path.name, 5e-6)
    numeric_check("Incremental value", f"{analysis}: CV AUC both", r["CV_AUC_combined_mean"], cv_both, inc_path.name, 5e-6)

apet_cv = row(performance.loc[performance["Analysis"].eq("T status, A+ subset, Tau PET only")], Biomarker="MTBR-Tau243")
numeric_check("A+ tau PET-only (n=52)", "MTBR-Tau243, CV accuracy (Results 3.4: 84.7%)", round(100 * apet_cv["CV accuracy"], 1), 84.7, performance_path.name, 1e-9)

# Two-step A/T comparison (tau PET-defined groups; code/Fig5_AT_two_step_comparison.py), 2026-11-09
two_path = OUTPUT / "AT_two_step_comparison.csv"
two = pd.read_csv(two_path)
for marker, r_a, q_a, r_t, q_t, auc in [
    ("eMTBR-tau243", 0.5729166667, None, 0.8823529412, None, 0.7864583333),
    ("p-tau217", 0.3333333333, None, 0.8398692810, None, 0.6666666667),
    ("GFAP", 0.5972222222, None, 0.0816993464, None, 0.7986111111),
    ("BD-tau", 0.1944444444, None, 0.7549019608, None, 0.5972222222),
]:
    r = row(two, Marker=marker)
    numeric_check("Two-step A/T", f"{marker}: amyloid-step r", r["Amyloid step r (A+T− vs A−T−)"], r_a, two_path.name, 5e-3)
    numeric_check("Two-step A/T", f"{marker}: tau-step r", r["Tau step r (A+T+ vs A+T−)"], r_t, two_path.name, 5e-3)
    numeric_check("Two-step A/T", f"{marker}: T− A-status AUC", r["T− A-status AUC"], auc, two_path.name, 5e-4)

# Age-adjusted amyloid step (code/Text_age_adjusted_amyloid_step.py), 2026-11-09
age_path = OUTPUT / "age_adjusted_amyloid_step.csv"
age = pd.read_csv(age_path)
for marker, lrt in [("eMTBR-tau243", 0.039), ("GFAP", 0.134), ("p-tau217", 0.141)]:
    r = row(age, Marker=marker)
    numeric_check("Age-adjusted amyloid step", f"{marker}: LRT p over age", round(r["LRT p over age"], 3), lrt, age_path.name, 1e-9)

# Manuscript rendering checks: key values must be present in the manuscript prose itself
# (2026-11-09: no longer also searched in the outline, which had let removed values pass).
docs = [MANUSCRIPT]
for item, pattern in [
    ("Primary tau PET AUCs", r"0\.960 \(95% CI, 0\.909[–-]0\.996.*p-tau217 \(AUC 0\.948 \[0\.888[–-]0\.989\]"),
    ("Operational PET/CSF AUCs", r"0\.943 \[95% CI, 0\.887[–-]0\.984\].*0\.934 \[0\.880[–-]0\.977\]"),
    ("A+ tau PET CV accuracy", r"Cross-validated accuracy was 84\.7%"),
    ("Tau PET-defined GMM", r"ΔBIC=16\.16; parametric-bootstrap p=0\.003"),
    ("Pooled GMM", r"ΔBIC=21\.38; p<0\.001"),
    ("Unimodality n=18", r"dip and Silverman tests both rejected unimodality \(p=0\.002 for each\).*dip test did not reject unimodality \(p=0\.170\), whereas the Silverman test did \(p=0\.037\)"),
    ("Plasma profile n=18", r"p-tau217 \(median 0\.976.*q=0\.001; rank-biserial r=1\.00\) and p-tau231.*q=0\.001; r=1\.00\)"),
    ("Incremental value", r"adding eMTBR-tau243 to p-tau217: likelihood-ratio p=0\.50 and 0\.73"),
    ("Two-step results", r"eMTBR-tau243 changed at both steps.*amyloid step r=0\.57, q=0\.016; tau step r=0\.88, q<0\.001"),
    ("T− A-status AUC", r"AUC of 0\.786 \(95% CI, 0\.619–0\.925\)"),
    ("Age-adjusted", r"eMTBR-tau243 remained associated with A status \(likelihood-ratio p=0\.039 over age alone\), whereas GFAP \(p=0\.134\) and p-tau217 \(p=0\.141\) did not"),
    ("A−T− tau PET n", r"A−T− \(n=16\), A\+T− \(n=18\), and A\+T\+ \(n=34\) participants, all classified by tau PET"),
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

