#!/usr/bin/env python3
"""Supplementary Table S9: bimodality validation and Low/High subgroup comparisons
for eMTBR-tau243 within A+T− defined by tau PET or, if unavailable, CSF (n=33).

v3.62 (2026-11-11). Pure formatting of the locked outputs written by
SI_PETCSF_GMM_full_sample.py, SI_PETCSF_unimodality_tests.py,
SI_PETCSF_GMM_by_reference.py and SI_PETCSF_subgroup_comparison.py
(outputs/Aplus_Tminus_PET_CSF/). No statistic is computed here.
"""
import os
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(SCRIPT_DIR.parent / "outputs")))
D = ROOT / "Aplus_Tminus_PET_CSF"
OUT = ROOT / "TableS9_PETCSF_Aplus_Tminus.md"


def fp(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def rng(s):
    return f"{s.min():.2f} to {s.max():.2f}"


full = pd.read_csv(D / "bimodality_validation_summary.csv").set_index("Marker")
loo = pd.read_csv(D / "bimodality_leave_one_out.csv")
ext = pd.read_csv(D / "bimodality_extreme_value_sensitivity.csv")
uni = pd.read_csv(D / "unimodality_test_results.csv").set_index("Source column")
strat = pd.read_csv(D / "GMM_reference_stratified_summary.csv").set_index("Reference stratum")
sloo = pd.read_csv(D / "GMM_reference_stratified_leave_one_out.csv")
sext = pd.read_csv(D / "GMM_reference_stratified_extreme_sensitivity.csv")
bio = pd.read_csv(D / "A+T-Subgroup_biomarkers.csv").set_index("Variable")
cont = pd.read_csv(D / "A+T-Subgroup_continuous.csv").set_index("Variable")
cat = pd.read_csv(D / "A+T-Subgroup_categorical.csv").set_index("Variable")
csf = pd.read_csv(D / "A+T-_CSF_Lumipulse_comparison.csv").set_index("Variable")
# The tau PET-defined stratum p value is taken from the primary n=18 run (Supplementary Table S7,
# Supplementary Figure S6b), as in FigS6_PETCSF_Aplus_Tminus_bimodality.py.
primary18 = pd.read_csv(ROOT / "Aplus_Tminus_tauPET_only" / "bimodality_validation_summary.csv").set_index("Marker")

assert int(full.loc["MTBR-Tau243", "N"]) == 33

lines = [
    "# Supplementary Table S9. eMTBR-tau243 heterogeneity within A+T− defined by tau PET or CSF (n=33)",
    "",
    "**(a) Gaussian-mixture and unimodality tests**",
    "",
    "| Population | Marker | N | Low n | High n | ΔBIC (2 vs 1 component) | Parametric-bootstrap LRT p | Decision boundary (pg/mL) | Leave-one-out ΔBIC range | Extreme-value ΔBIC range | Hartigan dip p | Silverman p |",
    "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
]
for marker, label, ucol in [("MTBR-Tau243", "eMTBR-tau243", "MTBR-Tau243"), ("p-Tau217", "p-tau217 (comparison)", "p-Tau217(≤0.758)")]:
    r = full.loc[marker]
    lines.append(
        f"| A+T− (tau PET or CSF) | {label} | {int(r['N'])} | {int(r['Observed low n'])} | {int(r['Observed high n'])} | "
        f"{r['Observed delta BIC']:.2f} | {fp(r['Parametric-bootstrap LRT p under 1-component null'])} | "
        f"{r['Observed decision boundary']:.3f} | {rng(loo.loc[loo['Marker'].eq(marker), 'Delta BIC'])} | "
        f"{rng(ext.loc[ext['Marker'].eq(marker) & ~ext['Scenario'].eq('All observations'), 'Delta BIC'])} | "
        f"{fp(uni.loc[ucol, 'Hartigan dip bootstrap p'])} | {fp(uni.loc[ucol, 'Silverman bootstrap p'])} |"
    )
for stratum, label in [("Tau PET-negative", "Tau PET-defined A+T−"), ("CSF-defined T-negative", "CSF-defined A+T−")]:
    r = strat.loc[stratum]
    lines.append(
        f"| {label} | eMTBR-tau243 | {int(r['N'])} | {int(r['Observed low n'])} | {int(r['Observed high n'])} | "
        f"{r['Observed delta BIC']:.2f} | {fp(primary18.loc['MTBR-Tau243', 'Parametric-bootstrap LRT p under 1-component null'] if stratum == 'Tau PET-negative' else r['Parametric-bootstrap LRT p'])} | {r['Observed decision boundary']:.3f} | "
        f"{rng(sloo.loc[sloo['Reference stratum'].eq(stratum), 'Delta BIC'])} | "
        f"{rng(sext.loc[sext['Reference stratum'].eq(stratum) & ~sext['Scenario'].eq('All observations'), 'Delta BIC'])} | — | — |"
    )

lines += [
    "",
    "**(b) Low versus high eMTBR-tau243 groups (A+T− by tau PET or CSF)**",
    "",
    "| Family | Variable | Low n | Low median [IQR] | High n | High median [IQR] | Rank-biserial r | FDR q |",
    "| --- | --- | --- | --- | --- | --- | --- | --- |",
]
names = {"AGE": "Age, years", "MMSE_Total": "MMSE", "MOCA_P": "MoCA-P"}
for v in ["AGE", "MMSE_Total", "MOCA_P"]:
    r = cont.loc[v]
    lines.append(
        f"| Clinical | {names[v]} | {int(r['Low n'])} | {r['Low median [IQR]']} | {int(r['High n'])} | {r['High median [IQR]']} | "
        f"{r['Rank-biserial r (positive means higher in High subgroup)']:.2f} | {fp(r['FDR q (all clinical comparisons)'])} |"
    )


def pct(s):
    return s.replace("女", "Female").replace("男", "Male")


for v, label in [("SEX", "Sex, n (%)"), ("APOE ε4 carrier", "APOE ε4 carrier, n (%)")]:
    r = cat.loc[v]
    lines.append(
        f"| Clinical | {label} | {int(r['Low n available'])} | {pct(r['Low distribution'])} | {int(r['High n available'])} | "
        f"{pct(r['High distribution'])} | — | {fp(r['FDR q (all clinical comparisons)'])} |"
    )
bnames = {"p-Tau217": "p-tau217", "p-Tau231": "p-tau231", "Aβ1-42": "Aβ1-42", "BD Tau": "BD-tau", "GFAP": "GFAP", "NFL": "NfL"}
for v in ["p-Tau217", "p-Tau231", "Aβ1-42", "BD Tau", "GFAP", "NFL"]:
    r = bio.loc[v]
    lines.append(
        f"| Plasma (pg/mL) | {bnames[v]} | {int(r['Low n'])} | {r['Low median [IQR]']} | {int(r['High n'])} | {r['High median [IQR]']} | "
        f"{r['Rank-biserial r (positive means higher in High subgroup)']:.2f} | {fp(r['FDR q (6 biomarkers)'])} |"
    )
for v in csf.index:
    r = csf.loc[v]
    lines.append(
        f"| CSF (contemporaneous) | {v} | {int(r['Low n'])} | {r['Low median [IQR]']} | {int(r['High n'])} | {r['High median [IQR]']} | "
        f"{r['Rank-biserial r (positive means higher in High group)']:.2f} | {fp(r['FDR q (evaluable CSF comparisons)'])} |"
    )

lines += [
    "",
    "A+T− includes the 18 tau PET-defined and 15 CSF-defined A+T− participants. ΔBIC = BIC for the one-component model minus BIC for the "
    "two-component model; positive values favor two components. Parametric-bootstrap p values use 1,000 replicates under the one-component "
    "null and are not corrected for multiplicity; the tau PET-defined row repeats Supplementary Table S7. Leave-one-out and extreme-value ranges describe the influence of individual or extreme "
    "observations, not confidence intervals. Hartigan dip and Silverman tests were applied to log10 concentrations. In (b), groups are the "
    "two-component assignments in the combined population; eMTBR-tau243 defines the groups and is not compared. The clinical (5), plasma (6), "
    "and CSF (5) comparisons are separate Benjamini–Hochberg FDR families; Mann–Whitney U tests (exact permutation for CSF) and Fisher's exact "
    "tests. CSF comparisons are limited to the 14 participants with CSF collected within 30 days of blood collection (6 low, 8 high). This combined-reference analysis is shown for "
    "comparison with the tau PET-defined analysis (Supplementary Tables S7 and S8) and is descriptive.",
]
OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(OUT.read_text(encoding="utf-8"))
