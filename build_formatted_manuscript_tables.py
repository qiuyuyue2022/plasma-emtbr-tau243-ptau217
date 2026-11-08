#!/usr/bin/env python3
"""Build reviewable, formatted manuscript-table drafts from locked output CSVs.

This script never recomputes any statistic. It only reformats existing,
already-locked CSV/MD outputs from the upstream analysis scripts into
submission-style Markdown tables. It must be run after:
  - Table1.py                                       (Table 1 draft, Supplementary Table S2)
  - Fig3_ROC.py                                      (biomarker_performance_results.csv)
  - 01_diagnostic_validation.py      (diagnostic_* CSVs)
  - 02a_bimodality_GMM_full_sample.py                 (bimodality_* CSVs)
  - 02b_bimodality_GMM_by_reference.py    (GMM_reference_stratified_* CSVs)
  - 03_Aplus_Tminus_subgroup_comparison.py (A+T-Subgroup_* CSVs; also writes
    TableS6_Aplus_Tminus_CSF_comparison.md directly, not via this script)

Renamed 2026-11-01 (author direction) from Tables1to4_build_revised_results_tables.py:
the old name and its build_table2()/build_table3() functions described a table
design (main-text 4-table package: participant characteristics, full 11-marker
performance grid, prespecified paired-AUC contrasts, A+T- low/high groups) that
has since been fully superseded -- see 00_story_outline_and_analysis_lock.md
Section 10 and changelog entries 3.4-3.6 for the full history. The dead
build_table2()/build_table3() functions have been removed; the best-marker-vs-
others Supplementary Table S1 (TableS1_best_marker_vs_others.md, moved from main
Table 2 on 2026-11-01) remains hand-maintained from the now-deleted
TableS7_best_vs_others_AUC_comparisons.csv, as before.

Renumbered again 2026-11-01 (author direction) to match true first-citation
order in methods_and_results_draft.md, after the previously-missing citations
for the full 11-marker grid, the biomarker-concentration table, and the GMM
sensitivity table were added to the Results text: old S6 (best marker vs
others) -> S1, old S1 (full 11-marker grid) -> S3, old S3 (A+T- low/high
groups) -> S5, old S5 (CSF comparison) -> S6; S2 and S4 were already cited in
the right relative position and are unchanged. Old Supplementary Table S7
(`TableS7_best_vs_others_AUC_comparisons.csv`) was deleted outright rather than
renumbered: it was found to be byte-for-byte the same 60 comparisons as the
best-marker-vs-others table (new S1), just as an unformatted CSV instead of a
submission-style Markdown table, so it carried no information S1 did not
already have. See changelog entry 3.10 for the full rationale and the complete
list of files touched by this renumbering.

Current outputs of this script:
  - Table1_participant_characteristics.md  (main Table 1, from Table1.py's draft)
  - TableS3_biomarker_performance_full.md  (Supplementary Table S3: full 11-marker
    x 6-population performance grid, replacing the orphaned biomarker_performance_tables.md
    whose generating script, 03_export_copyable_tables.py, no longer exists)
  - TableS5_Aplus_Tminus_groups.md         (Supplementary Table S5)
  - TableS4_GMM_sensitivity_analyses.md    (Supplementary Table S4: full-sample and
    reference-stratified GMM diagnostics, bootstrap, leave-one-out, extreme-value,
    and subgroup-membership sensitivity -- previously scattered across several
    unformatted analysis-passport files with no single manuscript-style table)
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def fmt(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def qfmt(value: float) -> str:
    return "<0.001" if value < 0.001 else fmt(value)


def pfmt(value: float) -> str:
    return "<0.001" if value < 0.001 else fmt(value)


def table(headers: list[str], rows: list[list[str]]) -> str:
    return "| " + " | ".join(headers) + " |\n| " + " | ".join(["---"] * len(headers)) + " |\n" + "\n".join(
        "| " + " | ".join(str(x) for x in row) + " |" for row in rows
    )


def build_table1() -> str:
    original = (OUT / "Table1_participant_characteristics_draft.md").read_text(encoding="utf-8")
    first = original.split("# Supplementary table.")[0]
    first = first.replace("A−T0", "A−, T reference unavailable")
    first = first.replace("  Not available (T0)", "  T reference unavailable")
    first = first.replace("T0 indicates that no tau PET or CSF T-reference result was available.",
                          "T reference unavailable indicates that no tau PET or eligible CSF T classification was available; it is not an AT(N) category.")
    if "All T+ participants had tau PET" not in first:
        first += "\n- All T+ participants had tau PET as the T reference; CSF contributed T− classifications only.\n"
    return first


def build_tableS3() -> str:
    d = pd.read_csv(OUT / "biomarker_performance_results.csv").set_index(["Analysis", "Biomarker"])
    analyses = [
        ("A status, non-control cohort", "A status, healthy controls excluded (n=101; primary)"),
        ("A status", "A status, all participants (n=135)"),
        ("T status, Tau PET subset", "T status, tau PET (n=68; primary)"),
        ("T status", "T status, tau PET/CSF (n=114)"),
        ("T status, A+ subset, Tau PET only", "T status within A+, tau PET (n=52; primary)"),
        ("T status, A+ subset", "T status within A+, tau PET/CSF (n=67)"),
    ]
    markers = [
        ("MTBR-Tau243", "eMTBR-tau243"),
        ("p-Tau217", "p-tau217"),
        ("p-Tau231", "p-tau231"),
        ("BD Tau", "BD-tau"),
        ("GFAP", "GFAP"),
        ("NFL", "NfL"),
        ("Aβ1-42", "Aβ1-42"),
        ("MTBR-Tau243/Aβ42", "eMTBR-tau243/Aβ42"),
        ("p-Tau217/Aβ42", "p-tau217/Aβ42"),
        ("p-Tau231/Aβ42", "p-tau231/Aβ42"),
        ("BD-Tau/Aβ42", "BD-tau/Aβ42"),
    ]
    rows = []
    for key, label in analyses:
        for marker, mlabel in markers:
            x = d.loc[(key, marker)]
            rows.append([
                label,
                mlabel,
                str(int(x["N"])),
                f"{fmt(x['AUC'])} ({fmt(x['AUC 95% CI lower'])}–{fmt(x['AUC 95% CI upper'])})",
                f"{x['Cutoff direction']}{fmt(x['Cutoff'], 2)}",
                f"{100*x['CV accuracy']:.1f}%",
                f"{100*x['CV balanced accuracy']:.1f}%",
            ])
    return (
        "# Supplementary Table S4. Full 11-marker diagnostic performance across six analysis populations\n\n"
        + table(
            ["Analysis population", "Biomarker", "n", "Apparent AUC (95% bootstrap CI)",
             "Sample-derived Youden cutoff", "Repeated-CV accuracy", "Repeated-CV balanced accuracy"],
            rows,
        )
        + (
            "\n\nApparent AUC 95% CIs use 3,000 nonparametric participant-level bootstrap samples. "
            "Cutoffs maximize the Youden index in the full analysis sample and are sample-derived and exploratory, "
            "not validated clinical thresholds. Repeated-CV metrics use stratified 5-fold cross-validation repeated 50 times "
            "and describe split-to-split variability, not independent-sample confidence intervals. "
            "The main text reports only eMTBR-tau243 and p-tau217 (the two leading markers) and the complementary "
            "best-marker-versus-other-ten comparison (Supplementary Table S2); this table is the only place the full 11-marker grid, "
            "including cross-validated accuracy and balanced accuracy for every marker and population, is shown.\n"
        )
    )


def build_tableS5() -> str:
    clinical = pd.read_csv(OUT / "A+T-Subgroup_continuous.csv")
    categorical = pd.read_csv(OUT / "A+T-Subgroup_categorical.csv")
    biomarkers = pd.read_csv(OUT / "A+T-Subgroup_biomarkers.csv")
    rows = []
    for _, x in clinical.iterrows():
        label = {"AGE": "Age, years", "MMSE_Total": "MMSE", "MOCA_P": "MoCA-P"}.get(x["Variable"], x["Variable"])
        rows.append([label, f"{x['Low median [IQR]']}; n={int(x['Low n'])}",
                     f"{x['High median [IQR]']}; n={int(x['High n'])}",
                     fmt(x["Rank-biserial r (positive means higher in High subgroup)"], 2),
                     qfmt(float(x["FDR q (all clinical comparisons)"]))])
    for _, x in categorical.iterrows():
        label = "Female sex" if x["Variable"] == "SEX" else x["Variable"]
        low_dist = str(x["Low distribution"]).replace("女:", "Female:").replace("男:", "Male:")
        high_dist = str(x["High distribution"]).replace("女:", "Female:").replace("男:", "Male:")
        rows.append([label, low_dist, high_dist,
                     "—", qfmt(float(x["FDR q (all clinical comparisons)"]))])
    for _, x in biomarkers.iterrows():
        label = {"p-Tau217": "p-tau217", "p-Tau231": "p-tau231", "BD Tau": "BD-tau", "NFL": "NfL"}.get(x["Variable"], x["Variable"])
        rows.append([label, f"{x['Low median [IQR]']}; n={int(x['Low n'])}",
                     f"{x['High median [IQR]']}; n={int(x['High n'])}",
                     fmt(x["Rank-biserial r (positive means higher in High subgroup)"], 2),
                     qfmt(float(x["FDR q (6 biomarkers)"]))])
    return "# Supplementary Table S7. Characteristics of eMTBR-tau243-defined A+T− groups\n\n" + table(
        ["Characteristic", "Low (n=12)", "High (n=21)", "Rank-biserial r", "FDR q"], rows
    ) + ("\n\nThe group-defining eMTBR-tau243 concentration is omitted as an independent comparison. "
         "Clinical/demographic and six non-defining plasma-biomarker comparisons form separate FDR families. "
         "Reference composition: tau PET-confirmed T−, 9 Low and 9 High; CSF-defined T−, 3 Low and 12 High. "
         "The contemporaneous CSF subset (Low 6, High 8) is exploratory and should be presented separately. "
         "Plasma concentrations are in pg/mL.\n")


def _extreme_value_range(path: Path, group_col: str, group_value: str) -> str:
    d = pd.read_csv(path)
    d = d.loc[(d[group_col] == group_value) & (d["Scenario"] != "All observations"), "Delta BIC"]
    return f"{fmt(d.min(), 2)} to {fmt(d.max(), 2)}"


def build_tableS4() -> str:
    full = pd.read_csv(OUT / "bimodality_validation_summary.csv").set_index("Marker")
    ref = pd.read_csv(OUT / "GMM_reference_stratified_summary.csv").set_index("Reference stratum")

    rows = []

    m = full.loc["MTBR-Tau243"]
    rows.append([
        "Full A+T− sample, eMTBR-tau243 (primary)",
        str(int(m["N"])), str(int(m["Observed low n"])), str(int(m["Observed high n"])),
        fmt(m["Observed delta BIC"], 2), pfmt(m["Parametric-bootstrap LRT p under 1-component null"]), "—",
        fmt(m["Observed decision boundary"], 3),
        f"{fmt(m['LOO minimum delta BIC'], 2)} to {fmt(m['LOO maximum delta BIC'], 2)}",
        _extreme_value_range(OUT / "bimodality_extreme_value_sensitivity.csv", "Marker", "MTBR-Tau243"),
    ])

    m = full.loc["p-Tau217"]
    rows.append([
        "Full A+T− sample, p-tau217 (negative comparison)",
        str(int(m["N"])), str(int(m["Observed low n"])), str(int(m["Observed high n"])),
        fmt(m["Observed delta BIC"], 2), pfmt(m["Parametric-bootstrap LRT p under 1-component null"]), "—",
        fmt(m["Observed decision boundary"], 3),
        f"{fmt(m['LOO minimum delta BIC'], 2)} to {fmt(m['LOO maximum delta BIC'], 2)}",
        _extreme_value_range(OUT / "bimodality_extreme_value_sensitivity.csv", "Marker", "p-Tau217"),
    ])

    for stratum, label in [("Tau PET-negative", "Reference-stratified, tau PET-confirmed T−"),
                            ("CSF-defined T-negative", "Reference-stratified, CSF-defined T−")]:
        m = ref.loc[stratum]
        rows.append([
            label,
            str(int(m["N"])), str(int(m["Observed low n"])), str(int(m["Observed high n"])),
            fmt(m["Observed delta BIC"], 2), pfmt(m["Parametric-bootstrap LRT p"]), qfmt(float(m["Bootstrap LRT BH-FDR q"])),
            fmt(m["Observed decision boundary"], 3),
            f"{fmt(m['LOO minimum delta BIC'], 2)} to {fmt(m['LOO maximum delta BIC'], 2)}",
            _extreme_value_range(OUT / "GMM_reference_stratified_extreme_sensitivity.csv", "Reference stratum", stratum),
        ])

    stability = pd.read_csv(OUT / "MTBR_subgroup_assignment_stability.csv")
    low = stability.loc[stability["Observed subgroup"] == "Low MTBR-Tau243", "Parametric-bootstrap high-assignment frequency"]
    high = stability.loc[stability["Observed subgroup"] == "High MTBR-Tau243", "Parametric-bootstrap high-assignment frequency"]

    return (
        "# Supplementary Table S6. Gaussian-mixture-model diagnostics and sensitivity analyses for eMTBR-tau243 bimodality within A+T−\n\n"
        + table(
            ["Model", "N", "Low n", "High n", "ΔBIC (2 vs 1 component)", "Parametric-bootstrap LRT p",
             "BH-FDR q", "Decision boundary", "Leave-one-out ΔBIC range", "Extreme-value ΔBIC range"],
            rows,
        )
        + (
            "\n\nΔBIC = BIC for the one-component model minus BIC for the two-component model; positive values favor "
            "two components. Parametric-bootstrap p values use 1,000 replicates under the one-component null. BH-FDR q is "
            "shown only where multiplicity correction was applied (jointly across the two reference-stratified tests); the "
            "two full-sample tests were each a single prespecified comparison and were not corrected against each other. "
            "Leave-one-out and extreme-value (remove-minimum/remove-maximum/remove-both/trim-two-per-tail) ranges describe "
            "the influence of individual or extreme observations on ΔBIC, not confidence intervals. Subgroup-assignment "
            "stability (parametric-bootstrap posterior high-component assignment frequency, full-sample model): Low-eMTBR-tau243 "
            f"participants (n={len(low)}) had frequency ≤{fmt(low.max(), 3)} (mean {fmt(low.mean(), 4)}); "
            f"High-eMTBR-tau243 participants (n={len(high)}) had frequency ≥{fmt(high.min(), 3)} (mean {fmt(high.mean(), 4)}), "
            "indicating highly stable component assignment. Source files: bimodality_validation_summary.csv, "
            "bimodality_leave_one_out.csv, bimodality_extreme_value_sensitivity.csv, MTBR_subgroup_assignment_stability.csv, "
            "GMM_reference_stratified_summary.csv, GMM_reference_stratified_leave_one_out.csv, "
            "GMM_reference_stratified_extreme_sensitivity.csv.\n"
        )
    )


# Supplementary Table S2 generator (added 2026-11-05; previously hand-maintained).
# Source: best_vs_others_AUC_comparisons.csv written by 01_diagnostic_validation.py.
S2_ANALYSES = [  # (analysis key in CSV, display label) in Results order, primary first
    ("A status, non-control cohort", "A status, healthy controls excluded (n=101; primary)"),
    ("A status", "A status, all participants (n=135)"),
    ("T status, Tau PET subset", "T status, tau PET (n=68; primary)"),
    ("T status", "T status, tau PET/CSF (n=114)"),
    ("T status, A+ subset, Tau PET only", "T status within A+, tau PET (n=52; primary)"),
    ("T status, A+ subset", "T status within A+, tau PET/CSF (n=67)"),
]
S2_DISPLAY = {
    "MTBR-Tau243": "eMTBR-tau243", "MTBR-Tau243/Aβ42": "eMTBR-tau243/Aβ42",
    "p-Tau217": "p-tau217", "p-Tau217/Aβ42": "p-tau217/Aβ42",
    "p-Tau231": "p-tau231", "p-Tau231/Aβ42": "p-tau231/Aβ42",
    "BD Tau": "BD-tau", "BD-Tau/Aβ42": "BD-tau/Aβ42", "NFL": "NfL",
    "GFAP": "GFAP", "Aβ1-42": "Aβ1-42",
}


def build_tableS2() -> str:
    data = pd.read_csv(OUT / "best_vs_others_AUC_comparisons.csv")
    rows = []
    for key, label in S2_ANALYSES:
        frame = data[data["Analysis"] == key].sort_values("delong_z", ascending=False, kind="mergesort")
        assert len(frame) == 10, (key, len(frame))
        best = S2_DISPLAY[frame["Best marker"].iloc[0]]
        for _, r in frame.iterrows():
            other = S2_DISPLAY[r["Other marker"]]
            rows.append([
                f"{label} ({best})", f"{best} − {other}", fmt(r["delta_auc"]),
                pfmt(r["delong_p"]), qfmt(r["DeLong FDR q (10 comparisons within analysis)"]),
                f"{fmt(r['bootstrap_ci_lower'])} to {fmt(r['bootstrap_ci_upper'])}",
                qfmt(r["Bootstrap FDR q (10 comparisons within analysis)"]),
            ])
    return (
        "# Supplementary Table S2. Leading-measure screening: best-performing measure versus each of the other ten measures\n\n"
        + table(["Analysis (best marker)", "Comparison (best − other)", "ΔAUC", "DeLong p",
                 "DeLong q (10 comparisons)", "Paired-bootstrap 95% CI", "Bootstrap q"], rows)
        + "\n\nWithin each analysis population, the marker with the highest apparent (full-sample) AUC was identified first, "
        "then compared against each of the other ten markers in the panel using paired DeLong tests and 5,000 paired "
        "participant-level bootstrap resamples (same participants/draw for both markers). Positive ΔAUC favors the best marker. "
        "BH-FDR was applied to the ten DeLong tests and, separately, the ten bootstrap tests within each analysis. "
        "Rows within each analysis are ordered by DeLong z statistic. This is the first step of the two-step comparison (Methods 2.4): "
        "measures not significantly lower than the best measure by DeLong test (q≥0.05) form the leading group. "
        "The head-to-head eMTBR-tau243 versus p-tau217 comparison (second step) is reported in the Results text.\n"
    )


def main() -> None:
    outputs = {
        "Table1_participant_characteristics.md": build_table1(),
        "TableS2_best_marker_vs_others.md": build_tableS2(),  # generated since 2026-11-05
        "TableS4_biomarker_performance_full.md": build_tableS3(),  # Supplementary Table S4 since 2026-11-05
        "TableS7_Aplus_Tminus_groups.md": build_tableS5(),  # Supplementary Table S7 since v3.28
        "TableS6_GMM_sensitivity_analyses.md": build_tableS4(),  # Supplementary Table S6 since v3.28
    }
    # NOTE: Supplementary Table S1 (best-marker-vs-others comparison,
    # TableS1_best_marker_vs_others.md -- moved from main Table 2 on 2026-11-01, renumbered
    # from S6 to S1 the same day to match true first-citation order; see
    # 00_story_outline_and_analysis_lock.md changelog 3.9-3.10) and Supplementary Table S6
    # (contemporaneous CSF comparison, TableS6_Aplus_Tminus_CSF_comparison.md, renumbered
    # from S5) are intentionally not built here -- the former is hand-maintained (its old
    # source CSV, TableS7_best_vs_others_AUC_comparisons.csv, was deleted as a redundant
    # duplicate), and the latter is written directly by
    # 03_Aplus_Tminus_subgroup_comparison.py's write_manuscript_table_s5().
    # See 00_story_outline_and_analysis_lock.md Section 10 for the authoritative table list.
    for filename, content in outputs.items():
        (OUT / filename).write_text(content, encoding="utf-8")
    print(f"Wrote {len(outputs)} formatted manuscript table drafts to outputs/: " + ", ".join(outputs))


if __name__ == "__main__":
    main()
