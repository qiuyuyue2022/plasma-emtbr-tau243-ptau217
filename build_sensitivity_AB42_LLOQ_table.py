"""Draft supplementary table: plasma Aβ1-42 LLOQ-censoring sensitivity analysis (post hoc).

Compares the locked main outputs (outputs/) with a rerun of the same scripts
(Fig3_ROC.py, 01_diagnostic_validation.py, 03_Aplus_Tminus_subgroup_comparison.py,
Fig5_AT_3group_distributions.py) on a copy of biomarker_result.xlsx in which plasma
Aβ1-42 values below the LLOQ (1.546 pg/mL; 23/135) were set to the LLOQ. Rerun outputs
and the censored input live in outputs/sensitivity_AB42_LLOQ_censoring/.
Numbered Supplementary Table S5 (first cited in Results 3.3); former S5-S7 became S6-S8 (v3.28).
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "outputs"
SENS = MAIN / "sensitivity_AB42_LLOQ_censoring"
OUTFILE = MAIN / "TableS5_AB42_LLOQ_sensitivity.md"

ANALYSES = [("A status, non-control cohort", "A, healthy controls excluded (n=101; primary)"),
            ("A status", "A, all participants (n=135)"),
            ("T status, Tau PET subset", "T, tau PET (n=68; primary)"),
            ("T status", "T, tau PET/CSF (n=114)"),
            ("T status, A+ subset, Tau PET only", "T within A+, tau PET (n=52; primary)"),
            ("T status, A+ subset", "T within A+, tau PET/CSF (n=67)")]
NAME = {"MTBR-Tau243": "eMTBR-tau243", "MTBR-Tau243/Aβ42": "eMTBR-tau243/Aβ42", "p-Tau217": "p-tau217",
        "p-Tau217/Aβ42": "p-tau217/Aβ42", "p-Tau231/Aβ42": "p-tau231/Aβ42", "BD-Tau/Aβ42": "BD-tau/Aβ42",
        "Aβ1-42": "Aβ1-42", "NFL": "NfL", "BD Tau": "BD-tau", "p-Tau231": "p-tau231", "GFAP": "GFAP"}
CONTRAST_UNUSED = {"MTBR ratio vs absolute": "eMTBR-tau243/Aβ42 − eMTBR-tau243",
            "p-Tau217 ratio vs absolute": "p-tau217/Aβ42 − p-tau217",
            "MTBR vs p-Tau217 absolute": "eMTBR-tau243 − p-tau217",
            "MTBR vs p-Tau217 ratio": "eMTBR-tau243/Aβ42 − p-tau217/Aβ42"}


def f3(x): return f"{x:.3f}"
def q3(x): return "<0.001" if x < 0.001 else f"{x:.3f}"
def table(h, rows): return "| " + " | ".join(h) + " |\n|" + "---|" * len(h) + "\n" + "\n".join("| " + " | ".join(r) + " |" for r in rows)


def panel_a():
    m, c = (pd.read_csv(d / "biomarker_performance_results.csv") for d in (MAIN, SENS))
    rows = []
    for key, lab in ANALYSES:
        for b in ["Aβ1-42", "MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "p-Tau231/Aβ42", "BD-Tau/Aβ42"]:
            x = m[(m.Analysis == key) & (m.Biomarker == b)].iloc[0]; y = c[(c.Analysis == key) & (c.Biomarker == b)].iloc[0]
            ci = lambda r: f"{f3(r.AUC)} ({f3(r['AUC 95% CI lower'])}–{f3(r['AUC 95% CI upper'])})"
            rows.append([lab, NAME[b], ci(x), ci(y), f"{y.AUC - x.AUC:+.3f}"])
    return table(["Analysis", "Measure", "Main AUC (95% CI)", "Censored AUC (95% CI)", "Change"], rows)


def panel_b():
    """Leading group (step 1 of the two-step comparison) before and after censoring."""
    m, c = (pd.read_csv(d / "best_vs_others_AUC_comparisons.csv") for d in (MAIN, SENS))
    qd, qb = "DeLong FDR q (10 comparisons within analysis)", "Bootstrap FDR q (10 comparisons within analysis)"
    def tier(d, key, col):
        x = d[d.Analysis == key]
        return [NAME[x["Best marker"].iloc[0]]] + [NAME[o] for o in x.loc[x[col] >= 0.05, "Other marker"]]
    rows = []
    for key, lab in ANALYSES:
        rows.append([lab, ", ".join(tier(m, key, qd)), ", ".join(tier(c, key, qd)), ", ".join(tier(c, key, qb))])
    return table(["Analysis", "Main: leading group (DeLong)", "Censored: leading group (DeLong)", "Censored: leading group (bootstrap)"], rows)


def panel_c():
    rows = []
    m, c = (pd.read_csv(d / "A+T-Subgroup_biomarkers.csv") for d in (MAIN, SENS))
    rc = "Rank-biserial r (positive means higher in High subgroup)"; qc = "FDR q (6 biomarkers)"
    x = m[m.Variable == "Aβ1-42"].iloc[0]; y = c[c.Variable == "Aβ1-42"].iloc[0]
    rows.append(["A+T− High vs Low eMTBR-tau243 (Results 3.7)", "Aβ1-42", f"High {x['High median [IQR]']} vs Low {x['Low median [IQR]']}; r={x[rc]:.2f}; q={q3(x[qc])}",
                 f"High {y['High median [IQR]']} vs Low {y['Low median [IQR]']}; r={y[rc]:.2f}; q={q3(y[qc])}"])
    m, c = (pd.read_csv(d / "Fig5_AT_3group_distributions_statistics.csv") for d in (MAIN, SENS))
    rc = "Rank-biserial r (positive means higher in A+T+)"; qc = "FDR q (11 biomarkers)"
    for b in ["Aβ1-42", "MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "p-Tau231/Aβ42", "BD-Tau/Aβ42"]:
        x = m[m.Biomarker == b].iloc[0]; y = c[c.Biomarker == b].iloc[0]
        rows.append(["A+T+ vs A+T− (Results 3.5; Figure 5)", NAME[b], f"r={x[rc]:.2f}; q={q3(x[qc])}", f"r={y[rc]:.2f}; q={q3(y[qc])}"])
    return table(["Analysis", "Measure", "Main", "Censored"], rows)


def main():
    text = ("# Supplementary Table S5. Sensitivity analysis censoring plasma Aβ1-42 at the lower limit of quantification\n\n"
            "**A. AUCs of plasma Aβ1-42 and the four Aβ42-normalized ratios**\n\n" + panel_a() +
            "\n\n**B. Leading group in each analysis population (first step of the two-step comparison; Supplementary Table S2)**\n\n" + panel_b() +
            "\n\n**C. Other analyses involving plasma Aβ1-42**\n\n" + panel_c() +
            "\n\nPost hoc sensitivity analysis. Plasma Aβ1-42 values below the LLOQ (1.546 pg/mL; 23/135 participants, 9 below the LOD) were set to the LLOQ, "
            "the four Aβ42-normalized ratios were recomputed, and all analyses involving Aβ1-42 were repeated with unchanged code, seeds, populations, and multiplicity families. "
            "Measures not involving Aβ1-42 were unchanged (maximum AUC difference 0), so the head-to-head eMTBR-tau243 versus p-tau217 comparison was unaffected. The leading group comprises the best measure and all measures not significantly lower than it (q≥0.05). "
            "Main analyses use values as reported by the laboratory; for research samples, values below the LOD were extrapolated from the calibration curve.\n")
    OUTFILE.write_text(text, encoding="utf-8")
    print(f"wrote {OUTFILE}")


if __name__ == "__main__":
    main()
