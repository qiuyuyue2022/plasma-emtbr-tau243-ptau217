#!/usr/bin/env python3
"""Compare CSF-defined and tau PET-defined A+T- participants (Stage 4' NEW-9, 2026-11-09).

Mann-Whitney U tests with rank-biserial r (positive = higher in CSF-defined) for
age, MMSE, MoCA-P, plasma p-tau217 and MTBR-tau243; BH-FDR across the five tests.
Read-only with respect to the raw data.
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(ROOT / "biomarker_result.xlsx")))
OUT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(ROOT / "outputs")))
VARS = {"Age": "AGE", "MMSE": "MMSE_Total", "MoCA-P": "MOCA_P",
        "Plasma p-tau217 (pg/mL)": "p-Tau217(≤0.758)", "Plasma MTBR-tau243 (pg/mL)": "MTBR-Tau243"}


def clean(s):
    return s.astype("string").str.strip().str.replace("−", "-", regex=False)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


d = pd.read_excel(INPUT, sheet_name="Sheet1")
a, t = clean(d["A"]), clean(d["T"])
pet = pd.to_numeric(d["taupet"], errors="coerce").eq(1)
csf_g = a.eq("A+") & t.eq("T-") & ~pet
pet_g = a.eq("A+") & t.eq("T-") & pet
assert (csf_g.sum(), pet_g.sum()) == (15, 18)
rows = []
for name, col in VARS.items():
    v = pd.to_numeric(d[col], errors="coerce")
    x, y = v[csf_g].dropna(), v[pet_g].dropna()
    u = mannwhitneyu(x, y, alternative="two-sided")
    rows.append({"Variable": name, "CSF-defined n": len(x), "CSF-defined median [IQR]": f"{x.median():.3g} [{x.quantile(.25):.3g}–{x.quantile(.75):.3g}]",
                 "Tau PET-defined n": len(y), "Tau PET-defined median [IQR]": f"{y.median():.3g} [{y.quantile(.25):.3g}–{y.quantile(.75):.3g}]",
                 "U": u.statistic, "p": u.pvalue, "rank-biserial r": 2 * u.statistic / (len(x) * len(y)) - 1})
r = pd.DataFrame(rows); r["FDR q"] = bh(r["p"])
r.to_csv(OUT / "Aplus_Tminus_reference_comparison.csv", index=False)
md = ["# CSF-defined vs tau PET-defined A+T−", "", "| Variable | CSF-defined (n) | median [IQR] | Tau PET-defined (n) | median [IQR] | r | p | q |", "|---|---|---|---|---|---|---|---|"]
for _, x in r.iterrows():
    md.append(f"| {x.Variable} | {x['CSF-defined n']} | {x['CSF-defined median [IQR]']} | {x['Tau PET-defined n']} | {x['Tau PET-defined median [IQR]']} | {x['rank-biserial r']:.2f} | {x.p:.3g} | {x['FDR q']:.3g} |")
(OUT / "Aplus_Tminus_reference_comparison.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print("\n".join(md))
