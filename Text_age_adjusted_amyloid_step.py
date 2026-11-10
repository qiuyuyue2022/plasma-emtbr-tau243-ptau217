#!/usr/bin/env python3
"""Age-adjusted amyloid-step sensitivity analysis (Stage 4.5 check, 2026-11-09).

Within tau PET-defined T− participants (A+T− n=18 vs A−T− n=16), logistic regression of
A status on each standardized log10 marker plus standardized age; the likelihood-ratio test
compares marker+age with age alone. Also reports the age-only AUC and, within all A−T−
participants (n=47), Spearman correlations of age with each marker. Read-only.
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(ROOT / "biomarker_result.xlsx")))
OUT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(ROOT / "outputs")))
MARKERS = {"eMTBR-tau243": "MTBR-Tau243", "p-tau217": "p-Tau217(≤0.758)", "p-tau231": "p-Tau231(≤5.793)",
           "BD-tau": "BD Tau（≤52.301）", "GFAP": "GFAP(≤131.880)", "NfL": "NFL(≤9.624 )", "Aβ1-42": "Aβ 1-42"}


def clean(s):
    return s.astype("string").str.strip().str.replace("−", "-", regex=False)


d = pd.read_excel(INPUT, sheet_name="Sheet1")
a, t = clean(d["A"]), clean(d["T"])
pet = pd.to_numeric(d["taupet"], errors="coerce").eq(1)
m = t.eq("T-") & pet
s = d[m]
y = a[m].eq("A+").astype(int).to_numpy()
age = pd.to_numeric(s["AGE"]).to_numpy()
za = (age - age.mean()) / age.std()
base = sm.Logit(y, sm.add_constant(za)).fit(disp=0)
rows = []
for name, col in MARKERS.items():
    x = np.log10(pd.to_numeric(s[col]).to_numpy()); z = (x - x.mean()) / x.std()
    r = sm.Logit(y, sm.add_constant(np.column_stack([z, za]))).fit(disp=0)
    lr = 2 * (r.llf - base.llf)
    aa = a.eq("A-") & t.eq("T-")
    rho, p_rho = spearmanr(pd.to_numeric(d.loc[aa, "AGE"]), pd.to_numeric(d.loc[aa, col]))
    rows.append({"Marker": name, "OR per SD (age-adjusted)": float(np.exp(r.params[1])), "Wald p": float(r.pvalues[1]),
                 "LRT p over age": float(chi2.sf(lr, 1)), "Age p": float(r.pvalues[2]),
                 "A−T− (n=47) Spearman rho with age": float(rho), "rho p": float(p_rho)})
res = pd.DataFrame(rows)
res.to_csv(OUT / "age_adjusted_amyloid_step.csv", index=False)
md = ["# Age-adjusted amyloid step (tau PET-defined T−: A+T− 18 vs A−T− 16)", "",
      f"Age-only AUC for A status: {roc_auc_score(y, age):.3f}; median age A−T− {np.median(age[y == 0]):.1f}, A+T− {np.median(age[y == 1]):.1f}.", "",
      "| Marker | OR per SD (age-adjusted) | Wald p | LRT p over age | Age p | A−T− rho with age (p) |", "|---|---|---|---|---|---|"]
for _, r in res.iterrows():
    md.append(f"| {r.Marker} | {r['OR per SD (age-adjusted)']:.2f} | {r['Wald p']:.3f} | {r['LRT p over age']:.3f} | {r['Age p']:.3f} | {r['A−T− (n=47) Spearman rho with age']:.2f} ({r['rho p']:.3f}) |")
(OUT / "age_adjusted_amyloid_step.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print("\n".join(md))
