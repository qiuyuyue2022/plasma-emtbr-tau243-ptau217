#!/usr/bin/env python3
"""Plasma storage duration and blood-to-CSF intervals (Revision Roadmap REV-19/REV-20).

Added 2026-11-09. Storage time = assay date (2026-09-22, single batch) minus
sample_date. Reports its distribution and Spearman correlation with each plasma
biomarker (all participants), and the blood-to-CSF interval for participants
whose reference standard was CSF. Tau/amyloid PET dates are not in the dataset;
the authors report that all PET scans were within 30 days of blood collection.
"""
from pathlib import Path
import os
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(ROOT / "biomarker_result.xlsx")))
OUT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(ROOT / "outputs")))
ASSAY_DATE = pd.Timestamp("2026-09-22")
MARKERS = {"MTBR-tau243": "MTBR-Tau243", "p-tau217": "p-Tau217(≤0.758)", "p-tau231": "p-Tau231(≤5.793)",
           "BD-tau": "BD Tau（≤52.301）", "GFAP": "GFAP(≤131.880)", "NfL": "NFL(≤9.624 )", "Aβ1-42": "Aβ 1-42"}

d = pd.read_excel(INPUT, sheet_name="Sheet1")
years = (ASSAY_DATE - pd.to_datetime(d["sample_date"], errors="coerce")).dt.days / 365.25
lines = ["# Plasma storage duration and blood-to-reference intervals", "",
         f"Storage (n={years.notna().sum()}): median {years.median():.1f} years (IQR {years.quantile(.25):.1f}–{years.quantile(.75):.1f}; "
         f"range {years.min():.1f}–{years.max():.1f}).", "", "| Biomarker | Spearman ρ with storage years | p |", "|---|---|---|"]
for name, col in MARKERS.items():
    v = pd.to_numeric(d[col], errors="coerce")
    ok = v.notna() & years.notna()
    rho, p = spearmanr(years[ok], v[ok])
    lines.append(f"| {name} | {rho:.2f} | {p:.3f} |")
csf_ref = d["T_level"].eq("CSF") | d["A_level"].eq("CSF")
gap = pd.to_numeric(d["CSF_Lumipulse_绝对时间差_天"], errors="coerce")
g = gap[csf_ref & gap.notna()]
lines += ["", f"Blood-to-CSF interval, participants with a CSF-based A or T reference and a recorded interval (n={len(g)}): "
          f"median {g.median():.0f} days (IQR {g.quantile(.25):.0f}–{g.quantile(.75):.0f}; maximum {g.max():.0f})."]
(OUT / "storage_time_and_intervals.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))

# --- Added 2026-11-09 (REV-20/REV-16, before re-review): storage by A/T group and A+T- component,
# within-group associations, and an A-T--referenced elevation proportion in tau PET-defined A+T-.
import numpy as np
from scipy.stats import kruskal, mannwhitneyu

def _clean(s):
    return s.astype("string").str.strip().str.replace("−", "-", regex=False)

A, T = _clean(d["A"]), _clean(d["T"])
pet = pd.to_numeric(d["taupet"], errors="coerce").eq(1)
grp = pd.Series(np.select([A.eq("A-") & T.eq("T-"), A.eq("A+") & T.eq("T-"), A.eq("A+") & T.eq("T+")],
                          ["A−T−", "A+T−", "A+T+"], "T unavailable"), index=d.index)
extra = ["", "## Storage time by group", "", "| Group | n | Median years (IQR) |", "|---|---|---|"]
for k in ["A−T−", "A+T−", "A+T+", "T unavailable"]:
    y = years[grp.eq(k)]
    extra.append(f"| {k} | {len(y)} | {y.median():.1f} ({y.quantile(.25):.1f}–{y.quantile(.75):.1f}) |")
kw = kruskal(*[years[grp.eq(k)] for k in ["A−T−", "A+T−", "A+T+"]]).pvalue
extra.append(f"\nKruskal–Wallis across A−T−, A+T−, A+T+: p={kw:.3f}.")
mem_path = OUT / "Aplus_Tminus_tauPET_only" / "A+T-Subgroup_membership.csv"
if mem_path.exists():
    mem = pd.read_csv(mem_path)[["SID", "Plasma tau subgroup"]]
    x = d.assign(years=years).merge(mem, on="SID")
    parts = {k: s["years"] for k, s in x.groupby("Plasma tau subgroup")}
    for k, y in parts.items():
        extra.append(f"- Tau PET-defined A+T−, {k}: n={len(y)}, median {y.median():.1f} years (IQR {y.quantile(.25):.1f}–{y.quantile(.75):.1f}).")
    extra.append(f"- Low vs High Mann–Whitney p={mannwhitneyu(*parts.values()).pvalue:.3f}.")
mt = pd.to_numeric(d["MTBR-Tau243"], errors="coerce")
p217 = pd.to_numeric(d["p-Tau217(≤0.758)"], errors="coerce")
for lab, mask in [("Tau PET-defined A+T−", A.eq("A+") & T.eq("T-") & pet), ("A−T−", A.eq("A-") & T.eq("T-"))]:
    r1, r2 = spearmanr(years[mask], mt[mask]), spearmanr(years[mask], p217[mask])
    extra.append(f"- {lab} (n={int(mask.sum())}): Spearman ρ storage vs MTBR-tau243 {r1[0]:.2f} (p={r1[1]:.2f}); vs p-tau217 {r2[0]:.2f} (p={r2[1]:.2f}).")
ref = mt[A.eq("A-") & T.eq("T-")]
p95 = float(np.percentile(ref, 95))
b = mt[A.eq("A+") & T.eq("T-") & pet]
extra += ["", "## A−T−-referenced elevation", "",
          f"A−T− 95th percentile of MTBR-tau243: {p95:.2f} pg/mL; tau PET-defined A+T− above it: {int((b > p95).sum())}/{len(b)}."]
with open(OUT / "storage_time_and_intervals.md", "a", encoding="utf-8") as fh:
    fh.write("\n".join(extra) + "\n")
print("\n".join(extra))
