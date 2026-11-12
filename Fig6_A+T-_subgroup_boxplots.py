#!/usr/bin/env python3
"""Main-text figure: eMTBR-tau243-defined Low vs High A+T- subgroup
comparisons within the 18 tau PET-defined A+T− patients -- clinical/demographic
and plasma biomarkers. Pure plotting from already-locked statistics (medians/IQR/FDR q read
from the existing CSVs) and raw per-participant values (joined from
biomarker_result.xlsx by SID) -- no new statistical test is computed here;
all q values are the already-locked FDR q values from
A+T-Subgroup_continuous.csv, A+T-Subgroup_categorical.csv,
and A+T-Subgroup_biomarkers.csv.

Two independently FDR-corrected families are shown as two rows:
  a) 5 clinical/demographic variables (n=18)
  b) 6 non-defining plasma biomarkers (n=18)

Panels carry no in-image captions or descriptive titles -- only bare
a/b panel letters and per-panel y-axis labels naming each measure.
Full descriptive text lives in the manuscript figure legend.
"""
import json
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl
mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

DEFAULT_INPUT = SCRIPT_DIR.parent / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(SCRIPT_DIR.parent / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# v3.60 (2026-11-11): tau PET-defined A+T− patients only (n=18); the pooled mode and the CSF
# row (only 6 of 18 with contemporaneous CSF; not evaluable) were removed. Statistics are read
# from the A+T− analysis folder; the figure is written to OUTPUT_DIR.
STATS_DIR = Path(os.environ.get("BIOMARKER_APTM_DIR", str(OUTPUT_DIR / "Aplus_Tminus_tauPET_only")))
EXPECTED_N = 18

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)

LOW_COLOR = "#4C78A8"
HIGH_COLOR = "#C44E52"
LOW = "Low MTBR-Tau243"
HIGH = "High MTBR-Tau243"

PLASMA_COL = {
    "p-Tau217": "p-Tau217(≤0.758)",
    "BD Tau": "BD Tau（≤52.301）",
    "p-Tau231": "p-Tau231(≤5.793)",
    "GFAP": "GFAP(≤131.880)",
    "Aβ1-42": "Aβ 1-42",
    "NFL": "NFL(≤9.624 )",
}
PLASMA_DISPLAY = {
    "p-Tau217": "p-tau217", "BD Tau": "BD-tau", "p-Tau231": "p-tau231",
    "GFAP": "GFAP", "Aβ1-42": "Aβ1-42", "NFL": "NfL",
}
PLASMA_ORDER = ["p-Tau217", "p-Tau231", "Aβ1-42", "BD Tau", "GFAP", "NFL"]


membership = pd.read_csv(STATS_DIR / "A+T-Subgroup_membership.csv")
df = pd.read_excel(INPUT, sheet_name="Sheet1")


def apoe4_carrier(value):
    if pd.isna(value):
        return np.nan
    return "ε4" in str(value)


df["APOE4_carrier"] = df["apoe"].map(apoe4_carrier)
merged = membership.merge(df, on="SID", how="left", validate="one_to_one")
assert len(merged) == EXPECTED_N, len(merged)

q_clin_cont = pd.read_csv(STATS_DIR / "A+T-Subgroup_continuous.csv").set_index("Variable")["FDR q (all clinical comparisons)"]
q_clin_cat = pd.read_csv(STATS_DIR / "A+T-Subgroup_categorical.csv").set_index("Variable")["FDR q (all clinical comparisons)"]
q_plasma = pd.read_csv(STATS_DIR / "A+T-Subgroup_biomarkers.csv").set_index("Variable")["FDR q (6 biomarkers)"]



def qfmt(v: float) -> str:
    return "q<0.001" if v < 0.001 else f"q={v:.3f}"


def sig_bracket(ax, y0, y1, q_value):
    ax.plot([1, 1, 2, 2], [y0, y1, y1, y0], color="black", linewidth=0.7, clip_on=False)
    sig = q_value < 0.05
    label = qfmt(q_value) + ("*" if sig else "")
    ax.text(1.5, y1, label, ha="center", va="bottom", fontsize=6.2,
             fontweight="bold" if sig else "normal",
             color="#1a1a1a" if sig else "#777777", clip_on=False)


def box_panel(ax, low_vals, high_vals, label, q_value, log=False):
    rng = np.random.default_rng(abs(hash(label)) % (2**32))
    bp = ax.boxplot(
        [low_vals, high_vals], positions=[1, 2], widths=0.5, showfliers=False,
        patch_artist=True,
        medianprops={"color": "black", "linewidth": 1.1},
        whiskerprops={"linewidth": 0.9, "color": "#555555"},
        capprops={"linewidth": 0.9, "color": "#555555"},
        boxprops={"linewidth": 0.9, "edgecolor": "#555555"},
    )
    for patch, color in zip(bp["boxes"], [LOW_COLOR, HIGH_COLOR]):
        patch.set_facecolor(color)
        patch.set_alpha(0.28)
    for vals, pos, color in zip([low_vals, high_vals], [1, 2], [LOW_COLOR, HIGH_COLOR]):
        jitter = rng.uniform(-0.09, 0.09, size=len(vals))
        ax.scatter(np.full(len(vals), pos) + jitter, vals, s=11, color=color,
                   edgecolor="white", linewidth=0.35, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Low\n(n={len(low_vals)})", f"High\n(n={len(high_vals)})"], fontsize=6.3)
    if log:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f"{y:g}"))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylabel(label, fontsize=6.8, fontweight="bold")
    lo = min(np.min(low_vals), np.min(high_vals))
    hi = max(np.max(low_vals), np.max(high_vals))
    if log:
        y0 = hi * 1.12
        y1 = hi * 1.28
        ax.set_ylim(top=hi * 1.9)
    else:
        span = hi - lo if hi > lo else (abs(hi) if hi else 1.0)
        y0 = hi + 0.06 * span
        y1 = y0 + 0.03 * span
        ax.set_ylim(top=hi + 0.32 * span)
    sig_bracket(ax, y0, y1, q_value)
    ax.tick_params(axis="y", labelsize=6.0)


def bar_panel(ax, low_counts, high_counts, categories, label, q_value, bar_colors):
    low_total = sum(low_counts)
    high_total = sum(high_counts)
    low_pct = [100 * c / low_total for c in low_counts]
    high_pct = [100 * c / high_total for c in high_counts]
    x = np.arange(len(categories))
    width = 0.32
    ax.bar(x - width / 2, low_pct, width=width, color=LOW_COLOR, alpha=0.75,
           edgecolor="#555555", linewidth=0.7, label=f"Low (n={low_total})")
    ax.bar(x + width / 2, high_pct, width=width, color=HIGH_COLOR, alpha=0.75,
           edgecolor="#555555", linewidth=0.7, label=f"High (n={high_total})")
    for xi, (lp, hp) in enumerate(zip(low_pct, high_pct)):
        ax.text(xi - width / 2, lp + 2, f"{lp:.0f}%", ha="center", va="bottom", fontsize=5.6)
        ax.text(xi + width / 2, hp + 2, f"{hp:.0f}%", ha="center", va="bottom", fontsize=5.6)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=6.2)
    ax.set_ylim(0, 100)
    ax.set_ylabel(f"{label}, % of group", fontsize=6.4, fontweight="bold")
    sig = q_value < 0.05
    qlabel = qfmt(q_value) + ("*" if sig else "")
    ax.text(0.5, 1.07, qlabel, ha="center", va="bottom", fontsize=6.2, transform=ax.transAxes,
             fontweight="bold" if sig else "normal", color="#1a1a1a" if sig else "#777777")
    ax.tick_params(axis="y", labelsize=6.0)


N_ROWS = 2
fig = plt.figure(figsize=(9.8, 6.4))
gs = fig.add_gridspec(N_ROWS, 6, left=0.07, right=0.98, top=0.90, bottom=0.07, hspace=1.05, wspace=0.6)

# --- Row a: clinical / demographic (joint FDR family of 5) ---
fig.text(0.015, 0.975, "a", fontsize=11, fontweight="bold")

low_m = merged["Plasma tau subgroup"].eq(LOW)
high_m = merged["Plasma tau subgroup"].eq(HIGH)

ax = fig.add_subplot(gs[0, 0])
box_panel(ax, merged.loc[low_m, "AGE"].dropna().to_numpy(), merged.loc[high_m, "AGE"].dropna().to_numpy(),
          "Age, years", q_clin_cont["AGE"])

ax = fig.add_subplot(gs[0, 1])
box_panel(ax, merged.loc[low_m, "MMSE_Total"].dropna().to_numpy(), merged.loc[high_m, "MMSE_Total"].dropna().to_numpy(),
          "MMSE", q_clin_cont["MMSE_Total"])

ax = fig.add_subplot(gs[0, 2])
box_panel(ax, merged.loc[low_m, "MOCA_P"].dropna().to_numpy(), merged.loc[high_m, "MOCA_P"].dropna().to_numpy(),
          "MoCA-P", q_clin_cont["MOCA_P"])

ax = fig.add_subplot(gs[0, 3])
sex_low = merged.loc[low_m, "SEX"].dropna()
sex_high = merged.loc[high_m, "SEX"].dropna()
bar_panel(ax, [int((sex_low == "女").sum()), int((sex_low == "男").sum())],
          [int((sex_high == "女").sum()), int((sex_high == "男").sum())],
          ["Female", "Male"], "Female sex", q_clin_cat["SEX"], None)

ax = fig.add_subplot(gs[0, 4])
apoe_low = merged.loc[low_m, "APOE4_carrier"].dropna()
apoe_high = merged.loc[high_m, "APOE4_carrier"].dropna()
bar_panel(ax, [int((apoe_low == False).sum()), int((apoe_low == True).sum())],
          [int((apoe_high == False).sum()), int((apoe_high == True).sum())],
          ["No ε4", "ε4+"], "APOE ε4 carrier", q_clin_cat["APOE ε4 carrier"], None)

fig.add_subplot(gs[0, 5]).axis("off")

# --- Row b: plasma biomarkers (linear scale; joint FDR family of 6) ---
fig.text(0.015, 0.47, "b", fontsize=11, fontweight="bold")

for i, var in enumerate(PLASMA_ORDER):
    ax = fig.add_subplot(gs[1, i])
    col = PLASMA_COL[var]
    low_vals = merged.loc[low_m, col].dropna().to_numpy()
    high_vals = merged.loc[high_m, col].dropna().to_numpy()
    box_panel(ax, low_vals, high_vals, f"{PLASMA_DISPLAY[var]} (pg/mL)", q_plasma[var], log=False)

from matplotlib.patches import Patch
legend_elems = [
    Patch(facecolor=LOW_COLOR, alpha=0.6, edgecolor="#555555", label="Low eMTBR-tau243"),
    Patch(facecolor=HIGH_COLOR, alpha=0.6, edgecolor="#555555", label="High eMTBR-tau243"),
]
fig.legend(handles=legend_elems, loc="upper right", bbox_to_anchor=(0.99, 0.995), fontsize=7.2, ncol=2)

stem = "Fig6_A+T-_subgroup_boxplots"
for ext, kwargs in [
    ("svg", {}),
    ("pdf", {}),
    ("png", {"dpi": 300}),
    ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}),
]:
    fig.savefig(OUTPUT_DIR / f"{stem}.{ext}", bbox_inches="tight", facecolor="white", **kwargs)

# Source data: per-participant values underlying every panel, for traceability.
source_cols = {
    "SID": merged["SID"],
    "Plasma tau subgroup": merged["Plasma tau subgroup"],
    "AGE": merged["AGE"],
    "MMSE_Total": merged["MMSE_Total"],
    "MOCA_P": merged["MOCA_P"],
    "SEX": merged["SEX"],
    "APOE4_carrier": merged["APOE4_carrier"],
}
for var, col in PLASMA_COL.items():
    source_cols[f"plasma_{var}"] = merged[col]
source_df = pd.DataFrame(source_cols)
source_df.to_csv(OUTPUT_DIR / f"{stem}_source_data.csv", index=False)

metadata = {
    "figure": "eMTBR-tau243-defined Low vs High A+T- subgroup comparisons",
    "n_main": int(len(merged)),
    "fdr_families": {
        "clinical_demographic": ["AGE", "SEX", "APOE4_carrier", "MMSE_Total", "MOCA_P"],
        "plasma_biomarkers": PLASMA_ORDER,
    },
    "scales": {"row_a": "linear", "row_b": "linear"},
    "source_stats": [
        "A+T-Subgroup_continuous.csv", "A+T-Subgroup_categorical.csv",
        "A+T-Subgroup_biomarkers.csv",
        "A+T-Subgroup_membership.csv",
    ],
}
(OUTPUT_DIR / f"{stem}_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

print(f"wrote {stem} (svg/pdf/png/tiff) + source data + metadata to {OUTPUT_DIR}")
