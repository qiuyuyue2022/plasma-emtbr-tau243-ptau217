import json
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl
mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("/Users/qiuyuyue/Documents/My-Papers/Blood biomarkers/biomarker_result.xlsx")
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
WORK_DIR = SCRIPT_DIR
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(SCRIPT_DIR.parent / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
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
})

LABELS = {
    "GFAP(≤131.880)": "GFAP",
    "NFL(≤9.624 )": "NFL",
    "Aβ 1-42": "Aβ1-42",
    "p-Tau217(≤0.758)": "p-Tau217",
    "p-Tau231(≤5.793)": "p-Tau231",
    "BD Tau（≤52.301）": "BD Tau",
    "MTBR-Tau243": "MTBR-Tau243",
    "p-Tau217/Aβ42": "p-Tau217/Aβ42",
    "p-Tau231/Aβ42": "p-Tau231/Aβ42",
    "MTBR-Tau243/Aβ42": "MTBR-Tau243/Aβ42",
    "BD-Tau/Aβ42": "BD-Tau/Aβ42",
}

COLORS = {
    # Ratios use saturated hues; the corresponding absolute biomarkers use
    # clearly different hues so paired measures remain separable at a glance.
    "MTBR-Tau243/Aβ42": "#0072B2",
    "p-Tau217/Aβ42": "#CC79A7",
    "BD-Tau/Aβ42": "#009E73",
    "p-Tau231/Aβ42": "#E69F00",
    "MTBR-Tau243": "#1A1A1A",
    "p-Tau217": "#D55E00",
    "BD Tau": "#56B4E9",
    "GFAP": "#7B3294",
    "p-Tau231": "#A6761D",
    "Aβ1-42": "#4D9221",
    "NFL": "#8C8C8C",
}

LINESTYLES = {
    "MTBR-Tau243/Aβ42": (0, (7, 2)),
    "p-Tau217/Aβ42": (0, (5, 1.5, 1.2, 1.5)),
    "BD-Tau/Aβ42": (0, (2.5, 1.5)),
    "p-Tau231/Aβ42": (0, (6, 1.5, 1.2, 1.5, 1.2, 1.5)),
    "MTBR-Tau243": "-",
    "p-Tau217": "-",
    "BD Tau": "--",
    "GFAP": "-.",
    "p-Tau231": ":",
    "Aβ1-42": "--",
    "NFL": ":",
}

ORDER = [
    "MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "BD-Tau/Aβ42", "p-Tau231/Aβ42",
    "MTBR-Tau243", "p-Tau217", "BD Tau", "GFAP", "p-Tau231", "Aβ1-42", "NFL",
]

# Presentation-only: the manufacturer-confirmed analyte name is eMTBR-tau243
# (see Methods 2.3). LABELS/COLORS/LINESTYLES/ORDER above, and the "Biomarker"
# column written to biomarker_performance_results.csv, deliberately keep the
# source-data label "MTBR-Tau243" unchanged, because Fig2Fig4/Fig5's cutoff_lookup
# and QC_audit_manuscript_numbers.py index that CSV by exact string match.
# This override is applied only where text is actually drawn on the figure.
DISPLAY_OVERRIDE = {  # display names follow the nomenclature lock (presentation only)
    "MTBR-Tau243": "eMTBR-tau243", "MTBR-Tau243/Aβ42": "eMTBR-tau243/Aβ42",
    "p-Tau217": "p-tau217", "p-Tau217/Aβ42": "p-tau217/Aβ42",
    "p-Tau231": "p-tau231", "p-Tau231/Aβ42": "p-tau231/Aβ42",
    "BD Tau": "BD-tau", "BD-Tau/Aβ42": "BD-tau/Aβ42", "NFL": "NfL",
}


def auc_score(y, score):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    pos = score[y == 1]
    neg = score[y == 0]
    diff = pos[:, None] - neg[None, :]
    return float((np.sum(diff > 0) + 0.5 * np.sum(diff == 0)) / diff.size)


def confusion(y, pred):
    y = np.asarray(y, dtype=int)
    pred = np.asarray(pred, dtype=int)
    tp = int(np.sum((y == 1) & (pred == 1)))
    tn = int(np.sum((y == 0) & (pred == 0)))
    fp = int(np.sum((y == 0) & (pred == 1)))
    fn = int(np.sum((y == 1) & (pred == 0)))
    return tn, fp, fn, tp


def oriented_auc(y, x):
    raw = auc_score(y, x)
    return (raw, 1) if raw >= 0.5 else (1 - raw, -1)


def roc_points(y, oriented_score):
    thresholds = np.r_[np.inf, np.sort(np.unique(oriented_score))[::-1], -np.inf]
    fpr, tpr = [], []
    for threshold in thresholds:
        pred = (oriented_score >= threshold).astype(int)
        tn, fp, fn, tp = confusion(y, pred)
        fpr.append(fp / (fp + tn))
        tpr.append(tp / (tp + fn))
    return np.asarray(fpr), np.asarray(tpr)


def optimal_metrics(y, x):
    auc, direction = oriented_auc(y, x)
    score = x * direction
    unique = np.unique(score)
    thresholds = np.r_[np.inf, (unique[:-1] + unique[1:]) / 2, -np.inf]
    best = None
    for threshold in thresholds:
        pred = (score >= threshold).astype(int)
        tn, fp, fn, tp = confusion(y, pred)
        sensitivity = tp / (tp + fn)
        specificity = tn / (tn + fp)
        candidate = (sensitivity + specificity - 1, sensitivity + specificity, threshold, pred)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    threshold = best[2]
    pred = best[3]
    tn, fp, fn, tp = confusion(y, pred)
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    return {
        "auc": auc,
        "direction_sign": direction,
        "direction": ">" if direction == 1 else "<",
        "cutoff": float(threshold / direction),
        "apparent_accuracy": (tp + tn) / len(y),
        "apparent_balanced_accuracy": (sensitivity + specificity) / 2,
        "apparent_sensitivity": sensitivity,
        "apparent_specificity": specificity,
    }


def bootstrap_ci(y, x, direction, seed, reps=3000):
    rng = np.random.default_rng(seed)
    oriented_score = x * direction
    values = []
    for _ in range(reps):
        idx = rng.integers(0, len(y), len(y))
        if np.unique(y[idx]).size < 2:
            continue
        values.append(auc_score(y[idx], oriented_score[idx]))
    return tuple(np.quantile(values, [0.025, 0.975]))


def cv_metrics(y, x, seed, repeats=50, folds=5):
    rng = np.random.default_rng(seed)
    pos_idx = np.flatnonzero(y == 1)
    neg_idx = np.flatnonzero(y == 0)
    accuracy, balanced, sensitivity, specificity = [], [], [], []
    for _ in range(repeats):
        pos_folds = np.array_split(rng.permutation(pos_idx), folds)
        neg_folds = np.array_split(rng.permutation(neg_idx), folds)
        for fold in range(folds):
            test = np.concatenate([pos_folds[fold], neg_folds[fold]])
            train = np.setdiff1d(np.arange(len(y)), test)
            trained = optimal_metrics(y[train], x[train])
            direction = trained["direction_sign"]
            pred = ((x[test] * direction) >= (trained["cutoff"] * direction)).astype(int)
            tn, fp, fn, tp = confusion(y[test], pred)
            se = tp / (tp + fn)
            sp = tn / (tn + fp)
            accuracy.append((tp + tn) / len(test))
            balanced.append((se + sp) / 2)
            sensitivity.append(se)
            specificity.append(sp)
    return {
        "cv_accuracy": float(np.mean(accuracy)),
        "cv_balanced_accuracy": float(np.mean(balanced)),
        "cv_sensitivity": float(np.mean(sensitivity)),
        "cv_specificity": float(np.mean(specificity)),
    }


df = pd.read_excel(INPUT, sheet_name="Sheet1")
base_biomarkers = [
    "GFAP(≤131.880)",
    "NFL(≤9.624 )",
    "Aβ 1-42",
    "p-Tau217(≤0.758)",
    "p-Tau231(≤5.793)",
    "BD Tau（≤52.301）",
    "MTBR-Tau243",
]
abeta42 = pd.to_numeric(df["Aβ 1-42"], errors="coerce")
ratio_numerators = {
    "p-Tau217/Aβ42": "p-Tau217(≤0.758)",
    "p-Tau231/Aβ42": "p-Tau231(≤5.793)",
    "MTBR-Tau243/Aβ42": "MTBR-Tau243",
    "BD-Tau/Aβ42": "BD Tau（≤52.301）",
}
for ratio_name, numerator in ratio_numerators.items():
    df[ratio_name] = pd.to_numeric(df[numerator], errors="coerce") / abeta42
biomarkers = base_biomarkers + list(ratio_numerators)
A_y = df["A"].eq("A+").astype(int)
T_y = df["T"].eq("T+").astype(int)

outcomes = [
    {
        "analysis": "A status",
        "title": "A+ versus A−",
        "subtitle": "All participants",
        "mask": df["A"].isin(["A+", "A-"]),
        "y": A_y,
    },
    {
        "analysis": "A status, non-control cohort",
        "title": "A+ versus A−",
        "subtitle": "HC excluded",
        "mask": df["A"].isin(["A+", "A-"]) & df["DIAG"].ne("健康对照"),
        "y": A_y,
    },
    {
        "analysis": "T status",
        "title": "T+ versus T−",
        "subtitle": "PET/CSF defined T",
        "mask": df["T"].isin(["T+", "T-"]),
        "y": T_y,
    },
    {
        "analysis": "T status, Tau PET subset",
        "title": "T+ versus T−",
        "subtitle": "PET defined T",
        "mask": df["T"].isin(["T+", "T-"]) & df["taupet"].eq(1),
        "y": T_y,
    },
    {
        "analysis": "T status, A+ subset",
        "title": "T+ versus T−",
        "subtitle": "A+ participants, PET/CSF defined T",
        "mask": df["A"].eq("A+") & df["T"].isin(["T+", "T-"]),
        "y": T_y,
    },
    {
        "analysis": "T status, A+ subset, Tau PET only",
        "title": "T+ versus T−",
        "subtitle": "A+ participants, PET defined T",
        "mask": df["A"].eq("A+") & df["T"].isin(["T+", "T-"]) & df["taupet"].eq(1),
        "y": T_y,
    },
]

# Preserve the original computation order (and therefore the established random
# seeds and numerical results), but use the manuscript narrative order for plots:
# A non-control (primary, n=101; v3.28) -> A overall -> primary tau PET T -> operational PET/CSF T -> A+ T.
PLOT_ANALYSIS_ORDER = [
    "A status, non-control cohort",
    "A status",
    "T status, Tau PET subset",
    "T status",
    "T status, A+ subset, Tau PET only",
    "T status, A+ subset",
]

rows = []
roc_data = {}
for outcome_idx, outcome in enumerate(outcomes):
    mask = outcome["mask"]
    y_all = outcome["y"][mask].to_numpy(dtype=int)
    outcome["n"] = int(len(y_all))
    outcome["positive"] = int(y_all.sum())
    outcome["negative"] = int(len(y_all) - y_all.sum())
    roc_data[outcome["analysis"]] = {}
    for biomarker_idx, biomarker in enumerate(biomarkers):
        local = mask & df[biomarker].notna()
        y = outcome["y"][local].to_numpy(dtype=int)
        x = pd.to_numeric(df.loc[local, biomarker], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(x)
        y, x = y[valid], x[valid]
        positive_n = int(y.sum())
        negative_n = int(len(y) - positive_n)
        if min(positive_n, negative_n) < 5:
            raise ValueError(
                f"{outcome['analysis']}: {biomarker} has {positive_n} positive and "
                f"{negative_n} negative participants with finite values; "
                "5-fold cross-validation requires at least 5 in each class."
            )
        apparent = optimal_metrics(y, x)
        ci_low, ci_high = bootstrap_ci(
            y, x, apparent["direction_sign"], 20261002 + outcome_idx * 100 + biomarker_idx
        )
        cv = cv_metrics(y, x, 20262002 + outcome_idx * 100 + biomarker_idx)
        short = LABELS[biomarker]
        fpr, tpr = roc_points(y, x * apparent["direction_sign"])
        roc_data[outcome["analysis"]][short] = {"fpr": fpr, "tpr": tpr}
        rows.append({
            "Analysis": outcome["analysis"],
            "N": int(len(y)),
            "Positive n": int(y.sum()),
            "Negative n": int(len(y) - y.sum()),
            "Biomarker": short,
            "AUC": apparent["auc"],
            "AUC 95% CI lower": ci_low,
            "AUC 95% CI upper": ci_high,
            "CV accuracy": cv["cv_accuracy"],
            "CV balanced accuracy": cv["cv_balanced_accuracy"],
            "CV sensitivity": cv["cv_sensitivity"],
            "CV specificity": cv["cv_specificity"],
            "Cutoff direction": apparent["direction"],
            "Cutoff": apparent["cutoff"],
        })

results = pd.DataFrame(rows)
results["rank"] = results.groupby("Analysis")["AUC"].rank(method="first", ascending=False)
results = results.sort_values(["Analysis", "rank"]).drop(columns="rank")
results.to_csv(OUTPUT_DIR / "biomarker_performance_results.csv", index=False)
(OUTPUT_DIR / "biomarker_performance_results.json").write_text(
    results.to_json(orient="records", force_ascii=False, indent=2), encoding="utf-8"
)
metadata = {
    "outcomes": [
        {key: value for key, value in outcome.items() if key not in {"mask", "y"}}
        for outcome in outcomes
    ],
    "bootstrap_replicates": 3000,
    "cv_folds": 5,
    "cv_repeats": 50,
}
(OUTPUT_DIR / "biomarker_performance_metadata.json").write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
)


def save_figure(fig, stem):
    fig.savefig(OUTPUT_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


# Manuscript Figure 2 (main figure): six-panel ROC comparison. The biological
# narrative remains A -> T; within T, the primary tau PET-only analysis
# precedes the operational PET/CSF analysis, and the A+-conditioned
# tau-PET-only analysis (primary within A+) closes out the grid.
plot_outcomes = [
    next(outcome for outcome in outcomes if outcome["analysis"] == analysis)
    for analysis in PLOT_ANALYSIS_ORDER
]
fig = plt.figure(figsize=(9.6, 5.7))
grid = fig.add_gridspec(2, 3, left=0.065, right=0.99, bottom=0.09, top=0.91, hspace=0.50, wspace=0.25)
axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in range(3)]
panel_letters = ["a", "b", "c", "d", "e", "f"]
for ax, outcome, letter in zip(axes, plot_outcomes, panel_letters):
    subset = results[results["Analysis"] == outcome["analysis"]].set_index("Biomarker")
    for marker in ORDER:
        curve = roc_data[outcome["analysis"]][marker]
        auc = subset.loc[marker, "AUC"]
        ax.plot(
            curve["fpr"], curve["tpr"],
            color=COLORS[marker], linestyle=LINESTYLES[marker],
            linewidth=1.6 if marker in {"MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "MTBR-Tau243", "p-Tau217"} else 1.15,
            alpha=0.98, label=f"{DISPLAY_OVERRIDE.get(marker, marker)}  {auc:.3f}",
        )
    ax.plot([0, 1], [0, 1], color="#B8B8B8", linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    # In-image panel title/subtitle removed at author direction (2026-11-01); the
    # population and n for each lettered panel are described in the Figure 3 legend
    # in manuscript/full_manuscript_v1.md.
    ax.text(-0.17, 1.15, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    ax.legend(title="Biomarker   AUC", title_fontsize=5.2, fontsize=4.75, loc="lower right", handlelength=2.3, labelspacing=0.20)
    ax.tick_params(labelsize=6.0, length=3)
axes[0].set_ylabel("Sensitivity", fontsize=7)
axes[3].set_ylabel("Sensitivity", fontsize=7)
axes[3].set_xlabel("1 − specificity", fontsize=7)
axes[4].set_xlabel("1 − specificity", fontsize=7)
axes[5].set_xlabel("1 − specificity", fontsize=7)
save_figure(fig, "Fig3_full_ROC_11marker")  # main Figure 3 since 2026-11-05 (was Figure 4)


# Supplementary Figure S2 (renumbered 2026-11-01 from S5, to match true first-citation
# order -- see changelog 3.11): matching six-panel AUC forest summary.
fig = plt.figure(figsize=(9.6, 5.9))
grid = fig.add_gridspec(2, 3, left=0.12, right=0.99, bottom=0.09, top=0.91, hspace=0.48, wspace=0.32)
axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in range(3)]
for ax, outcome, letter in zip(axes, plot_outcomes, panel_letters):
    subset = results[results["Analysis"] == outcome["analysis"]].set_index("Biomarker").loc[ORDER]
    y_positions = np.arange(len(ORDER))[::-1]
    for y_pos, marker in zip(y_positions, ORDER):
        row = subset.loc[marker]
        ax.errorbar(
            row["AUC"], y_pos,
            xerr=[[row["AUC"] - row["AUC 95% CI lower"]], [row["AUC 95% CI upper"] - row["AUC"]]],
            fmt="o", markersize=4.3 if marker in {"MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "MTBR-Tau243", "p-Tau217"} else 3.7,
            color=COLORS[marker], ecolor=COLORS[marker], elinewidth=1.1, capsize=2.2,
            markeredgecolor="white", markeredgewidth=0.45,
        )
    ax.axvline(0.5, color="#A8A8A8", linewidth=0.8, linestyle=(0, (3, 3)))
    ax.set_xlim(0.35, 1.01)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_yticks(y_positions)
    # Biomarker order is identical within each row.  Showing the category
    # labels only on the left-most panel prevents long labels from colliding
    # with the adjacent confidence intervals in the three-panel lower row.
    ax.set_yticklabels(
        [DISPLAY_OVERRIDE.get(m, m) for m in ORDER] if letter in {"a", "d"} else [],
        fontsize=6.4,
    )
    ax.set_xlabel("AUC (95% CI)", fontsize=7)
    # In-image panel title/subtitle removed at author direction (2026-11-01); the
    # population and n for each lettered panel match Figure 4 (same panel order),
    # to be described in the forthcoming supplementary figure legend for S2.
    ax.text(-0.18, 1.15, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    ax.tick_params(axis="x", labelsize=6.2, length=3)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
save_figure(fig, "FigS2_full_AUC_forest_11marker")

print(results.to_string(index=False))
