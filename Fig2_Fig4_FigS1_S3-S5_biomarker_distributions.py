import json
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl
mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("/Users/qiuyuyue/Documents/My-Papers/Blood biomarkers/biomarker_result.xlsx")
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(SCRIPT_DIR.parent / "outputs")))
RESULTS_FILE = OUTPUT_DIR / "biomarker_performance_results.csv"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
MAIN_REVISION = os.environ.get("DISTRIBUTION_MAIN_REVISION") == "1"
MAIN_ANALYSES = {"A status", "T status, Tau PET subset", "T status"}

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
})

NEGATIVE_COLOR = "#3775BA"
POSITIVE_COLOR = "#B64342"

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

DISPLAY_NAME = {
    "MTBR-Tau243": "eMTBR-tau243", "MTBR-Tau243/Aβ42": "eMTBR-tau243/Aβ42",
    "p-Tau217": "p-tau217", "p-Tau217/Aβ42": "p-tau217/Aβ42",
    "p-Tau231": "p-tau231", "p-Tau231/Aβ42": "p-tau231/Aβ42",
    "BD Tau": "BD-tau", "BD-Tau/Aβ42": "BD-tau/Aβ42", "NFL": "NfL",
}


def axis_label(marker, log_scale):
    """Display name plus unit (pg/mL for absolute concentrations; unitless ratio)."""
    name = DISPLAY_NAME.get(marker, marker)
    unit = "ratio" if "/" in marker else "pg/mL"
    return f"{name} ({unit}, log)" if log_scale else f"{name} ({unit})"


PANEL_ORDER = [
    "MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "BD-Tau/Aβ42", "p-Tau231/Aβ42",
    "MTBR-Tau243", "p-Tau217", "BD Tau", "p-Tau231",
    "GFAP", "Aβ1-42", "NFL",
]
LOG_SCALE_MARKERS = set(PANEL_ORDER)


def safe_values(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return values


def kde_density(values, x):
    """Gaussian-kernel density estimate without assuming a normal population."""
    values = safe_values(values)
    n = len(values)
    sd = float(np.std(values, ddof=1)) if n > 1 else 0.0
    q25, q75 = np.quantile(values, [0.25, 0.75])
    robust_sd = float((q75 - q25) / 1.349)
    scale_candidates = [v for v in (sd, robust_sd) if np.isfinite(v) and v > 0]
    scale = min(scale_candidates) if scale_candidates else max(abs(float(np.mean(values))) * 0.01, 1e-6)
    bandwidth = max(0.9 * scale * n ** (-1 / 5), 1e-9)
    z = (x[:, None] - values[None, :]) / bandwidth
    density = np.mean(np.exp(-0.5 * z ** 2), axis=1) / (bandwidth * np.sqrt(2 * np.pi))
    return density


def original_value_on_log_axis(value, _position):
    raw = 10 ** value
    if raw >= 100:
        return f"{raw:.0f}"
    if raw >= 10:
        return f"{raw:.1f}".rstrip("0").rstrip(".")
    if raw >= 1:
        return f"{raw:.2f}".rstrip("0").rstrip(".")
    return f"{raw:.3f}".rstrip("0").rstrip(".")


def save_figure(fig, stem):
    fig.savefig(OUTPUT_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(
        OUTPUT_DIR / f"{stem}.tiff",
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


df = pd.read_excel(INPUT, sheet_name="Sheet1")
base_columns = [
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
biomarker_columns = {LABELS[column]: column for column in base_columns + list(ratio_numerators)}
if not RESULTS_FILE.exists():
    raise FileNotFoundError(
        f"Missing {RESULTS_FILE}. Run Fig3_FigS2_TableS4_ROC.py first to generate the locked "
        "six-analysis performance and Youden-cutoff results."
    )
performance = pd.read_csv(RESULTS_FILE)
cutoff_lookup = {
    (row["Analysis"], row["Biomarker"]): (row["Cutoff direction"], float(row["Cutoff"]))
    for _, row in performance.iterrows()
}

analyses = [
    {
        "analysis": "A status",
        "stem": "FigS1_distribution_A_all",  # secondary A population (n=135) since v3.28
        "title": "Biomarker distributions for A status",
        "subtitle": "A+ versus A−; all participants",
        "mask": df["A"].isin(["A+", "A-"]),
        "group": df["A"],
        "negative": "A-",
        "positive": "A+",
        "negative_label": "A−",
        "positive_label": "A+",
    },
    {
        "analysis": "A status, non-control cohort",
        "stem": "Fig2_distribution_A_noncontrol",  # primary A population (n=101) since v3.28
        "title": "Biomarker distributions for A status",
        "subtitle": "A+ versus A−; healthy controls excluded",
        "mask": df["A"].isin(["A+", "A-"]) & df["DIAG"].ne("健康对照"),
        "group": df["A"],
        "negative": "A-",
        "positive": "A+",
        "negative_label": "A−",
        "positive_label": "A+",
    },
    {
        "analysis": "T status, Tau PET subset",
        "stem": "Fig4_distribution_T_tauPET",  # main Figure 4 since 2026-11-05 (was Figure 3)
        "title": "Biomarker distributions for T status",
        "subtitle": "T+ versus T−; tau PET-defined T",
        "mask": df["T"].isin(["T+", "T-"]) & df["taupet"].eq(1),
        "group": df["T"],
        "negative": "T-",
        "positive": "T+",
        "negative_label": "T−",
        "positive_label": "T+",
        "expected": (68, 34, 34),  # corrected 2026-11-01: locked cohort is n=68 (34 T+/34 T-), not the earlier n=67 snapshot
    },
    {
        "analysis": "T status",
        "stem": "FigS3_distribution_T_PET_CSF",
        "title": "Biomarker distributions for T status",
        "subtitle": "T+ versus T−; tau PET/CSF-defined T",
        "mask": df["T"].isin(["T+", "T-"]),
        "group": df["T"],
        "negative": "T-",
        "positive": "T+",
        "negative_label": "T−",
        "positive_label": "T+",
        "expected": (114, 80, 34),
    },
    {
        "analysis": "T status, A+ subset, Tau PET only",
        "stem": "FigS4_distribution_T_within_Aplus_tauPET",
        "title": "Biomarker distributions for T status",
        "subtitle": "T+ versus T−; A+ participants with tau PET-defined T",
        "mask": df["A"].eq("A+") & df["T"].isin(["T+", "T-"]) & df["taupet"].eq(1),
        "group": df["T"],
        "negative": "T-",
        "positive": "T+",
        "negative_label": "T−",
        "positive_label": "T+",
        "expected": (52, 18, 34),
    },
    {
        "analysis": "T status, A+ subset",
        "stem": "FigS5_distribution_T_within_Aplus_PET_CSF",
        "title": "Biomarker distributions for T status",
        "subtitle": "T+ versus T−; A+ participants with tau PET/CSF-defined T",
        "mask": df["A"].eq("A+") & df["T"].isin(["T+", "T-"]),
        "group": df["T"],
        "negative": "T-",
        "positive": "T+",
        "negative_label": "T−",
        "positive_label": "T+",
        "expected": (67, 33, 34),
    },
]

# Locked expected sizes for the two A-status analyses.
analyses[0]["expected"] = (135, 68, 67)
analyses[1]["expected"] = (101, 34, 67)

source_rows = []

for analysis in analyses:
    mask = analysis["mask"]
    group = analysis["group"]
    neg_n = int(np.sum(mask & group.eq(analysis["negative"])))
    pos_n = int(np.sum(mask & group.eq(analysis["positive"])))
    expected_total, expected_neg, expected_pos = analysis["expected"]
    observed = (neg_n + pos_n, neg_n, pos_n)
    if observed != (expected_total, expected_neg, expected_pos):
        raise AssertionError(
            f"{analysis['analysis']} population drift: observed total/negative/positive "
            f"{observed}, expected {(expected_total, expected_neg, expected_pos)}"
        )

    # Main-text figures (stem "Fig...") omit sample-derived Youden cutoffs (Methods 2.5);
    # supplementary figures (stem "FigS...") show them. Source data always include cutoffs.
    show_cutoff = analysis["stem"].startswith("FigS")
    fig, axes_grid = plt.subplots(3, 4, figsize=(7.5, 7.0))
    axes = axes_grid.ravel()

    for panel_index, marker in enumerate(PANEL_ORDER):
        ax = axes[panel_index]
        column = biomarker_columns[marker]
        neg_values = safe_values(
            pd.to_numeric(df.loc[mask & group.eq(analysis["negative"]), column], errors="coerce")
        )
        pos_values = safe_values(
            pd.to_numeric(df.loc[mask & group.eq(analysis["positive"]), column], errors="coerce")
        )

        use_log_scale = marker in LOG_SCALE_MARKERS
        if use_log_scale:
            neg_plot_values = np.log10(neg_values)
            pos_plot_values = np.log10(pos_values)
        else:
            neg_plot_values = neg_values
            pos_plot_values = pos_values

        combined = np.concatenate([neg_plot_values, pos_plot_values])
        data_low = float(np.min(combined))
        data_high = float(np.max(combined))
        data_span = max(data_high - data_low, max(abs(data_low), abs(data_high), 1.0) * 0.02)
        x_low = data_low - 0.08 * data_span
        x_high = data_high + 0.08 * data_span
        if not use_log_scale and data_low >= 0:
            x_low = max(0.0, x_low)
        x = np.linspace(x_low, x_high, 600)

        neg_density = kde_density(neg_plot_values, x)
        pos_density = kde_density(pos_plot_values, x)
        ax.fill_between(x, 0, neg_density, color=NEGATIVE_COLOR, alpha=0.15, linewidth=0)
        ax.fill_between(x, 0, pos_density, color=POSITIVE_COLOR, alpha=0.15, linewidth=0)
        ax.plot(x, neg_density, color=NEGATIVE_COLOR, linewidth=1.7)
        ax.plot(x, pos_density, color=POSITIVE_COLOR, linewidth=1.7)
        neg_median = float(np.median(neg_values))
        pos_median = float(np.median(pos_values))
        neg_median_plot = np.log10(neg_median) if use_log_scale else neg_median
        pos_median_plot = np.log10(pos_median) if use_log_scale else pos_median
        ax.axvline(neg_median_plot, color=NEGATIVE_COLOR, linewidth=0.8, linestyle=(0, (3, 2)), alpha=0.85)
        ax.axvline(pos_median_plot, color=POSITIVE_COLOR, linewidth=0.8, linestyle=(0, (3, 2)), alpha=0.85)
        cutoff_direction, cutoff = cutoff_lookup[(analysis["analysis"], marker)]
        if show_cutoff:
            cutoff_plot = np.log10(cutoff) if use_log_scale else cutoff
            ax.axvline(cutoff_plot, color="#222222", linewidth=0.95, linestyle=(0, (1.5, 1.5)), alpha=0.95)
        ymax = max(float(np.max(neg_density)), float(np.max(pos_density)))
        ax.vlines(neg_plot_values, -0.055 * ymax, -0.015 * ymax, color=NEGATIVE_COLOR, linewidth=0.35, alpha=0.38)
        ax.vlines(pos_plot_values, -0.115 * ymax, -0.075 * ymax, color=POSITIVE_COLOR, linewidth=0.35, alpha=0.38)
        ax.set_ylim(-0.13 * ymax, 1.08 * ymax)

        # Presentation-only: always show the manufacturer-confirmed analyte name
        # eMTBR-tau243 (Methods 2.3), regardless of MAIN_REVISION (which only
        # controls whether sample-derived cutoff lines/source-data export are shown).
        ax.set_xlabel(axis_label(marker, use_log_scale), fontsize=6.6, fontweight="bold")
        if panel_index in {0, 4, 8}:
            ax.set_ylabel("Kernel density", fontsize=6.5)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=6, length=3)
        if use_log_scale:
            log_ticks = np.arange(np.ceil(x_low), np.floor(x_high) + 1, dtype=float)
            ax.set_xticks(log_ticks)
            ax.xaxis.set_major_formatter(FuncFormatter(original_value_on_log_axis))
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        if show_cutoff:
            ax.text(
                0.03, 0.95,
                f"Cutoff {cutoff_direction} {cutoff:.3f}",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=5.4, color="#222222",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.6},
            )
        ax.text(
            0.98, 0.875,
            f"{analysis['negative_label']}: median={neg_median:.3g}",
            transform=ax.transAxes,
            ha="right", va="top", fontsize=5.3, color=NEGATIVE_COLOR,
        )
        ax.text(
            0.98, 0.80,
            f"{analysis['positive_label']}: median={pos_median:.3g}",
            transform=ax.transAxes,
            ha="right", va="top", fontsize=5.3, color=POSITIVE_COLOR,
        )

        if True:  # source data always exported (includes cutoffs for every analysis)
            for group_label, values in (
                (analysis["negative_label"], neg_values),
                (analysis["positive_label"], pos_values),
            ):
                group_median = float(np.median(values))
                for observation_index, value in enumerate(values, start=1):
                    source_rows.append({
                        "Analysis": analysis["analysis"],
                        "Figure stem": analysis["stem"],
                        "Biomarker": marker,
                        "Group": group_label,
                        "Observation within group": observation_index,
                        "Biomarker value": float(value),
                        "log10 biomarker value": float(np.log10(value)),
                        "Group median": group_median,
                        "Youden cutoff direction": cutoff_direction,
                        "Youden cutoff": float(cutoff),
                    })

    legend_ax = axes[11]
    legend_ax.set_axis_off()
    handles = [
        Line2D([0], [0], color=NEGATIVE_COLOR, lw=2.0, label=f"{analysis['negative_label']} (n={neg_n})"),
        Line2D([0], [0], color=POSITIVE_COLOR, lw=2.0, label=f"{analysis['positive_label']} (n={pos_n})"),
        Line2D([0], [0], color="#555555", lw=0.8, linestyle=(0, (3, 2)), label="Group median"),
    ]
    if show_cutoff:
        handles.append(Line2D([0], [0], color="#222222", lw=0.95, linestyle=(0, (1.5, 1.5)), label="Youden cutoff"))
    legend_ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.96),
        fontsize=7.5, handlelength=3.2, labelspacing=0.45,
    )
    # In-image caption/title text removed at author direction (2026-11-01); the full
    # methodological description and panel content live in the figure legend
    # (manuscript/full_manuscript_v1.md for Figures 2/4; a supplementary legends file for
    # the S1/S3/S4/S5 variants this same function also produces).
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.065, top=0.89, wspace=0.28, hspace=0.47)
    save_figure(fig, analysis["stem"])

pd.DataFrame(source_rows).to_csv(
    OUTPUT_DIR / "FigS_distributions_source_data.csv", index=False
)
(OUTPUT_DIR / "FigS_distributions_metadata.json").write_text(
    json.dumps(
        {
            "input": INPUT.name,
            "performance_input": RESULTS_FILE.name,
            "backend": "Python/matplotlib",
            "analyses": [
                {
                    "analysis": item["analysis"],
                    "figure_stem": item["stem"],
                    "subtitle": item["subtitle"],
                    "expected_total_negative_positive": list(item["expected"]),
                }
                for item in analyses
            ],
            "density_method": "Gaussian KDE on log10-transformed biomarker values",
            "cutoff_method": "Full-analysis-sample Youden index from Fig3_FigS2_TableS4_ROC.py",
            "exports": ["SVG", "PDF", "PNG 300 dpi", "TIFF 600 dpi LZW"],
            "interpretive_guardrail": (
                "KDE curves are descriptive; modality claims require the separate formal tests."
            ),
        },
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)

print("Created distribution figures: main Figures 2 (A, controls excluded, primary) and 4 (T, tau PET); supplementary S1 (A, all participants), S3, S4, S5:")
for analysis in analyses:
    print(analysis["stem"])
