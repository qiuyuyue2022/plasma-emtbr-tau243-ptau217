import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(SCRIPT_DIR / ".mplconfig"),
)
os.environ.setdefault(
    "XDG_CACHE_HOME",
    str(SCRIPT_DIR / ".cache"),
)

import matplotlib as mpl

mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


DEFAULT_INPUT = Path("/Users/qiuyuyue/Documents/My-Papers/Blood biomarkers/biomarker_result.xlsx")
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
STEM = "Fig5_AT_3group_distributions"

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
    }
)

DISPLAY_TITLE = {
    "MTBR-Tau243": "eMTBR-tau243", "MTBR-Tau243/Aβ42": "eMTBR-tau243/Aβ42",
    "p-Tau217": "p-tau217", "p-Tau217/Aβ42": "p-tau217/Aβ42",
    "p-Tau231": "p-tau231", "p-Tau231/Aβ42": "p-tau231/Aβ42",
    "BD Tau": "BD-tau", "BD-Tau/Aβ42": "BD-tau/Aβ42", "NFL": "NfL",
}

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

PANEL_ORDER = [
    "MTBR-Tau243/Aβ42",
    "p-Tau217/Aβ42",
    "BD-Tau/Aβ42",
    "p-Tau231/Aβ42",
    "MTBR-Tau243",
    "p-Tau217",
    "BD Tau",
    "p-Tau231",
    "GFAP",
    "Aβ1-42",
    "NFL",
]

GROUPS = [
    {
        "key": "A−T−",
        "mask": lambda data: data["A"].eq("A-") & data["T"].eq("T-"),
        "color": "#4C78A8",
        "linestyle": "-",
    },
    {
        "key": "A+T−",
        "mask": lambda data: data["A"].eq("A+") & data["T"].eq("T-"),
        "color": "#E69F00",
        "linestyle": (0, (5, 2)),
    },
    {
        "key": "A+T+",
        "mask": lambda data: data["A"].eq("A+") & data["T"].eq("T+"),
        "color": "#C23B53",
        "linestyle": (0, (5, 1.5, 1.2, 1.5)),
    },
]


def safe_values(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if np.any(values <= 0):
        raise ValueError("All biomarkers must be strictly positive for log10 plotting.")
    return values


def pooled_bandwidth(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = len(values)
    sd = float(np.std(values, ddof=1)) if n > 1 else 0.0
    q25, q75 = np.quantile(values, [0.25, 0.75])
    robust_sd = float((q75 - q25) / 1.349)
    scales = [v for v in (sd, robust_sd) if np.isfinite(v) and v > 0]
    scale = min(scales) if scales else max(abs(float(np.mean(values))) * 0.01, 1e-6)
    return max(0.9 * scale * n ** (-1 / 5), 1e-9)


def kde_density(values, x, bandwidth):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    z = (x[:, None] - values[None, :]) / bandwidth
    return np.mean(np.exp(-0.5 * z**2), axis=1) / (bandwidth * np.sqrt(2 * np.pi))


def raw_tick_label(log_value, _position):
    raw = 10**log_value
    if raw >= 100:
        return f"{raw:.0f}"
    if raw >= 10:
        return f"{raw:.1f}".rstrip("0").rstrip(".")
    if raw >= 1:
        return f"{raw:.2f}".rstrip("0").rstrip(".")
    return f"{raw:.3f}".rstrip("0").rstrip(".")


def bh_adjust(p_values):
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 1.0
    count = len(p)
    for reverse_rank in range(count - 1, -1, -1):
        idx = order[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, p[idx] * count / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted


def format_q(value):
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def save_figure(fig):
    fig.savefig(OUTPUT_DIR / f"{STEM}.svg", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{STEM}.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{STEM}.png", dpi=300, bbox_inches="tight")
    fig.savefig(
        OUTPUT_DIR / f"{STEM}.tiff",
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )


df = pd.read_excel(INPUT, sheet_name="Sheet1")
abeta42 = pd.to_numeric(df["Aβ 1-42"], errors="coerce")
ratio_numerators = {
    "p-Tau217/Aβ42": "p-Tau217(≤0.758)",
    "p-Tau231/Aβ42": "p-Tau231(≤5.793)",
    "MTBR-Tau243/Aβ42": "MTBR-Tau243",
    "BD-Tau/Aβ42": "BD Tau（≤52.301）",
}
for ratio_name, numerator in ratio_numerators.items():
    df[ratio_name] = pd.to_numeric(df[numerator], errors="coerce") / abeta42

base_columns = [
    "GFAP(≤131.880)",
    "NFL(≤9.624 )",
    "Aβ 1-42",
    "p-Tau217(≤0.758)",
    "p-Tau231(≤5.793)",
    "BD Tau（≤52.301）",
    "MTBR-Tau243",
]
biomarker_columns = {LABELS[column]: column for column in base_columns + list(ratio_numerators)}
group_masks = {group["key"]: group["mask"](df) for group in GROUPS}
group_counts = {key: int(mask.sum()) for key, mask in group_masks.items()}

stats_rows = []
for marker in PANEL_ORDER:
    column = biomarker_columns[marker]
    values_by_group = {
        key: safe_values(pd.to_numeric(df.loc[mask, column], errors="coerce"))
        for key, mask in group_masks.items()
    }
    t_minus = values_by_group["A+T−"]
    t_plus = values_by_group["A+T+"]
    u_statistic, p_value = mannwhitneyu(t_plus, t_minus, alternative="two-sided", method="asymptotic")
    rank_biserial = 2.0 * u_statistic / (len(t_plus) * len(t_minus)) - 1.0
    stats_rows.append(
        {
            "Biomarker": marker,
            "A−T− n": len(values_by_group["A−T−"]),
            "A−T− median": float(np.median(values_by_group["A−T−"])),
            "A+T− n": len(t_minus),
            "A+T− median": float(np.median(t_minus)),
            "A+T+ n": len(t_plus),
            "A+T+ median": float(np.median(t_plus)),
            "Mann-Whitney U (A+T+ vs A+T−)": float(u_statistic),
            "Raw p": float(p_value),
            "Rank-biserial r (positive means higher in A+T+)": float(rank_biserial),
        }
    )

statistics = pd.DataFrame(stats_rows)
statistics["FDR q (11 biomarkers)"] = bh_adjust(statistics["Raw p"].to_numpy())
statistics.to_csv(OUTPUT_DIR / f"{STEM}_statistics.csv", index=False)
stats_by_marker = statistics.set_index("Biomarker")

fig, axes_grid = plt.subplots(3, 4, figsize=(7.6, 7.05))
axes = axes_grid.ravel()

for panel_index, marker in enumerate(PANEL_ORDER):
    ax = axes[panel_index]
    column = biomarker_columns[marker]
    raw_by_group = {
        key: safe_values(pd.to_numeric(df.loc[mask, column], errors="coerce"))
        for key, mask in group_masks.items()
    }
    log_by_group = {key: np.log10(values) for key, values in raw_by_group.items()}
    combined = np.concatenate(list(log_by_group.values()))
    data_low = float(np.min(combined))
    data_high = float(np.max(combined))
    span = max(data_high - data_low, max(abs(data_low), abs(data_high), 1.0) * 0.02)
    x_low = data_low - 0.08 * span
    x_high = data_high + 0.08 * span
    x = np.linspace(x_low, x_high, 650)
    bandwidth = pooled_bandwidth(combined)

    densities = {}
    for group in GROUPS:
        key = group["key"]
        density = kde_density(log_by_group[key], x, bandwidth)
        densities[key] = density
        ax.fill_between(x, 0, density, color=group["color"], alpha=0.075, linewidth=0)
        ax.plot(
            x,
            density,
            color=group["color"],
            linestyle=group["linestyle"],
            linewidth=1.65,
            label=key,
        )
        median_log = np.log10(np.median(raw_by_group[key]))
        ax.axvline(
            median_log,
            color=group["color"],
            linewidth=0.75,
            linestyle=(0, (2, 2)),
            alpha=0.8,
        )

    ymax = max(float(np.max(density)) for density in densities.values())
    rug_levels = [-0.032, -0.073, -0.114]
    for group, rug_level in zip(GROUPS, rug_levels):
        key = group["key"]
        ax.vlines(
            log_by_group[key],
            rug_level * ymax,
            (rug_level + 0.025) * ymax,
            color=group["color"],
            linewidth=0.35,
            alpha=0.42,
        )
    ax.set_ylim(-0.13 * ymax, 1.10 * ymax)

    row = stats_by_marker.loc[marker]
    ax.text(
        0.03,
        0.96,
        f"A+T− vs A+T+: q={format_q(row['FDR q (11 biomarkers)'])}, "
        f"r={row['Rank-biserial r (positive means higher in A+T+)']:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.2,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 0.6},
    )
    unit = "ratio" if "/" in marker else "pg/mL"
    ax.set_xlabel(f"{DISPLAY_TITLE.get(marker, marker)} ({unit}, log)", fontsize=6.6, fontweight="bold")
    if panel_index in {0, 4, 8}:
        ax.set_ylabel("Kernel density", fontsize=6.5)
    else:
        ax.set_yticklabels([])
    ax.tick_params(labelsize=6, length=3)
    log_ticks = np.arange(np.ceil(x_low), np.floor(x_high) + 1, dtype=float)
    ax.set_xticks(log_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(raw_tick_label))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))

legend_ax = axes[11]
legend_ax.set_axis_off()
handles = [
    Line2D(
        [0],
        [0],
        color=group["color"],
        linestyle=group["linestyle"],
        lw=2.0,
        label=f"{group['key']} (n={group_counts[group['key']]})",
    )
    for group in GROUPS
]
handles.append(
    Line2D([0], [0], color="#555555", lw=0.8, linestyle=(0, (2, 2)), label="Group median")
)
legend_ax.legend(
    handles=handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.96),
    fontsize=7.5,
    handlelength=3.2,
    labelspacing=0.5,
)
# In-image caption/title text removed at author direction (2026-11-01); the full
# methodological description, group n, and q/r annotation key live in the
# Figure 5 legend (main_figure_legends_v2.md).
fig.subplots_adjust(left=0.075, right=0.985, bottom=0.065, top=0.89, wspace=0.28, hspace=0.47)
save_figure(fig)
plt.close(fig)

print(statistics.to_string(index=False))
print(f"Created {STEM} in PNG, PDF, SVG, and TIFF formats.")
