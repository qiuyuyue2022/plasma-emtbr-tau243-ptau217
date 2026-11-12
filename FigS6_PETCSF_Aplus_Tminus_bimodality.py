"""Create the integrated A+T- eMTBR-tau243 heterogeneity figure.

Scientific claim
----------------
The eMTBR-tau243-defined groups delineate a cross-sectional gradient of
AD-related molecular abnormality within A+T-. The full-group mixture is the
primary distributional evidence, recurrence in tau PET-confirmed T- is the key
robustness evidence, non-defining plasma biomarkers provide the main molecular
support, and contemporaneous CSF profiles provide exploratory support only.

This script does not redefine subgroups or perform new hypothesis tests. It
draws established results from the locked analysis outputs and checks that a
deterministic refit of the full-group GMM reproduces the locked diagnostics.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
# SI (v3.62, 2026-11-11): combined tau PET/CSF-defined A+T− analysis (n=33), restored from the
# pre-v3.60 code for the Supplementary Information. Statistics go to outputs/Aplus_Tminus_PET_CSF/.
ROOT_OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", PROJECT_DIR / "outputs"))
OUTPUT_DIR = ROOT_OUTPUT_DIR / "Aplus_Tminus_PET_CSF"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl

mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
import sklearn
from sklearn.mixture import GaussianMixture


# Mandatory publication/export settings: editable SVG/PDF text and compact
# double-column typography.
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "svg.hashsalt": "Aplus_Tminus_molecular_gradient_20261030",
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
NEUTRAL = "#8C8C8C"
TEXT = "#2B2B2B"

MEMBERSHIP_PATH = OUTPUT_DIR / "A+T-Subgroup_membership.csv"
FULL_SUMMARY_PATH = OUTPUT_DIR / "bimodality_validation_summary.csv"
REFERENCE_SUMMARY_PATH = OUTPUT_DIR / "GMM_reference_stratified_summary.csv"
REFERENCE_DENSITY_PATH = OUTPUT_DIR / "GMM_reference_stratified_figure_source_data.csv"
REFERENCE_ASSIGNMENTS_PATH = OUTPUT_DIR / "GMM_reference_stratified_assignments.csv"
STEM = ROOT_OUTPUT_DIR / "FigS6_PETCSF_Aplus_Tminus_bimodality"


def require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise KeyError(f"{name} is missing required columns: {missing}")


def fit_full_gmm(values: np.ndarray) -> dict[str, object]:
    one = GaussianMixture(
        n_components=1,
        covariance_type="full",
        n_init=100,
        max_iter=1000,
        reg_covar=1e-6,
        random_state=20261026,
    ).fit(values.reshape(-1, 1))
    two = GaussianMixture(
        n_components=2,
        covariance_type="full",
        n_init=100,
        max_iter=1000,
        reg_covar=1e-6,
        random_state=20261027,
    ).fit(values.reshape(-1, 1))
    order = np.argsort(two.means_.ravel())
    low_component, high_component = int(order[0]), int(order[1])
    x = values.reshape(-1, 1)
    return {
        "one": one,
        "two": two,
        "low_component": low_component,
        "high_component": high_component,
        "delta_bic": float(one.bic(x) - two.bic(x)),
    }


def normal_density(grid: np.ndarray, mean: float, variance: float) -> np.ndarray:
    return norm.pdf(grid, loc=mean, scale=np.sqrt(variance))


def build_full_density(membership: pd.DataFrame, full_summary: pd.DataFrame) -> pd.DataFrame:
    values = np.log10(membership["MTBR-Tau243"].to_numpy(dtype=float))
    fitted = fit_full_gmm(values)
    locked = full_summary.loc[full_summary["Marker"].eq("MTBR-Tau243")].iloc[0]
    if not np.isclose(fitted["delta_bic"], locked["Observed delta BIC"], atol=1e-6):
        raise AssertionError(
            f"Full-group GMM refit mismatch: {fitted['delta_bic']} vs "
            f"{locked['Observed delta BIC']}"
        )

    grid = np.linspace(values.min() - 0.20, values.max() + 0.20, 700)
    one = fitted["one"]
    two = fitted["two"]
    component_density = []
    for component in range(2):
        component_density.append(
            two.weights_.ravel()[component]
            * normal_density(
                grid,
                two.means_.ravel()[component],
                two.covariances_.reshape(-1)[component],
            )
        )
    return pd.DataFrame(
        {
            "log10 eMTBR-tau243": grid,
            "eMTBR-tau243": 10**grid,
            "One-component density": normal_density(
                grid,
                one.means_.ravel()[0],
                one.covariances_.reshape(-1)[0],
            ),
            "Two-component total density": component_density[0] + component_density[1],
            "Low-component weighted density": component_density[fitted["low_component"]],
            "High-component weighted density": component_density[fitted["high_component"]],
        }
    )


def original_scale_ticks(x_min: float, x_max: float) -> tuple[np.ndarray, list[str]]:
    candidates = np.array([0.1, 0.3, 1, 3, 10, 30, 100], dtype=float)
    log_candidates = np.log10(candidates)
    keep = (log_candidates >= x_min) & (log_candidates <= x_max)
    return log_candidates[keep], [f"{value:g}" for value in candidates[keep]]


def add_panel_label(axis: plt.Axes, label: str, x: float = -0.13, y: float = 1.06) -> None:
    axis.text(
        x,
        y,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="top",
    )


def plot_full_gmm(
    axis: plt.Axes,
    density: pd.DataFrame,
    membership: pd.DataFrame,
    locked: pd.Series,
) -> None:
    x = density["log10 eMTBR-tau243"].to_numpy()
    axis.fill_between(
        x,
        density["Low-component weighted density"],
        color=LOW_COLOR,
        alpha=0.18,
        linewidth=0,
    )
    axis.fill_between(
        x,
        density["High-component weighted density"],
        color=HIGH_COLOR,
        alpha=0.15,
        linewidth=0,
    )
    axis.plot(
        x,
        density["One-component density"],
        color=NEUTRAL,
        linewidth=1.1,
        linestyle=(0, (4, 2)),
        label="One-component fit",
    )
    axis.plot(
        x,
        density["Two-component total density"],
        color="#222222",
        linewidth=1.8,
        label="Two-component fit",
    )
    boundary = float(locked["Observed decision boundary"])
    axis.axvline(np.log10(boundary), color="#222222", linestyle=":", linewidth=1.0)

    ymax = max(
        density["One-component density"].max(),
        density["Two-component total density"].max(),
    )
    rug_low = -0.04 * ymax
    for _, row in membership.sort_values("MTBR-Tau243").iterrows():
        color = HIGH_COLOR if str(row["Plasma tau subgroup"]).startswith("High") else LOW_COLOR
        value = np.log10(float(row["MTBR-Tau243"]))
        axis.plot([value, value], [rug_low, 0], color=color, linewidth=0.75, alpha=0.75, clip_on=False)

    ticks, ticklabels = original_scale_ticks(x.min(), x.max())
    axis.set_xticks(ticks)
    axis.set_xticklabels(ticklabels)
    axis.set_xlim(x.min(), x.max())
    axis.set_ylim(rug_low * 1.7, ymax * 1.15)
    axis.set_ylabel("Density")
    # Panel title removed at author direction (2026-11-01); described in the
    # supplementary figure legend for S6 instead.
    axis.text(
        0.02,
        0.95,
        f"ΔBIC = {locked['Observed delta BIC']:.2f}\n"
        "p < 0.001",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=6.3,
        color=TEXT,
    )
    add_panel_label(axis, "a")


def plot_reference_stratum(
    axis: plt.Axes,
    stratum: str,
    density: pd.DataFrame,
    assignments: pd.DataFrame,
    summary: pd.Series,
) -> None:
    local_density = density.loc[density["Reference stratum"].eq(stratum)]
    local_people = assignments.loc[assignments["Reference stratum"].eq(stratum)]
    x = local_density["log10 eMTBR-tau243"].to_numpy()
    axis.fill_between(
        x,
        local_density["Low-component weighted density"],
        color=LOW_COLOR,
        alpha=0.13,
        linewidth=0,
    )
    axis.fill_between(
        x,
        local_density["High-component weighted density"],
        color=HIGH_COLOR,
        alpha=0.11,
        linewidth=0,
    )
    axis.plot(
        x,
        local_density["One-component density"],
        color="#A0A0A0",
        linewidth=0.9,
        linestyle=(0, (4, 2)),
    )
    axis.plot(
        x,
        local_density["Two-component total density"],
        color=TEXT,
        linewidth=1.55,
    )
    axis.axvline(
        np.log10(float(summary["Observed decision boundary"])),
        color="#333333",
        linewidth=0.9,
        linestyle=":",
    )
    ymax = max(
        local_density["One-component density"].max(),
        local_density["Two-component total density"].max(),
    )
    rug_low = -0.045 * ymax
    for _, person in local_people.iterrows():
        rug_color = HIGH_COLOR if person["Observed subgroup"] == "High" else LOW_COLOR
        value = np.log10(float(person["eMTBR-tau243"]))
        axis.plot([value, value], [rug_low, 0], color=rug_color, linewidth=0.65, alpha=0.7, clip_on=False)
    ticks, ticklabels = original_scale_ticks(x.min(), x.max())
    axis.set_xticks(ticks)
    axis.set_xticklabels(ticklabels, fontsize=5.8)
    axis.set_xlim(x.min(), x.max())
    axis.set_ylim(rug_low * 1.8, ymax * 1.17)
    title = "Tau PET-confirmed T−" if stratum == "Tau PET-negative" else "CSF-defined T−"
    # Panel title removed at author direction (2026-11-01); described in the
    # supplementary figure legend for S6 instead.
    p_value = float(summary["Parametric-bootstrap LRT p"])
    if stratum == "Tau PET-negative":
        # 2026-11-09 revision: the tau PET-defined group (n=18) is now the primary analysis; show the
        # p value from its primary run (02a with BIOMARKER_APTM_REFERENCE=tauPET), which the text and
        # Supplementary Table S7 report. No BH-FDR q: the two strata are no longer a correction family.
        primary = pd.read_csv(ROOT_OUTPUT_DIR / "Aplus_Tminus_tauPET_only" / "bimodality_validation_summary.csv").set_index("Marker")
        p_value = float(primary.loc["MTBR-Tau243", "Parametric-bootstrap LRT p under 1-component null"])
    p_text = f"p={p_value:.3f}"
    axis.text(
        0.03,
        0.95,
        f"ΔBIC {summary['Observed delta BIC']:.2f}\n{p_text}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=5.9,
        color=TEXT,
    )
    axis.tick_params(axis="y", labelsize=5.8, length=2.5)


def save_figure(figure: plt.Figure) -> None:
    figure.savefig(
        STEM.with_suffix(".svg"),
        bbox_inches="tight",
        metadata={"Date": "2026-10-30"},
    )
    figure.savefig(
        STEM.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    figure.savefig(
        STEM.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    figure.savefig(STEM.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    membership = pd.read_csv(MEMBERSHIP_PATH)
    full_summary = pd.read_csv(FULL_SUMMARY_PATH)
    reference_summary = pd.read_csv(REFERENCE_SUMMARY_PATH)
    reference_density = pd.read_csv(REFERENCE_DENSITY_PATH)
    reference_assignments = pd.read_csv(REFERENCE_ASSIGNMENTS_PATH)

    require_columns(
        membership,
        {"SID", "Plasma tau subgroup", "MTBR-Tau243"},
        "Subgroup membership",
    )
    if len(membership) != 33:
        raise AssertionError(f"Expected 33 A+T- participants, found {len(membership)}.")
    if membership["Plasma tau subgroup"].value_counts().to_dict() != {
        "High MTBR-Tau243": 21,
        "Low MTBR-Tau243": 12,
    }:
        raise AssertionError("Locked full-group subgroup counts changed.")

    full_density = build_full_density(membership, full_summary)
    locked_full = full_summary.loc[full_summary["Marker"].eq("MTBR-Tau243")].iloc[0]

    figure = plt.figure(figsize=(7.2, 2.75))
    outer = figure.add_gridspec(
        1,
        2,
        left=0.085,
        right=0.985,
        bottom=0.265,
        top=0.895,
        width_ratios=[1.04, 1.0],
        wspace=0.36,
    )

    axis_a = figure.add_subplot(outer[0, 0])
    plot_full_gmm(axis_a, full_density, membership, locked_full)

    reference_grid = outer[0, 1].subgridspec(1, 2, wspace=0.30)
    axis_b1 = figure.add_subplot(reference_grid[0, 0])
    axis_b2 = figure.add_subplot(reference_grid[0, 1], sharey=axis_b1)
    for axis, stratum in [
        (axis_b1, "Tau PET-negative"),
        (axis_b2, "CSF-defined T-negative"),
    ]:
        row = reference_summary.loc[reference_summary["Reference stratum"].eq(stratum)].iloc[0]
        plot_reference_stratum(
            axis,
            stratum,
            reference_density,
            reference_assignments,
            row,
        )
    axis_b1.set_ylabel("Density", fontsize=6.2)
    axis_b2.tick_params(labelleft=False)
    figure.text(0.545, 0.955, "b", fontsize=9, fontweight="bold", ha="left", va="top")
    figure.supxlabel("eMTBR-tau243 (pg/mL, log scale)", x=0.53, y=0.145, fontsize=7)

    legend_handles = [
        Line2D([0], [0], color=NEUTRAL, linewidth=1.1, linestyle=(0, (4, 2)), label="One-component fit"),
        Line2D([0], [0], color="#222222", linewidth=1.8, label="Two-component fit"),
        Line2D([0], [0], color="#222222", linewidth=1.0, linestyle=":", label="Component boundary"),
        Line2D([0], [0], marker="|", color=LOW_COLOR, linestyle="None", markersize=8, label="Low component"),
        Line2D([0], [0], marker="|", color=HIGH_COLOR, linestyle="None", markersize=8, label="High component"),
    ]
    figure.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=5,
        fontsize=6.0,
        handlelength=2.4,
        columnspacing=1.2,
    )
    # In-image suptitle/caption text removed at author direction (2026-11-01); the
    # full methodological description (effect-size convention, FDR families, CSF
    # exploratory caveat, reference-standard imbalance) lives in the supplementary
    # figure legend for S6.

    save_figure(figure)

    full_density.to_csv(
        OUTPUT_DIR / "FigS6_Aplus_Tminus_heterogeneity_full_GMM_source.csv",
        index=False,
    )
    metadata = {
        "core_conclusion": (
            "eMTBR-tau243-defined groups delineate a cross-sectional gradient of "
            "AD-related molecular abnormality within A+T-."
        ),
        "archetype": "asymmetric mixed-modality quantitative figure",
        "backend": "Python/matplotlib",
        "final_size_inches": [7.2, 2.75],
        "source_files": [
            str(MEMBERSHIP_PATH),
            str(FULL_SUMMARY_PATH),
            str(REFERENCE_SUMMARY_PATH),
            str(REFERENCE_DENSITY_PATH),
            str(REFERENCE_ASSIGNMENTS_PATH),
        ],
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "review_risk": (
            "Cross-sectional exploratory subgroups; no subtype, stage, or progression "
            "inference."
        ),
    }
    (OUTPUT_DIR / "FigS6_Aplus_Tminus_heterogeneity_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Figure and source data written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
