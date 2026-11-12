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
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", PROJECT_DIR / "outputs"))
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

APTM_DIR = OUTPUT_DIR / "Aplus_Tminus_tauPET_only"
STEM = OUTPUT_DIR / "FigS5_Aplus_Tminus_heterogeneity"  # FigS6 before v3.60; FigS4 in v3.60


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
        f"p = {float(locked['Parametric-bootstrap LRT p under 1-component null']):.3f}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=6.3,
        color=TEXT,
    )
    # single-panel figure: no panel letter


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
    # v3.60 (2026-11-11): Supplementary Figure S5 (S4 in v3.60) shows the mixture-model fit within the
    # 18 tau PET-defined A+T− patients, from the outputs of TableS7_FigS6_GMM_full_sample.py
    # (summary) and TableS8_Aplus_Tminus_subgroup_comparison.py (membership) in
    # outputs/Aplus_Tminus_tauPET_only/. The pooled (n=33) and CSF-defined (n=15) panels
    # were removed; the GMM is refitted here only to draw the curves and is checked
    # against the locked ΔBIC.
    membership = pd.read_csv(APTM_DIR / "A+T-Subgroup_membership.csv")
    summary = pd.read_csv(APTM_DIR / "bimodality_validation_summary.csv")
    require_columns(membership, {"SID", "Plasma tau subgroup", "MTBR-Tau243"}, "Subgroup membership")
    if len(membership) != 18:
        raise AssertionError(f"Expected 18 tau PET-defined A+T− patients, found {len(membership)}.")
    density = build_full_density(membership, summary)
    locked = summary.loc[summary["Marker"].eq("MTBR-Tau243")].iloc[0]

    figure = plt.figure(figsize=(3.6, 2.9))
    axis = figure.add_axes([0.15, 0.32, 0.80, 0.62])
    plot_full_gmm(axis, density, membership, locked)
    axis.set_ylabel("Density", fontsize=6.2)
    axis.set_xlabel("eMTBR-tau243 (pg/mL, log scale)", fontsize=7)
    axis.tick_params(labelsize=5.8, length=2.5)

    legend_handles = [
        Line2D([0], [0], color=NEUTRAL, linewidth=1.1, linestyle=(0, (4, 2)), label="One-component fit"),
        Line2D([0], [0], color="#222222", linewidth=1.8, label="Two-component fit"),
        Line2D([0], [0], color="#222222", linewidth=1.0, linestyle=":", label="Component boundary"),
        Line2D([0], [0], marker="|", color=LOW_COLOR, linestyle="None", markersize=8, label="Low component"),
        Line2D([0], [0], marker="|", color=HIGH_COLOR, linestyle="None", markersize=8, label="High component"),
    ]
    figure.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                  ncol=3, fontsize=5.8, handlelength=2.2, columnspacing=1.0)
    save_figure(figure)

    density.to_csv(OUTPUT_DIR / "FigS5_Aplus_Tminus_heterogeneity_source.csv", index=False)
    metadata = {
        "core_conclusion": (
            "Within 18 tau PET-defined A+T− patients, a two-component mixture fit "
            "log10 eMTBR-tau243 better than a one-component model (exploratory)."
        ),
        "backend": "Python/matplotlib",
        "final_size_inches": [3.6, 2.9],
        "source_files": [
            str(APTM_DIR / "A+T-Subgroup_membership.csv"),
            str(APTM_DIR / "bimodality_validation_summary.csv"),
        ],
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    (OUTPUT_DIR / "FigS5_Aplus_Tminus_heterogeneity_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Figure written to {STEM}")


if __name__ == "__main__":
    main()
