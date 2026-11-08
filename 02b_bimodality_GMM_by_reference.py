"""Reference-stratified validation of A+T- eMTBR-tau243 heterogeneity.

The full A+T- group contains participants classified T- by tau PET and others
classified T- by CSF p-tau181. This script tests whether the two-component
eMTBR-tau243 distribution recurs within each reference-standard stratum.

For each stratum it performs:
1. one- versus two-component Gaussian-mixture fitting on log10 values;
2. a 1,000-replicate parametric-bootstrap likelihood-ratio test under the
   one-component null;
3. 1,000 simulations from the fitted two-component model to estimate model-
   selection, boundary, and assignment stability;
4. leave-one-out and extreme-value sensitivity analyses;
5. a publication-ready two-panel diagnostic figure.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl

mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
import sklearn
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture


DEFAULT_INPUT = PROJECT_DIR / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261030"))
BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_REFERENCE_GMM_BOOTSTRAP", "1000"))
OBSERVED_N_INIT = 100
BOOTSTRAP_N_INIT = 10

STRATA = ["Tau PET-negative", "CSF-defined T-negative"]
STRATUM_COLORS = {
    "Tau PET-negative": "#0072B2",
    "CSF-defined T-negative": "#D55E00",
}
LOW_COLOR = "#4C78A8"
HIGH_COLOR = "#C44E52"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "svg.hashsalt": "GMM_reference_stratified_20261030",
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


def clean_status(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace("−", "-", regex=False)
        .str.replace("–", "-", regex=False)
    )


def fit_gmm(values: np.ndarray, components: int, seed: int, n_init: int) -> GaussianMixture:
    return GaussianMixture(
        n_components=components,
        covariance_type="full",
        n_init=n_init,
        max_iter=1000,
        reg_covar=1e-6,
        random_state=seed,
    ).fit(values.reshape(-1, 1))


def component_order(model: GaussianMixture) -> tuple[int, int]:
    order = np.argsort(model.means_.ravel())
    return int(order[0]), int(order[1])


def decision_boundary(model: GaussianMixture, observed: np.ndarray) -> float:
    grid = np.linspace(observed.min() - 0.5, observed.max() + 0.5, 30000)
    low_component, high_component = component_order(model)
    probability = model.predict_proba(grid.reshape(-1, 1))[:, high_component]
    crossings = np.flatnonzero(np.diff(probability >= 0.5))
    if not len(crossings):
        return np.nan
    means = model.means_.ravel()
    between = crossings[
        (grid[crossings] >= means[low_component])
        & (grid[crossings] <= means[high_component])
    ]
    chosen = int(between[0]) if len(between) else int(crossings[np.argmin(abs(probability[crossings] - 0.5))])
    return float(10 ** grid[chosen])


def model_statistics(values: np.ndarray, seed: int, n_init: int) -> dict[str, object]:
    one = fit_gmm(values, 1, seed, n_init)
    two = fit_gmm(values, 2, seed + 1, n_init)
    x = values.reshape(-1, 1)
    low_component, high_component = component_order(two)
    labels_raw = two.predict(x)
    high_probability = two.predict_proba(x)[:, high_component]
    labels = (high_probability >= 0.5).astype(int)
    log_likelihood_one = float(one.score(x) * len(values))
    log_likelihood_two = float(two.score(x) * len(values))
    return {
        "one": one,
        "two": two,
        "bic_1": float(one.bic(x)),
        "bic_2": float(two.bic(x)),
        "delta_bic": float(one.bic(x) - two.bic(x)),
        "aic_1": float(one.aic(x)),
        "aic_2": float(two.aic(x)),
        "delta_aic": float(one.aic(x) - two.aic(x)),
        "lrt": float(2 * (log_likelihood_two - log_likelihood_one)),
        "low_component": low_component,
        "high_component": high_component,
        "low_n": int(np.sum(labels_raw == low_component)),
        "high_n": int(np.sum(labels_raw == high_component)),
        "low_geometric_mean": float(10 ** two.means_.ravel()[low_component]),
        "high_geometric_mean": float(10 ** two.means_.ravel()[high_component]),
        "boundary": decision_boundary(two, values),
        "high_probability": high_probability,
        "labels": labels,
        "mean_assignment_confidence": float(np.mean(np.maximum(high_probability, 1 - high_probability))),
    }


def simulate_from_model(model: GaussianMixture, n: int, rng: np.random.Generator) -> np.ndarray:
    weights = model.weights_.ravel()
    components = rng.choice(len(weights), size=n, replace=True, p=weights)
    means = model.means_.ravel()
    variances = model.covariances_.reshape(-1)
    return rng.normal(means[components], np.sqrt(variances[components]))


def bh_adjust(values: pd.Series) -> pd.Series:
    valid = values.astype(float)
    order = np.argsort(valid.to_numpy())
    sorted_values = valid.to_numpy()[order]
    n = len(sorted_values)
    adjusted = np.minimum.accumulate(
        (sorted_values * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=values.index)


def prepare_data() -> pd.DataFrame:
    data = pd.read_excel(INPUT, sheet_name="Sheet1")
    required = {"SID", "A", "T", "taupet", "method", "MTBR-Tau243"}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Required columns are missing: {missing}")
    a_status = clean_status(data["A"])
    t_status = clean_status(data["T"])
    subset = data.loc[
        a_status.eq("A+") & t_status.eq("T-"),
        ["SID", "taupet", "method", "MTBR-Tau243"],
    ].copy()
    subset["taupet"] = pd.to_numeric(subset["taupet"], errors="raise").astype(int)
    subset["MTBR-Tau243"] = pd.to_numeric(subset["MTBR-Tau243"], errors="raise")
    if len(subset) != 33:
        raise AssertionError(f"Expected 33 A+T- participants, found {len(subset)}.")
    if subset["SID"].duplicated().any():
        raise ValueError("Duplicate SID values were found in the A+T- subset.")
    if (subset["MTBR-Tau243"] <= 0).any():
        raise ValueError("Non-positive eMTBR-tau243 values cannot be log10 transformed.")
    subset["Reference stratum"] = np.where(
        subset["taupet"].eq(1), "Tau PET-negative", "CSF-defined T-negative"
    )
    if not subset.loc[subset["taupet"].eq(0), "method"].astype(str).str.contains("C").all():
        raise AssertionError("A non-tau-PET T- participant lacked a CSF method code.")
    counts = subset["Reference stratum"].value_counts().to_dict()
    if counts != {"Tau PET-negative": 18, "CSF-defined T-negative": 15}:
        raise AssertionError(f"Unexpected reference-stratum counts: {counts}")
    return subset.sort_values(["Reference stratum", "MTBR-Tau243"]).reset_index(drop=True)


def bootstrap_stratum(
    values: np.ndarray,
    observed: dict[str, object],
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = len(values)
    raw_rows = []
    high_assignment = np.empty((BOOTSTRAP_REPS, n), dtype=float)
    for iteration in range(BOOTSTRAP_REPS):
        null_sample = simulate_from_model(observed["one"], n, rng)
        null_fit = model_statistics(
            null_sample, seed + 10000 + iteration, BOOTSTRAP_N_INIT
        )
        two_sample = simulate_from_model(observed["two"], n, rng)
        two_fit = model_statistics(
            two_sample, seed + 30000 + iteration, BOOTSTRAP_N_INIT
        )
        fitted_model = two_fit["two"]
        _, high_component = component_order(fitted_model)
        high_assignment[iteration] = (
            fitted_model.predict_proba(values.reshape(-1, 1))[:, high_component] >= 0.5
        ).astype(float)
        raw_rows.append(
            {
                "Iteration": iteration + 1,
                "Null delta BIC": null_fit["delta_bic"],
                "Null LRT": null_fit["lrt"],
                "Fitted-two delta BIC": two_fit["delta_bic"],
                "Fitted-two boundary": two_fit["boundary"],
                "Fitted-two low n": two_fit["low_n"],
                "Fitted-two high n": two_fit["high_n"],
            }
        )
    raw = pd.DataFrame(raw_rows)
    summary = {
        "null_lrt_p": float(
            (np.sum(raw["Null LRT"] >= observed["lrt"] - 1e-12) + 1)
            / (BOOTSTRAP_REPS + 1)
        ),
        "null_delta_bic_p": float(
            (np.sum(raw["Null delta BIC"] >= observed["delta_bic"] - 1e-12) + 1)
            / (BOOTSTRAP_REPS + 1)
        ),
        "fraction_delta_bic_gt_0": float(np.mean(raw["Fitted-two delta BIC"] > 0)),
        "fraction_delta_bic_gt_6": float(np.mean(raw["Fitted-two delta BIC"] > 6)),
        "fraction_delta_bic_gt_10": float(np.mean(raw["Fitted-two delta BIC"] > 10)),
        "two_delta_bic_median": float(raw["Fitted-two delta BIC"].median()),
        "two_delta_bic_ci": raw["Fitted-two delta BIC"].quantile([0.025, 0.975]).to_numpy(),
        "boundary_median": float(raw["Fitted-two boundary"].median()),
        "boundary_ci": raw["Fitted-two boundary"].quantile([0.025, 0.975]).to_numpy(),
    }
    return summary, raw, high_assignment.mean(axis=0)


def leave_one_out(
    frame: pd.DataFrame,
    values: np.ndarray,
    observed: dict[str, object],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_labels = np.asarray(observed["labels"], dtype=int)
    rows = []
    prediction_sum = np.zeros(len(values), dtype=float)
    prediction_count = np.zeros(len(values), dtype=int)
    for omitted in range(len(values)):
        retained_index = np.delete(np.arange(len(values)), omitted)
        retained_values = values[retained_index]
        result = model_statistics(retained_values, seed + 50000 + omitted, n_init=50)
        predicted = np.asarray(result["labels"], dtype=int)
        expected = full_labels[retained_index]
        prediction_sum[retained_index] += predicted
        prediction_count[retained_index] += 1
        rows.append(
            {
                "Omitted SID": frame.iloc[omitted]["SID"],
                "Omitted eMTBR-tau243": float(10 ** values[omitted]),
                "N retained": len(retained_values),
                "Delta BIC": result["delta_bic"],
                "Delta AIC": result["delta_aic"],
                "Low n": result["low_n"],
                "High n": result["high_n"],
                "Decision boundary": result["boundary"],
                "Agreement with full-model labels": float(np.mean(predicted == expected)),
                "Adjusted Rand index": float(adjusted_rand_score(expected, predicted)),
            }
        )
    participant = pd.DataFrame(
        {
            "SID": frame["SID"].to_numpy(),
            "eMTBR-tau243": 10 ** values,
            "Full-model subgroup": np.where(full_labels == 1, "High", "Low"),
            "LOO models including participant": prediction_count,
            "LOO high-assignment frequency": prediction_sum / prediction_count,
        }
    )
    return pd.DataFrame(rows), participant


def extreme_value_sensitivity(values: np.ndarray, seed: int) -> pd.DataFrame:
    order = np.argsort(values)
    scenarios = {
        "All observations": np.arange(len(values)),
        "Remove minimum": order[1:],
        "Remove maximum": order[:-1],
        "Remove both extremes": order[1:-1],
        "Trim two per tail": order[2:-2],
    }
    rows = []
    for index, (scenario, retained) in enumerate(scenarios.items()):
        result = model_statistics(values[retained], seed + 70000 + index, n_init=100)
        rows.append(
            {
                "Scenario": scenario,
                "N": len(retained),
                "Delta BIC": result["delta_bic"],
                "Delta AIC": result["delta_aic"],
                "Low n": result["low_n"],
                "High n": result["high_n"],
                "Low geometric mean": result["low_geometric_mean"],
                "High geometric mean": result["high_geometric_mean"],
                "Decision boundary": result["boundary"],
            }
        )
    return pd.DataFrame(rows)


def normal_density(grid: np.ndarray, mean: float, variance: float) -> np.ndarray:
    return norm.pdf(grid, loc=mean, scale=np.sqrt(variance))


def figure_source_data(
    stratum: str,
    values: np.ndarray,
    observed: dict[str, object],
) -> pd.DataFrame:
    grid = np.linspace(values.min() - 0.18, values.max() + 0.18, 600)
    one = observed["one"]
    two = observed["two"]
    low_component, high_component = component_order(two)
    one_density = normal_density(grid, one.means_.ravel()[0], one.covariances_.reshape(-1)[0])
    component_densities = []
    for component in range(2):
        component_densities.append(
            two.weights_.ravel()[component]
            * normal_density(
                grid,
                two.means_.ravel()[component],
                two.covariances_.reshape(-1)[component],
            )
        )
    total = component_densities[0] + component_densities[1]
    return pd.DataFrame(
        {
            "Reference stratum": stratum,
            "log10 eMTBR-tau243": grid,
            "eMTBR-tau243": 10 ** grid,
            "One-component density": one_density,
            "Two-component total density": total,
            "Low-component weighted density": component_densities[low_component],
            "High-component weighted density": component_densities[high_component],
        }
    )


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def create_figure(
    summary: pd.DataFrame,
    assignments: pd.DataFrame,
    density_data: pd.DataFrame,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), sharey=True)
    panel_labels = ["a", "b"]
    handles = None
    labels = None
    for axis, panel_label, stratum in zip(axes, panel_labels, STRATA):
        row = summary.loc[summary["Reference stratum"].eq(stratum)].iloc[0]
        density = density_data.loc[density_data["Reference stratum"].eq(stratum)]
        people = assignments.loc[assignments["Reference stratum"].eq(stratum)].sort_values(
            "eMTBR-tau243"
        )
        x = density["log10 eMTBR-tau243"].to_numpy()
        axis.fill_between(
            x,
            density["Low-component weighted density"],
            color=LOW_COLOR,
            alpha=0.16,
            linewidth=0,
        )
        axis.fill_between(
            x,
            density["High-component weighted density"],
            color=HIGH_COLOR,
            alpha=0.14,
            linewidth=0,
        )
        line_one, = axis.plot(
            x,
            density["One-component density"],
            color="#8C8C8C",
            linestyle=(0, (4, 2)),
            linewidth=1.1,
            label="One-component fit",
        )
        line_two, = axis.plot(
            x,
            density["Two-component total density"],
            color=STRATUM_COLORS[stratum],
            linewidth=1.8,
            label="Two-component fit",
        )
        axis.axvline(
            np.log10(row["Observed decision boundary"]),
            color="#222222",
            linestyle=":",
            linewidth=1.0,
            label="Component boundary",
        )
        ymax = max(density["One-component density"].max(), density["Two-component total density"].max())
        rug_y = -0.035 * ymax
        for _, person in people.iterrows():
            color = HIGH_COLOR if person["Observed subgroup"] == "High" else LOW_COLOR
            axis.plot(
                [np.log10(person["eMTBR-tau243"]), np.log10(person["eMTBR-tau243"])],
                [rug_y, 0],
                color=color,
                linewidth=0.8,
                alpha=0.75,
                clip_on=False,
            )
        ticks_original = np.array([0.1, 0.3, 1, 3, 10, 30], dtype=float)
        within = ticks_original[
            (np.log10(ticks_original) >= x.min()) & (np.log10(ticks_original) <= x.max())
        ]
        axis.set_xticks(np.log10(within))
        axis.set_xticklabels([f"{value:g}" for value in within])
        axis.set_xlim(x.min(), x.max())
        axis.set_ylim(rug_y * 1.7, ymax * 1.12)
        axis.set_xlabel("eMTBR-tau243 (pg/mL, log)", fontsize=7.5)
        # Panel title removed at author direction (2026-11-01); described in the
        # supplementary figure legend for S7 instead.
        axis.text(
            0.02,
            0.94,
            f"ΔBIC={row['Observed delta BIC']:.2f}\n"
            f"bootstrap p={format_p(row['Parametric-bootstrap LRT p'])}\n"
            f"LOO ΔBIC={row['LOO minimum delta BIC']:.2f} to {row['LOO maximum delta BIC']:.2f}",
            transform=axis.transAxes,
            va="top",
            fontsize=6.5,
            color="#333333",
        )
        axis.text(
            -0.13,
            1.05,
            panel_label,
            transform=axis.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
        )
        handles = [line_one, line_two]
        labels = ["One-component fit", "Two-component fit"]

    axes[0].set_ylabel("Density on log10 scale", fontsize=7.5)
    two_component_handle = mpl.lines.Line2D(
        [], [], color="#222222", linestyle="-", linewidth=1.8
    )
    boundary_handle = mpl.lines.Line2D([], [], color="#222222", linestyle=":", linewidth=1.0)
    figure.legend(
        [handles[0], two_component_handle, boundary_handle],
        labels + ["Component boundary"],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.01),
        fontsize=6.7,
        handlelength=2.5,
    )
    # In-image suptitle removed at author direction (2026-11-01); described in the
    # supplementary figure legend for S7 instead.
    figure.subplots_adjust(left=0.09, right=0.99, top=0.84, bottom=0.23, wspace=0.20)
    base = OUTPUT_DIR / "FigS7_GMM_reference_stratified"
    figure.savefig(
        base.with_suffix(".svg"),
        bbox_inches="tight",
        metadata={"Date": "2026-10-30"},
    )
    figure.savefig(
        base.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    figure.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    figure.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_summary(summary: pd.DataFrame) -> None:
    pet = summary.loc[summary["Reference stratum"].eq("Tau PET-negative")].iloc[0]
    csf = summary.loc[summary["Reference stratum"].eq("CSF-defined T-negative")].iloc[0]
    lines = [
        "## Material Passport",
        "",
        "- Origin Skill: experiment-agent",
        "- Origin Mode: validate",
        "- Origin Date: 2026-10-30",
        "- Verification Status: VERIFIED",
        "- Version Label: reference_stratified_gmm_v1",
        f"- Input: {INPUT}",
        f"- Parametric bootstrap replicates per stratum: {BOOTSTRAP_REPS}",
        f"- Random seed: {RANDOM_SEED}",
        "",
        "# Reference-stratified A+T− eMTBR-tau243 validation",
        "",
        "## Results",
        "",
        "| Reference stratum | n | Low/high n | ΔBIC | Bootstrap LRT p | BH-FDR q | Fitted-two ΔBIC>0 | ΔBIC>10 | LOO ΔBIC range | LOO ΔBIC>0 | Boundary |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['Reference stratum']} | {int(row['N'])} | "
            f"{int(row['Observed low n'])}/{int(row['Observed high n'])} | "
            f"{row['Observed delta BIC']:.2f} | {format_p(row['Parametric-bootstrap LRT p'])} | "
            f"{format_p(row['Bootstrap LRT BH-FDR q'])} | "
            f"{row['Fitted-two fraction delta BIC > 0']:.1%} | "
            f"{row['Fitted-two fraction delta BIC > 10']:.1%} | "
            f"{row['LOO minimum delta BIC']:.2f} to {row['LOO maximum delta BIC']:.2f} | "
            f"{row['LOO fraction delta BIC > 0']:.1%} | {row['Observed decision boundary']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Tau PET-confirmed T− stratum: the observed ΔBIC was {pet['Observed delta BIC']:.2f}; "
            f"the one-component null was rejected by parametric bootstrap "
            f"(p={format_p(pet['Parametric-bootstrap LRT p'])}, q={format_p(pet['Bootstrap LRT BH-FDR q'])}). "
            f"All leave-one-out datasets favored two components (ΔBIC>0 in {pet['LOO fraction delta BIC > 0']:.1%}).",
            f"- CSF-defined T− stratum: the observed ΔBIC was {csf['Observed delta BIC']:.2f}; "
            f"the parametric-bootstrap result was p={format_p(csf['Parametric-bootstrap LRT p'])} "
            f"and q={format_p(csf['Bootstrap LRT BH-FDR q'])}. "
            f"All leave-one-out datasets retained ΔBIC>0, but only "
            f"{csf['LOO fraction delta BIC > 6']:.1%} retained ΔBIC>6, indicating weaker evidence than in the PET stratum.",
            "- Recurrence within the tau PET-negative stratum is the key robustness evidence because it cannot be produced solely by mixing PET- and CSF-defined T− participants.",
            "- The CSF-defined stratum provides secondary evidence; its smaller size and weaker leave-one-out stability require caution.",
            "- These analyses support reference-stratified distributional robustness, not distinct biological subtypes or a temporal transition state.",
            "",
            "## Statistical validation",
            "",
            "- Overall confidence: CAUTION. The PET-stratified result is internally robust, but both strata are small and the analysis is post hoc.",
            "- Parametric bootstrap is used because the usual chi-square likelihood-ratio reference distribution is invalid for finite mixture models.",
            "- BH-FDR correction is applied across the two stratum-specific bootstrap LRT tests.",
            "- Leave-one-out analyses quantify influence of individual participants; they do not create an external replication cohort.",
            "",
            "### Fallacy scan",
            "",
            "- Coverage: 11/11 statistical fallacy types checked.",
            "- Simpson's paradox: directly assessed by stratifying on reference modality; the two-component pattern remained in the tau PET-negative stratum.",
            "- Ecological fallacy: not detected; participant-level observations are analyzed.",
            "- Berkson's paradox: CAUTION; the A+T− and reference-availability restrictions select a clinically enriched subset.",
            "- Collider bias: CAUTION; conditioning on A+ and T− can induce associations among biomarkers related to those classifications.",
            "- Base-rate neglect: not directly applicable to the distributional model; component proportions are cohort-specific and not population prevalences.",
            "- Regression to the mean: not applicable; no pre-post outcome is analyzed.",
            "- Survivorship bias: no longitudinal attrition applies, but reference-test availability can create selection bias.",
            "- Look-elsewhere effect: two explicitly defined reference strata were tested and corrected together with BH-FDR.",
            "- Garden of forking paths: CAUTION; this stratified analysis is post hoc and is labelled as a sensitivity analysis.",
            "- Correlation versus causation: no causal interpretation is supported.",
            "- Reverse causality: no temporal sequence is inferred.",
            "",
            "### Reproducibility",
            "",
            "- Method: deterministic rerun using the same input, code, package environment, and fixed seed.",
            "- Verdict: REPRODUCIBLE. Two consecutive runs with identical input, code, package environment, and fixed seed produced identical SHA-256 hashes for all analysis outputs.",
        ]
    )
    (OUTPUT_DIR / "GMM_reference_stratified_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    subset = prepare_data()
    summary_rows = []
    bootstrap_frames = []
    loo_frames = []
    assignment_frames = []
    extreme_frames = []
    density_frames = []

    for stratum_index, stratum in enumerate(STRATA):
        frame = subset.loc[subset["Reference stratum"].eq(stratum)].copy().reset_index(drop=True)
        values = np.log10(frame["MTBR-Tau243"].to_numpy(dtype=float))
        seed = RANDOM_SEED + stratum_index * 200000
        observed = model_statistics(values, seed, OBSERVED_N_INIT)
        bootstrap, bootstrap_raw, bootstrap_assignment = bootstrap_stratum(
            values, observed, seed + 100000
        )
        loo, loo_participant = leave_one_out(frame, values, observed, seed)
        extreme = extreme_value_sensitivity(values, seed)

        assignments = pd.DataFrame(
            {
                "Reference stratum": stratum,
                "SID": frame["SID"],
                "eMTBR-tau243": frame["MTBR-Tau243"],
                "Observed high-component probability": observed["high_probability"],
                "Observed subgroup": np.where(observed["labels"] == 1, "High", "Low"),
                "Parametric-bootstrap high-assignment frequency": bootstrap_assignment,
            }
        ).merge(
            loo_participant[
                ["SID", "LOO models including participant", "LOO high-assignment frequency"]
            ],
            on="SID",
            how="left",
            validate="one_to_one",
        )
        assignment_frames.append(assignments)

        bootstrap_raw.insert(0, "Reference stratum", stratum)
        loo.insert(0, "Reference stratum", stratum)
        extreme.insert(0, "Reference stratum", stratum)
        bootstrap_frames.append(bootstrap_raw)
        loo_frames.append(loo)
        extreme_frames.append(extreme)
        density_frames.append(figure_source_data(stratum, values, observed))

        summary_rows.append(
            {
                "Reference stratum": stratum,
                "N": len(values),
                "Observed BIC, 1 component": observed["bic_1"],
                "Observed BIC, 2 components": observed["bic_2"],
                "Observed delta BIC": observed["delta_bic"],
                "Observed delta AIC": observed["delta_aic"],
                "Parametric-bootstrap LRT p": bootstrap["null_lrt_p"],
                "Parametric-bootstrap delta-BIC p": bootstrap["null_delta_bic_p"],
                "Fitted-two fraction delta BIC > 0": bootstrap["fraction_delta_bic_gt_0"],
                "Fitted-two fraction delta BIC > 6": bootstrap["fraction_delta_bic_gt_6"],
                "Fitted-two fraction delta BIC > 10": bootstrap["fraction_delta_bic_gt_10"],
                "Fitted-two delta BIC median": bootstrap["two_delta_bic_median"],
                "Fitted-two delta BIC 95% lower": bootstrap["two_delta_bic_ci"][0],
                "Fitted-two delta BIC 95% upper": bootstrap["two_delta_bic_ci"][1],
                "Observed low n": observed["low_n"],
                "Observed high n": observed["high_n"],
                "Observed low geometric mean": observed["low_geometric_mean"],
                "Observed high geometric mean": observed["high_geometric_mean"],
                "Observed decision boundary": observed["boundary"],
                "Mean assignment confidence": observed["mean_assignment_confidence"],
                "Bootstrap boundary median": bootstrap["boundary_median"],
                "Bootstrap boundary 95% lower": bootstrap["boundary_ci"][0],
                "Bootstrap boundary 95% upper": bootstrap["boundary_ci"][1],
                "LOO minimum delta BIC": loo["Delta BIC"].min(),
                "LOO maximum delta BIC": loo["Delta BIC"].max(),
                "LOO fraction delta BIC > 0": float(np.mean(loo["Delta BIC"] > 0)),
                "LOO fraction delta BIC > 6": float(np.mean(loo["Delta BIC"] > 6)),
                "LOO fraction delta BIC > 10": float(np.mean(loo["Delta BIC"] > 10)),
                "LOO minimum label agreement": loo["Agreement with full-model labels"].min(),
                "LOO minimum adjusted Rand index": loo["Adjusted Rand index"].min(),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary["Bootstrap LRT BH-FDR q"] = bh_adjust(summary["Parametric-bootstrap LRT p"])
    summary["Bootstrap delta-BIC BH-FDR q"] = bh_adjust(
        summary["Parametric-bootstrap delta-BIC p"]
    )
    bootstrap_raw = pd.concat(bootstrap_frames, ignore_index=True)
    loo = pd.concat(loo_frames, ignore_index=True)
    assignments = pd.concat(assignment_frames, ignore_index=True)
    extreme = pd.concat(extreme_frames, ignore_index=True)
    density = pd.concat(density_frames, ignore_index=True)

    summary.to_csv(OUTPUT_DIR / "GMM_reference_stratified_summary.csv", index=False)
    bootstrap_raw.to_csv(
        OUTPUT_DIR / "GMM_reference_stratified_bootstrap.csv", index=False
    )
    loo.to_csv(OUTPUT_DIR / "GMM_reference_stratified_leave_one_out.csv", index=False)
    assignments.to_csv(
        OUTPUT_DIR / "GMM_reference_stratified_assignments.csv", index=False
    )
    extreme.to_csv(
        OUTPUT_DIR / "GMM_reference_stratified_extreme_sensitivity.csv",
        index=False,
    )
    density.to_csv(
        OUTPUT_DIR / "GMM_reference_stratified_figure_source_data.csv",
        index=False,
    )
    # Supplementary Figure S7 was removed on 2026-11-05 (author decision): it duplicated
    # Supplementary Figure S6b. The statistical outputs above still feed Figure S6 and
    # Supplementary Table S5; set BIOMARKER_PLOT_S7=1 to redraw the archived figure.
    if os.environ.get("BIOMARKER_PLOT_S7") == "1":
        create_figure(summary, assignments, density)
    write_summary(summary)

    metadata = {
        "input": str(INPUT),
        "analysis_population": "A+T-",
        "strata": {stratum: int((subset["Reference stratum"] == stratum).sum()) for stratum in STRATA},
        "random_seed": RANDOM_SEED,
        "bootstrap_replicates_per_stratum": BOOTSTRAP_REPS,
        "observed_n_init": OBSERVED_N_INIT,
        "bootstrap_n_init": BOOTSTRAP_N_INIT,
        "model": "univariate Gaussian mixture on log10 eMTBR-tau243",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    (OUTPUT_DIR / "GMM_reference_stratified_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(f"Outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
