"""Robustness analyses for MTBR-Tau243 bimodality within A+T-.

Primary variable: log10(MTBR-Tau243).
Negative-comparison variable: log10(p-Tau217), which is not used to define the
subgroups. A parametric bootstrap is used because the usual chi-square
likelihood-ratio reference distribution is invalid for finite mixture models.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import sklearn
from sklearn.mixture import GaussianMixture


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = PROJECT_DIR / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_MIXTURE_BOOTSTRAP", "1000"))

MARKERS = {
    "MTBR-Tau243": "MTBR-Tau243",
    "p-Tau217": "p-Tau217(≤0.758)",
}


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


def model_statistics(values: np.ndarray, seed: int, n_init: int = 100) -> dict[str, object]:
    one = fit_gmm(values, 1, seed, n_init)
    two = fit_gmm(values, 2, seed + 1, n_init)
    x = values.reshape(-1, 1)
    order = np.argsort(two.means_.ravel())
    low_component, high_component = int(order[0]), int(order[1])
    labels = two.predict(x)
    high_probability = two.predict_proba(x)[:, high_component]
    return {
        "one": one,
        "two": two,
        "bic_1": float(one.bic(x)),
        "bic_2": float(two.bic(x)),
        "delta_bic": float(one.bic(x) - two.bic(x)),
        "aic_1": float(one.aic(x)),
        "aic_2": float(two.aic(x)),
        "delta_aic": float(one.aic(x) - two.aic(x)),
        "lrt": float(2 * (two.score(x) * len(values) - one.score(x) * len(values))),
        "low_component": low_component,
        "high_component": high_component,
        "low_n": int(np.sum(labels == low_component)),
        "high_n": int(np.sum(labels == high_component)),
        "low_geometric_mean": float(10 ** two.means_.ravel()[low_component]),
        "high_geometric_mean": float(10 ** two.means_.ravel()[high_component]),
        "mean_assignment_confidence": float(np.mean(np.maximum(high_probability, 1 - high_probability))),
        "high_probability": high_probability,
    }


def simulate_from_model(model: GaussianMixture, n: int, rng: np.random.Generator) -> np.ndarray:
    weights = model.weights_.ravel()
    components = rng.choice(len(weights), size=n, replace=True, p=weights)
    means = model.means_.ravel()
    variances = model.covariances_.reshape(-1)
    return rng.normal(means[components], np.sqrt(variances[components]))


def decision_boundary(two: GaussianMixture, observed: np.ndarray) -> float:
    grid = np.linspace(observed.min() - 0.5, observed.max() + 0.5, 20000)
    order = np.argsort(two.means_.ravel())
    high = int(order[1])
    probability = two.predict_proba(grid.reshape(-1, 1))[:, high]
    crossing = np.flatnonzero(np.diff(probability >= 0.5))
    if not len(crossing):
        return np.nan
    between_means = crossing[
        (grid[crossing] >= two.means_.ravel()[order[0]])
        & (grid[crossing] <= two.means_.ravel()[order[1]])
    ]
    selected = between_means[0] if len(between_means) else crossing[np.argmin(abs(probability[crossing] - 0.5))]
    return float(10 ** grid[selected])


def bootstrap_validation(values: np.ndarray, observed: dict[str, object], seed: int) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    n = len(values)
    null_delta_bic = np.empty(BOOTSTRAP_REPS)
    null_lrt = np.empty(BOOTSTRAP_REPS)
    fitted_two_delta_bic = np.empty(BOOTSTRAP_REPS)
    fitted_two_boundary = np.empty(BOOTSTRAP_REPS)
    fitted_two_low_n = np.empty(BOOTSTRAP_REPS)
    original_high_assignments = np.empty((BOOTSTRAP_REPS, n), dtype=float)

    for i in range(BOOTSTRAP_REPS):
        null_sample = simulate_from_model(observed["one"], n, rng)
        null_fit = model_statistics(null_sample, seed + 10000 + i, n_init=10)
        null_delta_bic[i] = null_fit["delta_bic"]
        null_lrt[i] = null_fit["lrt"]

        two_sample = simulate_from_model(observed["two"], n, rng)
        two_fit = model_statistics(two_sample, seed + 30000 + i, n_init=10)
        fitted_two_delta_bic[i] = two_fit["delta_bic"]
        fitted_two_boundary[i] = decision_boundary(two_fit["two"], two_sample)
        fitted_two_low_n[i] = two_fit["low_n"]

        model = two_fit["two"]
        high_component = int(np.argmax(model.means_.ravel()))
        original_high_assignments[i] = (
            model.predict_proba(values.reshape(-1, 1))[:, high_component] >= 0.5
        ).astype(float)

    null_p_lrt = (np.sum(null_lrt >= observed["lrt"] - 1e-12) + 1) / (BOOTSTRAP_REPS + 1)
    null_p_bic = (np.sum(null_delta_bic >= observed["delta_bic"] - 1e-12) + 1) / (BOOTSTRAP_REPS + 1)
    return {
        "null_lrt_p": float(null_p_lrt),
        "null_delta_bic_p": float(null_p_bic),
        "two_model_delta_bic_median": float(np.median(fitted_two_delta_bic)),
        "two_model_delta_bic_ci": np.quantile(fitted_two_delta_bic, [0.025, 0.975]),
        "two_model_fraction_delta_bic_gt_0": float(np.mean(fitted_two_delta_bic > 0)),
        "two_model_fraction_delta_bic_gt_6": float(np.mean(fitted_two_delta_bic > 6)),
        "two_model_fraction_delta_bic_gt_10": float(np.mean(fitted_two_delta_bic > 10)),
        "boundary_median": float(np.nanmedian(fitted_two_boundary)),
        "boundary_ci": np.nanquantile(fitted_two_boundary, [0.025, 0.975]),
        "low_n_median": float(np.median(fitted_two_low_n)),
        "low_n_ci": np.quantile(fitted_two_low_n, [0.025, 0.975]),
        "participant_high_assignment_frequency": original_high_assignments.mean(axis=0),
    }


def sensitivity_analyses(values: np.ndarray, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    loo_rows = []
    for omitted in range(len(values)):
        retained = np.delete(values, omitted)
        result = model_statistics(retained, seed + 50000 + omitted, n_init=50)
        loo_rows.append(
            {
                "Omitted row position": omitted,
                "Omitted original value": 10 ** values[omitted],
                "Delta BIC": result["delta_bic"],
                "Delta AIC": result["delta_aic"],
                "Low n": result["low_n"],
                "High n": result["high_n"],
            }
        )

    order = np.argsort(values)
    scenarios = {
        "All observations": np.arange(len(values)),
        "Remove minimum": order[1:],
        "Remove maximum": order[:-1],
        "Remove both extremes": order[1:-1],
        "Trim two per tail": order[2:-2],
    }
    sensitivity_rows = []
    for i, (scenario, retained) in enumerate(scenarios.items()):
        result = model_statistics(values[retained], seed + 60000 + i, n_init=100)
        sensitivity_rows.append(
            {
                "Scenario": scenario,
                "N": len(retained),
                "Delta BIC": result["delta_bic"],
                "Delta AIC": result["delta_aic"],
                "Low n": result["low_n"],
                "High n": result["high_n"],
                "Low geometric mean": result["low_geometric_mean"],
                "High geometric mean": result["high_geometric_mean"],
                "Decision boundary": decision_boundary(result["two"], values[retained]),
            }
        )
    return pd.DataFrame(loo_rows), pd.DataFrame(sensitivity_rows)


def prepare_subset() -> pd.DataFrame:
    data = pd.read_excel(INPUT, sheet_name="Sheet1")
    required = {"SID", "A", "T", *MARKERS.values()}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Required columns are missing: {missing}")
    mask = clean_status(data["A"]).eq("A+") & clean_status(data["T"]).eq("T-")
    subset = data.loc[mask, ["SID", *MARKERS.values()]].copy()
    for column in MARKERS.values():
        subset[column] = pd.to_numeric(subset[column], errors="coerce")
    if subset[list(MARKERS.values())].isna().any().any():
        raise ValueError("Missing clustering-marker values were found in A+T- participants.")
    if (subset[list(MARKERS.values())] <= 0).any().any():
        raise ValueError("Non-positive clustering-marker values cannot be log10 transformed.")
    return subset.reset_index(drop=True)


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def main() -> None:
    subset = prepare_subset()
    summary_rows = []
    all_loo = []
    all_sensitivity = []
    primary_assignment = None

    for marker_index, (marker_name, source_column) in enumerate(MARKERS.items()):
        values = np.log10(subset[source_column].to_numpy(dtype=float))
        observed = model_statistics(values, RANDOM_SEED + marker_index * 1000, n_init=100)
        bootstrap = bootstrap_validation(values, observed, RANDOM_SEED + 100000 + marker_index * 100000)
        loo, sensitivity = sensitivity_analyses(values, RANDOM_SEED + marker_index * 1000)
        loo.insert(0, "Marker", marker_name)
        sensitivity.insert(0, "Marker", marker_name)
        all_loo.append(loo)
        all_sensitivity.append(sensitivity)

        observed_boundary = decision_boundary(observed["two"], values)
        summary_rows.append(
            {
                "Marker": marker_name,
                "N": len(values),
                "Observed BIC, 1 component": observed["bic_1"],
                "Observed BIC, 2 components": observed["bic_2"],
                "Observed delta BIC": observed["delta_bic"],
                "Observed delta AIC": observed["delta_aic"],
                "Parametric-bootstrap LRT p under 1-component null": bootstrap["null_lrt_p"],
                "Parametric-bootstrap delta-BIC p under 1-component null": bootstrap["null_delta_bic_p"],
                "Fraction selecting 2 components under fitted 2-component model (delta BIC > 0)": bootstrap["two_model_fraction_delta_bic_gt_0"],
                "Fraction with strong 2-component evidence (delta BIC > 6)": bootstrap["two_model_fraction_delta_bic_gt_6"],
                "Fraction with very strong 2-component evidence (delta BIC > 10)": bootstrap["two_model_fraction_delta_bic_gt_10"],
                "Two-model simulated delta BIC median": bootstrap["two_model_delta_bic_median"],
                "Two-model simulated delta BIC 95% lower": bootstrap["two_model_delta_bic_ci"][0],
                "Two-model simulated delta BIC 95% upper": bootstrap["two_model_delta_bic_ci"][1],
                "Observed low n": observed["low_n"],
                "Observed high n": observed["high_n"],
                "Observed low geometric mean": observed["low_geometric_mean"],
                "Observed high geometric mean": observed["high_geometric_mean"],
                "Observed decision boundary": observed_boundary,
                "Simulated boundary median": bootstrap["boundary_median"],
                "Simulated boundary 95% lower": bootstrap["boundary_ci"][0],
                "Simulated boundary 95% upper": bootstrap["boundary_ci"][1],
                "LOO minimum delta BIC": loo["Delta BIC"].min(),
                "LOO maximum delta BIC": loo["Delta BIC"].max(),
                "LOO fraction delta BIC > 10": float(np.mean(loo["Delta BIC"] > 10)),
            }
        )

        if marker_name == "MTBR-Tau243":
            primary_assignment = pd.DataFrame(
                {
                    "SID": subset["SID"],
                    "MTBR-Tau243": subset[source_column],
                    "Observed high-component probability": observed["high_probability"],
                    "Observed subgroup": np.where(
                        observed["high_probability"] >= 0.5,
                        "High MTBR-Tau243",
                        "Low MTBR-Tau243",
                    ),
                    "Parametric-bootstrap high-assignment frequency": bootstrap[
                        "participant_high_assignment_frequency"
                    ],
                }
            ).sort_values("MTBR-Tau243")

    summary = pd.DataFrame(summary_rows)
    loo = pd.concat(all_loo, ignore_index=True)
    sensitivity = pd.concat(all_sensitivity, ignore_index=True)
    summary.to_csv(OUTPUT_DIR / "bimodality_validation_summary.csv", index=False)
    loo.to_csv(OUTPUT_DIR / "bimodality_leave_one_out.csv", index=False)
    sensitivity.to_csv(OUTPUT_DIR / "bimodality_extreme_value_sensitivity.csv", index=False)
    primary_assignment.to_csv(OUTPUT_DIR / "MTBR_subgroup_assignment_stability.csv", index=False)

    mtbr = summary.loc[summary["Marker"] == "MTBR-Tau243"].iloc[0]
    ptau = summary.loc[summary["Marker"] == "p-Tau217"].iloc[0]
    lines = [
        "## Material Passport",
        "",
        "- artifact_type: mixture_model_validation",
        "- verification_status: VERIFIED",
        f"- input: {INPUT}",
        f"- seed: {RANDOM_SEED}",
        f"- parametric_bootstrap_replicates: {BOOTSTRAP_REPS}",
        "",
        "# A+T− bimodality validation summary",
        "",
        "## MTBR-Tau243",
        "",
        f"- Observed delta BIC: {mtbr['Observed delta BIC']:.2f}.",
        f"- Parametric-bootstrap LRT p under the one-component null: {format_p(mtbr['Parametric-bootstrap LRT p under 1-component null'])}.",
        f"- Fraction of simulations from the fitted two-component model retaining delta BIC > 0: {mtbr['Fraction selecting 2 components under fitted 2-component model (delta BIC > 0)']:.1%}.",
        f"- Fraction retaining delta BIC > 10: {mtbr['Fraction with very strong 2-component evidence (delta BIC > 10)']:.1%}.",
        f"- Leave-one-out delta BIC range: {mtbr['LOO minimum delta BIC']:.2f} to {mtbr['LOO maximum delta BIC']:.2f}.",
        f"- Model-implied decision boundary: {mtbr['Observed decision boundary']:.3f}.",
        "",
        "## p-Tau217 negative comparison",
        "",
        f"- Observed delta BIC: {ptau['Observed delta BIC']:.2f}.",
        f"- Parametric-bootstrap LRT p under the one-component null: {format_p(ptau['Parametric-bootstrap LRT p under 1-component null'])}.",
        "- p-Tau217 remains external to the subgroup definition.",
        "",
        "## Interpretation limits",
        "",
        "- Parametric-bootstrap results assess model stability under the fitted Gaussian-mixture assumptions.",
        "- They do not establish that the components are distinct biological subtypes.",
        "- Cross-sectional data cannot establish temporal progression from Low to High MTBR-Tau243 or from T− to T+.",
    ]
    (OUTPUT_DIR / "bimodality_validation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    metadata = {
        "input": str(INPUT),
        "random_seed": RANDOM_SEED,
        "bootstrap_replicates": BOOTSTRAP_REPS,
        "model": "univariate Gaussian mixture on log10 biomarker values",
        "observed_n_init": 100,
        "bootstrap_n_init": 10,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
    }
    (OUTPUT_DIR / "bimodality_validation_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
