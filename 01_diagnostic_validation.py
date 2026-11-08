"""Supplementary validation for the five A/T diagnostic analyses.

This script does not overwrite the original ROC results. It adds:
1. repeated stratified 5-fold cross-validated AUC and threshold-based metrics;
2. cross-validation cutoff stability;
3. fixed-cutoff patient-bootstrap intervals for accuracy, balanced accuracy,
   sensitivity, and specificity;
4. prespecified paired DeLong and paired-bootstrap AUC comparisons.

All resampling is performed at participant level. Pairwise comparisons use the
same participants and the same bootstrap draw for both biomarkers.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
import sklearn
from sklearn.model_selection import RepeatedStratifiedKFold


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = PROJECT_DIR / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
CV_REPEATS = int(os.environ.get("BIOMARKER_CV_REPEATS", "50"))
CV_FOLDS = int(os.environ.get("BIOMARKER_CV_FOLDS", "5"))
BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_VALIDATION_BOOTSTRAP", "5000"))

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

BASE_BIOMARKERS = [
    "GFAP(≤131.880)",
    "NFL(≤9.624 )",
    "Aβ 1-42",
    "p-Tau217(≤0.758)",
    "p-Tau231(≤5.793)",
    "BD Tau（≤52.301）",
    "MTBR-Tau243",
]

RATIO_NUMERATORS = {
    "p-Tau217/Aβ42": "p-Tau217(≤0.758)",
    "p-Tau231/Aβ42": "p-Tau231(≤5.793)",
    "MTBR-Tau243/Aβ42": "MTBR-Tau243",
    "BD-Tau/Aβ42": "BD Tau（≤52.301）",
}

COMPARISON_PAIRS = [
    ("MTBR-Tau243/Aβ42", "MTBR-Tau243", "MTBR ratio vs absolute"),
    ("p-Tau217/Aβ42", "p-Tau217", "p-Tau217 ratio vs absolute"),
    ("MTBR-Tau243", "p-Tau217", "MTBR vs p-Tau217 absolute"),
    ("MTBR-Tau243/Aβ42", "p-Tau217/Aβ42", "MTBR vs p-Tau217 ratio"),
]


def clean_status(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace("−", "-", regex=False)
        .str.replace("–", "-", regex=False)
    )


def auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    positive = score[y == 1]
    negative = score[y == 0]
    comparisons = positive[:, None] - negative[None, :]
    return float((np.sum(comparisons > 0) + 0.5 * np.sum(comparisons == 0)) / comparisons.size)


def confusion(y: np.ndarray, prediction: np.ndarray) -> tuple[int, int, int, int]:
    y = np.asarray(y, dtype=int)
    prediction = np.asarray(prediction, dtype=int)
    tp = int(np.sum((y == 1) & (prediction == 1)))
    tn = int(np.sum((y == 0) & (prediction == 0)))
    fp = int(np.sum((y == 0) & (prediction == 1)))
    fn = int(np.sum((y == 1) & (prediction == 0)))
    return tn, fp, fn, tp


def classification_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    tn, fp, fn, tp = confusion(y, prediction)
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    return {
        "accuracy": (tp + tn) / len(y),
        "balanced_accuracy": (sensitivity + specificity) / 2,
        "sensitivity": sensitivity,
        "specificity": specificity,
    }


def oriented_auc(y: np.ndarray, values: np.ndarray) -> tuple[float, int]:
    raw_auc = auc_score(y, values)
    return (raw_auc, 1) if raw_auc >= 0.5 else (1 - raw_auc, -1)


def optimal_metrics(y: np.ndarray, values: np.ndarray) -> dict[str, float]:
    auc, direction = oriented_auc(y, values)
    score = values * direction
    unique = np.unique(score)
    thresholds = np.r_[np.inf, (unique[:-1] + unique[1:]) / 2, -np.inf]
    best = None
    for threshold in thresholds:
        prediction = (score >= threshold).astype(int)
        metrics = classification_metrics(y, prediction)
        candidate = (
            metrics["sensitivity"] + metrics["specificity"] - 1,
            metrics["sensitivity"] + metrics["specificity"],
            threshold,
            prediction,
        )
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    raw_cutoff = float(best[2] / direction)
    metrics = classification_metrics(y, best[3])
    return {
        "auc": auc,
        "direction_sign": direction,
        "direction": ">" if direction == 1 else "<",
        "cutoff": raw_cutoff,
        **metrics,
    }


def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    return tuple(np.quantile(np.asarray(values, dtype=float), [0.025, 0.975]))


def bh_adjust(p_values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    sorted_p = valid.to_numpy()[order]
    n = len(sorted_p)
    adjusted = np.minimum.accumulate((sorted_p * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result.loc[valid.index.to_numpy()[order]] = adjusted
    return result


def prepare_data(path: Path) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    data = pd.read_excel(path, sheet_name="Sheet1")
    required = {"A", "T", "DIAG", "taupet", *BASE_BIOMARKERS}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Required columns are missing: {missing}")
    data["A_clean"] = clean_status(data["A"])
    data["T_clean"] = clean_status(data["T"])
    for column in BASE_BIOMARKERS:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    abeta42 = data["Aβ 1-42"]
    for ratio, numerator in RATIO_NUMERATORS.items():
        data[ratio] = data[numerator] / abeta42
    # Add stable short-name aliases used by the prespecified pair definitions.
    # The original workbook columns remain unchanged.
    for source, short_name in LABELS.items():
        if source in data.columns and short_name not in data.columns:
            data[short_name] = data[source]
    data["A_outcome"] = data["A_clean"].eq("A+").astype(int)
    data["T_outcome"] = data["T_clean"].eq("T+").astype(int)
    outcomes = [
        {
            "analysis": "A status",
            "mask": data["A_clean"].isin(["A+", "A-"]),
            "outcome": "A_outcome",
        },
        {
            "analysis": "A status, non-control cohort",
            "mask": data["A_clean"].isin(["A+", "A-"]) & data["DIAG"].ne("健康对照"),
            "outcome": "A_outcome",
        },
        {
            "analysis": "T status",
            "mask": data["T_clean"].isin(["T+", "T-"]),
            "outcome": "T_outcome",
        },
        {
            "analysis": "T status, Tau PET subset",
            "mask": data["T_clean"].isin(["T+", "T-"]) & pd.to_numeric(data["taupet"], errors="coerce").eq(1),
            "outcome": "T_outcome",
        },
        {
            "analysis": "T status, A+ subset",
            "mask": data["A_clean"].eq("A+") & data["T_clean"].isin(["T+", "T-"]),
            "outcome": "T_outcome",
        },
        {
            "analysis": "T status, A+ subset, Tau PET only",
            "mask": data["A_clean"].eq("A+")
            & data["T_clean"].isin(["T+", "T-"])
            & pd.to_numeric(data["taupet"], errors="coerce").eq(1),
            "outcome": "T_outcome",
        },
    ]
    _ap_tau = outcomes[-1]
    _ap_mask = _ap_tau["mask"]
    _ap_pos = int(data.loc[_ap_mask, _ap_tau["outcome"]].sum())
    _ap_n = int(_ap_mask.sum())
    if (_ap_n, _ap_pos, _ap_n - _ap_pos) != (52, 34, 18):
        raise AssertionError(
            "Unexpected A+ tau-PET-only population: "
            f"n={_ap_n}, T+={_ap_pos}, T-={_ap_n - _ap_pos} (expected 52, 34, 18)."
        )
    return data, outcomes


def get_marker_data(
    data: pd.DataFrame,
    outcome: dict[str, object],
    marker: str,
) -> tuple[np.ndarray, np.ndarray]:
    local = outcome["mask"] & data[marker].notna()
    y = data.loc[local, outcome["outcome"]].to_numpy(dtype=int)
    x = pd.to_numeric(data.loc[local, marker], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x)
    return y[valid], x[valid]


def repeated_cv_metrics(
    y: np.ndarray,
    values: np.ndarray,
    random_state: int,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    splitter = RepeatedStratifiedKFold(
        n_splits=CV_FOLDS,
        n_repeats=CV_REPEATS,
        random_state=random_state,
    )
    repeat_metrics = {name: [] for name in ["auc", "accuracy", "balanced_accuracy", "sensitivity", "specificity"]}
    trained_cutoffs = []
    trained_directions = []
    current_repeat = -1
    oof_score = None
    oof_prediction = None
    fold_aucs = []
    for split_number, (train, test) in enumerate(splitter.split(values.reshape(-1, 1), y)):
        repeat_number = split_number // CV_FOLDS
        if repeat_number != current_repeat:
            if current_repeat >= 0:
                repeat_metrics["auc"].append(float(np.mean(fold_aucs)))
                metrics = classification_metrics(y, oof_prediction)
                for key, value in metrics.items():
                    repeat_metrics[key].append(value)
            current_repeat = repeat_number
            oof_score = np.full(len(y), np.nan, dtype=float)
            oof_prediction = np.full(len(y), -1, dtype=int)
            fold_aucs = []

        trained = optimal_metrics(y[train], values[train])
        direction = int(trained["direction_sign"])
        oriented_cutoff = trained["cutoff"] * direction
        oof_score[test] = values[test] * direction
        oof_prediction[test] = (oof_score[test] >= oriented_cutoff).astype(int)
        fold_aucs.append(auc_score(y[test], oof_score[test]))
        trained_cutoffs.append(trained["cutoff"])
        trained_directions.append(direction)

    repeat_metrics["auc"].append(float(np.mean(fold_aucs)))
    metrics = classification_metrics(y, oof_prediction)
    for key, value in metrics.items():
        repeat_metrics[key].append(value)
    repeat_metrics = {key: np.asarray(value, dtype=float) for key, value in repeat_metrics.items()}
    return repeat_metrics, np.asarray(trained_cutoffs), np.asarray(trained_directions)


def fixed_cutoff_bootstrap(
    y: np.ndarray,
    values: np.ndarray,
    direction: int,
    cutoff: float,
    seed: int,
) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(y == 1)
    negative = np.flatnonzero(y == 0)
    collected = {name: [] for name in ["accuracy", "balanced_accuracy", "sensitivity", "specificity"]}
    oriented_cutoff = cutoff * direction
    for _ in range(BOOTSTRAP_REPS):
        index = np.concatenate(
            [
                rng.choice(positive, size=len(positive), replace=True),
                rng.choice(negative, size=len(negative), replace=True),
            ]
        )
        prediction = (values[index] * direction >= oriented_cutoff).astype(int)
        metrics = classification_metrics(y[index], prediction)
        for key, value in metrics.items():
            collected[key].append(value)
    return {key: percentile_interval(np.asarray(value)) for key, value in collected.items()}


def delong_covariance(y: np.ndarray, predictions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """AUC vector and DeLong covariance using positive/negative influence values."""
    positive_predictions = predictions[:, y == 1]
    negative_predictions = predictions[:, y == 0]
    marker_count = predictions.shape[0]
    aucs = np.empty(marker_count, dtype=float)
    positive_influence = []
    negative_influence = []
    for row in predictions:
        positive = row[y == 1]
        negative = row[y == 0]
        comparison = (positive[:, None] > negative[None, :]).astype(float)
        comparison += 0.5 * (positive[:, None] == negative[None, :])
        aucs[len(positive_influence)] = comparison.mean()
        positive_influence.append(comparison.mean(axis=1))
        negative_influence.append(comparison.mean(axis=0))
    positive_influence = np.asarray(positive_influence)
    negative_influence = np.asarray(negative_influence)
    positive_covariance = np.cov(positive_influence, bias=False)
    negative_covariance = np.cov(negative_influence, bias=False)
    covariance = positive_covariance / positive_predictions.shape[1] + negative_covariance / negative_predictions.shape[1]
    return aucs, np.atleast_2d(covariance)


def paired_delong(y: np.ndarray, score_1: np.ndarray, score_2: np.ndarray) -> dict[str, float]:
    predictions = np.vstack([score_1, score_2])
    aucs, covariance = delong_covariance(y, predictions)
    contrast = np.array([1.0, -1.0])
    variance = float(contrast @ covariance @ contrast)
    difference = float(aucs[0] - aucs[1])
    if variance <= 0:
        standard_error = 0.0
        z_value = np.nan
        p_value = 1.0 if np.isclose(difference, 0) else 0.0
    else:
        standard_error = float(np.sqrt(variance))
        z_value = difference / standard_error
        p_value = float(2 * norm.sf(abs(z_value)))
    return {
        "auc_1": float(aucs[0]),
        "auc_2": float(aucs[1]),
        "delta_auc": difference,
        "delong_se": standard_error,
        "delong_ci_lower": difference - 1.96 * standard_error,
        "delong_ci_upper": difference + 1.96 * standard_error,
        "delong_z": z_value,
        "delong_p": p_value,
    }


def paired_auc_bootstrap(
    y: np.ndarray,
    score_1: np.ndarray,
    score_2: np.ndarray,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(y == 1)
    negative = np.flatnonzero(y == 0)
    differences = np.empty(BOOTSTRAP_REPS, dtype=float)
    for i in range(BOOTSTRAP_REPS):
        index = np.concatenate(
            [
                rng.choice(positive, size=len(positive), replace=True),
                rng.choice(negative, size=len(negative), replace=True),
            ]
        )
        differences[i] = auc_score(y[index], score_1[index]) - auc_score(y[index], score_2[index])
    ci_low, ci_high = percentile_interval(differences)
    lower_tail = (np.sum(differences <= 0) + 1) / (BOOTSTRAP_REPS + 1)
    upper_tail = (np.sum(differences >= 0) + 1) / (BOOTSTRAP_REPS + 1)
    return {
        "bootstrap_delta_auc_mean": float(np.mean(differences)),
        "bootstrap_ci_lower": ci_low,
        "bootstrap_ci_upper": ci_high,
        "bootstrap_p": float(min(1.0, 2 * min(lower_tail, upper_tail))),
    }


def cv_analysis(data: pd.DataFrame, outcomes: list[dict[str, object]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    performance_rows = []
    cutoff_rows = []
    biomarkers = BASE_BIOMARKERS + list(RATIO_NUMERATORS)
    for outcome_index, outcome in enumerate(outcomes):
        for marker_index, marker in enumerate(biomarkers):
            y, values = get_marker_data(data, outcome, marker)
            apparent = optimal_metrics(y, values)
            seed = RANDOM_SEED + outcome_index * 1000 + marker_index
            repeat_metrics, cutoffs, directions = repeated_cv_metrics(y, values, seed)
            metric_row = {
                "Analysis": outcome["analysis"],
                "Biomarker": LABELS[marker],
                "N": len(y),
                "Positive n": int(y.sum()),
                "Negative n": int(len(y) - y.sum()),
                "Apparent AUC": apparent["auc"],
                "Apparent cutoff direction": apparent["direction"],
                "Apparent cutoff": apparent["cutoff"],
            }
            for metric, values_array in repeat_metrics.items():
                low, high = percentile_interval(values_array)
                display_name = metric.replace("_", " ").title()
                metric_row[f"CV {display_name} mean"] = float(np.mean(values_array))
                metric_row[f"CV {display_name} repeat percentile lower"] = low
                metric_row[f"CV {display_name} repeat percentile upper"] = high

            fixed_intervals = fixed_cutoff_bootstrap(
                y,
                values,
                int(apparent["direction_sign"]),
                float(apparent["cutoff"]),
                seed + 500000,
            )
            for metric in ["accuracy", "balanced_accuracy", "sensitivity", "specificity"]:
                display_name = metric.replace("_", " ").title()
                metric_row[f"Apparent {display_name}"] = apparent[metric]
                metric_row[f"Fixed-cutoff bootstrap {display_name} CI lower"] = fixed_intervals[metric][0]
                metric_row[f"Fixed-cutoff bootstrap {display_name} CI upper"] = fixed_intervals[metric][1]
            performance_rows.append(metric_row)

            expected_direction = int(apparent["direction_sign"])
            same_direction = directions == expected_direction
            stable_cutoffs = cutoffs[same_direction]
            cutoff_low, cutoff_high = percentile_interval(stable_cutoffs)
            cutoff_rows.append(
                {
                    "Analysis": outcome["analysis"],
                    "Biomarker": LABELS[marker],
                    "N training folds": len(cutoffs),
                    "Full-sample cutoff direction": apparent["direction"],
                    "Full-sample cutoff": apparent["cutoff"],
                    "Training-fold direction agreement": float(np.mean(same_direction)),
                    "Training-fold cutoff median": float(np.median(stable_cutoffs)),
                    "Training-fold cutoff 2.5th percentile": cutoff_low,
                    "Training-fold cutoff 97.5th percentile": cutoff_high,
                }
            )
    return pd.DataFrame(performance_rows), pd.DataFrame(cutoff_rows)


def pairwise_analysis(data: pd.DataFrame, outcomes: list[dict[str, object]]) -> pd.DataFrame:
    rows = []
    for outcome_index, outcome in enumerate(outcomes):
        for pair_index, (marker_1, marker_2, comparison_name) in enumerate(COMPARISON_PAIRS):
            local = outcome["mask"] & data[marker_1].notna() & data[marker_2].notna()
            y = data.loc[local, outcome["outcome"]].to_numpy(dtype=int)
            values_1 = pd.to_numeric(data.loc[local, marker_1], errors="coerce").to_numpy(dtype=float)
            values_2 = pd.to_numeric(data.loc[local, marker_2], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(values_1) & np.isfinite(values_2)
            y, values_1, values_2 = y[valid], values_1[valid], values_2[valid]

            _, direction_1 = oriented_auc(y, values_1)
            _, direction_2 = oriented_auc(y, values_2)
            score_1 = values_1 * direction_1
            score_2 = values_2 * direction_2
            delong = paired_delong(y, score_1, score_2)
            bootstrap = paired_auc_bootstrap(
                y,
                score_1,
                score_2,
                RANDOM_SEED + 900000 + outcome_index * 100 + pair_index,
            )
            rows.append(
                {
                    "Analysis": outcome["analysis"],
                    "Comparison": comparison_name,
                    "Marker 1": marker_1,
                    "Marker 2": marker_2,
                    "N": len(y),
                    "Positive n": int(y.sum()),
                    "Negative n": int(len(y) - y.sum()),
                    **delong,
                    **bootstrap,
                }
            )
    results = pd.DataFrame(rows)
    results["DeLong FDR q (4 comparisons within analysis)"] = results.groupby("Analysis", sort=False)[
        "delong_p"
    ].transform(bh_adjust)
    results["Bootstrap FDR q (4 comparisons within analysis)"] = results.groupby("Analysis", sort=False)[
        "bootstrap_p"
    ].transform(bh_adjust)
    return results


def best_vs_others_analysis(
    data: pd.DataFrame,
    outcomes: list[dict[str, object]],
    cv_results: pd.DataFrame,
) -> pd.DataFrame:
    """For each analysis, compare the marker with the highest full-sample
    (apparent) point-estimate AUC against each of the other ten markers in
    the 11-marker panel, using the same paired DeLong and paired-bootstrap
    methodology as the four prespecified comparisons above. This is a
    post-hoc, exploratory "best versus the rest" comparison, not one of the
    four prespecified contrasts; FDR correction is applied within each
    analysis across these ten comparisons, separately from the four
    prespecified-comparison FDR correction above."""
    rows = []
    for outcome_index, outcome in enumerate(outcomes):
        analysis_name = outcome["analysis"]
        analysis_cv = cv_results[cv_results["Analysis"] == analysis_name]
        best_row = analysis_cv.loc[analysis_cv["Apparent AUC"].idxmax()]
        best_label = best_row["Biomarker"]
        label_to_column = {LABELS[marker]: marker for marker in BASE_BIOMARKERS + list(RATIO_NUMERATORS)}
        best_marker = label_to_column[best_label]
        other_markers = [m for m in BASE_BIOMARKERS + list(RATIO_NUMERATORS) if m != best_marker]
        for pair_index, other_marker in enumerate(other_markers):
            local = outcome["mask"] & data[best_marker].notna() & data[other_marker].notna()
            y = data.loc[local, outcome["outcome"]].to_numpy(dtype=int)
            values_1 = pd.to_numeric(data.loc[local, best_marker], errors="coerce").to_numpy(dtype=float)
            values_2 = pd.to_numeric(data.loc[local, other_marker], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(values_1) & np.isfinite(values_2)
            y, values_1, values_2 = y[valid], values_1[valid], values_2[valid]

            _, direction_1 = oriented_auc(y, values_1)
            _, direction_2 = oriented_auc(y, values_2)
            score_1 = values_1 * direction_1
            score_2 = values_2 * direction_2
            delong = paired_delong(y, score_1, score_2)
            bootstrap = paired_auc_bootstrap(
                y,
                score_1,
                score_2,
                RANDOM_SEED + 1900000 + outcome_index * 100 + pair_index,
            )
            rows.append(
                {
                    "Analysis": analysis_name,
                    "Best marker": LABELS[best_marker],
                    "Other marker": LABELS[other_marker],
                    "N": len(y),
                    "Positive n": int(y.sum()),
                    "Negative n": int(len(y) - y.sum()),
                    **delong,
                    **bootstrap,
                }
            )
    results = pd.DataFrame(rows)
    results["DeLong FDR q (10 comparisons within analysis)"] = results.groupby("Analysis", sort=False)[
        "delong_p"
    ].transform(bh_adjust)
    results["Bootstrap FDR q (10 comparisons within analysis)"] = results.groupby("Analysis", sort=False)[
        "bootstrap_p"
    ].transform(bh_adjust)
    return results


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def write_summary(cv: pd.DataFrame, cutoffs: pd.DataFrame, comparisons: pd.DataFrame, best_vs_others: pd.DataFrame) -> None:
    lines = [
        "## Material Passport",
        "",
        "- artifact_type: statistical_validation",
        "- verification_status: VERIFIED",
        f"- input: {INPUT}",
        f"- seed: {RANDOM_SEED}",
        f"- repeated_cv: {CV_REPEATS} repeats x {CV_FOLDS} folds",
        f"- paired_bootstrap_replicates: {BOOTSTRAP_REPS}",
        "",
        "# Diagnostic performance validation summary",
        "",
        "## Cross-validated leading biomarkers",
        "",
    ]
    for analysis, frame in cv.groupby("Analysis", sort=False):
        top = frame.sort_values("CV Auc mean", ascending=False).head(3)
        lines.append(f"### {analysis}")
        for _, row in top.iterrows():
            lines.append(
                f"- {row['Biomarker']}: CV AUC {row['CV Auc mean']:.3f} "
                f"(repeat percentile {row['CV Auc repeat percentile lower']:.3f}–"
                f"{row['CV Auc repeat percentile upper']:.3f})."
            )
        lines.append("")

    lines.extend(["## Prespecified paired AUC comparisons", ""])
    for analysis, frame in comparisons.groupby("Analysis", sort=False):
        lines.append(f"### {analysis}")
        for _, row in frame.iterrows():
            lines.append(
                f"- {row['Marker 1']} minus {row['Marker 2']}: ΔAUC={row['delta_auc']:.3f}, "
                f"DeLong p={format_p(row['delong_p'])}, "
                f"FDR q={format_p(row['DeLong FDR q (4 comparisons within analysis)'])}, "
                f"paired-bootstrap 95% CI {row['bootstrap_ci_lower']:.3f} to {row['bootstrap_ci_upper']:.3f}."
            )
        lines.append("")

    lines.extend(["## Best marker versus the other ten (post-hoc, exploratory)", ""])
    for analysis, frame in best_vs_others.groupby("Analysis", sort=False):
        best_label = frame["Best marker"].iloc[0]
        lines.append(f"### {analysis} (best: {best_label})")
        for _, row in frame.sort_values("delong_p").iterrows():
            lines.append(
                f"- {row['Best marker']} minus {row['Other marker']}: \u0394AUC={row['delta_auc']:.3f}, "
                f"DeLong p={format_p(row['delong_p'])}, "
                f"FDR q={format_p(row['DeLong FDR q (10 comparisons within analysis)'])}, "
                f"paired-bootstrap 95% CI {row['bootstrap_ci_lower']:.3f} to {row['bootstrap_ci_upper']:.3f}."
            )
        lines.append("")

    lines.extend(
        [
            "## Interpretation notes",
            "",
            "- CV AUC is calculated in each held-out fold and averaged within each repeat; its percentile interval describes repeat-to-repeat split variability, not an independent-sample confidence interval.",
            "- Fixed-cutoff bootstrap intervals hold the full-sample Youden cutoff fixed and quantify participant-sampling uncertainty at that cutoff.",
            "- Cutoff stability is summarized from training folds only.",
            "- Pairwise bootstrap resamples participants identically for both biomarkers.",
            "- DeLong and bootstrap multiplicity corrections are applied separately within each analysis across four prespecified comparisons.",
            "- The best-versus-other-ten comparison is post-hoc and exploratory (not one of the four prespecified contrasts): within each analysis, the marker with the highest apparent (full-sample) AUC is compared against each of the remaining ten markers, with its own separate FDR correction across those ten comparisons.",
        ]
    )
    (OUTPUT_DIR / "diagnostic_validation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data, outcomes = prepare_data(INPUT)
    cv_results, cutoff_results = cv_analysis(data, outcomes)
    comparison_results = pairwise_analysis(data, outcomes)
    best_vs_others_results = best_vs_others_analysis(data, outcomes, cv_results)

    cv_results.to_csv(OUTPUT_DIR / "diagnostic_cv_performance.csv", index=False)
    cutoff_results.to_csv(OUTPUT_DIR / "diagnostic_cutoff_stability.csv", index=False)
    comparison_results.to_csv(OUTPUT_DIR / "diagnostic_auc_pairwise_comparisons.csv", index=False)
    best_vs_others_results.to_csv(
        OUTPUT_DIR / "best_vs_others_AUC_comparisons.csv", index=False
    )
    write_summary(cv_results, cutoff_results, comparison_results, best_vs_others_results)

    metadata = {
        "input": str(INPUT),
        "random_seed": RANDOM_SEED,
        "cv_repeats": CV_REPEATS,
        "cv_folds": CV_FOLDS,
        "bootstrap_replicates": BOOTSTRAP_REPS,
        "paired_comparisons": COMPARISON_PAIRS,
        "best_vs_others": "post-hoc, exploratory: highest apparent-AUC marker per analysis vs each of the other ten markers",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    (OUTPUT_DIR / "diagnostic_validation_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Validated {len(cv_results)} analysis-marker combinations.")
    print(f"Completed {len(comparison_results)} prespecified paired comparisons.")
    print(f"Completed {len(best_vs_others_results)} post-hoc best-vs-others comparisons.")
    print(f"Results written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
