"""Compare data-driven subgroups within the A+T- participants, on both the
biomarkers (Supplementary Table S8 and Figure 6).

Part 1 fits a two-component Gaussian mixture model to log10(MTBR-Tau243),
the biomarker with evidence of bimodality, to define Low/High MTBR-Tau243
subgroups, then compares clinical/demographic variables and the six
secondary plasma biomarkers between those subgroups. p-Tau217 is an
external comparison variable and does not contribute to subgroup
assignment. Component numbers are converted to stable, interpretable
Low/High labels, so the output does not depend on scikit-learn's arbitrary
component IDs.

v3.60 (2026-11-11): restricted to the 18 tau PET-defined A+T− patients; the
pooled (n=33) mode and the contemporaneous CSF comparison (Part 2; only 6 of
18 had CSF, 5 Low and 1 High, not evaluable) were removed.

Outputs are deliberately de-identified: the participant-level membership
file contains SID but not name or hospital ID.
"""

from __future__ import annotations

from itertools import combinations
import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu, norm, rankdata
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = PROJECT_DIR / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs" / "Aplus_Tminus_tauPET_only")))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
N_PERMUTATIONS = int(os.environ.get("BIOMARKER_SUBGROUP_PERMUTATIONS", "10000"))
MIN_N_PER_GROUP = 3

LOW = "Low MTBR-Tau243"
HIGH = "High MTBR-Tau243"

COL = {
    "abeta42": "Aβ 1-42",
    "ptau217": "p-Tau217(≤0.758)",
    "ptau231": "p-Tau231(≤5.793)",
    "bd_tau": "BD Tau（≤52.301）",
    "mtbr243": "MTBR-Tau243",
    "gfap": "GFAP(≤131.880)",
    "nfl": "NFL(≤9.624 )",
}

BIOMARKERS = [
    "p-Tau217",
    "BD Tau",
    "p-Tau231",
    "GFAP",
    "Aβ1-42",
    "NFL",
]



# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def clean_status(series: pd.Series) -> pd.Series:
    """Normalize ASCII and Unicode minus signs in A/T status columns."""
    return (
        series.astype("string")
        .str.strip()
        .str.replace("−", "-", regex=False)
        .str.replace("–", "-", regex=False)
    )


def bh_adjust(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p-values, preserving missing values."""
    result = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    sorted_p = valid.to_numpy()[order]
    n = len(sorted_p)
    adjusted = np.minimum.accumulate((sorted_p * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    original_positions = valid.index.to_numpy()[order]
    result.loc[original_positions] = adjusted
    return result


def hodges_lehmann_estimate(high: np.ndarray, low: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Hodges-Lehmann estimate of the location shift (High minus Low): the
    median of all pairwise differences, with the Mann-Whitney-based
    confidence interval (tie-corrected normal approximation; Hollander &
    Wolfe, "Nonparametric Statistical Methods", and Conover, "Practical
    Nonparametric Statistics"). This is the estimator and interval that
    pairs with the Mann-Whitney U test used for the p-value, rather than a
    bootstrap interval around the difference of the two sample medians."""
    m, n = len(high), len(low)
    diffs = np.subtract.outer(high, low).ravel()
    diffs.sort()
    n_pairs = m * n
    estimate = float(np.median(diffs))

    pooled = np.concatenate([high, low])
    ranks = rankdata(pooled, method="average")
    _, tie_counts = np.unique(ranks, return_counts=True)
    total_n = m + n
    tie_term = np.sum(tie_counts**3 - tie_counts) / (total_n * (total_n - 1)) if total_n > 1 else 0.0
    variance = m * n / 12 * ((total_n + 1) - tie_term)
    if variance <= 0:
        return estimate, estimate, estimate

    z = norm.ppf(1 - alpha / 2)
    k = int(np.floor(n_pairs / 2 - z * np.sqrt(variance) - 0.5))
    k = min(max(k, 0), n_pairs - 1)
    ci_lower = float(diffs[k])
    ci_upper = float(diffs[n_pairs - 1 - k])
    return estimate, ci_lower, ci_upper


def fold_change_estimate(high: np.ndarray, low: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Hodges-Lehmann location shift computed on log10-transformed values,
    then back-transformed to a fold-change (High / Low) with its 95% CI.
    Appropriate for concentration-type biomarkers, which are right-skewed
    on the raw scale; median(log(X)) = log(median(X)), so this exponentiated
    estimate is exactly the ratio of the two groups' typical levels, not an
    approximation. Does not change the Mann-Whitney p-value or rank-biserial
    r, which are already invariant to monotonic transforms such as log."""
    if np.any(high <= 0) or np.any(low <= 0):
        raise ValueError("fold_change_estimate requires strictly positive values for a log10 transform.")
    estimate, ci_lower, ci_upper = hodges_lehmann_estimate(np.log10(high), np.log10(low), alpha)
    return 10.0**estimate, 10.0**ci_lower, 10.0**ci_upper


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def format_number(value: float, digits: int = 4) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}g}"


# ---------------------------------------------------------------------------
# Part 1: define Low/High MTBR-Tau243 subgroups; feeds Table 4
# ---------------------------------------------------------------------------

def prepare_data(path: Path) -> pd.DataFrame:
    data = pd.read_excel(path, sheet_name="Sheet1")
    required = {
        "SID", "A", "T", "AGE", "SEX", "MMSE_Total", "MOCA_P", "apoe",
        *COL.values(),
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Required columns are missing from the workbook: {missing}")

    for source in COL.values():
        data[source] = pd.to_numeric(data[source], errors="coerce")
    data["MTBR-Tau243"] = data[COL["mtbr243"]]
    data["p-Tau217"] = data[COL["ptau217"]]
    data["BD Tau"] = data[COL["bd_tau"]]
    data["p-Tau231"] = data[COL["ptau231"]]
    data["GFAP"] = data[COL["gfap"]]
    data["Aβ1-42"] = data[COL["abeta42"]]
    data["NFL"] = data[COL["nfl"]]

    mask = clean_status(data["A"]).eq("A+") & clean_status(data["T"]).eq("T-")
    mask &= pd.to_numeric(data["taupet"], errors="coerce").eq(1)  # tau PET-defined T- only
    subset = data.loc[mask].copy()
    if len(subset) != 18:
        raise ValueError(f"Expected 18 tau PET-defined A+T- patients; found {len(subset)}.")
    if len(subset) < 10:
        raise ValueError(f"Only {len(subset)} A+T- participants were found; check A/T coding.")

    apoe_text = subset["apoe"].astype("string").str.lower()
    subset["APOE ε4 carrier"] = np.where(
        subset["apoe"].isna(), pd.NA,
        np.where(apoe_text.str.contains("ε4|e4", regex=True, na=False), "Yes", "No"),
    )
    return subset


def fit_gmm(values: np.ndarray, n_components: int) -> GaussianMixture:
    return GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        n_init=100,
        max_iter=1000,
        reg_covar=1e-6,
        random_state=RANDOM_SEED,
    ).fit(values)


def stable_binary_labels(model: GaussianMixture, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    raw = model.predict(values)
    # Inputs are log10 biomarker concentrations. The component with the larger
    # average log concentration is called High.
    high_component = int(np.argmax(model.means_.mean(axis=1)))
    labels = np.where(raw == high_component, HIGH, LOW)
    probabilities = model.predict_proba(values)[:, high_component]
    return labels, probabilities


def define_subgroups(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_columns = ["MTBR-Tau243"]
    diagnostic_columns = ["MTBR-Tau243", "p-Tau217"]
    valid = data[primary_columns].notna().all(axis=1) & (data[primary_columns] > 0).all(axis=1)
    if not valid.all():
        warnings.warn(f"Excluding {(~valid).sum()} rows with missing/non-positive MTBR-Tau243 values.")
    analysis = data.loc[valid].copy()
    x_primary = np.log10(analysis[primary_columns].to_numpy(dtype=float))

    gmm_1 = fit_gmm(x_primary, 1)
    gmm_2 = fit_gmm(x_primary, 2)
    labels, high_probability = stable_binary_labels(gmm_2, x_primary)
    analysis["Plasma tau subgroup"] = pd.Categorical(
        labels, categories=[LOW, HIGH], ordered=True
    )
    analysis["Probability of high plasma tau subgroup"] = high_probability
    analysis["Subgroup assignment confidence"] = np.maximum(high_probability, 1 - high_probability)

    diagnostics = []
    label_sets = {}
    for model_name, columns in [
        ("MTBR-Tau243 only", ["MTBR-Tau243"]),
        ("p-Tau217 only", ["p-Tau217"]),
        ("MTBR-Tau243 + p-Tau217", diagnostic_columns),
    ]:
        diagnostic_valid = analysis[columns].notna().all(axis=1) & (analysis[columns] > 0).all(axis=1)
        diagnostic_data = analysis.loc[diagnostic_valid]
        x = np.log10(diagnostic_data[columns].to_numpy(dtype=float))
        one = fit_gmm(x, 1)
        two = fit_gmm(x, 2)
        model_labels, probability = stable_binary_labels(two, x)
        label_sets[model_name] = pd.Series(model_labels, index=diagnostic_data.index)
        counts = pd.Series(model_labels).value_counts()
        diagnostics.append(
            {
                "Model": model_name,
                "Variables": " + ".join(columns),
                "N": len(x),
                "BIC, 1 component": one.bic(x),
                "BIC, 2 components": two.bic(x),
                "Delta BIC (1 minus 2; positive favors 2)": one.bic(x) - two.bic(x),
                "AIC, 1 component": one.aic(x),
                "AIC, 2 components": two.aic(x),
                "Delta AIC (1 minus 2; positive favors 2)": one.aic(x) - two.aic(x),
                "Low subgroup n": int(counts.get(LOW, 0)),
                "High subgroup n": int(counts.get(HIGH, 0)),
                "Mean assignment confidence": float(np.mean(np.maximum(probability, 1 - probability))),
            }
        )

    mtbr_labels = label_sets["MTBR-Tau243 only"]
    for i in range(1, len(diagnostics)):
        comparison_labels = label_sets[diagnostics[i]["Model"]]
        diagnostics[i]["Agreement with primary MTBR-only labels"] = float(
            np.mean(mtbr_labels.loc[comparison_labels.index] == comparison_labels)
        )
        diagnostics[i]["Adjusted Rand index vs primary MTBR-only labels"] = float(
            adjusted_rand_score(mtbr_labels.loc[comparison_labels.index], comparison_labels)
        )
    return analysis, pd.DataFrame(diagnostics)


def median_iqr(values: pd.Series) -> str:
    values = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if not len(values):
        return "NA"
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return f"{median:.3g} [{q1:.3g}, {q3:.3g}]"


def continuous_comparisons(data: pd.DataFrame, variables: list[str], family: str, effect_scale: str = "raw") -> pd.DataFrame:
    rows = []
    for variable in variables:
        low = pd.to_numeric(
            data.loc[data["Plasma tau subgroup"] == LOW, variable], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        high = pd.to_numeric(
            data.loc[data["Plasma tau subgroup"] == HIGH, variable], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        if not len(low) or not len(high):
            continue
        test = mannwhitneyu(high, low, alternative="two-sided", method="auto")
        rank_biserial = 2 * float(test.statistic) / (len(high) * len(low)) - 1
        row = {
            "Variable": variable,
            "Family": family,
            "Low n": len(low),
            "Low median [IQR]": median_iqr(pd.Series(low)),
            "High n": len(high),
            "High median [IQR]": median_iqr(pd.Series(high)),
            "Mann-Whitney U": float(test.statistic),
            "Raw p": float(test.pvalue),
            "Rank-biserial r (positive means higher in High subgroup)": rank_biserial,
        }
        if effect_scale == "log":
            fc_estimate, fc_low, fc_high = fold_change_estimate(high, low)
            row["Fold-change (High / Low), Hodges-Lehmann"] = fc_estimate
            row["Fold-change 95% CI lower"] = fc_low
            row["Fold-change 95% CI upper"] = fc_high
        else:
            hl_estimate, hl_low, hl_high = hodges_lehmann_estimate(high, low)
            row["Hodges-Lehmann shift (High minus Low)"] = hl_estimate
            row["Hodges-Lehmann 95% CI lower"] = hl_low
            row["Hodges-Lehmann 95% CI upper"] = hl_high
        rows.append(row)
    return pd.DataFrame(rows)


def pearson_chi2(table: np.ndarray) -> float:
    return float(chi2_contingency(table, correction=False)[0])


def permutation_chi2_pvalue(groups: np.ndarray, categories: np.ndarray) -> tuple[float, float]:
    group_levels = np.array([LOW, HIGH], dtype=object)
    category_levels = pd.unique(categories)

    def make_table(current_groups: np.ndarray) -> np.ndarray:
        return np.array(
            [[np.sum((current_groups == group) & (categories == category)) for category in category_levels]
             for group in group_levels],
            dtype=int,
        )

    observed = make_table(groups)
    observed_chi2 = pearson_chi2(observed)
    rng = np.random.default_rng(RANDOM_SEED)
    at_least_as_extreme = 0
    for _ in range(N_PERMUTATIONS):
        permuted = rng.permutation(groups)
        if pearson_chi2(make_table(permuted)) >= observed_chi2 - 1e-12:
            at_least_as_extreme += 1
    p_value = (at_least_as_extreme + 1) / (N_PERMUTATIONS + 1)
    return observed_chi2, p_value


def categorical_summary(values: pd.Series) -> str:
    observed = values.dropna().astype(str)
    n = len(observed)
    if n == 0:
        return "NA"
    counts = observed.value_counts(sort=False)
    return "; ".join(f"{level}: {count} ({100 * count / n:.1f}%)" for level, count in counts.items())


def categorical_comparisons(data: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    rows = []
    for variable in variables:
        frame = data[["Plasma tau subgroup", variable]].dropna().copy()
        frame[variable] = frame[variable].astype(str)
        if frame[variable].nunique() < 2:
            continue
        table = pd.crosstab(frame["Plasma tau subgroup"], frame[variable]).reindex(
            index=[LOW, HIGH], fill_value=0
        )
        if table.shape == (2, 2):
            odds_ratio, p_value = fisher_exact(table.to_numpy(), alternative="two-sided")
            test_name = "Fisher exact"
            statistic = float(odds_ratio)
        else:
            statistic, p_value = permutation_chi2_pvalue(
                frame["Plasma tau subgroup"].astype(str).to_numpy(),
                frame[variable].to_numpy(),
            )
            test_name = f"Permutation Pearson chi-square ({N_PERMUTATIONS:,} permutations)"
        chi2 = pearson_chi2(table.to_numpy())
        denominator = len(frame) * min(table.shape[0] - 1, table.shape[1] - 1)
        cramers_v = np.sqrt(chi2 / denominator) if denominator > 0 else np.nan
        rows.append(
            {
                "Variable": variable,
                "Low n available": int((frame["Plasma tau subgroup"] == LOW).sum()),
                "Low distribution": categorical_summary(
                    frame.loc[frame["Plasma tau subgroup"] == LOW, variable]
                ),
                "High n available": int((frame["Plasma tau subgroup"] == HIGH).sum()),
                "High distribution": categorical_summary(
                    frame.loc[frame["Plasma tau subgroup"] == HIGH, variable]
                ),
                "Test": test_name,
                "Test statistic (odds ratio for 2x2; chi-square otherwise)": statistic,
                "Raw p": float(p_value),
                "Cramer's V": float(cramers_v),
            }
        )
    return pd.DataFrame(rows)


def apply_clinical_fdr(continuous: pd.DataFrame, categorical: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    combined = pd.concat(
        [
            continuous[["Variable", "Raw p"]].assign(Source="continuous"),
            categorical[["Variable", "Raw p"]].assign(Source="categorical"),
        ],
        ignore_index=True,
    )
    combined["FDR q (all clinical comparisons)"] = bh_adjust(combined["Raw p"])
    q_map = combined.set_index(["Source", "Variable"])["FDR q (all clinical comparisons)"]
    continuous["FDR q (all clinical comparisons)"] = [q_map[("continuous", x)] for x in continuous["Variable"]]
    categorical["FDR q (all clinical comparisons)"] = [q_map[("categorical", x)] for x in categorical["Variable"]]
    return continuous, categorical


def p_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def write_summary(
    data: pd.DataFrame,
    diagnostics: pd.DataFrame,
    clinical_continuous: pd.DataFrame,
    categorical: pd.DataFrame,
    biomarkers: pd.DataFrame,
) -> None:
    low_n = int((data["Plasma tau subgroup"] == LOW).sum())
    high_n = int((data["Plasma tau subgroup"] == HIGH).sum())
    primary = diagnostics.loc[diagnostics["Model"] == "MTBR-Tau243 only"].iloc[0]
    clinical_all = pd.concat(
        [
            clinical_continuous[["Variable", "Raw p", "FDR q (all clinical comparisons)"]],
            categorical[["Variable", "Raw p", "FDR q (all clinical comparisons)"]],
        ], ignore_index=True,
    ).sort_values("Raw p")
    lines = [
        "# A+T− internal subgroup comparison",
        "",
        "## Subgroup definition",
        "",
        "- Two-component Gaussian mixture model on log10(MTBR-Tau243) only.",
        "- p-Tau217 is compared after subgroup assignment and does not define the groups.",
        f"- Low MTBR-Tau243: n={low_n}; High MTBR-Tau243: n={high_n}.",
        f"- Delta BIC (one minus two components): {primary['Delta BIC (1 minus 2; positive favors 2)']:.2f}.",
        f"- Mean posterior assignment confidence: {primary['Mean assignment confidence']:.3f}.",
        "",
        "## Clinical and demographic comparisons",
        "",
    ]
    for _, row in clinical_all.iterrows():
        lines.append(
            f"- {row['Variable']}: raw p={p_text(row['Raw p'])}, "
            f"FDR q={p_text(row['FDR q (all clinical comparisons)'])}."
        )
    lines.extend(
        [
            "",
            "## Secondary biomarker profiles",
            "",
        ]
    )
    for _, row in biomarkers.sort_values("Raw p").iterrows():
        lines.append(
            f"- {row['Variable']}: raw p={p_text(row['Raw p'])}, "
            f"FDR q={p_text(row['FDR q (6 biomarkers)'])}, "
            f"rank-biserial r={row['Rank-biserial r (positive means higher in High subgroup)']:.2f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "- MTBR-Tau243 defines the subgroups; its between-group p-value is descriptive, not independent validation.",
            "- p-Tau217 is a non-defining variable and provides an external comparison with the MTBR-defined groups.",
            "- This is an exploratory, cross-sectional analysis with a small A+T− sample.",
            "- Subgroups must not be interpreted as temporal stages without longitudinal conversion or continuous tau-reference evidence.",
        ]
    )
    (OUTPUT_DIR / "A+T-Subgroup_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    # --- Part 1: define Low/High MTBR-Tau243 subgroups ---
    data = prepare_data(INPUT)
    data, diagnostics = define_subgroups(data)

    clinical_variables = [
        "AGE",
        "MMSE_Total",
        "MOCA_P",
    ]
    categorical_variables = [
        "SEX",
        "APOE ε4 carrier",
    ]

    clinical_continuous = continuous_comparisons(data, clinical_variables, "Clinical/demographic")
    categorical = categorical_comparisons(data, categorical_variables)
    clinical_continuous, categorical = apply_clinical_fdr(clinical_continuous, categorical)

    biomarkers = continuous_comparisons(data, BIOMARKERS, "Biomarker", effect_scale="log")
    biomarkers["FDR q (6 biomarkers)"] = bh_adjust(biomarkers["Raw p"])

    membership_columns = [
        "SID",
        "Plasma tau subgroup",
        "Probability of high plasma tau subgroup",
        "Subgroup assignment confidence",
        "MTBR-Tau243",
        "p-Tau217",
    ]
    membership = data[membership_columns].sort_values(
        ["Plasma tau subgroup", "Probability of high plasma tau subgroup"]
    )

    diagnostics.to_csv(OUTPUT_DIR / "A+T-Subgroup_ModelDiagnostics.csv", index=False)
    clinical_continuous.to_csv(OUTPUT_DIR / "A+T-Subgroup_continuous.csv", index=False)
    categorical.to_csv(OUTPUT_DIR / "A+T-Subgroup_categorical.csv", index=False)
    biomarkers.to_csv(OUTPUT_DIR / "A+T-Subgroup_biomarkers.csv", index=False)
    membership.to_csv(OUTPUT_DIR / "A+T-Subgroup_membership.csv", index=False)
    write_summary(data, diagnostics, clinical_continuous, categorical, biomarkers)

    print(f"Input: {INPUT}")
    print(f"A+T- participants analyzed: {len(data)}")
    print(data["Plasma tau subgroup"].value_counts(sort=False).to_string())
    print(f"Results written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
