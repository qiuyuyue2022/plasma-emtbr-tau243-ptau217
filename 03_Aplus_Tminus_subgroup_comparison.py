"""Compare data-driven subgroups within the A+T- participants, on both the
biomarkers that define Table 4 and the contemporaneous CSF Lumipulse data
that feeds Figure 4.

Part 1 fits a two-component Gaussian mixture model to log10(MTBR-Tau243),
the biomarker with evidence of bimodality, to define Low/High MTBR-Tau243
subgroups, then compares clinical/demographic variables and the six
secondary plasma biomarkers between those subgroups. p-Tau217 is an
external comparison variable and does not contribute to subgroup
assignment. Component numbers are converted to stable, interpretable
Low/High labels, so the output does not depend on scikit-learn's arbitrary
component IDs.

Part 2 reuses those same Low/High subgroup assignments -- the GMM is
deliberately not refitted -- to compare contemporaneous CSF Lumipulse
biomarkers in the smaller subset of A+T- participants who have CSF data.

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
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs")))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
N_PERMUTATIONS = int(os.environ.get("BIOMARKER_SUBGROUP_PERMUTATIONS", "10000"))
MIN_N_PER_GROUP = 3
CSF_T_PTAU181_CUTOFF = 61.6

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

CSF_VARIABLES = [
    ("CSF Aβ42/Aβ40", "CSF_Lumipulse_Aβ42/Aβ40"),
    ("CSF p-tau181", "CSF_Lumipulse_p-tau"),
    ("CSF t-tau", "CSF_Lumipulse_t-tau"),
    ("CSF p-tau181/Aβ42", "CSF_Lumipulse_p-tau/Aβ42"),
    ("CSF t-tau/Aβ42", "CSF_Lumipulse_t-tau/Aβ42"),
]


# ---------------------------------------------------------------------------
# Shared helpers (used by both the subgroup-definition part and the CSF part)
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
    subset = data.loc[mask].copy()
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


# ---------------------------------------------------------------------------
# Part 2: contemporaneous CSF Lumipulse comparison in the same subgroups; feeds Figure 4
# ---------------------------------------------------------------------------

def exact_permutation_mann_whitney(
    high: np.ndarray, low: np.ndarray
) -> tuple[float, float, int]:
    """Tie-aware exact randomization p-value for the Mann-Whitney statistic."""
    pooled = np.concatenate([high, low])
    n_high = len(high)
    n_low = len(low)
    ranks = rankdata(pooled, method="average")
    observed_u = float(ranks[:n_high].sum() - n_high * (n_high + 1) / 2)
    expected_u = n_high * n_low / 2
    observed_distance = abs(observed_u - expected_u)

    extreme = 0
    total = 0
    for high_indices in combinations(range(len(pooled)), n_high):
        u_value = float(ranks[list(high_indices)].sum() - n_high * (n_high + 1) / 2)
        if abs(u_value - expected_u) >= observed_distance - 1e-12:
            extreme += 1
        total += 1
    return observed_u, extreme / total, total


def median_iqr_csf(values: np.ndarray) -> str:
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return f"{median:.4g} [{q1:.4g}, {q3:.4g}]"


def prepare_csf_data(membership: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the raw workbook again for the CSF-specific columns and merge in
    the Low/High subgroup assignment already computed in Part 1 (`membership`
    has SID + Plasma tau subgroup; the GMM is not refitted here)."""
    data = pd.read_excel(INPUT, sheet_name="Sheet1")
    required = {
        "SID", "A", "T", "method", "taupet", "CSF_Lumipulse_采样时间",
        "CSF_Lumipulse_绝对时间差_天", *[column for _, column in CSF_VARIABLES],
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Missing required columns in {INPUT}: {missing}")
    if not {"SID", "Plasma tau subgroup"}.issubset(membership.columns):
        raise KeyError("The locked membership table lacks SID or Plasma tau subgroup.")
    if membership["SID"].duplicated().any():
        raise ValueError("The locked membership table contains duplicate SID values.")

    mask = clean_status(data["A"]).eq("A+") & clean_status(data["T"]).eq("T-")
    all_at = data.loc[mask].copy()
    all_at = all_at.merge(
        membership[["SID", "Plasma tau subgroup"]],
        on="SID",
        how="left",
        validate="one_to_one",
    )
    if len(all_at) != 33:
        raise ValueError(f"Expected 33 A+T- participants; found {len(all_at)}.")
    if all_at["Plasma tau subgroup"].isna().any():
        raise ValueError("Some A+T- participants lack locked subgroup assignments.")

    for _, column in CSF_VARIABLES:
        all_at[column] = pd.to_numeric(all_at[column], errors="coerce")
    all_at["CSF_Lumipulse_绝对时间差_天"] = pd.to_numeric(
        all_at["CSF_Lumipulse_绝对时间差_天"], errors="coerce"
    )
    csf_subset = all_at.loc[all_at["CSF_Lumipulse_采样时间"].notna()].copy()
    if len(csf_subset) != 14:
        raise ValueError(f"Expected 14 contemporaneous CSF matches; found {len(csf_subset)}.")
    return all_at, csf_subset


def compare_csf(csf_subset: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for label, column in CSF_VARIABLES:
        low = csf_subset.loc[
            csf_subset["Plasma tau subgroup"] == LOW, column
        ].dropna().to_numpy(dtype=float)
        high = csf_subset.loc[
            csf_subset["Plasma tau subgroup"] == HIGH, column
        ].dropna().to_numpy(dtype=float)
        if len(low) == 0 or len(high) == 0:
            rows.append(
                {
                    "Variable": label,
                    "Source column": column,
                    "Low n": len(low),
                    "High n": len(high),
                    "Test status": "Not tested: one group has no observations",
                }
            )
            continue

        base = {
            "Variable": label,
            "Source column": column,
            "Low n": len(low),
            "Low median [IQR]": median_iqr_csf(low),
            "High n": len(high),
            "High median [IQR]": median_iqr_csf(high),
        }
        if len(low) < MIN_N_PER_GROUP or len(high) < MIN_N_PER_GROUP:
            rows.append(
                {
                    **base,
                    "Test status": f"Descriptive only: fewer than {MIN_N_PER_GROUP} observations in a group",
                }
            )
            continue

        u_statistic, p_value, permutations = exact_permutation_mann_whitney(high, low)
        rank_biserial = 2 * u_statistic / (len(high) * len(low)) - 1
        fc_estimate, fc_low, fc_high = fold_change_estimate(high, low)
        rows.append(
            {
                **base,
                "Fold-change (High / Low), Hodges-Lehmann": fc_estimate,
                "Fold-change 95% CI lower": fc_low,
                "Fold-change 95% CI upper": fc_high,
                "Mann-Whitney U (High group)": u_statistic,
                "Exact permutation count": permutations,
                "Raw p": p_value,
                "Rank-biserial r (positive means higher in High group)": rank_biserial,
                "Test status": "Tested",
            }
        )

    result = pd.DataFrame(rows)
    result["FDR q (evaluable CSF comparisons)"] = bh_adjust(result["Raw p"])
    return result


def build_coverage(all_at: pd.DataFrame, csf_subset: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in [LOW, HIGH]:
        full_n = int((all_at["Plasma tau subgroup"] == group).sum())
        matched = csf_subset.loc[csf_subset["Plasma tau subgroup"] == group]
        rows.append(
            {
                "Plasma tau subgroup": group,
                "A+T- total n": full_n,
                "Contemporaneous CSF n": len(matched),
                "CSF availability (%)": 100 * len(matched) / full_n,
                "Absolute blood-CSF interval, median [range], days": (
                    f"{matched['CSF_Lumipulse_绝对时间差_天'].median():.1f} "
                    f"[{matched['CSF_Lumipulse_绝对时间差_天'].min():.0f}, "
                    f"{matched['CSF_Lumipulse_绝对时间差_天'].max():.0f}]"
                ),
                "Tau PET reference n": int((matched["taupet"] == 1).sum()),
                "Non-Tau-PET reference n": int((matched["taupet"] != 1).sum()),
            }
        )
    return pd.DataFrame(rows)


def availability_fisher(all_at: pd.DataFrame) -> tuple[float, float]:
    table = []
    for group in [LOW, HIGH]:
        group_data = all_at.loc[all_at["Plasma tau subgroup"] == group]
        available = int(group_data["CSF_Lumipulse_采样时间"].notna().sum())
        table.append([available, len(group_data) - available])
    odds_ratio, p_value = fisher_exact(np.asarray(table), alternative="two-sided")
    return float(odds_ratio), float(p_value)


def write_csf_report(
    results: pd.DataFrame,
    coverage: pd.DataFrame,
    availability_or: float,
    availability_p: float,
) -> None:
    lines = [
        "## Material Passport",
        "",
        "- Origin Skill: experiment-agent",
        "- Origin Mode: validate",
        "- Origin Date: 2026-10-29",
        "- Verification Status: ANALYZED",
        "- Version Label: validation_v1",
        "",
        "# A+T− low versus high plasma eMTBR-tau243: CSF Lumipulse comparison",
        "",
        "## Analysis set and method",
        "",
        "- Subgroups were inherited from the locked 33-participant log10(eMTBR-tau243) GMM; the model was not refitted in the CSF subset.",
        "- Contemporaneous CSF was defined by the existing merge rule: nearest same-ID Lumipulse sample within 30 days of blood collection.",
        "- Continuous outcomes used a tie-aware exact permutation Mann-Whitney test, rank-biserial effect size, and "
        "a log10-scale Hodges-Lehmann estimate back-transformed to a fold-change (High / Low) with its 95% CI, "
        "appropriate for these right-skewed concentration measures.",
        f"- Benjamini-Hochberg FDR correction was applied jointly to all {int(results['Raw p'].notna().sum())} CSF comparisons.",
        "",
        "## Coverage",
        "",
        coverage.to_markdown(index=False),
        "",
        f"CSF availability comparison: Fisher exact odds ratio={availability_or:.3g}, p={availability_p:.3f}.",
        "",
        "## Statistical findings",
        "",
        "| Variable | Low median [IQR] | High median [IQR] | Fold-change, High / Low (95% CI) | Exact p | FDR q | Rank-biserial r |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in results.iterrows():
        lines.append(
            f"| {row['Variable']} | {row['Low median [IQR]']} (n={int(row['Low n'])}) | "
            f"{row['High median [IQR]']} (n={int(row['High n'])}) | "
            f"{format_number(row.get('Fold-change (High / Low), Hodges-Lehmann'))}x "
            f"({format_number(row.get('Fold-change 95% CI lower'))}, {format_number(row.get('Fold-change 95% CI upper'))}) | "
            f"{format_p(row.get('Raw p'))} | {format_p(row.get('FDR q (evaluable CSF comparisons)'))} | "
            f"{format_number(row.get('Rank-biserial r (positive means higher in High group)'), 2)} |"
        )

    lines.extend(
        [
            "",
            "## Warnings",
            "",
            "- Overall confidence: CAUTION. Only 14 participants had contemporaneous CSF (6 low, 8 high), and some analytes had additional missingness.",
            "- The T-reference composition was imbalanced: 5/6 low-group participants but only 1/8 high-group participants had Tau PET available.",
            f"- For participants whose T− classification used CSF p-tau181 ≤{CSF_T_PTAU181_CUTOFF:g}, comparison of absolute CSF p-tau181 is conditioned on the classification cutoff and is not an independent validation.",
            "- CSF p-tau181/Aβ42 and t-tau/Aβ42 were not used to define T status, but their subgroup differences remain exploratory and may reflect the Aβ42 denominator as well as the tau numerator.",
            "- The analysis is exploratory and cross-sectional. It cannot establish that plasma eMTBR-tau243 changes before CSF or PET tau pathology.",
            "",
            "## Fallacy scan",
            "",
            "- Coverage: 11/11 checked.",
            "- Simpson's paradox: CAUTION; reference-standard composition differs between subgroups and cannot be meaningfully stratified at this sample size.",
            "- Ecological fallacy: not applicable; participant-level data were analyzed.",
            "- Berkson's paradox: CAUTION; analysis is restricted to A+T− participants with contemporaneous CSF.",
            "- Collider bias: CAUTION; conditioning on A+T− status may alter biomarker associations.",
            "- Base-rate neglect: not applicable; no diagnostic predictive values are reported.",
            "- Regression to the mean: no repeated-change claim was made.",
            "- Survivorship/availability bias: CAUTION; only 14/33 participants contributed CSF data.",
            f"- Look-elsewhere effect: addressed with BH-FDR across all {int(results['Raw p'].notna().sum())} CSF tests.",
            "- Garden of forking paths: CAUTION; exploratory outcome family and 30-day window should be locked before manuscript reporting.",
            "- Correlation does not imply causation: CAUTION; group differences are associational.",
            "- Reverse causality: CAUTION; cross-sectional timing cannot establish temporal order.",
            "",
            "## Reproducibility",
            "",
            "- Method: deterministic group membership, exhaustive label permutations, and a closed-form Hodges-Lehmann interval.",
            "- Verdict: CANNOT_VERIFY from this single execution; the statistical output is ANALYZED.",
        ]
    )
    (OUTPUT_DIR / "A+T-_CSF_Lumipulse_validation.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def write_manuscript_table_s6(results: pd.DataFrame, coverage: pd.DataFrame) -> None:
    """Format the same CSF comparison as a manuscript-ready Supplementary
    Table S6, matching the TableS5/TableS1 markdown-table style, rather than
    the raw analysis-passport report written by write_csf_report(). This was
    added 2026-11-01 at author direction: the manuscript text (Results
    Section 3.7) and 00_story_outline_and_analysis_lock.md had both already
    referenced "Supplementary Table S5" as available, but no formatted table
    file had actually been generated -- only the passport-style
    A+T-_CSF_Lumipulse_validation.md existed. Renumbered S5 -> S6 later the
    same day to match true first-citation order (changelog 3.10)."""
    low_n = int(coverage.loc[coverage["Plasma tau subgroup"] == LOW, "Contemporaneous CSF n"].iloc[0])
    high_n = int(coverage.loc[coverage["Plasma tau subgroup"] == HIGH, "Contemporaneous CSF n"].iloc[0])
    lines = [
        "# Supplementary Table S8. Contemporaneous CSF comparison between eMTBR-tau243-defined Low and High A+T− groups",
        "",
        f"| Characteristic | Low (n={low_n}) | High (n={high_n}) | Fold-change, High/Low (95% CI) | Exact p | FDR q | Rank-biserial r |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in results.iterrows():
        if row.get("Test status") != "Tested":
            continue
        fc = format_number(row.get("Fold-change (High / Low), Hodges-Lehmann"), 3)
        fc_lo = format_number(row.get("Fold-change 95% CI lower"), 2)
        fc_hi = format_number(row.get("Fold-change 95% CI upper"), 2)
        r_value = row.get("Rank-biserial r (positive means higher in High group)")
        lines.append(
            f"| {row['Variable']} | {row['Low median [IQR]']} | {row['High median [IQR]']} | "
            f"{fc}x ({fc_lo}, {fc_hi}) | {format_p(row.get('Raw p'))} | "
            f"{format_p(row.get('FDR q (evaluable CSF comparisons)'))} | {r_value:.2f} |"
        )
    lines.append("")
    lines.append(
        "Subgroup membership (Low/High eMTBR-tau243) is inherited from the full-sample 33-participant "
        "Gaussian-mixture model (Supplementary Table S7) and was not refit in this CSF subset. Contemporaneous "
        "CSF was defined as the nearest same-participant LUMIPULSE sample within 30 days of blood collection; "
        f"only {low_n + high_n} of 33 A+T− participants ({low_n} Low, {high_n} High) had contemporaneous CSF, "
        "and CSF availability did not differ significantly between subgroups (Fisher exact OR=1.62, p=0.716). "
        "The T-reference composition of this CSF subset was imbalanced (5/6 Low participants but only 1/8 High "
        "participants had tau PET), so these ratio findings are exploratory profile differences rather than "
        "independent validation of the plasma-defined groups. CSF p-tau181/Aβ42 and t-tau/Aβ42 were not used to "
        "define T status; for CSF-classified participants, the operational T rule used absolute CSF p-tau181 at "
        "a cutoff of 61.6, so comparisons of absolute CSF p-tau181 are partly conditioned on that cutoff. "
        "Benjamini-Hochberg FDR correction was applied jointly across these five comparisons.\n"
    )
    (OUTPUT_DIR / "TableS8_Aplus_Tminus_CSF_comparison.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Entry point: Part 1 (subgroup definition, feeds Table 4), then Part 2
# (CSF comparison in the same subgroups, feeds Figure 4) -- in one pass, no
# round trip through a membership CSV between the two parts.
# ---------------------------------------------------------------------------

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

    # --- Part 2: compare contemporaneous CSF Lumipulse biomarkers in the ---
    # --- same subgroups (the GMM above is not refitted here) ---
    all_at, csf_subset = prepare_csf_data(membership)
    csf_results = compare_csf(csf_subset)
    coverage = build_coverage(all_at, csf_subset)
    availability_or, availability_p = availability_fisher(all_at)

    csf_results.to_csv(OUTPUT_DIR / "A+T-_CSF_Lumipulse_comparison.csv", index=False)
    coverage.to_csv(OUTPUT_DIR / "A+T-_CSF_Lumipulse_coverage.csv", index=False)
    write_csf_report(csf_results, coverage, availability_or, availability_p)
    write_manuscript_table_s6(csf_results, coverage)

    print(f"\nA+T- participants with contemporaneous CSF: {len(csf_subset)}")
    print(csf_subset["Plasma tau subgroup"].value_counts(sort=False).to_string())
    print("\nCSF comparison:")
    print(
        csf_results[
            [
                "Variable", "Low n", "Low median [IQR]", "High n",
                "High median [IQR]", "Raw p", "FDR q (evaluable CSF comparisons)",
                "Rank-biserial r (positive means higher in High group)",
            ]
        ].to_string(index=False)
    )
    print(f"\nResults written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
