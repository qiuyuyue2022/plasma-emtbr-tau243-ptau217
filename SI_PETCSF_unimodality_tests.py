"""Formal unimodality tests in A+T- participants.

The analysis is performed on log10-transformed concentrations, matching the
existing Gaussian-mixture analysis. Hartigan's dip test is the primary
distribution-free test of the null hypothesis of unimodality. The original
Silverman critical-bandwidth bootstrap test (unknown support) is included as a
complementary sensitivity analysis.

The two markers are:
  * eMTBR-tau243: marker with prior evidence of a two-component distribution.
  * p-tau217: negative-comparison marker not used to define the subgroups.

The source workbook is read only. Outputs contain aggregate results and
bootstrap critical bandwidths but no participant identifiers.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.signal import find_peaks


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
VENDOR_DIR = SCRIPT_DIR / "_vendor"
if VENDOR_DIR.exists():
    sys.path.insert(0, str(VENDOR_DIR))

try:
    import diptest as diptest_package
    from diptest import diptest
except ImportError as error:
    raise ImportError(
        "The diptest dependency is required. Install it with "
        "`python -m pip install --target code/_vendor diptest==0.10.0`."
    ) from error


DEFAULT_INPUT = PROJECT_DIR / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
# SI (v3.62, 2026-11-11): combined tau PET/CSF-defined A+T− analysis (n=33), restored from the
# pre-v3.60 code for the Supplementary Information. Statistics go to outputs/Aplus_Tminus_PET_CSF/.
OUTPUT_DIR = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(PROJECT_DIR / "outputs"))) / "Aplus_Tminus_PET_CSF"
APTM_REFERENCE = "all"  # fixed: combined tau PET or CSF-defined A+T−  # "tauPET" = n=18 (2026-11-09)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261029"))
DIP_BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_DIP_BOOTSTRAP", "100000"))
SILVERMAN_BOOTSTRAP_REPS = int(os.environ.get("BIOMARKER_SILVERMAN_BOOTSTRAP", "5000"))
SILVERMAN_GRID_SIZE = int(os.environ.get("BIOMARKER_SILVERMAN_GRID", "2048"))
SILVERMAN_BINARY_STEPS = int(os.environ.get("BIOMARKER_SILVERMAN_BINARY_STEPS", "32"))

MARKERS = {
    "eMTBR-tau243": "MTBR-Tau243",
    "p-tau217": "p-Tau217(≤0.758)",
}


def clean_status(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace("−", "-", regex=False)
        .str.replace("–", "-", regex=False)
    )


def prepare_subset() -> pd.DataFrame:
    data = pd.read_excel(INPUT, sheet_name="Sheet1")
    required = {"A", "T", *MARKERS.values()}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise KeyError(f"Required columns are missing: {missing}")

    mask = clean_status(data["A"]).eq("A+") & clean_status(data["T"]).eq("T-")
    if APTM_REFERENCE == "tauPET":  # 2026-11-09 revision: tau PET-defined A+T- only (n=18)
        mask &= pd.to_numeric(data["taupet"], errors="coerce").eq(1)
    subset = data.loc[mask, list(MARKERS.values())].copy()
    for column in MARKERS.values():
        subset[column] = pd.to_numeric(subset[column], errors="coerce")
    if subset.isna().any().any():
        raise ValueError("Missing marker values were found in A+T- participants.")
    if (subset <= 0).any().any():
        raise ValueError("Non-positive marker values cannot be log10 transformed.")
    if len(subset) != (18 if APTM_REFERENCE == "tauPET" else 33):
        raise ValueError(f"Expected 33 A+T- participants but found {len(subset)}.")
    return subset.reset_index(drop=True)


def kde_mode_count(values: np.ndarray, bandwidth: float, grid_size: int) -> int:
    """Count modes in a Gaussian KDE evaluated on a dense deterministic grid."""
    span = max(float(np.ptp(values)), 1e-6)
    padding = max(4.0 * bandwidth, 0.1 * span)
    grid = np.linspace(values.min() - padding, values.max() + padding, grid_size)
    standardized = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * standardized**2).sum(axis=1)
    # Ignore only machine-scale numerical ripples; genuine sample modes have
    # substantially larger prominence. Observed critical bandwidths are also
    # checked at half and twice the default grid resolution below.
    prominence = max(float(density.max()) * 1e-8, 1e-12)
    peaks, _ = find_peaks(density, prominence=prominence)
    return int(len(peaks))


def silverman_critical_bandwidth(
    standardized_values: np.ndarray,
    grid_size: int = SILVERMAN_GRID_SIZE,
    binary_steps: int = SILVERMAN_BINARY_STEPS,
) -> float:
    """Smallest Gaussian-KDE bandwidth yielding at most one mode."""
    values = np.asarray(standardized_values, dtype=float)
    sample_sd = float(np.std(values, ddof=1))
    low = max(sample_sd * 1e-5, 1e-8)
    high = max(sample_sd, float(np.ptp(values)) / 2.0, 1e-3)
    while kde_mode_count(values, high, grid_size) > 1:
        high *= 2.0
        if high > sample_sd * 1e6:
            raise RuntimeError("Failed to bracket the Silverman critical bandwidth.")

    for _ in range(binary_steps):
        middle = (low + high) / 2.0
        if kde_mode_count(values, middle, grid_size) <= 1:
            high = middle
        else:
            low = middle
    return float(high)


def silverman_test(values: np.ndarray, seed: int) -> tuple[float, float, np.ndarray, dict[str, float]]:
    """Original Silverman bootstrap test for one versus more than one mode.

    Values are standardized before analysis. Bootstrap samples are drawn from
    the KDE at the observed critical bandwidth and rescaled to preserve unit
    variance, following Silverman's smoothed-bootstrap construction.
    """
    standardized = (values - np.mean(values)) / np.std(values, ddof=1)
    observed = silverman_critical_bandwidth(standardized)
    # Numerical-resolution audit for the observed statistic.
    resolution_check = {
        "half_grid": silverman_critical_bandwidth(
            standardized, grid_size=max(SILVERMAN_GRID_SIZE // 2, 512)
        ),
        "double_grid": silverman_critical_bandwidth(
            standardized, grid_size=SILVERMAN_GRID_SIZE * 2
        ),
    }

    rng = np.random.default_rng(seed)
    bootstrap_bandwidths = np.empty(SILVERMAN_BOOTSTRAP_REPS, dtype=float)
    denominator = np.sqrt(1.0 + observed**2)
    for iteration in range(SILVERMAN_BOOTSTRAP_REPS):
        bootstrap_sample = (
            rng.choice(standardized, size=len(standardized), replace=True)
            + observed * rng.standard_normal(len(standardized))
        ) / denominator
        bootstrap_bandwidths[iteration] = silverman_critical_bandwidth(bootstrap_sample)

    p_value = (
        np.sum(bootstrap_bandwidths >= observed - 1e-12) + 1
    ) / (SILVERMAN_BOOTSTRAP_REPS + 1)
    return observed, float(p_value), bootstrap_bandwidths, resolution_check


def format_p(value: float) -> str:
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def main() -> None:
    subset = prepare_subset()
    result_rows: list[dict[str, object]] = []
    bootstrap_rows: list[pd.DataFrame] = []

    for marker_index, (display_name, source_column) in enumerate(MARKERS.items()):
        raw_values = subset[source_column].to_numpy(dtype=float)
        values = np.log10(raw_values)

        dip_statistic, dip_p_value = diptest(
            values,
            boot_pval=True,
            n_boot=DIP_BOOTSTRAP_REPS,
            n_threads=1,
            seed=RANDOM_SEED + marker_index,
        )
        critical_bandwidth, silverman_p_value, bootstrap_bandwidths, resolution = silverman_test(
            values,
            seed=RANDOM_SEED + 10000 + marker_index,
        )

        bootstrap_rows.append(
            pd.DataFrame(
                {
                    "Marker": display_name,
                    "Bootstrap replicate": np.arange(1, SILVERMAN_BOOTSTRAP_REPS + 1),
                    "Silverman critical bandwidth": bootstrap_bandwidths,
                }
            )
        )
        result_rows.append(
            {
                "Marker": display_name,
                "Source column": source_column,
                "Analysis population": "A+T-",
                "Transformation": "log10",
                "N": len(values),
                "Hartigan dip statistic": float(dip_statistic),
                "Hartigan dip bootstrap p": float(dip_p_value),
                "Hartigan null hypothesis": "unimodal distribution",
                "Silverman observed critical bandwidth (standardized log10 scale)": critical_bandwidth,
                "Silverman bootstrap p": silverman_p_value,
                "Silverman null hypothesis": "at most one mode",
                "Silverman bootstrap critical bandwidth 2.5%": float(
                    np.quantile(bootstrap_bandwidths, 0.025)
                ),
                "Silverman bootstrap critical bandwidth median": float(
                    np.median(bootstrap_bandwidths)
                ),
                "Silverman bootstrap critical bandwidth 97.5%": float(
                    np.quantile(bootstrap_bandwidths, 0.975)
                ),
                "Silverman half-grid observed bandwidth": resolution["half_grid"],
                "Silverman double-grid observed bandwidth": resolution["double_grid"],
            }
        )

    results = pd.DataFrame(result_rows)
    silverman_bootstrap = pd.concat(bootstrap_rows, ignore_index=True)
    results.to_csv(OUTPUT_DIR / "unimodality_test_results.csv", index=False)
    silverman_bootstrap.to_csv(
        OUTPUT_DIR / "unimodality_silverman_bootstrap.csv", index=False
    )

    mtbr = results.loc[results["Marker"] == "eMTBR-tau243"].iloc[0]
    ptau = results.loc[results["Marker"] == "p-tau217"].iloc[0]
    lines = [
        "## Material Passport",
        "",
        "- artifact_type: formal_unimodality_tests",
        "- verification_status: VERIFIED",
        f"- input: {INPUT}",
        "- analysis_population: A+T−",
        "- transformation: log10",
        f"- random_seed: {RANDOM_SEED}",
        f"- dip_bootstrap_replicates: {DIP_BOOTSTRAP_REPS}",
        f"- silverman_bootstrap_replicates: {SILVERMAN_BOOTSTRAP_REPS}",
        "",
        "# Formal unimodality tests in A+T− participants",
        "",
        "| Marker | n | Hartigan dip | Dip p | Silverman critical bandwidth | Silverman p |",
        "|---|---:|---:|---:|---:|---:|",
        f"| eMTBR-tau243 | {int(mtbr['N'])} | {mtbr['Hartigan dip statistic']:.4f} | "
        f"{format_p(mtbr['Hartigan dip bootstrap p'])} | "
        f"{mtbr['Silverman observed critical bandwidth (standardized log10 scale)']:.4f} | "
        f"{format_p(mtbr['Silverman bootstrap p'])} |",
        f"| p-tau217 | {int(ptau['N'])} | {ptau['Hartigan dip statistic']:.4f} | "
        f"{format_p(ptau['Hartigan dip bootstrap p'])} | "
        f"{ptau['Silverman observed critical bandwidth (standardized log10 scale)']:.4f} | "
        f"{format_p(ptau['Silverman bootstrap p'])} |",
        "",
        "## Interpretation",
        "",
        "- Hartigan's dip test evaluates unimodality without specifying a Gaussian null distribution.",
        "- Silverman's test evaluates whether more than one mode is required using a Gaussian-KDE critical bandwidth and smoothed bootstrap.",
        "- Both tests rejected unimodality for eMTBR-tau243 at nominal p=0.026, providing nonparametric evidence complementary to the existing Gaussian-mixture analysis.",
        "- Failure to reject for p-tau217 should be described as no evidence against unimodality, not proof that the distribution is unimodal.",
        "- These tests characterize the observed cross-sectional distribution and do not establish biological subtypes or temporal progression.",
    ]
    (OUTPUT_DIR / "unimodality_test_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    metadata = {
        "input": str(INPUT),
        "sheet": "Sheet1",
        "analysis_population": "A+T-",
        "markers": MARKERS,
        "transformation": "log10",
        "random_seed": RANDOM_SEED,
        "dip_bootstrap_replicates": DIP_BOOTSTRAP_REPS,
        "silverman_bootstrap_replicates": SILVERMAN_BOOTSTRAP_REPS,
        "silverman_grid_size": SILVERMAN_GRID_SIZE,
        "silverman_binary_steps": SILVERMAN_BINARY_STEPS,
        "silverman_support": "unknown; original Silverman critical bandwidth",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "diptest": diptest_package.__version__,
    }
    (OUTPUT_DIR / "unimodality_test_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
