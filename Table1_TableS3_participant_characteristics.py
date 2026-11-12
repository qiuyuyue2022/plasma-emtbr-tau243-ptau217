#!/usr/bin/env python3
"""Prepare manuscript Table 1 draft and Supplementary Table S3 (renumbered from S2 on 2026-11-05).

Outputs: Table1_participant_characteristics_draft.{csv,md}, Table1_metadata.json,
Table1_workbook_data.json (consumed by Table1_build_workbook_excel.mjs), and
TableS3_biomarkers_by_AT_group.{csv,md} (absolute biomarker concentrations and
full A/T distribution statistics; Supplementary Table S3).

The raw workbook is never modified. Table1_participant_characteristics_draft.md
is further revised (A-T0 label fix etc.) by
build_formatted_manuscript_tables.py into the final Table1_participant_characteristics.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "biomarker_result.xlsx"
OUTPUT_DIR = ROOT / "outputs"


# v3.60 (2026-11-11): T is defined by tau PET only, in patients only. Columns are the
# tau PET-defined A/T groups of patients, patients without tau PET, and healthy controls.
GROUP_ORDER = ["Overall", "A−T−", "A+T−", "A+T+", "Patients without tau PET", "Healthy controls"]


def assign_group(row: pd.Series) -> str:
    if row["DIAG"] == "健康对照":
        return "Healthy controls"
    if row["taupet"] != 1 or row["T"] not in ("T+", "T-"):
        return "Patients without tau PET"
    if row["A"] == "A-" and row["T"] == "T-":
        return "A−T−"
    if row["A"] == "A+" and row["T"] == "T-":
        return "A+T−"
    if row["A"] == "A+" and row["T"] == "T+":
        return "A+T+"
    raise ValueError(f"Unexpected A/T combination: A={row['A']!r}, T={row['T']!r}")


def subset(df: pd.DataFrame, group: str) -> pd.DataFrame:
    return df if group == "Overall" else df.loc[df["AT_group"] == group]


def fmt_mean_sd(values: pd.Series) -> str:
    x = pd.to_numeric(values, errors="coerce").dropna()
    return "—" if x.empty else f"{x.mean():.1f} ({x.std(ddof=1):.1f})"


def fmt_median_iqr(values: pd.Series, decimals: int = 1) -> str:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if x.empty:
        return "—"
    q1, med, q3 = x.quantile([0.25, 0.5, 0.75])
    return f"{med:.{decimals}f} [{q1:.{decimals}f}, {q3:.{decimals}f}]"


def fmt_median_iqr_n(values: pd.Series, decimals: int = 1) -> str:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if x.empty:
        return "—"
    return f"{fmt_median_iqr(x, decimals)}; n={len(x)}"


def fmt_n_pct(mask: pd.Series, denominator: int | None = None) -> str:
    valid = mask.dropna()
    den = len(valid) if denominator is None else denominator
    num = int(valid.astype(bool).sum())
    return "—" if den == 0 else f"{num} ({100 * num / den:.1f}%)"


def fmt_n_over_n_pct(mask: pd.Series) -> str:
    valid = mask.dropna()
    den = len(valid)
    num = int(valid.astype(bool).sum())
    return "—" if den == 0 else f"{num}/{den} ({100 * num / den:.1f}%)"


def apoe4_carrier(value: object) -> object:
    if pd.isna(value):
        return np.nan
    return "ε4" in str(value)


def diagnosis_group(value: object) -> str:
    value = str(value)
    # Clinical syndrome only (v3.33): the source labels AD / ADlike / ADmimic all denote a
    # clinically probable AD syndrome and differ only by A/T biomarker status (A+T+ / A+T- /
    # A-T-), so they are pooled here to keep the clinical dimension independent of A/T.
    probable_ad_syndrome = {
        "AD", "ADlike", "ADmimic",
        "ADlike+CAA", "LBD+ADlike", "AD+LBD", "AD+CAA", "ADlike+PSP",
    }
    if value == "健康对照":
        return "Healthy control"
    if value in probable_ad_syndrome:
        return "Probable AD syndrome"
    return "Other clinical syndrome"


def build_characteristics(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    group_ns = {g: len(subset(df, g)) for g in GROUP_ORDER}
    columns = ["Characteristic"] + [f"{g} (n={group_ns[g]})" for g in GROUP_ORDER]
    rows: list[list[str]] = []

    def add(label: str, formatter) -> None:
        rows.append([label] + [formatter(subset(df, g)) for g in GROUP_ORDER])

    add("Age, years", lambda d: fmt_mean_sd(d["AGE"]))
    add("Female sex", lambda d: fmt_n_pct(d["SEX"].eq("女"), len(d)))
    add("APOE ε4 carrier", lambda d: fmt_n_over_n_pct(d["APOE4_carrier"]))
    add("MMSE score", lambda d: fmt_median_iqr_n(d["MMSE_Total"], 1))
    add("MoCA-P score", lambda d: fmt_median_iqr_n(d["MOCA_P"], 1))
    add("Amyloid positive", lambda d: fmt_n_pct(d["A"].eq("A+"), len(d)))
    rows.append(["A reference"] + [""] * len(GROUP_ORDER))
    add("  PET", lambda d: fmt_n_pct(d["A_level"].eq("PET"), len(d)))
    add("  CSF", lambda d: fmt_n_pct(d["A_level"].eq("CSF"), len(d)))
    add("  Plasma", lambda d: fmt_n_pct(d["A_level"].eq("Plasma"), len(d)))

    out = pd.DataFrame(rows, columns=columns)
    metadata = {
        "input": str(INPUT),
        "n_total": len(df),
        "group_n": group_ns,
        "apoe_missing": int(df["apoe"].isna().sum()),
        "mmse_missing": int(df["MMSE_Total"].isna().sum()),
        "moca_p_missing": int(df["MOCA_P"].isna().sum()),
        "notes": [
            "Values are mean (SD), median [IQR], or n (%) unless otherwise indicated.",
            "APOE ε4 carriage is shown as n/N (%) because genotype was unavailable for some participants.",
            "A reference indicates the method by which amyloid status was established: PET or CSF for patients; PET, CSF, or (when neither was available) a plasma p-tau217 reference assay on the LUMIPULSE platform for healthy controls (Methods 2.2). This plasma reference assay is analytically distinct from the plasma p-tau217 measured on the Inspirer DX platform evaluated as a candidate biomarker in this study.",
            "A−T−, A+T−, and A+T+ are patients with tau PET; T status was defined by tau PET only. These 67 patients form the T-status population and the three groups of the two-step comparison.",
            "Patients without tau PET and healthy controls contributed to A-status and descriptive analyses only.",
        ],
    }
    return out, metadata


def build_biomarkers(df: pd.DataFrame) -> pd.DataFrame:
    markers = [
        ("GFAP", "GFAP(≤131.880)"),
        ("NfL", "NFL(≤9.624 )"),
        ("Aβ1-42", "Aβ 1-42"),
        ("p-tau217", "p-Tau217(≤0.758)"),
        ("p-tau231", "p-Tau231(≤5.793)"),
        ("BD-tau", "BD Tau（≤52.301）"),
        ("eMTBR-tau243", "MTBR-Tau243"),
    ]
    group_ns = {g: len(subset(df, g)) for g in GROUP_ORDER}
    columns = ["Biomarker"] + [f"{g} (n={group_ns[g]})" for g in GROUP_ORDER]
    rows = []
    for label, source_col in markers:
        rows.append(
            [label]
            + [fmt_median_iqr(subset(df, g)[source_col], 3) for g in GROUP_ORDER]
        )
    return pd.DataFrame(rows, columns=columns)


def markdown_table(df: pd.DataFrame) -> str:
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in df.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(INPUT, sheet_name="Sheet1")
    df["AT_group"] = df.apply(assign_group, axis=1)
    df["APOE4_carrier"] = df["apoe"].map(apoe4_carrier)
    df["diagnosis_group"] = df["DIAG"].map(diagnosis_group)
    df["T_reference"] = np.where(
        df["T"].eq("T0"),
        "Not available",
        np.where(df["taupet"].eq(1), "Tau PET", "CSF"),
    )

    table1, metadata = build_characteristics(df)
    biomarkers = build_biomarkers(df)

    table1.to_csv(OUTPUT_DIR / "Table1_participant_characteristics_draft.csv", index=False)
    biomarkers.to_csv(OUTPUT_DIR / "TableS3_biomarkers_by_AT_group.csv", index=False)
    (OUTPUT_DIR / "Table1_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "Table1_workbook_data.json").write_text(
        json.dumps(
            {
                "table1": [list(table1.columns)] + table1.astype(object).values.tolist(),
                "biomarkers": [list(biomarkers.columns)]
                + biomarkers.astype(object).values.tolist(),
                "metadata": metadata,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    md = [
        "# Table 1. Participant characteristics by tau PET-defined A/T group",
        "",
        markdown_table(table1),
        "",
        "## Notes",
        "",
    ]
    md.extend(f"- {note}" for note in metadata["notes"])
    md.extend(
        [
            "",
            "# Supplementary table. Absolute biomarker concentrations by A/T group",
            "",
            markdown_table(biomarkers),
            "",
            "Values are median [IQR]. Plasma biomarker concentrations are in pg/mL.",
        ]
    )
    (OUTPUT_DIR / "Table1_participant_characteristics_draft.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8"
    )

    # Amyloid-step and tau-step comparisons for the seven absolute biomarkers, read from
    # AT_two_step_comparison.csv (Fig5_AT_two_step_comparison.py; BH-FDR across the seven
    # markers within each step). No statistic is recomputed here.
    s3 = biomarkers.copy()
    stats_path = OUTPUT_DIR / "AT_two_step_comparison.csv"
    st = pd.read_csv(stats_path).set_index("Marker")
    name_map = {"NfL": "NfL"}
    def qfmt(q):
        return "<0.001" if q < 0.001 else f"{q:.3f}"
    for step, rcol, qcol in [("Amyloid step", "Amyloid step r (A+T− vs A−T−)", "Amyloid step q"),
                             ("Tau step", "Tau step r (A+T+ vs A+T−)", "Tau step q")]:
        s3[f"{step}, r"] = s3["Biomarker"].map(lambda b: f"{st.loc[name_map.get(b, b), rcol]:.2f}")
        s3[f"{step}, q"] = s3["Biomarker"].map(lambda b: qfmt(st.loc[name_map.get(b, b), qcol]))
    gn = metadata["group_n"]
    stats_note = (" The amyloid step compares A+T− with A−T−, and the tau step compares A+T+ with A+T−, "
                  "among patients with tau PET-defined T status (Mann–Whitney U tests; positive rank-biserial r "
                  "indicates higher concentrations in the later group; Benjamini–Hochberg q across the seven "
                  "markers within each step; same statistics as Figure 5).")
    tables2_md = [
        "# Supplementary Table S3. Plasma biomarker concentrations by tau PET-defined A/T group",
        "",
        markdown_table(s3),
        "",
        (f"Values are median [interquartile range] of absolute plasma concentrations (Overall, n={gn['Overall']}; "
         f"tau PET-defined A−T−, n={gn['A−T−']}; A+T−, n={gn['A+T−']}; A+T+, n={gn['A+T+']}; "
         f"patients without tau PET, n={gn['Patients without tau PET']}; healthy controls, n={gn['Healthy controls']}). "
         "All concentrations are in pg/mL." + stats_note),
    ]
    (OUTPUT_DIR / "TableS3_biomarkers_by_AT_group.md").write_text(
        "\n".join(tables2_md) + "\n", encoding="utf-8"
    )

    print(table1.to_string(index=False))
    print("\nGroup counts:", metadata["group_n"])


if __name__ == "__main__":
    main()
