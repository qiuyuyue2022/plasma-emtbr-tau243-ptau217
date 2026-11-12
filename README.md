# Plasma eMTBR-tau243 immunoassay versus p-tau217

Analysis code for the manuscript **"Plasma eMTBR-tau243 varies with both amyloid and tau status"** (PUMCH Dementia Cohort).

This repository contains the code that produced every statistic, table, and figure in the manuscript and its Supporting Information. It does not contain participant data (see [Data availability](#data-availability)).

---

## Contents

1. [The study in brief](#the-study-in-brief)
2. [Key definitions](#key-definitions)
3. [Analysis populations](#analysis-populations)
4. [Where each result comes from](#where-each-result-comes-from)
5. [Data availability and input format](#data-availability)
6. [Running the code](#running-the-code)
7. [Statistical methods at a glance](#statistical-methods-at-a-glance)
8. [Reproducibility notes](#reproducibility-notes)
9. [How to read the results](#how-to-read-the-results)

---

## The study in brief

**Question.** Plasma eMTBR-tau243 is a candidate marker of tau tangles. The existing evidence comes from mass spectrometry. We asked how a single-molecule **immunoassay** for eMTBR-tau243 relates to amyloid (A) and tau (T) status, and how it compares with plasma p-tau217.

**Design.** A cross-sectional diagnostic-accuracy study of 135 participants: 101 memory-clinic patients and 34 healthy controls. Seven plasma biomarkers were measured: eMTBR-tau243, p-tau217, p-tau231, BD-tau, GFAP, NfL, and Aβ1-42. Four Aβ42-normalized ratios were also derived, giving 11 measures in total.

**Main findings.**

- **T status.** eMTBR-tau243 discriminated tau PET-defined T status as well as p-tau217, with no significant paired difference (AUC 0.959 vs 0.947; n=67). Results were similar within amyloid-positive patients (AUC 0.941 vs 0.920; n=52). Adding either marker to the other did not clearly improve discrimination.
- **Amyloid and tau steps.** Across tau PET-defined A−T−, A+T−, and A+T+ groups, eMTBR-tau243 increased at both the amyloid step and the tau step. This differs from the tau-stage-specific profile reported for mass-spectrometry eMTBR-tau243; the immunoassay signal behaved more like a T1-type fluid tau marker.
- **Within A+T−.** eMTBR-tau243 values in tau PET-defined A+T− participants favored a two-component distribution, and the components tracked p-tau217 and p-tau231. This analysis is exploratory, and the two groups are descriptive.

> The source data label the analyte `MTBR-Tau243`. The manuscript calls it **eMTBR-tau243**. The immunoassay has not been compared with mass spectrometry in the same samples, so results from the two methods should not be assumed to be interchangeable.

---

## Key definitions

| Term | Definition used in the code |
| --- | --- |
| **A status** | Amyloid PET visual read. If amyloid PET was unavailable, CSF Aβ42/Aβ40 (LUMIPULSE; A+ if <0.07) was used. For 12 healthy controls with neither, a separate LUMIPULSE plasma p-tau217 reference assay was used; this is not the Inspirer DX p-tau217 evaluated in the study. Source column: `A_level`. |
| **T status** | **Tau PET only** ([18F]MK-6240, visual read; any elevated uptake = T+). Defined only for participants with tau PET within 30 days of blood collection (`taupet == 1`). |
| **Combined T reference** | Tau PET or, if tau PET was unavailable, CSF p-tau181 (T+ if >61.6 pg/mL). Used **only** in Supplementary Figures S4, S6, and S7 and Supplementary Table S9. These analyses are descriptive comparisons with the main tau PET-defined results. |
| **Patients / healthy controls** | Healthy controls are identified by the diagnosis column (`DIAG == "健康对照"`, Chinese for "healthy control"). All other participants are patients. |
| **Amyloid step / tau step** | A−T− → A+T− and A+T− → A+T+, respectively. |

---

## Analysis populations

| Analysis | Population | n | Role |
| --- | --- | --- | --- |
| A status | Patients | 101 (A− 34, A+ 67) | Primary A-status analysis |
| A status | All participants | 135 | Secondary |
| T status | Patients with tau PET | 67 (T− 33, T+ 34) | **Primary analysis of the study** |
| T status within A+ | A+ patients with tau PET | 52 (T− 18, T+ 34) | Primary A+-conditioned analysis |
| Two-step comparison | Tau PET-defined A−T− / A+T− / A+T+ | 15 / 18 / 34 | Secondary |
| A+T− heterogeneity | Tau PET-defined A+T− | 18 | Exploratory |
| Combined-reference A/T groups | A−T− / A+T− / A+T+ (T by tau PET or CSF) | 46 / 33 / 34 | Supplementary, descriptive |
| Combined-reference A+T− heterogeneity | A+T− (T by tau PET or CSF) | 33 | Supplementary, descriptive |

---

## Where each result comes from

Script names start with the display item they produce: `Fig#`, `FigS#`, `Table#`, or `TableS#`. Scripts whose results appear only in the text start with `Text_`. Scripts for the combined tau PET/CSF analyses in the Supporting Information start with `SI_PETCSF_`.

### Main text

| Manuscript item | Script(s) |
| --- | --- |
| Figure 1 — cohort and analysis-population flow | `Fig1_study_flow.py` |
| Figure 2 — biomarker distributions by A status (patients) | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure 3 — ROC curves, four analysis populations | `Fig3_FigS2_TableS4_ROC.py` |
| Figure 4 — biomarker distributions by tau PET-defined T status | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure 5 — two-step changes across tau PET-defined A/T groups | `Fig5_AT_two_step_comparison.py` |
| Figure 6 — profiles of the eMTBR-tau243 groups within tau PET-defined A+T− | `TableS8_Aplus_Tminus_subgroup_comparison.py` → `Fig6_A+T-_subgroup_boxplots.py` |
| Table 1 — participant characteristics | `Table1_TableS3_participant_characteristics.py` → `build_formatted_manuscript_tables.py` |
| Results reported in the text only (storage time, blood–CSF interval, CSF- vs tau PET-defined A+T−) | `Text_storage_time_and_intervals.py`, `Text_Aplus_Tminus_reference_comparison.py` |

### Supporting Information

| Item | Script(s) |
| --- | --- |
| Table S1 — assay parameters | Compiled from the manufacturer's kit specifications (no code) |
| Table S2 — best measure vs the other ten (paired DeLong and bootstrap) | `TableS2_diagnostic_validation.py` |
| Table S3 — concentrations by A/T group | `Table1_TableS3_participant_characteristics.py`, `TableS3_AT_group_statistics.py` |
| Table S4 — full diagnostic performance, 11 measures × 4 populations | `Fig3_FigS2_TableS4_ROC.py` |
| Table S5 — Aβ1-42 below-LLOQ sensitivity analysis | `TableS5_AB42_LLOQ_sensitivity.py` (see [Sensitivity analysis](#sensitivity-analysis-ab1-42-below-the-lloq)) |
| Table S6 — incremental value of eMTBR-tau243 and p-tau217 | `TableS6_incremental_value.py` |
| Table S7 — mixture-model and unimodality diagnostics (n=18) | `TableS7_FigS6_GMM_full_sample.py`, `TableS7_unimodality_tests.py` |
| Table S8 — characteristics of the A+T− groups (n=18) | `TableS8_Aplus_Tminus_subgroup_comparison.py` |
| Table S9 — A+T− heterogeneity, combined T reference (n=33) | `SI_PETCSF_*.py` → `TableS9_PETCSF_Aplus_Tminus.py` |
| Figure S1 — distributions by A status, all participants | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure S2 — AUC forest plot, four populations | `Fig3_FigS2_TableS4_ROC.py` |
| Figure S3 — distributions by T status within A+ | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure S4 — A/T groups with the combined T reference | `FigS4_AT_3group_PET_CSF.py` |
| Figure S5 — mixture model within tau PET-defined A+T− (n=18) | `FigS5_Aplus_Tminus_heterogeneity.py` |
| Figure S6 — mixture models within A+T−, combined T reference (n=33) | `SI_PETCSF_GMM_full_sample.py`, `SI_PETCSF_GMM_by_reference.py` → `FigS6_PETCSF_Aplus_Tminus_bimodality.py` |
| Figure S7 — clinical, plasma, and CSF profiles, combined T reference (n=33) | `SI_PETCSF_subgroup_comparison.py` → `FigS7_PETCSF_Aplus_Tminus_subgroup_boxplots.py` |

### Supporting tools (not needed to reproduce the statistics)

| Script | Purpose |
| --- | --- |
| `QC_audit_manuscript_numbers.py` | Checks that every key number in the manuscript text matches the analysis outputs |
| `assemble_manuscript.py`, `assemble_supplementary.py`, `renumber_citations.py` | Build the Word manuscript and Supporting Information from Markdown. They need the manuscript files, which are not in this repository. |
| `Table1_build_workbook_excel.mjs` | Optional Excel export of Table 1. It depends on an internal Node package and is not required. |

---

<a id="data-availability"></a>
## Data availability and input format

Individual-level clinical, imaging, and biomarker data are **not included** because of patient-privacy requirements and the terms of the ethical approval from Peking Union Medical College Hospital. Anonymized data are available from the corresponding author on reasonable request, as described in the manuscript. Without the study data, the code cannot reproduce the reported numbers.

The scripts read one Excel workbook, `biomarker_result.xlsx`, worksheet `Sheet1`, with one row per participant. Required columns:

| Group | Columns | Coding |
| --- | --- | --- |
| Identifier | `SID` | Study ID |
| Reference standards | `A`, `T`, `A_level`, `T_level`, `taupet`, `method` | `A`: `A+`/`A-`; `T`: `T+`/`T-`/`T0` (no T reference); `A_level`: `PET`/`CSF`/`Plasma`; `T_level`: `PET`/`CSF`; `taupet`: 1 if tau PET was within 30 days of blood collection |
| Diagnosis | `DIAG`, `N` | `DIAG == "健康对照"` marks healthy controls |
| Demographics and cognition | `AGE`, `SEX`, `apoe`, `MMSE_Total`, `MOCA_P` | `SEX`: `男` (male) / `女` (female); `apoe`: genotype such as `ε3/ε4` |
| Plasma biomarkers (pg/mL) | `MTBR-Tau243`, `p-Tau217(≤0.758)`, `p-Tau231(≤5.793)`, `BD Tau（≤52.301）`, `GFAP(≤131.880)`, `NFL(≤9.624 )`, `Aβ 1-42` | Column names are kept exactly as exported by the laboratory. The values in parentheses are the laboratory's reference limits; they are not used as cutoffs. |
| Sample timing | `sample_date` | Blood collection date (used for storage time) |
| CSF (Lumipulse) | `CSF_Lumipulse_Aβ42/Aβ40`, `CSF_Lumipulse_p-tau`, `CSF_Lumipulse_t-tau`, `CSF_Lumipulse_p-tau/Aβ42`, `CSF_Lumipulse_t-tau/Aβ42`, `CSF_Lumipulse_采样时间` (collection date), `CSF_Lumipulse_绝对时间差_天` (absolute blood–CSF interval, days) | Used for the combined T reference and the CSF profiles |

Many scripts include **assertions that lock the cohort counts** to the published values (for example, 101 patients and 67 with tau PET). On any other dataset these checks stop the script on purpose. Edit or remove them to use the code on new data.

---

## Running the code

### 1. Environment

The analyses were run with **Python 3.9.24**.

```bash
python -m pip install -r requirements.txt
```

`diptest` (Hartigan's dip test) is used by the unimodality scripts. If a system-wide install is not possible, install it locally with `python -m pip install --target _vendor diptest==0.10.0`; the scripts look for it there.

### 2. Folder layout

```text
project/
├── biomarker_result.xlsx          ← study data (not provided)
├── code/                          ← this repository (any folder name)
└── outputs/                       ← created automatically
```

By default, the scripts read `../biomarker_result.xlsx` and write to `../outputs/`. Most scripts accept the environment variables `BIOMARKER_INPUT` (path to the workbook) and `BIOMARKER_OUTPUT_DIR` (output folder).

### 3. Run order

Run from inside `code/`. Later steps read outputs from earlier ones.

```bash
# Diagnostic performance (Figure 3, Tables S2 and S4, Figure S2)
python Fig3_FigS2_TableS4_ROC.py
python TableS2_diagnostic_validation.py
python Fig3_FigS2_TableS4_ROC.py          # rerun so the Figure 3 legend uses the Table S2 comparisons
python TableS6_incremental_value.py

# Two-step A/T comparison and group tables (Figure 5, Table S3, Table 1)
python Fig5_AT_two_step_comparison.py
python TableS3_AT_group_statistics.py
python Table1_TableS3_participant_characteristics.py

# Exploratory A+T− analysis, tau PET-defined (n=18): Tables S7 and S8
python TableS7_FigS6_GMM_full_sample.py
python TableS7_unimodality_tests.py
python TableS8_Aplus_Tminus_subgroup_comparison.py

# Combined tau PET/CSF reference (n=33): Table S9
python SI_PETCSF_GMM_full_sample.py
python SI_PETCSF_unimodality_tests.py
python SI_PETCSF_GMM_by_reference.py
python SI_PETCSF_subgroup_comparison.py
python TableS9_PETCSF_Aplus_Tminus.py

# Figures
python Fig1_study_flow.py
python Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py
python Fig6_A+T-_subgroup_boxplots.py
python FigS4_AT_3group_PET_CSF.py
python FigS5_Aplus_Tminus_heterogeneity.py
python FigS6_PETCSF_Aplus_Tminus_bimodality.py
python FigS7_PETCSF_Aplus_Tminus_subgroup_boxplots.py

# Text-only results and formatted tables
python Text_storage_time_and_intervals.py
python Text_Aplus_Tminus_reference_comparison.py
python build_formatted_manuscript_tables.py
```

Figures are saved as SVG and PDF (editable text), PNG, and TIFF. Most figures also have a source-data CSV and a metadata JSON.

<a id="sensitivity-analysis-ab1-42-below-the-lloq"></a>
### 4. Sensitivity analysis: Aβ1-42 below the LLOQ

Supplementary Table S5 tests whether plasma Aβ1-42 values below the lower limit of quantification (1.546 pg/mL) affect the results.

1. Make a copy of the workbook in which every `Aβ 1-42` value below 1.546 is set to 1.546.
2. Rerun the main scripts on that copy, writing to a separate folder:

```bash
export BIOMARKER_INPUT=../outputs/sensitivity_AB42_LLOQ_censoring/biomarker_result_AB42_censored_at_LLOQ.xlsx
export BIOMARKER_OUTPUT_DIR=../outputs/sensitivity_AB42_LLOQ_censoring
python Fig3_FigS2_TableS4_ROC.py
python TableS2_diagnostic_validation.py
python Fig3_FigS2_TableS4_ROC.py
python TableS3_AT_group_statistics.py
BIOMARKER_OUTPUT_DIR=../outputs/sensitivity_AB42_LLOQ_censoring/Aplus_Tminus_tauPET_only \
  python TableS8_Aplus_Tminus_subgroup_comparison.py
unset BIOMARKER_INPUT BIOMARKER_OUTPUT_DIR
python TableS5_AB42_LLOQ_sensitivity.py   # compares the main and censored outputs
```

---

## Statistical methods at a glance

| Analysis | Method |
| --- | --- |
| Discrimination | ROC AUC with 95% CI from participant-level bootstrap. The direction of each marker is chosen in training folds for cross-validated AUC. |
| Comparing AUCs | Paired DeLong test and paired bootstrap ΔAUC. Benjamini–Hochberg (BH) FDR across the ten comparisons with the best measure in each population. The "leading group" is the best measure plus all measures not significantly lower than it. |
| Incremental value | Logistic regression on log10 concentrations; likelihood-ratio tests; bootstrap ΔAUC; repeated stratified 5-fold cross-validation |
| Two-step comparison | Mann–Whitney U tests with rank-biserial r; BH-FDR across the seven absolute markers within each step |
| A+T− heterogeneity | One- vs two-component Gaussian mixture models on log10 eMTBR-tau243, compared by BIC. Parametric-bootstrap likelihood-ratio test (1,000 replicates). Leave-one-out and extreme-value checks. Hartigan's dip test and Silverman's critical-bandwidth test. |
| Group profiles | Mann–Whitney U tests, or Fisher's exact test for categorical variables. Separate BH-FDR families for clinical, plasma, and CSF variables; exact permutation tests for the small CSF subset. |
| Distribution figures | Gaussian kernel densities fitted on the log10 scale; descriptive only |

All tests are two-sided. FDR q<0.05 is considered significant.

---

## Reproducibility notes

- **Fixed random seeds.** All bootstrap, cross-validation, and mixture-model steps use fixed seeds, set at the top of each script (most can be overridden with `BIOMARKER_RANDOM_SEED`). Rerunning on the same data gives identical outputs.
- **Read-only input.** No script modifies the input workbook.
- **Numeric audit.** `QC_audit_manuscript_numbers.py` compares the key numbers in the manuscript with the outputs. All 111 checks pass for the submitted version.
- **Retired code.** Analyses removed during revision are not part of this repository.

---

## How to read the results

- **Single-center, cross-sectional data.** Cutoffs were derived in this sample and have not been externally validated.
- **Research-use-only assay.** The eMTBR-tau243 kit is for research use only and should not be used for clinical T staging.
- **Exploratory A+T− analysis.** The two eMTBR-tau243 groups within A+T− are model-based, descriptive groups from small samples (n=18; n=33 with the combined reference). They are not validated clinical categories or temporal stages.
- **Open questions.** Same-sample comparison with mass spectrometry, quantitative PET, longitudinal follow-up, and independent replication are needed to define the biological and clinical meaning of the immunoassay signal.

For questions about the analyses, please contact the corresponding author listed in the manuscript.
