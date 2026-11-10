# Plasma eMTBR-tau243 immunoassay versus p-tau217

Analysis code accompanying the manuscript **“Plasma eMTBR-tau243 measured by immunoassay discriminates tau PET status but rises with amyloid before tau PET positivity.”**

## Study overview

This cross-sectional diagnostic-accuracy study measured eMTBR-tau243 (research-use-only assay) and six other plasma biomarkers by single-molecule immunoassay in 135 participants of the PUMCH Dementia Cohort: 101 memory-clinic patients and 34 healthy controls. The primary comparison was the ability of eMTBR-tau243 and p-tau217 to discriminate tau PET-defined T status, overall (n=68) and within amyloid-positive participants (n=52). Secondary analyses examined amyloid (A) status, marker changes from A−T− to A+T− (amyloid step) and from A+T− to A+T+ (tau step) in tau PET-defined groups (n=16, 18 and 34), and the distribution of eMTBR-tau243 within the 18 tau PET-defined A+T− participants.

eMTBR-tau243 discriminated tau PET-defined T status similarly to p-tau217, with no significant paired difference, and adding either marker to the other did not clearly improve discrimination. Unlike reports for mass-spectrometry eMTBR-tau243, the immunoassay signal changed at both the amyloid step and the tau step, behaving more like a T1-type fluid tau marker. Within tau PET-defined A+T−, eMTBR-tau243 values favored a two-component distribution that tracked p-tau217 and p-tau231; this analysis is exploratory, and the components are descriptive model-based groups, not validated clinical categories.

The source-data label `MTBR-Tau243` refers to the analyte described as eMTBR-tau243 in the manuscript. The immunoassay has not been compared with mass-spectrometry measurement in the same samples; results from the two methods should not be assumed interchangeable.

## Code and manuscript figures

Each script name begins with the manuscript display item(s) it produces or supplies (`Fig#`, `FigS#`, `Table#`, `TableS#`). Scripts whose results are reported only in the text begin with `Text_`; `QC_` marks the numeric audit.

| Manuscript item | Script |
| --- | --- |
| Figure 1 — study cohort and analysis populations | `Fig1_study_flow.py` |
| Figure 2 — biomarker distributions by A status in patients | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure 3 — ROC curves across six analysis populations | `Fig3_FigS2_TableS4_ROC.py` |
| Figure 4 — biomarker distributions by tau PET-defined T status | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Figure 5 — two-step changes across tau PET-defined A−T−, A+T−, and A+T+ | `Fig5_AT_two_step_comparison.py` |
| Figure 6 — clinical and molecular profiles of the exploratory A+T− groups | `Fig6_A+T-_subgroup_boxplots.py` |
| Table 1 — participant characteristics | `Table1_TableS3_participant_characteristics.py`, `build_formatted_manuscript_tables.py`, `Table1_build_workbook_excel.mjs` |
| Supplementary Figures S1, S3–S5 — biomarker distributions in secondary populations | `Fig2_Fig4_FigS1_S3-S5_biomarker_distributions.py` |
| Supplementary Figure S2 — AUC forest plot | `Fig3_FigS2_TableS4_ROC.py` |
| Supplementary Figure S6 — A+T− mixture-model analyses | `FigS6_Aplus_Tminus_heterogeneity.py` |
| Supplementary Table S2 — best measure versus the other ten | `TableS2_diagnostic_validation.py` |
| Supplementary Table S3 — concentrations by A/T group | `Table1_TableS3_participant_characteristics.py`, `TableS3_AT_group_statistics.py` |
| Supplementary Table S4 — full diagnostic performance | `Fig3_FigS2_TableS4_ROC.py` |
| Supplementary Table S5 — Aβ1-42 LLOQ sensitivity analysis | `TableS5_AB42_LLOQ_sensitivity.py` |
| Supplementary Table S6 — incremental value of eMTBR-tau243 and p-tau217 | `TableS6_incremental_value.py` |
| Supplementary Table S7 — mixture-model and unimodality diagnostics | `TableS7_FigS6_GMM_full_sample.py`, `TableS7_FigS6_GMM_by_reference.py`, `TableS7_unimodality_tests.py` |
| Supplementary Table S8 — characteristics of the A+T− groups | `TableS8_Aplus_Tminus_subgroup_comparison.py` |
| Results reported in the text only | `Text_storage_time_and_intervals.py`, `Text_Aplus_Tminus_reference_comparison.py`, `Text_age_adjusted_amyloid_step.py` |

Supplementary Table S1 (assay parameters) was compiled from the manufacturer's kit specifications. `build_formatted_manuscript_tables.py` formats Table 1 and Supplementary Tables S2, S4 and S6–S8 without recomputing statistics; `assemble_supplementary.py` assembles the Supporting Information. `TableS3_AT_group_statistics.py` also draws an unnumbered reference figure (`FigS_AT_3group_distributions`).

## Data availability

Individual-level clinical, imaging, and biomarker data are not included in this repository because of patient privacy requirements and the terms of the study’s ethical approval. Anonymized data supporting the findings are available from the corresponding author upon reasonable request, as described in the manuscript.

The code cannot reproduce the reported numerical results without the study data.

## Running the code

The analyses were performed using Python 3.9.24. Python dependencies are listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

The scripts expect an Excel workbook named `biomarker_result.xlsx` with a worksheet named `Sheet1`. Place it alongside the cloned repository, with generated files written to a sibling `outputs/` directory:

```text
project/
├── biomarker_result.xlsx
├── plasma-emtbr-tau243-ptau217/
└── outputs/
```

Some scripts contain defaults specific to the original analysis computer. For scripts that support them, set `BIOMARKER_INPUT` and `BIOMARKER_OUTPUT_DIR` to the appropriate local paths before running. `Table1_TableS3_participant_characteristics.py` reads `biomarker_result.xlsx` from the repository’s parent directory.

Run the principal analysis scripts before the figure and table scripts that use their outputs:

```bash
python Fig3_FigS2_TableS4_ROC.py
python TableS2_diagnostic_validation.py
python TableS7_FigS6_GMM_full_sample.py
python TableS7_FigS6_GMM_by_reference.py
python TableS7_unimodality_tests.py
python TableS8_Aplus_Tminus_subgroup_comparison.py
python TableS6_incremental_value.py
python Fig5_AT_two_step_comparison.py
python TableS3_AT_group_statistics.py
python Table1_TableS3_participant_characteristics.py
python build_formatted_manuscript_tables.py
```

The figure scripts can then be run individually using the manuscript-to-script table above.

## Interpretation

The analyses assess discrimination within this single-center cohort. Sample-derived cutoffs were not externally validated, and the research-use-only assay should not be used for clinical T staging. Same-sample comparison with mass spectrometry, quantitative PET, longitudinal follow-up and independent validation are needed to define the biological and clinical meaning of the immunoassay signal.
