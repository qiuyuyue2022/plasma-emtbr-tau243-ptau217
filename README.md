# Analysis code for plasma eMTBR-tau243 and p-tau217 in a memory-clinic cohort

This repository contains the analysis code for the study **“A plasma eMTBR-tau243 immunoassay versus p-tau217 for tau status in a memory clinic.”**

The study compares plasma eMTBR-tau243 and p-tau217 for discrimination of amyloid (A) and tau (T) status. It also explores variation in eMTBR-tau243 measurements among amyloid-positive, tau-negative (A+T−) participants. The scripts generate the diagnostic analyses, figures, and tables reported in the manuscript.

## Repository contents

| Code | Purpose |
| --- | --- |
| `Fig3_ROC.py` | Calculates biomarker performance and generates Figure 3 and Supplementary Figure S2. |
| `01_diagnostic_validation.py` | Performs cross-validation, cutoff stability analyses, and paired AUC comparisons. |
| `02a_bimodality_GMM_full_sample.py` | Fits Gaussian mixture models in the A+T− group and assesses model stability. |
| `02b_bimodality_GMM_by_reference.py` | Repeats mixture-model analyses within tau-reference subgroups. |
| `02c_unimodality_tests.py` | Performs additional tests of unimodality in the A+T− group. |
| `03_Aplus_Tminus_subgroup_comparison.py` | Compares clinical, plasma, and available CSF measurements between the A+T− subgroups. |
| `Table1.py` and `build_formatted_manuscript_tables.py` | Generate the participant-characteristics table and formatted supplementary tables. |
| `Fig1_study_flow.py` | Generates the study-flow diagram (Figure 1). |
| `Fig2Fig4_biomarker_distributions.py` | Generates the biomarker distributions in Figures 2 and 4 and Supplementary Figures S1 and S3–S5. |
| `Fig5_AT_3group_distributions.py` | Generates Figure 5, showing biomarker distributions across A/T groups. |
| `Fig6_A+T-_subgroup_boxplots.py` | Generates Figure 6, comparing the A+T− subgroups. |
| `FigS6_plot_A+T-_integrated_figure.py` | Generates Supplementary Figure S6. |

`requirements.txt` lists the Python dependencies. The `_obsolete/` directory contains historical scripts and is not part of the current analysis workflow.

## Data availability

Individual-level clinical, imaging, and biomarker data are not included in this repository because of patient privacy requirements and the terms of the study’s ethical approval. Anonymized data supporting the findings may be requested from the corresponding author, subject to applicable approvals.

The code alone cannot regenerate the reported numerical results without the study data. The manuscript describes the cohort, reference standards, assays, and statistical methods.

## Running the analyses

The scripts expect an Excel workbook named `biomarker_result.xlsx` in the parent directory of this repository. They read the worksheet named `Sheet1` and write results to a sibling `outputs/` directory:

```text
project/
├── biomarker_result.xlsx
├── code/                 # this repository
└── outputs/              # generated files
```

Install the Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

The main analyses can then be run from the repository directory in this order:

```bash
python Fig3_ROC.py
python 01_diagnostic_validation.py
python 02a_bimodality_GMM_full_sample.py
python 02b_bimodality_GMM_by_reference.py
python 02c_unimodality_tests.py
python 03_Aplus_Tminus_subgroup_comparison.py
python Table1.py
python build_formatted_manuscript_tables.py
```

The figure scripts can be run after the analyses that supply their input files:

```bash
python Fig1_study_flow.py
python Fig2Fig4_biomarker_distributions.py
python Fig5_AT_3group_distributions.py
python Fig6_A+T-_subgroup_boxplots.py
python FigS6_plot_A+T-_integrated_figure.py
```

Several scripts also accept `BIOMARKER_INPUT` and `BIOMARKER_OUTPUT_DIR` environment variables to specify alternative input and output paths. `Table1.py` uses the parent-directory layout shown above.

## Notes on interpretation

The reported cutoffs were estimated within the study cohort for research analyses. They have not been established as clinical decision thresholds. The exploratory A+T− subgroup analyses describe variation within this cohort and require independent validation.
