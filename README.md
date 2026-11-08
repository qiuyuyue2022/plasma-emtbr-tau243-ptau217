# Plasma eMTBR-tau243 immunoassay versus p-tau217

Analysis code accompanying the manuscript **“A plasma eMTBR-tau243 immunoassay versus p-tau217 for tau status in a memory clinic.”**

## Study overview

This cross-sectional study examined plasma biomarkers in 135 participants: 101 memory-clinic patients and 34 healthy controls. Its primary comparison was the ability of immunoassay-measured eMTBR-tau243 and p-tau217 to discriminate tau PET-defined T status, including within amyloid-positive participants. The analyses also examined amyloid (A) status and explored the distribution of eMTBR-tau243 among amyloid-positive, tau-negative (A+T−) participants.

The two markers had similar observed AUCs for the main comparisons, with no statistically significant paired differences detected. Within the 33 A+T− participants, eMTBR-tau243 measurements favored a two-component distribution. This subgroup finding is exploratory. The components are descriptive model-based groups, not validated clinical categories or diagnostic thresholds.

The source-data label `MTBR-Tau243` refers to the analyte described as eMTBR-tau243 in the manuscript. The immunoassay has not been directly compared with mass-spectrometry measurement in the same samples; results from the two methods should not be assumed interchangeable.

## Code and manuscript figures

| Manuscript item | Script |
| --- | --- |
| Figure 1 — study cohort and analysis populations | `Fig1_study_flow.py` |
| Figure 2 — biomarker distributions by A status in patients | `Fig2Fig4_biomarker_distributions.py` |
| Figure 3 — ROC curves across six analysis populations | `Fig3_ROC.py` |
| Figure 4 — biomarker distributions by tau PET-defined T status | `Fig2Fig4_biomarker_distributions.py` |
| Figure 5 — distributions across A−T−, A+T−, and A+T+ | `Fig5_AT_3group_distributions.py` |
| Figure 6 — clinical and molecular profiles of the exploratory A+T− groups | `Fig6_A+T-_subgroup_boxplots.py` |
| Supplementary Figure S6 — A+T− mixture-model analyses | `FigS6_plot_A+T-_integrated_figure.py` |
| Table 1 and formatted supplementary tables | `Table1.py`, `build_formatted_manuscript_tables.py` |

`Fig3_ROC.py` also generates Supplementary Figure S2. `Fig2Fig4_biomarker_distributions.py` generates Supplementary Figures S1 and S3–S5.

The supporting analysis scripts are:

- `01_diagnostic_validation.py`: cross-validation, cutoff stability, and paired AUC comparisons.
- `02a_bimodality_GMM_full_sample.py`: full-sample A+T− Gaussian mixture-model analysis and stability checks.
- `02b_bimodality_GMM_by_reference.py`: mixture-model sensitivity analyses stratified by T-reference method.
- `02c_unimodality_tests.py`: additional tests of unimodality.
- `03_Aplus_Tminus_subgroup_comparison.py`: comparisons of clinical, plasma, and available CSF measurements between the model-defined A+T− groups.

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
├── emtbr-tau243-immunoassay-vs-ptau217/
└── outputs/
```

Some scripts contain defaults specific to the original analysis computer. For scripts that support them, set `BIOMARKER_INPUT` and `BIOMARKER_OUTPUT_DIR` to the appropriate local paths before running. `Table1.py` reads `biomarker_result.xlsx` from the repository’s parent directory.

Run the principal analysis scripts before the figure and table scripts that use their outputs:

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

The figure scripts can then be run individually using the manuscript-to-script table above.

## Interpretation

The analyses assess discrimination within this study cohort. Sample-derived cutoffs were not externally validated and should not be used for clinical decisions. The A+T− mixture-model groups require independent and longitudinal validation.
