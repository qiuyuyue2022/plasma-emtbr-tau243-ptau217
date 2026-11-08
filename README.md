# Biomarker A/T分析代码

本文件夹保存了本次分析使用的ROC曲线、AUC森林图、非参数分布图、GMM双峰性分析和结果表格代码。

| 手稿编号 | outputs/ 文件名前缀 | 生成脚本 |
| --- | --- | --- |
| Figure 1 | `Fig1_study_flow` | `Fig1_study_flow.py` |
| Figure 2 | `Fig2_distribution_A_noncontrol`（A状态，排除健康对照n=101，主分析；不显示切点） | `Fig2Fig4_biomarker_distributions.py` |
| Figure 3 | `Fig3_full_ROC_11marker`（6面板ROC） | `Fig3_ROC.py` |
| Figure 4 | `Fig4_distribution_T_tauPET`（T状态，tau PET-only n=68；不显示切点） | `Fig2Fig4_biomarker_distributions.py` |
| Figure 5 | `Fig5_AT_3group_distributions` | `Fig5_AT_3group_distributions.py` |
| Figure 6 | `Fig6_A+T-_subgroup_boxplots` | `Fig6_A+T-_subgroup_boxplots.py`（读取`03_compare_…py`已锁定的q值，不重新计算统计检验） |
| Table 1 | `Table1_participant_characteristics.md` | `Table1.py`（草稿）→ `build_formatted_manuscript_tables.py`（定版）→ `05_build_Table1_workbook_excel.mjs`（Excel） |
| Supp. Figure S1 | `FigS1_distribution_A_all`（全体n=135，次要分析；显示切点） | `Fig2Fig4_biomarker_distributions.py` |
| Supp. Figure S2 | `FigS2_full_AUC_forest_11marker` | `Fig3_ROC.py` |
| Supp. Figure S3 | `FigS3_distribution_T_PET_CSF` | `Fig2Fig4_biomarker_distributions.py` |
| Supp. Figure S4/S5 | `FigS4/5_distribution_T_within_Aplus_*` | `Fig2Fig4_biomarker_distributions.py` |
| Supp. Figure S6 | `FigS6_Aplus_Tminus_heterogeneity` | `FigS6_plot_A+T-_integrated_figure.py` |
| ~~Supp. Figure S7~~ | 已删除（2026-11-05，与S6b重复） | `02b_bimodality_GMM_by_reference.py`默认不再出图；设`BIOMARKER_PLOT_S7=1`可重绘 |
| Supp. Table S1 | `TableS1_plasma_assay_parameters.md` | 人工维护（iomics试剂参数表，2026-11-04） |
| Supp. Table S2 | `TableS2_best_marker_vs_others.md` | `build_formatted_manuscript_tables.py`的`build_tableS2()`（2026-11-05起由脚本生成），数据来自`01_diagnostic_validation.py`的`best_vs_others_AUC_comparisons.csv` |
| Supp. Table S3 | `TableS3_biomarkers_by_AT_group.{md,csv}` | `Table1.py` |
| Supp. Table S4 | `TableS4_biomarker_performance_full.md` | `build_formatted_manuscript_tables.py`（数据来自`Fig3_ROC.py`） |
| Supp. Table S5 | `TableS5_AB42_LLOQ_sensitivity.md` | `build_sensitivity_AB42_LLOQ_table.py`（Aβ1-42截断敏感性分析；输入在`outputs/sensitivity_AB42_LLOQ_censoring/`） |
| Supp. Table S6 | `TableS6_GMM_sensitivity_analyses.md` | `build_formatted_manuscript_tables.py` |
| Supp. Table S7 | `TableS7_Aplus_Tminus_groups.md` | `build_formatted_manuscript_tables.py`（数据来自`03_compare_…py`） |
| Supp. Table S8 | `TableS8_Aplus_Tminus_CSF_comparison.md` | `03_Aplus_Tminus_subgroup_comparison.py`（`write_manuscript_table_s6()`） |

