# Biomarker A/T分析代码

本文件夹保存了本次分析使用的ROC曲线、AUC森林图、非参数分布图、GMM双峰性分析和结果表格代码。

**2026-11-01 全面重写说明：** 本README此前严重过时，列出的脚本文件名（如`Fig2_ROC.py`、`FigS3_plot_biomarker_distributions.py`、`01_validate_diagnostic_performance_feeds_FigTable.py`等）经过多轮改名/删除后已全部不存在。本次重写按`code/`目录当前的真实文件名、以及`manuscript/00_story_outline_and_analysis_lock.md`第10/11节当前锁定的Table/Figure编号重新核对。

## 定版图表编号对照（outputs/ 文件名 → 手稿编号，2026-11-05）

主文本 **6张图 + 1张主表 = 7个主显示项**。编号按正文首次引用顺序（Methods → Results → Discussion），2026-11-05已用脚本核对。编号的权威来源是`00_story_outline_and_analysis_lock.md`第10–11节；历次改号历史见其changelog（3.9–3.11、3.22–3.24）。被替换的旧图文件归档于`outputs/_superseded_20261105/`。

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

命名原则（2026-11-05）：直接产出某张图的出图脚本按图号命名（`Fig3_ROC.py`等）；分析脚本按内容命名、不带图表编号（`01_diagnostic_validation.py`、`02a/02b/02c_*`、`03_Aplus_Tminus_subgroup_comparison.py`），这样图表重新编号时脚本名不会过时。文件名中不使用空格。已过时的`04_sync_results_display_references_QC.py`移至`_obsolete/`，请勿运行。

## 文件说明（按依赖/运行顺序排列，而非文件名字母顺序）

**读原始数据、必须最先运行、后续脚本依赖其输出的脚本：**

1. `Fig3_ROC.py`
   - 读取原始Excel数据（`biomarker_result.xlsx`）。
   - 计算4个Aβ42标准化比值，对11个标志物在六组分析人群中计算AUC、3,000次患者水平Bootstrap 95% CI、Youden最佳切点，以及50次重复分层5折交叉验证的准确率/平衡准确率/敏感度/特异度。
   - **六组分析**：A全体(n=135)、A排除对照(n=101)、T-tauPET(n=68)、T-联合PET/CSF(n=114)、T-within-A+(n=67)、T-within-A+-tauPET-only(n=52)。
   - 输出`biomarker_performance_results.csv/json`（**`Fig2Fig4_biomarker_distributions.py`、`build_formatted_manuscript_tables.py`（TableS4，原S1）都依赖这份文件，必须最先运行**）、方法元数据，以及正文Figure 3（6面板ROC图）和补充Figure S2（AUC森林图，原编号S5，2026-11-01重新编号）。

2. `01_diagnostic_validation.py`
   - 独立读取原始Excel，对六组分析补充：训练折切点稳定性、固定切点的患者Bootstrap区间、4组预设配对AUC对比（DeLong + 配对Bootstrap，组内BH-FDR校正），以及post-hoc"最佳标志物 vs 其余十个"对比（各自单独做10重比较的BH-FDR校正）。
   - 输出`diagnostic_cv_performance.csv`、`diagnostic_cutoff_stability.csv`、`diagnostic_auc_pairwise_comparisons.csv`、`best_vs_others_AUC_comparisons.csv`（原名`TableS7_best_vs_others_AUC_comparisons.csv`，2026-11-01改名，是Supp. Table S2的原始数据而非独立编号表格）、`diagnostic_validation_summary.md/json`。

3. `02a_bimodality_GMM_full_sample.py`
   - 在A+T−人群中验证log10(eMTBR-tau243)的一/二成分Gaussian mixture model，p-tau217作为不参与分组的单峰对照。
   - 1,000次参数Bootstrap似然比检验、逐例剔除（LOO）、极端值/尾部敏感性分析、去标识化亚组归属稳定性。
   - 输出`bimodality_validation_summary.csv/md`、`bimodality_leave_one_out.csv`、`bimodality_extreme_value_sensitivity.csv`、`MTBR_subgroup_assignment_stability.csv`、`bimodality_validation_metadata.json`。**`Fig1_study_flow.py`、`FigS6_plot_A+T-_integrated_figure.py`、`build_formatted_manuscript_tables.py`（TableS6）依赖这些文件。**

4. `02b_bimodality_GMM_by_reference.py`
   - 在A+T−内按T参考标准分层：tau PET确认T−（n=18）和CSF定义T−（n=15），每层独立重复一/二成分GMM比较、1,000次参数Bootstrap、组间BH-FDR校正、LOO和极端值敏感性分析。
   - 输出`GMM_reference_stratified_*.csv/md/json`（**`FigS6_plot_A+T-_integrated_figure.py`、`build_formatted_manuscript_tables.py`（TableS6）依赖这些文件**）（原补充Figure S7已于2026-11-05删除，默认不再出图）。

5. `03_Aplus_Tminus_subgroup_comparison.py`
   - **第一部分**：在A+T−人群中，以呈双峰证据的log10(eMTBR-tau243)二成分GMM定义低/高eMTBR-tau243亚组（不重新拟合，复用第3号脚本的分组）；比较年龄、性别、APOE ε4、MMSE、MoCA-P，以及6项不参与分组的绝对值血浆生物标志物。
   - **第二部分**：在具有同期LUMIPULSE CSF结果的A+T−参与者（n=14）中，复用同一套低/高亚组，比较Aβ42/Aβ40、p-tau181、t-tau、p-tau181/Aβ42、t-tau/Aβ42，直接写出格式化的**Supplementary Table S8**（`write_manuscript_table_s6()`函数，2026-11-01新增，原编号S5，同日重新编号为S6）。
   - 输出模型诊断、临床连续/分类变量、血浆生物标志物（`A+T-Subgroup_*.csv`，**`build_formatted_manuscript_tables.py`构建TableS7（原S3）依赖这些文件**）、去标识化亚组成员、CSF比较结果（`A+T-_CSF_Lumipulse_*.csv/md`）和`TableS6_Aplus_Tminus_CSF_comparison.md`（原名`TableS5_Aplus_Tminus_CSF_comparison.md`）。脚本名中的"Table4"是历史遗留（该分析此前曾计划做主文本Table 4，后确认一直只在正文中以Supplementary Table引用，从未真正提升为主文本表，2026-11-01已将其产出文件改名为`TableS5_*`——先改名为`TableS3_*`，同日又因正文引用顺序核对改为`TableS5_*`，脚本文件名未再改动）。
   - **独立运行**，`FigS6_plot_A+T-_integrated_figure.py`依赖它的输出。

**不产出表格/图形、仅做格式化整理或校验的脚本：**

6. `Table1.py`
   - 读取原始Excel，按A/T分组计算参与者特征（年龄、性别、APOE、MMSE、MoCA-P、A/T参考方式构成）和7个绝对值生物标志物浓度分布。
   - 输出`Table1_participant_characteristics_draft.csv/md`、`Table1_metadata.json`、`Table1_workbook_data.json`（**`build_formatted_manuscript_tables.py`和`05_build_Table1_workbook_excel.mjs`依赖这份草稿/JSON**），以及**Supplementary Table S3**（`TableS2_biomarkers_by_AT_group.csv/md`，绝对浓度及完整A/T分布统计）。

7. `build_formatted_manuscript_tables.py`（2026-11-01由`Tables1to4_build_revised_results_tables.py`改名而来）
   - 读取上述脚本产出的草稿/锁定CSV，只做格式化整理，不重新计算任何统计量。
   - 产出定版**Table 1**（`Table1_participant_characteristics.md`，修正"A−T0"标签等文字问题）、**Supplementary Table S4**（原编号S1；11标志物×6人群完整表现网格，替换此前的孤儿文件`biomarker_performance_tables.md`——其原生成脚本`03_export_copyable_tables.py`早已被删除）、**Supplementary Table S7**（原编号S3）、**Supplementary Table S6**（编号未变；整合GMM全样本/分层诊断、Bootstrap、LOO、极端值和亚组归属稳定性，此前散落在多个未格式化的"分析护照"文件中，没有统一表格）。2026-11-01：`build_tableS1()`/`build_tableS3()`函数相应改名为`build_tableS3()`/`build_tableS5()`，匹配正文真实首次引用顺序重新编号后的表格。
   - 旧版本（`Tables1to4_build_revised_results_tables.py`）里的`build_table2()`/`build_table3()`对应已废弃的旧表格设计（完整11标志物主文本表、4组预设配对AUC对比表），已随改名一并删除；历史沿革见`00_story_outline_and_analysis_lock.md`第10节和changelog 3.4–3.7。

8. `05_build_Table1_workbook_excel.mjs`
   - 读取`Table1.py`产出的`Table1_workbook_data.json`，生成格式化的Table 1 Excel工作簿。

9. ~~`05_build_Table2_workbook_excel.mjs`~~ — 已于2026-11-05移至`_obsolete/`（生成的是已删除的旧主文本Table 2工作簿，且只覆盖5个分析人群）。
10. `02c_unimodality_tests.py`
    - 在A+T−的33名参与者中，对log10(eMTBR-tau243)和log10(p-tau217)做Hartigan's dip test（100,000次Bootstrap）和Silverman临界带宽检验（5,000次平滑Bootstrap）。
    - 输出`unimodality_test_results.csv`、`unimodality_silverman_bootstrap.csv`、`unimodality_test_summary.md/json`。**独立运行，结果以数值形式直接写入正文Section 3.6，不被其他脚本读取，也未被制成单独的表格。**

11. `Fig1_study_flow.py`
    - 从锁定数据程序化核对记忆门诊/健康对照来源、A/T状态、T参考标准和各分析人群人数，依赖第3号脚本（`02a_bimodality_GMM_full_sample.py`）的`bimodality_validation_summary.csv`。
    - 生成正文Figure 1研究队列与分析人群流程图：可编辑SVG/PDF、600 dpi TIFF、PNG、源数据CSV和元数据JSON。

12. `Fig2Fig4_biomarker_distributions.py`
    - 读取`Fig3_ROC.py`的`biomarker_performance_results.csv`。
    - 生成六组分析人群的11标志物非参数核密度分布图：正文Figure 2（A全体n=135，主分析）、Figure 4（T-tau PET-only），以及补充Figure S1（A排除对照n=101）/S3/S4/S5。是否显示Youden切点按文件名自动决定：正文图（stem以`Fig`开头）不显示，补充图（`FigS`）显示；源数据始终包含切点。原全局开关`DISTRIBUTION_MAIN_REVISION`已停用（2026-11-05，见lock文档changelog 3.22）。坐标轴标注单位（绝对浓度pg/mL，Aβ42比值为ratio）。

13. `Fig5_AT_3group_distributions.py`
    - 绘制A−T−、A+T−、A+T+三组的11项生物标志物分布图（正文Figure 5），A+T−与A+T+的Mann–Whitney U检验BH-FDR q值和秩二列效应量。
    - 导出PNG、PDF、SVG、600 dpi TIFF和统计结果CSV。

13b. `Fig6_A+T-_subgroup_boxplots.py`（2026-11-01新增）
    - 读取`A+T-Subgroup_membership.csv`（eMTBR-tau243 Low/High分组）并按SID关联`biomarker_result.xlsx`原始数据；读取`03_Aplus_Tminus_subgroup_comparison.py`产出的`A+T-Subgroup_continuous/categorical/biomarkers.csv`和`A+T-_CSF_Lumipulse_comparison.csv`中已锁定的FDR q值（本脚本不重新计算任何统计检验，scipy在部分环境下不可用）。
    - 绘制正文Figure 6三行箱型/柱状图：(a) 5项临床/人口学变量、(b) 6项非分组用血浆标志物、(c) n=14同期CSF子集的5项指标，均为线性坐标轴（经作者确认不用log坐标，因各指标Low/High合并取值范围均未超过约14倍）。图内只保留裸的a/b/c面板字母和各子图y轴标签，不包含说明文字或图注，完整说明放在`main_figure_legends_v2.md`的Figure 6图注中。
    - 导出SVG、PDF、300 dpi PNG、600 dpi TIFF，以及source data CSV（两份：主样本n=33、CSF子集n=14）和元数据JSON。

14. `FigS6_plot_A+T-_integrated_figure.py`
    - 将A+T−全样本GMM（第3号脚本）、按T参考标准分层的GMM稳健性（第4号脚本）、不参与分组的血浆指标效应量和同期CSF指标效应量（第5号脚本）整合为一张投稿级多面板图。
    - 输出补充Figure S6，附可编辑SVG/PDF、600 dpi TIFF、PNG、源数据CSV和元数据JSON。

**QC / 历史脚本：**

15. `04_audit_manuscript_results_QC.py`
    - 只读校验：核对锁定的AUC、交叉验证、GMM、血浆/CSF效应量和FDR q值是否与`methods_and_results_draft.md`、`00_story_outline_and_analysis_lock.md`中陈述的数值一致。任何不一致以非零状态退出。**每次改动正文数值或表格后都应重跑一遍。**

16. `04_sync_results_display_references_QC.py`
    - **历史性一次性脚本，不属于常规流水线，不要重新运行。** 这是2026-10-31前后用于把正文草稿从某个旧版本（硬编码"Version 1.5 → 1.6"）批量替换图表引用文字的脚本，字符串替换对象已是过时版本，现在的正文版本号是1.8，重新运行不会匹配到任何字符串，也不会产生任何效果。保留仅供查阅当时做了哪些文字替换。

## 其他说明

- `requirements.txt`：Python依赖。
- 所有脚本默认从`../biomarker_result.xlsx`读取原始数据、写入`../outputs/`，可通过环境变量`BIOMARKER_INPUT`/`BIOMARKER_OUTPUT_DIR`覆盖。
- 表格/图形的权威编号、状态和修订历史以`manuscript/00_story_outline_and_analysis_lock.md`第10/11节和changelog为准；本README只负责说明"每个编号对应哪个脚本/文件"，如有冲突以lock文档为准。
