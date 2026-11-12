"""Assemble the Supporting Information file from the current outputs.

Inputs:
  outputs/TableS1_*.md ... outputs/TableS9_*.md            supplementary tables (generated or hand-built)
  outputs/FigS1_*.png ... outputs/FigS7_*.png              supplementary figures
  manuscript/ars/phase4_drafting/supplementary_figure_legends_v1.md   figure titles and legends
  manuscript/full_manuscript_v1.md                         title (first line)

Nothing is recomputed: tables are copied verbatim (headings demoted one level), figures are embedded
with their legends. Outputs: manuscript/supporting_information_v1.md and, if pandoc is available,
supporting_information_v1.docx.
"""
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
MS = ROOT / "manuscript"
LEGENDS = MS / "ars" / "phase4_drafting" / "supplementary_figure_legends_v1.md"
OUT_MD = MS / "supporting_information_v1.md"
OUT_DOCX = MS / "supporting_information_v1.docx"

TABLES = [
    "TableS1_plasma_assay_parameters.md",
    "TableS2_best_marker_vs_others.md",
    "TableS3_biomarkers_by_AT_group.md",
    "TableS4_biomarker_performance_full.md",
    "TableS5_AB42_LLOQ_sensitivity.md",
    "TableS6_incremental_value.md",  # 2026-11-09 revision
    "TableS7_GMM_sensitivity_analyses.md",
    "TableS8_Aplus_Tminus_groups.md",  # tau PET-defined A+T- (n=18); CSF table removed
    "TableS9_PETCSF_Aplus_Tminus.md",  # v3.62: combined tau PET/CSF A+T- (n=33), restored
]
FIGURES = [
    "FigS1_distribution_A_all.png",
    "FigS2_full_AUC_forest_11marker.png",
    "FigS3_distribution_T_within_Aplus_tauPET.png",  # was S4 before v3.60
    "FigS4_AT_3group_PET_CSF.png",  # v3.61: combined tau PET/CSF A/T groups (restored)
    "FigS5_Aplus_Tminus_heterogeneity.png",  # S6 before v3.60; S4 in v3.60
    "FigS6_PETCSF_Aplus_Tminus_bimodality.png",  # v3.62: combined tau PET/CSF A+T- mixture models (restored)
    "FigS7_PETCSF_Aplus_Tminus_subgroup_boxplots.png",  # v3.62: combined tau PET/CSF A+T- subgroup profiles (restored)
]


def read(p):
    return Path(p).read_text(encoding="utf-8")


PAGE = '```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```'

title = read(MS / "full_manuscript_v1.md").splitlines()[0].lstrip("# ").strip()

# ---- tables
table_blocks, table_titles = [], []
for i, name in enumerate(TABLES, 1):
    text = read(OUT / name).strip()
    first, rest = text.split("\n", 1)
    if not first.startswith(f"# Supplementary Table S{i}."):
        raise ValueError(f"{name}: expected heading 'Supplementary Table S{i}.', found {first!r}")
    table_titles.append(first.lstrip("# ").strip())
    table_blocks.append("### " + first.lstrip("# ").strip() + "\n" + rest.strip())

# ---- figures
leg = read(LEGENDS)
fig_blocks, fig_titles = [], []
for i, img in enumerate(FIGURES, 1):
    m = re.search(rf"^## (Supplementary Figure S{i}\..*?)\n\n(.*?)(?=\n## |\Z)", leg, flags=re.S | re.M)
    if not m:
        raise ValueError(f"legend for Supplementary Figure S{i} not found")
    if not (OUT / img).exists():
        raise FileNotFoundError(OUT / img)
    head, body = m.group(1).strip(), m.group(2).strip()
    fig_titles.append(head)
    fig_blocks.append(f"### {head}\n\n![](../outputs/{img}){{width=100%}}\n\n{body}")

contents = "\n".join(f"- {t}" for t in table_titles + fig_titles)

doc = "\n\n".join([
    "# Supporting Information",
    f"**{title}**",
    "## Contents\n\n" + contents,
    PAGE,
    "## Supplementary Tables",
    ("\n\n" + PAGE + "\n\n").join(table_blocks),
    PAGE,
    "## Supplementary Figures",
    ("\n\n" + PAGE + "\n\n").join(fig_blocks),
]) + "\n"

# ---- typographic normalization (2026-11-09; formatting only, no values change)
# 1) hyphen-minus before a number -> true minus sign, as in the main text; negative zero -> 0
doc = re.sub(r"(?<=[\s|(\[;/])-(?=\d)", "−", doc)
doc = re.sub(r"−(0\.0+)(?![\d])", r"\1", doc)


# 2) Table S8 plasma rows: restore trailing zeros so values have 3 significant figures (6.7 -> 6.70)
def _sig3(num: str) -> str:
    digits = re.sub(r"[^\d]", "", num).lstrip("0")
    if len(digits) >= 3:
        return num
    pad = 3 - len(digits)
    return num + ("" if "." in num else ".") + "0" * pad


S8_PLASMA = ("p-tau217", "BD-tau", "p-tau231", "GFAP", "Aβ1-42", "NfL")


def _fix_s8(block: str) -> str:
    out = []
    for line in block.split("\n"):
        if any(line.startswith(f"| {m} |") for m in S8_PLASMA):
            cells = line.split(" | ")
            for k in (1, 2):
                cells[k] = re.sub(r"\d+(?:\.\d+)?", lambda mm: _sig3(mm.group(0)), cells[k].split("; n=")[0]) + \
                    "; n=" + cells[k].split("; n=")[1]
            line = " | ".join(cells)
        out.append(line)
    return "\n".join(out)


s8_start = doc.index("### Supplementary Table S8.")
s8_end = doc.index("## Supplementary Figures")
doc = doc[:s8_start] + _fix_s8(doc[s8_start:s8_end]) + doc[s8_end:]
doc = doc.replace("Contemporaneous CSF was available for 6 participants (5 Low, 1 High)",
                  "Values are median [interquartile range] or n (%). "
                  "Contemporaneous CSF was available for 6 participants (5 low, 1 high)")

# The table footnotes refer to the main-text tables/figures; flag any leftover internal file names.
leftover = re.findall(r"[\w\-+]+\.(?:csv|json|py)\b", doc)
if leftover:
    print("WARNING: internal file names in Supporting Information:", sorted(set(leftover)))

OUT_MD.write_text(doc, encoding="utf-8")
print(f"wrote {OUT_MD}  ({len(TABLES)} tables, {len(FIGURES)} figures)")

if shutil.which("pandoc"):
    subprocess.run(["pandoc", str(OUT_MD), "--resource-path", str(MS), "-o", str(OUT_DOCX)], check=True)
    print(f"wrote {OUT_DOCX}")
