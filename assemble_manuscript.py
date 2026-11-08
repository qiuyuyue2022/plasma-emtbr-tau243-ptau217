"""Assemble the full manuscript from its section sources and number all citations.

Inputs (manuscript/ars/phase4_drafting/ unless noted):
  abstract_conclusion_v1.md                     Abstract, keywords, Conclusion
  introduction_discussion_v4_source_with_keys.md Introduction and Discussion with [@Key] citations
  manuscript.md                                  Methods and Results ([@Key] citations; stops before the author checklist)
  main_figure_legends_v2.md                      main figure legends
  ../../outputs/Table1_participant_characteristics.md
  declarations_v1.md                            title-page authors/affiliations and Declarations
  ../../references_db.json                       verified reference records keyed by citation key

Citations are numbered by first appearance in the order Introduction -> Methods -> Results ->
Discussion -> Conclusion -> legends; consecutive runs of three or more are collapsed (e.g., 2–4).
Outputs: manuscript/full_manuscript_v1.md and, if pandoc is available, full_manuscript_v1.docx.
"""
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "manuscript"
D = MS / "ars" / "phase4_drafting"
OUT_MD = MS / "full_manuscript_v1.md"
OUT_DOCX = MS / "full_manuscript_v1.docx"

TITLE = "A plasma eMTBR-tau243 immunoassay versus p-tau217 for tau status in a memory clinic"


def read(p):
    return Path(p).read_text(encoding="utf-8")


def between(text, start, end=None):
    i = text.index(start)
    j = text.index(end, i) if end else len(text)
    return text[i:j].strip()


def words(text):
    t = re.sub(r"^#+.*$", "", text, flags=re.M)
    t = re.sub(r"\[[0-9,–]+\]", "", t)
    return len(t.split())


ac = read(D / "abstract_conclusion_v1.md")
abstract = between(ac, "**INTRODUCTION:**", "**Keywords:**")
keywords = between(ac, "**Keywords:**", "---")
conclusion = between(ac, "## 5. Conclusion")

idsrc = read(D / "introduction_discussion_v4_source_with_keys.md")
intro = between(idsrc, "## 1. Introduction", "---")
discussion = between(idsrc, "## 4. Discussion")
discussion = discussion.split("\n---")[0].strip()

mr = read(D / "manuscript.md")
methods = between(mr, "## 2. Methods", "## 3. Results")
results = between(mr, "## 3. Results", "## Author-input checklist")

leg = read(D / "main_figure_legends_v2.md")
legends = leg[leg.index("## Figure 1."):].strip()
table1 = read(MS.parent / "outputs" / "Table1_participant_characteristics.md").strip()

refs = json.loads(read(MS / "references_db.json"))

decl_src = read(D / "declarations_v1.md")
title_block = between(decl_src, "<!-- TITLE-PAGE -->", "<!-- DECLARATIONS -->").replace("<!-- TITLE-PAGE -->", "").strip()
decl_block = between(decl_src, "<!-- DECLARATIONS -->").replace("<!-- DECLARATIONS -->", "").strip()

# ---- number citations by first appearance
body_parts = [intro, methods, results, discussion, conclusion, legends]
order = []


def fmt(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ",".join(map(str, nums[i:j + 1])))
        i = j + 1
    return ",".join(out)


def rep(m):
    keys = [k.strip().lstrip("@") for k in m.group(1).split(";")]
    nums = []
    for k in keys:
        if k not in refs:
            raise KeyError(f"citation key not in references_db.json: {k}")
        if k not in order:
            order.append(k)
        nums.append(order.index(k) + 1)
    return "[" + fmt(nums) + "]"


numbered = [re.sub(r"\[(@[^\]]+)\]", rep, part) for part in body_parts]
intro_n, methods_n, results_n, discussion_n, conclusion_n, legends_n = numbered


def flag_gaps(t):
    # make placeholders conspicuous in Word: `[MATERIAL GAP: ...]` -> **[MATERIAL GAP: ...]**
    return re.sub(r"`(\[MATERIAL GAP[^`]*\])`", r"**\1**", t)


methods_n = flag_gaps(methods_n)

counts = {
    "Abstract": words(re.sub(r"\*\*[A-Z]+:\*\*", "", abstract.split("**Trial registration:**")[0])),
    "Introduction": words(intro_n),
    "Methods": words(methods_n),
    "Results": words(results_n),
    "Discussion": words(discussion_n),
    "Conclusion": words(conclusion_n),
}
main_text = sum(counts[k] for k in ["Introduction", "Methods", "Results", "Discussion", "Conclusion"])
gaps = len(re.findall(r"MATERIAL GAP", methods_n + results_n))
confirms = len(re.findall(r"\[(?:CONFIRM|TO BE SUPPLIED)", decl_src.split("<!-- TITLE-PAGE -->")[1]))

title_page = f"""# {TITLE}

{title_block}

**Article type:** Research Article (*Alzheimer's & Dementia*)

**Word counts:** Abstract {counts['Abstract']} (limit 250); Introduction {counts['Introduction']} (limit 650); Methods {counts['Methods']}; Results {counts['Results']}; Discussion {counts['Discussion']} + Conclusion {counts['Conclusion']} = {counts['Discussion'] + counts['Conclusion']} (limit 1,500); main text {main_text:,}. References: {len(order)} (limit 75). Display items: 6 figures + 1 table (limit 8).

**Draft status:** assembled automatically by `code/assemble_manuscript.py`; {gaps} MATERIAL GAP placeholder(s) remain in Methods; declarations drafted from the ChP template with {confirms} [CONFIRM]/[TO BE SUPPLIED] item(s) to resolve.
"""

declarations = "## Declarations\n\n" + decl_block

ref_list = "## References\n\n" + "\n".join(f"{i + 1}. {refs[k]}" for i, k in enumerate(order))

doc = "\n\n".join([
    title_page,
    "## Abstract\n\n" + abstract,
    keywords,
    "---",
    intro_n, methods_n, results_n, discussion_n, conclusion_n,
    "---",
    declarations,
    ref_list,
    "---",
    "## Tables\n\n" + table1.replace("# Table 1.", "### Table 1.").replace("\n## Notes\n", "\n**Notes**\n"),
    "---",
    "## Figure legends\n\n" + legends_n.replace("## Figure", "### Figure"),
]) + "\n"
OUT_MD.write_text(doc, encoding="utf-8")
print(f"wrote {OUT_MD}")
print("counts:", counts, "| references:", len(order), "| material gaps:", gaps)

if shutil.which("pandoc"):
    subprocess.run(["pandoc", str(OUT_MD), "-o", str(OUT_DOCX)], check=True)
    print(f"wrote {OUT_DOCX}")
