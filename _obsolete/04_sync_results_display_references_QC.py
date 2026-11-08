#!/usr/bin/env python3
"""Synchronize Results figure/table citations with the author-directed display map."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/ars/phase4_drafting/methods_and_results_draft.md"
text = MANUSCRIPT.read_text(encoding="utf-8")

replacements = [
    ("**Version:** 1.5", "**Version:** 1.6"),
    ("summarized in Figure 1. Participant characteristics", "summarized in Supplementary Figure S1. Participant characteristics"),
    ("Among all 135 participants, eMTBR-tau243/Aβ42", "Among all 135 participants (Figure 1; Table 2), eMTBR-tau243/Aβ42"),
    ("After healthy controls were excluded, 101 participants remained", "After healthy controls were excluded, 101 participants remained"),
    ("In the primary tau PET-only cohort (n=67; 34 T+ and 33 T−),", "In the primary tau PET-only cohort (n=67; 34 T+ and 33 T−; Figures 2 and 4c; Tables 2 and 3),"),
    ("The broader combined tau PET/CSF analysis included 114 participants (34 T+ and 80 T−).", "The broader combined tau PET/CSF analysis included 114 participants (34 T+ and 80 T−; Figures 3 and 4d; Tables 2 and 3)."),
    ("Among the 67 A+ participants, including 34 T+ and 33 T−,", "Among the 67 A+ participants, including 34 T+ and 33 T− (Figure 4e; Tables 2 and 3),"),
    ("Visual inspection of the nonparametric density curves suggested that the A+T− eMTBR-tau243 distribution contained separated low and high concentration regions. This pattern motivated the exploratory, result-driven mixture-model analysis.",
     "Within A+T−, the individual eMTBR-tau243 values showed separated low and high concentration regions (Figure 5a). This descriptive pattern motivated the exploratory, result-driven mixture-model analysis."),
    ("Within the 33 A+T− participants, the two-component Gaussian mixture model", "Within the 33 A+T− participants (Figure 5a), the two-component Gaussian mixture model"),
    ("The distributional pattern recurred within the 18 tau PET-confirmed T− participants.", "The distributional pattern recurred within the 18 tau PET-confirmed T− participants (Figure 5b)."),
    ("Compared with the low-eMTBR-tau243 group, the high-eMTBR-tau243 group had higher p-tau217 concentrations", "As shown in Figure 5c and Supplementary Table S3, the high-eMTBR-tau243 group had higher p-tau217 concentrations than the low group"),
    ("Fourteen A+T− participants had contemporaneous LUMIPULSE CSF measurements,", "Fourteen A+T− participants had contemporaneous LUMIPULSE CSF measurements (Figure 5d; Supplementary Table S5),"),
]

for before, after in replacements:
    count = text.count(before)
    if count != 1:
        raise AssertionError(f"Expected one occurrence of {before!r}, found {count}")
    text = text.replace(before, after, 1)

MANUSCRIPT.write_text(text, encoding="utf-8")
print(f"Updated {len(replacements)} display references in {MANUSCRIPT}")
