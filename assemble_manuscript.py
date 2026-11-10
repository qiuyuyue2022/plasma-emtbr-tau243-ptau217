"""Export the canonical manuscript Markdown to Word.

Edit manuscript/full_manuscript_v1.md directly. This script reads that file and
never rewrites it. The historical section drafts are not inputs.
"""

from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript" / "full_manuscript_v1.md"
DOCX = MANUSCRIPT.with_suffix(".docx")


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    if "[@" in text:
        raise ValueError("Unresolved citation keys in the canonical manuscript")

    references = text.split("## References\n", 1)
    if len(references) != 2:
        raise ValueError("Missing References section")
    numbered = re.findall(r"^(\d+)\. ", references[1].split("\n---\n", 1)[0], re.M)
    if [int(n) for n in numbered] != list(range(1, len(numbered) + 1)):
        raise ValueError("Reference list numbering is not consecutive")

    cited = set()
    for group in re.findall(r"\[(\d+(?:[,–-]\d+)*)\]", references[0]):
        for part in group.split(","):
            bounds = re.split(r"[–-]", part)
            if len(bounds) == 1:
                cited.add(int(bounds[0]))
            else:
                cited.update(range(int(bounds[0]), int(bounds[1]) + 1))
    expected = set(range(1, len(numbered) + 1))
    if cited != expected:
        raise ValueError(f"Citation/reference mismatch: missing {sorted(expected - cited)}, unknown {sorted(cited - expected)}")

    pandoc = shutil.which("pandoc")
    if pandoc is None:
        raise RuntimeError("pandoc is required to export the Word file")
    subprocess.run([pandoc, str(MANUSCRIPT), "-o", str(DOCX)], check=True)
    print(f"Source: {MANUSCRIPT}")
    print(f"Exported: {DOCX}")
    print(f"References: {len(numbered)}")


if __name__ == "__main__":
    main()
