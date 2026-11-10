#!/usr/bin/env python3
"""Renumber numbered citations in the canonical manuscript by first appearance.

Added 2026-11-09 (Stage 4'). Reads manuscript/full_manuscript_v1.md, maps each
existing [n] to a references_db.json key via the current numbered reference list,
accepts newly inserted [@Key; @Key] citations, renumbers all citations in order of
first appearance in the body (text before "## References"), and rebuilds the list
from references_db.json. Writes the file in place; keep a backup before running.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "manuscript" / "full_manuscript_v1.md"
DB = json.loads((ROOT / "manuscript" / "references_db.json").read_text(encoding="utf-8"))

text = MS.read_text(encoding="utf-8")
body, rest = text.split("## References\n", 1)
ref_block, tail = rest.split("\n---\n", 1)
entries = dict((int(n), e.strip()) for n, e in re.findall(r"^(\d+)\. (.+)$", ref_block, re.M))


def key_for(entry: str) -> str:
    doi = entry.split("doi:")[-1].strip() if "doi:" in entry else None
    hits = [k for k, v in DB.items() if (doi and v.endswith(doi)) or v == entry]
    if len(hits) != 1:
        raise ValueError(f"cannot map reference to one key: {entry[:80]} -> {hits}")
    return hits[0]


num2key = {n: key_for(e) for n, e in entries.items()}


def expand(group: str):
    out = []
    for part in group.split(","):
        b = re.split(r"[–-]", part.strip())
        out += [int(b[0])] if len(b) == 1 else list(range(int(b[0]), int(b[1]) + 1))
    return out


# 1) numbers -> keys
# A citation bracket follows text or punctuation directly; "26 [21–27]" (space before) is an IQR/CI, not a citation.
body = re.sub(r"(?<!\s)\[(\d+(?:\s*[,–-]\s*\d+)*)\]",
              lambda m: "[" + "; ".join("@" + num2key[n] for n in expand(m.group(1))) + "]", body)
# 2) first-appearance order
order = []
for grp in re.findall(r"\[(@[^\]]+)\]", body):
    for k in [x.strip().lstrip("@") for x in grp.split(";")]:
        if k not in DB:
            raise KeyError(f"unknown citation key {k}")
        if k not in order:
            order.append(k)
new = {k: i + 1 for i, k in enumerate(order)}


def fmt(nums):
    nums = sorted(set(nums)); runs = []; s = p = nums[0]
    for x in nums[1:]:
        if x == p + 1:
            p = x; continue
        runs.append((s, p)); s = p = x
    runs.append((s, p))
    return ",".join(f"{a}" if a == b else (f"{a},{b}" if b == a + 1 else f"{a}–{b}") for a, b in runs)


body = re.sub(r"\[(@[^\]]+)\]", lambda m: "[" + fmt([new[x.strip().lstrip("@")] for x in m.group(1).split(";")]) + "]", body)
refs = "\n".join(f"{i + 1}. {DB[k]}" for i, k in enumerate(order))
MS.write_text(body + "## References\n\n" + refs + "\n\n---\n" + tail, encoding="utf-8")
print(f"{len(order)} references; previously {len(entries)}")
