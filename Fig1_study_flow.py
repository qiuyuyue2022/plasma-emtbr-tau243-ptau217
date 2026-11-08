#!/usr/bin/env python3
"""Create the manuscript study-cohort and analysis-population flow figure."""

from __future__ import annotations

import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))

import matplotlib as mpl
mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd


DEFAULT_INPUT = ROOT / "biomarker_result.xlsx"
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(DEFAULT_INPUT)))
OUTPUT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(ROOT / "outputs")))
OUTPUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.8,
})


PITCH = 0.036  # absolute line pitch (data units), independent of box height


def flow_box(ax, x, y, w, h, label, number, *, detail=None, detail2=None, lw=0.8):
    """Rounded box whose text block (title / n / optional detail lines) is
    vertically centered using a fixed absolute line pitch, so centering
    holds regardless of box height. No explanatory/caption text goes in a
    box; only flow labels and counts (figure legend carries the prose)."""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.002,rounding_size=0.012",
        facecolor="#FFFFFF", edgecolor="#000000", linewidth=lw, zorder=2))
    lines = [(label, 7.6, False), (f"n = {number}", 7.4, True)]
    if detail:
        lines.append((detail, 6.6, False))
    if detail2:
        lines.append((detail2, 6.6, False))
    k = len(lines)
    y_mid = y + h / 2
    for i, (text, fs, bold) in enumerate(lines):
        offset = (k - 1) / 2 - i
        ax.text(x + w / 2, y_mid + offset * PITCH, text, ha="center", va="center",
                fontsize=fs, fontweight=("bold" if bold else "normal"),
                color="#000000", zorder=3)


def line(ax, *points, lw=0.8):
    xs, ys = zip(*points)
    ax.plot(xs, ys, color="#000000", linewidth=lw, solid_capstyle="butt", zorder=1)


def count_data():
    data = pd.read_excel(INPUT)
    required = {"SID", "DIAG", "A", "T", "taupet", "A_level"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if len(data) != data["SID"].nunique(dropna=True):
        raise AssertionError("Rows are not one-to-one with nonmissing SID values")

    health = data["DIAG"].astype(str).eq("健康对照")
    t_available = data["T"].isin(["T+", "T-"])
    tau_pet = data["taupet"].eq(1)
    a_plus = data["A"].eq("A+")
    t_minus = data["T"].eq("T-")

    counts = {
        "included": len(data),
        "memory_clinic": int((~health).sum()),
        "healthy_controls": int(health.sum()),
        "a_plus": int(a_plus.sum()),
        "a_minus": int(data["A"].eq("A-").sum()),
        "a_noncontrol_n": int((~health).sum()),
        "a_noncontrol_plus": int((~health & a_plus).sum()),
        "a_noncontrol_minus": int((~health & data["A"].eq("A-")).sum()),
        "t_available": int(t_available.sum()),
        "t_plus": int(data["T"].eq("T+").sum()),
        "t_minus": int(t_minus.sum()),
        "t0": int(data["T"].eq("T0").sum()),
        "t0_a_plus": int((data["T"].eq("T0") & data["A"].eq("A+")).sum()),
        "t0_a_minus": int((data["T"].eq("T0") & data["A"].eq("A-")).sum()),
        "tau_pet_n": int(tau_pet.sum()),
        "tau_pet_plus": int((tau_pet & data["T"].eq("T+")).sum()),
        "tau_pet_minus": int((tau_pet & t_minus).sum()),
        "csf_defined_tminus": int((~tau_pet & t_minus).sum()),
        "a_minus_tminus": int((data["A"].eq("A-") & t_minus).sum()),
        "a_plus_tminus": int((a_plus & t_minus).sum()),
        "a_plus_tplus": int((a_plus & data["T"].eq("T+")).sum()),
        "a_plus_t_available": int((a_plus & t_available).sum()),
        "a_plus_tau_pet": int((a_plus & tau_pet).sum()),
        "a_plus_tau_pet_tminus": int((a_plus & tau_pet & t_minus).sum()),
        "a_plus_csf_tminus": int((a_plus & ~tau_pet & t_minus).sum()),
        "patient_a_pet": int((~health & data["A_level"].eq("PET")).sum()),
        "patient_a_csf": int((~health & data["A_level"].eq("CSF")).sum()),
        "control_a_pet": int((health & data["A_level"].eq("PET")).sum()),
        "control_a_csf": int((health & data["A_level"].eq("CSF")).sum()),
        "control_a_plasma": int((health & data["A_level"].eq("Plasma")).sum()),
    }

    expected = {
        "included": 135, "memory_clinic": 101, "healthy_controls": 34,
        "a_plus": 67, "a_minus": 68, "t_available": 114,
        "t_plus": 34, "t_minus": 80, "t0": 21,
        "t0_a_plus": 0, "t0_a_minus": 21,
        "tau_pet_n": 68, "tau_pet_plus": 34, "tau_pet_minus": 34,
        "csf_defined_tminus": 46, "a_minus_tminus": 47,
        "a_plus_tminus": 33, "a_plus_tplus": 34,
        "a_plus_tau_pet": 52, "a_plus_tau_pet_tminus": 18,
        "a_plus_csf_tminus": 15,
        "patient_a_pet": 91, "patient_a_csf": 10,
        "control_a_pet": 15, "control_a_csf": 7, "control_a_plasma": 12,
    }
    drift = {key: (counts[key], value) for key, value in expected.items()
             if counts[key] != value}
    if drift:
        raise AssertionError(f"Locked cohort counts changed: {drift}")

    gmm = pd.read_csv(OUTPUT / "bimodality_validation_summary.csv")
    mtbr = gmm.loc[gmm["Marker"].eq("MTBR-Tau243")]
    if len(mtbr) != 1:
        raise AssertionError("Expected one locked MTBR-Tau243 GMM row")
    counts["mtbr_low"] = int(mtbr.iloc[0]["Observed low n"])
    counts["mtbr_high"] = int(mtbr.iloc[0]["Observed high n"])
    if (counts["mtbr_low"], counts["mtbr_high"]) != (12, 21):
        raise AssertionError("Locked full-sample eMTBR group sizes changed")
    return counts


def draw(counts):
    """Left-aligned recruitment-to-A/T spine with a right-hand exclusion
    branch for participants without an available T status, filling the
    full canvas width rather than a narrow centered column.

    Patients and controls are kept as two parallel columns through the
    A-reference step (rather than merging immediately), each carrying its
    own A-status split and A-reference-method breakdown, so the
    patient/control asymmetry described in Methods 2.1-2.2 (patients:
    PET/CSF; controls: PET/CSF/plasma, always A-negative) is visible in
    the figure itself. All explanatory prose (which group is all
    controls, the plasma-assay distinction) lives in the Figure 1 legend
    in the manuscript text, not inside the figure.
    """
    fig = plt.figure(figsize=(7.4, 4.55), facecolor="white")
    ax = fig.add_axes([0.03, 0.03, 0.95, 0.95])
    ax.set_xlim(0, 1)
    ax.set_ylim(0.03, 1)
    ax.axis("off")

    main_x, main_w = 0.05, 0.60
    main_cx = main_x + main_w / 2  # spine centerline

    # Row 1: recruitment sources, each box carrying its own A-status split
    # and A-reference-method breakdown (unified "X reference: … | …" form).
    rec_w = 0.26
    rec1_x = main_x
    rec2_x = main_x + main_w - rec_w
    row1_y, row1_h = 0.745, 0.22
    flow_box(ax, rec1_x, row1_y, rec_w, row1_h, "Memory-clinic participants",
             counts["memory_clinic"],
             detail=f'A+ {counts["a_noncontrol_plus"]}  |  A− {counts["a_noncontrol_minus"]}',
             detail2=f'A reference: PET {counts["patient_a_pet"]} | CSF {counts["patient_a_csf"]}',
             lw=1.0)
    flow_box(ax, rec2_x, row1_y, rec_w, row1_h, "HC", counts["healthy_controls"],
             detail=f'A− {counts["healthy_controls"]}',
             detail2=(f'A reference: PET {counts["control_a_pet"]} | '
                      f'CSF {counts["control_a_csf"]} | plasma {counts["control_a_plasma"]}'),
             lw=1.0)
    line(ax, (rec1_x + rec_w / 2, row1_y), (rec1_x + rec_w / 2, 0.68),
         (rec2_x + rec_w / 2, 0.68), (rec2_x + rec_w / 2, row1_y))
    line(ax, (main_cx, 0.68), (main_cx, 0.64))

    # Row 2: final included cohort, carrying the A-status summary.
    row2_y = 0.505
    flow_box(ax, main_x, row2_y, main_w, 0.135, "Final included cohort",
             counts["included"],
             detail=f'A+ {counts["a_plus"]}  |  A− {counts["a_minus"]}', lw=1.0)
    line(ax, (main_cx, row2_y), (main_cx, 0.455))

    # Branch point: continue down the spine to T-status assessment, and
    # branch right to the participants excluded from T-status analyses.
    branch_y = 0.455
    excl_x, excl_w = main_x + main_w + 0.07, 1.0 - (main_x + main_w + 0.07) - 0.02
    excl_cy = 0.41
    line(ax, (main_cx, branch_y), (main_cx, 0.43))
    line(ax, (main_cx, 0.43), (excl_x + excl_w / 2, 0.43),
         (excl_x + excl_w / 2, excl_cy + 0.0475))
    flow_box(ax, excl_x, excl_cy, excl_w, 0.095,
             "T status unavailable", counts["t0"])
    line(ax, (main_cx, 0.43), (main_cx, 0.39))

    # Row 3: T-status available, continuing the main spine. Same unified
    # "T reference: … | …" form as the A-reference lines above.
    row3_y = 0.255
    flow_box(ax, main_x, row3_y, main_w, 0.135, "T status available",
             counts["t_available"],
             detail=f'T reference: PET {counts["tau_pet_n"]} | CSF {counts["csf_defined_tminus"]}',
             lw=1.0)

    # Row 4: the three A/T groups, outer edges aligned with the main spine.
    split_y = row3_y
    line(ax, (main_cx, split_y), (main_cx, 0.205))
    group_y = 0.045
    group_w = (main_w - 2 * 0.025) / 3
    xs = [main_x, main_x + group_w + 0.025, main_x + 2 * (group_w + 0.025)]
    for gx in xs:
        line(ax, (gx + group_w / 2, 0.205), (gx + group_w / 2, group_y + 0.085))
    line(ax, (xs[0] + group_w / 2, 0.205), (xs[2] + group_w / 2, 0.205))
    for x, label, n in [
        (xs[0], "A−T−", counts["a_minus_tminus"]),
        (xs[1], "A+T−", counts["a_plus_tminus"]),
        (xs[2], "A+T+", counts["a_plus_tplus"]),
    ]:
        flow_box(ax, x, group_y, group_w, 0.085, label, n)

    return fig


def export_source_data(counts):
    rows = [
        ("Recruitment", "Memory-clinic participants", counts["memory_clinic"]),
        ("Recruitment", "Healthy controls", counts["healthy_controls"]),
        ("Patient A status", "A+", counts["a_noncontrol_plus"]),
        ("Patient A status", "A−", counts["a_noncontrol_minus"]),
        ("Patient A reference", "PET", counts["patient_a_pet"]),
        ("Patient A reference", "CSF", counts["patient_a_csf"]),
        ("Control A status", "A−", counts["healthy_controls"]),
        ("Control A reference", "PET", counts["control_a_pet"]),
        ("Control A reference", "CSF", counts["control_a_csf"]),
        ("Control A reference", "Plasma", counts["control_a_plasma"]),
        ("Included cohort", "Final included cohort", counts["included"]),
        ("A status", "A+", counts["a_plus"]),
        ("A status", "A−", counts["a_minus"]),
        ("T status", "T+", counts["t_plus"]),
        ("T status", "T−", counts["t_minus"]),
        ("T status", "T0 (unavailable)", counts["t0"]),
        ("T reference", "Tau PET", counts["tau_pet_n"]),
        ("T reference", "CSF-defined T−", counts["csf_defined_tminus"]),
        ("AT group", "A−T−", counts["a_minus_tminus"]),
        ("AT group", "A+T−", counts["a_plus_tminus"]),
        ("AT group", "A+T+", counts["a_plus_tplus"]),
        ("A+T− reference stratum", "Tau PET-confirmed T−", counts["a_plus_tau_pet_tminus"]),
        ("A+T− reference stratum", "CSF-defined T−", counts["a_plus_csf_tminus"]),
        ("A+T− eMTBR component", "Low", counts["mtbr_low"]),
        ("A+T− eMTBR component", "High", counts["mtbr_high"]),
    ]
    pd.DataFrame(rows, columns=["Section", "Population", "N"]).to_csv(
        OUTPUT / "Fig1_study_flow_SourceData.csv", index=False)

def main():
    counts = count_data()
    fig = draw(counts)
    stem = OUTPUT / "Fig1_study_flow"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    export_source_data(counts)
    print(f"Study-flow figure and source data written to {OUTPUT}")


if __name__ == "__main__":
    main()

