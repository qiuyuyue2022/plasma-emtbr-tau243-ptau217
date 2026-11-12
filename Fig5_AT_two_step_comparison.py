#!/usr/bin/env python3
"""Two-step A/T comparison in tau PET-defined groups (Stage 4' storyline, 2026-11-09).

Groups (tau PET-classified T status only): A−T−, A+T−, A+T+.
For each of the seven absolute plasma markers:
  * amyloid step  A−T− vs A+T−  (Mann–Whitney U, rank-biserial r, BH-FDR across 7 markers)
  * tau step      A+T− vs A+T+  (same)
  * A-status AUC within tau PET-defined T− (A+T− vs A−T−) with a 3,000-sample participant
    bootstrap 95% CI, and the paired-bootstrap ΔAUC for eMTBR-tau243 minus p-tau217.
Outputs: outputs/AT_two_step_comparison.{csv,md} and the figure Fig5_AT_two_step.{png,pdf,svg,tiff}
(panel a: three-group distributions per marker; panel b: amyloid-step |r| vs tau-step |r|).
Read-only with respect to the raw data.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
INPUT = Path(os.environ.get("BIOMARKER_INPUT", str(ROOT / "biomarker_result.xlsx")))
OUT = Path(os.environ.get("BIOMARKER_OUTPUT_DIR", str(ROOT / "outputs")))
SEED = int(os.environ.get("BIOMARKER_RANDOM_SEED", "20261026"))
N_BOOT = 3000
STEM = "Fig5_AT_two_step"

MARKERS = [  # display name, source column
    ("eMTBR-tau243", "MTBR-Tau243"),
    ("p-tau217", "p-Tau217(≤0.758)"),
    ("p-tau231", "p-Tau231(≤5.793)"),
    ("BD-tau", "BD Tau（≤52.301）"),
    ("GFAP", "GFAP(≤131.880)"),
    ("NfL", "NFL(≤9.624 )"),
    ("Aβ1-42", "Aβ 1-42"),
]
GROUPS = [("A−T−", "#2F6DB5"), ("A+T−", "#E69F00"), ("A+T+", "#C23B53")]
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#d9d9d9"


def clean(s: pd.Series) -> pd.Series:
    return s.astype("string").str.strip().str.replace("−", "-", regex=False)


def bh(p) -> np.ndarray:
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def rb(x, y):
    """Mann–Whitney p and rank-biserial r (positive = higher in y)."""
    u = mannwhitneyu(y, x, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1, u.pvalue


def nice_log_ticks(lo: float, hi: float):
    """1-2-5 ticks (raw scale) inside [10**lo, 10**hi]; thinned to at most 5."""
    cand = [m * 10.0 ** e for e in range(-3, 5) for m in (1, 2, 5)]
    ticks = [c for c in cand if 10 ** lo <= c <= 10 ** hi]
    if len(ticks) > 5:
        ticks = [c for c in ticks if str(c).lstrip("0.").startswith(("1", "1e"))] or ticks[::2]
    return ticks


def tick_label(v: float) -> str:
    return f"{v:g}" if v < 1000 else f"{int(v):,}"


def qtext(q: float) -> str:
    return "q<0.001" if q < 0.001 else f"q={q:.3f}"


def main() -> None:
    d = pd.read_excel(INPUT, sheet_name="Sheet1")
    a, t = clean(d["A"]), clean(d["T"])
    pet = pd.to_numeric(d["taupet"], errors="coerce").eq(1) & d["DIAG"].ne("健康对照")  # v3.60: patients only
    masks = {"A−T−": a.eq("A-") & t.eq("T-") & pet, "A+T−": a.eq("A+") & t.eq("T-") & pet,
             "A+T+": a.eq("A+") & t.eq("T+") & pet}
    counts = {k: int(m.sum()) for k, m in masks.items()}
    if counts != {"A−T−": 15, "A+T−": 18, "A+T+": 34}:
        raise AssertionError(f"Unexpected tau PET-defined group sizes: {counts}")

    vals = {name: {g: pd.to_numeric(d.loc[m, col], errors="coerce").dropna().to_numpy() for g, m in masks.items()}
            for name, col in MARKERS}
    rows = []
    for name, _ in MARKERS:
        v = vals[name]
        r1, p1 = rb(v["A−T−"], v["A+T−"])
        r2, p2 = rb(v["A+T−"], v["A+T+"])
        rows.append({"Marker": name,
                     **{f"{g} median": float(np.median(v[g])) for g in masks},
                     "Amyloid step r (A+T− vs A−T−)": r1, "Amyloid step p": p1,
                     "Tau step r (A+T+ vs A+T−)": r2, "Tau step p": p2})
    res = pd.DataFrame(rows)
    res["Amyloid step q"] = bh(res["Amyloid step p"])
    res["Tau step q"] = bh(res["Tau step p"])

    # A-status AUC within tau PET-defined T− (orientation fixed: higher = A+, except Aβ1-42 lower = A+)
    tm = masks["A−T−"] | masks["A+T−"]
    y = a[tm].eq("A+").astype(int).to_numpy()
    rng = np.random.default_rng(SEED)
    idx = [rng.integers(0, len(y), len(y)) for _ in range(N_BOOT)]
    idx = [i for i in idx if y[i].min() != y[i].max()]
    scores = {}
    for name, col in MARKERS:
        x = pd.to_numeric(d.loc[tm, col], errors="coerce").to_numpy()
        x = -x if name == "Aβ1-42" else x
        scores[name] = x
        auc = roc_auc_score(y, x)
        bs = [roc_auc_score(y[i], x[i]) for i in idx]
        res.loc[res.Marker == name, "T− A-status AUC"] = auc
        res.loc[res.Marker == name, "T− A-status AUC CI lower"] = np.percentile(bs, 2.5)
        res.loc[res.Marker == name, "T− A-status AUC CI upper"] = np.percentile(bs, 97.5)
    dd = [roc_auc_score(y[i], scores["eMTBR-tau243"][i]) - roc_auc_score(y[i], scores["p-tau217"][i]) for i in idx]
    d_obs = roc_auc_score(y, scores["eMTBR-tau243"]) - roc_auc_score(y, scores["p-tau217"])
    d_lo, d_hi = np.percentile(dd, [2.5, 97.5])
    res.to_csv(OUT / "AT_two_step_comparison.csv", index=False)

    md = ["# Two-step A/T comparison, tau PET-defined groups",
          f"Groups: A−T− n={counts['A−T−']}, A+T− n={counts['A+T−']}, A+T+ n={counts['A+T+']}. "
          "Rank-biserial r: positive = higher in the later group. BH-FDR across 7 markers within each step. "
          f"T− A-status AUC: A+T− vs A−T−, {len(idx)} participant bootstrap samples (seed {SEED}).", "",
          "| Marker | A−T− | A+T− | A+T+ | Amyloid step r (q) | Tau step r (q) | T− A-status AUC (95% CI) |",
          "|---|---|---|---|---|---|---|"]
    for _, r in res.iterrows():
        md.append(f"| {r.Marker} | {r['A−T− median']:.3g} | {r['A+T− median']:.3g} | {r['A+T+ median']:.3g} | "
                  f"{r['Amyloid step r (A+T− vs A−T−)']:.2f} ({qtext(r['Amyloid step q'])}) | "
                  f"{r['Tau step r (A+T+ vs A+T−)']:.2f} ({qtext(r['Tau step q'])}) | "
                  f"{r['T− A-status AUC']:.3f} ({r['T− A-status AUC CI lower']:.3f}–{r['T− A-status AUC CI upper']:.3f}) |")
    md += ["", f"Paired-bootstrap ΔAUC (eMTBR-tau243 − p-tau217) for A status within tau PET-defined T−: "
               f"{d_obs:.3f} (95% CI {d_lo:.3f} to {d_hi:.3f})."]
    (OUT / "AT_two_step_comparison.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))

    # ---------------- figure ----------------
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
    fig = plt.figure(figsize=(7.2, 4.3))
    gs = fig.add_gridspec(2, 6, width_ratios=[1, 1, 1, 1, 0.35, 3.3], left=0.07, right=0.985,
                          top=0.9, bottom=0.12, wspace=0.85, hspace=0.55)
    cells = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)]
    for (name, _), (rr, cc) in zip(MARKERS, cells):
        ax = fig.add_subplot(gs[rr, cc])
        row = res[res.Marker == name].iloc[0]
        data = [np.log10(vals[name][g]) for g, _ in GROUPS]
        bp = ax.boxplot(data, positions=[0, 1, 2], widths=0.55, showfliers=False, patch_artist=True,
                        medianprops={"color": INK, "linewidth": 1.1},
                        whiskerprops={"color": MUTED, "linewidth": 0.8}, capprops={"color": MUTED, "linewidth": 0.8},
                        boxprops={"linewidth": 0.8})
        for patch, (_, col) in zip(bp["boxes"], GROUPS):
            patch.set_facecolor(col); patch.set_alpha(0.28); patch.set_edgecolor(col)
        jr = np.random.default_rng(7)
        for k, ((g, col), arr) in enumerate(zip(GROUPS, data)):
            ax.scatter(k + jr.uniform(-0.13, 0.13, len(arr)), arr, s=9, color=col, edgecolor="#333333",
                       linewidth=0.3, zorder=3)
        lo, hi = min(map(np.min, data)), max(map(np.max, data))
        span = hi - lo
        ax.set_ylim(lo - 0.05 * span, hi + 0.42 * span)
        for (x0, x1, level, r, q) in [(0, 1, hi + 0.10 * span, row["Amyloid step r (A+T− vs A−T−)"], row["Amyloid step q"]),
                                      (1, 2, hi + 0.26 * span, row["Tau step r (A+T+ vs A+T−)"], row["Tau step q"])]:
            sig = q < 0.05
            ax.plot([x0, x0, x1, x1], [level - 0.03 * span, level, level, level - 0.03 * span],
                    color=INK if sig else MUTED, linewidth=0.7)
            ax.text((x0 + x1) / 2, level + 0.01 * span, f"r={r:.2f}{'*' if sig else ''}", ha="center", va="bottom",
                    fontsize=5.8, color=INK if sig else MUTED, fontweight="bold" if sig else "normal")
        ticks = nice_log_ticks(lo - 0.05 * span, hi)
        ax.set_yticks(np.log10(ticks))
        ax.set_yticklabels([tick_label(v) for v in ticks], fontsize=5.6)
        ax.set_ylim(lo - 0.05 * span, hi + 0.42 * span)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels([g[:2] + "\n" + g[2:] for g, _ in GROUPS], fontsize=5.6, linespacing=0.95)
        ax.set_title(name, fontsize=7, color=INK, pad=3)
        ax.set_ylabel("pg/mL (log)", fontsize=5.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    key = fig.add_subplot(gs[1, 3]); key.axis("off")
    key.text(0, 0.95, "Tau PET-defined groups", fontsize=6.2, color=INK, va="top", fontweight="bold")
    for k, (g, col) in enumerate(GROUPS):
        key.scatter([0.05], [0.72 - k * 0.17], s=22, color=col, edgecolor="#333333", linewidth=0.3)
        key.text(0.16, 0.72 - k * 0.17, f"{g} (n={counts[g]})", fontsize=6, va="center", color=INK)
    key.set_xlim(0, 1); key.set_ylim(0, 1)

    axb = fig.add_subplot(gs[:, 5])
    xa = res["Amyloid step r (A+T− vs A−T−)"].abs(); yt = res["Tau step r (A+T+ vs A+T−)"].abs()
    axb.plot([0, 1], [0, 1], color=GRID, linewidth=0.8, zorder=0)
    for (_, r), x, yv in zip(res.iterrows(), xa, yt):
        hl = r.Marker == "eMTBR-tau243"
        axb.scatter(x, yv, s=46 if hl else 30, color="#C23B53" if hl else "#555555", edgecolor="white",
                    linewidth=0.8, zorder=3)
        off = {"eMTBR-tau243": (0.035, 0.0, "left"), "p-tau217": (0.035, -0.02, "left"), "BD-tau": (0.035, 0.0, "left"),
               "GFAP": (0.0, 0.055, "center"), "Aβ1-42": (0.0, 0.055, "center"), "p-tau231": (0.035, 0.0, "left"),
               "NfL": (0.035, 0.0, "left")}[r.Marker]
        axb.text(x + off[0], yv + off[1], r.Marker, fontsize=6.4, ha=off[2], va="center",
                 color=INK, fontweight="bold" if hl else "normal")
    axb.set_xlim(0, 1); axb.set_ylim(0, 1)
    axb.set_aspect("equal", adjustable="box")
    axb.set_xlabel("Amyloid step |r|\n(A−T− → A+T−)", fontsize=6.6)
    axb.set_ylabel("Tau step |r|  (A+T− → A+T+)", fontsize=6.6)
    axb.grid(color=GRID, linewidth=0.5); axb.set_axisbelow(True)
    for s in ("top", "right"):
        axb.spines[s].set_visible(False)
    fig.text(0.01, 0.955, "a", fontsize=10, fontweight="bold")
    fig.text(0.575, 0.955, "b", fontsize=10, fontweight="bold")
    for ext, kw in [("png", {"dpi": 300}), ("pdf", {}), ("svg", {}),
                    ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(OUT / f"{STEM}.{ext}", facecolor="white", **kw)
    print(f"wrote {STEM}.*")


if __name__ == "__main__":
    main()
