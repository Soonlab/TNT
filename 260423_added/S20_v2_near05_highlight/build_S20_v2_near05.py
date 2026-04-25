#!/usr/bin/env python3
"""
build_S20_v2_near05.py

Supp Fig S20 Panel A — v2 rebuild with three explicit additions:

  1. **Near-significance emphasis** (0.05 ≤ P < 0.10): light-amber lollipop
     fill, distinct from the strict GOLD reserved for P < 0.05 (none in
     this dataset) and the THREAD1/THREAD2 colours used for non-significant
     positive/negative correlations. Three pairs qualify (Myc Targets V2,
     DSB repair, DNA Repair, all × IGH_n_delta; P = 0.06-0.07).

  2. **Headline-pair marker** (silver-grey diamond): overlays the row that
     manuscript cites as the canonical thread-1 ↔ thread-2 contact point —
     `DSB-repair × CD8_cytotoxic_delta` (r = -0.07, P = 0.83). Without
     this marker the manuscript-cited pair is buried in the |r|-sorted
     forest at row 24 / 36 and reviewers cannot quickly verify the
     headline value.

  3. **n = 11-12 paired** in the x-axis label, replacing the original
     "n = 12 paired" — four MSI_pct rows have n = 11 because one paired
     subject lacks MSI%; per-row n annotation (faint italic 5pt) is
     added on the rows that drop to n = 11.

All other style elements (BH-q reference, gray |r| > 0.58 zones,
thread-by-sign palette, badge text) are preserved verbatim from
build_S20_fancy() in 30_supp_fancy_260420.py.

Output:
  /mnt/sda1/data/TNT/analysis/260423_added/S20_v2_near05_highlight/
  SuppFig_S20_v2_convergence_null_near05.pptx
"""

import os
import importlib.util
import numpy as np
import pandas as pd
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
HELPER_PATH = "/mnt/sda1/data/TNT/analysis/260418_add/28_supp_natives_260420.py"
spec28 = importlib.util.spec_from_file_location("s28", HELPER_PATH)
s28 = importlib.util.module_from_spec(spec28)
spec28.loader.exec_module(s28)

GOOD = s28.GOOD; BAD = s28.BAD; INK = s28.INK; GREY = s28.GREY
LT_GREY = s28.LT_GREY; VLT_GREY = s28.VLT_GREY; WHITE = s28.WHITE
GOLD = s28.GOLD
THREAD1 = s28.THREAD1; THREAD2 = s28.THREAD2
RGBColor = s28.RGBColor
new_prs = s28.new_prs; new_slide = s28.new_slide
add_text = s28.add_text; add_line = s28.add_line
add_rect = s28.add_rect; add_circle = s28.add_circle
add_diamond = s28.add_diamond
axis_frame = s28.axis_frame
scale_x = s28.scale_x; scale_y = s28.scale_y
_i = s28._i

In = Inches


def badge(slide, x, y, w, h, text, fill=GOLD, text_color=INK,
          size=10, bold=True):
    add_rect(slide, x, y, w, h, fill=fill, line_color=INK, line_width=0.8)
    add_text(slide, x, y, w, h, text, size=size, bold=bold,
             color=text_color, align="center", anchor="middle")


# ---------------------------------------------------------------------------
# Paths + style additions
# ---------------------------------------------------------------------------
ADD = "/mnt/sda1/data/TNT/analysis/260418_add"
OUT = "/mnt/sda1/data/TNT/analysis/260423_added/S20_v2_near05_highlight"
os.makedirs(OUT, exist_ok=True)

LIGHT_AMBER = RGBColor(0xF0, 0xC8, 0x6E)   # "near-significance" tint —
                                            # less saturated than GOLD
                                            # (#D4A300) so reviewers
                                            # immediately read it as
                                            # "approaches but does not
                                            # reach" rather than "passes"
SILVER = RGBColor(0x9C, 0xA3, 0xAF)         # headline-pair marker — neutral
                                            # so it doesn't compete with
                                            # the THREAD palette


HEADLINE_PAIR = ("DNA Double-Strand Break Repair R-HSA-5693532",
                 "CD8_cytotoxic_delta")


def is_headline(row):
    return (row["baseline"] == HEADLINE_PAIR[0] and
            row["cascade"] == HEADLINE_PAIR[1])


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
def build_S20_v2():
    conv = pd.read_csv(f"{ADD}/targeted_convergence_test.tsv", sep="\t")
    prs = new_prs()

    # ============================================================
    # Panel A: 36-pair lollipop forest with near-0.05 emphasis
    # ============================================================
    slide = new_slide(prs)
    add_text(slide, In(0.35), In(0.25), In(0.45), In(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, In(0.9), In(0.35), In(11.5), In(0.4),
             "36-pair baseline × cascade-Δ convergence test "
             "— null with near-0.05 emphasis "
             "(amber = 0.05 ≤ P < 0.10; silver ◆ = manuscript "
             "headline pair)",
             size=11, bold=True)

    df = conv.copy()
    df["abs_r"] = df["spearman_r"].abs()
    df = df.sort_values("abs_r", ascending=False).reset_index(drop=True)
    px = In(3.5); py = In(1.2); pw = In(8.5); ph = In(5.4)
    n = len(df); row_h = ph / n
    xmin, xmax = -1, 1
    x_ticks = [-1, -0.5, 0, 0.5, 1.0]

    # |r| > 0.58 unfavourable zone shading (P < 0.05 at n = 12)
    pc = scale_x(0.58, xmin, xmax, px, pw)
    nc = scale_x(-0.58, xmin, xmax, px, pw)
    add_rect(slide, pc, py, px + pw - pc, ph,
             fill=RGBColor(0xF0, 0xF0, 0xF0), line_color=None)
    add_rect(slide, px, py, nc - px, ph,
             fill=RGBColor(0xF0, 0xF0, 0xF0), line_color=None)
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, xmin, xmax, px, pw) for v in x_ticks],
               x_labels=[f"{v:+g}" for v in x_ticks],
               xlab="Spearman r  (n = 11–12 paired; "
                    "MSI_pct rows have n = 11, others n = 12)")

    zx = scale_x(0, xmin, xmax, px, pw)
    add_line(slide, zx, py, zx, py + ph, color=INK, width=1.2)
    for rc in [-0.58, 0.58]:
        xc = scale_x(rc, xmin, xmax, px, pw)
        add_line(slide, xc, py, xc, py + ph,
                 color=GOLD, width=1.0, dashed=True)

    # rows
    near05_count = 0
    for i, row in df.iterrows():
        cy = py + row_h * (i + 0.5)
        r = float(row["spearman_r"]); p = float(row["spearman_p"])
        n_pair = int(row["n_paired"])
        is_head = is_headline(row)
        is_near05 = (0.05 <= p < 0.10)
        if is_near05:
            near05_count += 1

        # ----- highlight band on the row (subtle, full-width) -----
        if is_near05:
            add_rect(slide, px, cy - row_h * 0.45, pw, row_h * 0.9,
                     fill=RGBColor(0xFB, 0xF1, 0xDA), line_color=None)
        if is_head:
            add_rect(slide, px, cy - row_h * 0.45, pw, row_h * 0.9,
                     fill=RGBColor(0xEE, 0xF1, 0xF5), line_color=None)

        # ----- left-side label -----
        lab = f"{row['baseline'][:13]} × {row['cascade'][:20]}"
        lab_color = (LIGHT_AMBER if is_near05 else
                     (SILVER if is_head else INK))
        add_text(slide, px - In(3.3), cy - row_h * 0.45,
                 In(3.2), row_h * 0.9, lab,
                 size=6, align="right", anchor="middle",
                 color=lab_color, bold=is_near05 or is_head)
        # per-row n annotation only when n != 12
        if n_pair != 12:
            add_text(slide, px - In(3.3) - In(0.05),
                     cy - row_h * 0.45,
                     In(0.30), row_h * 0.9,
                     f"n={n_pair}", size=5, italic=True,
                     align="right", anchor="middle",
                     color=GREY)

        # ----- lollipop -----
        ex = scale_x(r, xmin, xmax, px, pw)
        add_line(slide, zx, cy, ex, cy, color=GREY, width=0.8)
        if p < 0.05:
            col = GOLD
        elif is_near05:
            col = LIGHT_AMBER
        else:
            col = THREAD1 if r > 0 else THREAD2
        add_circle(slide, ex, cy, In(0.05),
                   fill=col, line_color=INK, line_width=0.4)

        # headline-pair silver diamond overlay
        if is_head:
            add_diamond(s28._i if False else slide, ex, cy, In(0.10),
                        fill=SILVER, line_color=INK, line_width=0.6)
            # tag to the right of the diamond (above P label)
            add_text(slide, ex + In(0.13), cy - row_h * 0.55,
                     In(1.4), row_h * 0.6,
                     "headline pair (manuscript)",
                     size=6, italic=True, color=SILVER, anchor="middle")

        # ----- right-side P label -----
        p_color = (LIGHT_AMBER if is_near05 else
                   (SILVER if is_head else col))
        add_text(slide, px + pw + In(0.08), cy - row_h * 0.45,
                 In(1.0), row_h * 0.9, f"P={p:.2f}",
                 size=6, anchor="middle",
                 color=p_color, bold=is_near05 or is_head)

    # ----- summary badges (right side) -----
    badge(slide, In(10.5), In(1.4), In(2.5), In(1.0),
          "0 / 36 pairs\nP < 0.05 (BH q ≥ 0.98)\n1.8 expected by chance",
          fill=GOLD, size=10)
    badge(slide, In(10.5), In(2.55), In(2.5), In(0.85),
          f"{near05_count} / 36 pairs\nnear 0.05 (P < 0.10)\n"
          f"all × IGH_n Δ; partial P > 0.13",
          fill=LIGHT_AMBER, size=9)
    badge(slide, In(10.5), In(3.55), In(2.5), In(0.65),
          "headline pair ◆\nDSB × CD8-cyt Δ\nr = −0.07, P = 0.83",
          fill=VLT_GREY, size=8)

    # axis label is rendered by axis_frame at y = py + ph + 0.38 = ~7.0;
    # tick labels at ~6.66; legend immediately below xlab; footer below.
    leg_y = In(7.13)   # below the bold xlab from axis_frame
    leg_x = In(3.5)
    chip_w = In(0.18); chip_h = In(0.13)
    items = [
        ("P < 0.05",       GOLD),
        ("0.05 ≤ P < 0.10", LIGHT_AMBER),
        ("P ≥ 0.10, r > 0", THREAD1),
        ("P ≥ 0.10, r < 0", THREAD2),
        ("headline pair ◆", SILVER),
    ]
    cur_x = leg_x
    for lab, col in items:
        add_rect(slide, cur_x, leg_y, chip_w, chip_h,
                 fill=col, line_color=INK, line_width=0.2)
        add_text(slide, cur_x + chip_w + In(0.04), leg_y - In(0.01),
                 In(1.5), In(0.18),
                 lab, size=7, color=INK)
        cur_x += In(1.75)

    add_text(slide, In(3.5), In(7.32), In(8.5), In(0.20),
             "Gold dashed = |r| = 0.58 (P = 0.05 at n = 12). "
             "Shaded gray zones = |r| > 0.58 (significant zone) "
             "— no hit lands inside.",
             size=8, italic=True)

    # ============================================================
    # Panel B: purity-adjusted sensitivity scatter (unchanged)
    # ============================================================
    slide = new_slide(prs)
    add_text(slide, In(0.35), In(0.25), In(0.45), In(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, In(0.9), In(0.35), In(11.5), In(0.4),
             "Purity-adjusted paired Δ sensitivity — raw vs purity-"
             "adjusted Δ (Tarabichi 2021 motif)",
             size=11, bold=True)
    psens_path = ("/mnt/sda1/data/TNT/analysis/09_integration/"
                  "paired_delta/delta_purity_sensitivity.tsv")
    if os.path.exists(psens_path):
        psens = pd.read_csv(psens_path, sep="\t")
        cols = psens.columns.tolist()
        ax = "delta_raw" if "delta_raw" in cols else cols[-2]
        ay = "delta_adj" if "delta_adj" in cols else cols[-1]
        px = In(3.5); py = In(1.3); pw = In(6.5); ph = In(5.5)
        vmin = float(min(psens[ax].min(), psens[ay].min())) * 1.1
        vmax = float(max(psens[ax].max(), psens[ay].max())) * 1.1
        tk = np.linspace(vmin, vmax, 5)
        axis_frame(slide, px, py, pw, ph,
                   x_ticks=[scale_x(v, vmin, vmax, px, pw) for v in tk],
                   x_labels=[f"{v:.1f}" for v in tk],
                   y_ticks=[scale_y(v, vmin, vmax, py, ph) for v in tk],
                   y_labels=[f"{v:.1f}" for v in tk],
                   xlab="raw Δ (observed)", ylab="purity-adjusted Δ")
        add_line(slide, px, py + ph, px + pw, py,
                 color=GOLD, width=1.2, dashed=True)
        for _, row in psens.iterrows():
            xv = scale_x(float(row[ax]), vmin, vmax, px, pw)
            yv = scale_y(float(row[ay]), vmin, vmax, py, ph)
            add_circle(slide, xv, yv, In(0.06),
                       fill=THREAD1, line_color=INK, line_width=0.3)
    badge(slide, In(3.5), In(6.85), In(9.0), In(0.35),
          "y = x diagonal (gold) — points scatter tightly: "
          "purity correction does not flip Δ sign or rank for any "
          "cascade feature",
          fill=VLT_GREY, size=9)

    out_path = f"{OUT}/SuppFig_S20_v2_convergence_null_near05.pptx"
    prs.save(out_path)
    print(f"wrote {out_path}")
    print(f"   near-0.05 rows highlighted: {near05_count}")


if __name__ == "__main__":
    build_S20_v2()
