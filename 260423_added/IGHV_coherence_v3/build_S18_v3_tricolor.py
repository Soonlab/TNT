#!/usr/bin/env python3
"""
build_S18_v3_tricolor.py

Supp Fig S18 (IGHV directional coherence) — v3 rebuild with **three-class
colour-coded highlights** based on which group is statistically coherent:

    1. IGHV6-1     — GOLD  (#D4A300)
       between-group polar-opposite (Fisher P=0.061) AND within-good
       sign test P=0.031 (6/6 unanimous down). Strongest signal.

    2. IGHV3-38-3  — GOOD  (#0A7D6E, deep teal)
       within-good sign test P=0.063 (5/5 down) — good responders only
       are coherent; bad responders are mixed (3/3).

    3. IGHV3-66    — BAD   (#C53E1F, deep coral)
       within-bad sign test P=0.031 (6/6 down) — bad responders only
       are coherent; good responders are mixed.

This replaces the original three-gene focus set
{IGHV6-1, IGHV3-7, IGHV3-74} which highlighted IGHV3-7 / IGHV3-74 as
illustrative quadrant examples without individual statistical support
(both Fisher P=1.000).

Style is unchanged from 30_supp_fancy_260420.py build_S18_fancy() —
this script imports the same helpers from 28_supp_natives_260420.py to
guarantee identical look-and-feel.

Output: /mnt/sda1/data/TNT/analysis/260423_added/IGHV_coherence_v3/
        SuppFig_S18_v3_IGHV_coherence_tricolor.pptx
"""

import os
import importlib.util
import numpy as np
import pandas as pd
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE


# ---------------------------------------------------------------------------
# Shared helpers from 28_supp_natives_260420
# ---------------------------------------------------------------------------
HELPER_PATH = "/mnt/sda1/data/TNT/analysis/260418_add/28_supp_natives_260420.py"
spec28 = importlib.util.spec_from_file_location("s28", HELPER_PATH)
s28 = importlib.util.module_from_spec(spec28)
spec28.loader.exec_module(s28)

GOOD = s28.GOOD; BAD = s28.BAD; INK = s28.INK; GREY = s28.GREY
LT_GREY = s28.LT_GREY; VLT_GREY = s28.VLT_GREY; WHITE = s28.WHITE
GOLD = s28.GOLD
RGBColor = s28.RGBColor
new_prs = s28.new_prs; new_slide = s28.new_slide
add_text = s28.add_text; add_line = s28.add_line
add_rect = s28.add_rect; add_circle = s28.add_circle
add_diamond = s28.add_diamond
axis_frame = s28.axis_frame
scale_x = s28.scale_x; scale_y = s28.scale_y
_i = s28._i

In = Inches


def add_star(slide, cx, cy, r, fill=GOLD, line_color=None):
    r = max(_i(r), 1)
    shp = slide.shapes.add_shape(MSO_SHAPE.STAR_5_POINT,
                                 _i(cx) - r, _i(cy) - r, 2 * r, 2 * r)
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    if line_color is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line_color
        shp.line.width = Pt(0.5)
    shp.shadow.inherit = False
    s28.kill_shadow(shp)
    return shp


def badge(slide, x, y, w, h, text, fill=GOLD, text_color=INK,
          size=10, bold=True):
    add_rect(slide, x, y, w, h, fill=fill, line_color=INK, line_width=0.8)
    add_text(slide, x, y, w, h, text,
             size=size, bold=bold, color=text_color,
             align="center", anchor="middle")


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ADD = "/mnt/sda1/data/TNT/analysis/260418_add"
OUT = "/mnt/sda1/data/TNT/analysis/260423_added/IGHV_coherence_v3"
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------------------
# Highlight scheme
# ---------------------------------------------------------------------------
HIGHLIGHT = {
    "IGHV6-1":    {"color": GOLD,
                   "category": "polar-opposite",
                   "tag":      "★ polar-opposite\n(Fisher P=0.061;\ngood 6/6 down,\nbad 4/2 up)"},
    "IGHV3-38-3": {"color": GOOD,
                   "category": "good-coherent",
                   "tag":      "● good-coherent\n(sign P=0.063;\ngood 5/5 down,\nbad mixed 3/3)"},
    "IGHV3-66":   {"color": BAD,
                   "category": "bad-coherent",
                   "tag":      "● bad-coherent\n(sign P=0.031;\nbad 6/6 down,\ngood mixed)"},
}


def cat_color(v_gene):
    h = HIGHLIGHT.get(v_gene)
    return h["color"] if h else None


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
def build_S18_v3():
    ig = pd.read_csv(f"{ADD}/trust4_ighv_directional_stats.tsv", sep="\t")
    prs = new_prs()

    # ============================================================
    # Panel A: 53 V-gene forest with tri-colour highlight scheme
    # ============================================================
    slide = new_slide(prs)
    add_text(slide, In(0.35), In(0.25), In(0.45), In(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, In(0.9), In(0.35), In(11.5), In(0.4),
             "IGH V-gene directional-coherence forest "
             "— tri-colour highlight by group-level coherence "
             "(GOLD=polar-opposite, TEAL=good-only, CORAL=bad-only)",
             size=11, bold=True)

    df = (ig.copy()
          .sort_values("coherence_gap", ascending=False)
          .reset_index(drop=True))
    px = In(3.0); py = In(1.2); pw = In(8.5); ph = In(5.5)
    n_g = len(df)
    row_h = ph / n_g
    zx = px + pw / 2

    # gray zero-band |majority| < 0.55 (≈ 3/6)
    band_half = pw / 2 * ((3.3 - 3.0) / 6.0)
    add_rect(slide, zx - band_half, py, band_half * 2, ph,
             fill=RGBColor(0xEE, 0xEE, 0xEE), line_color=None)
    add_line(slide, zx, py, zx, py + ph, color=INK, width=1.2)
    add_text(slide, px, py - In(0.25), pw / 2, In(0.22),
             "bad (n up / 6) ←", size=9, align="right", color=BAD)
    add_text(slide, px + pw / 2, py - In(0.25), pw / 2, In(0.22),
             "→ good (n down / 6)", size=9, align="left", color=GOOD)

    # --- explicit x-axis at the bottom ---
    axis_y = py + ph
    add_line(slide, px, axis_y, px + pw, axis_y, color=INK, width=0.8)
    tick_positions = [
        (px,                      "-6",  BAD),
        (px + pw / 4,             "-3",  BAD),
        (zx,                      "0",   INK),
        (px + 3 * pw / 4,         "+3",  GOOD),
        (px + pw,                 "+6",  GOOD),
    ]
    for tx, lab, col in tick_positions:
        add_line(slide, tx, axis_y, tx, axis_y + In(0.06),
                 color=INK, width=0.6)
        add_text(slide, tx - In(0.18), axis_y + In(0.07),
                 In(0.36), In(0.16),
                 lab, size=8, color=col, align="center", bold=True)
    add_text(slide, px, axis_y + In(0.26), pw, In(0.20),
             "subject count (bad ↑ on left, good ↓ on right; max ±6)",
             size=8, italic=True, color=INK, align="center")

    for i, row in df.iterrows():
        cy = py + row_h * (i + 0.5)
        v = row["v_gene"]
        hi_color = cat_color(v)
        is_focus = hi_color is not None
        g_down = int(row["good_n_down"])
        b_up = int(row["bad_n_up"])
        right_len = pw / 2 * (g_down / 6.0)
        left_len = pw / 2 * (b_up / 6.0)

        # --- Bar coloring ---
        # Right (good_n_down) bar: keep GOOD/teal as group identity,
        # but if this row is highlighted, override to its category colour
        right_fill = hi_color if is_focus else GOOD
        # Left (bad_n_up) bar: keep BAD as group identity. Override
        # only for the bad-coherent highlight (IGHV3-66) so its identity
        # dominates visually.
        left_fill = hi_color if (is_focus and v == "IGHV3-66") else BAD

        add_rect(slide, zx, cy - row_h * 0.3, right_len, row_h * 0.6,
                 fill=right_fill, line_color=INK, line_width=0.2)
        add_rect(slide, zx - left_len, cy - row_h * 0.3,
                 left_len, row_h * 0.6, fill=left_fill,
                 line_color=INK, line_width=0.2)

        # Highlight ring: thin coloured frame around the entire row
        if is_focus:
            add_rect(slide, px, cy - row_h * 0.45, pw, row_h * 0.9,
                     fill=None, line_color=hi_color, line_width=1.2)

        # --- Label (left of plot) ---
        if is_focus:
            sym = {"polar-opposite": "★ ",
                   "good-coherent":  "● ",
                   "bad-coherent":   "● "}[HIGHLIGHT[v]["category"]]
        else:
            sym = "  "
        add_text(slide, px - In(1.3), cy - row_h * 0.4,
                 In(1.25), row_h * 0.8,
                 sym + v, size=6,
                 align="right", bold=is_focus,
                 color=(hi_color if is_focus else INK))

        # --- Per-row Fisher P (right of plot) ---
        fp = row.get("fisher_P_updown", 1)
        p_col = hi_color if is_focus else (GOLD if fp <= 0.10 else GREY)
        add_text(slide, px + pw + In(0.05), cy - row_h * 0.4,
                 In(0.9), row_h * 0.8, f"P={fp:.2f}",
                 size=6, color=p_col, anchor="middle",
                 bold=is_focus)

    # --- callout badges to the right, one per highlighted V-gene ---
    legend_x = In(11.9)
    for v, info in HIGHLIGHT.items():
        if v not in df["v_gene"].values:
            continue
        idx = df.index[df["v_gene"] == v].tolist()[0]
        cy_v = py + row_h * (idx + 0.5)
        # connector line from row to badge
        add_line(slide, px + pw + In(1.0), cy_v,
                 legend_x - In(0.05), cy_v,
                 color=info["color"], width=0.8, dashed=True)
        badge(slide, legend_x, cy_v - In(0.32),
              In(1.4), In(0.65), info["tag"],
              fill=info["color"],
              text_color=WHITE if info["color"] in (GOOD, BAD, GOLD)
                                    else INK,
              size=6, bold=True)

    # footnote
    add_text(slide, In(3.0), In(7.20), In(8.5), In(0.45),
             "Sorted by coherence_gap descending. Gray band = "
             "|majority − 0.5| < 3.3/6 (indistinguishable from chance). "
             "Tri-colour highlights mark the three V-genes with "
             "individual-level statistical support: IGHV6-1 (between-"
             "group Fisher P=0.061 + within-good sign P=0.031, polar-"
             "opposite), IGHV3-38-3 (within-good sign P=0.063, good-"
             "coherent only), IGHV3-66 (within-bad sign P=0.031, bad-"
             "coherent only). Repertoire-level aggregate Wilcoxon "
             "P = 0.035.",
             size=8, italic=True)

    # legend block (left side, before the plot, to avoid colliding with
    # right-hand per-row callouts)
    leg_x = In(0.45); leg_y = In(0.85)
    add_text(slide, leg_x, leg_y, In(2.4), In(0.20),
             "Highlight category",
             size=10, bold=True, color=INK)
    for k, (lab, col) in enumerate([
            ("polar-opposite", GOLD),
            ("good-coherent only", GOOD),
            ("bad-coherent only", BAD)]):
        ry = leg_y + In(0.26 + 0.20 * k)
        add_rect(slide, leg_x, ry, In(0.22), In(0.16), fill=col,
                 line_color=INK, line_width=0.2)
        add_text(slide, leg_x + In(0.28), ry - In(0.01),
                 In(2.2), In(0.18),
                 lab, size=8, color=INK)

    # ============================================================
    # Panel B: pattern-class scatter, same tri-colour highlight
    # ============================================================
    slide = new_slide(prs)
    add_text(slide, In(0.35), In(0.25), In(0.45), In(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, In(0.9), In(0.35), In(11.5), In(0.4),
             "Pattern-class scatter (good vs bad majority fraction) "
             "— tri-colour highlight on three statistically supported "
             "V-genes",
             size=11, bold=True)
    px = In(2.5); py = In(1.2); pw = In(5.5); ph = In(5.5)
    zx_ = scale_x(0.75, 0.5, 1.0, px, pw)
    zy_ = scale_y(0.75, 0.5, 1.0, py, ph)
    # quadrant background shading
    add_rect(slide, zx_, py, px + pw - zx_, zy_ - py,
             fill=RGBColor(0xEE, 0xF7, 0xF4), line_color=None)
    add_rect(slide, zx_, zy_, px + pw - zx_, py + ph - zy_,
             fill=RGBColor(0xE5, 0xEF, 0xE8), line_color=None)
    add_rect(slide, px, py, zx_ - px, zy_ - py,
             fill=RGBColor(0xFA, 0xF0, 0xEC), line_color=None)
    add_rect(slide, px, zy_, zx_ - px, py + ph - zy_,
             fill=RGBColor(0xF5, 0xF5, 0xF5), line_color=None)
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, 0.5, 1.0, px, pw)
                        for v in [0.5, 0.625, 0.75, 0.875, 1.0]],
               x_labels=["0.5", "0.625", "0.75", "0.875", "1.0"],
               y_ticks=[scale_y(v, 0.5, 1.0, py, ph)
                        for v in [0.5, 0.625, 0.75, 0.875, 1.0]],
               y_labels=["0.5", "0.625", "0.75", "0.875", "1.0"],
               xlab="good majority fraction",
               ylab="bad majority fraction")
    add_line(slide, zx_, py, zx_, py + ph,
             color=GREY, width=0.6, dashed=True)
    add_line(slide, px, zy_, px + pw, zy_,
             color=GREY, width=0.6, dashed=True)

    # plot all V-genes
    for _, row in df.iterrows():
        v = row["v_gene"]
        gx = scale_x(float(row["good_majority_frac"]), 0.5, 1.0, px, pw)
        bx_ = scale_y(float(row["bad_majority_frac"]), 0.5, 1.0, py, ph)
        hi_color = cat_color(v)
        if hi_color is not None:
            add_star(slide, gx, bx_, In(0.11),
                     fill=hi_color, line_color=INK)
            add_text(slide, gx + In(0.10), bx_ - In(0.08),
                     In(1.2), In(0.18), v,
                     size=8, color=hi_color, bold=True)
        else:
            add_circle(slide, gx, bx_, In(0.04),
                       fill=INK, line_color=INK, line_width=0.3)

    # quadrant labels
    add_text(slide, zx_ + In(0.1), py + In(0.1), pw, In(0.22),
             "↗ both-coherent", size=8, color=GOOD, bold=True)
    add_text(slide, zx_ + In(0.1), zy_ + In(0.1), pw, In(0.22),
             "→ good-coherent only", size=8, color=GOOD)
    add_text(slide, px + In(0.1), py + In(0.1), pw, In(0.22),
             "↑ bad-coherent only", size=8, color=BAD)
    add_text(slide, px + In(0.1), zy_ + In(0.1), pw, In(0.22),
             "↙ stochastic", size=8, color=GREY)

    # right-hand legend
    lx = In(8.5); ly = In(1.3)
    add_text(slide, lx, ly, In(4.8), In(0.22),
             "Highlight V-genes (statistical support)",
             size=10, bold=True, color=INK)
    for k, (v, info) in enumerate(HIGHLIGHT.items()):
        ry = ly + In(0.32 + 0.30 * k)
        add_star(slide, lx + In(0.1), ry + In(0.1), In(0.09),
                 fill=info["color"], line_color=INK)
        add_text(slide, lx + In(0.28), ry - In(0.02),
                 In(4.6), In(0.16),
                 v, size=10, bold=True, color=info["color"])
        # description line
        add_text(slide, lx + In(0.28), ry + In(0.13),
                 In(4.6), In(0.16),
                 info["tag"].replace("\n", "  "),
                 size=7, color=INK)

    # save
    out_path = f"{OUT}/SuppFig_S18_v3_IGHV_coherence_tricolor.pptx"
    prs.save(out_path)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    build_S18_v3()
