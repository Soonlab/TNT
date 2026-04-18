#!/usr/bin/env python3
"""
19_fig6_native_pptx.py

Build Figure 6 (§3.9 Target engagement of four Thread-1 baseline factors
in paired radiation-phase biopsies) and its supplementary companion
Supp Fig S17 (member-pathway detail + per-subject Δ heatmap) as
native-editable PowerPoint decks.

Rules (same as script 18 / Fig 5):
    * one panel per slide
    * every in-plot text element is a add_textbox TEXT_BOX (editable,
      not grouped)
    * every line / tick / marker / bar / ribbon / cell is a native
      connector or auto-shape
    * Arial throughout, no panel titles inside the plot
    * TNT project palette: GOOD=#0a7d6e (deep teal) / BAD=#c53e1f
      (deep coral)
    * all shapes get empty <a:effectLst/> via kill_shadow to override
      PowerPoint theme's default soft shadow

Main Figure 6 panels
--------------------
  A  2x2 grid of paired pre->post composite-z spaghetti plots for the
     four Thread-1 factors (DSB/HDR repair, Tumor cell-cycle, E2F/MYC,
     EMT); 12 subjects per sub-plot (6 good + 6 bad); annotation box
     per sub-plot with within-group sign counts and binomial P.
  B  Directional-concordance grouped bar chart: for each of the four
     factors, the fraction of subjects moving in the biologically
     predicted direction, good vs bad, with n_pred/n_total labels and
     a chance-level (0.5) reference line. EMT good 6/6 (P=0.016) is
     highlighted because it is the only composite-level nominal
     single-factor significance.

Supp Figure S17 panels
----------------------
  A  Member-pathway sign-count dotted stack: for each of the 4 factors
     and each of its member pathways, show how many of the 6 good
     responders and how many of the 6 bad responders moved in the
     predicted direction.
  B  Per-subject Δ heatmap across all 17 signatures (4 composites + 13
     member pathways) × 12 paired subjects grouped good | bad,
     diverging colour scale with composite rows annotated.

Style references consulted (paired-biopsy pharmacodynamic convention)
---------------------------------------------------------------------
  - Tumeh et al Nature 2014 (PD-1 pre/post CD8 density): hollow-circle
    paired points connected by thin line, group-coloured.
  - Ribas lab (Tumeh / Riaz) Cell 2017 (pre/on-treatment ICI): small
    per-subject spaghetti arranged in sub-panel grids by marker.
  - Cercek et al NEJM 2022 (dostarlimab rectal): per-subject response
    waterfall + paired lesion trajectories.
  - Petitprez et al Nature 2020 / Helmink et al Nature 2020 / Cabrita
    et al Nature 2020 (B-cell TLS trilogy): directional fraction bars
    with n/N annotations.
  - Thorsson et al Immunity 2018 (TCGA immune landscape): grouped bar
    charts with chance-reference lines and small-N transparency.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree


# ---------------------------------------------------------------------------
# Shared infrastructure --- same as script 18
# ---------------------------------------------------------------------------
def kill_shadow(shape):
    elem = shape._element
    spPr = elem.find(qn('p:spPr'))
    if spPr is None:
        return
    for el in spPr.findall(qn('a:effectLst')):
        spPr.remove(el)
    etree.SubElement(spPr, qn('a:effectLst'))


GOOD = RGBColor(0x0A, 0x7D, 0x6E)
BAD = RGBColor(0xC5, 0x3E, 0x1F)
INK = RGBColor(0x22, 0x22, 0x22)
LINE = RGBColor(0x33, 0x33, 0x33)
GREY = RGBColor(0xBB, 0xBB, 0xBB)
LT_GREY = RGBColor(0xDD, 0xDD, 0xDD)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HIGHLIGHT = RGBColor(0xD4, 0xA3, 0x00)   # gold for EMT 6/6 callout
GOOD_HEX = (0x0A, 0x7D, 0x6E)
BAD_HEX = (0xC5, 0x3E, 0x1F)

FONT = "Arial"
SLIDE_W = Inches(6.5)
SLIDE_H = Inches(4.5)

OUT = "/data/data/TNT/analysis/260418_add/ppt"
DATA = "/data/data/TNT/analysis/260418_add"
os.makedirs(OUT, exist_ok=True)


def new_slide(prs, w=SLIDE_W, h=SLIDE_H):
    prs.slide_width = w
    prs.slide_height = h
    return prs.slides.add_slide(prs.slide_layouts[6])


def add_text(slide, x, y, w, h, text, size=8, bold=False,
             color=INK, align="left", anchor="middle", font=FONT):
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    tf.margin_bottom = 0
    tf.word_wrap = True
    tf.vertical_anchor = {"top": MSO_ANCHOR.TOP, "bottom": MSO_ANCHOR.BOTTOM,
                          "middle": MSO_ANCHOR.MIDDLE}[anchor]
    p = tf.paragraphs[0]
    p.alignment = {"left": PP_ALIGN.LEFT, "center": PP_ALIGN.CENTER,
                   "right": PP_ALIGN.RIGHT}[align]
    r = p.add_run()
    r.text = str(text)
    r.font.name = font
    r.font.size = Pt(size)
    r.font.bold = bool(bold)
    r.font.color.rgb = color
    kill_shadow(tb)
    return tb


def _i(v):
    """Coerce to int --- PowerPoint rejects XML width='0.0' type strings,
    so every coordinate / size passed to a shape creator must be an
    integer number of EMU."""
    return int(round(float(v)))


def add_line(slide, x1, y1, x2, y2, color=LINE, width=0.5, dashed=False):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,
                                   _i(x1), _i(y1), _i(x2), _i(y2))
    c.line.color.rgb = color
    c.line.width = Pt(width)
    if dashed:
        from pptx.enum.dml import MSO_LINE_DASH_STYLE
        try:
            c.line.dash_style = MSO_LINE_DASH_STYLE.DASH
        except Exception:
            pass
    kill_shadow(c)
    return c


def add_rect(slide, x, y, w, h, fill=None, line_color=None, line_width=0.5):
    # sizes must be strictly positive for PowerPoint to render
    w = max(_i(w), 1)
    h = max(_i(h), 1)
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, _i(x), _i(y), w, h)
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line_color is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line_color
        shp.line.width = Pt(line_width)
    shp.shadow.inherit = False
    kill_shadow(shp)
    return shp


def add_circle(slide, cx, cy, r, fill=None, line_color=None, line_width=0.5):
    r = max(_i(r), 1)
    shp = slide.shapes.add_shape(MSO_SHAPE.OVAL,
                                 _i(cx) - r, _i(cy) - r, 2 * r, 2 * r)
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line_color is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line_color
        shp.line.width = Pt(line_width)
    shp.shadow.inherit = False
    kill_shadow(shp)
    return shp


def add_diamond(slide, cx, cy, r, fill=None, line_color=None, line_width=0.5):
    r = max(_i(r), 1)
    shp = slide.shapes.add_shape(MSO_SHAPE.DIAMOND,
                                 _i(cx) - r, _i(cy) - r, 2 * r, 2 * r)
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line_color is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line_color
        shp.line.width = Pt(line_width)
    shp.shadow.inherit = False
    kill_shadow(shp)
    return shp


def lighten(rgb_hex, factor=0.72):
    """Blend an RGB triple towards white by 'factor' (0=identity, 1=white)."""
    r0, g0, b0 = rgb_hex
    return RGBColor(int(r0 + (255 - r0) * factor),
                    int(g0 + (255 - g0) * factor),
                    int(b0 + (255 - b0) * factor))


def draw_panel_letter(slide, letter):
    add_text(slide, Inches(0.15), Inches(0.1), Inches(0.4), Inches(0.35),
             letter, size=14, bold=True, color=INK, align="left")


# ---------------------------------------------------------------------------
# Load paired baseline-factor data
# ---------------------------------------------------------------------------
delta_df = pd.read_csv(f"{DATA}/baseline_factor_per_subject_delta.tsv",
                       sep="\t")
sign_df = pd.read_csv(f"{DATA}/baseline_factor_sign_table.tsv", sep="\t")
stats_df = pd.read_csv(f"{DATA}/baseline_factor_pharmacodynamics_stats.tsv",
                       sep="\t")

FACTORS = [
    ("DSB_HDR_repair", "DSB / HDR repair", "down"),
    ("Tumor_cellcycle", "Tumor cell-cycle", "down"),
    ("E2F_MYC_cellcycle", "E2F / MYC", "down"),
    ("EMT", "EMT", "up"),
]

# subset to composite-level rows for Panel A
comp = delta_df[delta_df.member == "composite"].copy()


# ===========================================================================
# MAIN FIGURE 6
# ===========================================================================
prs_main = Presentation()


def build_A():
    """2x2 grid of paired pre->post composite-z spaghetti plots."""
    s = new_slide(prs_main)
    draw_panel_letter(s, "A")

    # common y range: pool composite deltas and pre/post values to get
    # a comfortable symmetric range per factor
    layout = {
        # (factor key, row, col)
        "DSB_HDR_repair": (0, 0),
        "Tumor_cellcycle": (0, 1),
        "E2F_MYC_cellcycle": (1, 0),
        "EMT": (1, 1),
    }

    grid_ox = Inches(0.75)
    grid_oy = Inches(0.50)
    sub_w = Inches(2.55)
    sub_h = Inches(1.55)
    col_gap = Inches(0.45)
    row_gap = Inches(0.52)

    for fname, pretty, pred in FACTORS:
        row, col = layout[fname]
        px = grid_ox + col * (sub_w + col_gap)
        py = grid_oy + row * (sub_h + row_gap)

        sub = comp[comp.factor == fname].sort_values(
            ["response_bin", "subject_id"])

        # per-factor y range (symmetric)
        vals = np.concatenate([sub["pre"].values, sub["post"].values])
        y_span = max(3.0, float(np.ceil(np.max(np.abs(vals)))))
        y_lo, y_hi = -y_span, y_span

        # mini axes: x from pre(0) to post(1), y from -y_span to +y_span
        # use inner plot area margins
        ax_x = px + Inches(0.42)
        ax_y = py + Inches(0.08)
        ax_w = sub_w - Inches(0.52)
        ax_h = sub_h - Inches(0.45)

        def tx(v): return int(ax_x + (v) * ax_w)
        def ty(v):
            return int(ax_y + ax_h - (v - y_lo) / (y_hi - y_lo) * ax_h)

        # spines
        add_line(s, ax_x, ax_y, ax_x, ax_y + ax_h, LINE, 0.6)       # y
        add_line(s, ax_x, ax_y + ax_h, ax_x + ax_w,
                 ax_y + ax_h, LINE, 0.6)                            # x

        # y-axis ticks
        for yv in [y_lo, y_lo / 2, 0, y_hi / 2, y_hi]:
            yy = ty(yv)
            add_line(s, int(ax_x - Inches(0.04)), yy, int(ax_x), yy,
                     LINE, 0.4)
            add_text(s, int(ax_x - Inches(0.36)), yy - Inches(0.07),
                     Inches(0.32), Inches(0.14),
                     f"{yv:+.0f}" if abs(yv) >= 1 else f"{yv:.1f}",
                     size=6, color=INK, align="right")
        # x-axis tick labels
        for xv, xlab in [(0, "pre"), (1, "post")]:
            xx = tx(xv)
            add_line(s, xx, int(ax_y + ax_h), xx,
                     int(ax_y + ax_h + Inches(0.04)), LINE, 0.4)
            add_text(s, xx - Inches(0.2), int(ax_y + ax_h + Inches(0.05)),
                     Inches(0.4), Inches(0.14),
                     xlab, size=7, color=INK, align="center")

        # zero reference line (neutral for all factors)
        add_line(s, tx(0), ty(0), tx(1), ty(0), GREY, 0.4, dashed=True)

        # sub-plot factor name (below the plot, below tick labels)
        add_text(s, px, py + sub_h - Inches(0.22),
                 sub_w, Inches(0.18),
                 f"{pretty}   (pred: {pred})",
                 size=8, bold=True, color=INK, align="center")

        # =================================================================
        # ~~~ FANCY spaghetti composition (Nature/Cell paired-pre/post
        # convention with group summaries):
        #   (i)   very light individual spaghetti slopes beneath
        #   (ii)  thin IQR rectangles at pre and post per group
        #   (iii) hollow circle markers at each individual timepoint
        #   (iv)  bold group-median slope on top, filled diamond markers
        #   (v)   small raincloud-style jitter dots at pre and post
        # =================================================================

        good = sub[sub.response_bin == "good"]
        bad = sub[sub.response_bin == "bad"]

        def iqr(arr):
            q25, q75 = np.percentile(arr, [25, 75])
            return q25, q75

        g_pre_q = iqr(good.pre.values); g_post_q = iqr(good.post.values)
        b_pre_q = iqr(bad.pre.values); b_post_q = iqr(bad.post.values)
        g_med_pre = float(np.median(good.pre.values))
        g_med_post = float(np.median(good.post.values))
        b_med_pre = float(np.median(bad.pre.values))
        b_med_post = float(np.median(bad.post.values))

        # (i) Individual slopes --- use LIGHT group colour so they form
        # a soft underlay; post median trace will sit on top in full
        # saturation.
        LIGHT_GOOD = lighten(GOOD_HEX, 0.62)
        LIGHT_BAD = lighten(BAD_HEX, 0.62)
        for _, row_ in sub.iterrows():
            color = LIGHT_GOOD if row_.response_bin == "good" else LIGHT_BAD
            add_line(s, tx(0), ty(row_.pre), tx(1), ty(row_.post),
                     color, 0.8)

        # (ii) IQR rectangles at pre and post, per group --- a small
        # vertical bar spanning 25th-75th percentile.
        iqr_w = Inches(0.10)
        half = iqr_w / 2
        # good --- slight offset left at each timepoint to separate from bad
        for xv, qlo, qhi in [(0, g_pre_q[0], g_pre_q[1]),
                             (1, g_post_q[0], g_post_q[1])]:
            cx = tx(xv) - int(iqr_w)
            add_rect(s, cx - half, ty(qhi), iqr_w,
                     max(ty(qlo) - ty(qhi), 2),
                     fill=LIGHT_GOOD, line_color=GOOD, line_width=0.6)
        for xv, qlo, qhi in [(0, b_pre_q[0], b_pre_q[1]),
                             (1, b_post_q[0], b_post_q[1])]:
            cx = tx(xv) + int(iqr_w)
            add_rect(s, cx - half, ty(qhi), iqr_w,
                     max(ty(qlo) - ty(qhi), 2),
                     fill=LIGHT_BAD, line_color=BAD, line_width=0.6)

        # (iii) hollow circle markers at each individual timepoint
        # (on top of light slopes, below IQR rects? -- keep them on top
        # so the reader can identify each subject)
        for _, row_ in sub.iterrows():
            color = GOOD if row_.response_bin == "good" else BAD
            add_circle(s, tx(0), ty(row_.pre), Emu(22000),
                       fill=WHITE, line_color=color, line_width=0.9)
            add_circle(s, tx(1), ty(row_.post), Emu(22000),
                       fill=WHITE, line_color=color, line_width=0.9)

        # (iv) Bold group-median slope (GOOD + BAD) with filled diamond
        # endpoints for strong visual emphasis.
        add_line(s, tx(0) - int(iqr_w), ty(g_med_pre),
                 tx(1) - int(iqr_w), ty(g_med_post), GOOD, 2.2)
        add_diamond(s, tx(0) - int(iqr_w), ty(g_med_pre), Emu(42000),
                    fill=GOOD, line_color=WHITE, line_width=1.2)
        add_diamond(s, tx(1) - int(iqr_w), ty(g_med_post), Emu(42000),
                    fill=GOOD, line_color=WHITE, line_width=1.2)
        add_line(s, tx(0) + int(iqr_w), ty(b_med_pre),
                 tx(1) + int(iqr_w), ty(b_med_post), BAD, 2.2)
        add_diamond(s, tx(0) + int(iqr_w), ty(b_med_pre), Emu(42000),
                    fill=BAD, line_color=WHITE, line_width=1.2)
        add_diamond(s, tx(1) + int(iqr_w), ty(b_med_post), Emu(42000),
                    fill=BAD, line_color=WHITE, line_width=1.2)

        # =================================================================
        # Annotation box with within-group sign counts + between-group MW
        # =================================================================
        g = sign_df[(sign_df.factor == fname) & (sign_df.group == "good")].iloc[0]
        b = sign_df[(sign_df.factor == fname) & (sign_df.group == "bad")].iloc[0]
        try:
            mw_p = float(stats_df[(stats_df.factor == fname)
                                  & (stats_df.level == "composite")]
                         .iloc[0]["mw_p"])
        except Exception:
            mw_p = float("nan")
        arrow = "↓" if pred == "down" else "↑"

        ann_x = ax_x + Inches(0.05)
        ann_y = ax_y + Inches(0.03)
        box_w = Inches(1.35); box_h = Inches(0.58)
        add_rect(s, ann_x, ann_y, box_w, box_h,
                 fill=WHITE, line_color=GREY, line_width=0.3)
        add_text(s, ann_x + Inches(0.04), ann_y + Emu(8000),
                 box_w - Inches(0.08), Inches(0.16),
                 f"good {arrow} : {int(g.n_predicted)}/{int(g.n_total)}  "
                 f"(P = {float(g.sign_binomial_one_sided_P):.3f})",
                 size=6, color=GOOD, bold=True, anchor="top", align="left")
        add_text(s, ann_x + Inches(0.04), ann_y + Inches(0.18),
                 box_w - Inches(0.08), Inches(0.16),
                 f"bad  {arrow} : {int(b.n_predicted)}/{int(b.n_total)}  "
                 f"(P = {float(b.sign_binomial_one_sided_P):.3f})",
                 size=6, color=BAD, bold=True, anchor="top", align="left")
        add_text(s, ann_x + Inches(0.04), ann_y + Inches(0.36),
                 box_w - Inches(0.08), Inches(0.16),
                 f"MW Δ (good vs bad) P = {mw_p:.2f}",
                 size=6, color=INK, anchor="top", align="left")

        # Small direction arrow in the top-right corner of each sub-plot
        arr_x = ax_x + ax_w - Inches(0.28)
        arr_y = ax_y + Inches(0.12)
        add_text(s, arr_x, arr_y - Inches(0.07),
                 Inches(0.26), Inches(0.20),
                 f"pred {arrow}", size=7, bold=True,
                 color=RGBColor(0x66, 0x66, 0x66), align="center")

    # Shared y-axis title (rotated, at far left)
    yt = add_text(s, Inches(0.18), Inches(1.35), Inches(0.45), Inches(1.7),
                  "Composite z-score (member ssGSEA mean)",
                  size=8, color=INK, align="center")
    yt.rotation = -90

    # Legend at bottom of slide
    leg_y = Inches(4.15)
    add_circle(s, Inches(1.35), leg_y + Inches(0.06),
               Emu(30000), fill=WHITE, line_color=GOOD, line_width=1.0)
    add_line(s, Inches(1.45), leg_y + Inches(0.06),
             Inches(1.75), leg_y + Inches(0.06), GOOD, 1.1)
    add_circle(s, Inches(1.80), leg_y + Inches(0.06),
               Emu(30000), fill=WHITE, line_color=GOOD, line_width=1.0)
    add_text(s, Inches(1.88), leg_y - Emu(5000), Inches(0.9), Inches(0.18),
             "good (n=6)", size=7, color=INK, align="left")
    add_circle(s, Inches(3.2), leg_y + Inches(0.06),
               Emu(30000), fill=WHITE, line_color=BAD, line_width=1.0)
    add_line(s, Inches(3.3), leg_y + Inches(0.06),
             Inches(3.6), leg_y + Inches(0.06), BAD, 1.1)
    add_circle(s, Inches(3.65), leg_y + Inches(0.06),
               Emu(30000), fill=WHITE, line_color=BAD, line_width=1.0)
    add_text(s, Inches(3.73), leg_y - Emu(5000), Inches(0.9), Inches(0.18),
             "bad (n=6)", size=7, color=INK, align="left")


def build_B():
    """Directional-concordance grouped bar chart."""
    s = new_slide(prs_main)
    draw_panel_letter(s, "B")

    # plot area
    px = Inches(0.95); py = Inches(0.55)
    pw = Inches(4.8); ph = Inches(3.10)

    def tx(i): return int(px + (i + 0.5) / len(FACTORS) * pw)

    y_lo, y_hi = 0.0, 1.12
    def ty(v):
        return int(py + ph - (v - y_lo) / (y_hi - y_lo) * ph)

    # spines
    add_line(s, px, py, px, py + ph, LINE, 0.6)
    add_line(s, px, py + ph, px + pw, py + ph, LINE, 0.6)

    # y-axis ticks (%) -- 0, 25, 50, 75, 100
    for yv, lab in [(0, "0 %"), (0.25, "25 %"), (0.5, "50 %"),
                    (0.75, "75 %"), (1.0, "100 %")]:
        yy = ty(yv)
        add_line(s, px - Inches(0.05), yy, px, yy, LINE, 0.4)
        add_text(s, px - Inches(0.55), yy - Inches(0.09),
                 Inches(0.5), Inches(0.18),
                 lab, size=7, color=INK, align="right")

    # y-axis title (rotated)
    yt = add_text(s, Inches(0.10), py + ph / 2 - Inches(0.8),
                  Inches(0.35), Inches(1.6),
                  "Fraction of subjects moving in predicted direction",
                  size=8, color=INK, align="center")
    yt.rotation = -90

    # chance line at 0.5
    add_line(s, px + Inches(0.02), ty(0.5), px + pw - Inches(0.02),
             ty(0.5), GREY, 0.5, dashed=True)
    add_text(s, px + pw + Inches(0.02), ty(0.5) - Inches(0.08),
             Inches(0.6), Inches(0.16),
             "chance", size=6, color=INK, align="left")

    # bars --- two per factor (good, bad)
    group_w = pw / len(FACTORS)
    bar_w = int(group_w * 0.28)
    for i, (fname, pretty, pred) in enumerate(FACTORS):
        center = tx(i)
        g = sign_df[(sign_df.factor == fname) & (sign_df.group == "good")].iloc[0]
        b = sign_df[(sign_df.factor == fname) & (sign_df.group == "bad")].iloc[0]

        # good bar (left)
        xg = center - int(bar_w * 1.15)
        top_g = ty(float(g.fraction_predicted))
        add_rect(s, xg, top_g, bar_w, int(ty(0) - top_g),
                 fill=GOOD, line_color=WHITE, line_width=0.3)
        add_text(s, xg - Inches(0.15), top_g - Inches(0.24),
                 Inches(0.6), Inches(0.18),
                 f"{int(g.n_predicted)}/{int(g.n_total)}",
                 size=7, bold=True, color=GOOD, align="center")
        # highlight EMT 6/6 P=0.016 with gold asterisk
        if fname == "EMT":
            add_text(s, xg - Inches(0.15), top_g - Inches(0.40),
                     Inches(0.6), Inches(0.16),
                     f"★ P = {float(g.sign_binomial_one_sided_P):.3f}",
                     size=6, bold=True, color=HIGHLIGHT, align="center")

        # bad bar (right)
        xb = center + int(bar_w * 0.15)
        top_b = ty(float(b.fraction_predicted))
        add_rect(s, xb, top_b, bar_w, int(ty(0) - top_b),
                 fill=BAD, line_color=WHITE, line_width=0.3)
        add_text(s, xb - Inches(0.15), top_b - Inches(0.24),
                 Inches(0.6), Inches(0.18),
                 f"{int(b.n_predicted)}/{int(b.n_total)}",
                 size=7, bold=True, color=BAD, align="center")

        # factor name below
        add_text(s, center - group_w / 2 + Inches(0.08),
                 py + ph + Inches(0.12),
                 group_w - Inches(0.15), Inches(0.22),
                 pretty, size=7, color=INK, align="center")
        add_text(s, center - group_w / 2 + Inches(0.08),
                 py + ph + Inches(0.30),
                 group_w - Inches(0.15), Inches(0.16),
                 f"(pred: {pred})", size=6, color=RGBColor(0x66, 0x66, 0x66),
                 align="center")

    # x-axis title
    add_text(s, px, py + ph + Inches(0.55), pw, Inches(0.22),
             "Thread-1 baseline factor",
             size=8, color=INK, align="center")

    # legend top-right
    leg_x = px + Inches(3.6); leg_y = py + Inches(0.08)
    add_rect(s, leg_x, leg_y, Inches(0.20), Inches(0.13), fill=GOOD)
    add_text(s, leg_x + Inches(0.25), leg_y - Emu(10000),
             Inches(0.9), Inches(0.16),
             "good (n=6)", size=7, color=INK, align="left")
    add_rect(s, leg_x, leg_y + Inches(0.20), Inches(0.20),
             Inches(0.13), fill=BAD)
    add_text(s, leg_x + Inches(0.25), leg_y + Inches(0.19),
             Inches(0.9), Inches(0.16),
             "bad (n=6)", size=7, color=INK, align="left")


build_A()
build_B()
deck_main = f"{OUT}/Fig6_target_engagement_native_editable.pptx"
prs_main.save(deck_main)
print(f"wrote {deck_main}")


# ===========================================================================
# SUPP FIGURE S17 --- member-pathway detail and per-subject Δ heatmap
# ===========================================================================
prs_supp = Presentation()


def build_S17A():
    """Member-level sign-count bar: 17 signature rows (4 composites + 13
    members), for each row show the good-6 and bad-6 predicted-direction
    counts as horizontal bars."""
    s = new_slide(prs_supp)
    draw_panel_letter(s, "A")

    # build row list
    rows = []
    for fname, pretty, pred in FACTORS:
        rows.append((fname, "composite", pretty, pred, True))
        mems = (delta_df[(delta_df.factor == fname)
                         & (delta_df.member != "composite")]
                .member.unique().tolist())
        for m in sorted(mems):
            rows.append((fname, m, m, pred, False))

    n = len(rows)
    px = Inches(2.8); py = Inches(0.40)
    pw = Inches(3.45); ph = Inches(3.75)

    def ty(i):
        return int(py + (i + 0.5) / n * ph)
    def tx(v):
        return int(px + (v + 6) / 12 * pw)   # range -6 to +6

    # spines
    add_line(s, px + pw / 2, py, px + pw / 2, py + ph, LINE, 0.8)  # 0 axis
    add_line(s, px, py + ph, px + pw, py + ph, LINE, 0.6)          # x spine
    # x ticks -6..+6 by 2
    for v in [-6, -4, -2, 0, 2, 4, 6]:
        xx = tx(v)
        add_line(s, xx, int(py + ph), xx,
                 int(py + ph + Inches(0.05)), LINE, 0.4)
        add_text(s, xx - Inches(0.15), int(py + ph + Inches(0.06)),
                 Inches(0.3), Inches(0.14),
                 f"{v:+d}", size=6, color=INK, align="center")
    add_text(s, px, py + ph + Inches(0.22), pw, Inches(0.18),
             "← bad (predicted-direction count)     "
             "good (predicted-direction count) →",
             size=7, color=INK, align="center")

    # composite rows get a bold left label, members indented
    for i, (fname, mname, pretty, pred, is_comp) in enumerate(rows):
        yy = ty(i)
        # row label
        label = pretty if is_comp else (
            "  " + (pretty[:28] + "…" if len(pretty) > 28 else pretty))
        if is_comp:
            label = f"▸ {pretty}"
        add_text(s, Inches(0.25), yy - Inches(0.09),
                 px - Inches(0.30), Inches(0.16),
                 label, size=6 if is_comp else 5,
                 bold=is_comp, color=INK, align="right")

        # compute counts for this row
        if is_comp:
            sub = delta_df[(delta_df.factor == fname)
                           & (delta_df.member == "composite")]
        else:
            sub = delta_df[(delta_df.factor == fname)
                           & (delta_df.member == mname)]
        good_deltas = sub[sub.response_bin == "good"]["delta"].values
        bad_deltas = sub[sub.response_bin == "bad"]["delta"].values

        def count_predicted(arr, pred_dir):
            if pred_dir == "down":
                return int((arr < 0).sum())
            return int((arr > 0).sum())

        g_pred = count_predicted(good_deltas, pred)
        b_pred = count_predicted(bad_deltas, pred)

        # good count as right-going bar
        bar_h = int(ph / n * 0.52)
        add_rect(s, tx(0), yy - bar_h // 2,
                 tx(g_pred) - tx(0), bar_h, fill=GOOD,
                 line_color=WHITE, line_width=0.3)
        add_text(s, tx(g_pred) + Inches(0.04), yy - Inches(0.08),
                 Inches(0.3), Inches(0.14),
                 f"{g_pred}/6", size=5, color=GOOD, bold=True, align="left")
        # bad count as left-going bar
        add_rect(s, tx(-b_pred), yy - bar_h // 2,
                 tx(0) - tx(-b_pred), bar_h, fill=BAD,
                 line_color=WHITE, line_width=0.3)
        add_text(s, tx(-b_pred) - Inches(0.28), yy - Inches(0.08),
                 Inches(0.28), Inches(0.14),
                 f"{b_pred}/6", size=5, color=BAD, bold=True, align="right")

    # legend
    leg_y = Inches(4.20)
    add_rect(s, Inches(1.9), leg_y, Inches(0.16), Inches(0.12), fill=GOOD)
    add_text(s, Inches(2.08), leg_y - Emu(10000), Inches(0.9),
             Inches(0.16), "good  →  predicted direction count",
             size=6, color=INK, align="left")
    add_rect(s, Inches(1.9), leg_y + Inches(0.14), Inches(0.16),
             Inches(0.12), fill=BAD)
    add_text(s, Inches(2.08), leg_y + Inches(0.13), Inches(0.9),
             Inches(0.16), "bad  ←  predicted direction count",
             size=6, color=INK, align="left")

    # caption
    add_text(s, Inches(0.15), Inches(4.35), SLIDE_W - Inches(0.3),
             Inches(0.12),
             "▸ = composite (4 Thread-1 factors); unprefixed rows = "
             "member ssGSEA pathways (13 total)",
             size=5, color=RGBColor(0x66, 0x66, 0x66), align="center")


def build_S17B():
    """Per-subject Δ heatmap for all 17 signatures × 12 subjects grouped
    good | bad."""
    s = new_slide(prs_supp)
    draw_panel_letter(s, "B")

    # row list = 17 signatures (same ordering as S17A)
    rows = []
    for fname, pretty, pred in FACTORS:
        rows.append((fname, "composite", pretty, pred, True))
        mems = (delta_df[(delta_df.factor == fname)
                         & (delta_df.member != "composite")]
                .member.unique().tolist())
        for m in sorted(mems):
            rows.append((fname, m, m, pred, False))

    # subject order: 6 good then 6 bad, each sorted by subject_id
    good_subj = sorted(delta_df[delta_df.response_bin == "good"]
                       .subject_id.unique().tolist())
    bad_subj = sorted(delta_df[delta_df.response_bin == "bad"]
                      .subject_id.unique().tolist())
    subj_order = good_subj + bad_subj

    n_rows = len(rows)
    n_cols = len(subj_order)

    px = Inches(2.2); py = Inches(0.45)
    pw = Inches(4.00); ph = Inches(3.55)
    cell_w = pw / n_cols
    cell_h = ph / n_rows

    # diverging ramp: teal for predicted direction, coral for opposite
    def ramp(v, pred_dir):
        v = float(v)
        # orient so positive intensity = predicted direction
        oriented = -v if pred_dir == "down" else v
        lim = 2.0
        t = max(-1, min(1, oriented / lim))
        if t >= 0:
            r = int(255 - t * (255 - GOOD_HEX[0]))
            g = int(255 - t * (255 - GOOD_HEX[1]))
            b = int(255 - t * (255 - GOOD_HEX[2]))
        else:
            t = -t
            r = int(255 - t * (255 - BAD_HEX[0]))
            g = int(255 - t * (255 - BAD_HEX[1]))
            b = int(255 - t * (255 - BAD_HEX[2]))
        return RGBColor(r, g, b)

    for i, (fname, mname, pretty, pred, is_comp) in enumerate(rows):
        for j, sid in enumerate(subj_order):
            sel = delta_df[(delta_df.factor == fname)
                           & (delta_df.member == ("composite" if is_comp else mname))
                           & (delta_df.subject_id == sid)]
            v = float(sel.delta.values[0]) if not sel.empty else 0.0
            x = int(px + j * cell_w)
            y = int(py + i * cell_h)
            add_rect(s, x, y, int(cell_w), int(cell_h),
                     fill=ramp(v, pred))

    # row labels (truncated)
    for i, (fname, mname, pretty, pred, is_comp) in enumerate(rows):
        label = f"▸ {pretty}" if is_comp else pretty
        if len(label) > 30:
            label = label[:29] + "…"
        yy = int(py + (i + 0.5) * cell_h)
        add_text(s, Inches(0.25), yy - Inches(0.07),
                 px - Inches(0.30), Inches(0.14),
                 label, size=5 if not is_comp else 6,
                 bold=is_comp, color=INK, align="right")

    # column labels (subject_id) and group divider
    for j, sid in enumerate(subj_order):
        xx = int(px + (j + 0.5) * cell_w)
        color = GOOD if sid in good_subj else BAD
        add_text(s, xx - Inches(0.18), py + ph + Inches(0.03),
                 Inches(0.36), Inches(0.14),
                 str(sid), size=5, bold=True, color=color, align="center")

    # vertical divider between good and bad columns
    div_x = int(px + len(good_subj) * cell_w)
    add_line(s, div_x, py, div_x, py + ph, INK, 0.8)

    # group labels above the column headers
    add_text(s, px, py - Inches(0.18), cell_w * len(good_subj),
             Inches(0.14), f"good (n={len(good_subj)})",
             size=7, bold=True, color=GOOD, align="center")
    add_text(s, div_x, py - Inches(0.18), cell_w * len(bad_subj),
             Inches(0.14), f"bad (n={len(bad_subj)})",
             size=7, bold=True, color=BAD, align="center")

    # x-axis title
    add_text(s, px, py + ph + Inches(0.20),
             pw, Inches(0.18),
             "Subject ID",
             size=7, color=INK, align="center")

    # colour bar
    cb_x = px + pw + Inches(0.18)
    cb_y = py + Inches(0.10)
    cb_w = Inches(0.18)
    cb_h = Inches(1.0)
    for k in range(40):
        t = 1 - 2 * k / 39
        # simulate predicted-direction-oriented ramp with pred='down'
        oriented = t
        if oriented >= 0:
            r = int(255 - oriented * (255 - GOOD_HEX[0]))
            g = int(255 - oriented * (255 - GOOD_HEX[1]))
            b = int(255 - oriented * (255 - GOOD_HEX[2]))
        else:
            o = -oriented
            r = int(255 - o * (255 - BAD_HEX[0]))
            g = int(255 - o * (255 - BAD_HEX[1]))
            b = int(255 - o * (255 - BAD_HEX[2]))
        add_rect(s, cb_x, int(cb_y + k * cb_h / 40),
                 cb_w, int(cb_h / 40) + Emu(2000),
                 fill=RGBColor(r, g, b))
    add_rect(s, cb_x, cb_y, cb_w, cb_h, line_color=LINE, line_width=0.4)
    add_text(s, cb_x - Inches(0.30), cb_y - Inches(0.12),
             cb_w + Inches(0.6), Inches(0.14),
             "Δ vs predicted", size=6, color=INK, align="center")
    add_text(s, cb_x + cb_w + Inches(0.02), cb_y - Inches(0.04),
             Inches(0.3), Inches(0.12),
             "pred", size=5, color=GOOD, bold=True, align="left")
    add_text(s, cb_x + cb_w + Inches(0.02), cb_y + cb_h - Inches(0.08),
             Inches(0.3), Inches(0.12),
             "opp.", size=5, color=BAD, bold=True, align="left")
    add_text(s, cb_x + cb_w + Inches(0.02), cb_y + cb_h / 2 - Inches(0.06),
             Inches(0.3), Inches(0.12),
             "0", size=5, color=INK, align="left")

    # caption
    add_text(s, Inches(0.15), Inches(4.30),
             SLIDE_W - Inches(0.3), Inches(0.15),
             "Cells coloured by Δ (post − pre) re-oriented by predicted "
             "direction: teal = moved in predicted direction, coral = opposite.",
             size=5, color=RGBColor(0x66, 0x66, 0x66), align="center")


build_S17A()
build_S17B()
deck_supp = f"{OUT}/SuppFig_S17_target_engagement_members_native_editable.pptx"
prs_supp.save(deck_supp)
print(f"wrote {deck_supp}")
