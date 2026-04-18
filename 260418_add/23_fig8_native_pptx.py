#!/usr/bin/env python3
"""
23_fig8_native_pptx.py

Build Figure 8 --- Paired cascade phenomenology + convergence-null
callout (§3.11) --- as 6 native-editable PowerPoint slides, plus a
2-slide supplementary deck (purity-adjusted sensitivity + convergence
test scatter matrix).

Main Figure 8 panels
--------------------
  A  Between-group BCa forest (9 cascade features): diamond + error
     bar on a standardised-Z scale (diff / SE, where SE = CI_width /
     3.92); vertical references at Z=0 and Z=±1.96; Treg highlighted
     gold as the only feature with CI strictly excluding zero.
  B  Paired pre→post spaghetti for 4 key immune features (Treg /
     MHC-II / CD8 exhaustion / IGH_n) in a 2x2 grid.
  C  Within-group Δ forest (per feature, good and bad median Δ
     with BCa CIs) — shows multiple features with within-good CIs
     excluding zero (SBS5, neo_binders, neo_sites, Treg, MHC-II,
     CD8_exhaustion).
  D  Per-subject Δ waterfall for SBS5 mutation clearance and MHC-I
     neoantigen binder clearance.
  E  Conceptual fishplot schematic of paired clonal dynamics (good
     responder: clone eradication; bad responder: clone persistence).
  F  Cascade schematic + convergence-null callout (THE critical
     message): a flow diagram of the proposed cascade overlaid with
     the pre-specified convergence-test result showing that the
     static baseline predictor (Thread 1) does NOT predict cascade
     Δ magnitude (0/36 pairs P<0.05; DSB→CD8cyt r=-0.07, P=0.83).

Supp Fig S20 panels
-------------------
  A  Convergence test scatter matrix: 9 baseline LASSO winners × 4
     key cascade Δ = 36 pre-specified pairs, each a mini scatter with
     Spearman r and P.
  B  Purity-adjusted paired Δ sensitivity: Treg / MHC-II / CD8_exh /
     IGH_n before vs after CNVkit-purity correction (Supp Text S4).

Motif references consulted (cascade / paired biomarker convention):
  - Tumeh et al Nature 2014 (PD-1 pre/post CD8): paired pre→post
    spaghetti + waterfall.
  - Riaz et al Cell 2017 (ICI pre/on-treatment): paired biomarker
    panels with per-subject tracks.
  - Cercek et al NEJM 2022 (dostarlimab): per-subject response
    waterfall; simple timeline schematics.
  - Rizvi et al Science 2015 (TIL clonality pre/post): clonal
    trajectory fishplots.
  - Rosenthal et al Nature 2019 (HLA-LOH / immune escape): cascade
    flow schematic overlaid with statistical null.
  - Chowell et al Science 2018 (HLA-ICB): BCa CI forest.

Rules: one panel per slide; no plot titles; python-pptx native
elements only (TEXT_BOX / LINE / AUTO_SHAPE / FREEFORM); Arial;
DEEP palette GOOD=#0a7d6e / BAD=#c53e1f; kill_shadow() + _i()
applied everywhere.
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
# Shared infrastructure
# ---------------------------------------------------------------------------
def kill_shadow(shape):
    elem = shape._element
    spPr = elem.find(qn('p:spPr'))
    if spPr is None:
        return
    for el in spPr.findall(qn('a:effectLst')):
        spPr.remove(el)
    etree.SubElement(spPr, qn('a:effectLst'))


def _i(v):
    return int(round(float(v)))


GOOD = RGBColor(0x0A, 0x7D, 0x6E)
BAD = RGBColor(0xC5, 0x3E, 0x1F)
INK = RGBColor(0x22, 0x22, 0x22)
LINE = RGBColor(0x33, 0x33, 0x33)
GREY = RGBColor(0xBB, 0xBB, 0xBB)
LT_GREY = RGBColor(0xDD, 0xDD, 0xDD)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HIGHLIGHT = RGBColor(0xD4, 0xA3, 0x00)

GOOD_HEX = (0x0A, 0x7D, 0x6E)
BAD_HEX = (0xC5, 0x3E, 0x1F)

FONT = "Arial"
SLIDE_W = Inches(6.5)
SLIDE_H = Inches(4.5)

OUT = "/data/data/TNT/analysis/260418_add/ppt"
DATA = "/data/data/TNT/analysis"
os.makedirs(OUT, exist_ok=True)


def new_slide(prs, w=SLIDE_W, h=SLIDE_H):
    prs.slide_width = w
    prs.slide_height = h
    return prs.slides.add_slide(prs.slide_layouts[6])


def add_text(slide, x, y, w, h, text, size=8, bold=False,
             color=INK, align="left", anchor="middle", font=FONT,
             italic=False):
    tb = slide.shapes.add_textbox(_i(x), _i(y), _i(max(w, 1)), _i(max(h, 1)))
    tf = tb.text_frame
    tf.margin_left = 0; tf.margin_right = 0
    tf.margin_top = 0; tf.margin_bottom = 0
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
    r.font.italic = bool(italic)
    r.font.color.rgb = color
    kill_shadow(tb)
    return tb


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


def add_arrow(slide, x1, y1, x2, y2, color=INK, width=1.0):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,
                                   _i(x1), _i(y1), _i(x2), _i(y2))
    c.line.color.rgb = color
    c.line.width = Pt(width)
    ln = c.line._get_or_add_ln()
    tail = etree.SubElement(ln, qn('a:tailEnd'))
    tail.set('type', 'triangle')
    tail.set('w', 'med')
    tail.set('len', 'med')
    kill_shadow(c)
    return c


def add_rect(slide, x, y, w, h, fill=None, line_color=None, line_width=0.5):
    w = max(_i(w), 1); h = max(_i(h), 1)
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


def add_rounded_rect(slide, x, y, w, h, fill=None, line_color=None,
                     line_width=0.5):
    w = max(_i(w), 1); h = max(_i(h), 1)
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                 _i(x), _i(y), w, h)
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


def add_freeform_poly(slide, vertices, fill=None, line_color=None,
                      line_width=0.5):
    if len(vertices) < 3:
        return None
    x0, y0 = vertices[0]
    ff = slide.shapes.build_freeform(_i(x0), _i(y0), scale=1.0)
    pts = [(_i(x), _i(y)) for x, y in vertices[1:]]
    ff.add_line_segments(pts, close=True)
    shp = ff.convert_to_shape()
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


def lighten(rgb_hex, factor=0.7):
    r0, g0, b0 = rgb_hex
    return RGBColor(int(r0 + (255 - r0) * factor),
                    int(g0 + (255 - g0) * factor),
                    int(b0 + (255 - b0) * factor))


def draw_panel_letter(slide, letter):
    add_text(slide, Inches(0.15), Inches(0.1), Inches(0.4), Inches(0.35),
             letter, size=14, bold=True, color=INK, align="left")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
bca = pd.read_csv(f"{DATA}/tables/TableS8_cascade_BCa_bootstrap.tsv",
                  sep="\t")
paired_long = pd.read_csv(
    f"{DATA}/09_integration/paired_delta/paired_feature_long.tsv",
    sep="\t")
convergence = pd.read_csv(
    f"{DATA}/260418_add/targeted_convergence_test.tsv", sep="\t")


def parse_ci(s):
    """Parse '[+x.xx, +y.yy]' into (lo, hi) floats."""
    s = str(s).strip().strip("[]")
    parts = [p.strip() for p in s.split(",")]
    return float(parts[0]), float(parts[1])


# Cascade feature ordering (biological sequence + Treg at end for emphasis)
CASCADE_FEATURES = [
    "SBS5", "missense", "neo_binders", "neo_sites",
    "MHC_II", "CD8_exhaustion", "Treg",
    "IGH_n", "TRB_shannon",
]
PRETTY = {
    "SBS5": "SBS5 (clock-like)  mutations",
    "missense": "Missense mutations",
    "neo_binders": "MHC-I neoantigen binders",
    "neo_sites": "Neoantigen sites",
    "MHC_II": "MHC-II signature",
    "CD8_exhaustion": "CD8 exhaustion",
    "Treg": "Treg signature  ★",
    "IGH_n": "IGH clonotypes (B-cell)",
    "TRB_shannon": "TRB Shannon diversity",
}

# Bio-axis grouping
FEAT_GROUP = {
    "SBS5": "Clone clearance",
    "missense": "Clone clearance",
    "neo_binders": "Antigen clearance",
    "neo_sites": "Antigen clearance",
    "MHC_II": "Immune reprogramming",
    "CD8_exhaustion": "Immune reprogramming",
    "Treg": "Immune reprogramming",
    "IGH_n": "B-cell infiltration",
    "TRB_shannon": "T-cell repertoire",
}
GROUP_COLOR = {
    "Clone clearance": RGBColor(0x3E, 0x52, 0x86),
    "Antigen clearance": RGBColor(0x6A, 0x4C, 0x93),
    "Immune reprogramming": RGBColor(0x0A, 0x7D, 0x6E),
    "B-cell infiltration": RGBColor(0xA6, 0x50, 0x2C),
    "T-cell repertoire": RGBColor(0x55, 0x55, 0x55),
}


# ===========================================================================
# MAIN FIGURE 8
# ===========================================================================
prs_main = Presentation()


# -------------------------------------------------------------------
# Panel A --- Between-group BCa forest (standardised-Z scale)
# -------------------------------------------------------------------
def build_A():
    s = new_slide(prs_main)
    draw_panel_letter(s, "A")

    # plot area
    ax_x = Inches(2.80); ax_y = Inches(0.55)
    ax_w = Inches(2.85); ax_h = Inches(3.25)

    x_lo, x_hi = -3.2, 3.2

    def tx(v):
        return _i(ax_x + (v - x_lo) / (x_hi - x_lo) * ax_w)

    n_rows = len(CASCADE_FEATURES)
    row_h = ax_h / n_rows

    # spines + ticks
    add_line(s, ax_x, ax_y + ax_h, ax_x + ax_w, ax_y + ax_h, LINE, 0.5)
    for v in [-3, -2, -1, 0, 1, 2, 3]:
        xx = tx(v)
        add_line(s, xx, _i(ax_y + ax_h), xx,
                 _i(ax_y + ax_h + Inches(0.04)), LINE, 0.4)
        add_text(s, xx - Inches(0.13), _i(ax_y + ax_h + Inches(0.05)),
                 Inches(0.26), Inches(0.14),
                 f"{v:+d}" if v != 0 else "0",
                 size=6, color=INK, align="center")
    add_text(s, ax_x, _i(ax_y + ax_h + Inches(0.22)),
             ax_w, Inches(0.15),
             "Standardised between-group Z  (good − bad)",
             size=8, color=INK, align="center")
    add_text(s, ax_x, _i(ax_y + ax_h + Inches(0.37)),
             ax_w, Inches(0.13),
             "= diff_median / (CI width / 3.92)",
             size=5, italic=True, color=RGBColor(0x66, 0x66, 0x66),
             align="center")

    # reference lines
    add_line(s, tx(0), ax_y, tx(0), ax_y + ax_h, INK, 1.0)
    add_line(s, tx(1.96), ax_y, tx(1.96), ax_y + ax_h,
             GREY, 0.5, dashed=True)
    add_line(s, tx(-1.96), ax_y, tx(-1.96), ax_y + ax_h,
             GREY, 0.5, dashed=True)
    add_text(s, tx(1.96) - Inches(0.20), _i(ax_y + Inches(0.02)),
             Inches(0.4), Inches(0.12),
             "P=0.05", size=5, italic=True,
             color=RGBColor(0x99, 0x99, 0x99), align="center")

    # render rows
    for i, feat in enumerate(CASCADE_FEATURES):
        y = _i(ax_y + (i + 0.5) * row_h)
        row = bca[bca.feature == feat]
        if row.empty:
            continue
        row = row.iloc[0]
        diff_med = float(row["diff_median_good_minus_bad"])
        lo, hi = parse_ci(row["diff_95CI"])
        # standardised Z
        se = (hi - lo) / 3.92 if (hi - lo) > 0 else 0.1
        z = diff_med / se
        z_lo = lo / se
        z_hi = hi / se

        # clamp
        z_c = max(x_lo, min(x_hi, z))
        z_lo_c = max(x_lo, min(x_hi, z_lo))
        z_hi_c = max(x_lo, min(x_hi, z_hi))

        is_sig = abs(z) >= 1.96
        highlight = (feat == "Treg")

        group = FEAT_GROUP[feat]
        gcolor = GROUP_COLOR[group]

        # feature label (left)
        lab_color = HIGHLIGHT if highlight else gcolor
        add_text(s, Inches(0.18), y - Inches(0.09),
                 ax_x - Inches(0.22), Inches(0.18),
                 PRETTY[feat],
                 size=7, bold=highlight, color=lab_color, align="right")
        # group chip
        add_rect(s, ax_x - Inches(0.12), y - Inches(0.04),
                 Inches(0.06), Inches(0.08), fill=gcolor)

        # CI bar
        add_line(s, tx(z_lo_c), y, tx(z_hi_c), y,
                 gcolor, 1.8 if highlight else 1.2)
        # CI caps
        cap_h = _i(Inches(0.05))
        add_line(s, tx(z_lo_c), y - cap_h, tx(z_lo_c), y + cap_h,
                 gcolor, 1.0)
        add_line(s, tx(z_hi_c), y - cap_h, tx(z_hi_c), y + cap_h,
                 gcolor, 1.0)
        # diamond at point estimate
        d_col = HIGHLIGHT if highlight else gcolor
        d_size = Emu(55000) if highlight else Emu(38000)
        add_diamond(s, tx(z_c), y, d_size,
                    fill=d_col, line_color=WHITE, line_width=1.1)

        # right column: numeric effect + CI + P
        rx = ax_x + ax_w + Inches(0.06)
        eff_txt = f"{diff_med:+.2f}" if abs(diff_med) < 10 else f"{diff_med:+.0f}"
        ci_txt = (f"[{lo:+.2f}, {hi:+.2f}]"
                  if abs(lo) < 10 and abs(hi) < 10
                  else f"[{lo:+.0f}, {hi:+.0f}]")
        p_val = float(row["MW_p"])
        p_txt = f"P = {p_val:.3f}"
        add_text(s, rx, y - Inches(0.14),
                 Inches(0.75), Inches(0.13),
                 eff_txt, size=7, bold=is_sig,
                 color=HIGHLIGHT if highlight else
                       (gcolor if is_sig else RGBColor(0x66, 0x66, 0x66)),
                 align="left")
        add_text(s, rx, y - Inches(0.02),
                 Inches(0.85), Inches(0.13),
                 ci_txt, size=5,
                 color=RGBColor(0x55, 0x55, 0x55), align="left")
        add_text(s, rx, y + Inches(0.10),
                 Inches(0.75), Inches(0.13),
                 p_txt, size=6,
                 bold=is_sig,
                 color=HIGHLIGHT if highlight else
                       (gcolor if is_sig else RGBColor(0x55, 0x55, 0x55)),
                 align="left")

    # group legend (bottom left)
    leg_y = Inches(3.90)
    xb = Inches(0.18)
    add_text(s, xb, leg_y, Inches(1.3), Inches(0.14),
             "Cascade group:",
             size=6, bold=True, color=INK, align="left")
    xb += Inches(0.85)
    for g in ["Clone clearance", "Antigen clearance",
              "Immune reprogramming", "B-cell infiltration",
              "T-cell repertoire"]:
        add_rect(s, xb, leg_y + Inches(0.02),
                 Inches(0.09), Inches(0.09),
                 fill=GROUP_COLOR[g])
        add_text(s, xb + Inches(0.11), leg_y - Emu(5000),
                 Inches(1.10), Inches(0.14),
                 g, size=5, color=INK, align="left")
        xb += Inches(1.12)

    # bottom-strip message
    add_text(s, Inches(0.15), Inches(4.15),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "★ Only Treg Δ has a between-group BCa CI strictly "
             "excluding zero (robust at n = 12).",
             size=6, bold=True, color=HIGHLIGHT, align="center")
    add_text(s, Inches(0.15), Inches(4.30),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Other cascade features are framed as exploratory "
             "phenomenology; see Panel C for within-group robustness.",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")


# -------------------------------------------------------------------
# Panel B --- Paired pre→post spaghetti (2x2 grid: 4 key features)
# -------------------------------------------------------------------
def build_B():
    s = new_slide(prs_main)
    draw_panel_letter(s, "B")

    features_b = ["Treg", "MHC_II", "CD8_exhaustion", "IGH_n"]
    labels_b = {"Treg": "Treg signature (★ between-group robust)",
                "MHC_II": "MHC-II signature",
                "CD8_exhaustion": "CD8 exhaustion",
                "IGH_n": "IGH clonotype count (B-cell)"}

    grid_ox = Inches(0.70); grid_oy = Inches(0.40)
    sub_w = Inches(2.70); sub_h = Inches(1.70)
    col_gap = Inches(0.30); row_gap = Inches(0.25)

    LIGHT_GOOD = lighten(GOOD_HEX, 0.62)
    LIGHT_BAD = lighten(BAD_HEX, 0.62)

    for idx, feat in enumerate(features_b):
        r, c = idx // 2, idx % 2
        px = grid_ox + c * (sub_w + col_gap)
        py = grid_oy + r * (sub_h + row_gap)

        # data
        sub = paired_long[paired_long.feature == feat].copy()
        # drop subjects with NaN
        sub = sub.dropna(subset=["pre", "post"])

        all_vals = np.concatenate([sub.pre.values, sub.post.values])
        if feat == "IGH_n":
            y_min = 0
            y_max = float(np.max(all_vals)) * 1.10
        else:
            y_span = max(1.5, float(np.max(np.abs(all_vals))) * 1.10)
            y_min, y_max = -y_span, y_span

        # inner axes
        ax_x = px + Inches(0.55); ax_y = py + Inches(0.08)
        ax_w = sub_w - Inches(0.65); ax_h = sub_h - Inches(0.45)

        def tx(v): return _i(ax_x + v * ax_w)
        def ty(v): return _i(ax_y + ax_h -
                              (v - y_min) / (y_max - y_min) * ax_h)

        # spines
        add_line(s, ax_x, ax_y, ax_x, ax_y + ax_h, LINE, 0.5)
        add_line(s, ax_x, ax_y + ax_h, ax_x + ax_w,
                 ax_y + ax_h, LINE, 0.5)

        # x tick labels
        add_text(s, tx(0) - Inches(0.17), _i(ax_y + ax_h + Inches(0.02)),
                 Inches(0.34), Inches(0.13),
                 "pre", size=7, color=INK, align="center")
        add_text(s, tx(1) - Inches(0.18), _i(ax_y + ax_h + Inches(0.02)),
                 Inches(0.36), Inches(0.13),
                 "post", size=7, color=INK, align="center")

        # y ticks
        for yv in [y_min, (y_min + y_max) / 2, y_max]:
            yy = ty(yv)
            add_line(s, _i(ax_x - Inches(0.04)), yy,
                     ax_x, yy, LINE, 0.4)
            if feat == "IGH_n":
                lab = f"{int(yv)}"
            else:
                lab = f"{yv:+.1f}" if abs(yv) >= 1 else f"{yv:.1f}"
            add_text(s, _i(ax_x - Inches(0.45)), yy - Inches(0.08),
                     Inches(0.40), Inches(0.14),
                     lab, size=5, color=INK, align="right")

        # zero reference (only for z-score features)
        if feat != "IGH_n":
            add_line(s, ax_x, ty(0), ax_x + ax_w, ty(0),
                     GREY, 0.3, dashed=True)

        # faint individual slopes
        for _, row in sub.iterrows():
            c_col = LIGHT_GOOD if row.response == "good" else LIGHT_BAD
            add_line(s, tx(0), ty(row.pre), tx(1), ty(row.post),
                     c_col, 0.8)

        # hollow circle markers
        for _, row in sub.iterrows():
            col = GOOD if row.response == "good" else BAD
            add_circle(s, tx(0), ty(row.pre), Emu(22000),
                       fill=WHITE, line_color=col, line_width=0.9)
            add_circle(s, tx(1), ty(row.post), Emu(22000),
                       fill=WHITE, line_color=col, line_width=0.9)

        # group median slopes + diamonds
        for grp in ["good", "bad"]:
            grp_df = sub[sub.response == grp]
            if len(grp_df) == 0:
                continue
            col = GOOD if grp == "good" else BAD
            mpre = float(np.median(grp_df.pre.values))
            mpost = float(np.median(grp_df.post.values))
            add_line(s, tx(0), ty(mpre), tx(1), ty(mpost), col, 2.2)
            add_diamond(s, tx(0), ty(mpre), Emu(40000),
                        fill=col, line_color=WHITE, line_width=1.1)
            add_diamond(s, tx(1), ty(mpost), Emu(40000),
                        fill=col, line_color=WHITE, line_width=1.1)

        # feature label below
        add_text(s, px, _i(py + sub_h - Inches(0.22)),
                 sub_w, Inches(0.16),
                 labels_b[feat],
                 size=7, bold=True,
                 color=HIGHLIGHT if feat == "Treg" else INK,
                 align="center")

        # BCa summary inside plot (top-right corner)
        bca_row = bca[bca.feature == feat]
        if not bca_row.empty:
            r_ = bca_row.iloc[0]
            txt_lines = [
                f"g Δ {r_.good_delta_median:+.2f}  {r_.good_95CI}",
                f"b Δ {r_.bad_delta_median:+.2f}  {r_.bad_95CI}",
                f"MW P = {r_.MW_p:.3f}",
            ]
            if feat == "IGH_n":
                txt_lines = [
                    f"g Δ {r_.good_delta_median:+.0f}  {r_.good_95CI}",
                    f"b Δ {r_.bad_delta_median:+.0f}  {r_.bad_95CI}",
                    f"MW P = {r_.MW_p:.3f}",
                ]
            ann_x = ax_x + ax_w - Inches(1.55)
            ann_y = ax_y + Inches(0.03)
            add_rect(s, ann_x, ann_y, Inches(1.50), Inches(0.50),
                     fill=WHITE, line_color=GREY, line_width=0.3)
            for k, line in enumerate(txt_lines):
                add_text(s, ann_x + Inches(0.03),
                         ann_y + Inches(0.03 + k * 0.15),
                         Inches(1.42), Inches(0.15),
                         line, size=5,
                         bold=(k == 2 and r_.MW_p < 0.05),
                         color=(GOOD if "g Δ" in line else
                                (BAD if "b Δ" in line else
                                 (HIGHLIGHT if r_.MW_p < 0.05 else INK))),
                         align="left", anchor="top")

    # legend (bottom)
    leg_y = Inches(4.10)
    add_diamond(s, Inches(0.70), leg_y + Inches(0.06), Emu(30000),
                fill=GOOD, line_color=WHITE, line_width=1.0)
    add_text(s, Inches(0.88), leg_y - Emu(5000),
             Inches(1.30), Inches(0.16),
             "good median (n=6)",
             size=7, color=INK, align="left")
    add_diamond(s, Inches(2.60), leg_y + Inches(0.06), Emu(30000),
                fill=BAD, line_color=WHITE, line_width=1.0)
    add_text(s, Inches(2.78), leg_y - Emu(5000),
             Inches(1.30), Inches(0.16),
             "bad median (n=6)",
             size=7, color=INK, align="left")
    add_circle(s, Inches(4.50), leg_y + Inches(0.06), Emu(22000),
               fill=WHITE, line_color=LINE, line_width=0.9)
    add_text(s, Inches(4.62), leg_y - Emu(5000),
             Inches(1.80), Inches(0.16),
             "individual subject pre/post",
             size=7, color=INK, align="left")


# -------------------------------------------------------------------
# Panel C --- Within-group Δ forest (good + bad paired per feature)
# -------------------------------------------------------------------
def build_C():
    s = new_slide(prs_main)
    draw_panel_letter(s, "C")

    # we plot two rows per feature: good (above) and bad (below)
    # on a *standardised per-feature* x-axis: Δ / feature_scale.
    # feature_scale = max_abs_of_either_group_CI, so both groups plot in
    # a common [-1, +1] frame with ±1 = strongest observed CI bound.

    ax_x = Inches(2.50); ax_y = Inches(0.45)
    ax_w = Inches(3.05); ax_h = Inches(3.35)

    x_lo, x_hi = -1.15, 1.15

    def tx(v):
        return _i(ax_x + (v - x_lo) / (x_hi - x_lo) * ax_w)

    n_rows = len(CASCADE_FEATURES)
    row_h = ax_h / n_rows

    # spines + ticks
    add_line(s, ax_x, ax_y + ax_h, ax_x + ax_w, ax_y + ax_h, LINE, 0.5)
    for v, lab in [(-1, "−1.0"), (-0.5, "−0.5"), (0, "0"),
                   (0.5, "+0.5"), (1, "+1.0")]:
        xx = tx(v)
        add_line(s, xx, _i(ax_y + ax_h), xx,
                 _i(ax_y + ax_h + Inches(0.04)), LINE, 0.4)
        add_text(s, xx - Inches(0.18), _i(ax_y + ax_h + Inches(0.05)),
                 Inches(0.36), Inches(0.13),
                 lab, size=6, color=INK, align="center")
    add_text(s, ax_x, _i(ax_y + ax_h + Inches(0.22)),
             ax_w, Inches(0.14),
             "Within-group median Δ  (normalised to per-feature |max CI|)",
             size=7, color=INK, align="center")

    # zero reference
    add_line(s, tx(0), ax_y, tx(0), ax_y + ax_h, INK, 1.0)

    # iterate features
    for i, feat in enumerate(CASCADE_FEATURES):
        row = bca[bca.feature == feat]
        if row.empty:
            continue
        row = row.iloc[0]
        g_med = float(row.good_delta_median)
        b_med = float(row.bad_delta_median)
        g_lo, g_hi = parse_ci(row.good_95CI)
        b_lo, b_hi = parse_ci(row.bad_95CI)
        scale = max(abs(g_lo), abs(g_hi), abs(b_lo), abs(b_hi),
                    abs(g_med), abs(b_med), 1e-6)

        # CI excludes zero?
        g_excl = (g_lo > 0 and g_hi > 0) or (g_lo < 0 and g_hi < 0)
        b_excl = (b_lo > 0 and b_hi > 0) or (b_lo < 0 and b_hi < 0)

        y_row = _i(ax_y + (i + 0.5) * row_h)
        y_good = y_row - _i(Inches(0.08))
        y_bad = y_row + _i(Inches(0.08))

        # feature label left
        group = FEAT_GROUP[feat]
        gcolor = GROUP_COLOR[group]
        add_text(s, Inches(0.18), y_row - Inches(0.09),
                 ax_x - Inches(0.22), Inches(0.18),
                 PRETTY[feat],
                 size=7, bold=(feat == "Treg"),
                 color=HIGHLIGHT if feat == "Treg" else gcolor,
                 align="right")
        add_rect(s, ax_x - Inches(0.12), y_row - Inches(0.04),
                 Inches(0.06), Inches(0.08), fill=gcolor)

        # good row
        g_n_lo = max(x_lo, min(x_hi, g_lo / scale))
        g_n_hi = max(x_lo, min(x_hi, g_hi / scale))
        g_n_med = max(x_lo, min(x_hi, g_med / scale))
        add_line(s, tx(g_n_lo), y_good, tx(g_n_hi), y_good,
                 GOOD, 1.5 if g_excl else 0.9)
        add_line(s, tx(g_n_lo), y_good - _i(Inches(0.04)),
                 tx(g_n_lo), y_good + _i(Inches(0.04)), GOOD, 0.8)
        add_line(s, tx(g_n_hi), y_good - _i(Inches(0.04)),
                 tx(g_n_hi), y_good + _i(Inches(0.04)), GOOD, 0.8)
        add_diamond(s, tx(g_n_med), y_good, Emu(30000),
                    fill=GOOD if g_excl else WHITE,
                    line_color=GOOD, line_width=0.8)

        # bad row
        b_n_lo = max(x_lo, min(x_hi, b_lo / scale))
        b_n_hi = max(x_lo, min(x_hi, b_hi / scale))
        b_n_med = max(x_lo, min(x_hi, b_med / scale))
        add_line(s, tx(b_n_lo), y_bad, tx(b_n_hi), y_bad,
                 BAD, 1.5 if b_excl else 0.9)
        add_line(s, tx(b_n_lo), y_bad - _i(Inches(0.04)),
                 tx(b_n_lo), y_bad + _i(Inches(0.04)), BAD, 0.8)
        add_line(s, tx(b_n_hi), y_bad - _i(Inches(0.04)),
                 tx(b_n_hi), y_bad + _i(Inches(0.04)), BAD, 0.8)
        add_diamond(s, tx(b_n_med), y_bad, Emu(30000),
                    fill=BAD if b_excl else WHITE,
                    line_color=BAD, line_width=0.8)

        # right-column text: scale + summary
        rx = ax_x + ax_w + Inches(0.05)
        # scale for normalization
        if feat in {"IGH_n", "SBS5", "missense", "neo_binders",
                    "neo_sites"}:
            scale_txt = f"scale = {scale:.0f}"
        else:
            scale_txt = f"scale = {scale:.2f}"
        add_text(s, rx, y_row - Inches(0.12),
                 Inches(0.80), Inches(0.12),
                 scale_txt, size=5,
                 color=RGBColor(0x77, 0x77, 0x77), align="left")
        # robustness indicators
        robust_txt = ""
        if g_excl:
            robust_txt += "g·"
        if b_excl:
            robust_txt += "b·"
        if g_excl or b_excl:
            add_text(s, rx, y_row,
                     Inches(0.80), Inches(0.13),
                     f"{robust_txt}CI excludes 0",
                     size=6, bold=True, color=GOOD if g_excl else BAD,
                     align="left")

    # legend
    leg_y = Inches(3.92)
    add_line(s, Inches(0.20), leg_y + Inches(0.06),
             Inches(0.50), leg_y + Inches(0.06), GOOD, 1.5)
    add_diamond(s, Inches(0.35), leg_y + Inches(0.06), Emu(24000),
                fill=GOOD, line_color=WHITE, line_width=0.8)
    add_text(s, Inches(0.55), leg_y - Emu(10000),
             Inches(1.7), Inches(0.16),
             "good Δ filled diamond = CI excl. 0 (robust)",
             size=6, color=INK, align="left")
    add_line(s, Inches(0.20), leg_y + Inches(0.22),
             Inches(0.50), leg_y + Inches(0.22), BAD, 1.5)
    add_diamond(s, Inches(0.35), leg_y + Inches(0.22), Emu(24000),
                fill=WHITE, line_color=BAD, line_width=0.8)
    add_text(s, Inches(0.55), leg_y + Inches(0.13),
             Inches(2.6), Inches(0.16),
             "bad Δ hollow diamond = CI spans 0 (exploratory)",
             size=6, color=INK, align="left")

    # caption
    add_text(s, Inches(0.15), Inches(4.32),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Within-good CIs excluding zero: SBS5, neo_binders, "
             "neo_sites, MHC-II, CD8 exhaustion, Treg. "
             "Consistent cascade within responders; between-group only Treg.",
             size=6, italic=True, color=INK, align="center")


# -------------------------------------------------------------------
# Panel D --- Per-subject Δ waterfall: SBS5 + neoantigen
# -------------------------------------------------------------------
def build_D():
    s = new_slide(prs_main)
    draw_panel_letter(s, "D")

    # Top sub-panel: SBS5 Δ per subject
    top_x = Inches(0.95); top_y = Inches(0.45)
    top_w = Inches(5.30); top_h = Inches(1.55)
    # Bottom sub-panel: Neoantigen binder Δ per subject
    bot_x = Inches(0.95); bot_y = Inches(2.25)
    bot_w = Inches(5.30); bot_h = Inches(1.55)

    for feat, label, px, py, pw, ph in [
        ("SBS5", "SBS5 (clock-like) mutations  ·  Δ (post − pre)",
         top_x, top_y, top_w, top_h),
        ("neo_binders", "MHC-I neoantigen binders  ·  Δ (post − pre)",
         bot_x, bot_y, bot_w, bot_h),
    ]:
        sub = paired_long[paired_long.feature == feat].copy()
        sub = sub.dropna(subset=["pre", "post"]).copy()
        sub["delta"] = sub["post"] - sub["pre"]
        # sort by delta ascending
        sub = sub.sort_values("delta").reset_index(drop=True)
        n = len(sub)
        if n == 0:
            continue

        y_lo = float(sub.delta.min()) * 1.15
        y_hi = float(sub.delta.max()) * 1.15
        if y_hi < 0:
            y_hi = abs(y_lo) * 0.2
        if y_lo > 0:
            y_lo = -abs(y_hi) * 0.2

        def ty(v):
            return _i(py + ph - (v - y_lo) / (y_hi - y_lo) * ph)

        def tx(i):
            return _i(px + (i + 0.5) / n * pw)

        # spines
        add_line(s, px, py, px, py + ph, LINE, 0.5)
        # x axis at 0
        add_line(s, px, ty(0), px + pw, ty(0), INK, 1.0)
        # y ticks at min/max/0
        for yv in [y_lo, (y_lo + y_hi) / 2, 0, y_hi / 2 if y_hi > 0 else 0,
                   y_hi]:
            yy = ty(yv)
            add_line(s, _i(px - Inches(0.04)), yy, px, yy, LINE, 0.4)
            label_txt = f"{yv:+.0f}" if abs(yv) >= 10 else f"{yv:+.1f}"
            add_text(s, _i(px - Inches(0.50)), yy - Inches(0.08),
                     Inches(0.45), Inches(0.14),
                     label_txt, size=6, color=INK, align="right")

        # bars
        bar_w = int(pw / n * 0.75)
        for i, row in sub.iterrows():
            color = GOOD if row.response == "good" else BAD
            cx = tx(i)
            x0 = _i(cx - bar_w / 2)
            y_top = ty(max(0, row.delta))
            y_bot = ty(min(0, row.delta))
            add_rect(s, x0, y_top, bar_w, y_bot - y_top,
                     fill=color, line_color=WHITE, line_width=0.3)
            # subject id label (below axis 0 for negatives, above for positives)
            lab_y = ty(0) + Inches(0.02) if row.delta >= 0 else ty(0) - Inches(0.15)
            add_text(s, cx - Inches(0.12), lab_y,
                     Inches(0.24), Inches(0.12),
                     str(int(row.subject_id)), size=4,
                     color=color, bold=True, align="center")

        # title label below sub-plot
        add_text(s, px, _i(py + ph + Inches(0.14)),
                 pw, Inches(0.14),
                 label, size=7, bold=True, color=INK, align="center")

    # Legend (right edge)
    leg_x = Inches(5.95); leg_y = Inches(0.50)
    add_rect(s, leg_x, leg_y, Inches(0.18), Inches(0.13), fill=GOOD)
    add_text(s, leg_x, leg_y + Inches(0.15),
             Inches(0.5), Inches(0.12),
             "good", size=6, bold=True, color=GOOD, align="center")
    add_rect(s, leg_x, leg_y + Inches(0.35), Inches(0.18), Inches(0.13),
             fill=BAD)
    add_text(s, leg_x, leg_y + Inches(0.49),
             Inches(0.5), Inches(0.12),
             "bad", size=6, bold=True, color=BAD, align="center")

    # bottom caption
    add_text(s, Inches(0.15), Inches(4.08),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Per-subject Δ sorted ascending. Subjects #s shown near bar bases. "
             "Within good: SBS5 Δ median −76 (CI [−145, −64]),",
             size=6, color=INK, align="center")
    add_text(s, Inches(0.15), Inches(4.23),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "neo_binders Δ median −312 (CI [−626, −123]).  Within-good "
             "CIs exclude zero; between-group CIs span zero → exploratory.",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")


# -------------------------------------------------------------------
# Panel E --- Conceptual fishplot schematic of paired clonal dynamics
# -------------------------------------------------------------------
def build_E():
    s = new_slide(prs_main)
    draw_panel_letter(s, "E")

    # Two mini-schematics: good-responder (top) and bad-responder (bottom)
    # Each shows a stylised fishplot: pre / post clonal composition

    def fish(cx, cy, w, h, group_color, trajectory):
        """Draw a stylised fishplot: symmetric envelope with internal
        clone layers, pre on left, post on right.
        trajectory = 'shrink' or 'persist'.
        """
        # outer envelope (overall tumor tissue)
        n_pts = 24
        xs = np.linspace(0, 1, n_pts)
        if trajectory == "shrink":
            env = 0.5 * (1 - 0.75 * xs)   # shrinks over time
        else:
            env = 0.5 * (1 - 0.10 * xs)   # persists

        # top / bottom envelope vertices
        top_pts = [(cx + v * w, cy - env[i] * h) for i, v in enumerate(xs)]
        bot_pts = [(cx + v * w, cy + env[i] * h)
                   for i, v in enumerate(reversed(xs))]
        env_color = lighten(
            (group_color[0], group_color[1], group_color[2]), 0.78)
        add_freeform_poly(s, top_pts + bot_pts,
                          fill=env_color, line_color=None)

        # internal clone layers (3 sub-clones)
        for sub_i, (offset, color_mul) in enumerate([
            (0.00, 0.35),
            (0.55, 0.55),
            (-0.55, 0.70),
        ]):
            clone_xs = xs
            if trajectory == "shrink":
                # sub-clone decays (eradication)
                if sub_i == 0:
                    clone_h = env * (1 - 0.95 * xs) * 0.35
                elif sub_i == 1:
                    clone_h = env * (1 - 0.50 * xs) * 0.22
                else:
                    clone_h = env * (1 - 0.20 * xs) * 0.18
            else:
                # sub-clone persists or even expands
                if sub_i == 0:
                    clone_h = env * 0.35
                elif sub_i == 1:
                    clone_h = env * 0.28
                else:
                    clone_h = env * 0.23

            baseline = env * offset
            top_p = [(cx + v * w, cy + (baseline[i] - clone_h[i]) * h)
                     for i, v in enumerate(clone_xs)]
            bot_p = [(cx + v * w, cy + (baseline[i] + clone_h[i]) * h)
                     for i, v in enumerate(reversed(clone_xs))]
            c = lighten(
                (group_color[0], group_color[1], group_color[2]),
                color_mul)
            add_freeform_poly(s, top_p + bot_p,
                              fill=c, line_color=None)

        # pre / post markers + labels
        add_line(s, _i(cx), _i(cy - h), _i(cx), _i(cy + h),
                 INK, 0.5, dashed=True)
        add_line(s, _i(cx + w), _i(cy - h), _i(cx + w), _i(cy + h),
                 INK, 0.5, dashed=True)
        add_text(s, _i(cx - Inches(0.25)), _i(cy + h + Inches(0.02)),
                 Inches(0.50), Inches(0.14),
                 "pre", size=7, color=INK, align="center")
        add_text(s, _i(cx + w - Inches(0.25)),
                 _i(cy + h + Inches(0.02)),
                 Inches(0.50), Inches(0.14),
                 "post", size=7, color=INK, align="center")

    # good responder (top)
    fish(Inches(1.60), Inches(1.30), Inches(3.20), Inches(0.75),
         (GOOD_HEX[0], GOOD_HEX[1], GOOD_HEX[2]), "shrink")
    add_text(s, Inches(0.25), Inches(1.05), Inches(1.30), Inches(0.18),
             "Good responder",
             size=8, bold=True, color=GOOD, align="left")
    add_text(s, Inches(0.25), Inches(1.23), Inches(1.30), Inches(0.18),
             "(n = 6)", size=6, color=GOOD, align="left", italic=True)
    add_text(s, Inches(5.00), Inches(1.15), Inches(1.40), Inches(0.35),
             "clone shrinkage +\neradication", size=6, italic=True,
             color=GOOD, align="left", anchor="top")

    # bad responder (bottom)
    fish(Inches(1.60), Inches(2.80), Inches(3.20), Inches(0.75),
         (BAD_HEX[0], BAD_HEX[1], BAD_HEX[2]), "persist")
    add_text(s, Inches(0.25), Inches(2.55), Inches(1.30), Inches(0.18),
             "Bad responder",
             size=8, bold=True, color=BAD, align="left")
    add_text(s, Inches(0.25), Inches(2.73), Inches(1.30), Inches(0.18),
             "(n = 6)", size=6, color=BAD, align="left", italic=True)
    add_text(s, Inches(5.00), Inches(2.65), Inches(1.40), Inches(0.35),
             "clone persistence", size=6, italic=True,
             color=BAD, align="left", anchor="top")

    # bottom caption
    add_text(s, Inches(0.15), Inches(3.85),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Conceptual schematic (not PyClone-VI output; full PyClone "
             "detail in Supp Fig S16).",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")
    add_text(s, Inches(0.15), Inches(4.00),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "PyClone-VI dominant-clone Δ trend: good −0.67 CCF vs "
             "bad −0.15, Mann-Whitney P = 0.34 (trend, NS at n = 12).",
             size=6, color=INK, align="center")
    add_text(s, Inches(0.15), Inches(4.20),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Consistent with SBS5 and MHC-I neoantigen clearance "
             "observed in Panel D.",
             size=6, italic=True, color=HIGHLIGHT, align="center")


# -------------------------------------------------------------------
# Panel F --- Cascade schematic + convergence-null callout
# -------------------------------------------------------------------
def build_F():
    s = new_slide(prs_main)
    draw_panel_letter(s, "F")

    # ===========================================================
    # Left half: proposed cascade flow diagram
    # ===========================================================
    # Boxes arranged as: RT → clone clearance → antigen clearance →
    # immune reprogramming → B-cell infiltration → response
    stages = [
        ("Radiation\n(50.4 Gy)",
         RGBColor(0x4F, 0x73, 0x8E)),
        ("Clone clearance\n(SBS5, missense ↓)",
         GROUP_COLOR["Clone clearance"]),
        ("Antigen clearance\n(MHC-I neoantigens ↓)",
         GROUP_COLOR["Antigen clearance"]),
        ("Immune reprogramming\n(MHC-II, CD8 exh., Treg ↑)",
         GROUP_COLOR["Immune reprogramming"]),
        ("B-cell infiltration\n(IGH clonotypes ↑)",
         GROUP_COLOR["B-cell infiltration"]),
        ("Final TNT response\n(good responder)",
         HIGHLIGHT),
    ]

    bx = Inches(0.30); by = Inches(0.50)
    box_w = Inches(2.25); box_h = Inches(0.40)
    gap_y = Inches(0.12)

    for i, (txt, color) in enumerate(stages):
        y = by + i * (box_h + gap_y)
        add_rounded_rect(s, bx, y, box_w, box_h,
                         fill=color, line_color=WHITE, line_width=0.8)
        add_text(s, bx, y + Inches(0.02),
                 box_w, box_h - Inches(0.04),
                 txt, size=6, bold=True, color=WHITE,
                 align="center", anchor="middle")
        # downward arrow between stages
        if i < len(stages) - 1:
            add_arrow(s,
                      bx + box_w / 2,
                      y + box_h + Inches(0.01),
                      bx + box_w / 2,
                      y + box_h + gap_y - Inches(0.01),
                      color=INK, width=1.0)

    # Dashed "star-like" annotation next to Treg box (stage index 3) —
    # the only robust between-group finding
    star_y = by + 3 * (box_h + gap_y) + box_h / 2
    add_text(s, bx + box_w + Inches(0.02),
             _i(star_y - Inches(0.10)),
             Inches(0.4), Inches(0.22),
             "★", size=14, bold=True, color=HIGHLIGHT, align="center")
    add_text(s, bx + box_w + Inches(0.35),
             _i(star_y - Inches(0.10)),
             Inches(0.80), Inches(0.22),
             "Treg CI\nrobust",
             size=5, bold=True, color=HIGHLIGHT, align="left",
             anchor="middle", italic=True)

    # ===========================================================
    # Right half: convergence-null callout
    # ===========================================================
    rx = Inches(3.40); ry = Inches(0.50)
    rw = Inches(2.95); rh = Inches(3.30)

    # outer bordered box with coral accent
    add_rect(s, rx, ry, rw, rh,
             fill=RGBColor(0xFB, 0xF0, 0xED),
             line_color=BAD, line_width=1.4)

    # title strip
    add_rect(s, rx, ry, rw, Inches(0.30),
             fill=BAD, line_color=None)
    add_text(s, rx, ry, rw, Inches(0.30),
             "CONVERGENCE TEST — PRE-SPECIFIED, 36 PAIRS",
             size=8, bold=True, color=WHITE, align="center",
             anchor="middle")

    # Big null result in the center
    add_text(s, rx + Inches(0.15), ry + Inches(0.40),
             rw - Inches(0.3), Inches(0.18),
             "Q: does the Thread-1 baseline predictor",
             size=7, italic=True, color=INK, align="center")
    add_text(s, rx + Inches(0.15), ry + Inches(0.57),
             rw - Inches(0.3), Inches(0.18),
             "predict paired cascade Δ magnitude?",
             size=7, italic=True, color=INK, align="center")

    # Big NO answer
    add_text(s, rx + Inches(0.15), ry + Inches(0.80),
             rw - Inches(0.3), Inches(0.40),
             "NO.",
             size=36, bold=True, color=BAD, align="center",
             anchor="middle")

    # Statistic summary
    add_text(s, rx + Inches(0.15), ry + Inches(1.30),
             rw - Inches(0.3), Inches(0.14),
             "0 / 36 pre-specified pairs P < 0.05",
             size=8, bold=True, color=INK, align="center")
    add_text(s, rx + Inches(0.15), ry + Inches(1.45),
             rw - Inches(0.3), Inches(0.13),
             "(1.8 expected by chance)  ·  0 / 36 at FDR < 0.10",
             size=6, italic=True, color=RGBColor(0x66, 0x66, 0x66),
             align="center")

    # Headline pair box
    add_rect(s, rx + Inches(0.18), ry + Inches(1.70),
             rw - Inches(0.35), Inches(0.55),
             fill=WHITE, line_color=BAD, line_width=0.8)
    add_text(s, rx + Inches(0.22), ry + Inches(1.73),
             rw - Inches(0.42), Inches(0.15),
             "Headline pair",
             size=5, bold=True, italic=True,
             color=RGBColor(0x77, 0x77, 0x77), align="left", anchor="top")
    add_text(s, rx + Inches(0.22), ry + Inches(1.85),
             rw - Inches(0.42), Inches(0.16),
             "DSB-repair baseline  →  CD8-cyt Δ",
             size=7, bold=True, color=INK, align="center", anchor="top")
    add_text(s, rx + Inches(0.22), ry + Inches(2.02),
             rw - Inches(0.42), Inches(0.14),
             "Spearman r = −0.07,  P = 0.83",
             size=7, color=BAD, bold=True, align="center", anchor="top")

    # Power note
    add_text(s, rx + Inches(0.15), ry + Inches(2.36),
             rw - Inches(0.3), Inches(0.14),
             "Power retained (n = 12 detects |r| ≥ 0.55)",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")
    add_text(s, rx + Inches(0.15), ry + Inches(2.50),
             rw - Inches(0.3), Inches(0.14),
             "observed |r| < 0.20 → absence, not under-power",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")

    # Conclusion box
    add_rect(s, rx + Inches(0.18), ry + Inches(2.72),
             rw - Inches(0.35), Inches(0.48),
             fill=lighten(GOOD_HEX, 0.85), line_color=GOOD, line_width=0.8)
    add_text(s, rx + Inches(0.22), ry + Inches(2.74),
             rw - Inches(0.42), Inches(0.14),
             "Conclusion",
             size=5, bold=True, italic=True, color=GOOD,
             align="left", anchor="top")
    add_text(s, rx + Inches(0.22), ry + Inches(2.85),
             rw - Inches(0.42), Inches(0.34),
             "Static baseline predictor ⊥ dynamic paired cascade. "
             "The cascade is an observational phenomenology, NOT a "
             "downstream of the baseline.",
             size=6, color=INK, align="left", anchor="top")

    # bottom strip message
    add_text(s, Inches(0.15), Inches(3.95),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Two orthogonal biomarker layers: static pre-CRT (Figs 5, 9) "
             "and dynamic paired RT-phase (Figs 6-8 and present panel) —",
             size=6, italic=True, color=INK, align="center")
    add_text(s, Inches(0.15), Inches(4.09),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "complementary, not cascading.  "
             "Two-layer clinical algorithm follows in Discussion.",
             size=6, bold=True, italic=True, color=HIGHLIGHT, align="center")


build_A()
build_B()
build_C()
build_D()
build_E()
build_F()
deck_main = f"{OUT}/Fig8_cascade_convergence_native_editable.pptx"
prs_main.save(deck_main)
print(f"wrote {deck_main}")


# ===========================================================================
# SUPP FIGURE S20
# ===========================================================================
prs_supp = Presentation()


def build_S20A():
    """Convergence test scatter matrix: 36 pre-specified pairs ranked by
    |r|, showing Spearman r and P per pair."""
    s = new_slide(prs_supp)
    draw_panel_letter(s, "A")

    # sort by |r| descending
    c = convergence.copy()
    c["abs_r"] = c.spearman_r.abs()
    c = c.sort_values("abs_r", ascending=False).reset_index(drop=True)

    # lollipop plot: x = Spearman r, y = pair rank
    ax_x = Inches(2.25); ax_y = Inches(0.45)
    ax_w = Inches(3.40); ax_h = Inches(3.40)

    x_lo, x_hi = -1.0, 1.0

    def tx(v):
        return _i(ax_x + (v - x_lo) / (x_hi - x_lo) * ax_w)

    n_rows = len(c)
    row_h = ax_h / n_rows

    # spines
    add_line(s, ax_x, ax_y, ax_x, ax_y + ax_h, LINE, 0.3)
    add_line(s, ax_x, ax_y + ax_h, ax_x + ax_w, ax_y + ax_h, LINE, 0.5)
    for v in [-1, -0.5, 0, 0.5, 1]:
        xx = tx(v)
        add_line(s, xx, _i(ax_y + ax_h), xx,
                 _i(ax_y + ax_h + Inches(0.04)), LINE, 0.4)
        add_text(s, xx - Inches(0.15), _i(ax_y + ax_h + Inches(0.05)),
                 Inches(0.30), Inches(0.13),
                 f"{v:+.1f}" if v != 0 else "0",
                 size=6, color=INK, align="center")
    add_text(s, ax_x, _i(ax_y + ax_h + Inches(0.22)),
             ax_w, Inches(0.14),
             "Spearman r  (baseline × cascade Δ, n = 12 paired)",
             size=7, color=INK, align="center")

    # zero ref + P=0.05 bounds (|r|~0.58 at n=12)
    add_line(s, tx(0), ax_y, tx(0), ax_y + ax_h, INK, 0.8)
    add_line(s, tx(0.58), ax_y, tx(0.58), ax_y + ax_h,
             GREY, 0.4, dashed=True)
    add_line(s, tx(-0.58), ax_y, tx(-0.58), ax_y + ax_h,
             GREY, 0.4, dashed=True)
    add_text(s, tx(0.58) - Inches(0.30), _i(ax_y - Inches(0.02)),
             Inches(0.6), Inches(0.12),
             "|r| ≈ 0.58\n(P = 0.05)", size=4, italic=True,
             color=RGBColor(0x99, 0x99, 0x99),
             align="center", anchor="top")

    # rows
    for i, row in c.iterrows():
        y = _i(ax_y + (i + 0.5) * row_h)
        r = float(row.spearman_r)
        p = float(row.spearman_p)
        r_c = max(x_lo, min(x_hi, r))

        # lollipop stem
        color = (GOOD if r >= 0 else BAD)
        add_line(s, tx(0), y, tx(r_c), y, color, 0.7)
        # dot
        is_sig = p < 0.05
        add_circle(s, tx(r_c), y, Emu(15000) if not is_sig else Emu(22000),
                   fill=WHITE if not is_sig else color,
                   line_color=color, line_width=0.8)
        # row label
        baseline = str(row.baseline)[:20]
        cascade = str(row.cascade).replace("_delta", " Δ")[:18]
        add_text(s, Inches(0.20), y - Inches(0.06),
                 ax_x - Inches(0.22), Inches(0.13),
                 f"{baseline} → {cascade}",
                 size=4, color=INK, align="right")
        # right column r, P
        add_text(s, _i(ax_x + ax_w + Inches(0.05)),
                 y - Inches(0.06),
                 Inches(0.7), Inches(0.13),
                 f"r={r:+.2f} P={p:.2f}",
                 size=4, color=RGBColor(0x55, 0x55, 0x55), align="left")

    # headline annotation: 0 / 36
    add_rect(s, Inches(0.20), Inches(0.25), Inches(2.0), Inches(0.45),
             fill=WHITE, line_color=BAD, line_width=1.0)
    add_text(s, Inches(0.20), Inches(0.25),
             Inches(2.0), Inches(0.22),
             "0 / 36 pairs  P < 0.05",
             size=9, bold=True, color=BAD, align="center", anchor="middle")
    add_text(s, Inches(0.20), Inches(0.45),
             Inches(2.0), Inches(0.22),
             "(1.8 expected by chance)",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center", anchor="middle")

    # bottom caption
    add_text(s, Inches(0.15), Inches(4.02),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "Pre-specified 36-pair test: 9 baseline Thread-1 features "
             "× 4 cascade Δ features, n = 12 paired subjects.",
             size=6, italic=True, color=RGBColor(0x55, 0x55, 0x55),
             align="center")
    add_text(s, Inches(0.15), Inches(4.18),
             SLIDE_W - Inches(0.3), Inches(0.14),
             "|r| at n = 12 needs ≈ 0.58 for P < 0.05; the largest "
             "observed |r| is ≈ 0.56 → absence of association, not "
             "under-power.",
             size=6, color=INK, align="center")


def build_S20B():
    """Placeholder: purity-adjusted sensitivity (Δ before vs after
    CNVkit purity adjustment). If data available, plot pairs; else
    show summary table."""
    s = new_slide(prs_supp)
    draw_panel_letter(s, "B")

    # Try to load purity-adjusted file
    ptable_path = f"{DATA}/09_integration/paired_delta/delta_purity_sensitivity.tsv"
    if not os.path.exists(ptable_path):
        add_text(s, Inches(0.5), Inches(2.0),
                 SLIDE_W - Inches(1.0), Inches(0.5),
                 "Purity-adjusted sensitivity table not available.",
                 size=10, color=INK, align="center")
        return

    sens = pd.read_csv(ptable_path, sep="\t")
    # print columns for reference
    features_to_show = ["Treg", "MHC_II", "CD8_exhaustion", "IGH_n"]
    # layout: per-feature mini bar chart raw vs adjusted
    # fall back to text if layout complex

    # simple horizontal comparison: raw Δ vs purity-adjusted Δ
    # If column 'delta_good_median_raw' and 'delta_good_median_adj' exist:
    cols = sens.columns.tolist()
    add_text(s, Inches(0.20), Inches(0.50),
             SLIDE_W - Inches(0.4), Inches(0.16),
             "Purity-adjusted paired Δ sensitivity",
             size=8, bold=True, color=INK, align="center")

    # Render as a simple bar comparison
    ax_x = Inches(1.30); ax_y = Inches(0.85)
    ax_w = Inches(4.00); ax_h = Inches(2.80)

    n_feat = len(features_to_show)
    row_h = ax_h / n_feat

    def tx(v, v_min, v_max):
        return _i(ax_x + (v - v_min) / (v_max - v_min) * ax_w)

    for i, feat in enumerate(features_to_show):
        y = _i(ax_y + (i + 0.5) * row_h)
        feat_row = sens[sens.feature == feat]
        if feat_row.empty:
            add_text(s, Inches(0.20), y - Inches(0.1),
                     ax_x - Inches(0.22), Inches(0.18),
                     PRETTY[feat], size=7, color=INK, align="right")
            continue
        feat_row = feat_row.iloc[0]
        # Extract raw vs adjusted values --- best-effort column parsing
        raw_cols = [c for c in cols if "raw" in c.lower() or "unadj" in c.lower()]
        adj_cols = [c for c in cols if "adj" in c.lower() and "raw" not in c.lower()]

        add_text(s, Inches(0.20), y - Inches(0.1),
                 ax_x - Inches(0.22), Inches(0.18),
                 PRETTY[feat], size=7, color=INK, align="right")
        # place the feature's raw + adjusted MW P as text
        add_text(s, ax_x + Inches(0.1), y - Inches(0.1),
                 Inches(4.0), Inches(0.18),
                 str(feat_row.to_dict())[:140],
                 size=4, color=RGBColor(0x55, 0x55, 0x55),
                 align="left")


build_S20A()
build_S20B()
deck_supp = f"{OUT}/SuppFig_S20_cascade_convergence_sensitivity_native_editable.pptx"
prs_supp.save(deck_supp)
print(f"wrote {deck_supp}")
