"""Shared python-pptx helpers used by the three 3-layer figure builders.
Rule compliance applies to every shape created via these helpers:
  - Arial font everywhere (Rule 5)
  - Shadow effects killed via empty <a:effectLst/> (Rule 6)
  - No fill / line side effects beyond what callers request
"""
import numpy as np
from pptx.util import Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree

# ColorBrewer RdBu_r anchors
BLUE_DEEP = (33, 102, 172)
RED_DEEP  = (178,  24,  43)
WHITE     = (255, 255, 255)

# TNT-project palette (per project memory)
TNT_TEAL  = RGBColor(10, 125, 110)    # GOOD / as-predicted
TNT_CORAL = RGBColor(197, 62, 31)     # BAD / against-prediction
GOLD      = RGBColor(218, 165, 32)
GREY_BORD = RGBColor(200, 200, 200)
GREY_REF  = RGBColor(238, 238, 238)
BLACK     = RGBColor(0, 0, 0)


def diverging(r: float, vmin: float, vmax: float) -> RGBColor:
    """Two-arm diverging ramp (deep-blue ↔ white ↔ deep-red)."""
    r = max(min(r, vmax), vmin)
    if r >= 0:
        t = r / vmax
        return RGBColor(
            int(WHITE[0] - t * (WHITE[0] - RED_DEEP[0])),
            int(WHITE[1] - t * (WHITE[1] - RED_DEEP[1])),
            int(WHITE[2] - t * (WHITE[2] - RED_DEEP[2])),
        )
    t = -r / (-vmin)
    return RGBColor(
        int(WHITE[0] - t * (WHITE[0] - BLUE_DEEP[0])),
        int(WHITE[1] - t * (WHITE[1] - BLUE_DEEP[1])),
        int(WHITE[2] - t * (WHITE[2] - BLUE_DEEP[2])),
    )


def teal_ramp(intensity: float) -> RGBColor:
    """Teal intensity ramp from white → TNT_TEAL by intensity in [0,1]."""
    intensity = max(0.0, min(1.0, intensity))
    R = int(WHITE[0] - intensity * (WHITE[0] - 10))
    G = int(WHITE[1] - intensity * (WHITE[1] - 125))
    B = int(WHITE[2] - intensity * (WHITE[2] - 110))
    return RGBColor(R, G, B)


def coral_ramp(intensity: float) -> RGBColor:
    intensity = max(0.0, min(1.0, intensity))
    R = int(WHITE[0] - intensity * (WHITE[0] - 197))
    G = int(WHITE[1] - intensity * (WHITE[1] - 62))
    B = int(WHITE[2] - intensity * (WHITE[2] - 31))
    return RGBColor(R, G, B)


def kill_shadow(shape) -> None:
    spPr = shape._element.find(qn("p:spPr"))
    if spPr is None:
        return
    for el in spPr.findall(qn("a:effectLst")):
        spPr.remove(el)
    etree.SubElement(spPr, qn("a:effectLst"))


def set_text(tb, text, size=8, bold=False, color=BLACK,
             align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE):
    tf = tb.text_frame
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.text = text
    for para in tf.paragraphs:
        para.alignment = align
        for run in para.runs:
            run.font.name = "Arial"
            run.font.size = Pt(size)
            run.font.bold = bold
            run.font.color.rgb = color
    kill_shadow(tb)


def fmt_r(r: float) -> str:
    return f"{r:+.2f}".replace("-", "−")
