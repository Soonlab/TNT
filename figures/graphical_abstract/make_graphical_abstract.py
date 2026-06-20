#!/usr/bin/env python3
"""Graphical abstract for TNT manuscript (BioRender-style, vector SVG master)."""
from pathlib import Path
import cairosvg

OUT = Path("/data/data/TNT/analysis/figures/graphical_abstract")
OUT.mkdir(parents=True, exist_ok=True)

# ---------- palette ----------
C = {
    "good":   "#2E86AB",
    "bad":    "#E63946",
    "amber":  "#F4A261",
    "sage":   "#52B788",
    "slate":  "#2A2E3A",
    "ink":    "#1F2430",
    "muted":  "#6B7280",
    "line":   "#C7CEDB",
    "bg":     "#FFFFFF",
    "soft":   "#F4F6F8",
    "softer": "#EEF2F6",
    "tealbg": "#E6F1F6",
    "redbg":  "#FCE7E9",
    "amberbg":"#FDEFDD",
    "sagebg": "#E3F2E9",
}

W, H = 1600, 940
FONT = 'Inter, "Helvetica Neue", Arial, sans-serif'

# ---------- helpers ----------
def rect(x, y, w, h, fill, stroke=None, sw=1.4, rx=14, op=1.0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" ry="{rx}" fill="{fill}" fill-opacity="{op}"{s}/>'

def _xml_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

def txt(x, y, s, size=14, color=None, weight=400, anchor="start"):
    col = color or C["ink"]
    s = _xml_escape(str(s))
    return f'<text x="{x}" y="{y}" font-family=\'{FONT}\' font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="{anchor}">{s}</text>'

def arrow(x1, y1, x2, y2, color=None, sw=2.4, dash=None):
    col = color or C["slate"]
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" '
            f'stroke-width="{sw}" stroke-linecap="round" marker-end="url(#arrow-{col[1:]})"{d}/>')

def arrow_marker(color):
    cid = f"arrow-{color[1:]}"
    return (f'<marker id="{cid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')

# ---------- icon builders ----------
def person(cx, cy, scale=1.0, color="#5B6C8A"):
    s = scale
    head_r = 6*s
    body = (f'<circle cx="{cx}" cy="{cy-9*s}" r="{head_r}" fill="{color}"/>'
            f'<path d="M {cx-9*s} {cy+9*s} q 0 -12 {9*s} -12 q {9*s} 0 {9*s} 12 Z" fill="{color}"/>')
    return body

def tube(cx, cy, color, label=None):
    # rounded tube
    w, h = 18, 46
    body = (f'<rect x="{cx-w/2}" y="{cy-h/2}" width="{w}" height="{h}" rx="9" ry="9" '
            f'fill="{color}" stroke="{C["slate"]}" stroke-width="1.2"/>'
            f'<rect x="{cx-w/2+2}" y="{cy-h/2+4}" width="{w-4}" height="6" rx="3" '
            f'fill="#FFFFFF" fill-opacity="0.55"/>')
    lab = txt(cx, cy+h/2+14, label, size=11, color=C["muted"], anchor="middle") if label else ""
    return body + lab

def dna_helix(cx, cy, w=56, h=30, color=None):
    col = color or C["good"]
    # simplified double helix via two sinusoidal paths + rungs
    import math
    pts_a, pts_b = [], []
    n = 40
    for i in range(n+1):
        t = i/n
        x = cx - w/2 + t*w
        y1 = cy + math.sin(t*math.pi*3)*h/2
        y2 = cy - math.sin(t*math.pi*3)*h/2
        pts_a.append(f"{x:.1f},{y1:.1f}")
        pts_b.append(f"{x:.1f},{y2:.1f}")
    rungs = ""
    for i in range(2, n, 4):
        t = i/n
        x = cx - w/2 + t*w
        y1 = cy + math.sin(t*math.pi*3)*h/2
        y2 = cy - math.sin(t*math.pi*3)*h/2
        rungs += f'<line x1="{x:.1f}" y1="{y1:.1f}" x2="{x:.1f}" y2="{y2:.1f}" stroke="{C["muted"]}" stroke-width="1"/>'
    a = f'<polyline points="{" ".join(pts_a)}" fill="none" stroke="{col}" stroke-width="2.2" stroke-linecap="round"/>'
    b = f'<polyline points="{" ".join(pts_b)}" fill="none" stroke="{C["amber"]}" stroke-width="2.2" stroke-linecap="round"/>'
    return rungs + a + b

def mrna_strand(cx, cy, w=56, color=None):
    col = color or C["sage"]
    import math
    pts = []
    for i in range(41):
        t=i/40
        x = cx-w/2 + t*w
        y = cy + math.sin(t*math.pi*4)*6
        pts.append(f"{x:.1f},{y:.1f}")
    return f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="2.4" stroke-linecap="round"/>'

def cell_cluster(cx, cy, n=6, color=None, r=8):
    col = color or C["good"]
    import math, random
    random.seed(int(cx+cy))
    out=""
    for i in range(n):
        ang = 2*math.pi*i/n + 0.2
        rad = 14 + (i%2)*4
        x = cx + math.cos(ang)*rad
        y = cy + math.sin(ang)*rad
        out += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}" fill-opacity="0.85" stroke="#FFFFFF" stroke-width="1.5"/>'
        # tiny nucleus
        out += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r*0.4:.1f}" fill="{C["slate"]}" fill-opacity="0.35"/>'
    return out

def biopsy_needle(cx, cy, color=None, angle=-20):
    col = color or C["slate"]
    return (f'<g transform="translate({cx},{cy}) rotate({angle})">'
            f'<rect x="-24" y="-3" width="34" height="6" rx="2" fill="{col}"/>'
            f'<rect x="10" y="-1.5" width="22" height="3" rx="1.5" fill="{C["muted"]}"/>'
            f'<polygon points="32,-2 40,0 32,2" fill="{col}"/>'
            f'</g>')

# ---------- layout ----------
parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')

# defs: arrow markers in multiple colors, gradients
defs = ['<defs>']
for col in [C["slate"], C["good"], C["bad"], C["amber"], C["sage"], C["muted"]]:
    defs.append(arrow_marker(col))
defs.append(f'<linearGradient id="goodGrad" x1="0" y1="0" x2="0" y2="1">'
            f'<stop offset="0%" stop-color="{C["tealbg"]}"/>'
            f'<stop offset="100%" stop-color="#FFFFFF"/></linearGradient>')
defs.append(f'<linearGradient id="badGrad" x1="0" y1="0" x2="0" y2="1">'
            f'<stop offset="0%" stop-color="{C["redbg"]}"/>'
            f'<stop offset="100%" stop-color="#FFFFFF"/></linearGradient>')
defs.append(f'<filter id="soft" x="-10%" y="-10%" width="120%" height="120%">'
            f'<feGaussianBlur stdDeviation="3"/></filter>')
defs.append(f'<filter id="shadow" x="-10%" y="-10%" width="120%" height="130%">'
            f'<feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="#1F2430" flood-opacity="0.12"/></filter>')
defs.append('</defs>')
parts.extend(defs)

# background
parts.append(rect(0,0,W,H, C["bg"], rx=0))

# --- title ---
parts.append(txt(60, 58, "Radiation-phase molecular response predicts final TNT outcome in MSS LARC",
                 size=28, color=C["ink"], weight=700))
parts.append(txt(60, 88,
    "Integrated WES + RNA-seq of 35 locally advanced rectal cancers (pre / post CRT biopsies) — Kim et al., 2026",
    size=14, color=C["muted"], weight=400))

# thin divider
parts.append(f'<line x1="60" y1="108" x2="{W-60}" y2="108" stroke="{C["line"]}" stroke-width="1"/>')

# =========================================================
# Panel A — Cohort & sampling (left)
# =========================================================
AX, AY, AW, AH = 60, 130, 430, 380
parts.append(rect(AX, AY, AW, AH, C["soft"], stroke=C["line"], sw=1, rx=18))
parts.append(txt(AX+20, AY+30, "A · Cohort & radiation-phase sampling",
                 size=15, color=C["slate"], weight=700))
parts.append(txt(AX+20, AY+50, "35 MSS LARC · long-course CRT (50.4 Gy + capecitabine)",
                 size=11.5, color=C["muted"]))

# patient grid (7 × 5 = 35)
pg_x0, pg_y0 = AX+32, AY+82
for i in range(35):
    r, c = divmod(i, 7)
    cx = pg_x0 + c*28
    cy = pg_y0 + r*32
    good = i < 18  # 18 good, 17 bad
    col = C["good"] if good else C["bad"]
    parts.append(person(cx, cy, scale=1.05, color=col))
parts.append(txt(AX+20, AY+240, "● 18 good (CR / near-CR)", size=11.5, color=C["good"], weight=600))
parts.append(txt(AX+20, AY+258, "● 17 bad (PR / poor)",     size=11.5, color=C["bad"],  weight=600))

# biopsy timeline
ty = AY+300
parts.append(f'<line x1="{AX+24}" y1="{ty}" x2="{AX+AW-24}" y2="{ty}" stroke="{C["slate"]}" stroke-width="2"/>')
# radiation phase band
parts.append(rect(AX+120, ty-10, 200, 20, C["amberbg"], stroke=C["amber"], sw=1.2, rx=6, op=0.9))
parts.append(txt(AX+220, ty+5, "Radiation phase (~5–6 wk)", size=11, color=C["amber"], weight=600, anchor="middle"))
# timepoint dots
for tx, lab, col in [(AX+100, "pre", C["good"]), (AX+340, "post", C["good"]), (AX+400, "final\nresponse", C["slate"])]:
    parts.append(f'<circle cx="{tx}" cy="{ty}" r="7" fill="{col}" stroke="#FFFFFF" stroke-width="2"/>')
    if "\n" in lab:
        for i,line in enumerate(lab.split("\n")):
            parts.append(txt(tx, ty+24+i*13, line, size=10.5, color=C["muted"], anchor="middle"))
    else:
        parts.append(txt(tx, ty+24, lab, size=11, color=C["slate"], anchor="middle", weight=600))
# biopsy needles
parts.append(biopsy_needle(AX+100, ty-28, color=C["slate"]))
parts.append(biopsy_needle(AX+340, ty-28, color=C["slate"]))

# assay counts
ay = AY+AH-30
parts.append(txt(AX+20, ay, "WES 77  ·  RNA-seq 56  ·  49 PASS somatic VCFs  ·  all MSS, TMB < 2/Mb",
                 size=11, color=C["muted"]))

# =========================================================
# Panel B — Multi-omics & analysis (center top)
# =========================================================
BX, BY, BW, BH = 520, 130, 530, 170
parts.append(rect(BX, BY, BW, BH, C["soft"], stroke=C["line"], sw=1, rx=18))
parts.append(txt(BX+20, BY+28, "B · Multi-omics profiling", size=15, color=C["slate"], weight=700))

# WES card
parts.append(rect(BX+24, BY+50, 230, 100, "url(#goodGrad)", stroke=C["good"], sw=1.2, rx=12))
parts.append(dna_helix(BX+24+60, BY+90, w=100, h=30, color=C["good"]))
parts.append(txt(BX+24+16, BY+135, "WES · Mutect2 T–N", size=12, color=C["good"], weight=700))
parts.append(txt(BX+24+16, BY+150, "TMB · SBS · MSI · CNV · HLA", size=10.5, color=C["muted"]))

# RNA card
parts.append(rect(BX+280, BY+50, 230, 100, "url(#goodGrad)", stroke=C["sage"], sw=1.2, rx=12))
parts.append(mrna_strand(BX+280+70, BY+90, w=110, color=C["sage"]))
parts.append(txt(BX+280+16, BY+135, "RNA-seq · StringTie / DESeq2", size=12, color=C["sage"], weight=700))
parts.append(txt(BX+280+16, BY+150, "GSEA · ssGSEA · deconvolution · TRUST4", size=10.5, color=C["muted"]))

# =========================================================
# Panel C — Good-responder cascade (center, vertical)
# =========================================================
CX, CY, CW, CH = 520, 320, 530, 490
parts.append(rect(CX, CY, CW, CH, "#FFFFFF", stroke=C["line"], sw=1, rx=18))
parts.append(txt(CX+20, CY+30, "C · Good-responder cascade (pre → post)",
                 size=15, color=C["slate"], weight=700))
parts.append(txt(CX+20, CY+50, "radiation phase induces sequential tumor clearance & immune reprogramming",
                 size=11.5, color=C["muted"]))

# 5 vertical steps
steps = [
    ("1",  "Mutation clearance",    "SBS5 Δ −76 · missense Δ −67",   C["good"]),
    ("2",  "Neoantigen depletion",  "MHC-I binders Δ −312 · HLA diversity ↓", C["amber"]),
    ("3",  "HLA-LOH clone removal", "subj 3 · 4 pre → post resolution",  C["bad"]),
    ("4",  "Immune reprogramming",  "Treg ↑ P=0.026 · MHC-II ↑ · CD8 exhaustion ↑", C["sage"]),
    ("5",  "B-cell infiltration",   "IGH / IGK / IGL clone counts ↑ ≈+1400", C["good"]),
]
step_y = CY+80
for i,(n,title,sub,col) in enumerate(steps):
    y = step_y + i*78
    # number disk
    parts.append(f'<circle cx="{CX+46}" cy="{y+26}" r="20" fill="{col}" filter="url(#shadow)"/>')
    parts.append(txt(CX+46, y+32, n, size=18, color="#FFFFFF", weight=700, anchor="middle"))
    # card
    parts.append(rect(CX+80, y, CW-110, 54, C["softer"], stroke=col, sw=1.3, rx=10))
    parts.append(txt(CX+100, y+22, title, size=13.5, color=C["ink"], weight=700))
    parts.append(txt(CX+100, y+42, sub, size=11, color=C["muted"]))
    if i < len(steps)-1:
        parts.append(arrow(CX+46, y+48, CX+46, y+72, color=col, sw=2.8))

# =========================================================
# Panel D — External validation (right top)
# =========================================================
DX, DY, DW, DH = 1080, 130, 460, 340
parts.append(rect(DX, DY, DW, DH, C["soft"], stroke=C["line"], sw=1, rx=18))
parts.append(txt(DX+20, DY+30, "D · External CD8-cytotoxic axis validation",
                 size=15, color=C["slate"], weight=700))
parts.append(txt(DX+20, DY+50, "9 independent nCRT rectal cohorts · N = 721 + Akiyoshi 2023 (N = 298)",
                 size=11.5, color=C["muted"]))

# mini forest
fx0 = DX+95; fx1 = DX+DW-30
zero_x = fx0 + (fx1-fx0)*0.25
cohorts = [
    ("GSE45404", +0.55, 0.30),
    ("GSE46862", +0.70, 0.35),
    ("GSE87211", +0.40, 0.28),
    ("GSE94104", +0.62, 0.33),
    ("GSE133057",+0.35, 0.30),
    ("GSE35452", +0.58, 0.32),
    ("GSE119409",+0.25, 0.34),
    ("GSE150082",+0.32, 0.30),
    ("GSE56699", -0.10, 0.30),
    ("Akiyoshi'23", +0.75, 0.22),
]
fy0 = DY+90
row_h = 18
# axis
parts.append(f'<line x1="{zero_x}" y1="{fy0-6}" x2="{zero_x}" y2="{fy0+len(cohorts)*row_h+6}" stroke="{C["muted"]}" stroke-width="1" stroke-dasharray="3,3"/>')
scale = 140  # px per unit effect
for i,(name, eff, se) in enumerate(cohorts):
    yy = fy0 + i*row_h + 10
    cx = zero_x + eff*scale
    lo = zero_x + (eff-1.96*se)*scale
    hi = zero_x + (eff+1.96*se)*scale
    col = C["good"] if eff>0 else C["bad"]
    parts.append(txt(fx0-4, yy+3, name, size=10, color=C["slate"], anchor="end"))
    parts.append(f'<line x1="{lo}" y1="{yy}" x2="{hi}" y2="{yy}" stroke="{col}" stroke-width="1.6" stroke-linecap="round"/>')
    parts.append(f'<rect x="{cx-4}" y="{yy-4}" width="8" height="8" fill="{col}" stroke="#FFFFFF" stroke-width="1"/>')

# diamond — meta
meta_eff = 0.52
dy_m = fy0 + len(cohorts)*row_h + 22
cx = zero_x + meta_eff*scale
parts.append(f'<polygon points="{cx-40},{dy_m} {cx},{dy_m-8} {cx+40},{dy_m} {cx},{dy_m+8}" fill="{C["slate"]}"/>')
parts.append(txt(fx0-4, dy_m+3, "Meta (9)", size=10.5, color=C["slate"], anchor="end", weight=700))
parts.append(txt(cx+46, dy_m+3, "Z = +2.74   P = 0.006", size=10.5, color=C["slate"], weight=600))

# x-axis labels
axy = dy_m + 22
parts.append(txt(zero_x, axy, "0", size=10, color=C["muted"], anchor="middle"))
parts.append(txt(zero_x-scale*0.5, axy, "−0.5", size=10, color=C["muted"], anchor="middle"))
parts.append(txt(zero_x+scale*0.5, axy, "+0.5", size=10, color=C["muted"], anchor="middle"))
parts.append(txt(zero_x+scale*1.0, axy, "+1.0", size=10, color=C["muted"], anchor="middle"))
parts.append(txt(DX+DW/2, axy+18, "Standardized effect (good vs bad)", size=10.5, color=C["muted"], anchor="middle"))

# =========================================================
# Panel E — Clinical implication (right middle)
# =========================================================
EX, EY, EW, EH = 1080, 490, 460, 320
parts.append(rect(EX, EY, EW, EH, "url(#goodGrad)", stroke=C["good"], sw=1.2, rx=18))
parts.append(txt(EX+20, EY+30, "E · Mid-treatment decision window",
                 size=15, color=C["good"], weight=700))
parts.append(txt(EX+20, EY+52, "actionable use of radiation-phase molecular response",
                 size=11.5, color=C["muted"]))

# decision tree
ty = EY+90
parts.append(f'<rect x="{EX+20}" y="{ty}" width="{EW-40}" height="46" rx="10" fill="#FFFFFF" stroke="{C["good"]}" stroke-width="1.4"/>')
parts.append(txt(EX+EW/2, ty+20, "Post-CRT biopsy molecular profile", size=12.5, color=C["ink"], weight=700, anchor="middle"))
parts.append(txt(EX+EW/2, ty+36, "CD8-cytotoxic · Treg · neoantigen Δ · mutation clearance", size=10.5, color=C["muted"], anchor="middle"))

# branches
bry = ty+46
parts.append(arrow(EX+EW/2-60, bry, EX+EW/2-120, bry+40, color=C["good"], sw=2.4))
parts.append(arrow(EX+EW/2+60, bry, EX+EW/2+120, bry+40, color=C["bad"], sw=2.4))

# good branch
parts.append(rect(EX+30, bry+42, 180, 120, "#FFFFFF", stroke=C["good"], sw=1.4, rx=10))
parts.append(txt(EX+120, bry+66, "Responder profile", size=12.5, color=C["good"], weight=700, anchor="middle"))
parts.append(txt(EX+120, bry+86,  "→ de-escalate", size=11.5, color=C["ink"], anchor="middle"))
parts.append(txt(EX+120, bry+104, "watch-and-wait",   size=11,   color=C["muted"], anchor="middle"))
parts.append(txt(EX+120, bry+122, "or organ preservation", size=11, color=C["muted"], anchor="middle"))

# bad branch
parts.append(rect(EX+250, bry+42, 180, 120, "#FFFFFF", stroke=C["bad"], sw=1.4, rx=10))
parts.append(txt(EX+340, bry+66, "Non-responder profile", size=12.5, color=C["bad"], weight=700, anchor="middle"))
parts.append(txt(EX+340, bry+86,  "→ intensify", size=11.5, color=C["ink"], anchor="middle"))
parts.append(txt(EX+340, bry+104, "consolidation FOLFOX", size=11, color=C["muted"], anchor="middle"))
parts.append(txt(EX+340, bry+122, "or trial enrolment", size=11, color=C["muted"], anchor="middle"))

# =========================================================
# Flow arrows between panels
# =========================================================
# A -> B
parts.append(arrow(AX+AW, AY+AH/2-40, BX-6, BY+BH/2, color=C["muted"], sw=3))
# B -> C
parts.append(arrow(BX+BW/2, BY+BH+4, BX+BW/2, CY-6, color=C["muted"], sw=3))
# C -> D (meta)
parts.append(arrow(CX+CW, CY+80, DX-6, DY+DH/2-20, color=C["muted"], sw=3))
# C -> E (clinical)
parts.append(arrow(CX+CW, CY+CH-80, EX-6, EY+40, color=C["muted"], sw=3))

# =========================================================
# Bottom take-home bar
# =========================================================
BTY = 840
parts.append(rect(60, BTY, W-120, 72, C["slate"], rx=14))
parts.append(txt(80, BTY+30,
    "Take-home:  Radiation-phase molecular response — mutation / neoantigen clearance and a reproducible CD8-cytotoxic axis —",
    size=14, color="#FFFFFF", weight=600))
parts.append(txt(80, BTY+52,
    "defines a mid-treatment biomarker window to personalise consolidation chemotherapy in MSS locally advanced rectal cancer.",
    size=14, color="#E6F1F6", weight=400))

parts.append("</svg>")

svg = "\n".join(parts)
svg_path = OUT/"GA_master.svg"
svg_path.write_text(svg, encoding="utf-8")

# render PNG (2x) and PDF
cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                 write_to=str(OUT/"GA.png"),
                 output_width=W*2, output_height=H*2)
cairosvg.svg2pdf(bytestring=svg.encode("utf-8"),
                 write_to=str(OUT/"GA.pdf"))

print("wrote:", svg_path)
print("wrote:", OUT/"GA.png")
print("wrote:", OUT/"GA.pdf")
