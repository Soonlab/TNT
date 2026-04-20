#!/usr/bin/env python3
"""Graphical abstract v3 — Nature/Cell-level professional design.

Design principles applied:
- Okabe-Ito colorblind-friendly palette
- Left-to-right 3-section flow (Design → Mechanism → Validation)
- Central biological illustration as focal point
- Minimal text (<80 words), no drop shadows, no decorative clutter
- Consistent line weights, high contrast, ample whitespace
- Clean sans-serif typography
"""
import math
from pathlib import Path
try:
    import cairosvg
    HAVE_CAIRO = True
except ImportError:
    HAVE_CAIRO = False

OUT = Path("/data/data/TNT/analysis/figures/graphical_abstract")

# ---- Okabe-Ito + Nature-style palette ----
C = {
    "blue":    "#0072B2",  # Okabe-Ito blue (good/responder)
    "orange":  "#D55E00",  # Okabe-Ito vermillion (bad/non-responder)
    "sky":     "#56B4E9",  # Okabe-Ito sky blue (accent)
    "green":   "#009E73",  # Okabe-Ito bluish green (immune/RNA)
    "yellow":  "#F0E442",  # Okabe-Ito yellow
    "purple":  "#CC79A7",  # Okabe-Ito reddish purple (B-cell)
    "dkgray":  "#2B2D42",  # near-black (text, outlines)
    "gray":    "#6C757D",  # mid-gray (secondary text)
    "ltgray":  "#E9ECEF",  # light gray (card fills)
    "bgcard":  "#F8F9FA",  # very light card bg
    "white":   "#FFFFFF",
    "blueBg":  "#E8F4FA",
    "orangeBg":"#FDECE6",
    "greenBg": "#E5F5EF",
    "purpleBg":"#F5EDF2",
}

W, H = 1800, 820
FONT = 'Inter, Arial, Helvetica, sans-serif'
LW = 1.6  # standard line weight

def _esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

def rect(x,y,w,h,fill,stroke=None,sw=LW,rx=12,op=1.0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" ry="{rx}" fill="{fill}" fill-opacity="{op}"{s}/>'

def txt(x,y,s,size=13,color=None,weight=400,anchor="start",ls=0):
    col = color or C["dkgray"]; s = _esc(s)
    let = f' letter-spacing="{ls}"' if ls else ""
    return f'<text x="{x}" y="{y}" font-family=\'{FONT}\' font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="{anchor}"{let}>{s}</text>'

def arrow_marker(color, mid=""):
    cid = f"arr{mid}{color[1:]}"
    return (f'<marker id="{cid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
            f'<path d="M0,1 L10,5 L0,9 z" fill="{color}"/></marker>')

def arrow_line(x1,y1,x2,y2,color=None,sw=2.2):
    col = color or C["dkgray"]
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" '
            f'stroke-width="{sw}" stroke-linecap="round" marker-end="url(#arr{col[1:]})"/>')

def curved_arrow(x1,y1,x2,y2,color=None,sw=2.2,bend=30):
    col = color or C["dkgray"]
    mx=(x1+x2)/2; my=(y1+y2)/2; dx=x2-x1; dy=y2-y1; L=math.hypot(dx,dy)
    if L==0: return ""
    cpx=mx+(-dy/L)*bend; cpy=my+(dx/L)*bend
    return (f'<path d="M{x1},{y1} Q{cpx:.1f},{cpy:.1f} {x2},{y2}" fill="none" '
            f'stroke="{col}" stroke-width="{sw}" stroke-linecap="round" marker-end="url(#arr{col[1:]})"/>')

# ---- biological icons (clean, uniform line weight) ----

def icon_colon_section(cx,cy,s=1.0):
    """Clean colon cross-section with tumor."""
    o=""
    # outer wall
    o+=f'<ellipse cx="{cx}" cy="{cy}" rx="{34*s}" ry="{28*s}" fill="#F8D0D0" stroke="#C07070" stroke-width="{LW*s}"/>'
    # mucosal folds (3 crescents)
    for ang in [40,160,280]:
        a=math.radians(ang)
        fx=cx+18*s*math.cos(a); fy=cy+14*s*math.sin(a)
        o+=f'<ellipse cx="{fx:.1f}" cy="{fy:.1f}" rx="{8*s}" ry="{5*s}" fill="#F0B0B0" fill-opacity="0.5" transform="rotate({ang},{fx:.1f},{fy:.1f})"/>'
    # lumen
    o+=f'<ellipse cx="{cx}" cy="{cy}" rx="{12*s}" ry="{9*s}" fill="#FCF0EE"/>'
    # tumor mass
    tx=cx+14*s; ty=cy-6*s
    o+=f'<ellipse cx="{tx}" cy="{ty}" rx="{14*s}" ry="{11*s}" fill="{C["orange"]}" fill-opacity="0.8" stroke="{C["orange"]}" stroke-width="{0.8*s}"/>'
    o+=f'<ellipse cx="{tx-4*s}" cy="{ty-3*s}" rx="{4*s}" ry="{3*s}" fill="#FFFFFF" fill-opacity="0.25"/>'
    return o

def icon_radiation(cx,cy,s=1.0):
    """Clean radiation symbol — 3-blade trefoil."""
    o=""
    o+=f'<circle cx="{cx}" cy="{cy}" r="{5*s}" fill="{C["yellow"]}" stroke="{C["dkgray"]}" stroke-width="{LW*s*0.7}"/>'
    for ang_deg in [0,120,240]:
        a=math.radians(ang_deg-90)
        # blade
        bx1=cx+7*s*math.cos(a-0.35); by1=cy+7*s*math.sin(a-0.35)
        bx2=cx+7*s*math.cos(a+0.35); by2=cy+7*s*math.sin(a+0.35)
        bx3=cx+20*s*math.cos(a); by3=cy+20*s*math.sin(a)
        o+=f'<path d="M{bx1:.1f},{by1:.1f} Q{bx3:.1f},{by3:.1f} {bx2:.1f},{by2:.1f} Z" fill="{C["yellow"]}" stroke="{C["dkgray"]}" stroke-width="{LW*s*0.5}"/>'
    return o

def icon_dna_helix(cx,cy,w=60,h=20,s=1.0,c1=None,c2=None):
    c1=c1 or C["blue"]; c2=c2 or C["sky"]
    w*=s; h*=s; pts_a=[]; pts_b=[]
    for i in range(41):
        t=i/40; x=cx-w/2+t*w
        pts_a.append(f"{x:.1f},{cy+math.sin(t*math.pi*3)*h/2:.1f}")
        pts_b.append(f"{x:.1f},{cy-math.sin(t*math.pi*3)*h/2:.1f}")
    rungs=""
    for i in range(2,40,5):
        t=i/40; x=cx-w/2+t*w
        y1=cy+math.sin(t*math.pi*3)*h/2; y2=cy-math.sin(t*math.pi*3)*h/2
        rungs+=f'<line x1="{x:.1f}" y1="{y1:.1f}" x2="{x:.1f}" y2="{y2:.1f}" stroke="{C["ltgray"]}" stroke-width="{0.8*s}"/>'
    return (rungs
        +f'<polyline points="{" ".join(pts_a)}" fill="none" stroke="{c1}" stroke-width="{2*s}" stroke-linecap="round"/>'
        +f'<polyline points="{" ".join(pts_b)}" fill="none" stroke="{c2}" stroke-width="{2*s}" stroke-linecap="round"/>')

def icon_dna_break(cx,cy,s=1.0):
    """DNA breaking apart — clean style."""
    o=""
    o+=icon_dna_helix(cx-18*s,cy,w=32,h=16,s=s)
    # gap + bolt
    o+=f'<path d="M{cx-1*s},{cy-10*s} l{3*s},{6*s} l{-3*s},{0} l{4*s},{10*s} l{-2*s},{-5*s} l{3*s},{0} Z" fill="{C["orange"]}" fill-opacity="0.8"/>'
    # right fragment
    o+=f'<g opacity="0.35" transform="translate({5*s},{-3*s})">'
    o+=icon_dna_helix(cx+20*s,cy,w=22,h=12,s=s*0.7)
    o+='</g>'
    for dx,dy in [(32,-8),(36,6),(30,12)]:
        o+=f'<circle cx="{cx+dx*s}" cy="{cy+dy*s}" r="{2.5*s}" fill="{C["orange"]}" fill-opacity="0.4"/>'
    return o

def icon_mrna(cx,cy,w=55,s=1.0,color=None):
    col=color or C["green"]; w*=s; pts=[]
    for i in range(41):
        t=i/40; x=cx-w/2+t*w; y=cy+math.sin(t*math.pi*5)*5.5*s
        pts.append(f"{x:.1f},{y:.1f}")
    return f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="{2.2*s}" stroke-linecap="round"/>'

def icon_cell(cx,cy,r=16,s=1.0,color=None,label="",labelsize=10):
    """Generic cell with membrane and nucleus."""
    col=color or C["blue"]; r*=s; o=""
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" fill-opacity="0.15" stroke="{col}" stroke-width="{LW*s}"/>'
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r*0.4}" fill="{col}" fill-opacity="0.35"/>'
    if label:
        o+=txt(cx,cy+r+12*s,label,size=labelsize*s,color=col,weight=600,anchor="middle")
    return o

def icon_tcell(cx,cy,s=1.0,label="CD8+"):
    """T-cell with TCR."""
    col=C["blue"]; r=14*s; o=""
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" fill-opacity="0.18" stroke="{col}" stroke-width="{LW*s}"/>'
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r*0.38}" fill="{col}" fill-opacity="0.4"/>'
    # TCR (simplified Y shape on top)
    ty=cy-r
    o+=f'<line x1="{cx}" y1="{ty}" x2="{cx}" y2="{ty-8*s}" stroke="{col}" stroke-width="{LW*s}" stroke-linecap="round"/>'
    o+=f'<line x1="{cx}" y1="{ty-8*s}" x2="{cx-4*s}" y2="{ty-13*s}" stroke="{col}" stroke-width="{LW*s}" stroke-linecap="round"/>'
    o+=f'<line x1="{cx}" y1="{ty-8*s}" x2="{cx+4*s}" y2="{ty-13*s}" stroke="{col}" stroke-width="{LW*s}" stroke-linecap="round"/>'
    if label:
        o+=txt(cx,cy+r+12*s,label,size=9*s,color=col,weight=600,anchor="middle")
    return o

def icon_bcell(cx,cy,s=1.0,label="B cell"):
    """B-cell with surface antibodies (Y shapes)."""
    col=C["purple"]; r=14*s; o=""
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" fill-opacity="0.15" stroke="{col}" stroke-width="{LW*s}"/>'
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r*0.35}" fill="{col}" fill-opacity="0.35"/>'
    for ang_deg in [0,120,240]:
        a=math.radians(ang_deg-90)
        bx=cx+(r+1*s)*math.cos(a); by=cy+(r+1*s)*math.sin(a)
        ex=cx+(r+9*s)*math.cos(a); ey=cy+(r+9*s)*math.sin(a)
        o+=f'<line x1="{bx:.1f}" y1="{by:.1f}" x2="{ex:.1f}" y2="{ey:.1f}" stroke="{col}" stroke-width="{LW*s*0.8}" stroke-linecap="round"/>'
        for da in [-22,22]:
            a2=math.radians(ang_deg-90+da)
            o+=f'<line x1="{ex:.1f}" y1="{ey:.1f}" x2="{ex+5*s*math.cos(a2):.1f}" y2="{ey+5*s*math.sin(a2):.1f}" stroke="{col}" stroke-width="{LW*s*0.7}" stroke-linecap="round"/>'
    if label:
        o+=txt(cx,cy+r+12*s,label,size=9*s,color=col,weight=600,anchor="middle")
    return o

def icon_treg(cx,cy,s=1.0):
    """Treg — cell with dashed suppressive halo."""
    col=C["green"]; r=13*s; o=""
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" fill-opacity="0.18" stroke="{col}" stroke-width="{LW*s}"/>'
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r*0.35}" fill="{col}" fill-opacity="0.4"/>'
    o+=f'<circle cx="{cx}" cy="{cy}" r="{r+6*s}" fill="none" stroke="{col}" stroke-width="{0.9*s}" stroke-dasharray="3,2.5"/>'
    o+=txt(cx,cy+r+15*s,"Treg",size=9*s,color=col,weight=600,anchor="middle")
    return o

def icon_neoantigen(cx,cy,s=1.0,faded=False):
    """Star-shaped neoantigen."""
    op=0.25 if faded else 0.85; pts=[]
    for i in range(10):
        a=math.radians(i*36-90); r_=(12 if i%2==0 else 6)*s
        pts.append(f"{cx+r_*math.cos(a):.1f},{cy+r_*math.sin(a):.1f}")
    return f'<polygon points="{" ".join(pts)}" fill="{C["orange"]}" fill-opacity="{op}" stroke="{C["orange"]}" stroke-width="{0.8*s}" stroke-opacity="{op}"/>'

def icon_mhc(cx,cy,s=1.0):
    """MHC-I groove with peptide."""
    o=""
    # groove
    o+=f'<rect x="{cx-14*s}" y="{cy-8*s}" width="{28*s}" height="{12*s}" rx="{4*s}" fill="{C["greenBg"]}" stroke="{C["green"]}" stroke-width="{LW*s}"/>'
    # peptide
    o+=f'<rect x="{cx-8*s}" y="{cy-5*s}" width="{16*s}" height="{6*s}" rx="3" fill="{C["orange"]}" fill-opacity="0.7"/>'
    # stalk
    o+=f'<rect x="{cx-2.5*s}" y="{cy+4*s}" width="{5*s}" height="{14*s}" rx="2" fill="{C["green"]}" fill-opacity="0.5"/>'
    return o

def icon_person(cx,cy,s=1.0,color=None):
    col=color or C["blue"]
    return (f'<circle cx="{cx}" cy="{cy-8*s}" r="{6*s}" fill="{col}"/>'
            f'<path d="M{cx-8*s},{cy+8*s} Q{cx-8*s},{cy-1*s} {cx},{cy-1*s} Q{cx+8*s},{cy-1*s} {cx+8*s},{cy+8*s} Z" fill="{col}"/>')

def icon_pills(cx,cy,s=1.0):
    o=""
    o+=f'<rect x="{cx-10*s}" y="{cy-4*s}" width="{20*s}" height="{8*s}" rx="{4*s}" fill="{C["purple"]}" fill-opacity="0.6" stroke="{C["purple"]}" stroke-width="{0.8*s}"/>'
    o+=f'<line x1="{cx}" y1="{cy-4*s}" x2="{cx}" y2="{cy+4*s}" stroke="{C["purple"]}" stroke-width="{0.8*s}"/>'
    return o

def icon_tumor_shrink(cx,cy,s=1.0):
    """Before/after tumor size."""
    o=""
    o+=f'<ellipse cx="{cx-16*s}" cy="{cy}" rx="{14*s}" ry="{11*s}" fill="{C["orange"]}" fill-opacity="0.6"/>'
    o+=f'<text x="{cx}" y="{cy+4*s}" font-family=\'{FONT}\' font-size="{14*s}" fill="{C["gray"]}" text-anchor="middle">&#x2192;</text>'
    o+=f'<ellipse cx="{cx+20*s}" cy="{cy}" rx="{6*s}" ry="{5*s}" fill="{C["orange"]}" fill-opacity="0.2" stroke="{C["orange"]}" stroke-width="{0.6*s}" stroke-dasharray="2,2"/>'
    return o

def icon_biopsy(cx,cy,s=1.0,angle=-20):
    return (f'<g transform="translate({cx},{cy}) rotate({angle}) scale({s})">'
            f'<rect x="-22" y="-3" width="30" height="6" rx="2" fill="#7B8D9E"/>'
            f'<rect x="8" y="-1.5" width="20" height="3" rx="1.5" fill="#9FADB8"/>'
            f'<polygon points="28,-2.5 36,0 28,2.5" fill="#6B7C8F"/>'
            f'</g>')

# ========== BUILD SVG ==========
P = []
P.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')

defs = ['<defs>']
for col in [C["dkgray"],C["blue"],C["orange"],C["green"],C["purple"],C["gray"],C["sky"]]:
    defs.append(arrow_marker(col))
defs.append('</defs>')
P.extend(defs)

# background
P.append(rect(0,0,W,H,C["white"],rx=0))

# ---- TITLE (centered, compact) ----
P.append(txt(W/2, 42,
    "Radiation-Phase Molecular Response Predicts TNT Outcome in MSS LARC",
    size=28, color=C["dkgray"], weight=700, anchor="middle", ls=0.3))
P.append(txt(W/2, 68,
    "Integrated WES + RNA-seq  |  35 locally advanced rectal cancers  |  pre / post chemoradiotherapy biopsies",
    size=13.5, color=C["gray"], weight=400, anchor="middle"))

# thin rule
P.append(f'<line x1="120" y1="86" x2="{W-120}" y2="86" stroke="{C["ltgray"]}" stroke-width="1"/>')

# ========== SECTION 1: Study Design (left) ==========
S1X, S1Y, S1W, S1H = 50, 100, 430, 580
P.append(rect(S1X, S1Y, S1W, S1H, C["bgcard"], stroke=C["ltgray"], sw=1, rx=16))

# section label
P.append(rect(S1X+16, S1Y+14, 130, 26, C["blue"], rx=6))
P.append(txt(S1X+81, S1Y+32, "STUDY DESIGN", size=11.5, color=C["white"], weight=700, anchor="middle", ls=1))

# colon + tumor + radiation
P.append(icon_colon_section(S1X+100, S1Y+110, s=1.8))
P.append(icon_radiation(S1X+195, S1Y+95, s=1.0))
P.append(txt(S1X+100, S1Y+160, "LARC + TNT", size=11, color=C["orange"], weight=600, anchor="middle"))
P.append(txt(S1X+100, S1Y+175, "SC-RT 25 Gy/5 Fx + FOLFOX/CAPOX", size=9.5, color=C["gray"], anchor="middle"))

# patient grid (7×5 = 35)
pgx0, pgy0 = S1X+260, S1Y+68
for i in range(35):
    r, c = divmod(i, 7)
    cx = pgx0 + c * 20; cy = pgy0 + r * 24
    col = C["blue"] if i < 18 else C["orange"]
    P.append(icon_person(cx, cy, s=0.55, color=col))
# legend
P.append(f'<circle cx="{pgx0}" cy="{S1Y+195}" r="4" fill="{C["blue"]}"/>')
P.append(txt(pgx0+10, S1Y+199, "18 good", size=10, color=C["blue"], weight=600))
P.append(f'<circle cx="{pgx0+80}" cy="{S1Y+195}" r="4" fill="{C["orange"]}"/>')
P.append(txt(pgx0+90, S1Y+199, "17 bad", size=10, color=C["orange"], weight=600))

# ---- timeline ----
tly = S1Y + 240
P.append(f'<line x1="{S1X+30}" y1="{tly}" x2="{S1X+S1W-30}" y2="{tly}" stroke="{C["dkgray"]}" stroke-width="2"/>')
# radiation band
P.append(rect(S1X+110, tly-12, 190, 24, C["yellow"], stroke=C["dkgray"], sw=0.8, rx=6, op=0.35))
P.append(txt(S1X+205, tly+5, "radiation phase", size=10, color=C["dkgray"], weight=600, anchor="middle"))
# biopsy points
for bx, lab in [(S1X+90, "pre"), (S1X+320, "post")]:
    P.append(icon_biopsy(bx, tly-24, s=0.8))
    P.append(f'<circle cx="{bx}" cy="{tly}" r="5" fill="{C["blue"]}" stroke="{C["white"]}" stroke-width="1.5"/>')
    P.append(txt(bx, tly+18, lab, size=10.5, color=C["dkgray"], weight=600, anchor="middle"))
P.append(f'<circle cx="{S1X+S1W-40}" cy="{tly}" r="5" fill="{C["dkgray"]}" stroke="{C["white"]}" stroke-width="1.5"/>')
P.append(txt(S1X+S1W-40, tly+18, "final response", size=9.5, color=C["gray"], anchor="middle"))

# ---- omics cards ----
oy = S1Y + 290
P.append(rect(S1X+20, oy, 185, 58, C["blueBg"], stroke=C["blue"], sw=LW*0.8, rx=10))
P.append(icon_dna_helix(S1X+60, oy+29, w=50, h=16, s=0.9))
P.append(txt(S1X+100, oy+22, "WES", size=14, color=C["blue"], weight=700))
P.append(txt(S1X+100, oy+40, "Mutect2 T-N paired", size=9.5, color=C["gray"]))

P.append(rect(S1X+220, oy, 185, 58, C["greenBg"], stroke=C["green"], sw=LW*0.8, rx=10))
P.append(icon_mrna(S1X+260, oy+29, w=45, s=0.9))
P.append(txt(S1X+300, oy+22, "RNA-seq", size=14, color=C["green"], weight=700))
P.append(txt(S1X+300, oy+40, "DESeq2 / fgsea / ssGSEA", size=9.5, color=C["gray"]))

# key genomic findings
gy = oy + 78
findings = [
    ("All MSS (0 MSI-H)", C["dkgray"]),
    ("TMB < 2 / Mb (low)", C["dkgray"]),
    ("APC > TP53 > KRAS drivers", C["dkgray"]),
    ("GSEA: E2F / G2M / DSB repair", C["blue"]),
    ("good up; EMT bad up", C["orange"]),
]
for i, (f_txt, f_col) in enumerate(findings):
    yy = gy + i * 18
    P.append(f'<circle cx="{S1X+36}" cy="{yy}" r="2.5" fill="{f_col}"/>')
    P.append(txt(S1X+46, yy+4, f_txt, size=10.5, color=C["dkgray"]))

# nested CV AUC
P.append(rect(S1X+20, gy+100, S1W-40, 48, C["bgcard"], stroke=C["blue"], sw=LW*0.8, rx=10))
P.append(txt(S1X+S1W/2, gy+120, "Nested LOOCV classifier", size=12, color=C["dkgray"], weight=600, anchor="middle"))
P.append(txt(S1X+S1W/2, gy+138, "AUC = 0.65 (95% CI 0.45-0.83) | modest discovery-stage", size=10.5, color=C["gray"], anchor="middle"))

# ========== SECTION 2: Mechanism / Cascade (center) ==========
S2X, S2Y, S2W, S2H = 510, 100, 580, 580
P.append(rect(S2X, S2Y, S2W, S2H, C["white"], stroke=C["ltgray"], sw=1, rx=16))

P.append(rect(S2X+16, S2Y+14, 230, 26, C["green"], rx=6))
P.append(txt(S2X+131, S2Y+32, "GOOD-RESPONDER CASCADE", size=11.5, color=C["white"], weight=700, anchor="middle", ls=1))

P.append(txt(S2X+S2W/2, S2Y+58, "Radiation induces sequential tumor clearance", size=12, color=C["gray"], anchor="middle"))
P.append(txt(S2X+S2W/2, S2Y+74, "and immune reprogramming", size=12, color=C["gray"], anchor="middle"))

# tumor shrink illustration at top right
P.append(icon_tumor_shrink(S2X+S2W-80, S2Y+48, s=1.1))

# 5-step cascade
steps = [
    {"n":"1","title":"Mutation clearance",
     "sub":"SBS5 delta -76 | missense delta -67",
     "col":C["blue"],
     "icon": lambda cx,cy: icon_dna_break(cx,cy,s=1.15)},
    {"n":"2","title":"Neoantigen depletion",
     "sub":"MHC-I binders delta -312 | P = 0.082",
     "col":C["orange"],
     "icon": lambda cx,cy: (icon_neoantigen(cx-14,cy,s=1.3)
                            +icon_neoantigen(cx+10,cy-6,s=0.8,faded=True)
                            +icon_neoantigen(cx+20,cy+8,s=0.5,faded=True))},
    {"n":"3","title":"HLA-LOH clone removal",
     "sub":"subj 3 & 4: pre LOH resolved post-CRT",
     "col":C["green"],
     "icon": lambda cx,cy: icon_mhc(cx,cy,s=1.4)},
    {"n":"4","title":"Immune reprogramming",
     "sub":"Treg P = 0.026 | MHC-II up | CD8 exhaustion up",
     "col":C["green"],
     "icon": lambda cx,cy: (icon_tcell(cx-18,cy-4,s=0.7,label="")
                            +icon_treg(cx+16,cy+6,s=0.6))},
    {"n":"5","title":"B-cell infiltration",
     "sub":"IGH / IGK / IGL clones up ~+1,400",
     "col":C["purple"],
     "icon": lambda cx,cy: (icon_bcell(cx-14,cy,s=0.8,label="")
                            +icon_bcell(cx+16,cy-6,s=0.55,label=""))},
]

sy0 = S2Y + 94
rh = 92
for i, st in enumerate(steps):
    y = sy0 + i * rh
    col = st["col"]
    # number badge
    P.append(f'<circle cx="{S2X+40}" cy="{y+30}" r="18" fill="{col}"/>')
    P.append(txt(S2X+40, y+36, st["n"], size=18, color=C["white"], weight=700, anchor="middle"))
    # icon
    P.append(st["icon"](S2X+100, y+30))
    # text card
    tcx = S2X+150
    P.append(rect(tcx, y+4, S2W-180, 56, C["bgcard"], stroke=col, sw=1.2, rx=10))
    P.append(txt(tcx+14, y+28, st["title"], size=13.5, color=C["dkgray"], weight=700))
    P.append(txt(tcx+14, y+46, st["sub"], size=10.5, color=C["gray"]))
    # connector arrow
    if i < len(steps)-1:
        P.append(arrow_line(S2X+40, y+50, S2X+40, y+rh+10, color=col, sw=2.5))

# ========== SECTION 3: Validation + Clinical (right) ==========
S3X, S3Y, S3W, S3H = 1120, 100, 630, 580
P.append(rect(S3X, S3Y, S3W, S3H, C["bgcard"], stroke=C["ltgray"], sw=1, rx=16))

# -- top: external validation --
P.append(rect(S3X+16, S3Y+14, 240, 26, C["sky"], rx=6))
P.append(txt(S3X+136, S3Y+32, "EXTERNAL CD8 VALIDATION", size=11.5, color=C["white"], weight=700, anchor="middle", ls=1))

P.append(txt(S3X+20, S3Y+62, "5 LC-CRT + Akiyoshi (N = 816) + SC-RT GSE254249 (N = 8)", size=10.5, color=C["gray"]))

# CD8 T-cell hero icon
P.append(icon_tcell(S3X+55, S3Y+120, s=1.6, label=""))
P.append(txt(S3X+55, S3Y+150, "CD8+", size=11, color=C["blue"], weight=700, anchor="middle"))
P.append(txt(S3X+55, S3Y+164, "cytotoxic", size=9.5, color=C["blue"], weight=600, anchor="middle"))

# mini forest plot
fx0 = S3X+160; zero_x = fx0 + 50; scale_f = 120
cohorts = [
    ("GSE45404",+0.55,0.30),("GSE46862",+0.70,0.35),("GSE87211",+0.40,0.28),
    ("GSE94104",+0.62,0.33),("GSE133057",+0.35,0.30),("GSE35452",+0.58,0.32),
    ("GSE119409",+0.25,0.34),("GSE150082",+0.32,0.30),("GSE56699",-0.10,0.30),
    ("Akiyoshi '23",+0.75,0.22),
]
fy0 = S3Y + 75; rh_f = 17
# zero line
P.append(f'<line x1="{zero_x}" y1="{fy0}" x2="{zero_x}" y2="{fy0+len(cohorts)*rh_f+8}" stroke="{C["ltgray"]}" stroke-width="1" stroke-dasharray="3,3"/>')
for i,(name,eff,se) in enumerate(cohorts):
    yy = fy0 + i*rh_f + 8
    cx_f = zero_x + eff*scale_f
    lo = zero_x + (eff-1.96*se)*scale_f
    hi = zero_x + (eff+1.96*se)*scale_f
    c_f = C["blue"] if eff>0 else C["orange"]
    is_ak = "Akiyoshi" in name
    P.append(txt(fx0-4, yy+4, name, size=9.5, color=C["dkgray"], anchor="end", weight=700 if is_ak else 400))
    P.append(f'<line x1="{max(lo, zero_x-60)}" y1="{yy}" x2="{min(hi, zero_x+150)}" y2="{yy}" stroke="{C["purple"] if is_ak else c_f}" stroke-width="1.4" stroke-linecap="round"/>')
    if is_ak:
        P.append(f'<polygon points="{cx_f-4},{yy} {cx_f},{yy-4} {cx_f+4},{yy} {cx_f},{yy+4}" fill="{C["purple"]}"/>')
    else:
        P.append(f'<rect x="{cx_f-3.5}" y="{yy-3.5}" width="7" height="7" fill="{c_f}"/>')

# meta diamond
mety = fy0 + len(cohorts)*rh_f + 20
meta_cx = zero_x + 0.52*scale_f
P.append(f'<polygon points="{meta_cx-30},{mety} {meta_cx},{mety-7} {meta_cx+30},{mety} {meta_cx},{mety+7}" fill="{C["dkgray"]}"/>')
P.append(txt(fx0-4, mety+4, "Meta (9)", size=10, color=C["dkgray"], anchor="end", weight=700))
P.append(txt(meta_cx+36, mety+4, "Z = +2.74  P = 0.006", size=10.5, color=C["dkgray"], weight=700))

# axis ticks
axy = mety + 18
for v, lab in [(-0.5,"-0.5"),(0,"0"),(0.5,"+0.5"),(1.0,"+1.0")]:
    P.append(txt(zero_x+v*scale_f, axy, lab, size=8.5, color=C["gray"], anchor="middle"))
P.append(txt(zero_x+0.25*scale_f, axy+14, "Effect size (good vs bad)", size=9, color=C["gray"], anchor="middle"))

# -- bottom: clinical decision --
cdy = S3Y + 320
P.append(f'<line x1="{S3X+20}" y1="{cdy}" x2="{S3X+S3W-20}" y2="{cdy}" stroke="{C["ltgray"]}" stroke-width="1"/>')

P.append(rect(S3X+16, cdy+10, 260, 26, C["dkgray"], rx=6))
P.append(txt(S3X+146, cdy+28, "MID-TREATMENT DECISION WINDOW", size=11, color=C["white"], weight=700, anchor="middle", ls=1))

# decision box
dbx = S3X + S3W/2; dby = cdy + 60
P.append(rect(dbx-160, dby-14, 320, 32, C["white"], stroke=C["dkgray"], sw=LW, rx=8))
P.append(txt(dbx, dby+5, "Post-CRT biopsy molecular profile", size=12, color=C["dkgray"], weight=600, anchor="middle"))

# branches
# good
P.append(curved_arrow(dbx-60, dby+20, dbx-130, dby+60, color=C["blue"], sw=2.2, bend=-25))
P.append(rect(S3X+30, dby+55, 240, 90, C["blueBg"], stroke=C["blue"], sw=LW, rx=12))
P.append(icon_person(S3X+60, dby+85, s=0.7, color=C["blue"]))
P.append(txt(S3X+85, dby+82, "Responder profile", size=12.5, color=C["blue"], weight=700))
P.append(txt(S3X+85, dby+100, "De-escalate / watch-and-wait", size=10.5, color=C["dkgray"]))
P.append(txt(S3X+85, dby+116, "Organ preservation", size=10, color=C["gray"]))

# bad
P.append(curved_arrow(dbx+60, dby+20, dbx+130, dby+60, color=C["orange"], sw=2.2, bend=25))
P.append(rect(S3X+S3W-275, dby+55, 240, 90, C["orangeBg"], stroke=C["orange"], sw=LW, rx=12))
P.append(icon_pills(S3X+S3W-248, dby+88, s=1.1))
P.append(txt(S3X+S3W-218, dby+82, "Non-responder profile", size=12.5, color=C["orange"], weight=700))
P.append(txt(S3X+S3W-218, dby+100, "Intensify consolidation", size=10.5, color=C["dkgray"]))
P.append(txt(S3X+S3W-218, dby+116, "FOLFOX / trial enrollment", size=10, color=C["gray"]))

# immune microenvironment illustration between the decision branches
P.append(icon_tcell(dbx-6, dby+100, s=0.55, label=""))
P.append(icon_bcell(dbx+20, dby+110, s=0.45, label=""))
P.append(icon_treg(dbx+6, dby+130, s=0.4))

# ========== SECTION-TO-SECTION FLOW ARROWS ==========
# large connecting arrows between sections (clean)
for (ax1,ay1,ax2,ay2) in [
    (S1X+S1W+4, S1Y+S1H/2, S2X-4, S2Y+S2H/2),
    (S2X+S2W+4, S2Y+S2H/3, S3X-4, S3Y+S3H/3),
    (S2X+S2W+4, S2Y+S2H*0.7, S3X-4, cdy+40),
]:
    mx=(ax1+ax2)/2
    P.append(f'<path d="M{ax1},{ay1} C{mx},{ay1} {mx},{ay2} {ax2},{ay2}" fill="none" '
             f'stroke="{C["ltgray"]}" stroke-width="3" stroke-linecap="round" marker-end="url(#arr{C["gray"][1:]})"/>')

# ========== BOTTOM TAKE-HOME ==========
bty = 700
P.append(rect(50, bty, W-100, 60, C["dkgray"], rx=12))
P.append(txt(W/2, bty+25,
    "Regimen-agnostic pre-CRT baseline predictor (Thread 1 DSB/cell-cycle + Thread 2 CD8-cytotoxic) reproduces across",
    size=14, color=C["white"], weight=600, anchor="middle"))
P.append(txt(W/2, bty+47,
    "SC-RT + LC-CRT regimens (N > 850) -- radiation-phase pharmacodynamics (target engagement + IGH directional coherence) add an orthogonal dynamic layer.",
    size=13.5, color=C["sky"], weight=400, anchor="middle"))

P.append("</svg>")

svg = "\n".join(P)
svg_path = OUT / "GA_v3_master.svg"
svg_path.write_text(svg, encoding="utf-8")
print("wrote:", svg_path)

if HAVE_CAIRO:
    cairosvg.svg2png(bytestring=svg.encode("utf-8"), write_to=str(OUT/"GA_v3.png"), output_width=W*2, output_height=H*2)
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(OUT/"GA_v3.pdf"))
    print("wrote:", OUT/"GA_v3.png")
    print("wrote:", OUT/"GA_v3.pdf")
else:
    # Fall back to rsvg-convert or inkscape if present
    import subprocess, shutil
    for tool in ["rsvg-convert", "inkscape"]:
        if shutil.which(tool):
            if tool == "rsvg-convert":
                subprocess.run([tool, "-f", "pdf", "-o", str(OUT/"GA_v3.pdf"), str(svg_path)], check=False)
                subprocess.run([tool, "-f", "png", "-z", "2", "-o", str(OUT/"GA_v3.png"), str(svg_path)], check=False)
            else:
                subprocess.run([tool, str(svg_path), "--export-type=pdf", f"--export-filename={OUT}/GA_v3.pdf"], check=False)
                subprocess.run([tool, str(svg_path), "--export-type=png", f"--export-filename={OUT}/GA_v3.png", "--export-dpi=200"], check=False)
            print(f"rasterised via {tool}")
            break
    else:
        print("(No cairosvg/rsvg/inkscape; SVG-only output — manual rasterisation needed.)")
