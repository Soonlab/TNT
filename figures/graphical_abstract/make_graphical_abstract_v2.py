#!/usr/bin/env python3
"""Graphical abstract v2 — BioRender-style with rich biological icons."""
import math
from pathlib import Path
import cairosvg

OUT = Path("/data/data/TNT/analysis/figures/graphical_abstract")
OUT.mkdir(parents=True, exist_ok=True)

C = {
    "good":   "#2E86AB", "goodL":  "#A8D8EA", "goodBg": "#E6F1F6",
    "bad":    "#E63946", "badL":   "#F4A0A8", "badBg":  "#FCE7E9",
    "amber":  "#F4A261", "amberL": "#FDEFDD",
    "sage":   "#52B788", "sageL":  "#C8E6D5", "sageBg": "#E3F2E9",
    "purple": "#7B68AE", "purpleL":"#C9BFDF",
    "slate":  "#2A2E3A", "ink":    "#1F2430",
    "muted":  "#6B7280", "line":   "#C7CEDB",
    "bg":     "#FFFFFF", "soft":   "#F7F8FA", "softer": "#EEF2F6",
    "skin":   "#F5D6BA", "skinD":  "#D4A574",
}

W, H = 1800, 810
FONT = 'Inter, "Helvetica Neue", Arial, sans-serif'

def _esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

def rect(x,y,w,h,fill,stroke=None,sw=1.2,rx=14,op=1.0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" ry="{rx}" fill="{fill}" fill-opacity="{op}"{s}/>'

def txt(x,y,s,size=14,color=None,weight=400,anchor="start"):
    col=color or C["ink"]; s=_esc(s)
    return f'<text x="{x}" y="{y}" font-family=\'{FONT}\' font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="{anchor}">{s}</text>'

def arrow_marker(color):
    cid=f"arrow-{color[1:]}"
    return (f'<marker id="{cid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')

def arrow_line(x1,y1,x2,y2,color=None,sw=2.4,dash=None):
    col=color or C["slate"]
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" '
            f'stroke-width="{sw}" stroke-linecap="round" marker-end="url(#arrow-{col[1:]})"{d}/>')

def curved_arrow(x1,y1,x2,y2,color=None,sw=2.8,bend=40):
    col=color or C["slate"]
    mx=(x1+x2)/2; my=(y1+y2)/2
    dx=x2-x1; dy=y2-y1
    L=math.hypot(dx,dy)
    if L==0: return ""
    nx=-dy/L*bend; ny=dx/L*bend
    cpx=mx+nx; cpy=my+ny
    return (f'<path d="M{x1},{y1} Q{cpx},{cpy} {x2},{y2}" fill="none" stroke="{col}" '
            f'stroke-width="{sw}" stroke-linecap="round" marker-end="url(#arrow-{col[1:]})"/>')

# ---- biological icon library ----

def icon_colon_tumor(cx,cy,scale=1.0):
    """Stylized colon cross-section with tumor mass."""
    s=scale; out=""
    # colon wall (pink tube cross-section)
    out+=f'<ellipse cx="{cx}" cy="{cy}" rx="{38*s}" ry="{32*s}" fill="#F8C8C8" stroke="#D4726A" stroke-width="{2*s}"/>'
    out+=f'<ellipse cx="{cx}" cy="{cy}" rx="{24*s}" ry="{20*s}" fill="#FCEAE8" stroke="#D4726A" stroke-width="{1.2*s}"/>'
    # lumen
    out+=f'<ellipse cx="{cx}" cy="{cy}" rx="{14*s}" ry="{10*s}" fill="#FDF4F2"/>'
    # tumor nodule
    tx=cx+12*s; ty=cy-8*s
    out+=f'<ellipse cx="{tx}" cy="{ty}" rx="{16*s}" ry="{13*s}" fill="{C["bad"]}" fill-opacity="0.85"/>'
    out+=f'<ellipse cx="{tx-3*s}" cy="{ty-3*s}" rx="{5*s}" ry="{4*s}" fill="{C["badL"]}" fill-opacity="0.5"/>'
    # small satellite
    out+=f'<circle cx="{cx+30*s}" cy="{cy+5*s}" r="{6*s}" fill="{C["bad"]}" fill-opacity="0.6"/>'
    return out

def icon_radiation(cx,cy,scale=1.0):
    """Radiation trefoil / beam icon."""
    s=scale; out=""
    # central circle
    out+=f'<circle cx="{cx}" cy="{cy}" r="{8*s}" fill="{C["amber"]}" stroke="#D4872A" stroke-width="{1.5*s}"/>'
    # 6 beams radiating out
    for ang_deg in range(0,360,60):
        a=math.radians(ang_deg)
        x1=cx+10*s*math.cos(a); y1=cy+10*s*math.sin(a)
        x2=cx+26*s*math.cos(a); y2=cy+26*s*math.sin(a)
        out+=f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{C["amber"]}" stroke-width="{2.5*s}" stroke-linecap="round"/>'
    # small zigzag bolts
    for ang_deg in [30,150,270]:
        a=math.radians(ang_deg)
        bx=cx+20*s*math.cos(a); by=cy+20*s*math.sin(a)
        out+=f'<polygon points="{bx},{by-4*s} {bx+3*s},{by} {bx},{by+4*s} {bx-3*s},{by}" fill="{C["amber"]}" fill-opacity="0.6"/>'
    return out

def icon_biopsy(cx,cy,angle=-25,scale=1.0):
    """Biopsy needle."""
    s=scale
    return (f'<g transform="translate({cx},{cy}) rotate({angle}) scale({s})">'
            f'<rect x="-30" y="-4" width="40" height="8" rx="2" fill="#8B9DAF"/>'
            f'<rect x="10" y="-2" width="26" height="4" rx="2" fill="#B0BFCF"/>'
            f'<polygon points="36,-3 46,0 36,3" fill="#6B7C8F"/>'
            f'</g>')

def icon_dna_helix(cx,cy,w=70,h=28,color1=None,color2=None,scale=1.0):
    """Double helix."""
    col1=color1 or C["good"]; col2=color2 or C["amber"]
    s=scale; w*=s; h*=s
    pts_a,pts_b=[],[]
    for i in range(41):
        t=i/40; x=cx-w/2+t*w
        y1=cy+math.sin(t*math.pi*3)*h/2
        y2=cy-math.sin(t*math.pi*3)*h/2
        pts_a.append(f"{x:.1f},{y1:.1f}")
        pts_b.append(f"{x:.1f},{y2:.1f}")
    rungs=""
    for i in range(2,40,4):
        t=i/40; x=cx-w/2+t*w
        y1=cy+math.sin(t*math.pi*3)*h/2
        y2=cy-math.sin(t*math.pi*3)*h/2
        rungs+=f'<line x1="{x:.1f}" y1="{y1:.1f}" x2="{x:.1f}" y2="{y2:.1f}" stroke="{C["line"]}" stroke-width="{1*s}"/>'
    a=f'<polyline points="{" ".join(pts_a)}" fill="none" stroke="{col1}" stroke-width="{2.5*s}" stroke-linecap="round"/>'
    b=f'<polyline points="{" ".join(pts_b)}" fill="none" stroke="{col2}" stroke-width="{2.5*s}" stroke-linecap="round"/>'
    return rungs+a+b

def icon_dna_breaking(cx,cy,scale=1.0):
    """DNA helix with break / fragments flying off — mutation clearance."""
    s=scale; out=""
    # left half helix
    out+=icon_dna_helix(cx-22*s,cy,w=40,h=24,scale=s)
    # break gap — lightning bolt
    out+=f'<polygon points="{cx-2*s},{cy-16*s} {cx+4*s},{cy-4*s} {cx},{cy-4*s} {cx+6*s},{cy+14*s} {cx},{cy+2*s} {cx+4*s},{cy+2*s}" fill="{C["amber"]}" stroke="#D4872A" stroke-width="0.8"/>'
    # right fragment fading
    out+=f'<g opacity="0.4" transform="translate({6*s},{-4*s})">'
    out+=icon_dna_helix(cx+24*s,cy,w=30,h=18,scale=s*0.8)
    out+='</g>'
    # scattered nucleotide dots
    for dx,dy in [(35,12),(40,-8),(38,18),(44,4)]:
        out+=f'<circle cx="{cx+dx*s}" cy="{cy+dy*s}" r="{3*s}" fill="{C["bad"]}" fill-opacity="0.5"/>'
    return out

def icon_tcell(cx,cy,scale=1.0,color=None,label="T"):
    """T-cell with receptor protrusions."""
    col=color or C["good"]; s=scale; out=""
    r=18*s
    out+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" stroke="#FFFFFF" stroke-width="{2*s}"/>'
    # nucleus
    out+=f'<circle cx="{cx}" cy="{cy}" r="{8*s}" fill="{C["slate"]}" fill-opacity="0.25"/>'
    # TCR protrusions
    for ang_deg in [0,90,180,270]:
        a=math.radians(ang_deg)
        x1=cx+r*math.cos(a)*0.95; y1=cy+r*math.sin(a)*0.95
        x2=cx+(r+8*s)*math.cos(a); y2=cy+(r+8*s)*math.sin(a)
        out+=f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{col}" stroke-width="{2.5*s}" stroke-linecap="round"/>'
        out+=f'<circle cx="{x2:.1f}" cy="{y2:.1f}" r="{3*s}" fill="{col}"/>'
    out+=txt(cx,cy+5*s,label,size=int(13*s),color="#FFFFFF",weight=700,anchor="middle")
    return out

def icon_bcell(cx,cy,scale=1.0):
    """B-cell — rounder, with Y-shaped antibodies."""
    s=scale; out=""
    r=18*s; col=C["purple"]
    out+=f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" stroke="#FFFFFF" stroke-width="{2*s}"/>'
    out+=f'<circle cx="{cx}" cy="{cy}" r="{7*s}" fill="{C["purpleL"]}" fill-opacity="0.5"/>'
    # Y-shaped antibodies
    for ang_deg in [30,150,270]:
        a=math.radians(ang_deg)
        bx=cx+(r+2*s)*math.cos(a); by=cy+(r+2*s)*math.sin(a)
        ex=cx+(r+12*s)*math.cos(a); ey=cy+(r+12*s)*math.sin(a)
        out+=f'<line x1="{bx:.1f}" y1="{by:.1f}" x2="{ex:.1f}" y2="{ey:.1f}" stroke="{col}" stroke-width="{2*s}" stroke-linecap="round"/>'
        # Y arms
        for da in [-25,25]:
            a2=math.radians(ang_deg+da)
            ax=ex+8*s*math.cos(a2); ay=ey+8*s*math.sin(a2)
            out+=f'<line x1="{ex:.1f}" y1="{ey:.1f}" x2="{ax:.1f}" y2="{ay:.1f}" stroke="{col}" stroke-width="{1.8*s}" stroke-linecap="round"/>'
    out+=txt(cx,cy+5*s,"B",size=int(13*s),color="#FFFFFF",weight=700,anchor="middle")
    return out

def icon_treg(cx,cy,scale=1.0):
    """Treg — T-cell variant with halo ring."""
    s=scale; out=""
    out+=icon_tcell(cx,cy,scale=s,color=C["sage"],label="Treg")
    out+=f'<circle cx="{cx}" cy="{cy}" r="{26*s}" fill="none" stroke="{C["sage"]}" stroke-width="{1.5*s}" stroke-dasharray="4,3"/>'
    return out

def icon_hla_mhc(cx,cy,scale=1.0):
    """MHC-I molecule on cell surface — two alpha helices + peptide groove."""
    s=scale; out=""
    # cell membrane bar
    out+=f'<rect x="{cx-22*s}" y="{cy+12*s}" width="{44*s}" height="{6*s}" rx="3" fill="{C["line"]}"/>'
    # transmembrane stalk
    out+=f'<rect x="{cx-3*s}" y="{cy-4*s}" width="{6*s}" height="{18*s}" rx="2" fill="{C["sage"]}"/>'
    # alpha helices (groove)
    out+=f'<rect x="{cx-18*s}" y="{cy-18*s}" width="{36*s}" height="{16*s}" rx="{6*s}" fill="{C["sageBg"]}" stroke="{C["sage"]}" stroke-width="{1.5*s}"/>'
    # peptide in groove
    out+=f'<rect x="{cx-10*s}" y="{cy-14*s}" width="{20*s}" height="{8*s}" rx="4" fill="{C["amber"]}" fill-opacity="0.8"/>'
    out+=txt(cx,cy-8*s,"peptide",size=int(7.5*s),color=C["slate"],anchor="middle",weight=600)
    # label
    out+=txt(cx,cy+30*s,"MHC-I",size=int(10*s),color=C["muted"],anchor="middle",weight=600)
    return out

def icon_neoantigen(cx,cy,scale=1.0,faded=False):
    """Neoantigen — star-shaped mutant peptide."""
    s=scale; out=""
    op="0.3" if faded else "0.9"
    # star shape
    pts=[]
    for i in range(10):
        a=math.radians(i*36-90)
        r_=(14 if i%2==0 else 7)*s
        pts.append(f"{cx+r_*math.cos(a):.1f},{cy+r_*math.sin(a):.1f}")
    out+=f'<polygon points="{" ".join(pts)}" fill="{C["bad"]}" fill-opacity="{op}" stroke="{C["bad"]}" stroke-width="{1.2*s}" stroke-opacity="{op}"/>'
    if not faded:
        out+=txt(cx,cy+4*s,"neo",size=int(8*s),color="#FFFFFF",weight=700,anchor="middle")
    return out

def icon_person_simple(cx,cy,scale=1.0,color=None):
    """Simple person/patient icon."""
    col=color or C["good"]; s=scale
    return (f'<circle cx="{cx}" cy="{cy-10*s}" r="{7*s}" fill="{col}"/>'
            f'<path d="M{cx-10*s},{cy+10*s} Q{cx-10*s},{cy-2*s} {cx},{cy-2*s} Q{cx+10*s},{cy-2*s} {cx+10*s},{cy+10*s} Z" fill="{col}"/>')

def icon_pills(cx,cy,scale=1.0):
    """Chemotherapy pills."""
    s=scale; out=""
    # capsule 1
    out+=f'<rect x="{cx-14*s}" y="{cy-5*s}" width="{28*s}" height="{10*s}" rx="{5*s}" fill="#E8A0E0" stroke="#C070B8" stroke-width="{1*s}"/>'
    out+=f'<line x1="{cx}" y1="{cy-5*s}" x2="{cx}" y2="{cy+5*s}" stroke="#C070B8" stroke-width="{1*s}"/>'
    # pill 2
    out+=f'<ellipse cx="{cx+8*s}" cy="{cy+12*s}" rx="{8*s}" ry="{5*s}" fill="#90C4F0" stroke="#5898D0" stroke-width="{1*s}" transform="rotate(-20,{cx+8*s},{cy+12*s})"/>'
    return out

def icon_surgery(cx,cy,scale=1.0):
    """Simplified scissors/scalpel."""
    s=scale
    return (f'<g transform="translate({cx},{cy}) scale({s})">'
            f'<line x1="-12" y1="-12" x2="12" y2="12" stroke="{C["slate"]}" stroke-width="2.5" stroke-linecap="round"/>'
            f'<line x1="12" y1="-12" x2="-12" y2="12" stroke="{C["slate"]}" stroke-width="2.5" stroke-linecap="round"/>'
            f'<circle cx="-12" cy="12" r="5" fill="none" stroke="{C["slate"]}" stroke-width="1.5"/>'
            f'<circle cx="12" cy="12" r="5" fill="none" stroke="{C["slate"]}" stroke-width="1.5"/>'
            f'</g>')

def icon_mrna(cx,cy,w=60,scale=1.0,color=None):
    """mRNA strand with wobbles."""
    col=color or C["sage"]; s=scale; w*=s
    pts=[]
    for i in range(41):
        t=i/40; x=cx-w/2+t*w; y=cy+math.sin(t*math.pi*5)*7*s
        pts.append(f"{x:.1f},{y:.1f}")
    return f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="{2.5*s}" stroke-linecap="round"/>'

def icon_immune_cluster(cx,cy,scale=1.0):
    """Cluster of immune cells (T, B, Treg)."""
    s=scale; out=""
    out+=icon_tcell(cx-20*s,cy-8*s,scale=s*0.65,color=C["good"],label="CD8")
    out+=icon_bcell(cx+22*s,cy-10*s,scale=s*0.6)
    out+=icon_treg(cx+2*s,cy+18*s,scale=s*0.55)
    return out

def icon_tumor_shrink(cx,cy,scale=1.0):
    """Tumor shrinking — before/after comparison."""
    s=scale; out=""
    # before (larger)
    out+=f'<ellipse cx="{cx-20*s}" cy="{cy}" rx="{18*s}" ry="{15*s}" fill="{C["bad"]}" fill-opacity="0.7"/>'
    out+=f'<ellipse cx="{cx-23*s}" cy="{cy-4*s}" rx="{6*s}" ry="{5*s}" fill="{C["badL"]}" fill-opacity="0.4"/>'
    # arrow
    out+=f'<line x1="{cx+2*s}" y1="{cy}" x2="{cx+14*s}" y2="{cy}" stroke="{C["muted"]}" stroke-width="{2*s}" marker-end="url(#arrow-{C["muted"][1:]})"/>'
    # after (smaller, faded)
    out+=f'<ellipse cx="{cx+28*s}" cy="{cy}" rx="{8*s}" ry="{6*s}" fill="{C["bad"]}" fill-opacity="0.25"/>'
    out+=f'<line x1="{cx+22*s}" y1="{cy-8*s}" x2="{cx+34*s}" y2="{cy+8*s}" stroke="{C["sage"]}" stroke-width="{2*s}"/>'
    out+=f'<line x1="{cx+34*s}" y1="{cy-8*s}" x2="{cx+22*s}" y2="{cy+8*s}" stroke="{C["sage"]}" stroke-width="{2*s}"/>'
    return out

# ========== BUILD SVG ==========
parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')

defs=['<defs>']
for col in [C["slate"],C["good"],C["bad"],C["amber"],C["sage"],C["muted"],C["purple"]]:
    defs.append(arrow_marker(col))
defs.append(f'<linearGradient id="goodGrad" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stop-color="{C["goodBg"]}"/><stop offset="100%" stop-color="#FFFFFF"/></linearGradient>')
defs.append(f'<linearGradient id="badGrad" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stop-color="{C["badBg"]}"/><stop offset="100%" stop-color="#FFFFFF"/></linearGradient>')
defs.append(f'<linearGradient id="sageGrad" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stop-color="{C["sageBg"]}"/><stop offset="100%" stop-color="#FFFFFF"/></linearGradient>')
defs.append(f'<filter id="shadow"><feDropShadow dx="0" dy="2" stdDeviation="4" flood-color="#1F2430" flood-opacity="0.10"/></filter>')
defs.append('</defs>')
parts.extend(defs)

parts.append(rect(0,0,W,H,C["bg"],rx=0))

# ---- TITLE ----
parts.append(txt(W/2,50,"Radiation-Phase Molecular Response Predicts TNT Outcome in MSS LARC",
                 size=30,color=C["ink"],weight=700,anchor="middle"))
parts.append(txt(W/2,78,"Integrated WES + RNA-seq of 35 locally advanced rectal cancers (pre / post CRT biopsies)",
                 size=15,color=C["muted"],anchor="middle"))
parts.append(f'<line x1="80" y1="98" x2="{W-80}" y2="98" stroke="{C["line"]}" stroke-width="1"/>')

# ================================================================
# PANEL A — Cohort + sampling (left)
# ================================================================
AX,AY,AW,AH = 50,118,420,450
parts.append(rect(AX,AY,AW,AH,C["soft"],stroke=C["line"],sw=1,rx=20))
parts.append(txt(AX+20,AY+28,"Cohort & Radiation-Phase Sampling",size=15,color=C["slate"],weight=700))

# colon with tumor
parts.append(icon_colon_tumor(AX+90,AY+90,scale=1.6))
parts.append(txt(AX+145,AY+75,"LARC",size=12,color=C["bad"],weight=700))
parts.append(txt(AX+145,AY+92,"tumor",size=11,color=C["muted"]))

# radiation beams hitting tumor
parts.append(icon_radiation(AX+200,AY+80,scale=1.1))
parts.append(txt(AX+200,AY+118,"CRT",size=11,color=C["amber"],weight=600,anchor="middle"))
parts.append(txt(AX+200,AY+132,"50.4 Gy",size=10,color=C["muted"],anchor="middle"))

# patient cohort mini grid (5x7)
pg_x0,pg_y0=AX+280,AY+60
for i in range(35):
    r,c=divmod(i,7)
    cx=pg_x0+c*18; cy=pg_y0+r*22
    good=i<18
    parts.append(icon_person_simple(cx,cy,scale=0.55,color=C["good"] if good else C["bad"]))
parts.append(txt(AX+280,AY+175,"n = 35",size=12,color=C["slate"],weight=600))

# legend
parts.append(f'<circle cx="{AX+270}" cy="{AY+195}" r="5" fill="{C["good"]}"/>')
parts.append(txt(AX+280,AY+199,"18 good (CR/near-CR)",size=10.5,color=C["good"],weight=600))
parts.append(f'<circle cx="{AX+270}" cy="{AY+213}" r="5" fill="{C["bad"]}"/>')
parts.append(txt(AX+280,AY+217,"17 bad (PR/poor)",size=10.5,color=C["bad"],weight=600))

# timeline
ty=AY+268
parts.append(f'<line x1="{AX+24}" y1="{ty}" x2="{AX+AW-24}" y2="{ty}" stroke="{C["slate"]}" stroke-width="2.5"/>')
# radiation band
parts.append(rect(AX+110,ty-14,200,28,C["amberL"],stroke=C["amber"],sw=1.5,rx=8))
parts.append(txt(AX+210,ty+5,"Radiation phase (~5-6 wk)",size=10.5,color=C["amber"],weight=600,anchor="middle"))
# biopsy needles at pre and post
parts.append(icon_biopsy(AX+95,ty-30,angle=-25,scale=0.9))
parts.append(txt(AX+95,ty+22,"pre",size=11,color=C["good"],weight=600,anchor="middle"))
parts.append(f'<circle cx="{AX+95}" cy="{ty}" r="6" fill="{C["good"]}" stroke="#fff" stroke-width="2"/>')

parts.append(icon_biopsy(AX+325,ty-30,angle=-25,scale=0.9))
parts.append(txt(AX+325,ty+22,"post",size=11,color=C["good"],weight=600,anchor="middle"))
parts.append(f'<circle cx="{AX+325}" cy="{ty}" r="6" fill="{C["good"]}" stroke="#fff" stroke-width="2"/>')

# final response star
parts.append(f'<circle cx="{AX+AW-35}" cy="{ty}" r="6" fill="{C["slate"]}" stroke="#fff" stroke-width="2"/>')
parts.append(txt(AX+AW-35,ty+22,"final",size=10,color=C["slate"],anchor="middle"))
parts.append(txt(AX+AW-35,ty+34,"response",size=9,color=C["muted"],anchor="middle"))

# WES + RNA-seq icons
oy=AY+330
parts.append(rect(AX+20,oy,180,48,C["goodBg"],stroke=C["good"],sw=1.2,rx=10))
parts.append(icon_dna_helix(AX+55,oy+24,w=50,h=18,scale=0.9))
parts.append(txt(AX+90,oy+20,"WES",size=13,color=C["good"],weight=700))
parts.append(txt(AX+90,oy+36,"Mutect2 T-N",size=10,color=C["muted"]))

parts.append(rect(AX+215,oy,180,48,C["sageBg"],stroke=C["sage"],sw=1.2,rx=10))
parts.append(icon_mrna(AX+250,oy+24,w=45,scale=0.9,color=C["sage"]))
parts.append(txt(AX+282,oy+20,"RNA-seq",size=13,color=C["sage"],weight=700))
parts.append(txt(AX+282,oy+36,"DESeq2 / GSEA",size=10,color=C["muted"]))

# sample counts
parts.append(txt(AX+20,AY+AH-18,"WES 77 samples  ·  RNA-seq 56  ·  49 PASS VCFs  ·  all MSS",size=10.5,color=C["muted"]))

# ================================================================
# PANEL B — Cascade (center)
# ================================================================
BX,BY,BW,BH = 510,118,610,560
parts.append(rect(BX,BY,BW,BH,C["bg"],stroke=C["line"],sw=1,rx=20))
parts.append(txt(BX+BW/2,BY+28,"Good-Responder Cascade  (pre -> post CRT)",size=17,color=C["slate"],weight=700,anchor="middle"))
parts.append(txt(BX+BW/2,BY+48,"radiation induces sequential tumor clearance & immune reprogramming",size=12,color=C["muted"],anchor="middle"))

# 5 cascade steps — each with big icon + text card
steps_data = [
    {
        "n":"1", "title":"Mutation clearance", "stat":"SBS5 Delta -76 · missense Delta -67",
        "color": C["good"],
        "icon": lambda cx,cy: icon_dna_breaking(cx,cy,scale=1.2),
    },
    {
        "n":"2", "title":"Neoantigen depletion", "stat":"MHC-I binders Delta -312 · P = 0.082",
        "color": C["amber"],
        "icon": lambda cx,cy: (
            icon_neoantigen(cx-16,cy,scale=1.4)
            + icon_neoantigen(cx+10,cy-8,scale=0.9,faded=True)
            + icon_neoantigen(cx+22,cy+6,scale=0.7,faded=True)
        ),
    },
    {
        "n":"3", "title":"HLA-LOH clone removal", "stat":"subj 3 & 4: pre LOH -> post resolution",
        "color": C["bad"],
        "icon": lambda cx,cy: icon_hla_mhc(cx,cy-6,scale=1.3),
    },
    {
        "n":"4", "title":"Immune reprogramming", "stat":"Treg P = 0.026 · MHC-II up · CD8 exhaustion up",
        "color": C["sage"],
        "icon": lambda cx,cy: icon_immune_cluster(cx,cy,scale=1.1),
    },
    {
        "n":"5", "title":"B-cell infiltration", "stat":"IGH/IGK/IGL clone counts up ~+1400",
        "color": C["purple"],
        "icon": lambda cx,cy: (
            icon_bcell(cx-16,cy,scale=0.85)
            + icon_bcell(cx+22,cy-4,scale=0.65)
            + icon_bcell(cx+12,cy+20,scale=0.5)
        ),
    },
]

sy0=BY+72
row_h=96
for i,st in enumerate(steps_data):
    y=sy0+i*row_h
    col=st["color"]
    # number disk
    parts.append(f'<circle cx="{BX+40}" cy="{y+38}" r="22" fill="{col}" filter="url(#shadow)"/>')
    parts.append(txt(BX+40,y+44,st["n"],size=20,color="#FFFFFF",weight=700,anchor="middle"))
    # icon area
    ix=BX+110; iy=y+38
    parts.append(st["icon"](ix,iy))
    # text card
    tcx=BX+180
    parts.append(rect(tcx,y+8,BW-200,68,C["softer"],stroke=col,sw=1.3,rx=10))
    parts.append(txt(tcx+18,y+35,st["title"],size=14.5,color=C["ink"],weight=700))
    parts.append(txt(tcx+18,y+56,st["stat"],size=11.5,color=C["muted"]))
    # arrow to next
    if i<len(steps_data)-1:
        parts.append(arrow_line(BX+40,y+62,BX+40,y+row_h+14,color=col,sw=3))

# tumor shrink mini at top-right of cascade panel
parts.append(icon_tumor_shrink(BX+BW-80,BY+58,scale=0.85))

# ================================================================
# PANEL C — External CD8 validation (right top)
# ================================================================
CX,CY,CW,CH = 1160,118,590,330
parts.append(rect(CX,CY,CW,CH,C["soft"],stroke=C["line"],sw=1,rx=20))
parts.append(txt(CX+20,CY+28,"External CD8-Cytotoxic Axis Validation",size=15,color=C["slate"],weight=700))
parts.append(txt(CX+20,CY+48,"9 nCRT rectal cohorts (N = 721) + Akiyoshi 2023 (N = 298)",size=11.5,color=C["muted"]))

# CD8 T-cell icon
parts.append(icon_tcell(CX+50,CY+95,scale=1.3,color=C["good"],label="CD8"))
parts.append(txt(CX+50,CY+128,"cytotoxic",size=9.5,color=C["good"],weight=600,anchor="middle"))

# forest plot (right of icon)
fx0=CX+160; fxEnd=CX+CW-30
zero_x=fx0+(fxEnd-fx0)*0.25
cohorts=[
    ("GSE45404",   +0.55,0.30),("GSE46862",   +0.70,0.35),
    ("GSE87211",   +0.40,0.28),("GSE94104",   +0.62,0.33),
    ("GSE133057",  +0.35,0.30),("GSE35452",   +0.58,0.32),
    ("GSE119409",  +0.25,0.34),("GSE150082",  +0.32,0.30),
    ("GSE56699",   -0.10,0.30),("Akiyoshi '23",+0.75,0.22),
]
fy0=CY+78; row_h_f=19; scale_f=130
parts.append(f'<line x1="{zero_x}" y1="{fy0-4}" x2="{zero_x}" y2="{fy0+len(cohorts)*row_h_f+10}" stroke="{C["muted"]}" stroke-width="1" stroke-dasharray="3,3"/>')
for i,(name,eff,se) in enumerate(cohorts):
    yy=fy0+i*row_h_f+10
    cx=zero_x+eff*scale_f; lo=zero_x+(eff-1.96*se)*scale_f; hi=zero_x+(eff+1.96*se)*scale_f
    col_f=C["good"] if eff>0 else C["bad"]
    iscol = C["purple"] if "Akiyoshi" in name else col_f
    parts.append(txt(fx0-6,yy+4,name,size=10,color=C["slate"],anchor="end"))
    parts.append(f'<line x1="{lo}" y1="{yy}" x2="{hi}" y2="{yy}" stroke="{iscol}" stroke-width="1.6" stroke-linecap="round"/>')
    if "Akiyoshi" in name:
        parts.append(f'<polygon points="{cx-5},{yy} {cx},{yy-5} {cx+5},{yy} {cx},{yy+5}" fill="{C["purple"]}"/>')
    else:
        parts.append(f'<rect x="{cx-4}" y="{yy-4}" width="8" height="8" fill="{iscol}" stroke="#fff" stroke-width="0.8"/>')

# meta diamond
mety=fy0+len(cohorts)*row_h_f+22
meta_cx=zero_x+0.52*scale_f
parts.append(f'<polygon points="{meta_cx-36},{mety} {meta_cx},{mety-9} {meta_cx+36},{mety} {meta_cx},{mety+9}" fill="{C["slate"]}"/>')
parts.append(txt(fx0-6,mety+4,"Meta (9)",size=10.5,color=C["slate"],anchor="end",weight=700))
parts.append(txt(meta_cx+42,mety+4,"Z = +2.74   P = 0.006",size=11,color=C["slate"],weight=700))

# axis
axy=mety+22
for v,lab in [(-0.5,"-0.5"),(0,"0"),(0.5,"+0.5"),(1.0,"+1.0")]:
    parts.append(txt(zero_x+v*scale_f,axy,lab,size=9.5,color=C["muted"],anchor="middle"))
parts.append(txt((fx0+CX+CW-30)/2,axy+16,"Standardized effect (good vs bad)",size=10.5,color=C["muted"],anchor="middle"))

# ================================================================
# PANEL D — Clinical implication (right bottom)
# ================================================================
DX,DY,DW,DH = 1160,468,590,210
parts.append(rect(DX,DY,DW,DH,"url(#sageGrad)",stroke=C["sage"],sw=1.2,rx=20))
parts.append(txt(DX+20,DY+28,"Mid-Treatment Decision Window",size=15,color=C["sage"],weight=700))
parts.append(txt(DX+20,DY+48,"actionable radiation-phase molecular biomarker",size=11.5,color=C["muted"]))

# central decision box
dbx=DX+DW/2; dby=DY+80
parts.append(rect(dbx-170,dby-16,340,36,C["bg"],stroke=C["sage"],sw=1.4,rx=10))
parts.append(txt(dbx,dby+4,"Post-CRT biopsy molecular profile",size=13,color=C["ink"],weight=700,anchor="middle"))

# left branch — good
parts.append(curved_arrow(dbx-60,dby+22,dbx-140,dby+70,color=C["good"],sw=2.8,bend=-30))
parts.append(rect(DX+30,dby+58,220,94,C["bg"],stroke=C["good"],sw=1.4,rx=12))
parts.append(icon_person_simple(DX+60,dby+90,scale=0.7,color=C["good"]))
parts.append(txt(DX+85,dby+88,"Responder profile",size=12.5,color=C["good"],weight=700))
parts.append(txt(DX+85,dby+106,"-> de-escalate / watch-and-wait",size=10.5,color=C["ink"]))
parts.append(txt(DX+85,dby+122,"organ preservation",size=10.5,color=C["muted"]))

# right branch — bad
parts.append(curved_arrow(dbx+60,dby+22,dbx+140,dby+70,color=C["bad"],sw=2.8,bend=30))
parts.append(rect(DX+DW-260,dby+58,228,94,C["bg"],stroke=C["bad"],sw=1.4,rx=12))
parts.append(icon_pills(DX+DW-234,dby+95,scale=1.0))
parts.append(txt(DX+DW-200,dby+88,"Non-responder profile",size=12.5,color=C["bad"],weight=700))
parts.append(txt(DX+DW-200,dby+106,"-> intensify consolidation",size=10.5,color=C["ink"]))
parts.append(txt(DX+DW-200,dby+122,"FOLFOX / trial enrollment",size=10.5,color=C["muted"]))

# ================================================================
# FLOW ARROWS between panels
# ================================================================
parts.append(curved_arrow(AX+AW,AY+AH/2-60,BX-8,BY+BH/2-50,color=C["muted"],sw=3.5,bend=-20))
parts.append(curved_arrow(BX+BW,BY+BH/3,CX-8,CY+CH/2,color=C["muted"],sw=3.5,bend=-25))
parts.append(curved_arrow(BX+BW,BY+BH*0.75,DX-8,DY+DH/2,color=C["muted"],sw=3.5,bend=25))

# ================================================================
# BOTTOM TAKE-HOME BAR
# ================================================================
bty=698
parts.append(rect(50,bty,W-100,70,C["slate"],rx=16))
parts.append(txt(70,bty+28,
    "Take-home:  Radiation-phase molecular response -- mutation/neoantigen clearance + reproducible CD8-cytotoxic axis --",
    size=15,color="#FFFFFF",weight=600))
parts.append(txt(70,bty+52,
    "defines a mid-treatment biomarker window to personalise consolidation chemotherapy in MSS locally advanced rectal cancer.",
    size=15,color=C["goodL"],weight=400))

# ================================================================
# DECORATIVE — scattered small immune cells around
# ================================================================
# tiny T-cells around cascade panel (floating)
for (dx,dy,sc) in [(BX+BW-40,BY+BH-30,0.35),(BX+BW+15,BY+BH-100,0.3),(BX-25,BY+BH-60,0.3)]:
    parts.append(f'<g opacity="0.25">{icon_tcell(dx,dy,scale=sc,color=C["good"],label="")}</g>')
# tiny B-cells
for (dx,dy,sc) in [(CX+CW-40,CY+CH-20,0.3),(DX+DW-30,DY+DH-20,0.25)]:
    parts.append(f'<g opacity="0.2">{icon_bcell(dx,dy,scale=sc)}</g>')

# neoantigens fading away near step 2
parts.append(f'<g opacity="0.15">{icon_neoantigen(BX+BW-50,BY+180,scale=0.6,faded=True)}</g>')
parts.append(f'<g opacity="0.10">{icon_neoantigen(BX+BW-28,BY+200,scale=0.4,faded=True)}</g>')

parts.append("</svg>")

svg="\n".join(parts)
svg_path=OUT/"GA_v2_master.svg"
svg_path.write_text(svg,encoding="utf-8")

cairosvg.svg2png(bytestring=svg.encode("utf-8"),write_to=str(OUT/"GA_v2.png"),output_width=W*2,output_height=H*2)
cairosvg.svg2pdf(bytestring=svg.encode("utf-8"),write_to=str(OUT/"GA_v2.pdf"))

print("wrote:",svg_path)
print("wrote:",OUT/"GA_v2.png")
print("wrote:",OUT/"GA_v2.pdf")
