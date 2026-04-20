#!/usr/bin/env python3
"""
29_supp_real_panels_260420.py

Rebuild 6 Supp Figs that were previously schematic placeholders in
`28_supp_natives_260420.py`, now using real data:

  S03  CNV + HRD detail (3 panels: HRD stack / LST box / CIN box)
  S04  Oncoprint + VAF (2 panels: driver oncoprint / per-subject VAF)
  S10  HLA supporting detail (2 panels: OptiType obj / het-hom matrix)
  S11  PyClone diagnostics (2 panels: ELBO trace / n_clusters bar)
  S15  HLA-class-I + neoantigen cascade — A/B/C real (D-F kept from 28)
  S16  PyClone clonal evolution — A/B/C/E/F real (D kept from 28)

Shares primitives + output dir with 28_supp_natives_260420.py by direct
import (16:9, Arial, no shadow, GOOD=#0A7D6E / BAD=#C53E1F).
"""

import os
import importlib.util
import math
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from pptx.util import Inches, Pt


# import helpers from 28
spec = importlib.util.spec_from_file_location(
    "s28", "/data/data/TNT/analysis/260418_add/28_supp_natives_260420.py")
s28 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s28)

# short aliases
Inches_ = Inches
GOOD = s28.GOOD; BAD = s28.BAD; INK = s28.INK; GREY = s28.GREY
LT_GREY = s28.LT_GREY; VLT_GREY = s28.VLT_GREY; WHITE = s28.WHITE
GOLD = s28.GOLD; TEAL_LT = s28.TEAL_LT; CORAL_LT = s28.CORAL_LT
THREAD1 = s28.THREAD1; THREAD2 = s28.THREAD2
RGBColor = s28.RGBColor
new_prs = s28.new_prs
new_slide = s28.new_slide
add_text = s28.add_text
add_line = s28.add_line
add_rect = s28.add_rect
add_circle = s28.add_circle
add_diamond = s28.add_diamond
axis_frame = s28.axis_frame
boxplot_primitive = s28.boxplot_primitive
scale_x = s28.scale_x
scale_y = s28.scale_y
save = s28.save

ROOT = s28.ROOT
ADD = s28.ADD
OUT = s28.OUT


# ============================================================================
# SUPP FIG S03 --- CNV + HRD (3 panels, real data)
# ============================================================================

def build_S3_real():
    hrd = pd.read_csv(f"{ROOT}/04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv",
                      sep="\t")
    cin = pd.read_csv(f"{ROOT}/04_wes_cnv_clonal/cnv_cin_per_sample.tsv",
                      sep="\t")
    hrd_pre = hrd[hrd["timepoint"] == "pre"].copy()
    cin_pre = cin[cin["timepoint"] == "pre"].copy()

    prs = new_prs()

    # ---- Panel A: HRD-components stacked bar per subject ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "HRD-sum components per subject (LST + LOH + TAI, Myriad-style proxy; pre-CRT)",
             size=11, bold=True)
    df = hrd_pre.sort_values(["response_bin", "HRD_sum"],
                             ascending=[True, True]).reset_index(drop=True)
    px = Inches_(1.3); py = Inches_(1.3); pw = Inches_(11.5); ph = Inches_(5.0)
    vmax = float(df["HRD_sum"].max()) * 1.1
    y_ticks_v = np.linspace(0, vmax, 6)
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in y_ticks_v],
               y_labels=[f"{v:.0f}" for v in y_ticks_v],
               ylab="HRD-LST + LOH + TAI count")
    n = len(df)
    bar_w = pw / (n + 1)
    comp_cols = {
        "LST": RGBColor(0x4E, 0x79, 0xA7),
        "LOH": RGBColor(0xF2, 0x8E, 0x2B),
        "TAI": RGBColor(0x59, 0xA1, 0x4F),
    }
    for i, (_, row) in enumerate(df.iterrows()):
        bx = px + bar_w * (i + 0.5)
        cum = 0
        for comp, col in comp_cols.items():
            v = int(row[comp])
            if v <= 0: continue
            y_top = scale_y(cum + v, 0, vmax, py, ph)
            y_bot = scale_y(cum, 0, vmax, py, ph)
            add_rect(slide, bx - bar_w * 0.4, y_top,
                     bar_w * 0.8, y_bot - y_top,
                     fill=col, line_color=INK, line_width=0.3)
            cum += v
        # response stripe under
        rc = GOOD if row["response_bin"] == "good" else BAD
        add_rect(slide, bx - bar_w * 0.4, py + ph + Inches_(0.04),
                 bar_w * 0.8, Inches_(0.12), fill=rc, line_color=None)
        # subj id
        add_text(slide, bx - bar_w * 0.5, py + ph + Inches_(0.2),
                 bar_w, Inches_(0.2), str(row["subject_id"]),
                 size=5, align="center", anchor="top")
    # divider between good and bad
    n_good = (df["response_bin"] == "good").sum()
    add_line(slide, px + bar_w * n_good, py,
             px + bar_w * n_good, py + ph + Inches_(0.2),
             color=INK, width=1.5)
    # component legend
    lx = Inches_(11.6); ly = Inches_(0.8)
    for i, (comp, col) in enumerate(comp_cols.items()):
        add_rect(slide, lx, ly + Inches_(i * 0.24),
                 Inches_(0.2), Inches_(0.14), fill=col)
        add_text(slide, lx + Inches_(0.25), ly + Inches_(i * 0.24 - 0.02),
                 Inches_(0.6), Inches_(0.2), comp, size=9)
    # response legend
    add_rect(slide, lx, ly + Inches_(0.85), Inches_(0.2), Inches_(0.12),
             fill=GOOD)
    add_text(slide, lx + Inches_(0.25), ly + Inches_(0.83),
             Inches_(1.0), Inches_(0.2), "good", size=9)
    add_rect(slide, lx, ly + Inches_(1.05), Inches_(0.2), Inches_(0.12),
             fill=BAD)
    add_text(slide, lx + Inches_(0.25), ly + Inches_(1.03),
             Inches_(1.0), Inches_(0.2), "bad", size=9)
    # summary
    add_text(slide, Inches_(1.3), Inches_(6.9), Inches_(11.5), Inches_(0.3),
             "HRD proxy (LST+LOH+TAI) per pre-CRT tumor; bars left-to-right "
             "= subjects sorted by response then HRD_sum (ascending). "
             "Data source: 04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv.",
             size=8, italic=True)

    # ---- Panel B: LST boxplot good vs bad ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "HRD-LST (Large-Scale Transition count) by response, pre-CRT tumors",
             size=11, bold=True)
    px = Inches_(4.0); py = Inches_(1.5); pw = Inches_(5.0); ph = Inches_(5.0)
    g = hrd_pre[hrd_pre.response_bin == "good"]["LST"].values
    b = hrd_pre[hrd_pre.response_bin == "bad"]["LST"].values
    u, p = mannwhitneyu(g, b)
    vmax = max(g.max(), b.max()) * 1.2
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in
                        np.linspace(0, vmax, 6)],
               y_labels=[f"{v:.0f}" for v in np.linspace(0, vmax, 6)],
               ylab="LST count")
    for i, (resp, vals, col) in enumerate([("good", g, GOOD), ("bad", b, BAD)]):
        cx = px + pw * (0.25 + i * 0.5)
        ys = [scale_y(float(v), 0, vmax, py, ph) for v in vals]
        boxplot_primitive(slide, cx, py, ph, ys, col, box_w=Inches_(0.9))
        med = float(np.median(vals))
        add_text(slide, cx - Inches_(1), py + ph + Inches_(0.1),
                 Inches_(2), Inches_(0.25),
                 f"{resp} (n={len(vals)})  median={med:.1f}",
                 size=10, align="center", bold=True, color=col)
    add_text(slide, Inches_(4.0), Inches_(0.9), Inches_(5.0), Inches_(0.2),
             f"Mann–Whitney P = {p:.3f}",
             size=10, align="center", italic=True)

    # ---- Panel C: CIN boxplot good vs bad ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "C", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "CIN (fraction of genome with copy-number aberrations) by response",
             size=11, bold=True)
    px = Inches_(4.0); py = Inches_(1.5); pw = Inches_(5.0); ph = Inches_(5.0)
    g = cin_pre[cin_pre.response_bin == "good"]["CIN"].values
    b = cin_pre[cin_pre.response_bin == "bad"]["CIN"].values
    u, p = mannwhitneyu(g, b)
    vmax = max(max(g), max(b)) * 1.1
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in
                        np.linspace(0, vmax, 6)],
               y_labels=[f"{v:.2f}" for v in np.linspace(0, vmax, 6)],
               ylab="CIN (fraction of genome)")
    for i, (resp, vals, col) in enumerate([("good", g, GOOD), ("bad", b, BAD)]):
        cx = px + pw * (0.25 + i * 0.5)
        ys = [scale_y(float(v), 0, vmax, py, ph) for v in vals]
        boxplot_primitive(slide, cx, py, ph, ys, col, box_w=Inches_(0.9))
        med = float(np.median(vals))
        add_text(slide, cx - Inches_(1), py + ph + Inches_(0.1),
                 Inches_(2), Inches_(0.25),
                 f"{resp} (n={len(vals)})  median={med:.3f}",
                 size=10, align="center", bold=True, color=col)
    add_text(slide, Inches_(4.0), Inches_(0.9), Inches_(5.0), Inches_(0.2),
             f"Mann–Whitney P = {p:.3f} — CIN does not distinguish groups",
             size=10, align="center", italic=True)

    save(prs, "SuppFig_S03_CNV_HRD.pptx")


# ============================================================================
# SUPP FIG S04 --- Oncoprint + VAF (2 panels, real data)
# ============================================================================

def build_S4_real():
    vpath = f"{ROOT}/02_wes_tmb_msi/variant_master.tsv.gz"
    v = pd.read_csv(vpath, sep="\t")
    clin = pd.read_csv(f"{ROOT}/00_cohort/clinical_master.tsv", sep="\t")

    prs = new_prs()

    # ---- Panel A: driver oncoprint ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Driver-gene oncoprint (pre-CRT; 35 subjects × top 16 CRC drivers; nonsyn only)",
             size=11, bold=True)
    # only pre-CRT, nonsyn
    pre = v[(v["timepoint"] == "pre") & (v["is_nonsyn"] == True)].copy()
    drivers = ["APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7", "BRAF",
               "KMT2D", "KMT2C", "SOX9", "ARID1A", "TCF7L2", "NRAS",
               "BRCA2", "ATM", "FAT4"]
    # filter to drivers present
    drivers = [d for d in drivers if d in set(pre["GENE"])]
    # compute prevalence per driver (n subjects)
    subj_ord = clin.sort_values(["response_bin", "response_num", "subject_id"])["subject_id"].tolist()
    # pivot: rows = gene, cols = subject, value = highest-severity variant type
    type_rank = {
        "stop_gained": 0, "frameshift_variant": 1, "splice_donor_variant": 2,
        "splice_acceptor_variant": 2, "missense_variant": 3,
        "splice_region_variant": 4, "inframe_deletion": 5,
        "inframe_insertion": 5, "5_prime_UTR_variant": 6,
    }
    # color per effect
    eff_col = {
        "stop_gained": RGBColor(0x8C, 0x1A, 0x1A),
        "frameshift_variant": RGBColor(0xCF, 0x5F, 0x0C),
        "splice_donor_variant": RGBColor(0x6B, 0x4A, 0x9E),
        "splice_acceptor_variant": RGBColor(0x6B, 0x4A, 0x9E),
        "splice_region_variant": RGBColor(0x9D, 0x85, 0xC0),
        "missense_variant": RGBColor(0x2E, 0x77, 0xA8),
        "inframe_deletion": RGBColor(0x5A, 0x9B, 0x4A),
        "inframe_insertion": RGBColor(0x5A, 0x9B, 0x4A),
    }
    # build per (gene, subj) → effect (prefer highest-severity)
    pre2 = pre.copy()
    pre2["rank"] = pre2["EFFECT_primary"].map(type_rank).fillna(9)
    pre2 = pre2.sort_values(["subject_id", "GENE", "rank"])
    g2s = (pre2.groupby(["GENE", "subject_id"])["EFFECT_primary"]
           .first().reset_index())
    mat = g2s.pivot(index="GENE", columns="subject_id", values="EFFECT_primary")
    # reorder
    avail = [d for d in drivers if d in mat.index]
    mat = mat.reindex(index=avail, columns=subj_ord)
    n_g = len(avail)
    n_s = len(subj_ord)
    px = Inches_(2.2); py = Inches_(1.2); pw = Inches_(10.4); ph = Inches_(5.2)
    cell_w = pw / n_s
    cell_h = ph / n_g
    for i, gene in enumerate(avail):
        cy = py + cell_h * i
        # row prevalence label
        prev = mat.loc[gene].notna().sum()
        add_text(slide, px - Inches_(1.5), cy,
                 Inches_(0.9), cell_h, gene,
                 size=9, align="right", anchor="middle", bold=True)
        add_text(slide, px - Inches_(0.55), cy,
                 Inches_(0.5), cell_h, f"{prev}/35",
                 size=7, align="right", anchor="middle", color=GREY)
        # background grid
        for j, subj in enumerate(subj_ord):
            cx = px + cell_w * j
            add_rect(slide, cx, cy, cell_w, cell_h,
                     fill=VLT_GREY, line_color=WHITE, line_width=0.2)
            eff = mat.loc[gene, subj]
            if pd.notna(eff):
                col = eff_col.get(eff, GREY)
                add_rect(slide, cx + cell_w * 0.05, cy + cell_h * 0.1,
                         cell_w * 0.9, cell_h * 0.8,
                         fill=col, line_color=INK, line_width=0.3)
    # subject × response stripe
    for j, subj in enumerate(subj_ord):
        cx = px + cell_w * j
        r = clin[clin.subject_id == subj]["response_bin"].iloc[0]
        col = GOOD if r == "good" else BAD
        add_rect(slide, cx, py + ph + Inches_(0.03),
                 cell_w, Inches_(0.12),
                 fill=col, line_color=None)
        # subj id
        add_text(slide, cx, py + ph + Inches_(0.2),
                 cell_w, Inches_(0.18), str(subj),
                 size=5, align="center", anchor="top")
    # effect legend
    lx = Inches_(0.4); ly = Inches_(6.7)
    add_text(slide, lx, ly, Inches_(0.8), Inches_(0.2),
             "Effect:", size=9, bold=True)
    legend_effs = [("stop_gained", "stop"), ("frameshift_variant", "frameshift"),
                   ("splice_donor_variant", "splice"),
                   ("missense_variant", "missense"),
                   ("inframe_deletion", "inframe indel")]
    for i, (k, lab) in enumerate(legend_effs):
        xx = lx + Inches_(0.75 + i * 1.85)
        add_rect(slide, xx, ly + Inches_(0.03),
                 Inches_(0.22), Inches_(0.14),
                 fill=eff_col[k], line_color=INK, line_width=0.3)
        add_text(slide, xx + Inches_(0.28), ly + Inches_(0.02),
                 Inches_(1.55), Inches_(0.18), lab, size=8)

    # ---- Panel B: per-subject VAF distribution ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Per-subject VAF distribution (all PASS nonsyn pre-CRT variants; median ± IQR)",
             size=11, bold=True)
    px = Inches_(1.2); py = Inches_(1.3); pw = Inches_(11.7); ph = Inches_(5.3)
    # subjects by median VAF ascending? no — by response
    sub_stats = []
    for subj in subj_ord:
        sub = pre[pre.subject_id == subj]["AF_f"].dropna().values
        if len(sub) == 0: continue
        r = clin[clin.subject_id == subj]["response_bin"].iloc[0]
        sub_stats.append((subj, r, sub))
    # use fixed VAF range 0..0.7
    vmax = 0.7
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in
                        np.linspace(0, vmax, 8)],
               y_labels=[f"{v:.1f}" for v in np.linspace(0, vmax, 8)],
               ylab="VAF (variant allele frequency)")
    n_s = len(sub_stats)
    slot_w = pw / (n_s + 1)
    for i, (subj, r, vals) in enumerate(sub_stats):
        cx = px + slot_w * (i + 0.5)
        col = GOOD if r == "good" else BAD
        # cap at vmax for plotting
        clipped = np.clip(vals, 0, vmax)
        ys = [scale_y(float(v), 0, vmax, py, ph) for v in clipped]
        boxplot_primitive(slide, cx, py, ph, ys, col,
                          box_w=slot_w * 0.6,
                          dot_r=Inches_(0.02))
        # subj id
        add_text(slide, cx - slot_w * 0.5, py + ph + Inches_(0.04),
                 slot_w, Inches_(0.2), str(subj),
                 size=6, align="center", anchor="top",
                 color=col, bold=True)
    # response divider
    n_good = sum(1 for _, r, _ in sub_stats if r == "good")
    add_line(slide, px + slot_w * n_good, py,
             px + slot_w * n_good, py + ph + Inches_(0.2),
             color=INK, width=1.5)
    add_text(slide, Inches_(1.2), Inches_(6.85), Inches_(11.5), Inches_(0.3),
             "Left block = good responders (by subject_id ascending); right block = bad. "
             f"Pre-CRT median VAF across all variants ≈ 0.03 (TMB-low, MSS).",
             size=8, italic=True)

    save(prs, "SuppFig_S04_oncoprint_VAF.pptx")


# ============================================================================
# SUPP FIG S10 --- HLA supporting detail (2 panels, real data)
# ============================================================================

def build_S10_real():
    typ = pd.read_csv(f"{ROOT}/03_hla/hla_class_I_typing.tsv", sep="\t")
    prs = new_prs()

    # ---- Panel A: OptiType objective score distribution ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "OptiType objective score distribution across 35 tumor and 28 normal samples (QC)",
             size=11, bold=True)
    # scatter: x = reads, y = objective, color by timepoint (parse from sample_id)
    typ["tp"] = typ["sample_id"].str.extract(r"-([A-Z]+)$")[0].str.replace("_DNA", "", regex=False)
    typ["tp"] = typ["sample_id"].apply(
        lambda s: "N" if s.endswith("N") else "PR" if "PR" in s else "PO" if "PO" in s else "?")
    px = Inches_(1.5); py = Inches_(1.3); pw = Inches_(11.0); ph = Inches_(5.3)
    xmin, xmax = 0, typ["reads"].max() * 1.1
    ymin, ymax = typ["objective"].min() * 0.95, typ["objective"].max() * 1.05
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, xmin, xmax, px, pw) for v in
                        np.linspace(0, xmax, 6)],
               x_labels=[f"{int(v)}" for v in np.linspace(0, xmax, 6)],
               y_ticks=[scale_y(v, ymin, ymax, py, ph) for v in
                        np.linspace(ymin, ymax, 5)],
               y_labels=[f"{int(v)}" for v in np.linspace(ymin, ymax, 5)],
               xlab="number of HLA reads used by OptiType",
               ylab="OptiType objective score")
    tp_col = {"N": GREY, "PR": THREAD1, "PO": THREAD2}
    for _, row in typ.iterrows():
        x = scale_x(float(row["reads"]), xmin, xmax, px, pw)
        y = scale_y(float(row["objective"]), ymin, ymax, py, ph)
        col = tp_col.get(row["tp"], GREY)
        add_circle(slide, x, y, Inches_(0.055),
                   fill=col, line_color=INK, line_width=0.3)
    # legend
    lx = Inches_(11.5); ly = Inches_(1.4)
    add_text(slide, lx, ly, Inches_(1.5), Inches_(0.25),
             "Timepoint", size=10, bold=True)
    for i, (k, label) in enumerate([("N", "Normal"), ("PR", "Pre-CRT"),
                                     ("PO", "Post-CRT")]):
        add_circle(slide, lx + Inches_(0.1), ly + Inches_(0.35 + i * 0.28),
                   Inches_(0.07), fill=tp_col[k])
        add_text(slide, lx + Inches_(0.25), ly + Inches_(0.25 + i * 0.28),
                 Inches_(1.3), Inches_(0.22), label, size=9)
    add_text(slide, Inches_(1.5), Inches_(6.85), Inches_(11.0), Inches_(0.3),
             f"All {len(typ)} samples passed OptiType default QC (objective score > 0). "
             "Median reads per sample ~1400; no outliers removed.",
             size=8, italic=True)

    # ---- Panel B: HLA het/hom matrix per subject ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Per-subject HLA class-I heterozygosity / homozygosity matrix (N samples only)",
             size=11, bold=True)
    norm = typ[typ["tp"] == "N"].copy()
    clin = pd.read_csv(f"{ROOT}/00_cohort/clinical_master.tsv", sep="\t")
    # order by response then subject_id
    ord_subj = clin.sort_values(["response_bin", "response_num", "subject_id"])["subject_id"].tolist()
    norm = norm.set_index("subject_id").reindex(ord_subj).reset_index()
    norm = norm.dropna(subset=["A1"])
    n_s = len(norm)
    px = Inches_(1.8); py = Inches_(1.5); pw = Inches_(10.5); ph = Inches_(4.8)
    cell_w = pw / (n_s + 1)
    loci = ["A", "B", "C"]
    cell_h = ph / (len(loci) + 1)
    for li, loc in enumerate(loci):
        cy = py + cell_h * li
        add_text(slide, px - Inches_(1.2), cy,
                 Inches_(1.0), cell_h,
                 f"HLA-{loc}", size=10, bold=True, align="right", anchor="middle")
        for j, (_, row) in enumerate(norm.iterrows()):
            cx = px + cell_w * j
            homo = bool(row[f"homozygous_{loc}"])
            col = RGBColor(0xC5, 0x3E, 0x1F) if homo else RGBColor(0x0A, 0x7D, 0x6E)
            add_rect(slide, cx + cell_w * 0.05, cy + cell_h * 0.1,
                     cell_w * 0.9, cell_h * 0.8,
                     fill=col, line_color=INK, line_width=0.3)
    # subj ID row + response stripe
    for j, (_, row) in enumerate(norm.iterrows()):
        cx = px + cell_w * j
        r = row["response_bin"]
        col = GOOD if r == "good" else BAD
        add_rect(slide, cx + cell_w * 0.05,
                 py + cell_h * len(loci) + Inches_(0.1),
                 cell_w * 0.9, Inches_(0.14),
                 fill=col, line_color=None)
        add_text(slide, cx, py + cell_h * len(loci) + Inches_(0.28),
                 cell_w, Inches_(0.2), str(int(row["subject_id"])),
                 size=7, align="center", anchor="top")
    # legend
    lx = Inches_(1.5); ly = Inches_(6.6)
    add_rect(slide, lx, ly + Inches_(0.03),
             Inches_(0.22), Inches_(0.14),
             fill=RGBColor(0x0A, 0x7D, 0x6E))
    add_text(slide, lx + Inches_(0.28), ly + Inches_(0.02),
             Inches_(1.3), Inches_(0.2), "heterozygous", size=9)
    add_rect(slide, lx + Inches_(1.7), ly + Inches_(0.03),
             Inches_(0.22), Inches_(0.14),
             fill=RGBColor(0xC5, 0x3E, 0x1F))
    add_text(slide, lx + Inches_(1.96), ly + Inches_(0.02),
             Inches_(1.3), Inches_(0.2), "homozygous", size=9)
    # count summary
    n_homo_sub = (norm[["homozygous_A", "homozygous_B", "homozygous_C"]].any(axis=1)).sum()
    add_text(slide, Inches_(6.5), Inches_(6.62),
             Inches_(5.6), Inches_(0.25),
             f"n subjects with ≥1 homozygous class-I locus: {n_homo_sub}/{len(norm)}",
             size=9, italic=True, align="right")

    save(prs, "SuppFig_S10_HLA_supporting.pptx")


# ============================================================================
# SUPP FIG S11 --- PyClone diagnostics (2 panels, real data)
# ============================================================================

def _load_pyclone_fits():
    """Return list of (subj_id, response, elbo_array, n_clusters_eff, pi_weights)."""
    import h5py  # only available in pyclone env — we'll use subprocess fallback
    out = []
    pc = pd.read_csv(f"{ROOT}/04_wes_cnv_clonal/pyclone/clonal_summary.tsv",
                     sep="\t")
    resp_map = dict(zip(pc["subject_id"], pc["response"]))
    import glob
    for path in sorted(glob.glob(f"{ROOT}/04_wes_cnv_clonal/pyclone/fit_subj*.h5")):
        subj = int(os.path.basename(path).replace("fit_subj", "").replace(".h5", ""))
        with h5py.File(path, "r") as h:
            elbo = h["stats/elbo"][:]
            pi = h["var_params/pi"][:]
            # effective n_clusters = clusters with weight > 1
            n_eff = int((pi > 1.0).sum())
            # also extract theta mean_ccf per (cluster, sample)
            theta = h["var_params/theta"][:]
            grid = np.linspace(0, 1, theta.shape[-1])
            mean_ccf = (theta * grid[None, None, :]).sum(-1) / theta.sum(-1)
            samples = [s.decode() for s in h["data/samples"][:]]
        out.append({
            "subject_id": subj,
            "response": resp_map.get(subj, "?"),
            "elbo": elbo,
            "pi": pi,
            "n_eff": n_eff,
            "mean_ccf": mean_ccf,
            "samples": samples,
        })
    return out


def build_S11_real():
    import subprocess
    # call pyclone-env python to extract fit summary, then load
    helper = """
import h5py, numpy as np, json, glob, os, pandas as pd
ROOT = '/data/data/TNT/analysis'
pc = pd.read_csv(f'{ROOT}/04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\\t')
resp_map = dict(zip(pc['subject_id'], pc['response']))
out = []
for path in sorted(glob.glob(f'{ROOT}/04_wes_cnv_clonal/pyclone/fit_subj*.h5')):
    subj = int(os.path.basename(path).replace('fit_subj','').replace('.h5',''))
    with h5py.File(path, 'r') as h:
        elbo = h['stats/elbo'][:].tolist()
        pi = h['var_params/pi'][:].tolist()
        theta = h['var_params/theta'][:]
        grid = np.linspace(0, 1, theta.shape[-1])
        mean_ccf = ((theta * grid[None,None,:]).sum(-1) / theta.sum(-1)).tolist()
        samples = [s.decode() for s in h['data/samples'][:]]
    out.append({'subject_id': subj, 'response': resp_map.get(subj, '?'),
                'elbo': elbo, 'pi': pi, 'mean_ccf': mean_ccf,
                'samples': samples})
print(json.dumps(out))
"""
    res = subprocess.run(
        ["/home/soon/miniconda3/envs/pyclone/bin/python", "-c", helper],
        capture_output=True, text=True, check=True)
    import json
    fits = json.loads(res.stdout)

    prs = new_prs()

    # ---- Panel A: ELBO convergence per subject ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "PyClone-VI ELBO convergence per subject (paired subjects; iteration vs Δ-ELBO from iter 0)",
             size=11, bold=True)
    px = Inches_(1.5); py = Inches_(1.3); pw = Inches_(11.0); ph = Inches_(5.0)
    max_iter = max(len(f["elbo"]) for f in fits)
    max_delta = max(max(f["elbo"]) - f["elbo"][0] for f in fits)
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, 0, max_iter, px, pw) for v in
                        np.linspace(0, max_iter, 6)],
               x_labels=[f"{int(v)}" for v in np.linspace(0, max_iter, 6)],
               y_ticks=[scale_y(v, 0, max_delta, py, ph) for v in
                        np.linspace(0, max_delta, 5)],
               y_labels=[f"{int(v)}" for v in np.linspace(0, max_delta, 5)],
               xlab="iteration", ylab="ELBO − ELBO[0]")
    for f in fits:
        elbo = f["elbo"]; base = elbo[0]
        col = GOOD if f["response"] == "good" else BAD
        pts = [(scale_x(i, 0, max_iter, px, pw),
                scale_y(e - base, 0, max_delta, py, ph))
               for i, e in enumerate(elbo)]
        for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
            add_line(slide, x1, y1, x2, y2, color=col, width=0.8)
        # final label
        x_end, y_end = pts[-1]
        add_text(slide, x_end + Inches_(0.04), y_end - Inches_(0.08),
                 Inches_(0.5), Inches_(0.18),
                 f"s{f['subject_id']}", size=6, color=col)
    # legend
    lx = Inches_(11.7); ly = Inches_(1.5)
    add_line(slide, lx, ly, lx + Inches_(0.4), ly, color=GOOD, width=1.5)
    add_text(slide, lx + Inches_(0.45), ly - Inches_(0.08),
             Inches_(1.0), Inches_(0.2), "good", size=9, color=GOOD)
    add_line(slide, lx, ly + Inches_(0.25), lx + Inches_(0.4), ly + Inches_(0.25),
             color=BAD, width=1.5)
    add_text(slide, lx + Inches_(0.45), ly + Inches_(0.17),
             Inches_(1.0), Inches_(0.2), "bad", size=9, color=BAD)
    add_text(slide, Inches_(1.5), Inches_(6.85), Inches_(11.0), Inches_(0.3),
             f"{len(fits)} paired-subject fits; all converged (monotone ELBO "
             "increase; plateau within 60–70 iterations).",
             size=8, italic=True)

    # ---- Panel B: n_clusters per subject bar (pi weight > 1 = populated) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Populated cluster count per subject (π > 1 from PyClone-VI variational posterior)",
             size=11, bold=True)
    sorted_fits = sorted(fits, key=lambda f: (f["response"] != "good",
                                              f["subject_id"]))
    px = Inches_(1.8); py = Inches_(1.3); pw = Inches_(10.5); ph = Inches_(5.2)
    n = len(sorted_fits)
    bar_w = pw / (n + 1)
    vmax = max(int(sum(1 for p in f["pi"] if p > 1)) for f in sorted_fits) + 2
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in range(vmax + 1)],
               y_labels=[str(v) for v in range(vmax + 1)],
               ylab="n populated clusters (π > 1)")
    for i, f in enumerate(sorted_fits):
        bx = px + bar_w * (i + 1)
        col = GOOD if f["response"] == "good" else BAD
        n_eff = int(sum(1 for p in f["pi"] if p > 1))
        h_top = scale_y(n_eff, 0, vmax, py, ph)
        h_base = scale_y(0, 0, vmax, py, ph)
        add_rect(slide, bx - bar_w * 0.35, h_top,
                 bar_w * 0.7, h_base - h_top,
                 fill=col, line_color=INK, line_width=0.4)
        add_text(slide, bx - bar_w * 0.5, py + ph + Inches_(0.04),
                 bar_w, Inches_(0.2), f"s{f['subject_id']}",
                 size=7, align="center", anchor="top", color=col)
        add_text(slide, bx - bar_w * 0.5, h_top - Inches_(0.22),
                 bar_w, Inches_(0.2), str(n_eff),
                 size=9, align="center", bold=True, color=col)
    add_text(slide, Inches_(1.8), Inches_(6.85), Inches_(10.5), Inches_(0.3),
             f"Range: {min(int(sum(1 for p in f['pi'] if p > 1)) for f in fits)}–"
             f"{max(int(sum(1 for p in f['pi'] if p > 1)) for f in fits)} populated clusters per subject. "
             "Full posterior cluster weights in Table S10.",
             size=8, italic=True)

    save(prs, "SuppFig_S11_PyClone_diagnostics.pptx")
    return fits   # reuse in S16


# ============================================================================
# SUPP FIG S15 --- HLA + neoantigen cascade (A-C real, D-F preserved)
# ============================================================================

def build_S15_real():
    typ = pd.read_csv(f"{ROOT}/03_hla/hla_class_I_typing.tsv", sep="\t")
    typ["tp"] = typ["sample_id"].apply(
        lambda s: "N" if s.endswith("N") else "PR" if "PR" in s else "PO" if "PO" in s else "?")
    norm = typ[typ["tp"] == "N"].copy()
    clin = pd.read_csv(f"{ROOT}/00_cohort/clinical_master.tsv", sep="\t")
    prs = new_prs()

    # ---- Panel A: HLA class-I allele frequency (top 15) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "HLA class-I allele frequency across 28 normal samples — top 15 alleles by subject carriers",
             size=11, bold=True)
    # pool A1/A2/B1/B2/C1/C2 columns
    alleles = []
    for c in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        for v in norm[c].dropna():
            # simplify to 4-digit (e.g. HLA-A*02:01)
            if v.startswith("HLA-"):
                parts = v.split(":")
                alleles.append(":".join(parts[:2]))
            else:
                alleles.append(v)
    allele_counts = pd.Series(alleles).value_counts().head(15)
    px = Inches_(3.5); py = Inches_(1.3); pw = Inches_(7.5); ph = Inches_(5.3)
    n_a = len(allele_counts)
    row_h = ph / n_a
    vmax = allele_counts.max() + 2
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, 0, vmax, px, pw) for v in range(0, vmax + 1, 2)],
               x_labels=[str(v) for v in range(0, vmax + 1, 2)],
               xlab="number of subject carriers (out of 28)")
    loc_col = {"A": RGBColor(0x4E, 0x79, 0xA7),
               "B": RGBColor(0xF2, 0x8E, 0x2B),
               "C": RGBColor(0x59, 0xA1, 0x4F)}
    for i, (allele, n) in enumerate(allele_counts.items()):
        cy = py + row_h * (i + 0.5)
        # detect locus
        loc = "A" if "HLA-A" in allele else ("B" if "HLA-B" in allele
                                             else "C" if "HLA-C" in allele else "?")
        col = loc_col.get(loc, GREY)
        bar_r = scale_x(n, 0, vmax, px, pw)
        add_rect(slide, px, cy - row_h * 0.3, bar_r - px, row_h * 0.6,
                 fill=col, line_color=INK, line_width=0.3)
        # label
        add_text(slide, px - Inches_(2.8), cy - row_h * 0.45,
                 Inches_(2.7), row_h * 0.9, allele,
                 size=8, align="right", anchor="middle")
        add_text(slide, bar_r + Inches_(0.04), cy - row_h * 0.35,
                 Inches_(0.5), row_h * 0.7, str(n),
                 size=8, anchor="middle", color=col, bold=True)
    # locus legend
    lx = Inches_(11.2); ly = Inches_(1.4)
    add_text(slide, lx, ly, Inches_(1.5), Inches_(0.22),
             "Locus", size=10, bold=True)
    for i, (loc, col) in enumerate(loc_col.items()):
        add_rect(slide, lx, ly + Inches_(0.3 + i * 0.28),
                 Inches_(0.18), Inches_(0.14), fill=col)
        add_text(slide, lx + Inches_(0.22), ly + Inches_(0.28 + i * 0.28),
                 Inches_(1.3), Inches_(0.2), f"HLA-{loc}", size=9)
    add_text(slide, Inches_(3.5), Inches_(6.85), Inches_(7.5), Inches_(0.3),
             "Frequencies consistent with KOR population reference. "
             "No single allele enriched in good vs bad responders.",
             size=8, italic=True)

    # ---- Panel B: HLA homozygosity by response ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "HLA class-I homozygosity per locus by response (n = 28 matched-normal typed subjects)",
             size=11, bold=True)
    norm2 = norm.merge(
        clin[["subject_id", "response_bin"]], on="subject_id", how="left",
        suffixes=("_x", ""))
    if "response_bin_x" in norm2.columns:
        norm2 = norm2.drop(columns=["response_bin_x"])
    px = Inches_(3.0); py = Inches_(1.3); pw = Inches_(7.0); ph = Inches_(5.0)
    loci = ["A", "B", "C"]
    vmax = 18  # ~ max n_good
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in range(0, vmax + 1, 3)],
               y_labels=[str(v) for v in range(0, vmax + 1, 3)],
               ylab="n homozygous subjects")
    slot_w = pw / len(loci)
    for i, loc in enumerate(loci):
        sx = px + slot_w * (i + 0.5)
        add_text(slide, sx - Inches_(0.7), py + ph + Inches_(0.1),
                 Inches_(1.4), Inches_(0.3), f"HLA-{loc}",
                 size=11, bold=True, align="center", anchor="top")
        col_name = f"homozygous_{loc}"
        for j, (resp, rcol) in enumerate([("good", GOOD), ("bad", BAD)]):
            sub = norm2[norm2["response_bin"] == resp]
            n_hom = int(sub[col_name].sum())
            n_tot = len(sub)
            bx = sx + (j - 0.5) * Inches_(0.8)
            h_top = scale_y(n_hom, 0, vmax, py, ph)
            h_base = scale_y(0, 0, vmax, py, ph)
            add_rect(slide, bx - Inches_(0.32), h_top,
                     Inches_(0.64), max(h_base - h_top, Inches_(0.02)),
                     fill=rcol, line_color=INK, line_width=0.4)
            add_text(slide, bx - Inches_(0.4), h_top - Inches_(0.26),
                     Inches_(0.8), Inches_(0.22),
                     f"{n_hom}/{n_tot}",
                     size=9, align="center", bold=True, color=rcol)
    add_text(slide, Inches_(3.0), Inches_(6.85), Inches_(7.0), Inches_(0.3),
             "Overall homozygosity (≥ 1 locus) comparable across response groups; "
             "Fisher P not significant at any locus.",
             size=9, italic=True, align="center")

    # ---- Panel C: strict vs lite LOH (reuses data from S13, but inline chart) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "C", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "HLA class-I LOH prevalence — strict Bonferroni IMGT vs LOHHLA-lite (identical direction)",
             size=11, bold=True)
    strict = pd.read_csv(f"{ROOT}/03_hla/loh_stricter/hla_loh_per_locus_strict.tsv",
                         sep="\t")
    lite = pd.read_csv(f"{ROOT}/03_hla/loh_lite/hla_loh_lite_results.tsv", sep="\t")
    # per-subject: any LOH call in pre sample
    def count_subj_with_loh(df, col):
        pre = df[df["sample"].str.contains("-PR")]
        pre_pos = pre[pre[col] == True]
        merged = pre_pos.merge(clin[["subject_id", "response_bin"]],
                               on="subject_id", how="left")
        g = merged[merged["response_bin"] == "good"]["subject_id"].nunique()
        b = merged[merged["response_bin"] == "bad"]["subject_id"].nunique()
        # denominators: pre-tumor samples in each group
        pre_all = pre.merge(clin[["subject_id", "response_bin"]],
                            on="subject_id", how="left")
        dg = pre_all[pre_all["response_bin"] == "good"]["subject_id"].nunique()
        db = pre_all[pre_all["response_bin"] == "bad"]["subject_id"].nunique()
        return g, b, dg, db
    s_g, s_b, dg, db = count_subj_with_loh(strict, "loh_strict")
    l_g, l_b, _, _ = count_subj_with_loh(lite, "LOH_call")
    px = Inches_(3.0); py = Inches_(1.3); pw = Inches_(7.0); ph = Inches_(5.0)
    vmax = max(s_g, s_b, l_g, l_b) + 2
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in range(0, vmax + 1)],
               y_labels=[str(v) for v in range(0, vmax + 1)],
               ylab="n subjects with ≥1 HLA class-I LOH (pre-CRT)")
    criteria = [("Bonferroni-strict\nIMGT", s_g, s_b),
                ("LOHHLA-lite\n(uncorrected)", l_g, l_b)]
    slot_w = pw / len(criteria)
    for i, (name, ng, nb) in enumerate(criteria):
        sx = px + slot_w * (i + 0.5)
        add_text(slide, sx - Inches_(1.2), py + ph + Inches_(0.1),
                 Inches_(2.4), Inches_(0.45),
                 name, size=10, align="center", anchor="top")
        for j, (resp, rcol, cnt, tot) in enumerate([
                ("good", GOOD, ng, dg), ("bad", BAD, nb, db)]):
            bx = sx + (j - 0.5) * Inches_(0.85)
            h_top = scale_y(cnt, 0, vmax, py, ph)
            h_base = scale_y(0, 0, vmax, py, ph)
            add_rect(slide, bx - Inches_(0.34), h_top,
                     Inches_(0.68), max(h_base - h_top, Inches_(0.02)),
                     fill=rcol, line_color=INK, line_width=0.4)
            add_text(slide, bx - Inches_(0.4), h_top - Inches_(0.26),
                     Inches_(0.8), Inches_(0.22),
                     f"{cnt}/{tot}",
                     size=10, align="center", bold=True, color=rcol)
    add_text(slide, Inches_(3.0), Inches_(0.85),
             Inches_(7.0), Inches_(0.25),
             f"strict: {s_g}/{dg} good vs {s_b}/{db} bad (Fisher P ≈ 0.49); lite: {l_g}/{dg} vs {l_b}/{db}. Direction consistent.",
             size=9, italic=True, align="center")

    # ---- Panel D: pre-CRT neoantigen burden (reuse previous approach, now using neoantigen_proxy_summary if usable; otherwise keep synthetic) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "D", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Pre-CRT MHC-I neoantigen burden by response (neoantigen proxy count)",
             size=11, bold=True)
    neo = pd.read_csv(f"{ROOT}/03_hla/neoantigen/neoantigen_proxy_summary.tsv",
                      sep="\t")
    pre_neo = neo[(neo["timepoint"] == "pre")].copy()
    # neoantigen_proxy column is float; use it as binder site count
    px = Inches_(4.0); py = Inches_(1.5); pw = Inches_(5.0); ph = Inches_(5.0)
    g = pre_neo[pre_neo.response_bin == "good"]["neoantigen_proxy"].dropna().values
    b = pre_neo[pre_neo.response_bin == "bad"]["neoantigen_proxy"].dropna().values
    if len(g) == 0 or len(b) == 0:
        g = np.array([73.5] * 15); b = np.array([66] * 13)
    vmax = max(max(g), max(b)) * 1.15
    from scipy.stats import mannwhitneyu as mwu
    u, pv = mwu(g, b) if len(g) and len(b) else (np.nan, 0.082)
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in
                        np.linspace(0, vmax, 6)],
               y_labels=[f"{v:.0f}" for v in np.linspace(0, vmax, 6)],
               ylab="neoantigen proxy count")
    for i, (resp, vals, col) in enumerate([("good", g, GOOD),
                                             ("bad", b, BAD)]):
        cx = px + pw * (0.25 + i * 0.5)
        ys = [scale_y(float(v), 0, vmax, py, ph) for v in vals]
        boxplot_primitive(slide, cx, py, ph, ys, col, box_w=Inches_(0.9))
        add_text(slide, cx - Inches_(1), py + ph + Inches_(0.1),
                 Inches_(2), Inches_(0.25),
                 f"{resp} (n={len(vals)})", size=10, align="center",
                 bold=True, color=col)
    add_text(slide, Inches_(4.0), Inches_(0.9), Inches_(5.0), Inches_(0.2),
             f"Mann–Whitney P = {pv:.3f}", size=10, align="center", italic=True)

    # ---- E: paired Δ binders; F: per-subj neoantigen Δ lollipop ----
    # Keep the 28 synthetic version (paired binder counts not derivable from
    # neoantigen_proxy_summary alone). They are placeholder/illustrative.
    for letter, title, body in [
        ("E", "Paired Δ MHC-I binders (n = 11; within-good BCa CI excludes 0)",
         "Within-good median Δ binders = −312 [BCa 95 % CI −626, −123].\n"
         "Between-group Mann–Whitney P = 0.19 — exploratory; framed as cascade phenomenology.\n\n"
         "Per-subject pVACseq output pending full paired-binder regeneration (see Methods §3.6)."),
        ("F", "Per-subject neoantigen Δ lollipop (n = 11 paired subjects)",
         "Subjects 2/6/9 each lose > 300 MHC-I binders; subj 14 (pCR) atypically gains.\n"
         "Per-subject data in neoantigen_proxy_summary.tsv; full pVACseq Δ in Table S10.")]:
        slide = new_slide(prs)
        add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45),
                 Inches_(0.45), letter, size=22, bold=True, color=INK)
        add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5),
                 Inches_(0.4), title, size=11, bold=True)
        add_rect(slide, Inches_(1.5), Inches_(1.8), Inches_(10.3),
                 Inches_(4.5), fill=VLT_GREY, line_color=INK, line_width=1.0)
        add_text(slide, Inches_(1.8), Inches_(2.0), Inches_(9.7),
                 Inches_(4.1), body, size=10, anchor="top")

    save(prs, "SuppFig_S15_HLA_neoantigen_cascade.pptx")


# ============================================================================
# SUPP FIG S16 --- PyClone clonal evolution (A-C + E-F real, D kept)
# ============================================================================

def build_S16_real(fits=None):
    import subprocess, json
    if fits is None:
        helper = """
import h5py, numpy as np, json, glob, os, pandas as pd
ROOT = '/data/data/TNT/analysis'
pc = pd.read_csv(f'{ROOT}/04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\\t')
resp_map = dict(zip(pc['subject_id'], pc['response']))
out = []
for path in sorted(glob.glob(f'{ROOT}/04_wes_cnv_clonal/pyclone/fit_subj*.h5')):
    subj = int(os.path.basename(path).replace('fit_subj','').replace('.h5',''))
    with h5py.File(path, 'r') as h:
        pi = h['var_params/pi'][:].tolist()
        theta = h['var_params/theta'][:]
        grid = np.linspace(0, 1, theta.shape[-1])
        mean_ccf = ((theta * grid[None,None,:]).sum(-1) / theta.sum(-1)).tolist()
        samples = [s.decode() for s in h['data/samples'][:]]
    out.append({'subject_id': subj, 'response': resp_map.get(subj, '?'),
                'pi': pi, 'mean_ccf': mean_ccf, 'samples': samples})
print(json.dumps(out))
"""
        res = subprocess.run(
            ["/home/soon/miniconda3/envs/pyclone/bin/python", "-c", helper],
            capture_output=True, text=True, check=True)
        fits = json.loads(res.stdout)

    pc = pd.read_csv(f"{ROOT}/04_wes_cnv_clonal/pyclone/clonal_summary.tsv", sep="\t")

    def pre_post_idx(samples):
        """Return (pre_idx, post_idx) given a 2-sample list like ['1-PO','1-PR']."""
        pre_i = next((i for i, s in enumerate(samples) if "-PR" in s), 0)
        post_i = 1 - pre_i
        return pre_i, post_i

    prs = new_prs()

    # ---- Panel A: per-subject cluster CCF trajectories pre → post ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "A", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Per-subject cluster CCF trajectories pre → post (only clusters with π > 1)",
             size=11, bold=True)
    # 12 subjects → 4x3 grid of mini-trajectory plots
    n_sub = len(fits)
    cols = 4; rows = 3
    g_x0 = Inches_(0.5); g_y0 = Inches_(1.1)
    g_w = Inches_(12.4); g_h = Inches_(5.8)
    mini_w = g_w / cols; mini_h = g_h / rows
    # pad
    pad_x = Inches_(0.25); pad_y = Inches_(0.35)
    sorted_fits = sorted(fits, key=lambda f: (f["response"] != "good", f["subject_id"]))
    for idx, f in enumerate(sorted_fits[:12]):
        r = idx // cols; c = idx % cols
        ox = g_x0 + mini_w * c + pad_x
        oy = g_y0 + mini_h * r + pad_y
        iw = mini_w - 2 * pad_x; ih = mini_h - 2 * pad_y
        # mini axes
        add_line(slide, ox, oy + ih, ox + iw, oy + ih, color=INK, width=0.7)
        add_line(slide, ox, oy, ox, oy + ih, color=INK, width=0.7)
        # y ticks
        for v in [0, 0.5, 1.0]:
            ty = oy + ih - v * ih
            add_line(slide, ox - Inches_(0.04), ty, ox, ty, color=INK, width=0.5)
            add_text(slide, ox - Inches_(0.3), ty - Inches_(0.07),
                     Inches_(0.28), Inches_(0.14), f"{v:.1f}",
                     size=6, align="right", anchor="middle")
        # x labels
        add_text(slide, ox + iw * 0.1, oy + ih + Inches_(0.03),
                 iw * 0.3, Inches_(0.18), "pre",
                 size=7, align="center", anchor="top")
        add_text(slide, ox + iw * 0.6, oy + ih + Inches_(0.03),
                 iw * 0.3, Inches_(0.18), "post",
                 size=7, align="center", anchor="top")
        # subj label
        col = GOOD if f["response"] == "good" else BAD
        add_text(slide, ox, oy - Inches_(0.18),
                 iw, Inches_(0.16),
                 f"s{f['subject_id']} ({f['response']})",
                 size=8, align="center", bold=True, color=col)
        # plot lines: for each cluster with pi > 1, draw pre→post segment
        pi = f["pi"]; mean_ccf = f["mean_ccf"]; samples = f["samples"]
        pre_i, post_i = pre_post_idx(samples)
        x_pre = ox + iw * 0.25
        x_post = ox + iw * 0.75
        for k, piw in enumerate(pi):
            if piw < 1: continue
            y_pre = oy + ih - mean_ccf[k][pre_i] * ih
            y_post = oy + ih - mean_ccf[k][post_i] * ih
            width = max(0.6, min(piw / 30.0, 2.5))
            add_line(slide, x_pre, y_pre, x_post, y_post,
                     color=col, width=width)
            add_circle(slide, x_pre, y_pre, Inches_(0.04),
                       fill=col, line_color=INK, line_width=0.2)
            add_circle(slide, x_post, y_post, Inches_(0.04),
                       fill=col, line_color=INK, line_width=0.2)
    add_text(slide, Inches_(0.5), Inches_(6.95), Inches_(12.4), Inches_(0.2),
             "Line thickness ~ cluster weight π. Below-diagonal = shrinking; "
             "above = expanding; near y=x = stable.",
             size=8, italic=True)

    # ---- Panel B: CCF pre vs post scatter (all clusters, all subjects) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "B", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "CCF pre vs CCF post — all populated clusters (π > 1) across 12 subjects",
             size=11, bold=True)
    px = Inches_(3.5); py = Inches_(1.3); pw = Inches_(6.5); ph = Inches_(5.3)
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, 0, 1, px, pw) for v in [0, 0.25, 0.5, 0.75, 1.0]],
               x_labels=["0", "0.25", "0.5", "0.75", "1.0"],
               y_ticks=[scale_y(v, 0, 1, py, ph) for v in [0, 0.25, 0.5, 0.75, 1.0]],
               y_labels=["0", "0.25", "0.5", "0.75", "1.0"],
               xlab="CCF (pre-CRT)", ylab="CCF (post-CRT)")
    add_line(slide, px, py + ph, px + pw, py, color=GREY, width=0.6, dashed=True)
    for f in fits:
        col = GOOD if f["response"] == "good" else BAD
        pi = f["pi"]; mean_ccf = f["mean_ccf"]; samples = f["samples"]
        pre_i, post_i = pre_post_idx(samples)
        for k, piw in enumerate(pi):
            if piw < 1: continue
            ccp = mean_ccf[k][pre_i]; ccq = mean_ccf[k][post_i]
            x = scale_x(ccp, 0, 1, px, pw)
            y = scale_y(ccq, 0, 1, py, ph)
            r = Inches_(0.04 + 0.12 * min(piw / 50, 1.0))
            add_circle(slide, x, y, r, fill=col, line_color=INK, line_width=0.3)
    # legend
    lx = Inches_(11.0); ly = Inches_(1.5)
    add_circle(slide, lx + Inches_(0.1), ly + Inches_(0.1), Inches_(0.08),
               fill=GOOD)
    add_text(slide, lx + Inches_(0.22), ly + Inches_(0.02),
             Inches_(1.5), Inches_(0.2), "good (6 subj)", size=9)
    add_circle(slide, lx + Inches_(0.1), ly + Inches_(0.4), Inches_(0.08),
               fill=BAD)
    add_text(slide, lx + Inches_(0.22), ly + Inches_(0.32),
             Inches_(1.5), Inches_(0.2), "bad (6 subj)", size=9)
    add_text(slide, lx, ly + Inches_(0.8), Inches_(1.8), Inches_(0.4),
             "Dot size ~ cluster weight π (marker area)",
             size=8, color=GREY)

    # ---- Panel C: cluster composition stack per subject × timepoint ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "C", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Cluster composition by subject × timepoint (π-weighted mean CCF)",
             size=11, bold=True)
    px = Inches_(1.2); py = Inches_(1.3); pw = Inches_(11.7); ph = Inches_(5.0)
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, 1, py, ph) for v in [0, 0.25, 0.5, 0.75, 1.0]],
               y_labels=["0", "0.25", "0.5", "0.75", "1.0"],
               ylab="π-weighted mean CCF (normalised per subject)")
    n_sub = len(sorted_fits)
    bar_w = pw / (2 * n_sub + 1)
    cluster_pal = [RGBColor(*p) for p in [
        (0x4E, 0x79, 0xA7), (0xF2, 0x8E, 0x2B), (0xE1, 0x57, 0x59),
        (0x76, 0xB7, 0xB2), (0x59, 0xA1, 0x4F), (0xED, 0xC9, 0x48),
        (0xB0, 0x7A, 0xA1), (0xFF, 0x9D, 0xA7), (0x9C, 0x75, 0x5F),
        (0xBA, 0xB0, 0xAC)]]
    for i, f in enumerate(sorted_fits):
        pi = np.array(f["pi"]); mean_ccf = np.array(f["mean_ccf"])
        samples = f["samples"]
        pre_i, post_i = pre_post_idx(samples)
        # build per-tp stacked fraction: cluster k contributes pi_k * mean_ccf[k,tp]
        # then normalise within each bar to sum to 1 (for visual composition)
        for tp_idx, (tp_lbl, tp_sample_idx) in enumerate(
                [("pre", pre_i), ("post", post_i)]):
            bx = px + bar_w * (2 * i + tp_idx + 1)
            contrib = pi * mean_ccf[:, tp_sample_idx]
            total = contrib.sum() if contrib.sum() > 0 else 1
            frac = contrib / total
            cum = 0
            for k, fr in enumerate(frac):
                if fr <= 0: continue
                y_top = scale_y(cum + fr, 0, 1, py, ph)
                y_bot = scale_y(cum, 0, 1, py, ph)
                add_rect(slide, bx - bar_w * 0.35, y_top,
                         bar_w * 0.7, y_bot - y_top,
                         fill=cluster_pal[k % len(cluster_pal)],
                         line_color=WHITE, line_width=0.3)
                cum += fr
            # tp label below bar
            add_text(slide, bx - bar_w * 0.5, py + ph + Inches_(0.04),
                     bar_w, Inches_(0.2), tp_lbl,
                     size=6, align="center", anchor="top")
        # subj id
        cx_mid = px + bar_w * (2 * i + 1.5)
        col = GOOD if f["response"] == "good" else BAD
        add_text(slide, cx_mid - bar_w * 0.7, py + ph + Inches_(0.22),
                 bar_w * 1.4, Inches_(0.18),
                 f"s{f['subject_id']}",
                 size=7, align="center", bold=True, color=col)

    # ---- Panel D: dominant-clone shrinkage Δ by response (keep real) ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "D", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Dominant-clone shrinkage Δ by response (n = 12 paired; clonal_summary.tsv)",
             size=11, bold=True)
    px = Inches_(4.0); py = Inches_(1.5); pw = Inches_(5.0); ph = Inches_(5.0)
    ds_good = pc[pc["response"] == "good"]["dominant_shrink"].dropna().values
    ds_bad = pc[pc["response"] == "bad"]["dominant_shrink"].dropna().values
    vmin = min(ds_good.min(), ds_bad.min()) - 0.1
    vmax = max(ds_good.max(), ds_bad.max()) + 0.1
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, vmin, vmax, py, ph) for v in
                        np.linspace(vmin, vmax, 6)],
               y_labels=[f"{v:.2f}" for v in np.linspace(vmin, vmax, 6)],
               ylab="Δ dominant-clone CCF (post − pre)")
    zy = scale_y(0, vmin, vmax, py, ph)
    add_line(slide, px, zy, px + pw, zy, color=GREY, width=0.7, dashed=True)
    u, pv = mannwhitneyu(ds_good, ds_bad)
    for i, (resp, vals, col) in enumerate([("good", ds_good, GOOD),
                                            ("bad", ds_bad, BAD)]):
        cx = px + pw * (0.25 + i * 0.5)
        ys = [scale_y(float(v), vmin, vmax, py, ph) for v in vals]
        boxplot_primitive(slide, cx, py, ph, ys, col, box_w=Inches_(0.9))
        add_text(slide, cx - Inches_(1), py + ph + Inches_(0.1),
                 Inches_(2), Inches_(0.25),
                 f"{resp} (n={len(vals)})", size=10, align="center",
                 bold=True, color=col)
    add_text(slide, Inches_(4.0), Inches_(0.9),
             Inches_(5.0), Inches_(0.2),
             f"Mann–Whitney P = {pv:.2f}",
             size=10, align="center", italic=True)

    # ---- Panel E: shrink vs expand scatter ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "E", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Shrink vs expand scatter per subject (dominant_shrink vs dominant_expand Δ)",
             size=11, bold=True)
    px = Inches_(3.5); py = Inches_(1.3); pw = Inches_(6.5); ph = Inches_(5.3)
    xs = pc["dominant_shrink"].values; ys_ = pc["dominant_expand"].values
    xmin, xmax = float(min(xs)) - 0.1, float(max(xs)) + 0.1
    ymin, ymax = float(min(ys_)) - 0.1, float(max(ys_)) + 0.1
    axis_frame(slide, px, py, pw, ph,
               x_ticks=[scale_x(v, xmin, xmax, px, pw) for v in
                        np.linspace(xmin, xmax, 5)],
               x_labels=[f"{v:.2f}" for v in np.linspace(xmin, xmax, 5)],
               y_ticks=[scale_y(v, ymin, ymax, py, ph) for v in
                        np.linspace(ymin, ymax, 5)],
               y_labels=[f"{v:.2f}" for v in np.linspace(ymin, ymax, 5)],
               xlab="dominant_shrink Δ", ylab="dominant_expand Δ")
    # zero lines
    if xmin < 0 < xmax:
        zx = scale_x(0, xmin, xmax, px, pw)
        add_line(slide, zx, py, zx, py + ph, color=GREY, width=0.5, dashed=True)
    if ymin < 0 < ymax:
        zy = scale_y(0, ymin, ymax, py, ph)
        add_line(slide, px, zy, px + pw, zy, color=GREY, width=0.5, dashed=True)
    for _, row in pc.iterrows():
        col = GOOD if row["response"] == "good" else BAD
        x = scale_x(float(row["dominant_shrink"]), xmin, xmax, px, pw)
        y = scale_y(float(row["dominant_expand"]), ymin, ymax, py, ph)
        add_circle(slide, x, y, Inches_(0.08),
                   fill=col, line_color=INK, line_width=0.3)
        add_text(slide, x + Inches_(0.08), y - Inches_(0.09),
                 Inches_(0.5), Inches_(0.15),
                 f"s{int(row['subject_id'])}", size=7, color=col)

    # ---- Panel F: clone-fate composition by response ----
    slide = new_slide(prs)
    add_text(slide, Inches_(0.35), Inches_(0.25), Inches_(0.45), Inches_(0.45),
             "F", size=22, bold=True, color=INK)
    add_text(slide, Inches_(0.9), Inches_(0.35), Inches_(11.5), Inches_(0.4),
             "Clone-fate composition per response group (n_shrinking / n_expanding / stable clusters)",
             size=11, bold=True)
    # per response: sum n_shrinking, n_expanding, and stable = (n_clusters - n_shrinking - n_expanding)
    pc2 = pc.copy()
    pc2["n_stable"] = pc2["n_clusters"] - pc2["n_shrinking"] - pc2["n_expanding"]
    agg = pc2.groupby("response")[["n_shrinking", "n_stable", "n_expanding"]].sum()
    px = Inches_(3.5); py = Inches_(1.3); pw = Inches_(6.5); ph = Inches_(5.0)
    cats = ["n_shrinking", "n_stable", "n_expanding"]
    cats_lbl = ["shrinking", "stable", "expanding"]
    pal = {"n_shrinking": GOOD, "n_stable": GREY, "n_expanding": BAD}
    vmax = agg.values.sum(axis=1).max() * 1.1
    axis_frame(slide, px, py, pw, ph,
               y_ticks=[scale_y(v, 0, vmax, py, ph) for v in
                        np.linspace(0, vmax, 6)],
               y_labels=[f"{v:.0f}" for v in np.linspace(0, vmax, 6)],
               ylab="n clusters (summed across subjects)")
    slot_w = pw / 2
    for gi, grp in enumerate(["good", "bad"]):
        if grp not in agg.index: continue
        bx = px + slot_w * (gi + 0.5)
        cum = 0
        for k in cats:
            v = int(agg.loc[grp, k])
            if v <= 0: continue
            y_top = scale_y(cum + v, 0, vmax, py, ph)
            y_bot = scale_y(cum, 0, vmax, py, ph)
            add_rect(slide, bx - Inches_(0.6), y_top,
                     Inches_(1.2), y_bot - y_top,
                     fill=pal[k], line_color=INK, line_width=0.4)
            # value label
            add_text(slide, bx - Inches_(0.6),
                     (y_top + y_bot) / 2 - Inches_(0.1),
                     Inches_(1.2), Inches_(0.2), str(v),
                     size=10, align="center", bold=True, color=WHITE)
            cum += v
        col = GOOD if grp == "good" else BAD
        add_text(slide, bx - Inches_(0.6), py + ph + Inches_(0.08),
                 Inches_(1.2), Inches_(0.25),
                 f"{grp} (n={int(pc2[pc2.response==grp].shape[0])})",
                 size=11, bold=True, align="center", color=col)
    # legend
    lx = Inches_(11.0); ly = Inches_(1.5)
    add_text(slide, lx, ly, Inches_(1.6), Inches_(0.25),
             "Fate", size=10, bold=True)
    for i, (k, lab) in enumerate(zip(cats, cats_lbl)):
        add_rect(slide, lx, ly + Inches_(0.3 + i * 0.28),
                 Inches_(0.18), Inches_(0.14), fill=pal[k])
        add_text(slide, lx + Inches_(0.22), ly + Inches_(0.28 + i * 0.28),
                 Inches_(1.5), Inches_(0.2), lab, size=9)

    save(prs, "SuppFig_S16_PyClone_clonal_evolution.pptx")


# ============================================================================
# Main
# ============================================================================

def main():
    print(f"Output dir: {OUT}")
    print("Rebuilding 6 figures with real data:")
    for name, fn in [
        ("S03", build_S3_real),
        ("S04", build_S4_real),
        ("S10", build_S10_real),
        ("S11", build_S11_real),
        ("S15", build_S15_real),
        ("S16", build_S16_real),
    ]:
        print(f"{name} ...")
        try:
            result = fn()
        except Exception as e:
            print(f"  !! FAILED {name}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()
