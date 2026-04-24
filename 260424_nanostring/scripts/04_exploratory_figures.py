#!/usr/bin/env python
"""Figures for exploratory round A-F+H."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

ROOT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring")
TABLES = ROOT / "tables"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

GOOD = "#0a7d6e"
BAD = "#c53e1f"
GREY = "#666666"
GOLD = "#D4A300"
BLUE = "#1f5b82"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6})

COHORT = [(2,"good"),(4,"good"),(14,"good"),(10,"bad"),(11,"bad"),(13,"bad")]
SUBJ_OF_COL = {"TNT RNA 5":2, "TNT RNA 6":2, "TNT RNA 11":4, "TNT RNA 12":4,
               "TNT RNA 41":14, "TNT RNA 42":14, "TNT RNA 29":10, "TNT RNA 30":10,
               "TNT RNA 32":11, "TNT RNA 33":11, "TNT RNA 38":13, "TNT RNA 39":13}


# ============== Fig A: Pre/Post/Δ composite heatmap ==============
def fig_preposdelta_heatmap():
    pre = pd.read_csv(TABLES / "v2_pre_MW.tsv", sep="\t")
    post = pd.read_csv(TABLES / "v2_post_MW.tsv", sep="\t")
    dlt = pd.read_csv(TABLES / "v2_delta_MW.tsv", sep="\t")

    # Standardize ordering by pre P_1s good>bad
    order = pre.sort_values("MW_P_1s_good_gt_bad")["composite"].tolist()
    def lookup(df, composite_order, col):
        return df.set_index("composite").loc[composite_order, col].values

    pre_delta = lookup(pre, order, "good_mean") - lookup(pre, order, "bad_mean")
    post_delta = lookup(post, order, "good_mean") - lookup(post, order, "bad_mean")
    dlt_delta = lookup(dlt, order, "good_mean") - lookup(dlt, order, "bad_mean")
    pre_p = lookup(pre, order, "MW_P_1s_good_gt_bad")
    post_p = lookup(post, order, "MW_P_1s_good_gt_bad")
    dlt_p = lookup(dlt, order, "MW_P_1s_good_gt_bad")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6), gridspec_kw={"width_ratios": [3, 2]})

    # LEFT: delta matrix as heatmap (good_mean - bad_mean)
    M = np.vstack([pre_delta, post_delta, dlt_delta]).T  # 23 x 3
    vmax = float(np.abs(M).max())
    cmap = mcolors.LinearSegmentedColormap.from_list("gb", [BAD, "white", GOOD])
    im = ax1.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    ax1.set_xticks([0, 1, 2]); ax1.set_xticklabels(["pre", "post", "Δ"], fontsize=10, fontweight="bold")
    ax1.set_yticks(range(len(order))); ax1.set_yticklabels(order, fontsize=8.5)
    ax1.set_title("good − bad (composite z-score mean)", fontsize=9, fontweight="bold")
    # annotate values + stars for P ≤ 0.05
    P = np.vstack([pre_p, post_p, dlt_p]).T
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]; p = P[i, j]
            star = "★" if p <= 0.05 else ("·" if p <= 0.10 else "")
            col = "white" if abs(v) > vmax*0.6 else "black"
            ax1.text(j, i, f"{v:+.2f}\n{star}", ha="center", va="center",
                     fontsize=6.5, color=col)
    fig.colorbar(im, ax=ax1, shrink=0.5, label="Δ z (good − bad)")

    # RIGHT: P-value bars (one-sided good>bad, log scale)
    y = np.arange(len(order))[::-1]
    ax2.barh(y - 0.25, -np.log10(pre_p + 1e-12), height=0.2, color=BLUE, label="pre")
    ax2.barh(y + 0.00, -np.log10(post_p + 1e-12), height=0.2, color=GREY, label="post")
    ax2.barh(y + 0.25, -np.log10(dlt_p + 1e-12), height=0.2, color=GOLD, label="Δ")
    ax2.axvline(-np.log10(0.05), color="red", lw=0.6, ls="--", alpha=0.7)
    ax2.text(-np.log10(0.05), len(order), " P=0.05", fontsize=7, color="red", va="bottom")
    ax2.set_yticks(y); ax2.set_yticklabels([])
    ax2.set_xlabel("-log10 one-sided P (good > bad)", fontsize=8)
    ax2.legend(loc="lower right", fontsize=7, frameon=False)
    ax2.set_title("Directional evidence per timepoint", fontsize=9, fontweight="bold")

    fig.suptitle("NanoString composites: pre vs post vs Δ (n=3 good pCR vs n=3 bad poor)",
                 fontsize=10, fontweight="bold", y=0.99)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_pre_post_delta_heatmap.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============== Fig B: Canonical signatures pre/post boxplot ==============
def fig_canonical():
    pre_df = pd.read_csv(TABLES / "v2_composite_pre.tsv", sep="\t", index_col="subject")
    post_df = pd.read_csv(TABLES / "v2_composite_post.tsv", sep="\t", index_col="subject")
    pre_mw = pd.read_csv(TABLES / "v2_pre_MW.tsv", sep="\t").set_index("composite")
    post_mw = pd.read_csv(TABLES / "v2_post_MW.tsv", sep="\t").set_index("composite")

    targets = ["Ayers_TIS", "IFNg_6", "IFNg_10_Ayers", "CD8_cytotoxic", "IMPRES_pos", "M1_macro"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))

    for ax, sig in zip(axes.ravel(), targets):
        for df_, col_offset, tag in [(pre_df, 0, "pre"), (post_df, 1, "post")]:
            for s, bn in COHORT:
                x = col_offset + (0 if bn == "good" else 0.35)
                y = df_.loc[s, sig]
                col = GOOD if bn == "good" else BAD
                ax.scatter(x + np.random.uniform(-0.06, 0.06), y, s=55, color=col,
                           edgecolor="black", lw=0.5, zorder=3)
            # group means
            good_m = df_.loc[[s for s,b in COHORT if b=="good"], sig].mean()
            bad_m = df_.loc[[s for s,b in COHORT if b=="bad"], sig].mean()
            ax.plot([col_offset-0.12, col_offset+0.12], [good_m]*2, color=GOOD, lw=2)
            ax.plot([col_offset+0.22, col_offset+0.48], [bad_m]*2, color=BAD, lw=2)
        # labels
        ax.set_xticks([0, 0.35, 1, 1.35])
        ax.set_xticklabels(["good", "bad", "good", "bad"], fontsize=8)
        # top bracket labels
        ax.set_title(sig.replace("_", " "), fontsize=9, fontweight="bold")
        ax.set_ylabel("composite z", fontsize=8)
        # annotate P
        p_pre = pre_mw.loc[sig, "MW_P_1s_good_gt_bad"]
        p_post = post_mw.loc[sig, "MW_P_1s_good_gt_bad"]
        ax.text(0.175, ax.get_ylim()[1]*0.95, f"P={p_pre:.2f}",
                ha="center", fontsize=8, color=(GOOD if p_pre<=0.10 else GREY),
                fontweight="bold" if p_pre<=0.05 else "normal")
        ax.text(1.175, ax.get_ylim()[1]*0.95, f"P={p_post:.2f}",
                ha="center", fontsize=8, color=(GOOD if p_post<=0.10 else GREY),
                fontweight="bold" if p_post<=0.05 else "normal")
        # dividers
        ax.axvline(0.675, color=GREY, lw=0.5, ls=":")
        # subtitle for pre/post groups
        ax.text(0.175, ax.get_ylim()[0] - 0.12*(ax.get_ylim()[1]-ax.get_ylim()[0]),
                "pre", ha="center", fontsize=8, color=BLUE, fontweight="bold")
        ax.text(1.175, ax.get_ylim()[0] - 0.12*(ax.get_ylim()[1]-ax.get_ylim()[0]),
                "post", ha="center", fontsize=8, color=BLUE, fontweight="bold")

    leg = [Line2D([0],[0], marker="o", color="w", markerfacecolor=GOOD, markeredgecolor="black",
                  label="good pCR (n=3)", markersize=7),
           Line2D([0],[0], marker="o", color="w", markerfacecolor=BAD, markeredgecolor="black",
                  label="bad poor (n=3)", markersize=7)]
    fig.legend(handles=leg, loc="upper center", ncol=2, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, 1.03))
    fig.suptitle("Regulatory-grade signatures: pre- and post-SC-RT comparison (NanoString)",
                 fontsize=11, fontweight="bold", y=1.07)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_canonical_signatures.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============== Fig C: IAE vs IBI fingerprint heatmap ==============
def fig_iae_ibi():
    iaa = pd.read_csv(TABLES / "v2_IAE_vs_IBI_descriptive.tsv", sep="\t")
    genes = pd.read_csv(TABLES / "v2_IAE_vs_IBI_gene_descriptive.tsv", sep="\t")

    top_comp = iaa.head(16)
    top_genes = genes.head(24)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6.5),
                                    gridspec_kw={"width_ratios": [1, 1]})

    # LEFT: composite fingerprint
    y = np.arange(len(top_comp))[::-1]
    vals = top_comp["IAE_minus_IBI"].values
    colors = [GOOD if v > 0 else BAD for v in vals]
    ax1.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    ax1.set_yticks(y); ax1.set_yticklabels(top_comp["feature"], fontsize=8)
    ax1.axvline(0, color="black", lw=0.5)
    ax1.set_xlabel("IAE − IBI mean Δ (composite z)", fontsize=8)
    ax1.set_title("Composite level: IAE (n=2) vs IBI (n=3)", fontsize=9, fontweight="bold")
    # annotate values
    for yi, v in zip(y, vals):
        ax1.text(v + (0.05 if v > 0 else -0.05), yi, f"{v:+.2f}",
                 ha="left" if v > 0 else "right", va="center", fontsize=7, color=GREY)

    # RIGHT: top genes
    y = np.arange(len(top_genes))[::-1]
    vals = top_genes["IAE_minus_IBI"].values
    colors = [GOOD if v > 0 else BAD for v in vals]
    ax2.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    ax2.set_yticks(y); ax2.set_yticklabels(top_genes["gene"], fontsize=7.5, family="monospace")
    ax2.axvline(0, color="black", lw=0.5)
    ax2.set_xlabel("IAE − IBI mean Δ (log2 z)", fontsize=8)
    ax2.set_title("Gene level: top 24 by |IAE − IBI|", fontsize=9, fontweight="bold")
    for yi, v in zip(y, vals):
        ax2.text(v + (0.1 if v > 0 else -0.1), yi, f"{v:+.1f}",
                 ha="left" if v > 0 else "right", va="center", fontsize=6.5, color=GREY)

    leg = [Line2D([0],[0], color=GOOD, lw=3, label="IAE > IBI (good inflamed has more)"),
           Line2D([0],[0], color=BAD, lw=3, label="IBI > IAE (bad inflamed has more)")]
    fig.legend(handles=leg, loc="upper center", ncol=2, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Inflamed-but-Ineffective (IBI, n=3 bad inflamed) vs "
                 "Inflamed-Active-Effective (IAE, n=2 good inflamed)",
                 fontsize=10, fontweight="bold", y=1.06)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_IAE_vs_IBI.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============== Fig D: Subject radar (subj 4 + 11) ==============
def fig_subject_radar():
    pre_df = pd.read_csv(TABLES / "v2_composite_pre.tsv", sep="\t", index_col="subject")
    post_df = pd.read_csv(TABLES / "v2_composite_post.tsv", sep="\t", index_col="subject")
    dlt_df = pd.read_csv(TABLES / "v2_composite_delta.tsv", sep="\t", index_col="subject")

    axes_cats = ["TLS_8","Plasma_proxy","GC_TF","Naive_B","Memory_B","Ayers_TIS",
                 "IFNg_6","CD8_cytotoxic","Teff_cytotoxic","Treg","CD8_exh",
                 "HLA_II","HLA_I_machinery_narrow","M1_macro","M2_macro",
                 "NK_activating","DC_mature"]
    n = len(axes_cats)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False).tolist()
    angles_closed = angles + angles[:1]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5),
                              subplot_kw=dict(projection="polar"))
    rmax = 2.8

    for ax, subj in zip(axes, [4, 2, 11]):
        bn = [b for s, b in COHORT if s == subj][0]
        col = GOOD if bn == "good" else BAD
        title = f"Subject {subj} ({bn}{' pCR' if bn=='good' else ' poor'})"
        # plot pre, post, delta
        pre_vals = pre_df.loc[subj, axes_cats].values.tolist() + [pre_df.loc[subj, axes_cats[0]]]
        post_vals = post_df.loc[subj, axes_cats].values.tolist() + [post_df.loc[subj, axes_cats[0]]]
        dlt_vals = dlt_df.loc[subj, axes_cats].values.tolist() + [dlt_df.loc[subj, axes_cats[0]]]

        ax.plot(angles_closed, pre_vals, color=BLUE, lw=1.2, label="pre", alpha=0.8)
        ax.fill(angles_closed, pre_vals, color=BLUE, alpha=0.1)
        ax.plot(angles_closed, post_vals, color=col, lw=1.8, label="post")
        ax.fill(angles_closed, post_vals, color=col, alpha=0.15)

        ax.set_xticks(angles)
        ax.set_xticklabels([c.replace("_", "\n").replace("-", "\n") for c in axes_cats],
                           fontsize=6.5)
        ax.set_ylim(-rmax, rmax)
        ax.set_yticks([-2, -1, 0, 1, 2])
        ax.set_yticklabels([], fontsize=6)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=18)
        ax.legend(loc="lower right", fontsize=7, frameon=False, bbox_to_anchor=(1.15, -0.12))

    fig.suptitle("Subject fingerprints: s4 (atypical good) vs s2 (textbook good) vs s11 (bad, RNA-seq gap)",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_subject_radar.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============== Fig E: Directional concordance waterfall ==============
def fig_direction_waterfall():
    pre = pd.read_csv(TABLES / "v2_pre_MW.tsv", sep="\t")
    pre = pre.sort_values("MW_P_1s_good_gt_bad")
    vals = (pre["good_mean"] - pre["bad_mean"]).values

    fig, ax = plt.subplots(figsize=(9, 4.5))
    y = np.arange(len(pre))[::-1]
    colors = [GOOD if v > 0 else BAD for v in vals]
    ax.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(pre["composite"], fontsize=8)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Pre-treatment good − bad (composite z mean)", fontsize=9)
    ax.set_title(f"23/23 composites show good > bad direction pre-treatment "
                 f"(sign test vs 50:50 ⇒ P = (1/2)^23 ≈ 1.2e-7 if independent; "
                 f"composites share genes, so this is an upper bound)",
                 fontsize=8.5, fontweight="bold")
    # annotate P
    for yi, (_, row) in zip(y, pre.iterrows()):
        delta = row["good_mean"] - row["bad_mean"]
        p = row["MW_P_1s_good_gt_bad"]
        star = "★" if p <= 0.05 else ("·" if p <= 0.10 else "")
        ax.text(delta + (0.02 if delta > 0 else -0.02), yi,
                f" P={p:.2f} {star}",
                ha="left" if delta > 0 else "right", va="center",
                fontsize=6.5, color=GREY,
                fontweight="bold" if p <= 0.05 else "normal")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_direction_waterfall_pre.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============== Fig F: HLA-I machinery per-subject heatmap ==============
def fig_hla_i_heatmap():
    zmat = pd.read_csv(TABLES / "logz_matrix.tsv", sep="\t", index_col=0)
    genes = ["NLRC5","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9",
             "TAPBP","HLA-E","HLA-F","HLA-G","CIITA"]
    genes = [g for g in genes if g in zmat.index]
    subj_cols = [(s, bn, pre, post) for (s, bn), (_, pre, post, _) in
                 zip(COHORT, [(c[0], c[1], c[2], c[3]) for c in
                              [(2,"TNT RNA 5","TNT RNA 6","good"),
                               (4,"TNT RNA 11","TNT RNA 12","good"),
                               (14,"TNT RNA 41","TNT RNA 42","good"),
                               (10,"TNT RNA 29","TNT RNA 30","bad"),
                               (11,"TNT RNA 32","TNT RNA 33","bad"),
                               (13,"TNT RNA 38","TNT RNA 39","bad")]])]
    cols_flat = []
    for s, bn, pre, post in subj_cols:
        cols_flat.append((f"s{s}-pre", pre, bn))
        cols_flat.append((f"s{s}-post", post, bn))
    mat = np.zeros((len(genes), len(cols_flat)))
    for i, g in enumerate(genes):
        for j, (_, col, _) in enumerate(cols_flat):
            mat[i, j] = zmat.loc[g, col]

    fig, ax = plt.subplots(figsize=(10, 5.5))
    vmax = float(np.abs(mat).max())
    cmap = mcolors.LinearSegmentedColormap.from_list("pw", ["#3a567e","white","#d15b28"])
    im = ax.imshow(mat, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols_flat)))
    ax.set_xticklabels([c[0] for c in cols_flat], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8, family="monospace")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            col = "white" if abs(v) > vmax*0.6 else "black"
            ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=5.5, color=col)
    # response bar on top
    for j, (_, _, bn) in enumerate(cols_flat):
        ax.add_patch(plt.Rectangle((j-0.5, -1.2), 1, 0.5,
                                     facecolor=GOOD if bn=="good" else BAD,
                                     edgecolor="none", clip_on=False))
    ax.set_ylim(len(genes)-0.5, -1.8)
    ax.set_title("HLA class I antigen presentation machinery per subject × timepoint (log2+z)",
                 fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.6, label="log2+z")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"FigEx_HLA_I_machinery_heatmap.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_preposdelta_heatmap()
    fig_canonical()
    fig_iae_ibi()
    fig_subject_radar()
    fig_direction_waterfall()
    fig_hla_i_heatmap()
    print("=== Figures saved ===")
    for p in sorted(FIG.glob("FigEx_*.pdf")):
        print(f"  {p.name}")
