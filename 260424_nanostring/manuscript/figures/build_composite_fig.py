#!/usr/bin/env python
"""Assemble manuscript-ready composite Supp Fig S22 from exploratory NanoString figures."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

SRC_TAB = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring/tables")
OUT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring/manuscript/figures")
OUT.mkdir(parents=True, exist_ok=True)

GOOD = "#0a7d6e"
BAD = "#c53e1f"
GREY = "#666666"
GOLD = "#D4A300"
BLUE = "#1f5b82"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6})

COHORT = [(2,"good"),(4,"good"),(14,"good"),(10,"bad"),(11,"bad"),(13,"bad")]


def build():
    pre_df = pd.read_csv(SRC_TAB / "v2_composite_pre.tsv", sep="\t", index_col="subject")
    post_df = pd.read_csv(SRC_TAB / "v2_composite_post.tsv", sep="\t", index_col="subject")
    dlt_df = pd.read_csv(SRC_TAB / "v2_composite_delta.tsv", sep="\t", index_col="subject")
    pre_mw = pd.read_csv(SRC_TAB / "v2_pre_MW.tsv", sep="\t").set_index("composite")
    post_mw = pd.read_csv(SRC_TAB / "v2_post_MW.tsv", sep="\t").set_index("composite")
    dlt_mw = pd.read_csv(SRC_TAB / "v2_delta_MW.tsv", sep="\t").set_index("composite")
    iae = pd.read_csv(SRC_TAB / "v2_IAE_vs_IBI_descriptive.tsv", sep="\t")

    fig = plt.figure(figsize=(13.5, 11.5))
    gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.35,
                           height_ratios=[1.2, 1.0, 1.0])

    # ============== Panel A: pre/post/Δ heatmap ==============
    axA = fig.add_subplot(gs[0, :2])
    order = pre_mw.sort_values("MW_P_1s_good_gt_bad").index.tolist()
    pre_delta = (pre_mw.loc[order, "good_mean"] - pre_mw.loc[order, "bad_mean"]).values
    post_delta = (post_mw.loc[order, "good_mean"] - post_mw.loc[order, "bad_mean"]).values
    dlt_delta = (dlt_mw.loc[order, "good_mean"] - dlt_mw.loc[order, "bad_mean"]).values
    M = np.vstack([pre_delta, post_delta, dlt_delta]).T
    P = np.vstack([pre_mw.loc[order, "MW_P_1s_good_gt_bad"].values,
                    post_mw.loc[order, "MW_P_1s_good_gt_bad"].values,
                    dlt_mw.loc[order, "MW_P_1s_good_gt_bad"].values]).T
    vmax = float(np.abs(M).max())
    cmap = mcolors.LinearSegmentedColormap.from_list("gb", [BAD, "white", GOOD])
    im = axA.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    axA.set_xticks([0, 1, 2]); axA.set_xticklabels(["pre", "post", "Δ"], fontsize=9, fontweight="bold")
    axA.set_yticks(range(len(order))); axA.set_yticklabels(order, fontsize=7.5)
    axA.set_title("A   Pre/post/Δ — good minus bad composite z (★ one-sided P ≤ 0.05, · ≤ 0.10)",
                   fontsize=9.5, fontweight="bold", loc="left")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]; p = P[i, j]
            star = "★" if p <= 0.05 else ("·" if p <= 0.10 else "")
            col = "white" if abs(v) > vmax*0.6 else "black"
            axA.text(j, i, f"{v:+.2f}\n{star}", ha="center", va="center", fontsize=5.5, color=col)
    fig.colorbar(im, ax=axA, shrink=0.55, label="good − bad (z)")

    # ============== Panel B: Pre-treatment waterfall ==============
    axB = fig.add_subplot(gs[0, 2])
    vals = pre_delta
    colors = [GOOD if v > 0 else BAD for v in vals]
    y = np.arange(len(order))[::-1]
    axB.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    axB.set_yticks(y); axB.set_yticklabels([o[:18] for o in order], fontsize=6.5)
    axB.axvline(0, color="black", lw=0.4)
    axB.set_xlabel("Pre good − bad (z)", fontsize=8)
    axB.set_title("B   23/23 good > bad pre-treatment",
                   fontsize=9.5, fontweight="bold", loc="left")
    for yi, v, p in zip(y, vals, pre_mw.loc[order, "MW_P_1s_good_gt_bad"].values):
        mark = "★" if p <= 0.05 else ""
        axB.text(v + (0.02 if v > 0 else -0.02), yi, f" {mark}",
                  ha="left" if v > 0 else "right", va="center",
                  fontsize=8, color=GOLD, fontweight="bold")

    # ============== Panel C: Canonical signatures 2x3 ==============
    targets = [("Ayers_TIS", "Ayers TIS"),
               ("IFNg_6", "IFN-γ 6-gene"),
               ("IFNg_10_Ayers", "IFN-γ 10-gene"),
               ("CD8_cytotoxic", "CD8 cytotoxic"),
               ("IMPRES_pos", "IMPRES"),
               ("M1_macro", "M1 macrophage")]
    for i, (key, title) in enumerate(targets):
        axC = fig.add_subplot(gs[1, 0]) if i == 0 else None
        if i == 0:
            # Create sub-gridspec for 2x3 canonical panel
            from matplotlib.gridspec import GridSpecFromSubplotSpec
            sub_gs = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[1, :], hspace=0.55, wspace=0.35)
            axC.remove()
            canonical_axes = [fig.add_subplot(sub_gs[r, c]) for r in range(2) for c in range(3)]

        ax = canonical_axes[i]
        for df_, x_off, tag in [(pre_df, 0, "pre"), (post_df, 1.0, "post")]:
            for subj, bn in COHORT:
                xpos = x_off + (0 if bn == "good" else 0.35)
                xpos += np.random.uniform(-0.05, 0.05)
                y = df_.loc[subj, key]
                col = GOOD if bn == "good" else BAD
                ax.scatter(xpos, y, s=42, color=col, edgecolor="black", lw=0.4, zorder=3)
            gm = df_.loc[[s for s,b in COHORT if b=="good"], key].mean()
            bm = df_.loc[[s for s,b in COHORT if b=="bad"], key].mean()
            ax.plot([x_off-0.12, x_off+0.12], [gm]*2, color=GOOD, lw=2)
            ax.plot([x_off+0.22, x_off+0.48], [bm]*2, color=BAD, lw=2)
        ax.set_xticks([0, 0.35, 1.0, 1.35])
        ax.set_xticklabels(["g","b","g","b"], fontsize=7)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_ylabel("z", fontsize=7)
        p_pre = pre_mw.loc[key, "MW_P_1s_good_gt_bad"]
        p_post = post_mw.loc[key, "MW_P_1s_good_gt_bad"]
        ax.axvline(0.68, color=GREY, lw=0.4, ls=":")
        ax.text(0.175, ax.get_ylim()[1]*0.92, f"pre P={p_pre:.2f}",
                ha="center", fontsize=7,
                color=(GOLD if p_pre<=0.05 else GREY),
                fontweight="bold" if p_pre<=0.05 else "normal")
        ax.text(1.175, ax.get_ylim()[1]*0.92, f"post P={p_post:.2f}",
                ha="center", fontsize=7,
                color=(GOLD if p_post<=0.05 else GREY),
                fontweight="bold" if p_post<=0.05 else "normal")
        if i == 0:
            ax.text(-0.22, 1.12, "C   Regulatory-grade signatures (left: pre, right: post)",
                    transform=ax.transAxes, fontsize=9.5, fontweight="bold")

    # ============== Panel D: IBI vs IAE composite ==============
    axD = fig.add_subplot(gs[2, 0])
    top = iae.head(12)
    y = np.arange(len(top))[::-1]
    vals = top["IAE_minus_IBI"].values
    colors = [GOOD if v > 0 else BAD for v in vals]
    axD.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    axD.set_yticks(y); axD.set_yticklabels(top["feature"], fontsize=7)
    axD.axvline(0, color="black", lw=0.4)
    axD.set_xlabel("IAE − IBI mean Δ (z)", fontsize=8)
    axD.set_title("D   IAE (n=2) vs IBI (n=3)  composite fingerprint",
                   fontsize=9.5, fontweight="bold", loc="left")

    # ============== Panel E: Top genes IAE vs IBI ==============
    axE = fig.add_subplot(gs[2, 1])
    genes = pd.read_csv(SRC_TAB / "v2_IAE_vs_IBI_gene_descriptive.tsv", sep="\t").head(16)
    y = np.arange(len(genes))[::-1]
    vals = genes["IAE_minus_IBI"].values
    colors = [GOOD if v > 0 else BAD for v in vals]
    axE.barh(y, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.7)
    axE.set_yticks(y); axE.set_yticklabels(genes["gene"], fontsize=7, family="monospace")
    axE.axvline(0, color="black", lw=0.4)
    axE.set_xlabel("IAE − IBI (gene log2 z)", fontsize=8)
    axE.set_title("E   Top discriminator genes",
                   fontsize=9.5, fontweight="bold", loc="left")

    # ============== Panel F: Cascade T2 scatter ==============
    axF = fig.add_subplot(gs[2, 2])
    cd = pd.read_csv(SRC_TAB / "composite_subject_delta.tsv", sep="\t", index_col=0)
    for subj, bn in COHORT:
        x = cd.loc[f"subj_{subj}", "CD8_exh"]
        y = cd.loc[f"subj_{subj}", "TLS_8"]
        col = GOOD if bn == "good" else BAD
        axF.scatter(x, y, s=55, color=col, edgecolor="black", lw=0.5, zorder=3)
        axF.annotate(f"s{subj}", (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7)
    # regression
    xs = cd.loc[[f"subj_{s}" for s,_ in COHORT], "CD8_exh"].values
    ys = cd.loc[[f"subj_{s}" for s,_ in COHORT], "TLS_8"].values
    m, c = np.polyfit(xs, ys, 1)
    xr = np.linspace(xs.min(), xs.max(), 30)
    axF.plot(xr, m*xr + c, color=GREY, lw=1, ls="--", alpha=0.6)
    axF.set_xlabel("Δ CD8 exhaustion (z)", fontsize=8)
    axF.set_ylabel("Δ TLS-8 (z)", fontsize=8)
    axF.set_title("F   Pre-spec T2 cascade arrow 4\n(Pearson r = +0.82, P = 0.046)",
                   fontsize=9.5, fontweight="bold", loc="left")
    axF.axhline(0, color=GREY, lw=0.3, ls=":"); axF.axvline(0, color=GREY, lw=0.3, ls=":")

    # Legend
    leg = [Line2D([0],[0], marker="o", color="w", markerfacecolor=GOOD, markeredgecolor="black",
                  label="good pCR (n=3)", markersize=7),
           Line2D([0],[0], marker="o", color="w", markerfacecolor=BAD, markeredgecolor="black",
                  label="bad poor (n=3)", markersize=7)]
    fig.legend(handles=leg, loc="upper center", ncol=2, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, 1.00))
    fig.suptitle("Supp Fig S22.  NanoString PanCancer Immune orthogonal validation — "
                 "pre-registered Arrow 5 Δ null + extreme-phenotype pre-treatment immune axis replication",
                 fontsize=10.5, fontweight="bold", y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"FigS22_NanoString_composite.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT}/FigS22_NanoString_composite.pdf")


if __name__ == "__main__":
    build()
