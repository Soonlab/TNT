#!/usr/bin/env python
"""Figure set for NanoString Arrow 5 rescue test."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring")
TABLES = ROOT / "tables"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

GOOD = "#0a7d6e"
BAD = "#c53e1f"
GREY = "#666666"
GOLD = "#D4A300"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6})

COHORT = [(2,"good"),(4,"good"),(14,"good"),(10,"bad"),(11,"bad"),(13,"bad")]


def fig_primary_paired():
    """Fig A: pre/post slope per subject for 4 primary composites."""
    comp = pd.read_csv(TABLES / "composite_scores.tsv", sep="\t")
    pri = ["TLS_8", "Plasma_proxy", "GC_TF"]
    all_comps = ["CXCL13"] + pri  # 4 targets

    # CXCL13 values from logz matrix
    zmat = pd.read_csv(TABLES / "logz_matrix.tsv", sep="\t", index_col=0)

    fig, axes = plt.subplots(1, 4, figsize=(11, 2.9), sharey=False)
    pairs = [(2,"TNT RNA 5","TNT RNA 6","good"),(4,"TNT RNA 11","TNT RNA 12","good"),
             (14,"TNT RNA 41","TNT RNA 42","good"),(10,"TNT RNA 29","TNT RNA 30","bad"),
             (11,"TNT RNA 32","TNT RNA 33","bad"),(13,"TNT RNA 38","TNT RNA 39","bad")]

    for ax, target in zip(axes, all_comps):
        if target == "CXCL13":
            pre_vals = {s: zmat.loc["CXCL13", p] for (s, p, q, _b) in pairs}
            post_vals = {s: zmat.loc["CXCL13", q] for (s, p, q, _b) in pairs}
            title = "P1  CXCL13"
        else:
            sub = comp[comp["composite"] == target]
            pre_vals = {s: float(sub[sub["sample"] == p]["score"].iloc[0]) for (s, p, q, _b) in pairs}
            post_vals = {s: float(sub[sub["sample"] == q]["score"].iloc[0]) for (s, p, q, _b) in pairs}
            title = {"TLS_8":"P2  TLS-8", "Plasma_proxy":"P3  Plasma-proxy", "GC_TF":"P4  GC-TF"}[target]

        for subj, _p, _q, bn in pairs:
            col = GOOD if bn == "good" else BAD
            ax.plot([0, 1], [pre_vals[subj], post_vals[subj]], "-o",
                    color=col, lw=1.2, ms=5, alpha=0.85)
            ax.annotate(f"s{subj}", (1, post_vals[subj]), xytext=(4, 0),
                        textcoords="offset points", fontsize=7, color=col, va="center")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["pre", "post"])
        ax.set_xlim(-0.25, 1.5)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.axhline(0, color=GREY, lw=0.4, ls="--", alpha=0.6)
        ax.set_ylabel("z-score" if target == "CXCL13" else "composite z", fontsize=8)

        # Δ summary
        g_delta = np.mean([post_vals[s] - pre_vals[s] for (s, *_) in pairs if _[-1] == "good"])
        b_delta = np.mean([post_vals[s] - pre_vals[s] for (s, *_) in pairs if _[-1] == "bad"])
        ax.text(0.02, 0.97, f"Δ good={g_delta:+.2f}\nΔ bad ={b_delta:+.2f}",
                transform=ax.transAxes, fontsize=7, va="top", ha="left",
                family="monospace", color=GREY)

    # legend
    from matplotlib.lines import Line2D
    leg = [Line2D([0],[0], marker="o", color=GOOD, label="good (pCR)", lw=1.2, ms=5),
           Line2D([0],[0], marker="o", color=BAD, label="bad (poor)", lw=1.2, ms=5)]
    fig.legend(handles=leg, loc="upper center", ncol=2, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Primary Arrow 5 rescue tests (P1–P4) — NanoString nCounter paired Δ",
                 y=1.10, fontsize=10, fontweight="bold")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"Fig_primary_paired.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def fig_mw_summary():
    """Fig B: horizontal bar of |Δ good - Δ bad| with direction arrow, P1-P4 + S1-S3."""
    p = pd.read_csv(TABLES / "P1_P4_primary.tsv", sep="\t")
    s = pd.read_csv(TABLES / "S1_S3_lineage.tsv", sep="\t")
    df = pd.concat([p, s], ignore_index=True)
    df["delta_diff"] = df["good_mean"] - df["bad_mean"]
    df["direction"] = np.where(df["delta_diff"] > 0, "good>bad", "bad>good")
    df["P_1s"] = df["MW_P_1s_good_gt_bad"]
    df["label"] = df["target"].str.replace("_", " ")

    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    y = np.arange(len(df))[::-1]
    colors = [GOOD if d > 0 else BAD for d in df["delta_diff"]]
    ax.barh(y, df["delta_diff"], color=colors, edgecolor="black", linewidth=0.4, height=0.62)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels(df["label"], fontsize=9)
    ax.set_xlabel("Δ good − Δ bad (composite z-score units)", fontsize=9)
    ax.set_title("Between-group Δ difference (positive = good > bad, pre-registered direction)",
                 fontsize=9, fontweight="bold")
    # P annotations
    for i, (yi, row) in enumerate(zip(y, df.itertuples())):
        xt = row.delta_diff + (0.05 if row.delta_diff > 0 else -0.05)
        ha = "left" if row.delta_diff > 0 else "right"
        ax.text(xt, yi, f"  1s P={row.P_1s:.2f}", va="center", ha=ha, fontsize=7.5,
                color=GREY, family="monospace")
    # pre-spec direction shading
    ax.axvspan(0, ax.get_xlim()[1]*1.02, alpha=0.05, color=GOOD, zorder=-5)
    ax.axvspan(ax.get_xlim()[0]*1.02, 0, alpha=0.05, color=BAD, zorder=-5)
    ax.text(0.99, 0.02, "pre-registered alternative: good > bad",
            transform=ax.transAxes, fontsize=7, color=GREY, ha="right", va="bottom", style="italic")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"Fig_primary_secondary_bar.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def fig_cascade_scatter():
    """Fig C: T1–T3 scatter (Δ composite A × Δ composite B)."""
    t = pd.read_csv(TABLES / "T1_T3_cascade.tsv", sep="\t")
    cd = pd.read_csv(TABLES / "composite_subject_delta.tsv", sep="\t", index_col=0)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    defs = [
        ("T1  HLA-II × Plasma-proxy",  "HLA_II",     "Plasma_proxy"),
        ("T2  CD8-exh × TLS-8",        "CD8_exh",    "TLS_8"),
        ("T3  HLA-I × T-cell",         "HLA_I_axis", "T_cell"),
    ]
    for ax, (title, a, b) in zip(axes, defs):
        for subj, bn in COHORT:
            col = GOOD if bn == "good" else BAD
            ax.scatter(cd.loc[f"subj_{subj}", a], cd.loc[f"subj_{subj}", b],
                       color=col, s=60, edgecolor="black", lw=0.6, zorder=3)
            ax.annotate(f"s{subj}", (cd.loc[f'subj_{subj}', a], cd.loc[f'subj_{subj}', b]),
                        xytext=(5, 3), textcoords="offset points", fontsize=7)
        # regression line
        x_ = cd[a].values; y_ = cd[b].values
        m, b_ = np.polyfit(x_, y_, 1)
        xs = np.linspace(x_.min(), x_.max(), 50)
        ax.plot(xs, m*xs + b_, color=GREY, lw=1, ls="--", alpha=0.7, zorder=1)
        row = t[t["test"].str.startswith(title.split()[0])].iloc[0]
        ax.text(0.02, 0.97, f"r={row.pearson_r:+.2f}  P={row.pearson_P:.3f}\n"
                            f"ρ={row.spearman_r:+.2f}  P={row.spearman_P:.3f}",
                transform=ax.transAxes, fontsize=7, va="top",
                family="monospace", color=GREY,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GREY, lw=0.4))
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_xlabel(f"Δ {a}", fontsize=8)
        ax.set_ylabel(f"Δ {b}", fontsize=8)
        ax.axhline(0, color=GREY, lw=0.3, ls=":"); ax.axvline(0, color=GREY, lw=0.3, ls=":")
    fig.suptitle("Cascade internal coherence (T1–T3, n=6 paired)",
                 fontsize=10, fontweight="bold", y=1.03)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"Fig_cascade_scatter.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def fig_platform_concordance():
    """Fig D: NanoString Δ × RNA-seq Δ per-gene Pearson r histogram + top genes."""
    s4 = pd.read_csv(TABLES / "S4_platform_concordance.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))
    ax = axes[0]
    ax.hist(s4["pearson_r"], bins=40, color=GOOD, edgecolor="black", linewidth=0.3, alpha=0.75)
    med = s4["pearson_r"].median()
    ax.axvline(med, color=GOLD, lw=1.6, label=f"median r={med:.2f}")
    ax.axvline(0, color="black", lw=0.4, ls=":")
    pos = (s4["pearson_r"] > 0).mean() * 100
    ax.set_xlabel("Pearson r (NanoString Δ × RNA-seq Δ, n=5 subjects)", fontsize=8)
    ax.set_ylabel("Genes", fontsize=8)
    ax.set_title(f"Platform concordance: {pos:.1f}% genes r>0", fontsize=9, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    ax = axes[1]
    top = s4.sort_values("pearson_r", ascending=False).head(15)
    ax.barh(np.arange(len(top))[::-1], top["pearson_r"], color=GOOD, edgecolor="black", linewidth=0.3)
    ax.set_yticks(np.arange(len(top))[::-1]); ax.set_yticklabels(top["gene"], fontsize=7.5)
    ax.set_xlabel("Pearson r", fontsize=8)
    ax.set_title("Top 15 concordant genes", fontsize=9, fontweight="bold")
    ax.set_xlim(0.9, 1.0)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"Fig_platform_concordance.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def fig_subject_fingerprint():
    """Fig E: heatmap of 10 composites × 6 subjects (Δ z)."""
    cd = pd.read_csv(TABLES / "composite_subject_delta.tsv", sep="\t", index_col=0)
    order = ["TLS_8", "Plasma_proxy", "GC_TF", "Naive_B", "Memory_B", "BAFF_APRIL",
             "HLA_II", "CD8_exh", "HLA_I_axis", "T_cell"]
    subj_order = [f"subj_{s}" for s, _ in COHORT]
    mat = cd.loc[subj_order, order].T  # composites × subjects

    fig, ax = plt.subplots(figsize=(6, 4))
    import matplotlib.colors as mcolors
    vmax = float(np.max(np.abs(mat.values)))
    cmap = mcolors.LinearSegmentedColormap.from_list("ig", [BAD, "white", GOOD])
    im = ax.imshow(mat.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(subj_order)))
    ax.set_xticklabels([f"s{s}\n({bn})" for s, bn in COHORT], fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    # annotate values
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            col = "white" if abs(v) > vmax*0.6 else "black"
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.5, color=col)
    ax.set_title("Subject-level composite Δ fingerprint (post − pre)", fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Δ z-score")
    # group annotation bar
    for j, (s, bn) in enumerate(COHORT):
        ax.add_patch(plt.Rectangle((j-0.5, -1.1), 1, 0.5, facecolor=GOOD if bn=="good" else BAD,
                                   edgecolor="none", clip_on=False))
    ax.set_ylim(len(order)-0.5, -1.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"Fig_subject_fingerprint.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_primary_paired()
    fig_mw_summary()
    fig_cascade_scatter()
    fig_platform_concordance()
    fig_subject_fingerprint()
    print("=== Figures written ===")
    for p in sorted(FIG.glob("*.pdf")):
        print(f"  {p.name}")
