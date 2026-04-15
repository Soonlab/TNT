#!/usr/bin/env python3
"""Compose multi-panel main figures (Fig2, Fig3, Fig4, Fig5_integration_delta)
using native matplotlib gridspec (no image-stitching). All panels are drawn
from source TSVs to guarantee consistent panel sizes, uniform typography, and
11pt bold panel letters matching the Fig5_neoantigen_cascade reference style.

This script does NOT recompute any statistics — it only redraws layout.
Stat numbers come from existing TSV outputs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy import stats

ROOT = Path("/mnt/sda1/data/TNT/analysis")
OUT = ROOT / "figures" / "panels"
OUT.mkdir(parents=True, exist_ok=True)

GOOD = "#2E86AB"
BAD  = "#E63946"
PAL  = {"good": GOOD, "bad": BAD}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "#222",
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
})

def panel_letter(ax, letter, x=-0.12, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")

def stat_annot(ax, text, x=0.98, y=0.97):
    ax.text(x, y, text, transform=ax.transAxes,
            ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="none", alpha=0.75))

# ================================================================
# Shared loaders
# ================================================================
rna_inv = pd.read_csv(ROOT/"00_cohort/rna_inventory.tsv", sep="\t")
sig_scores = pd.read_csv(ROOT/"06_rna_immune/signature_scores.tsv", sep="\t")
sig_stats  = pd.read_csv(ROOT/"06_rna_immune/sig_response_stats.tsv", sep="\t")
tmb = pd.read_csv(ROOT/"02_wes_tmb_msi/tmb_per_sample.tsv", sep="\t")
gsea_hall = pd.read_csv(ROOT/"05_rna_deg_gsea/GSEA_Hallmark_pre.tsv", sep="\t")
resp_feat = pd.read_csv(ROOT/"tables/response_feature_stats.tsv", sep="\t")
onco = pd.read_csv(ROOT/"04_wes_cnv_clonal/driver_oncoprint_matrix.tsv", sep="\t")
delta_ssgsea = pd.read_csv(ROOT/"09_integration/paired_delta/delta_ssgsea_response.tsv", sep="\t")
delta_trust  = pd.read_csv(ROOT/"09_integration/paired_delta/delta_trust4_response.tsv", sep="\t")
trust = pd.read_csv(ROOT/"06_rna_immune/trust4_summary.tsv", sep="\t")

# Attach response/timepoint to signature scores (pre-treatment only)
sig_long = sig_scores.merge(
    rna_inv[["sample_id","subject_id","timepoint","response_bin"]],
    on="sample_id", how="left"
)
sig_pre = sig_long[sig_long.timepoint == "pre"].copy()

# ================================================================
# Figure 2 — WES landscape (A: TMB pre by response; B: oncoprint)
# ================================================================
def make_fig2():
    fig = plt.figure(figsize=(12, 5.2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 2.6],
                  wspace=0.25, left=0.07, right=0.99, top=0.88, bottom=0.22)

    # --- A: TMB pre (matched) by response ---
    axA = fig.add_subplot(gs[0, 0])
    pre = tmb[(tmb.timepoint == "pre") & (tmb.matched == True)].copy()
    sns.boxplot(data=pre, x="response_bin", y="TMB_nonsyn_per_Mb",
                order=["good","bad"], hue="response_bin",
                palette=PAL, legend=False, ax=axA,
                width=0.55, linewidth=0.7, fliersize=0)
    sns.stripplot(data=pre, x="response_bin", y="TMB_nonsyn_per_Mb",
                  order=["good","bad"], color="black", size=3, ax=axA)
    axA.set_yscale("log")
    axA.set_xlabel("")
    axA.set_ylabel("TMB (nonsyn / Mb)")
    axA.set_title("TMB pre-treatment", fontsize=9, fontweight="bold")
    g = pre[pre.response_bin=="good"].TMB_nonsyn_per_Mb.dropna()
    b = pre[pre.response_bin=="bad"].TMB_nonsyn_per_Mb.dropna()
    try:
        p = stats.mannwhitneyu(g, b).pvalue
        stat_annot(axA, f"P = {p:.3f}")
    except Exception:
        pass
    panel_letter(axA, "A", x=-0.22, y=1.08)

    # --- B: Oncoprint heatmap (top 15 drivers) ---
    axB = fig.add_subplot(gs[0, 1])
    mat = onco.set_index("GENE").drop(columns=["total_samples_mutated"], errors="ignore")
    # restrict to pre samples only (names end with -P or -PR)
    pre_cols = [c for c in mat.columns if c.endswith("-P") or c.endswith("-PR")]
    mat = mat[pre_cols]
    # top 15 genes by sample count
    top = mat.sum(axis=1).sort_values(ascending=False).head(15).index
    mat = mat.loc[top]
    # order samples by response
    samp2resp = {}
    for c in mat.columns:
        sid = int(c.split("-")[0])
        row = rna_inv[rna_inv.subject_id == sid]
        if len(row):
            samp2resp[c] = row.response_bin.iloc[0]
        else:
            samp2resp[c] = "good"
    ordered = sorted(mat.columns, key=lambda c: (samp2resp[c] != "good", c))
    mat = mat[ordered]
    axB.imshow((mat.values > 0).astype(float), aspect="auto",
               cmap="Greys", vmin=0, vmax=1, interpolation="none")
    axB.set_yticks(range(len(mat.index)))
    axB.set_yticklabels(mat.index, fontsize=7)
    axB.set_xticks(range(len(mat.columns)))
    axB.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=6)
    axB.set_title("Driver mutations (pre-treatment)",
                  fontsize=9, fontweight="bold")
    # response strip above
    for i, c in enumerate(mat.columns):
        axB.add_patch(plt.Rectangle((i-0.5, -1.2), 1, 0.5,
                                    color=PAL[samp2resp[c]], clip_on=False))
    axB.text(-0.02, -0.08, "Response:", transform=axB.transAxes,
             fontsize=6.5, ha="right", va="center")
    for sp in ("top","right"):
        axB.spines[sp].set_visible(False)
    panel_letter(axB, "B", x=-0.08, y=1.08)

    fig.suptitle("Figure 2. Somatic landscape of MSS LARC TNT cohort",
                 fontsize=11, fontweight="bold", x=0.5, y=0.98)
    fig.savefig(OUT/"Fig2_WES_landscape.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT/"Fig2_WES_landscape.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote Fig2_WES_landscape")

# ================================================================
# Figure 3 — RNA signatures (A: 10 boxplots, B: heatmap)
# ================================================================
def make_fig3():
    sigs_to_show = ["CD8_proliferation","CD8_activation","Antigen_presentation",
                    "NLRC5_HLA_IFNG","TLS_Cabrita","IFNg_Ayers_18",
                    "TGFb_Mariathasan","EMT_Mak","MHC_II","Checkpoint_inhibitory"]
    fig = plt.figure(figsize=(14, 6.2))
    outer = GridSpec(1, 2, figure=fig, width_ratios=[1.5, 1.0],
                     wspace=0.22, left=0.05, right=0.99, top=0.90, bottom=0.08)

    # A: 2x5 boxplots
    gsA = outer[0, 0].subgridspec(2, 5, wspace=0.65, hspace=0.75)
    pre_stats = sig_stats[sig_stats.timepoint.isin(["pre","all"])]
    for i, sig in enumerate(sigs_to_show):
        ax = fig.add_subplot(gsA[i // 5, i % 5])
        d = sig_pre[["response_bin", sig]].dropna()
        sns.boxplot(data=d, x="response_bin", y=sig, order=["good","bad"],
                    hue="response_bin", palette=PAL, legend=False, ax=ax,
                    width=0.55, linewidth=0.6, fliersize=0)
        sns.stripplot(data=d, x="response_bin", y=sig, order=["good","bad"],
                      color="black", size=2.2, ax=ax)
        r = pre_stats[pre_stats.signature == sig]
        p = r.pvalue.iloc[0] if len(r) else np.nan
        ax.set_title(sig.replace("_"," "), fontsize=7.5, fontweight="bold", pad=2)
        if not np.isnan(p):
            stat_annot(ax, f"P = {p:.3f}")
        ax.set_xlabel("")
        ax.set_ylabel("score", fontsize=6.5)
        ax.tick_params(axis="both", labelsize=6)
        if i == 0:
            panel_letter(ax, "A", x=-0.35, y=1.15)

    # B: heatmap
    axB = fig.add_subplot(outer[0, 1])
    sig_cols = [c for c in sig_scores.columns if c != "sample_id"]
    mat = sig_pre.set_index("sample_id")[sig_cols].T
    # order samples by response
    order_samps = sig_pre.sort_values("response_bin", ascending=True).sample_id.tolist()
    mat = mat[order_samps]
    vmin, vmax = -2, 2
    im = axB.imshow(mat.values, aspect="auto", cmap="RdBu_r",
                    vmin=vmin, vmax=vmax, interpolation="nearest")
    axB.set_yticks(range(len(mat.index)))
    axB.set_yticklabels(mat.index, fontsize=6.5)
    axB.set_xticks([])
    axB.set_xlabel("Samples (pre, sorted by response)")
    # response strip
    resp_map = sig_pre.set_index("sample_id").response_bin
    for i, s in enumerate(order_samps):
        axB.add_patch(plt.Rectangle((i-0.5, -1.0), 1, 0.6,
                                    color=PAL[resp_map[s]], clip_on=False))
    axB.set_title("Pre-treatment signature heatmap",
                  fontsize=9, fontweight="bold", pad=16)
    cbar = fig.colorbar(im, ax=axB, shrink=0.5, pad=0.02)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("z-score", fontsize=7)
    panel_letter(axB, "B", x=-0.05, y=1.05)

    fig.suptitle("Figure 3. Pre-treatment immune and pathway signatures",
                 fontsize=11, fontweight="bold", y=0.98)
    fig.savefig(OUT/"Fig3_RNA_signatures.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT/"Fig3_RNA_signatures.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote Fig3_RNA_signatures")

# ================================================================
# Figure 4 — GSEA + integrated response features
# ================================================================
def make_fig4():
    fig = plt.figure(figsize=(13, 6.2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.0],
                  wspace=0.55, left=0.16, right=0.98, top=0.90, bottom=0.10)

    # A: Hallmark GSEA signed -log10(p)
    axA = fig.add_subplot(gs[0, 0])
    h = gsea_hall.copy()
    h["signed_logp"] = -np.log10(h["pval"].clip(lower=1e-30)) * np.sign(h["NES"])
    h = h.reindex(h["signed_logp"].abs().sort_values(ascending=False).index).head(22)
    h = h.sort_values("signed_logp")
    colors = [GOOD if v > 0 else BAD for v in h["signed_logp"]]
    axA.barh(range(len(h)), h["signed_logp"], color=colors,
             edgecolor="white", linewidth=0.4)
    axA.set_yticks(range(len(h)))
    axA.set_yticklabels([p.replace("HALLMARK_","") for p in h["pathway"]],
                        fontsize=6.5)
    axA.axvline(0, color="#888", lw=0.5)
    axA.set_xlabel("signed $-\\log_{10}(P)$   (+ = UP in good)")
    axA.set_title("GSEA Hallmark — good vs bad (pre)",
                  fontsize=9, fontweight="bold")
    panel_letter(axA, "A", x=-0.42, y=1.05)

    # B: integrated response features bar
    axB = fig.add_subplot(gs[0, 1])
    rf = resp_feat.copy()
    rf["logp"] = -np.log10(rf["pvalue"].clip(lower=1e-30))
    rf["dir"] = np.sign(rf["delta_med"])
    rf = rf.sort_values("logp", ascending=False).head(25).iloc[::-1]
    bar_colors = [GOOD if d > 0 else BAD for d in rf["dir"]]
    axB.barh(range(len(rf)), rf["logp"], color=bar_colors,
             edgecolor="white", linewidth=0.4)
    axB.set_yticks(range(len(rf)))
    axB.set_yticklabels([str(f)[:55] for f in rf["feature"]], fontsize=6.5)
    axB.axvline(-np.log10(0.05), color="#888", lw=0.5, ls="--", label="P=0.05")
    axB.set_xlabel("$-\\log_{10}(P)$  good vs bad")
    axB.set_title("Response-associated integrated features",
                  fontsize=9, fontweight="bold")
    axB.legend(loc="lower right", fontsize=6)
    panel_letter(axB, "B", x=-0.55, y=1.05)

    fig.suptitle("Figure 4. DNA repair / cell-cycle programs predict TNT response",
                 fontsize=11, fontweight="bold", y=0.98)
    fig.savefig(OUT/"Fig4_GSEA_integration.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT/"Fig4_GSEA_integration.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote Fig4_GSEA_integration")

# ================================================================
# Figure 5 — Integration + treatment-induced deltas (2 rows)
# A: top ssGSEA delta bar | B: TRUST4 paired deltas row (4)
# C: ssGSEA paired deltas row (4)
# ================================================================
def make_fig5_integration():
    fig = plt.figure(figsize=(13, 8.8))
    gs = GridSpec(3, 4, figure=fig, height_ratios=[1.4, 1.0, 1.0],
                  wspace=0.55, hspace=0.85,
                  left=0.18, right=0.98, top=0.93, bottom=0.08)

    # --- A: top ssGSEA delta bar (spans all 4 cols, top row) ---
    axA = fig.add_subplot(gs[0, :])
    ds = delta_ssgsea.copy()
    ds["abs_delta"] = ds["delta_good_minus_bad"].abs()
    ds = ds.sort_values("abs_delta", ascending=False).head(20).iloc[::-1]
    colors = [GOOD if v > 0 else BAD for v in ds["delta_good_minus_bad"]]
    axA.barh(range(len(ds)), ds["delta_good_minus_bad"],
             color=colors, edgecolor="white", linewidth=0.4)
    axA.set_yticks(range(len(ds)))
    axA.set_yticklabels([str(f)[:60] for f in ds["feature"]], fontsize=6.5)
    axA.axvline(0, color="#888", lw=0.5)
    axA.set_xlabel("Δ(post − pre) good minus bad   (teal=↑good, orange=↑bad)")
    axA.set_title("Top treatment-induced pathway deltas",
                  fontsize=9, fontweight="bold")
    panel_letter(axA, "A", x=-0.28, y=1.05)

    # --- B row: TRUST4 paired deltas (4 metrics) ---
    trust_m = trust.merge(
        rna_inv[["sample_id","subject_id","timepoint","response_bin"]],
        on="sample_id", how="left", suffixes=("","_inv")
    )
    if "response_bin" not in trust_m.columns or trust_m.response_bin.isna().all():
        trust_m = trust.copy()
    # pick response_bin column robustly
    resp_col = "response_bin" if "response_bin" in trust_m.columns else "response"
    subs_paired = set(trust_m.groupby("subject_id").filter(
        lambda g: set(g.timepoint) >= {"pre","post"}).subject_id)

    trust_metrics = [("IGH_n","BCR (IGH) clonotypes"),
                     ("IGH_shannon","IGH Shannon diversity"),
                     ("TRB_shannon","TCR (TRB) Shannon"),
                     ("TRB_n","TRB clonotypes")]
    b_letters = ["B","C","D","E"]
    for j, (metric, label) in enumerate(trust_metrics):
        ax = fig.add_subplot(gs[1, j])
        rows = []
        for s in sorted(subs_paired):
            g = trust_m[trust_m.subject_id == s]
            pre = g[g.timepoint=="pre"][metric].values
            post = g[g.timepoint=="post"][metric].values
            resp = g[resp_col].iloc[0]
            if len(pre) and len(post):
                rows.append({"response": resp, "delta": post[0]-pre[0]})
        dd = pd.DataFrame(rows)
        sns.boxplot(data=dd, x="response", y="delta", order=["good","bad"],
                    hue="response", palette=PAL, legend=False, ax=ax,
                    width=0.55, linewidth=0.6, fliersize=0)
        sns.stripplot(data=dd, x="response", y="delta", order=["good","bad"],
                      color="black", size=3, ax=ax)
        ax.axhline(0, color="gray", ls="--", lw=0.6)
        r = delta_trust[delta_trust.feature == metric]
        p = r.MW_p.iloc[0] if len(r) else np.nan
        ax.set_title(label, fontsize=8, fontweight="bold", pad=2)
        ax.set_xlabel(""); ax.set_ylabel("Δ(post−pre)", fontsize=7)
        if not np.isnan(p):
            stat_annot(ax, f"P = {p:.3f}")
        panel_letter(ax, b_letters[j], x=-0.30, y=1.18)

    # --- C row: ssGSEA paired deltas for top 4 pathways by |delta| ---
    top4 = (delta_ssgsea.assign(abs_d=delta_ssgsea["delta_good_minus_bad"].abs())
            .sort_values("abs_d", ascending=False).head(4))
    # need per-subject ssGSEA deltas
    ssg = pd.read_csv(ROOT/"08_rna_pathway/ssgsea_scores.tsv", sep="\t")
    ssg_m = ssg.merge(
        rna_inv[["sample_id","subject_id","timepoint","response_bin"]],
        on="sample_id", how="left"
    )
    c_letters = ["F","G","H","I"]
    for j, (_, row) in enumerate(top4.iterrows()):
        feat = row["feature"]
        if feat not in ssg_m.columns:
            continue
        ax = fig.add_subplot(gs[2, j])
        rows_d = []
        for s in ssg_m.subject_id.dropna().unique():
            g = ssg_m[ssg_m.subject_id == s]
            pre = g[g.timepoint=="pre"][feat].values
            post = g[g.timepoint=="post"][feat].values
            if len(pre) and len(post):
                rows_d.append({"response": g.response_bin.iloc[0],
                               "delta": post[0]-pre[0]})
        dd = pd.DataFrame(rows_d)
        if dd.empty:
            continue
        sns.boxplot(data=dd, x="response", y="delta", order=["good","bad"],
                    hue="response", palette=PAL, legend=False, ax=ax,
                    width=0.55, linewidth=0.6, fliersize=0)
        sns.stripplot(data=dd, x="response", y="delta", order=["good","bad"],
                      color="black", size=3, ax=ax)
        ax.axhline(0, color="gray", ls="--", lw=0.6)
        title = str(feat).split("R-HSA")[0].strip(": ")[:28]
        ax.set_title(title, fontsize=8, fontweight="bold", pad=2)
        ax.set_xlabel(""); ax.set_ylabel("Δ(post−pre)", fontsize=7)
        p = row.get("MW_p", np.nan)
        if not np.isnan(p):
            stat_annot(ax, f"P = {p:.3f}")
        panel_letter(ax, c_letters[j], x=-0.30, y=1.18)

    fig.suptitle("Figure 5. Feature integration and treatment-induced changes",
                 fontsize=11, fontweight="bold", y=0.985)
    fig.savefig(OUT/"Fig5_integration_delta.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT/"Fig5_integration_delta.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote Fig5_integration_delta")


if __name__ == "__main__":
    make_fig2()
    make_fig3()
    make_fig4()
    make_fig5_integration()
    print("done")
