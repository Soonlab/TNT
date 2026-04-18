#!/usr/bin/env python3
"""
17_figures_pharmacodynamics.py

Figures for the pharmacodynamics (script 15) and IGHV directional-consistency
(script 16) analyses.

    Fig A: baseline 4-factor spaghetti plot (pre->post per subject, colored by
           response), one panel per factor, with annotated sign counts.
    Fig B: composite sign bar (factor x group, showing n_predicted / n_total).
    Fig C: IGHV repertoire coherence summary -- paired scatter of
           good_majority_frac vs bad_majority_frac for all 53 V-genes, plus
           pattern pie.
    Fig D: focus V-gene panel (IGHV6-1, IGHV3-7, IGHV3-74, and top 3 other
           good_coherent_bad_mixed genes) -- spaghetti in fraction space.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

OUT = "/data/data/TNT/analysis/260418_add"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 120,
})
GOOD = "#2E86AB"
BAD = "#E63946"

# ---------------------------------------------------------------------------
# Fig A: spaghetti for 4 baseline factors (composite level)
# ---------------------------------------------------------------------------
w = pd.read_csv(f"{OUT}/baseline_factor_per_subject_delta.tsv", sep="\t")
sign = pd.read_csv(f"{OUT}/baseline_factor_sign_table.tsv", sep="\t")

factors = ["DSB_HDR_repair", "Tumor_cellcycle", "E2F_MYC_cellcycle", "EMT"]
preds = {"DSB_HDR_repair": "down", "Tumor_cellcycle": "down",
         "E2F_MYC_cellcycle": "down", "EMT": "up"}

fig, axes = plt.subplots(1, 4, figsize=(11.8, 3.3), sharey=False)
for ax, fname in zip(axes, factors):
    sub = w[(w.factor == fname) & (w.member == "composite")]
    for _, row in sub.iterrows():
        c = GOOD if row.response_bin == "good" else BAD
        ax.plot([0, 1], [row.pre, row.post], color=c, alpha=0.75, lw=1.2,
                marker="o", mfc="white", mec=c, mew=1.0, ms=4.5)
    # annotate sign counts
    sg = sign[(sign.factor == fname) & (sign.group == "good")].iloc[0]
    sb = sign[(sign.factor == fname) & (sign.group == "bad")].iloc[0]
    pred_arrow = "↓" if preds[fname] == "down" else "↑"
    txt = (f"good {pred_arrow}: {int(sg.n_predicted)}/{int(sg.n_total)}  "
           f"(P={sg.sign_binomial_one_sided_P:.3f})\n"
           f"bad  {pred_arrow}: {int(sb.n_predicted)}/{int(sb.n_total)}  "
           f"(P={sb.sign_binomial_one_sided_P:.3f})")
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=7.5,
            ha="left", va="bottom",
            bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5,
                      boxstyle="round,pad=0.3"))
    ax.set_title(f"{fname}  (pred: {preds[fname]})", fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["pre", "post"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel("composite z-score" if fname == factors[0] else "")
    ax.axhline(0, color="0.8", lw=0.5, zorder=0)

# legend
good_line = plt.Line2D([], [], color=GOOD, lw=1.2, marker="o", mfc="white", mec=GOOD,
                       label="good (n=6)")
bad_line = plt.Line2D([], [], color=BAD, lw=1.2, marker="o", mfc="white", mec=BAD,
                      label="bad (n=6)")
fig.legend(handles=[good_line, bad_line], loc="upper center",
           ncol=2, bbox_to_anchor=(0.5, 1.02), frameon=False)
fig.suptitle("Paired pre -> post RT-phase biopsy: 4 baseline factor trajectories (n=12)",
             y=1.07, fontsize=10)
fig.tight_layout()
for ext in ["pdf", "png"]:
    fig.savefig(f"{OUT}/FigE_baseline_spaghetti.{ext}",
                dpi=400 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("wrote FigE_baseline_spaghetti.{pdf,png}")

# ---------------------------------------------------------------------------
# Fig B: composite sign bar (stacked)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 3.2))
x = np.arange(len(factors))
width = 0.36

for i, grp in enumerate(["good", "bad"]):
    color = GOOD if grp == "good" else BAD
    y = [sign[(sign.factor == f) & (sign.group == grp)]
         ["fraction_predicted"].values[0] for f in factors]
    labels = [f"{int(sign[(sign.factor == f) & (sign.group == grp)]['n_predicted'].values[0])}/6"
              for f in factors]
    bars = ax.bar(x + (i - 0.5) * width, y, width, color=color, alpha=0.85,
                  label=f"{grp} (n=6)")
    for b, lab in zip(bars, labels):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02,
                lab, ha="center", fontsize=7.5)

ax.axhline(0.5, color="0.5", ls="--", lw=0.6, zorder=0, label="chance (0.5)")
ax.set_xticks(x); ax.set_xticklabels(factors, rotation=20, ha="right")
ax.set_ylim(0, 1.12)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, decimals=0))
ax.set_ylabel("Fraction of subjects moving in predicted direction")
ax.set_title("Baseline factor directional concordance (paired n=6+6)", fontsize=9.5)
ax.legend(loc="upper right", frameon=False)
fig.tight_layout()
for ext in ["pdf", "png"]:
    fig.savefig(f"{OUT}/FigF_baseline_sign_bar.{ext}",
                dpi=400 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("wrote FigF_baseline_sign_bar.{pdf,png}")

# ---------------------------------------------------------------------------
# Fig C: IGHV repertoire coherence -- paired scatter + pattern pie
# ---------------------------------------------------------------------------
vstats = pd.read_csv(f"{OUT}/trust4_ighv_directional_stats.tsv", sep="\t")

fig = plt.figure(figsize=(10.5, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.1, 0.9], wspace=0.4)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

# (1) scatter good_majority_frac vs bad_majority_frac
for _, r in vstats.iterrows():
    color = "0.7"
    marker = "o"
    if r.pattern == "good_coherent_bad_mixed":
        color = GOOD; marker = "o"
    elif r.pattern == "bad_coherent_good_mixed":
        color = BAD; marker = "o"
    ax1.scatter(r.bad_majority_frac, r.good_majority_frac,
                color=color, alpha=0.85, s=28, lw=0.5, edgecolor="white",
                zorder=3)
lims = (0.45, 1.03)
ax1.plot(lims, lims, "--", color="0.5", lw=0.6, zorder=0)
ax1.set_xlim(lims); ax1.set_ylim(lims)
ax1.set_xlabel("bad responder V-gene majority fraction")
ax1.set_ylabel("good responder V-gene majority fraction")
ax1.set_title("V-gene directional coherence, good vs bad", fontsize=9.5)
# annotate the IGHV6-1 and user-focus genes
for vg in ["IGHV6-1", "IGHV3-7", "IGHV3-74", "IGHV1-45", "IGHV3-38-3"]:
    r = vstats[vstats.v_gene == vg]
    if r.empty:
        continue
    r = r.iloc[0]
    ax1.annotate(vg, (r.bad_majority_frac, r.good_majority_frac),
                 xytext=(4, 4), textcoords="offset points", fontsize=7)

# (2) Wilcoxon summary --- boxplot of paired differences
diffs = (vstats.good_majority_frac - vstats.bad_majority_frac).dropna()
ax2.axhline(0, color="0.5", lw=0.6)
ax2.boxplot([diffs], widths=0.4, patch_artist=True,
            boxprops=dict(facecolor="#E3E9F2", edgecolor="0.3", lw=0.8),
            medianprops=dict(color="#222", lw=1.2),
            whiskerprops=dict(color="0.3", lw=0.6),
            capprops=dict(color="0.3", lw=0.6),
            flierprops=dict(marker="o", mfc="white", mec="0.4",
                            ms=3.5, lw=0.4))
# scatter
jitter = np.random.default_rng(0).normal(0, 0.035, len(diffs))
ax2.scatter(1 + jitter, diffs, color="#4C6DA6", alpha=0.65, s=18, lw=0, zorder=4)
from scipy import stats as sst
w_stat, w_p = sst.wilcoxon(diffs, alternative="greater")
ax2.text(1, diffs.max() + 0.03,
         f"Wilcoxon one-sided P = {w_p:.3f}\n(53 V-genes)",
         ha="center", fontsize=8,
         bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5,
                   boxstyle="round,pad=0.3"))
ax2.set_xticks([1]); ax2.set_xticklabels(["good - bad\n(majority_frac)"])
ax2.set_ylabel("Δ majority fraction (paired per V-gene)")
ax2.set_title("Aggregate coherence gap", fontsize=9.5)

# (3) pattern pie
pc = vstats.pattern.value_counts()
# order for readability
order = ["good_coherent_bad_mixed", "bad_coherent_good_mixed",
         "both_coherent_same", "both_coherent_opposite", "both_mixed"]
pc = pc.reindex([k for k in order if k in pc.index])
palette = {
    "good_coherent_bad_mixed": GOOD,
    "bad_coherent_good_mixed": BAD,
    "both_coherent_same": "#7FB069",
    "both_coherent_opposite": "#C5A572",
    "both_mixed": "0.8",
}
colors = [palette[k] for k in pc.index]
wedges, texts, autotexts = ax3.pie(pc.values, colors=colors,
                                   autopct=lambda p: f"{p*sum(pc.values)/100:.0f}",
                                   startangle=90,
                                   wedgeprops=dict(edgecolor="white", lw=1))
for t in autotexts:
    t.set_fontsize(8)
ax3.set_title("Pattern breakdown\n(53 V-genes)", fontsize=9.5)
ax3.legend(wedges, pc.index, loc="center left",
           bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)

fig.suptitle("IGHV repertoire: good responders show more coherent RT-phase changes",
             y=1.02, fontsize=10)
fig.tight_layout()
for ext in ["pdf", "png"]:
    fig.savefig(f"{OUT}/FigG_ighv_coherence_summary.{ext}",
                dpi=400 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("wrote FigG_ighv_coherence_summary.{pdf,png}")

# ---------------------------------------------------------------------------
# Fig D: focus V-gene spaghetti
# ---------------------------------------------------------------------------
wide_ighv = pd.read_csv(f"{OUT}/trust4_ighv_per_subject_delta.tsv", sep="\t")

# IGHV6-1 (strongest); user-prior IGHV3-7, IGHV3-74; top 3 by Fisher P
focus_user = ["IGHV6-1", "IGHV3-7", "IGHV3-74"]
top_by_fisher = vstats[vstats.pattern == "good_coherent_bad_mixed"].sort_values(
    "fisher_P_updown").head(6)["v_gene"].tolist()
focus_list = []
for vg in focus_user + top_by_fisher:
    if vg not in focus_list:
        focus_list.append(vg)
focus_list = focus_list[:9]  # cap to 9 for 3x3

ncols = 3
nrows = int(np.ceil(len(focus_list) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 2.8 * nrows),
                         sharey=False)
axes_flat = axes.flat if nrows * ncols > 1 else [axes]

for ax, vg in zip(axes_flat, focus_list):
    sub = wide_ighv[wide_ighv.v_gene == vg]
    for _, row in sub.iterrows():
        c = GOOD if row.response_bin == "good" else BAD
        ax.plot([0, 1], [row.pre, row.post], color=c, alpha=0.75, lw=1.1,
                marker="o", mfc="white", mec=c, mew=0.9, ms=4.2)
    stat = vstats[vstats.v_gene == vg].iloc[0]
    txt = (f"good ↑/↓: {int(stat.good_n_up)}/{int(stat.good_n_down)}\n"
           f"bad ↑/↓: {int(stat.bad_n_up)}/{int(stat.bad_n_down)}\n"
           f"Fisher P={stat.fisher_P_updown:.3f}  MW P={stat.mw_P_delta:.3f}")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, fontsize=7,
            ha="left", va="top",
            bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5,
                      boxstyle="round,pad=0.25"))
    ax.set_title(vg, fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["pre", "post"])
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylabel("fraction of IGH repertoire")

# hide empty panels
for ax in list(axes_flat)[len(focus_list):]:
    ax.set_visible(False)

good_line = plt.Line2D([], [], color=GOOD, lw=1.2, marker="o", mfc="white", mec=GOOD,
                       label="good (n=6)")
bad_line = plt.Line2D([], [], color=BAD, lw=1.2, marker="o", mfc="white", mec=BAD,
                      label="bad (n=6)")
fig.legend(handles=[good_line, bad_line], loc="upper center",
           ncol=2, bbox_to_anchor=(0.5, 1.01), frameon=False)
fig.suptitle("Focus IGHV genes: pre -> post fraction trajectories",
             y=1.03, fontsize=10)
fig.tight_layout()
for ext in ["pdf", "png"]:
    fig.savefig(f"{OUT}/FigH_ighv_focus_spaghetti.{ext}",
                dpi=400 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("wrote FigH_ighv_focus_spaghetti.{pdf,png}")

print("\nAll figures written to 260418_add/")
