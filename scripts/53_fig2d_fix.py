"""Fix Fig 2D TMB/MSI waterfall — good/bad color visibility issue.

Root cause: all 41 matched tumors are MSS (max MSI 0.19 %, 35 exactly 0.0 %).
The original plot set y-axis 0–25 % to display the MSI-H 20 % cutoff line, which
made every bar near-zero and visually indistinguishable — so the legend said
good/poor but the bars were too small to show the color difference.

Fix: rescale y-axis to the actual data range (0–0.25 %) with a broken-axis
visual note that MSI-H cutoff is at 20 % (out of frame). Add a compact response
color strip above the bars so every sample is explicitly good/bad regardless of
bar height, and set bar heights to a minimum visible floor.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# text-as-text for PDF/SVG → downstream LibreOffice PPT keeps strings editable
plt.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.linewidth': 0.6,
})

BASE = Path('/mnt/sda1/data/TNT/analysis')
OUT  = Path('/data/data/TNT/analysis/figures/panels_v3')

GOOD = '#2E86AB'; BAD = '#E63946'; TMB_COLOR = '#9467bd'

# Load MSI + TMB (match the original script's merging logic)
msi_files = list((BASE/'02_wes_tmb_msi/msi').glob('*.tsv'))
msi = pd.concat([pd.read_csv(f, sep='\t') for f in msi_files])
tmb = pd.read_csv(BASE/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
mer = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id')

# matched tumors only — pre-CRT + post-CRT (drop normals and unmatched)
# UNMATCHED list from original script: subjects in tumor-only set
UNMATCHED = {13, 15, 16, 17, 18, 19, 33}
mer = mer[(mer.timepoint != 'normal') & (~mer.subject_id.isin(UNMATCHED))].copy()
mer = mer.sort_values('MSI_pct', ascending=False).reset_index(drop=True)
print(f'Matched tumor samples plotted: {len(mer)}')
print(mer['response_bin'].value_counts())

fig = plt.figure(figsize=(13, 5.2))
gs = fig.add_gridspec(2, 1, height_ratios=[0.04, 1.0], hspace=0.02)

# -------- Response color strip (top) --------
ax_strip = fig.add_subplot(gs[0])
strip_colors = [GOOD if r == 'good' else BAD for r in mer.response_bin]
for i, c in enumerate(strip_colors):
    ax_strip.add_patch(Rectangle((i - 0.5, 0), 1, 1, color=c, edgecolor='white', linewidth=0.3))
ax_strip.set_xlim(-0.6, len(mer) - 0.4)
ax_strip.set_ylim(0, 1)
ax_strip.set_xticks([]); ax_strip.set_yticks([])
for sp in ax_strip.spines.values(): sp.set_visible(False)
ax_strip.text(-1.8, 0.5, 'Response', ha='right', va='center', fontsize=9,
              fontweight='bold', color='#1d3557')

# -------- Main panel (bars + TMB overlay) --------
ax1 = fig.add_subplot(gs[1])

x = np.arange(len(mer))
# Give a visible minimum height so the color shows even when MSI ≈ 0
MIN_BAR = 0.003  # % (visual floor; real MSI values range 0–0.19)
plotted = np.maximum(mer.MSI_pct.values, MIN_BAR)
colors = [GOOD if r == 'good' else BAD for r in mer.response_bin]

bars = ax1.bar(x, plotted, color=colors, edgecolor='white', linewidth=0.6, alpha=0.95)

# Annotate actual MSI value above any non-zero bar
for i, (v, real) in enumerate(zip(plotted, mer.MSI_pct.values)):
    if real > 0.005:
        ax1.annotate(f'{real:.2f}', (i, v), xytext=(0, 3),
                     textcoords='offset points',
                     ha='center', va='bottom', fontsize=7.5, color='#333')

# Cutoff note (20% MSI-H is way above our y-max; show it as a callout box)
ax1.text(0.98, 0.96,
    'MSI-H clinical cutoff = 20 %\n(all samples shown here are ≤ 0.19 %; axis rescaled)',
    transform=ax1.transAxes, ha='right', va='top',
    fontsize=9, color='#d62828', fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
              edgecolor='#d62828', linewidth=0.8))

ax1.set_ylabel('MSI percentage (%)', color='#1d3557', fontsize=10)
ax1.set_ylim(0, 0.25)
ax1.set_xlim(-0.6, len(mer) - 0.4)
ax1.set_xticks(x)
ax1.set_xticklabels(mer.sample_id, rotation=90, fontsize=6.5)
ax1.tick_params(labelsize=8)
ax1.spines[['top','right']].set_visible(False)

# Secondary axis: TMB
ax2 = ax1.twinx()
ax2.scatter(x, mer.TMB_nonsyn_per_Mb, color=TMB_COLOR, s=38, edgecolor='white',
            linewidth=0.5, alpha=0.88, zorder=3)
ax2.axhline(10, color=TMB_COLOR, ls=':', lw=1, alpha=0.6)
ax2.text(len(mer)*0.01, 10.2, 'TMB 10/Mb', color=TMB_COLOR, fontsize=8.5,
         va='bottom', ha='left')
ax2.set_ylabel('TMB (nonsynonymous mutations / Mb)', color=TMB_COLOR, fontsize=10)
ax2.tick_params(labelcolor=TMB_COLOR, labelsize=8)
ax2.spines[['top']].set_visible(False)
ax2.spines['right'].set_color(TMB_COLOR)

# Legend (outside, bottom)
legend_handles = [
    mpatches.Patch(color=GOOD, label='Good (TRG 0–1)'),
    mpatches.Patch(color=BAD, label='Poor (TRG 2–3)'),
    matplotlib.lines.Line2D([0], [0], marker='o', color='w',
                            markerfacecolor=TMB_COLOR, markersize=8,
                            label='TMB (right axis)')
]
fig.legend(handles=legend_handles, loc='lower center',
           bbox_to_anchor=(0.5, -0.08),
           fontsize=10, frameon=False, ncol=3)

fig.suptitle('Sample-level MSI percentage (waterfall) with TMB overlay — all matched tumors are MSS (MSI < 1 %) and TMB-low (< 10 / Mb)',
             fontsize=11, fontweight='bold', color='#1d3557', y=1.00)

for ext in ('pdf', 'png'):
    plt.savefig(OUT/f'Fig2D_TMB_MSI_waterfall.{ext}', dpi=600, bbox_inches='tight')
print(f'Saved {OUT}/Fig2D_TMB_MSI_waterfall.pdf/png')
plt.close()
