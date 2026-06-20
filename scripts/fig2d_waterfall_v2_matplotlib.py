"""Fig 2D v2 — TMB + MSI waterfall, response-grouped with paired adjacency.

Redesign rationale:
  - Previous ordering (MSI desc) clustered all MSS samples indistinguishably at
    near-zero and broke response-group legibility.
  - New ordering: Poor group (left) | vertical separator | Good group (right).
    Within each group, paired subjects come first (PR immediately adjacent to
    PO for visual pre→post comparison), then single-pre unpaired subjects.
    Sort by subject_id ascending; within subject PR→PO.

Output: panels_v3/Fig2D_TMB_MSI_waterfall.{pdf,png}  (overwrites old)
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

plt.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype':'none',
    'font.family':'DejaVu Sans', 'font.size':9, 'axes.linewidth':0.6,
})

BASE = Path('/mnt/sda1/data/TNT/analysis')
OUT  = BASE/'figures/panels_v3'
OUT_SUB = BASE/'genome_medicine_submission/main_figures'

GOOD='#2E86AB'; BAD='#E63946'; TMB_COLOR='#6a4c93'
PRE='#457b9d'; POST='#f4a261'
SEP='#1d3557'

UNMATCHED = {13,15,16,17,18,19,33}

msi_files = list((BASE/'02_wes_tmb_msi/msi').glob('*.tsv'))
msi = pd.concat([pd.read_csv(f, sep='\t') for f in msi_files])
tmb = pd.read_csv(BASE/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
mer = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id')
mer = mer[(mer.timepoint!='normal') & (~mer.subject_id.isin(UNMATCHED))].copy()

def build_order(sub):
    """For one response group: paired subjects (PR→PO) asc, then single subjects asc."""
    paired = sub.groupby('subject_id').filter(lambda x: len(x)==2)
    paired_subj = sorted(paired.subject_id.unique())
    single_subj = sorted([s for s in sub.subject_id.unique() if s not in paired_subj])
    order = []
    for s in paired_subj:
        rows = sub[sub.subject_id==s].sort_values('timepoint',
                    key=lambda x: x.map({'pre':0,'post':1}))
        order.extend(rows.sample_id.tolist())
    for s in single_subj:
        order.extend(sub[sub.subject_id==s].sample_id.tolist())
    return order, len(paired_subj)*2  # count of paired samples (for separator position)

bad_order, bad_paired_n = build_order(mer[mer.response_bin=='bad'])
good_order, good_paired_n = build_order(mer[mer.response_bin=='good'])
order = bad_order + good_order
sep_x = len(bad_order) - 0.5  # between bad and good
bad_sub_sep = bad_paired_n - 0.5
good_sub_sep = len(bad_order) + good_paired_n - 0.5

df = mer.set_index('sample_id').loc[order].reset_index()
print(f'Plot order: {len(df)} samples')
print(f'  Bad n={len(bad_order)} (paired {bad_paired_n}, single {len(bad_order)-bad_paired_n})')
print(f'  Good n={len(good_order)} (paired {good_paired_n}, single {len(good_order)-good_paired_n})')

x = np.arange(len(df))
MIN_BAR = 0.003
plotted = np.maximum(df.MSI_pct.values, MIN_BAR)
resp_colors = [BAD if r=='bad' else GOOD for r in df.response_bin]
tp_colors   = [PRE if t=='pre' else POST for t in df.timepoint]

fig = plt.figure(figsize=(14, 5.6))
gs = fig.add_gridspec(3, 1, height_ratios=[0.05, 0.05, 1.0], hspace=0.04)

# Response strip
ax_r = fig.add_subplot(gs[0])
for i, c in enumerate(resp_colors):
    ax_r.add_patch(Rectangle((i-0.5, 0), 1, 1, color=c, edgecolor='white', linewidth=0.3))
ax_r.set_xlim(-0.6, len(df)-0.4); ax_r.set_ylim(0,1); ax_r.set_xticks([]); ax_r.set_yticks([])
for sp in ax_r.spines.values(): sp.set_visible(False)
ax_r.text(-1.4, 0.5, 'Response', ha='right', va='center', fontsize=9, fontweight='bold', color=SEP)

# Timepoint strip
ax_t = fig.add_subplot(gs[1])
for i, c in enumerate(tp_colors):
    ax_t.add_patch(Rectangle((i-0.5, 0), 1, 1, color=c, edgecolor='white', linewidth=0.3))
ax_t.set_xlim(-0.6, len(df)-0.4); ax_t.set_ylim(0,1); ax_t.set_xticks([]); ax_t.set_yticks([])
for sp in ax_t.spines.values(): sp.set_visible(False)
ax_t.text(-1.4, 0.5, 'Timepoint', ha='right', va='center', fontsize=9, fontweight='bold', color=SEP)

# Main MSI bars + TMB
ax1 = fig.add_subplot(gs[2])
bars = ax1.bar(x, plotted, color=resp_colors, edgecolor='white', linewidth=0.5, alpha=0.95)
for i, real in enumerate(df.MSI_pct.values):
    if real > 0.005:
        ax1.annotate(f'{real:.2f}', (i, max(real, MIN_BAR)), xytext=(0,3),
                     textcoords='offset points', ha='center', va='bottom',
                     fontsize=7.5, color='#333')

# Group separator (bad vs good)
ax1.axvline(sep_x, color=SEP, ls='-', lw=1.6, alpha=0.8)
# Sub-separators (paired vs single within each group)
ax1.axvline(bad_sub_sep,  color=SEP, ls=':', lw=0.8, alpha=0.5)
ax1.axvline(good_sub_sep, color=SEP, ls=':', lw=0.8, alpha=0.5)

# Group labels just below response strip, above all bars
ax1.text((sep_x)/2, 0.235, 'Poor (TRG 2–3)', ha='center', va='top',
         fontsize=11, fontweight='bold', color=BAD,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                   edgecolor=BAD, linewidth=0.8))
ax1.text((sep_x + len(df)-0.5)/2, 0.235, 'Good (TRG 0–1)', ha='center', va='top',
         fontsize=11, fontweight='bold', color=GOOD,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                   edgecolor=GOOD, linewidth=0.8))

# MSI cutoff note — placed in upper area where all bars are well below y-axis top
ax1.text(0.50, 0.62,
    'MSI-H clinical cutoff = 20 %   (all samples shown ≤ 0.19 %; axis rescaled)',
    transform=ax1.transAxes, ha='center', va='center',
    fontsize=9, color='#d62828', fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
              edgecolor='#d62828', linewidth=0.8))

ax1.set_ylabel('MSI percentage (%)', color=SEP, fontsize=10)
ax1.set_ylim(0, 0.25)
ax1.set_xlim(-0.6, len(df)-0.4)
ax1.set_xticks(x); ax1.set_xticklabels(df.sample_id, rotation=90, fontsize=6.8)
ax1.tick_params(labelsize=8)
ax1.spines[['top','right']].set_visible(False)

# TMB secondary axis
ax2 = ax1.twinx()
ax2.scatter(x, df.TMB_nonsyn_per_Mb, color=TMB_COLOR, s=38, edgecolor='white',
            linewidth=0.5, alpha=0.88, zorder=3)
ax2.axhline(10, color=TMB_COLOR, ls=':', lw=1, alpha=0.6)
ax2.text(0.2, 10.2, 'TMB 10/Mb cutoff', color=TMB_COLOR, fontsize=8.5, va='bottom', ha='left')
ax2.set_ylabel('TMB (nonsynonymous / Mb)', color=TMB_COLOR, fontsize=10)
ax2.tick_params(labelcolor=TMB_COLOR, labelsize=8)
ax2.spines[['top']].set_visible(False)
ax2.spines['right'].set_color(TMB_COLOR)

# Legend
legend_handles = [
    mpatches.Patch(color=GOOD, label='Good (TRG 0–1)'),
    mpatches.Patch(color=BAD,  label='Poor (TRG 2–3)'),
    mpatches.Patch(color=PRE,  label='Pre-CRT'),
    mpatches.Patch(color=POST, label='Post-CRT'),
    matplotlib.lines.Line2D([0],[0], marker='o', color='w',
                            markerfacecolor=TMB_COLOR, markersize=8,
                            label='TMB (right axis)'),
    matplotlib.lines.Line2D([0],[0], color=SEP, lw=1.6, label='Poor | Good group boundary'),
    matplotlib.lines.Line2D([0],[0], color=SEP, ls=':', lw=0.9,
                            label='Paired | single boundary'),
]
fig.legend(handles=legend_handles, loc='lower center',
           bbox_to_anchor=(0.5, -0.14), fontsize=9.5, frameon=False, ncol=4)

fig.suptitle('Sample-level MSI and TMB waterfall — Poor (left) vs Good (right); paired PR/PO adjacent',
             fontsize=11, fontweight='bold', color=SEP, y=1.00)

for d in (OUT, OUT_SUB):
    d.mkdir(parents=True, exist_ok=True)
    plt.savefig(d/f'Fig2D_TMB_MSI_waterfall.pdf', bbox_inches='tight')
    plt.savefig(d/f'Fig2D_TMB_MSI_waterfall.png', dpi=600, bbox_inches='tight')
print(f'Saved: {OUT/"Fig2D_TMB_MSI_waterfall.pdf/png"}  (+submission mirror)')
plt.close()
