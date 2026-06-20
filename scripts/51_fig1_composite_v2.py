"""Compose the new Fig 1 (Option 2) = 4 panels:
  A Study design (wide, top)
  B Sankey (sex→cT→response)
  C Sample matrix (WES × RNA × timepoint)
  D Headline 3-narrative preview forest (wide, bottom)
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PV3  = Path('/data/data/TNT/analysis/figures/panels_v3')
MAIN = Path('/data/data/TNT/analysis/genome_medicine_submission/main_figures')

plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9, 'axes.linewidth':0.6})

def place(ax, path, letter=None):
    if not Path(path).exists():
        ax.text(0.5, 0.5, f'missing: {Path(path).name}', ha='center', va='center',
                transform=ax.transAxes, color='red'); ax.axis('off'); return
    img = mpimg.imread(path)
    ax.imshow(img); ax.axis('off')
    if letter:
        ax.text(-0.02, 1.03, letter, transform=ax.transAxes,
                fontsize=16, fontweight='bold', color='#1F3B5C', va='top')

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0],
                      width_ratios=[1.0, 1.0],
                      hspace=0.08, wspace=0.08)

# A — Study design (top-left, spans left column)
ax = fig.add_subplot(gs[0, :])
place(ax, PV3/'Fig1A_study_design_v2.png', letter='A')

# B — Sankey (bottom-left)
ax = fig.add_subplot(gs[1, 0])
place(ax, PV3/'Fig1A_sankey.png', letter='B')

# C — Sample matrix (bottom-right would be covered by D preview — rearrange):
# Use 3-row layout:  row 1 = A design (wide)  /  row 2 = B sankey + C matrix  /  row 3 = D preview
plt.close(fig)

fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.20, 0.85], hspace=0.12, wspace=0.06)

ax = fig.add_subplot(gs[0, :])
place(ax, PV3/'Fig1A_study_design_v2.png', letter='A')

ax = fig.add_subplot(gs[1, 0])
place(ax, PV3/'Fig1A_sankey.png', letter='B')

ax = fig.add_subplot(gs[1, 1])
place(ax, PV3/'Fig1D_sample_matrix.png', letter='C')

ax = fig.add_subplot(gs[2, :])
place(ax, PV3/'Fig1D_preview_forest.png', letter='D')

# NO suptitle (user rule: no figure title)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig1_cohort.{ext}', dpi=600, bbox_inches='tight')
plt.close(fig)
print('saved new Fig1_cohort')
