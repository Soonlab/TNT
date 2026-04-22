"""
Rebuild Fig 3C / Fig 4D — DEG volcano with standard volcano-plot colouring.

Fix: prior version (fig3_v32_polished.py:297-315) coloured every gene by
pathway-category membership regardless of P-value, so genes BELOW the
P=0.01 line appeared vivid while significant genes in the "Other"
bucket (most of the top-ranked DEGs) were rendered grey. This inverts
the standard volcano convention (significant = coloured,
non-significant = grey).

Further empirical finding (2026-04-22 recompute): of the 228 genes with
P < 0.01 in the discovery DEG, 0 intersect the curated pathway gene
lists (cell cycle / DNA repair / immune / EMT). Pathway-category
colouring therefore carries no information at this significance
threshold; individual-gene signal in this cohort is carried by
non-canonical markers. We drop the pathway colour scheme and instead
use the standard EnhancedVolcano-style direction colouring:

  - Non-significant (P >= 0.01)     → soft grey hexbin + faint dots
  - Significant + up in good        → teal (GOOD_DEEP)
  - Significant + up in poor        → coral (BAD_DEEP)

Pathway-level biology is carried by Fig 3 (GSEA / ssGSEA) — the volcano
now honestly shows "many nominally-significant non-canonical genes,
but none are pathway-specific at this threshold", which matches the
manuscript's own stance (no FDR-significant DEGs; pathway signal is
aggregate).

Output (overwrites):
  figures/panels_v3/Fig3C_volcano_journal.{pdf,png}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from adjustText import adjust_text

sys.path.insert(0, '/data/data/TNT/analysis/scripts')
from _fig_style import setup_style, save_panel, add_axis_spines  # noqa: E402
setup_style()

GOOD_DEEP = '#0a7d6e'
BAD_DEEP  = '#c53e1f'
INK       = '#0e2a47'

ROOT = Path('/data/data/TNT/analysis')
OUT  = ROOT/'figures/panels_v3'

deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')

# -------------- pathway category gene lists (identical to fig3_v32) --------
CC_GENES  = {'MKI67','TOP2A','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7',
             'PCNA','CCNB1','CCNB2','CCNE1','CCNE2','CCND1','CDK1','CDK2',
             'CDK4','CDK6','AURKA','AURKB','BUB1','BUB1B','CDC20','PLK1',
             'E2F1','E2F2','MYBL2','FOXM1','TTK','TPX2','KIF2C'}
DR_GENES  = {'BRCA1','BRCA2','PALB2','ATM','ATR','CHEK1','CHEK2','MSH2','MSH6',
             'MLH1','PMS2','RAD51','RAD51B','RAD51C','RAD51D','RAD52','RAD54L',
             'XRCC2','XRCC3','FANCA','FANCD2','FANCI','BRIP1','BARD1','NBN',
             'MRE11','EXO1','BLM','WRN','POLQ','LIG4','XRCC4','XRCC6','PRKDC'}
EMT_GENES = {'VIM','CDH2','SNAI1','SNAI2','TWIST1','TWIST2','FN1','MMP2','MMP3',
             'MMP9','MMP14','ZEB1','ZEB2','TGFB1','TGFB2','TGFBR1','COL1A1',
             'COL3A1','COL4A1','FAP','ACTA2','S100A4','POSTN','SPARC'}
IM_GENES  = {'CD8A','CD8B','GZMA','GZMB','GZMK','GZMH','PRF1','GNLY','IFNG',
             'CD3D','CD3E','HLA-A','HLA-B','HLA-C','B2M','NLRC5','TAP1','TAP2',
             'CIITA','HLA-DRA','HLA-DRB1','MS4A1','CD79A','CD79B','BANK1',
             'FOXP3','CTLA4','PDCD1','LAG3','HAVCR2','TIGIT','CD274','CXCL9',
             'CXCL10','CCL19','CCL21','CXCL13'}

def cat(g):
    if g in CC_GENES: return 'Cell cycle'
    if g in DR_GENES: return 'DNA repair'
    if g in IM_GENES: return 'Immune'
    if g in EMT_GENES: return 'EMT/Stromal'
    return 'Other'

deg_p = deg.copy()
deg_p['cat']     = deg_p.gene.apply(cat)
deg_p['-log10p'] = -np.log10(deg_p.pvalue.replace(0, 1e-300))

NS_COLOR = '#c0c6cf'             # soft grey for non-significant
P_THRESH = 0.01
FC_HIGHLIGHT = 1.0               # extra emphasis for |log2FC| >= 1

sig_mask   = deg_p.pvalue < P_THRESH
sig_df     = deg_p[sig_mask].copy()
ns_df      = deg_p[~sig_mask].copy()

up_good = sig_df[sig_df.log2FoldChange > 0]
up_poor = sig_df[sig_df.log2FoldChange < 0]

print(f'DEG volcano: {len(deg_p)} genes total; '
      f'{sig_mask.sum()} above P<{P_THRESH}; '
      f'{(~sig_mask).sum()} below (non-significant)')
print(f'  sig-direction: up-in-good={len(up_good)}, up-in-poor={len(up_poor)}')
print(f'  sig-by-category (empirical): '
      + ', '.join(f"{c}={(sig_df.cat==c).sum()}" for c in
                  ['Cell cycle','DNA repair','Immune','EMT/Stromal','Other']))

fig, ax = plt.subplots(figsize=(9, 7))

# ----- non-significant: hexbin density (soft grey, behind) -----
ax.hexbin(ns_df.log2FoldChange.clip(-5, 5),
          ns_df['-log10p'].clip(0, 8),
          gridsize=55, mincnt=2,
          cmap=LinearSegmentedColormap.from_list('h', ['#ffffff', '#dde2e8']),
          zorder=1)
# any straggler non-sig points outside hexbin density  → faint dots
ax.scatter(ns_df.log2FoldChange.clip(-5, 5),
           ns_df['-log10p'].clip(0, 8),
           s=8, alpha=0.35, color=NS_COLOR,
           edgecolor='none', zorder=2)

# ----- significant: coloured by direction (good-up = teal, poor-up = coral) -----
ax.scatter(up_good.log2FoldChange.clip(-5, 5),
           up_good['-log10p'].clip(0, 8),
           s=30, alpha=0.88, color=GOOD_DEEP,
           edgecolor='#0e2a47', linewidth=0.35,
           label=f'Up in Good (n={len(up_good)})', zorder=4)
ax.scatter(up_poor.log2FoldChange.clip(-5, 5),
           up_poor['-log10p'].clip(0, 8),
           s=30, alpha=0.88, color=BAD_DEEP,
           edgecolor='#0e2a47', linewidth=0.35,
           label=f'Up in Poor (n={len(up_poor)})', zorder=4)

# extra emphasis ring for |log2FC| >= 1
big_fc = sig_df[sig_df.log2FoldChange.abs() >= FC_HIGHLIGHT]
ax.scatter(big_fc.log2FoldChange.clip(-5, 5),
           big_fc['-log10p'].clip(0, 8),
           s=80, facecolors='none',
           edgecolor='#D4A300', linewidth=1.1, zorder=5)

# ----- top-labeled genes (most extreme by -log10p × |log2FC| among sig) -----
to_label = (sig_df.assign(score=sig_df['-log10p'] * sig_df.log2FoldChange.abs())
                   .sort_values('score', ascending=False)
                   .head(35))
texts = []
for _, r in to_label.iterrows():
    color_ = GOOD_DEEP if r.log2FoldChange > 0 else BAD_DEEP
    t = ax.text(np.clip(r.log2FoldChange, -5, 5),
                np.clip(r['-log10p'], 0, 8),
                r.gene, fontsize=8.5,
                color=color_,
                fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.3),
                expand_points=(1.4, 1.5), force_points=0.5)
except Exception:
    pass

# axis decorations
ax.axvline(0, color=INK, lw=0.6, alpha=0.5)
ax.axhline(-np.log10(P_THRESH), color='#5a6772', lw=0.7, ls='--', alpha=0.7)
ax.text(-4.95, -np.log10(P_THRESH) + 0.08,
        f'P = {P_THRESH}', fontsize=8.5, color='#5a6772', ha='left')

ax.set_xlim(-5.2, 5.2)
ax.set_ylim(0, 8.5)

# direction arrows
y_arrow = 7.6
ax.annotate('', xy=(4.7, y_arrow), xytext=(0.7, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=GOOD_DEEP,
                            lw=2.5, mutation_scale=18))
ax.text(2.7, y_arrow + 0.18, 'Up in Good responders',
        ha='center', va='bottom', fontsize=10,
        color=GOOD_DEEP, fontweight='bold')
ax.annotate('', xy=(-4.7, y_arrow), xytext=(-0.7, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=BAD_DEEP,
                            lw=2.5, mutation_scale=18))
ax.text(-2.7, y_arrow + 0.18, 'Up in Poor responders',
        ha='center', va='bottom', fontsize=10,
        color=BAD_DEEP, fontweight='bold')

ax.set_xlabel('log2 fold change (good vs poor)',
              fontsize=11, fontweight='bold', color=INK)
ax.set_ylabel('−log10(p-value)',
              fontsize=11, fontweight='bold', color=INK)

# legend lower-left: direction colours + |log2FC|>=1 ring + non-sig chip
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker='o', linestyle='', markersize=8,
           markerfacecolor=GOOD_DEEP, markeredgecolor='#0e2a47',
           markeredgewidth=0.4,
           label=f'Up in Good, P < {P_THRESH}  (n = {len(up_good)})'),
    Line2D([0], [0], marker='o', linestyle='', markersize=8,
           markerfacecolor=BAD_DEEP, markeredgecolor='#0e2a47',
           markeredgewidth=0.4,
           label=f'Up in Poor, P < {P_THRESH}  (n = {len(up_poor)})'),
    Line2D([0], [0], marker='o', linestyle='', markersize=11,
           markerfacecolor='none', markeredgecolor='#D4A300',
           markeredgewidth=1.2,
           label=f'|log2 FC| ≥ {FC_HIGHLIGHT:.0f}'),
    Line2D([0], [0], marker='o', linestyle='', markersize=7,
           markerfacecolor=NS_COLOR, markeredgecolor='none',
           label=f'Non-significant (P ≥ {P_THRESH})'),
]
ax.legend(handles=legend_handles, loc='lower left',
          fontsize=9, title='Significance  /  direction',
          title_fontsize=10, frameon=True, framealpha=0.95,
          edgecolor=INK)

add_axis_spines(ax)
save_panel(fig, 'Fig3C_volcano_journal', OUT)
print(f'wrote {OUT}/Fig3C_volcano_journal.{{pdf,png}}')
