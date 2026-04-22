"""
Rebuild Fig 3C / Fig 4D — GENE-level DEG volcano with pathway-category
colouring (Cell cycle / DNA repair / Immune / EMT-Stromal / Other).

User directive (2026-04-22, second pass): Fig 4D (volcano) must be
GENE-level (each dot = gene) because Fig 3B already carries the
pathway-level NES × FDR bubble. The original fig3_v32 style was
"gene-level + pathway-category colouring" — that is the target.

Fix the original's inversion bug: previously coloured every gene in a
pathway list regardless of P-value, so non-significant pathway members
were vivid and significant non-pathway genes were grey. Here we:

  (a) use broad Hallmark gene sets (~200 gene each) so pathway-member
      coverage of the 228 sig genes is actually non-empty;
  (b) gate colouring by significance — only sig genes get category
      colour; ns genes stay grey (standard EnhancedVolcano convention);
  (c) for sig genes that are NOT in any of the four categories,
      colour by DIRECTION (good-up teal / poor-up coral) so the plot
      still reads as a proper volcano.

Data:
  - 05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv (22,838 genes, 228 at P<0.01)
  - Hallmark gene sets fetched via gseapy.get_library(MSigDB_Hallmark_2020)

Category → gene set (union across member Hallmarks):
  Cell cycle  = E2F_TARGETS ∪ G2M_CHECKPOINT ∪ MYC_TARGETS_V1/V2 ∪ MITOTIC_SPINDLE
  DNA repair  = DNA_REPAIR ∪ UV_RESPONSE_UP/DN
  Immune      = IFN α/γ ∪ INFLAMMATORY ∪ COMPLEMENT ∪ IL2/IL6 ∪ ALLOGRAFT ∪ TNFα ∪ COAGULATION
  EMT/Stromal = EMT ∪ MYOGENESIS ∪ APICAL_JUNCTION ∪ ANGIOGENESIS ∪ HEDGEHOG

A gene in multiple categories is assigned by priority
(Cell cycle > DNA repair > Immune > EMT > Other).

Output (overwrites):
  figures/panels_v3/Fig3C_volcano_journal.{pdf,png}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from adjustText import adjust_text
import gseapy

sys.path.insert(0, '/data/data/TNT/analysis/scripts')
from _fig_style import setup_style, save_panel, add_axis_spines  # noqa: E402
setup_style()

GOOD_DEEP = '#0a7d6e'
BAD_DEEP  = '#c53e1f'
INK       = '#0e2a47'

ROOT = Path('/data/data/TNT/analysis')
OUT  = ROOT/'figures/panels_v3'

deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
deg['-log10p'] = -np.log10(deg.pvalue.replace(0, 1e-300))
deg = deg.dropna(subset=['log2FoldChange', '-log10p']).copy()

# ---- fetch Hallmark gene sets (gseapy caches locally) ----
hm_sets = gseapy.get_library(name='MSigDB_Hallmark_2020', organism='Human')

CELL_CYCLE = (set(hm_sets.get('E2F Targets', []))
              | set(hm_sets.get('G2-M Checkpoint', []))
              | set(hm_sets.get('MYC Targets V1', []))
              | set(hm_sets.get('MYC Targets V2', []))
              | set(hm_sets.get('Mitotic Spindle', [])))
DNA_REPAIR = (set(hm_sets.get('DNA Repair', []))
              | set(hm_sets.get('UV Response Up', []))
              | set(hm_sets.get('UV Response Dn', [])))
IMMUNE = (set(hm_sets.get('Interferon Alpha Response', []))
          | set(hm_sets.get('Interferon Gamma Response', []))
          | set(hm_sets.get('Inflammatory Response', []))
          | set(hm_sets.get('Complement', []))
          | set(hm_sets.get('IL-2/STAT5 Signaling', []))
          | set(hm_sets.get('IL-6/JAK/STAT3 Signaling', []))
          | set(hm_sets.get('Allograft Rejection', []))
          | set(hm_sets.get('TNF-alpha Signaling via NF-kB', []))
          | set(hm_sets.get('Coagulation', [])))
EMT_STROMAL = (set(hm_sets.get('Epithelial Mesenchymal Transition', []))
               | set(hm_sets.get('Myogenesis', []))
               | set(hm_sets.get('Apical Junction', []))
               | set(hm_sets.get('Angiogenesis', []))
               | set(hm_sets.get('Hedgehog Signaling', [])))

print(f"Hallmark gene-set sizes:")
print(f"  Cell cycle  = {len(CELL_CYCLE)} genes")
print(f"  DNA repair  = {len(DNA_REPAIR)} genes")
print(f"  Immune      = {len(IMMUNE)} genes")
print(f"  EMT/Stromal = {len(EMT_STROMAL)} genes")

# priority-assigned category (Cell cycle > DNA repair > Immune > EMT > Other)
def cat(g):
    if g in CELL_CYCLE:    return 'Cell cycle'
    if g in DNA_REPAIR:    return 'DNA repair'
    if g in IMMUNE:        return 'Immune'
    if g in EMT_STROMAL:   return 'EMT/Stromal'
    return 'Other'

deg['cat'] = deg.gene.apply(cat)

P_THRESH = 0.01
sig_mask = deg.pvalue < P_THRESH
sig_df = deg[sig_mask].copy()
ns_df  = deg[~sig_mask].copy()

print(f"\nSig genes (P < {P_THRESH}): {len(sig_df)}")
print("By category:")
for c in ['Cell cycle', 'DNA repair', 'Immune', 'EMT/Stromal', 'Other']:
    n = (sig_df['cat'] == c).sum()
    n_up = (((sig_df['cat'] == c) & (sig_df.log2FoldChange > 0))).sum()
    n_dn = (((sig_df['cat'] == c) & (sig_df.log2FoldChange < 0))).sum()
    print(f"  {c:12s}  n={n:3d}  (up-good {n_up:3d}, up-poor {n_dn:3d})")

# ---- plotting ----
# Pathway-category palette matches Fig 5A GROUP_PAL (18_fig5_native_pptx.py
# line 422-428) for cross-figure consistency: plum / amber / cerulean / moss.
CAT_COLOR = {
    'Cell cycle':   '#E39A1E',   # amber  (matches Fig 5A Cell cycle)
    'DNA repair':   '#A24E78',   # plum   (matches Fig 5A DNA repair)
    'Immune':       '#118AB2',   # cerulean (matches Fig 5A Immune)
    'EMT/Stromal':  '#5A8261',   # moss   (matches Fig 5A EMT / hypoxia)
    'Other-good':   GOOD_DEEP,   # sig, no pathway, up-good (project teal)
    'Other-poor':   BAD_DEEP,    # sig, no pathway, up-poor (project coral)
}
NS_COLOR = '#c0c6cf'

fig, ax = plt.subplots(figsize=(9, 7))

# non-significant: hexbin density backdrop + faint dots
ax.hexbin(ns_df.log2FoldChange.clip(-5, 5),
          ns_df['-log10p'].clip(0, 8),
          gridsize=55, mincnt=2,
          cmap=LinearSegmentedColormap.from_list('h', ['#ffffff', '#dde2e8']),
          zorder=1)
ax.scatter(ns_df.log2FoldChange.clip(-5, 5),
           ns_df['-log10p'].clip(0, 8),
           s=8, alpha=0.30, color=NS_COLOR,
           edgecolor='none', zorder=2)

# significant genes: colour by pathway category (priority order), else direction
def marker_color(row):
    if row['cat'] != 'Other':
        return CAT_COLOR[row['cat']]
    return CAT_COLOR['Other-good'] if row.log2FoldChange > 0 else CAT_COLOR['Other-poor']

sig_df['marker_color'] = sig_df.apply(marker_color, axis=1)

# plot "Other" direction-coloured first so pathway dots land on top
for label, color, draw_first in [
        ('Other-good', CAT_COLOR['Other-good'], True),
        ('Other-poor', CAT_COLOR['Other-poor'], True),
        ('EMT/Stromal', CAT_COLOR['EMT/Stromal'], False),
        ('Immune', CAT_COLOR['Immune'], False),
        ('DNA repair', CAT_COLOR['DNA repair'], False),
        ('Cell cycle', CAT_COLOR['Cell cycle'], False)]:
    if label.startswith('Other'):
        dir_sign = 1 if label.endswith('good') else -1
        sub = sig_df[(sig_df['cat'] == 'Other') & (np.sign(sig_df.log2FoldChange) == dir_sign)]
    else:
        sub = sig_df[sig_df['cat'] == label]
    if not len(sub):
        continue
    z = 3 if draw_first else 5
    s = 22 if draw_first else 36
    ax.scatter(sub.log2FoldChange.clip(-5, 5),
               sub['-log10p'].clip(0, 8),
               s=s, alpha=0.85 if draw_first else 0.95,
               color=color, edgecolor='#0e2a47', linewidth=0.35,
               zorder=z)

# Gold emphasis ring removed (user feedback 2026-04-22: too many dots at
# |log2FC| >= 1 cluttered the visual; pathway-category colouring already
# carries the structural emphasis).

# labels: top pathway-coloured genes + top direction-only genes
pathway_sig = sig_df[sig_df['cat'] != 'Other']
to_label_path = (pathway_sig
                 .assign(score=pathway_sig['-log10p'] * pathway_sig.log2FoldChange.abs())
                 .sort_values('score', ascending=False)
                 .head(25))
other_sig = sig_df[sig_df['cat'] == 'Other']
to_label_other = (other_sig
                  .assign(score=other_sig['-log10p'] * other_sig.log2FoldChange.abs())
                  .sort_values('score', ascending=False)
                  .head(10))
to_label = pd.concat([to_label_path, to_label_other], ignore_index=True)

texts = []
for _, r in to_label.iterrows():
    color_ = r['marker_color']
    t = ax.text(np.clip(r.log2FoldChange, -5, 5),
                np.clip(r['-log10p'], 0, 8),
                r.gene, fontsize=8.5,
                color=color_, fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.35),
                expand_points=(1.4, 1.5), force_points=0.55)
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

# legend
legend_handles = []
for c in ['Cell cycle', 'DNA repair', 'Immune', 'EMT/Stromal']:
    n = ((sig_df['cat'] == c)).sum()
    if n == 0:
        continue
    legend_handles.append(
        Line2D([0], [0], marker='o', linestyle='', markersize=8,
               markerfacecolor=CAT_COLOR[c], markeredgecolor='#0e2a47',
               markeredgewidth=0.4,
               label=f'{c}  (n = {n})'))
n_other_up = ((sig_df['cat'] == 'Other') & (sig_df.log2FoldChange > 0)).sum()
n_other_dn = ((sig_df['cat'] == 'Other') & (sig_df.log2FoldChange < 0)).sum()
legend_handles.append(
    Line2D([0], [0], marker='o', linestyle='', markersize=7,
           markerfacecolor=GOOD_DEEP, markeredgecolor='#0e2a47',
           markeredgewidth=0.4,
           label=f'Other, up-good  (n = {n_other_up})'))
legend_handles.append(
    Line2D([0], [0], marker='o', linestyle='', markersize=7,
           markerfacecolor=BAD_DEEP, markeredgecolor='#0e2a47',
           markeredgewidth=0.4,
           label=f'Other, up-poor  (n = {n_other_dn})'))
legend_handles.append(
    Line2D([0], [0], marker='o', linestyle='', markersize=6,
           markerfacecolor=NS_COLOR, markeredgecolor='none',
           label=f'Non-significant (P ≥ {P_THRESH})'))

ax.legend(handles=legend_handles, loc='lower left',
          fontsize=8.5, title='Pathway category / direction',
          title_fontsize=9.5, frameon=True, framealpha=0.95,
          edgecolor=INK)

add_axis_spines(ax)
save_panel(fig, 'Fig3C_volcano_journal', OUT)
print(f'\nwrote {OUT}/Fig3C_volcano_journal.{{pdf,png}}')
