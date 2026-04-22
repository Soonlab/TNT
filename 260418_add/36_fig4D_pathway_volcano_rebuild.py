"""
Rebuild Fig 3C / Fig 4D — **pathway-level volcano** (Hallmark GSEA).

User directive (2026-04-22): the volcano should be pathway-level, not
single-gene-level. Pathway categories (cell cycle / DNA repair /
immune / EMT-Stromal / Other) are the intended colouring dimension —
matching the Thread 1 / Thread 2 narrative of the manuscript (cell-
cycle + DNA-repair dominant in good responders, EMT/stromal dominant
in poor; immune axis external-validated).

Data: `05_rna_deg_gsea/GSEA_Hallmark_pre.tsv` (fgsea Hallmark,
N = 50 pathways; 10 at padj < 0.05; 22 at padj < 0.25).

X-axis  : NES (normalized enrichment score)
Y-axis  : -log10(padj)
Colour  : 5 curated categories derived from Hallmark naming convention
Thresholds: padj < 0.05 solid dashed line, padj < 0.25 (Broad default)
Emphasis: |NES| >= 2 gold ring (large effect); top ~15 labelled by
          -log10(padj) * |NES|.

Output (overwrites):
  figures/panels_v3/Fig3C_volcano_journal.{pdf,png}
Also rebuilt Fig3_RNA_signatures composite via the inline runner
inside 35_rebuild_composite_fig3.py (called downstream).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/data/data/TNT/analysis/scripts')
from _fig_style import setup_style, save_panel, add_axis_spines  # noqa: E402
setup_style()

GOOD_DEEP = '#0a7d6e'
BAD_DEEP  = '#c53e1f'
INK       = '#0e2a47'

ROOT = Path('/data/data/TNT/analysis')
OUT  = ROOT/'figures/panels_v3'

hm = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
re = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')

# ---- Hallmark pathway → 5 categories ----
HM_CELL_CYCLE = {
    'HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT',
    'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2',
    'HALLMARK_MITOTIC_SPINDLE',
}
HM_DNA_REPAIR = {'HALLMARK_DNA_REPAIR'}
HM_IMMUNE = {
    'HALLMARK_INTERFERON_ALPHA_RESPONSE',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE',
    'HALLMARK_INFLAMMATORY_RESPONSE',
    'HALLMARK_COMPLEMENT',
    'HALLMARK_IL2_STAT5_SIGNALING',
    'HALLMARK_IL6_JAK_STAT3_SIGNALING',
    'HALLMARK_ALLOGRAFT_REJECTION',
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
    'HALLMARK_COAGULATION',
}
HM_EMT_STROMAL = {
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
    'HALLMARK_MYOGENESIS',
    'HALLMARK_APICAL_JUNCTION',
    'HALLMARK_APICAL_SURFACE',
    'HALLMARK_ANGIOGENESIS',
    'HALLMARK_HEDGEHOG_SIGNALING',
}

# ---- Curated Reactome additions (augment sparse Hallmark rows) ----
RE_CELL_CYCLE = {
    'REACTOME_MITOTIC_PROMETAPHASE',
    'REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE',
    'REACTOME_MITOTIC_SPINDLE_CHECKPOINT',
    'REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION',
    'REACTOME_MITOTIC_G2_G2_M_PHASES',
}
RE_DNA_REPAIR = {
    'REACTOME_HOMOLOGY_DIRECTED_REPAIR',
    'REACTOME_MISMATCH_REPAIR',
}
RE_IMMUNE = {
    'REACTOME_CD22_MEDIATED_BCR_REGULATION',
    'REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR',
    'REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS',
    'REACTOME_TCR_SIGNALING',
    'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',
}
RE_EMT_STROMAL = {
    'REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES',
    'REACTOME_ELASTIC_FIBRE_FORMATION',
    'REACTOME_ECM_PROTEOGLYCANS',
    'REACTOME_COLLAGEN_DEGRADATION',
    'REACTOME_COLLAGEN_FORMATION',
    'REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS',
}

def cat(name):
    if name in HM_CELL_CYCLE | RE_CELL_CYCLE: return 'Cell cycle'
    if name in HM_DNA_REPAIR | RE_DNA_REPAIR: return 'DNA repair'
    if name in HM_IMMUNE    | RE_IMMUNE:     return 'Immune'
    if name in HM_EMT_STROMAL | RE_EMT_STROMAL: return 'EMT/Stromal'
    return 'Other'

# assemble combined frame: all Hallmark + selected Reactome
re_keep = RE_CELL_CYCLE | RE_DNA_REPAIR | RE_IMMUNE | RE_EMT_STROMAL
re_sel = re[re.pathway.isin(re_keep)].copy()
# add a source column for readability
hm['source'] = 'Hallmark'
re_sel['source'] = 'Reactome'
hm = pd.concat([hm, re_sel], ignore_index=True)
hm['cat'] = hm.pathway.apply(cat)
hm['-log10padj'] = -np.log10(hm.padj.replace(0, 1e-300))
hm['short'] = (hm.pathway
               .str.replace('HALLMARK_', '', regex=False)
               .str.replace('REACTOME_', '', regex=False)
               .str.replace('_', ' ')
               .str.title())
# abbreviate long Reactome names
hm['short'] = hm['short'].replace({
    'Antigen Activates B Cell Receptor Bcr Leading To Generation Of Second Messengers': 'BCR → 2nd messengers',
    'Assembly Of Collagen Fibrils And Other Multimeric Structures': 'Collagen assembly',
    'Signaling By The B Cell Receptor Bcr': 'BCR signaling',
    'Cd22 Mediated Bcr Regulation': 'CD22 BCR reg',
    'Homology Directed Repair': 'HDR',
    'Mismatch Repair': 'MMR',
    'Mitotic Prometaphase': 'Mitotic prometaphase',
    'Mitotic Metaphase And Anaphase': 'Metaphase-anaphase',
    'Mitotic Spindle Checkpoint': 'Spindle checkpoint',
    'Mitotic G1 Phase And G1 S Transition': 'G1-S transition',
    'Mitotic G2 G2 M Phases': 'G2/M phase',
    'Interferon Alpha Beta Signaling': 'IFN α/β signaling',
    'Ecm Proteoglycans': 'ECM proteoglycans',
    'Elastic Fibre Formation': 'Elastic fibre',
    'Collagen Degradation': 'Collagen degradation',
    'Collagen Formation': 'Collagen formation',
    'Integrin Cell Surface Interactions': 'Integrin–cell surface',
    'Tcr Signaling': 'TCR signaling',
})

# simple rendering palette per category
PAL = {
    'Cell cycle':   '#057a64',
    'DNA repair':   '#00567d',
    'Immune':       '#c11456',
    'EMT/Stromal':  '#b03219',
    'Other':        '#c0c6cf',
}
NES_EMPHASIS = 2.0
PADJ_STRICT  = 0.05
PADJ_LENIENT = 0.25

n_sig = (hm.padj < PADJ_STRICT).sum()
n_lenient = (hm.padj < PADJ_LENIENT).sum()
print(f'Hallmark: N={len(hm)} pathways, padj<{PADJ_STRICT}: {n_sig}, padj<{PADJ_LENIENT}: {n_lenient}')
print(f'By category: '
      + ', '.join(f'{c}={(hm.cat==c).sum()}' for c in PAL))

# ---- plot ----
fig, ax = plt.subplots(figsize=(9, 7))

# NES reference lines
ax.axvline(0, color=INK, lw=0.6, alpha=0.5)
ax.axhline(-np.log10(PADJ_STRICT), color='#5a6772', lw=0.7, ls='--', alpha=0.8)
ax.axhline(-np.log10(PADJ_LENIENT), color='#a0a6ad', lw=0.5, ls=':', alpha=0.8)
ax.text(-3.1, -np.log10(PADJ_STRICT) + 0.15,
        f'BH padj = {PADJ_STRICT}', fontsize=8.5, color='#5a6772')
ax.text(-3.1, -np.log10(PADJ_LENIENT) + 0.10,
        f'BH padj = {PADJ_LENIENT}', fontsize=8, color='#8a8f98', style='italic')

# scatter by category (non-significant categories under, significant over)
cat_order = ['Other', 'EMT/Stromal', 'Immune', 'DNA repair', 'Cell cycle']
for c in cat_order:
    sub = hm[hm.cat == c]
    is_sig = sub.padj < PADJ_LENIENT
    # non-sig → faded
    ns = sub[~is_sig]
    if len(ns):
        ax.scatter(ns.NES, ns['-log10padj'], s=45,
                   color=PAL[c], alpha=0.28, edgecolor='#dde2e8',
                   linewidth=0.3, zorder=2)
    # sig → vivid
    sg = sub[is_sig]
    if len(sg):
        ax.scatter(sg.NES, sg['-log10padj'], s=110,
                   color=PAL[c], alpha=0.92, edgecolor=INK,
                   linewidth=0.55, label=c if c != 'Other' else None,
                   zorder=4)
# explicit "Other (sig)" entry in legend (sometimes absent)
other_sig = hm[(hm.cat == 'Other') & (hm.padj < PADJ_LENIENT)]
if len(other_sig):
    ax.scatter([], [], s=90, color=PAL['Other'], edgecolor=INK,
               linewidth=0.55, label='Other')

# |NES| >= 2 emphasis ring (large-effect)
big = hm[hm.NES.abs() >= NES_EMPHASIS]
ax.scatter(big.NES, big['-log10padj'], s=250, facecolors='none',
           edgecolor='#D4A300', linewidth=1.3, zorder=5,
           label=f'|NES| ≥ {NES_EMPHASIS:.0f}')

# pathway labels (top by -log10padj * |NES|, only among padj < lenient)
to_label = (hm[hm.padj < PADJ_LENIENT]
            .assign(score=lambda d: d['-log10padj'] * d.NES.abs())
            .sort_values('score', ascending=False)
            .head(15))
from adjustText import adjust_text
texts = []
for _, r in to_label.iterrows():
    t = ax.text(r.NES, r['-log10padj'], r['short'],
                fontsize=9, color=PAL[r['cat']], fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.4),
                expand_points=(1.5, 1.6), force_points=0.6)
except Exception:
    pass

# direction arrows
xlo, xhi = -3.2, 3.2
nlp_max = max(26, float(hm['-log10padj'].max()) * 1.03)
y_arrow = nlp_max * 0.95
ax.annotate('', xy=(xhi * 0.85, y_arrow), xytext=(0.6, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=GOOD_DEEP,
                            lw=2.5, mutation_scale=18))
ax.text(xhi * 0.45, y_arrow + nlp_max * 0.02, 'Enriched in Good',
        ha='center', va='bottom', fontsize=10,
        color=GOOD_DEEP, fontweight='bold')
ax.annotate('', xy=(-xhi * 0.85, y_arrow), xytext=(-0.6, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=BAD_DEEP,
                            lw=2.5, mutation_scale=18))
ax.text(-xhi * 0.45, y_arrow + nlp_max * 0.02, 'Enriched in Poor',
        ha='center', va='bottom', fontsize=10,
        color=BAD_DEEP, fontweight='bold')

ax.set_xlim(xlo, xhi)
ax.set_ylim(0, nlp_max)
ax.set_xlabel('NES (normalised enrichment score)',
              fontsize=11, fontweight='bold', color=INK)
ax.set_ylabel('−log₁₀(BH padj)',
              fontsize=11, fontweight='bold', color=INK)

# legend
ax.legend(loc='center right', fontsize=9, title='Pathway category',
          title_fontsize=10, frameon=True, framealpha=0.95,
          edgecolor=INK)

add_axis_spines(ax)
save_panel(fig, 'Fig3C_volcano_journal', OUT)
print(f'wrote {OUT}/Fig3C_volcano_journal.{{pdf,png}} '
      f'({len(hm)} Hallmark pathways)')
