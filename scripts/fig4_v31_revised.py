"""
Fig 4 v3.1 — revisions per user feedback:
  4A: redesigned with recent papers (Mariathasan Nature 2018, Bagaev Cancer Cell 2021,
      Krishna Cancer Cell 2024, Pan Nat Med 2024) — modern compact running ES with
      leading-edge density rather than barcode + ranked metric
  4B: separate size legend (gene set count) from category legend
  4C: clarify Hallmark / Reactome — use color-coded source bar instead of ◇/○ marker
  4D: edges visible (thicker, darker, lower threshold)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from adjustText import adjust_text
import gseapy as gp
import warnings; warnings.filterwarnings('ignore')

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}

HALLMARK_CAT = {
    'E2F_TARGETS':'Proliferation','G2M_CHECKPOINT':'Proliferation','MYC_TARGETS_V1':'Proliferation',
    'MYC_TARGETS_V2':'Proliferation','MITOTIC_SPINDLE':'Proliferation','MTORC1_SIGNALING':'Signaling',
    'DNA_REPAIR':'DNA repair','UV_RESPONSE_DN':'DNA repair','UNFOLDED_PROTEIN_RESPONSE':'Stress',
    'EPITHELIAL_MESENCHYMAL_TRANSITION':'Stromal/EMT','MYOGENESIS':'Stromal/EMT','APICAL_JUNCTION':'Stromal/EMT',
    'TGF_BETA_SIGNALING':'Stromal/EMT','HEDGEHOG_SIGNALING':'Signaling','WNT_BETA_CATENIN_SIGNALING':'Signaling',
    'NOTCH_SIGNALING':'Signaling','PI3K_AKT_MTOR_SIGNALING':'Signaling','KRAS_SIGNALING_UP':'Signaling',
    'KRAS_SIGNALING_DN':'Signaling','IL2_STAT5_SIGNALING':'Immune','IL6_JAK_STAT3_SIGNALING':'Immune',
    'INTERFERON_GAMMA_RESPONSE':'Immune','INTERFERON_ALPHA_RESPONSE':'Immune','INFLAMMATORY_RESPONSE':'Immune',
    'TNFA_SIGNALING_VIA_NFKB':'Immune','COMPLEMENT':'Immune','ALLOGRAFT_REJECTION':'Immune','COAGULATION':'Immune',
    'OXIDATIVE_PHOSPHORYLATION':'Metabolism','GLYCOLYSIS':'Metabolism','HYPOXIA':'Metabolism',
    'CHOLESTEROL_HOMEOSTASIS':'Metabolism','FATTY_ACID_METABOLISM':'Metabolism','XENOBIOTIC_METABOLISM':'Metabolism',
    'ADIPOGENESIS':'Metabolism','PEROXISOME':'Metabolism','HEME_METABOLISM':'Metabolism',
    'BILE_ACID_METABOLISM':'Metabolism','ESTROGEN_RESPONSE_EARLY':'Signaling','ESTROGEN_RESPONSE_LATE':'Signaling',
    'ANDROGEN_RESPONSE':'Signaling','P53_PATHWAY':'TP53','APOPTOSIS':'Apoptosis','PROTEIN_SECRETION':'Other',
    'UV_RESPONSE_UP':'Other','ANGIOGENESIS':'Other','PANCREAS_BETA_CELLS':'Other','SPERMATOGENESIS':'Other',
    'REACTIVE_OXYGEN_SPECIES_PATHWAY':'Stress',
}
GSEA_CAT_DEEP = {
    'Proliferation':'#057a64','DNA repair':'#00567d','Immune':'#c11456','Stromal/EMT':'#b03219',
    'Metabolism':'#d4a300','Signaling':'#7a3aad','Apoptosis':'#5a3582','Stress':'#0099b8',
    'TP53':'#ad2831','Other':'#5a6772',
}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
gsea_h = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
gsea_r = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
tpm = pd.read_csv(ROOT/'06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)

deg_p = deg.dropna(subset=['log2FoldChange','gene']).copy()
deg_p = deg_p.sort_values('log2FoldChange', ascending=False).reset_index(drop=True)
ranks = pd.Series(deg_p.log2FoldChange.values, index=deg_p.gene.values)

hm_sets = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')

def gseapy_pathway_key(name):
    if name.startswith('HALLMARK_'):
        return name.replace('HALLMARK_','').replace('_',' ').title()
    return name

def running_es(ranks, gene_set):
    n = len(ranks)
    in_set = ranks.index.isin(gene_set).astype(int)
    n_hits = in_set.sum()
    if n_hits == 0: return None, None
    weights = np.abs(ranks.values)
    sum_w_in = (weights * in_set).sum()
    if sum_w_in == 0: return None, None
    increment = (weights * in_set) / sum_w_in
    decrement = (1 - in_set) / (n - n_hits)
    cumsum = np.cumsum(increment - decrement)
    return cumsum, in_set

# ============================================================
# 4A — modern running-ES (Mariathasan Nature 2018 + Pan Nat Med 2024 + Krishna Cancer Cell 2024 motif)
# - 4 pathways in 2x2 grid
# - Each panel: smooth ES area + leading-edge gene density rug below
# - Sample-level expression heatmap (good vs bad) of leading-edge genes
# ============================================================
top_paths = [
    ('HALLMARK_E2F_TARGETS', 'E2F targets'),
    ('HALLMARK_G2M_CHECKPOINT', 'G2M checkpoint'),
    ('HALLMARK_MYC_TARGETS_V1', 'MYC targets V1'),
    ('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT (mesenchymal)'),
]

# Sample-level pre data
pre_samples = rna_inv[rna_inv.timepoint=='pre'].sample_id.tolist()
pre_resp = {s: rna_inv[rna_inv.sample_id==s].response_bin.iloc[0] for s in pre_samples}
sample_order_pre = sorted(pre_samples, key=lambda s: (pre_resp[s]!='good', s))

fig = plt.figure(figsize=(15, 10))
outer_gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20)

for idx, (pw, label) in enumerate(top_paths):
    row, col = divmod(idx, 2)
    inner = outer_gs[row, col].subgridspec(3, 1, height_ratios=[3.0, 0.35, 1.4], hspace=0.10)
    ax_es = fig.add_subplot(inner[0])
    ax_rug = fig.add_subplot(inner[1])
    ax_hm = fig.add_subplot(inner[2])

    key = gseapy_pathway_key(pw)
    gset = set(hm_sets.get(key, []))
    if not gset:
        # Try alternative key formats
        for k in hm_sets:
            if k.lower().replace(' ','').replace('-','') == pw.replace('HALLMARK_','').replace('_','').lower():
                gset = set(hm_sets[k]); break
    es, in_set = running_es(ranks, gset)
    if es is None:
        ax_es.text(0.5, 0.5, f'No data for\n{label}', transform=ax_es.transAxes, ha='center', va='center')
        ax_rug.set_visible(False); ax_hm.set_visible(False)
        continue
    row_g = gsea_h[gsea_h.pathway==pw]
    nes_val = row_g.NES.iloc[0] if len(row_g) else 0
    pval = row_g.pval.iloc[0] if len(row_g) else 1
    color = GOOD_DEEP if nes_val>0 else BAD_DEEP

    # ES curve with smooth gradient fill (modern look)
    x = np.arange(len(es))
    ax_es.plot(x, es, color=color, lw=2.6)
    ax_es.fill_between(x, 0, es, color=color, alpha=0.25)
    ax_es.axhline(0, color='#0e2a47', lw=0.6, alpha=0.5)
    # Peak
    if nes_val>0: peak = np.argmax(es); peak_y = es[peak]
    else: peak = np.argmin(es); peak_y = es[peak]
    ax_es.scatter([peak],[peak_y], color=color, s=110, edgecolor='white', linewidth=1.8, zorder=6)
    ax_es.set_ylabel('ES', fontsize=10, color='#0e2a47', fontweight='bold')
    ax_es.set_xlim(0, len(es)); ax_es.set_xticks([])
    add_axis_spines(ax_es)
    # Pathway label box
    ax_es.text(0.02, 0.96, label, transform=ax_es.transAxes, fontsize=11, fontweight='bold',
               color='#0e2a47', va='top', ha='left',
               bbox=dict(facecolor='white', edgecolor=color, alpha=0.92, boxstyle='round,pad=0.4'))
    ax_es.text(0.98, 0.04 if nes_val>0 else 0.96, f'NES = {nes_val:+.2f}\np = {pval:.2g}',
               transform=ax_es.transAxes, fontsize=10, color=color, fontweight='bold',
               ha='right', va='bottom' if nes_val>0 else 'top')

    # Leading-edge density rug
    hit_pos = np.where(in_set==1)[0]
    for hp in hit_pos:
        ax_rug.axvline(hp, color=color, lw=0.65, alpha=0.85)
    ax_rug.set_xlim(0, len(es)); ax_rug.set_yticks([]); ax_rug.set_xticks([])
    ax_rug.text(-0.005, 0.5, 'Leading\nedge', transform=ax_rug.transAxes, ha='right', va='center',
                fontsize=8.5, color='#0e2a47')
    for s in ['top','right','left','bottom']: ax_rug.spines[s].set_visible(False)

    # Sample-level expression heatmap of leading-edge genes
    if nes_val>0:
        leading = [g for i,g in enumerate(ranks.index) if in_set[i]==1 and i<=peak]
    else:
        leading = [g for i,g in enumerate(ranks.index) if in_set[i]==1 and i>=peak]
    leading = [g for g in leading if g in log_tpm.index][:20]  # cap at 20 for readability
    if leading:
        sub_mat = log_tpm.loc[leading, sample_order_pre]
        sub_z = sub_mat.sub(sub_mat.mean(axis=1), axis=0).div(sub_mat.std(axis=1), axis=0)
        ax_hm.imshow(sub_z.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto',
                     interpolation='nearest', extent=[0, len(sample_order_pre), 0, len(leading)])
        ax_hm.set_yticks([]); ax_hm.set_xticks([])
        # Response strip below heatmap
        for j, s in enumerate(sample_order_pre):
            ax_hm.add_patch(Rectangle((j, -0.45), 1, 0.4, color=PAL[pre_resp[s]],
                                       clip_on=False, transform=ax_hm.transData))
        ax_hm.set_ylim(-0.5, len(leading))
        ax_hm.text(-0.005, 0.5, f'Leading-edge\ngenes (n={len(leading)})',
                   transform=ax_hm.transAxes, ha='right', va='center', fontsize=8.5, color='#0e2a47')
        for s in ['top','right','left','bottom']: ax_hm.spines[s].set_visible(False)
        # X-axis label
        ax_hm.set_xlabel('Pre-treatment samples (good ← → poor)', fontsize=9, color='#0e2a47')

save_panel(fig, 'Fig4A_running_ES', OUT)

# ============================================================
# 4B — bubble plot with separate size legend
# ============================================================
fig, ax = plt.subplots(figsize=(11, 7.5))
h = gsea_h.copy()
h['name'] = h.pathway.str.replace('HALLMARK_','')
h['category'] = h.name.map(HALLMARK_CAT).fillna('Other')
h['neglogp'] = -np.log10(h.pval.replace(0, 1e-300))

cat_handles = []
for cat, color in GSEA_CAT_DEEP.items():
    sub = h[h.category==cat]
    if len(sub)==0: continue
    sizes = sub['size'] * 1.5
    ax.scatter(sub.NES, sub.neglogp, s=sizes, c=color, alpha=0.85,
               edgecolor='#0e2a47', linewidth=0.8, zorder=3)
    cat_handles.append(mpatches.Patch(color=color, label=cat))

texts = []
sig_h = h[h.pval < 1e-3]
for _, r in sig_h.iterrows():
    name_short = r['name'].replace('_',' ').title()[:30]
    t = ax.text(r.NES, r.neglogp, name_short, fontsize=8.5, color='#0e2a47', fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.4),
                expand_points=(1.5, 1.6), force_points=0.6)
except: pass

ax.axvline(0, color='#0e2a47', lw=0.7, alpha=0.5)
ax.axhline(-np.log10(0.05), color='#5a6772', lw=0.6, ls='--', alpha=0.5)
ax.text(ax.get_xlim()[1]*0.97, -np.log10(0.05)+0.5, 'p = 0.05', fontsize=8.5, color='#5a6772', ha='right')
ax.set_xlabel('Normalized enrichment score (NES)\n← UP in poor          UP in good →',
              fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('−log10(p-value)', fontsize=11, fontweight='bold', color='#0e2a47')

# Category legend
leg1 = ax.legend(handles=cat_handles, loc='upper left', fontsize=9, title='Category',
                 title_fontsize=10, ncol=1, bbox_to_anchor=(1.01, 1.0), frameon=False)
ax.add_artist(leg1)

# Size legend (separate, scatter handles with actual sizes)
size_examples = [25, 100, 200]
size_handles = []
for sz in size_examples:
    h_dot = ax.scatter([], [], s=sz*1.5, color='#5a6772', edgecolor='#0e2a47', linewidth=0.8, alpha=0.7,
                        label=f'{sz} genes')
    size_handles.append(h_dot)
ax.legend(handles=size_handles, loc='lower left', fontsize=9, title='Gene set size',
          title_fontsize=10, ncol=1, bbox_to_anchor=(1.01, 0.0), frameon=False,
          scatterpoints=1, labelspacing=1.8, handletextpad=2.0)

add_axis_spines(ax)
save_panel(fig, 'Fig4B_Hallmark_bubble', OUT)

# ============================================================
# 4C — top dotplot: source as colored bar (clearer than ◇/○)
# ============================================================
h_top = gsea_h.copy()
h_top['name'] = h_top.pathway.str.replace('HALLMARK_','')
h_top['source'] = 'Hallmark'
h_top['category'] = h_top.name.map(HALLMARK_CAT).fillna('Other')

r_top = gsea_r.copy()
r_top['name'] = r_top.pathway.str.replace('REACTOME_','').str.replace('_',' ').str.title()
r_top['source'] = 'Reactome'
def reactome_cat(n):
    n = n.lower()
    if 'cycle' in n or 'mitotic' in n or 'm phase' in n or 's phase' in n: return 'Proliferation'
    if 'repair' in n or 'recombination' in n or 'replication' in n or 'damage' in n: return 'DNA repair'
    if 'immune' in n or 'antigen' in n or 'mhc' in n or 'interferon' in n or 'inflammatory' in n: return 'Immune'
    if 'extracellular matrix' in n or 'ecm' in n or 'collagen' in n: return 'Stromal/EMT'
    if 'metabol' in n or 'glycolysis' in n or 'oxidative' in n: return 'Metabolism'
    if 'signal' in n: return 'Signaling'
    if 'apopt' in n: return 'Apoptosis'
    return 'Other'
r_top['category'] = r_top.name.apply(reactome_cat)

combined = pd.concat([h_top.assign(label=h_top.name.str.replace('_',' ').str.title()),
                      r_top.assign(label=r_top.name)], ignore_index=True)
combined = combined.sort_values('pval').head(25).reset_index(drop=True)
combined = combined.sort_values('NES').reset_index(drop=True)

# Use grid: source-color bar on left, then labels, then dot
fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(1, 4, width_ratios=[0.15, 4.5, 6.0, 0.6], wspace=0.02)

# Source color bar (left)
ax_src = fig.add_subplot(gs[0])
SRC_COLOR = {'Hallmark':'#264653', 'Reactome':'#e9c46a'}
src_arr = np.array([[to_rgb(SRC_COLOR[s]) for s in combined.source]])
ax_src.imshow(src_arr.transpose(1,0,2), aspect='auto', interpolation='nearest',
              extent=[0,1,0,len(combined)])
ax_src.set_xticks([]); ax_src.set_yticks([])
ax_src.invert_yaxis()
for s in ['top','right','left','bottom']: ax_src.spines[s].set_visible(False)

# Label column
ax_lbl = fig.add_subplot(gs[1])
ax_lbl.axis('off')
ax_lbl.set_xlim(0,1); ax_lbl.set_ylim(0, len(combined))
for i, r in combined.iterrows():
    ax_lbl.text(0.99, len(combined)-1-i+0.5, r.label[:55], ha='right', va='center',
                fontsize=9.5, color='#0e2a47')

# Dotplot
ax_dot = fig.add_subplot(gs[2])
y = np.arange(len(combined))
sizes = combined['size'] * 1.5
colors = [GSEA_CAT_DEEP[c] for c in combined.category]
ax_dot.axvline(0, color='#0e2a47', lw=0.9)
for i, r in combined.iterrows():
    color = GSEA_CAT_DEEP[r.category]
    ax_dot.plot([0, r.NES], [len(combined)-1-i, len(combined)-1-i], color=color, lw=1.2, alpha=0.45)
ax_dot.scatter(combined.NES, len(combined)-1-y, s=sizes, c=colors, alpha=0.92,
               edgecolor='#0e2a47', linewidth=0.7, zorder=3)
ax_dot.set_yticks([]); ax_dot.set_ylim(-0.5, len(combined)-0.5)
ax_dot.set_xlabel('Normalized enrichment score (NES)', fontsize=11, fontweight='bold', color='#0e2a47')
add_axis_spines(ax_dot, sides=('bottom',))
ax_dot.spines['left'].set_visible(False)

# Legend column (right)
ax_leg = fig.add_subplot(gs[3])
ax_leg.axis('off')
# Source legend
ax_leg.text(0.05, 0.99, 'Source', transform=ax_leg.transAxes, fontsize=10, fontweight='bold',
            color='#0e2a47', va='top', ha='left')
for i, (k, c) in enumerate(SRC_COLOR.items()):
    y_pos = 0.93 - i*0.05
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.018), 0.25, 0.035, color=c, transform=ax_leg.transAxes))
    ax_leg.text(0.34, y_pos, k, transform=ax_leg.transAxes, fontsize=9, va='center', color='#0e2a47')

# Category legend
ax_leg.text(0.05, 0.77, 'Category', transform=ax_leg.transAxes, fontsize=10, fontweight='bold',
            color='#0e2a47', va='top', ha='left')
cats_used = list(combined.category.unique())
for i, c in enumerate(cats_used):
    y_pos = 0.71 - i*0.055
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.018), 0.25, 0.035, color=GSEA_CAT_DEEP[c],
                                transform=ax_leg.transAxes))
    ax_leg.text(0.34, y_pos, c, transform=ax_leg.transAxes, fontsize=9, va='center', color='#0e2a47')

save_panel(fig, 'Fig4C_pathway_dotplot', OUT)

# ============================================================
# 4D — enrichment map with VISIBLE edges
# ============================================================
import networkx as nx

h_sel = gsea_h[gsea_h.pval < 0.05].copy()
h_sel['name'] = h_sel.pathway.str.replace('HALLMARK_','')
h_sel['category'] = h_sel.name.map(HALLMARK_CAT).fillna('Other')

G = nx.Graph()
for _, r in h_sel.iterrows():
    key = gseapy_pathway_key(r.pathway)
    gset = set(hm_sets.get(key, []))
    G.add_node(r['name'], nes=r.NES, pval=r.pval, size=len(gset),
               category=r.category, gene_set=gset)

# Edges with lower threshold for visibility
nodes = list(G.nodes())
JACCARD_THR = 0.06  # lowered from 0.10
for i, n1 in enumerate(nodes):
    for n2 in nodes[i+1:]:
        s1 = G.nodes[n1]['gene_set']; s2 = G.nodes[n2]['gene_set']
        if not s1 or not s2: continue
        jac = len(s1 & s2) / len(s1 | s2)
        if jac >= JACCARD_THR:
            G.add_edge(n1, n2, weight=jac)

print(f'  4D: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

fig, ax = plt.subplots(figsize=(12, 9))
pos = nx.spring_layout(G, seed=42, k=2.0, iterations=300, weight='weight')

# Draw EDGES first (visible darker color, thicker)
edge_widths = [G[u][v]['weight'] * 18 for u, v in G.edges()]  # scaled up
nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                        edge_color='#5a6772', alpha=0.55, style='solid')

# Draw NODES on top
node_sizes = [max(250, G.nodes[n]['size'] * 6) for n in G.nodes()]
node_colors = [GSEA_CAT_DEEP[G.nodes[n]['category']] for n in G.nodes()]
node_edgecolors = [GOOD_DEEP if G.nodes[n]['nes']>0 else BAD_DEEP for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
                        edgecolors=node_edgecolors, linewidths=2.5, alpha=0.95)

# Labels
texts = []
for n in G.nodes():
    label = n.replace('_',' ').title()[:25]
    t = ax.text(pos[n][0], pos[n][1], label, fontsize=8, color='#0e2a47', fontweight='bold',
                ha='center', va='center')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.3),
                expand_points=(1.3, 1.4), force_points=0.4)
except: pass

# Legend: category + NES sign + edge meaning
cat_handles = [mpatches.Patch(color=c, label=k) for k, c in GSEA_CAT_DEEP.items()
               if k in [G.nodes[n]['category'] for n in G.nodes()]]
nes_handles = [
    matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgrey',
                             markeredgecolor=GOOD_DEEP, markeredgewidth=2.5, markersize=14,
                             label='UP in good (NES>0)'),
    matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgrey',
                             markeredgecolor=BAD_DEEP, markeredgewidth=2.5, markersize=14,
                             label='UP in poor (NES<0)'),
]
edge_handles = [
    matplotlib.lines.Line2D([0],[0], color='#5a6772', lw=1.0, alpha=0.55, label='Jaccard ~ 0.06'),
    matplotlib.lines.Line2D([0],[0], color='#5a6772', lw=3.0, alpha=0.55, label='Jaccard ~ 0.20'),
    matplotlib.lines.Line2D([0],[0], color='#5a6772', lw=5.5, alpha=0.55, label='Jaccard ≥ 0.40'),
]
ax.legend(handles=cat_handles + nes_handles + edge_handles, loc='upper left', fontsize=8.5,
          bbox_to_anchor=(1.0, 1.0), frameon=False, title='Category / direction / overlap',
          title_fontsize=9.5)

ax.set_axis_off()
save_panel(fig, 'Fig4D_enrichment_map', OUT)

print('\n=== Fig 4 v3.1 (4 revised panels) saved ===')
