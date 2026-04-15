"""
Figure 4 v3 — journal-style pathway enrichment.
Motifs:
  4A: Subramanian PNAS 2005 + Mootha Nat Genet 2003 — running enrichment score (ES) curves with barcode for top 4 pathways
  4B: Litchfield Cell 2021 + Bagaev Cancer Cell 2021 — Hallmark bubble plot (NES × −log10p × size × category)
  4C: Wu Cell 2024 + clusterProfiler — top 30 pathway dotplot (categorized)
  4D: Reimand Nat Protoc 2019 + Subramanian — Enrichment map network (Jaccard overlap)
  4E: Hänzelmann GSVA 2013 + Bagaev TME — ssGSEA pathway × sample heatmap clustered
  4F: Mariathasan Nature 2018 — Hallmark category NES distribution box

Style: no titles, saturated colors, no overlapping text.
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

# Saturated palette
GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}
PAL_STAGE_DEEP = {'T2':'#7fb0c4', 'T2/T3':'#1c5d7e', 'T3':'#0e2a47', 'T4':'#a01b2b'}

# Hallmark categories (Liberzon 2015 + extension)
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
ssg = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce')
ssg_stats = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_response_stats.tsv', sep='\t')

# ============================================================
# 4A — Running enrichment score plots (Subramanian/Mootha iconic style)
# ============================================================
# Compute ranked stat = log2FC (descending)
deg_p = deg.dropna(subset=['log2FoldChange','gene']).copy()
deg_p = deg_p.sort_values('log2FoldChange', ascending=False).reset_index(drop=True)
ranks = pd.Series(deg_p.log2FoldChange.values, index=deg_p.gene.values)

# Get Hallmark gene sets via gseapy
hm_sets = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')

# Top 4 pathways: mix of UP (E2F, G2M) + DOWN (EMT, Myogenesis)
top_paths_for_es = [
    ('HALLMARK_E2F_TARGETS', 'Hallmark · E2F targets'),
    ('HALLMARK_G2M_CHECKPOINT', 'Hallmark · G2M checkpoint'),
    ('HALLMARK_MYC_TARGETS_V1', 'Hallmark · MYC targets V1'),
    ('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'Hallmark · EMT'),
]

def gseapy_pathway_key(name):
    """Map fgsea result name (with HALLMARK_ prefix) to gseapy library key."""
    if name.startswith('HALLMARK_'):
        # gseapy library uses Title Case format like 'E2F Targets'
        return name.replace('HALLMARK_','').replace('_',' ').title()
    return name

def running_es(ranks, gene_set):
    """Compute weighted Kolmogorov-Smirnov running ES (Subramanian PNAS 2005)."""
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

fig, axes = plt.subplots(4, 4, figsize=(15, 11), sharex=True,
                         gridspec_kw={'height_ratios':[3, 0.4, 0.5, 0.0001], 'hspace':0.05})
# We'll only use rows 0,1,2 per pathway-column
for col, (pw, label) in enumerate(top_paths_for_es):
    key = gseapy_pathway_key(pw)
    gene_set = set(hm_sets.get(key, []))
    if not gene_set:
        # Try alternative formats
        for k in hm_sets:
            if k.lower().replace(' ','_').replace('-','_') == pw.replace('HALLMARK_','').lower():
                gene_set = set(hm_sets[k]); break
    es, in_set = running_es(ranks, gene_set)
    if es is None:
        for r in range(3): axes[r, col].set_visible(False)
        continue
    # Compute NES from fgsea result for label
    row = gsea_h[gsea_h.pathway==pw]
    nes_val = row.NES.iloc[0] if len(row) else 0
    pval = row.pval.iloc[0] if len(row) else 1
    color = GOOD_DEEP if nes_val>0 else BAD_DEEP

    # Top: running ES curve
    ax = axes[0, col]
    ax.plot(np.arange(len(es)), es, color=color, lw=2.4)
    # Mark peak
    if nes_val>0:
        peak = np.argmax(es); peak_y = es[peak]
    else:
        peak = np.argmin(es); peak_y = es[peak]
    ax.axhline(0, color='#0e2a47', lw=0.6, alpha=0.5)
    ax.fill_between(np.arange(len(es)), 0, es, color=color, alpha=0.20)
    ax.scatter([peak],[peak_y], color=color, s=80, edgecolor='white', linewidth=1.5, zorder=5)
    ax.set_ylabel('Enrichment\nscore (ES)', fontsize=10, color='#0e2a47', fontweight='bold')
    ax.set_xlim(0, len(es))
    add_axis_spines(ax)
    # In-plot annotation: pathway name + NES + p
    ax.text(0.02, 0.96, label, transform=ax.transAxes, fontsize=10.5, fontweight='bold',
            color='#0e2a47', va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.35'))
    ax.text(0.98, 0.04 if nes_val>0 else 0.96, f'NES = {nes_val:+.2f}\np = {pval:.2g}',
            transform=ax.transAxes, fontsize=9.5, color=color, fontweight='bold',
            ha='right', va='bottom' if nes_val>0 else 'top')

    # Middle: barcode (gene set member positions)
    ax_b = axes[1, col]
    hit_positions = np.where(in_set==1)[0]
    for hp in hit_positions:
        ax_b.axvline(hp, color=color, lw=0.6, alpha=0.85)
    ax_b.set_xlim(0, len(es)); ax_b.set_yticks([]); ax_b.set_xticks([])
    for s in ['top','right','left','bottom']: ax_b.spines[s].set_visible(False)
    ax_b.text(-0.01, 0.5, 'Members', transform=ax_b.transAxes, ha='right', va='center',
              fontsize=8.5, color='#0e2a47')

    # Bottom: ranked metric heatmap (red→blue gradient)
    ax_r = axes[2, col]
    rank_arr = np.array([ranks.values])
    rank_norm = np.clip(rank_arr, -3, 3)
    cmap_rank = LinearSegmentedColormap.from_list('rk', [BAD_DEEP,'white',GOOD_DEEP])
    ax_r.imshow(rank_norm, aspect='auto', cmap=cmap_rank, extent=[0, len(es), 0, 1])
    ax_r.set_yticks([]); ax_r.set_xlim(0, len(es))
    ax_r.set_xlabel('Gene rank  (high LFC ← →  low LFC)', fontsize=9.5, color='#0e2a47')
    ax_r.tick_params(labelsize=8)
    for s in ['top','right','left']: ax_r.spines[s].set_visible(False)
    # Hide unused row
    axes[3, col].set_visible(False)

# Hide unused row 3
fig.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.06, hspace=0.05, wspace=0.18)
save_panel(fig, 'Fig4A_running_ES', OUT)

# ============================================================
# 4B — Hallmark bubble plot (Litchfield/Bagaev style)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 7.5))
h = gsea_h.copy()
h['name'] = h.pathway.str.replace('HALLMARK_','')
h['category'] = h.name.map(HALLMARK_CAT).fillna('Other')
h['neglogp'] = -np.log10(h.pval.replace(0, 1e-300))

for cat, color in GSEA_CAT_DEEP.items():
    sub = h[h.category==cat]
    if len(sub)==0: continue
    sizes = sub['size'] * 1.5
    ax.scatter(sub.NES, sub.neglogp, s=sizes, c=color, alpha=0.85,
               edgecolor='#0e2a47', linewidth=0.8, label=cat, zorder=3)

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
ax.legend(loc='upper left', fontsize=9, title='Category', title_fontsize=10, ncol=1,
          bbox_to_anchor=(1.01, 1), frameon=False)

# Size legend
size_examples = [50, 100, 200]
for i, sz in enumerate(size_examples):
    ax.scatter([], [], s=sz*1.5, color='#5a6772', edgecolor='#0e2a47', linewidth=0.5,
               label=f'{sz} genes', alpha=0.7)
add_axis_spines(ax)
save_panel(fig, 'Fig4B_Hallmark_bubble', OUT)

# ============================================================
# 4C — Top pathway dotplot (clusterProfiler/Wu Cell 2024 style)
# ============================================================
# Combine top Hallmark + Reactome
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
    if 'repair' in n or 'recombination' in n or 'replication' in n or 'damage' in n or 'hrd' in n.lower() or 'hdr' in n: return 'DNA repair'
    if 'immune' in n or 'antigen' in n or 'mhc' in n or 'interferon' in n or 'inflammatory' in n: return 'Immune'
    if 'extracellular matrix' in n or 'ecm' in n or 'collagen' in n: return 'Stromal/EMT'
    if 'metabol' in n or 'glycolysis' in n or 'oxidative' in n: return 'Metabolism'
    if 'signal' in n: return 'Signaling'
    if 'apopt' in n: return 'Apoptosis'
    return 'Other'
r_top['category'] = r_top.name.apply(reactome_cat)

combined = pd.concat([h_top.assign(label=h_top.name.str.replace('_',' ').str.title()),
                      r_top.assign(label=r_top.name)], ignore_index=True)

# Top 25 by significance, across both
combined = combined.sort_values('pval').head(25).reset_index(drop=True)
combined = combined.sort_values('NES').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, 8))
y = np.arange(len(combined))
x = combined.NES
sizes = combined['size'] * 1.5
colors = [GSEA_CAT_DEEP[c] for c in combined.category]

# Vertical reference at 0
ax.axvline(0, color='#0e2a47', lw=0.8)
# Stems
for i, r in combined.iterrows():
    color = GSEA_CAT_DEEP[r.category]
    ax.plot([0, r.NES], [i, i], color=color, lw=1.0, alpha=0.4)

scatter = ax.scatter(x, y, s=sizes, c=colors, alpha=0.92, edgecolor='#0e2a47', linewidth=0.7, zorder=3)

ax.set_yticks(y)
labels = []
for _, r in combined.iterrows():
    src_marker = '◇' if r.source=='Hallmark' else '○'
    labels.append(f'{src_marker} {r.label[:48]}')
ax.set_yticklabels(labels, fontsize=9, color='#0e2a47')
ax.set_xlabel('Normalized enrichment score (NES)', fontsize=11, fontweight='bold', color='#0e2a47')

# Legend: category + source
cat_handles = [mpatches.Patch(color=c, label=k) for k, c in GSEA_CAT_DEEP.items()
               if k in combined.category.unique()]
src_handles = [matplotlib.lines.Line2D([0],[0], marker='d', color='w', markerfacecolor='#5a6772', markersize=8, label='Hallmark'),
               matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='#5a6772', markersize=8, label='Reactome')]
leg1 = ax.legend(handles=cat_handles, loc='lower right', fontsize=9, frameon=False,
                 bbox_to_anchor=(0.99, 0.02), title='Category', title_fontsize=9.5)
ax.add_artist(leg1)

add_axis_spines(ax)
save_panel(fig, 'Fig4C_pathway_dotplot', OUT)

# ============================================================
# 4D — Pathway enrichment map (network, Reimand Nat Protoc 2019)
# ============================================================
import networkx as nx

# Use Hallmark only (manageable size)
h_sel = gsea_h[gsea_h.pval < 0.05].copy()
h_sel['name'] = h_sel.pathway.str.replace('HALLMARK_','')
h_sel['category'] = h_sel.name.map(HALLMARK_CAT).fillna('Other')

G = nx.Graph()
for _, r in h_sel.iterrows():
    key = gseapy_pathway_key(r.pathway)
    gset = set(hm_sets.get(key, []))
    G.add_node(r['name'], nes=r.NES, pval=r.pval, size=len(gset),
               category=r.category, gene_set=gset)

# Edges by Jaccard overlap
nodes = list(G.nodes())
for i, n1 in enumerate(nodes):
    for n2 in nodes[i+1:]:
        s1 = G.nodes[n1]['gene_set']; s2 = G.nodes[n2]['gene_set']
        if not s1 or not s2: continue
        jac = len(s1 & s2) / len(s1 | s2)
        if jac >= 0.10:
            G.add_edge(n1, n2, weight=jac)

fig, ax = plt.subplots(figsize=(11, 8.5))
pos = nx.spring_layout(G, seed=42, k=2.5, iterations=200, weight='weight')

# Draw edges first
for u, v, d in G.edges(data=True):
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color='#aab3bf', lw=d['weight']*5, alpha=0.55, zorder=1)

# Draw nodes: size ∝ gene set size, color by category, edge by NES sign
for n, d in G.nodes(data=True):
    color = GSEA_CAT_DEEP[d['category']]
    edge_color = GOOD_DEEP if d['nes']>0 else BAD_DEEP
    size = max(180, d['size'] * 6)
    ax.scatter(pos[n][0], pos[n][1], s=size, c=color, edgecolor=edge_color,
               linewidth=2.5, alpha=0.92, zorder=3)

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

# Legend in corner
cat_handles = [mpatches.Patch(color=c, label=k) for k, c in GSEA_CAT_DEEP.items()
               if k in [G.nodes[n]['category'] for n in G.nodes()]]
nes_handles = [matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgrey',
                                        markeredgecolor=GOOD_DEEP, markeredgewidth=2.5, markersize=12, label='UP in good (NES>0)'),
               matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgrey',
                                        markeredgecolor=BAD_DEEP, markeredgewidth=2.5, markersize=12, label='UP in poor (NES<0)')]
ax.legend(handles=cat_handles+nes_handles, loc='upper left', fontsize=8.5,
          bbox_to_anchor=(1.0, 1.0), frameon=False, title='Category / direction', title_fontsize=9.5)

ax.set_axis_off()
save_panel(fig, 'Fig4D_enrichment_map', OUT)

# ============================================================
# 4E — ssGSEA pathway × sample heatmap (GSVA / Bagaev style)
# ============================================================
ssg_m = ssg.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in ssg_m.columns:
    ssg_m = ssg.reset_index(); ssg_m.columns = ['sample_id'] + list(ssg_m.columns[1:])
ssg_m = ssg_m.merge(rna_inv[['sample_id','timepoint','response_bin']], on='sample_id')

pre = ssg_m[ssg_m.timepoint=='pre']
top30 = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(30).pathway.tolist()
present = [p for p in top30 if p in pre.columns]
mat = pre.set_index('sample_id')[present]
mat_z = mat.sub(mat.mean()).div(mat.std())
order_samples = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
mat_z = mat_z.loc[order_samples]
n_s = len(order_samples); n_p = len(present)

# Cluster pathways
from scipy.cluster.hierarchy import linkage, leaves_list
Z = linkage(mat_z.T.values, method='average', metric='correlation')
order_idx = leaves_list(Z)
present_ordered = [present[i] for i in order_idx]
mat_z = mat_z[present_ordered]

# Categorize pathway
def pw_cat(p):
    pl = p.lower()
    if 'cycle' in pl or 'mitotic' in pl or 'phase' in pl or 'cdk' in pl: return 'Proliferation'
    if 'repair' in pl or 'recombination' in pl or 'damage' in pl or 'hdr' in pl or 'replication' in pl: return 'DNA repair'
    if 'immune' in pl or 'antigen' in pl or 'mhc' in pl or 'interferon' in pl or 'inflammatory' in pl: return 'Immune'
    if 'extracellular' in pl or 'ecm' in pl or 'collagen' in pl or 'epithelial' in pl: return 'Stromal/EMT'
    if 'metabol' in pl or 'glycolysis' in pl or 'oxidative' in pl: return 'Metabolism'
    if 'signal' in pl: return 'Signaling'
    if 'myc' in pl or 'e2f' in pl or 'g2m' in pl: return 'Proliferation'
    return 'Other'

fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(3, 5, height_ratios=[0.18, 5.5, 0.7],
                      width_ratios=[0.05, 1.6, 0.18, 6.0, 1.5],
                      hspace=0.06, wspace=0.0)

# Top annotation (response only here)
def ann_strip(ax, vals, palette_map):
    arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[0, n_s, 0, 1])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, n_s)

resp_a = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_samples])
ax_ann = fig.add_subplot(gs[0, 3])
ann_strip(ax_ann, resp_a.tolist(), PAL)
ax_ann_lbl = fig.add_subplot(gs[0, 1])
ax_ann_lbl.axis('off')
ax_ann_lbl.text(1.0, 0.5, 'Response', ha='right', va='center', fontsize=11,
                color='#0e2a47', fontweight='bold', transform=ax_ann_lbl.transAxes)

# Pathway labels (col 1) + module color strip (col 2) + heatmap (col 3)
ax_pw_lbl = fig.add_subplot(gs[1, 1])
ax_pw_lbl.axis('off')
ax_pw_lbl.set_xlim(0, 1); ax_pw_lbl.set_ylim(0, n_p)
ax_pw_lbl.invert_yaxis()
for i, p in enumerate(present_ordered):
    short = p.split(' R-HSA')[0][:48]
    ax_pw_lbl.text(1.0, i+0.5, short, ha='right', va='center', fontsize=8.5, color='#0e2a47')

ax_mod = fig.add_subplot(gs[1, 2])
mod_colors = [GSEA_CAT_DEEP[pw_cat(p)] for p in present_ordered]
arr_mod = np.array([[to_rgb(c) for c in mod_colors]])
ax_mod.imshow(arr_mod.transpose(1,0,2), aspect='auto', interpolation='nearest', extent=[0,1,0,n_p])
ax_mod.set_xticks([]); ax_mod.set_yticks([])
for s in ['top','right','left','bottom']: ax_mod.spines[s].set_visible(False)
ax_mod.invert_yaxis()

ax_h = fig.add_subplot(gs[1, 3])
im = ax_h.imshow(mat_z.values.T, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                 aspect='auto', interpolation='nearest', extent=[0, n_s, 0, n_p])
ax_h.set_xticks([]); ax_h.set_yticks([])
ax_h.invert_yaxis()
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# Right column: colorbar (top) + module legend (below)
cax = fig.add_subplot(gs[1, 4])
# Sub-divide cax into top half (colorbar) + bottom half (legend)
# Use insert axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cb_ax = inset_axes(cax, width="40%", height="20%", loc='upper left', bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=cax.transAxes)
cb = fig.colorbar(im, cax=cb_ax, orientation='vertical')
cb.set_label('z-score', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.8)

# Module legend below
cax.axis('off')
mod_present = []
for p in present_ordered:
    c = pw_cat(p)
    if c not in mod_present: mod_present.append(c)
y_step = 0.07
y_start = 0.50
for i, mod in enumerate(mod_present):
    y = y_start - i*y_step
    cax.add_patch(Rectangle((0.05, y-0.02), 0.18, 0.04, color=GSEA_CAT_DEEP[mod],
                             transform=cax.transAxes))
    cax.text(0.27, y, mod, fontsize=9.5, va='center', transform=cax.transAxes, color='#0e2a47')
cax.text(0.05, y_start+0.06, 'Pathway category', fontsize=10, color='#0e2a47', fontweight='bold',
         transform=cax.transAxes)

# Bottom legend (response)
ax_clin = fig.add_subplot(gs[2, 1:4])
ax_clin.axis('off')
ax_clin.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good'),
                        mpatches.Patch(color=BAD_DEEP, label='Poor')],
               loc='center', fontsize=10, ncol=2, frameon=False)

save_panel(fig, 'Fig4E_ssGSEA_heatmap', OUT)

# ============================================================
# 4F — Hallmark category NES distribution (Mariathasan style)
# ============================================================
fig, ax = plt.subplots(figsize=(9, 5.5))
h_box = h.copy()
cat_order = ['Proliferation','DNA repair','Apoptosis','Immune','Signaling','Metabolism','Stress','Stromal/EMT','Other','TP53']
cat_present = [c for c in cat_order if c in h_box.category.values]

# Box + jitter per category
for i, cat in enumerate(cat_present):
    sub = h_box[h_box.category==cat]
    color = GSEA_CAT_DEEP[cat]
    bp = ax.boxplot([sub.NES], positions=[i], widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color='white', linewidth=1.6),
                    boxprops=dict(facecolor=color, alpha=0.9, edgecolor=color, linewidth=0.8),
                    whiskerprops=dict(color=color, linewidth=0.8),
                    capprops=dict(color=color, linewidth=0.8))
    # Jitter
    jx = i + np.random.uniform(-0.15, 0.15, len(sub))
    ax.scatter(jx, sub.NES, s=40, color=color, alpha=0.85, edgecolor='#0e2a47', linewidth=0.5, zorder=3)

ax.axhline(0, color='#0e2a47', lw=0.9)
ax.axhspan(-0.3, 0.3, color='#dee2e6', alpha=0.4, zorder=0)
ax.set_xticks(range(len(cat_present)))
ax.set_xticklabels(cat_present, rotation=30, ha='right', fontsize=10, color='#0e2a47')
ax.set_ylabel('Normalized enrichment score (NES)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.text(0.99, 0.97, 'NES > 0  ↑ Good\nNES < 0  ↑ Poor', transform=ax.transAxes,
        ha='right', va='top', fontsize=9.5, color='#0e2a47', fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.92, boxstyle='round,pad=0.4'))
add_axis_spines(ax)
save_panel(fig, 'Fig4F_category_NES_box', OUT)

print('\n=== Fig 4 v3 (6 journal-style panels) saved ===')
print(f'Output: {OUT}')
