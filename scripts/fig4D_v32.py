"""
Fig 4D v3.2 — enrichment map with denser edges (overlap coefficient threshold).
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import networkx as nx
import gseapy as gp
from adjustText import adjust_text
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle
import warnings; warnings.filterwarnings('ignore')

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
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
gsea_h = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
hm_sets = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')

def gseapy_pathway_key(name):
    if name.startswith('HALLMARK_'):
        return name.replace('HALLMARK_','').replace('_',' ').title()
    return name

# Use ALL Hallmark (50) so the network shows the full landscape
h_sel = gsea_h.copy()
h_sel['name'] = h_sel.pathway.str.replace('HALLMARK_','')
h_sel['category'] = h_sel.name.map(HALLMARK_CAT).fillna('Other')

G = nx.Graph()
for _, r in h_sel.iterrows():
    key = gseapy_pathway_key(r.pathway)
    gset = set(hm_sets.get(key, []))
    G.add_node(r['name'], nes=r.NES, pval=r.pval, size=len(gset),
               category=r.category, gene_set=gset)

# Use OVERLAP COEFFICIENT (Reimand Nat Protoc 2019 default — better for varied set sizes)
nodes = list(G.nodes())
OVERLAP_THR = 0.12
for i, n1 in enumerate(nodes):
    for n2 in nodes[i+1:]:
        s1 = G.nodes[n1]['gene_set']; s2 = G.nodes[n2]['gene_set']
        if not s1 or not s2: continue
        # Overlap coefficient: |A∩B| / min(|A|,|B|)
        overlap = len(s1 & s2) / min(len(s1), len(s2))
        if overlap >= OVERLAP_THR:
            G.add_edge(n1, n2, weight=overlap)

print(f'  4D revised: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (overlap coef ≥{OVERLAP_THR})')

fig, ax = plt.subplots(figsize=(13, 9.5))
pos = nx.spring_layout(G, seed=42, k=2.0, iterations=300, weight='weight')

# Edges first — VERY visible
edge_widths = [G[u][v]['weight'] * 9 for u, v in G.edges()]
edge_colors = ['#3a4a5c' for _ in G.edges()]
nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                        edge_color=edge_colors, alpha=0.65, style='solid')

# Nodes on top
node_sizes = [max(280, G.nodes[n]['size'] * 6) for n in G.nodes()]
node_colors = [GSEA_CAT_DEEP[G.nodes[n]['category']] for n in G.nodes()]
node_edgecolors = [GOOD_DEEP if G.nodes[n]['nes']>0 else BAD_DEEP for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
                        edgecolors=node_edgecolors, linewidths=2.5, alpha=0.95)

texts = []
for n in G.nodes():
    label = n.replace('_',' ').title()[:25]
    t = ax.text(pos[n][0], pos[n][1], label, fontsize=8.5, color='#0e2a47', fontweight='bold',
                ha='center', va='center')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.3),
                expand_points=(1.3, 1.4), force_points=0.4)
except: pass

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
    matplotlib.lines.Line2D([0],[0], color='#3a4a5c', lw=2.0, alpha=0.65, label=f'overlap = {OVERLAP_THR}'),
    matplotlib.lines.Line2D([0],[0], color='#3a4a5c', lw=4.0, alpha=0.65, label='overlap ≈ 0.45'),
    matplotlib.lines.Line2D([0],[0], color='#3a4a5c', lw=6.0, alpha=0.65, label='overlap ≥ 0.65'),
]
ax.legend(handles=cat_handles + nes_handles + edge_handles, loc='upper left', fontsize=8.5,
          bbox_to_anchor=(1.0, 1.0), frameon=False, title='Category / direction / overlap',
          title_fontsize=9.5)
ax.set_axis_off()
save_panel(fig, 'Fig4D_enrichment_map', OUT)

print('=== Fig 4D v3.2 saved ===')
