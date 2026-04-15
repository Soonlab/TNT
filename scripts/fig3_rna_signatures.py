"""
Figure 3 — RNA-seq signatures & DEG
  3A: Raincloud grid - 8 key signatures (pre, good vs bad)
  3B: Annotated volcano with pathway category coloring
  3C: Signature correlation network
  3D: 22-signature heatmap with multi-tier annotation
  3E: Signature ranking (lollipop with effect size)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import FancyBboxPatch
from adjustText import adjust_text

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
sig_stats = pd.read_csv(ROOT/'06_rna_immune/sig_response_stats.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')

sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index()
    sigs_m.columns = ['sample_id'] + list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

# ==========================================================
# 3A — Raincloud grid (8 key signatures)
# ==========================================================
key_sigs_list = [
    ('CD8_proliferation', 'CD8 proliferation'),
    ('CD8_activation', 'CD8 activation'),
    ('Antigen_presentation', 'Antigen presentation'),
    ('NLRC5_HLA_IFNG', 'NLRC5–HLA–IFNγ'),
    ('TLS_Cabrita', 'TLS (Cabrita)'),
    ('IFNg_Ayers_18', 'IFNγ (Ayers)'),
    ('TGFb_Mariathasan', 'TGFβ (Mariathasan)'),
    ('EMT_Mak', 'EMT (Mak)'),
]

fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharey=False)
pre = sigs_m[sigs_m.timepoint=='pre']
for ax, (sig, label) in zip(axes.flat, key_sigs_list):
    raincloud(ax, pre, 'response_bin', sig, ['good','bad'], PAL_RESP)
    g = pre[pre.response_bin=='good'][sig].dropna()
    b = pre[pre.response_bin=='bad'][sig].dropna()
    if len(g)>=3 and len(b)>=3:
        p = stats.mannwhitneyu(g, b).pvalue
        stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.15, p)
    ax.set_xlabel('')
    ax.set_ylabel('z-score' if ax in [axes[0,0], axes[1,0]] else '')
    ax.set_title(label, fontsize=10.5)
    annotate_count(ax, pre, 'response_bin', sig, ['good','bad'], dy=0.15)

fig.suptitle('Pre-treatment immune & stromal signatures (good vs poor responders)',
             fontsize=12, fontweight='bold', y=1.0, color='#1d3557')
fig.tight_layout()
save_panel(fig, 'Fig3A_raincloud_grid', OUT)

# ==========================================================
# 3B — Annotated volcano with pathway categories
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 6))
deg_p = deg.copy()
deg_p['-log10p'] = -np.log10(deg_p['pvalue'].replace(0, 1e-300))

# Manual category mapping
CC_GENES = {'MKI67','TOP2A','CCNB1','CCNB2','CDK1','CDC20','MCM2','MCM5','MCM7','BIRC5','CENPF','PLK1','AURKA','AURKB'}
DR_GENES = {'BRCA1','BRCA2','RAD51','RAD51B','PALB2','ATM','ATR','CHEK1','CHEK2','MRE11','RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2','BLM','BRIP1','EXO1','POLD1','POLE'}
EMT_GENES = {'VIM','CDH2','SNAI1','SNAI2','TWIST1','FN1','MMP2','MMP9','ZEB1','ZEB2','TGFB1','COL1A1','COL3A1','FAP','ACTA2','S100A4'}
IM_GENES = {'CD8A','CD8B','GZMA','GZMB','PRF1','IFNG','CD3D','CD3E','HLA-A','HLA-B','HLA-C','B2M','NLRC5','MS4A1','CD79A','CD79B','FOXP3'}

def cat(gene):
    if gene in CC_GENES: return 'Cell cycle'
    if gene in DR_GENES: return 'DNA repair'
    if gene in EMT_GENES: return 'EMT/Stromal'
    if gene in IM_GENES: return 'Immune'
    return 'Other'

deg_p['category'] = deg_p['gene'].apply(cat)
cat_colors = {'Cell cycle':'#2a9d8f','DNA repair':'#118ab2','EMT/Stromal':'#e76f51',
              'Immune':'#ef476f','Other':'#cccccc'}

# Background non-significant
ns = deg_p[(deg_p['pvalue']>=0.01) | (deg_p['log2FoldChange'].abs()<1)]
ax.scatter(ns.log2FoldChange, ns['-log10p'], s=4, alpha=0.25, color='#dee2e6')

# Significant
sig_df = deg_p[(deg_p['pvalue']<0.01) | (deg_p['log2FoldChange'].abs()>=1.5)]
for c in ['Cell cycle','DNA repair','Immune','EMT/Stromal','Other']:
    sub = sig_df[sig_df.category==c]
    if c=='Other':
        ax.scatter(sub.log2FoldChange, sub['-log10p'], s=10, alpha=0.55, color=cat_colors[c], edgecolor='none', label=c)
    else:
        ax.scatter(sub.log2FoldChange, sub['-log10p'], s=35, alpha=0.85, color=cat_colors[c],
                   edgecolor='white', linewidth=0.5, label=c, zorder=3)

# Label key genes
texts = []
to_label = sig_df[(sig_df.category != 'Other') & (sig_df['-log10p']>2.5)]
for _, r in to_label.head(25).iterrows():
    t = ax.text(r.log2FoldChange, r['-log10p'], r['gene'], fontsize=7.5,
                color='#1d3557', fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#6c757d', lw=0.4),
                expand_points=(1.2,1.2), force_points=0.3)
except: pass

ax.axvline(0, color='#1d3557', lw=0.5, alpha=0.3)
ax.axhline(-np.log10(0.01), color='#1d3557', lw=0.5, ls='--', alpha=0.3)
ax.text(ax.get_xlim()[1]*0.95, -np.log10(0.01)+0.05, 'p = 0.01', fontsize=8, color='#6c757d', ha='right')
ax.set_xlabel('log2 fold change (good vs bad)')
ax.set_ylabel('−log10(p-value)')
ax.set_title('Differential expression: pre-treatment good vs bad\n(genes colored by pathway category)')
ax.legend(loc='upper left', fontsize=8.5, title='Category', title_fontsize=9)
save_panel(fig, 'Fig3B_volcano_annotated', OUT)

# ==========================================================
# 3C — Signature correlation network
# ==========================================================
import networkx as nx
fig, ax = plt.subplots(figsize=(7, 6))
all_sig_cols = [c for c in sigs.columns if c in sigs_m.columns]
pre = sigs_m[sigs_m.timepoint=='pre'][all_sig_cols + ['response_bin']]
corr = pre[all_sig_cols].corr(method='spearman')

G = nx.Graph()
for s in all_sig_cols: G.add_node(s)
for i, s1 in enumerate(all_sig_cols):
    for j, s2 in enumerate(all_sig_cols):
        if i >= j: continue
        c = corr.loc[s1, s2]
        if abs(c) >= 0.4:
            G.add_edge(s1, s2, weight=abs(c), sign=np.sign(c))

# Layout
pos = nx.spring_layout(G, seed=42, k=1.6, iterations=100)

# Node color by response association sign (delta good - bad)
pre_stats = sig_stats[sig_stats.timepoint=='pre'].set_index('signature')
node_colors = []
node_sizes = []
for n in G.nodes():
    if n in pre_stats.index:
        d = pre_stats.loc[n, 'delta_good_minus_bad']
        p = pre_stats.loc[n, 'pvalue']
        node_colors.append(GOOD if d>0 else BAD)
        node_sizes.append(400 + 800*(-np.log10(p+1e-3))/3)
    else:
        node_colors.append('#cccccc'); node_sizes.append(300)

# Draw edges
edges = G.edges(data=True)
for u, v, d in edges:
    color = '#2a9d8f' if d['sign']>0 else '#e76f51'
    nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], edge_color=color, width=d['weight']*3, alpha=0.5, ax=ax)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, edgecolors='#1d3557',
                       linewidths=1, ax=ax, alpha=0.92)
labels = {n: n.replace('_',' ').replace(' Mariathasan','').replace(' Cabrita','').replace(' Ayers 18','')[:15] for n in G.nodes()}
nx.draw_networkx_labels(G, pos, labels, font_size=7.5, font_color='#1d3557', font_weight='bold', ax=ax)

# Legend
legend_items = [
    mpatches.Patch(color=GOOD, label='Up in good (pre)'),
    mpatches.Patch(color=BAD, label='Up in bad (pre)'),
    mpatches.Patch(color='#cccccc', label='No association'),
]
edge_legend = [matplotlib.lines.Line2D([0],[0], color='#2a9d8f', lw=2, label='Positive corr (ρ≥0.4)'),
               matplotlib.lines.Line2D([0],[0], color='#e76f51', lw=2, label='Negative corr (ρ≤-0.4)')]
ax.legend(handles=legend_items+edge_legend, loc='upper left', fontsize=8, bbox_to_anchor=(0.0, -0.02), ncol=2)
ax.set_axis_off()
ax.set_title('Signature correlation network (pre-treatment)\nnode size ∝ −log10(p) for response association', pad=8)
save_panel(fig, 'Fig3C_signature_network', OUT)

# ==========================================================
# 3D — 22-signature heatmap with multi-tier annotation
# ==========================================================
pre = sigs_m[sigs_m.timepoint=='pre']
sig_cols = [c for c in sigs.columns if c in sigs_m.columns]
pre_mat = pre.set_index('sample_id')[sig_cols]
pre_mat_z = pre_mat.sub(pre_mat.mean()).div(pre_mat.std())
order_samples = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
pre_mat_z = pre_mat_z.loc[order_samples]

# Annotation
resp_annot = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_samples], index=order_samples)
subj_annot = pd.Series([pre.set_index('sample_id').loc[s,'subject_id'] for s in order_samples], index=order_samples)
clin_idx = clin.set_index('subject_id')
ct_annot = pd.Series([clin_idx.loc[subj_annot[s], 'cT'] for s in order_samples], index=order_samples)
sex_annot = pd.Series([clin_idx.loc[subj_annot[s], 'sex'] for s in order_samples], index=order_samples)

col_colors_df = pd.DataFrame({
    'Response': resp_annot.map(PAL_RESP),
    'cT': ct_annot.map(PAL_STAGE),
    'Sex': sex_annot.map({'M':'#264653','F':'#e63946'}),
})

g = sns.clustermap(pre_mat_z.T, row_cluster=True, col_cluster=False, cmap='RdBu_r', center=0,
                    vmin=-2.5, vmax=2.5, col_colors=col_colors_df, figsize=(13, 8),
                    xticklabels=False, yticklabels=True, dendrogram_ratio=(0.06, 0.06),
                    cbar_kws={'label':'z-score','shrink':0.4},
                    linewidths=0.15, linecolor='#fafafa')
g.ax_heatmap.set_xlabel('Pre-treatment samples (n=33; sorted: good → bad)', fontsize=10)
g.ax_heatmap.set_ylabel('')
g.fig.suptitle('Pre-treatment 22-signature landscape', y=1.01, fontsize=12, fontweight='bold', color='#1d3557')
g.fig.savefig(OUT/'Fig3D_signature_heatmap.pdf', bbox_inches='tight')
g.fig.savefig(OUT/'Fig3D_signature_heatmap.png', dpi=400, bbox_inches='tight', facecolor='white')
plt.close('all')
print('  ✓ Fig3D_signature_heatmap')

# ==========================================================
# 3E — Lollipop chart of signature effect size
# ==========================================================
fig, ax = plt.subplots(figsize=(7.5, 7))
pre_st = sig_stats[sig_stats.timepoint=='pre'].sort_values('pvalue').head(20).iloc[::-1].reset_index(drop=True)
y = np.arange(len(pre_st))
colors = [GOOD if d>0 else BAD for d in pre_st.delta_good_minus_bad]
# Stems
for i, (_, r) in enumerate(pre_st.iterrows()):
    ax.plot([0, r.delta_good_minus_bad], [i, i], color='#1d3557', lw=1.2, alpha=0.5)
ax.scatter(pre_st.delta_good_minus_bad, y, s=[200*(-np.log10(p+1e-3))/3 for p in pre_st.pvalue],
           c=colors, edgecolor='#1d3557', linewidth=0.9, zorder=3, alpha=0.92)
# significance stars
for i, (_, r) in enumerate(pre_st.iterrows()):
    star = sig_symbol(r.pvalue)
    if star and star != 'ns':
        ax.text(r.delta_good_minus_bad + (0.05 if r.delta_good_minus_bad>0 else -0.05),
                i, star, va='center', ha='left' if r.delta_good_minus_bad>0 else 'right',
                fontsize=10, fontweight='bold', color='#1d3557')
ax.set_yticks(y)
ax.set_yticklabels([s.replace('_',' ') for s in pre_st.signature], fontsize=9)
ax.axvline(0, color='#1d3557', lw=0.8)
ax.set_xlabel('Δ z-score  (good − bad)')
ax.set_title('Signature ranking by response association (pre-treatment)\nbubble size ∝ −log10(p), color = direction')
add_axis_spines(ax, sides=('left','bottom'))
save_panel(fig, 'Fig3E_lollipop', OUT)

print('\n=== Fig 3 (5 panels) complete ===')
