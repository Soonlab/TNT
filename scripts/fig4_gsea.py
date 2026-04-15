"""
Figure 4 — GSEA & pathway enrichment
  4A: Hallmark GSEA bubble plot (NES vs -log10p, size=geneset, color=category)
  4B: Reactome top pathway lollipop with NES bars
  4C: ssGSEA top pathway grid (raincloud)
  4D: ssGSEA pathway heatmap (top 30 × 33 samples, multi-tier annotation)
  4E: Pathway category summary - radial plot
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from adjustText import adjust_text

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
gsea_h = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
gsea_r = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')
ssg = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce')
ssg_stats = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_response_stats.tsv', sep='\t')

# ===========================================================
# 4A — Hallmark GSEA bubble plot
# ===========================================================
HALLMARK_CAT = {
    'E2F_TARGETS':'Proliferation','G2M_CHECKPOINT':'Proliferation','MYC_TARGETS_V1':'Proliferation','MYC_TARGETS_V2':'Proliferation',
    'MITOTIC_SPINDLE':'Proliferation','MTORC1_SIGNALING':'Signaling',
    'DNA_REPAIR':'DNA repair','UV_RESPONSE_DN':'DNA repair','UNFOLDED_PROTEIN_RESPONSE':'Stress',
    'EPITHELIAL_MESENCHYMAL_TRANSITION':'Stromal/EMT','MYOGENESIS':'Stromal/EMT','APICAL_JUNCTION':'Stromal/EMT',
    'TGF_BETA_SIGNALING':'Stromal/EMT','HEDGEHOG_SIGNALING':'Signaling','WNT_BETA_CATENIN_SIGNALING':'Signaling',
    'NOTCH_SIGNALING':'Signaling','PI3K_AKT_MTOR_SIGNALING':'Signaling','KRAS_SIGNALING_UP':'Signaling',
    'IL2_STAT5_SIGNALING':'Immune','IL6_JAK_STAT3_SIGNALING':'Immune','INTERFERON_GAMMA_RESPONSE':'Immune',
    'INTERFERON_ALPHA_RESPONSE':'Immune','INFLAMMATORY_RESPONSE':'Immune','TNFA_SIGNALING_VIA_NFKB':'Immune',
    'COMPLEMENT':'Immune','ALLOGRAFT_REJECTION':'Immune','COAGULATION':'Immune',
    'OXIDATIVE_PHOSPHORYLATION':'Metabolism','GLYCOLYSIS':'Metabolism','HYPOXIA':'Metabolism',
    'CHOLESTEROL_HOMEOSTASIS':'Metabolism','FATTY_ACID_METABOLISM':'Metabolism','XENOBIOTIC_METABOLISM':'Metabolism',
    'ADIPOGENESIS':'Metabolism','PEROXISOME':'Metabolism','HEME_METABOLISM':'Metabolism','BILE_ACID_METABOLISM':'Metabolism',
    'ESTROGEN_RESPONSE_EARLY':'Signaling','ESTROGEN_RESPONSE_LATE':'Signaling','ANDROGEN_RESPONSE':'Signaling',
    'P53_PATHWAY':'Other','APOPTOSIS':'Apoptosis','PROTEIN_SECRETION':'Other','UV_RESPONSE_UP':'Other',
    'ANGIOGENESIS':'Other','PANCREAS_BETA_CELLS':'Other','SPERMATOGENESIS':'Other',
    'REACTIVE_OXYGEN_SPECIES_PATHWAY':'Stress','KRAS_SIGNALING_DN':'Signaling',
}

fig, ax = plt.subplots(figsize=(9, 7))
h = gsea_h.copy()
h['name'] = h.pathway.str.replace('HALLMARK_','')
h['category'] = h.name.map(HALLMARK_CAT).fillna('Other')
h['neglogp'] = -np.log10(h.pval.replace(0, 1e-300))

# Plot all
for cat, color in GSEA_CAT.items():
    sub = h[h.category==cat]
    if len(sub)==0: continue
    sizes = sub['size'] * 1.5
    ax.scatter(sub.NES, sub.neglogp, s=sizes, c=color, alpha=0.78,
               edgecolor='#1d3557', linewidth=0.8, label=cat, zorder=3)

# Label significant pathways
texts = []
sig_h = h[h.pval < 1e-3]
for _, r in sig_h.iterrows():
    name_short = r['name'].replace('_',' ').title()[:25]
    t = ax.text(r.NES, r.neglogp, name_short, fontsize=7.5, color='#1d3557', fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#6c757d', lw=0.4))
except: pass

ax.axvline(0, color='#1d3557', lw=0.6, alpha=0.5)
ax.axhline(-np.log10(0.05), color='#6c757d', lw=0.6, ls='--', alpha=0.5)
ax.text(ax.get_xlim()[1]*0.98, -np.log10(0.05)+0.2, 'p = 0.05', fontsize=8, color='#6c757d', ha='right')
ax.set_xlabel('Normalized enrichment score (NES)\n← UP in bad        UP in good →')
ax.set_ylabel('−log10(p-value)')
ax.set_title('Hallmark GSEA — pre-treatment good vs bad responders\n(bubble size ∝ gene set size)')
ax.legend(loc='upper left', fontsize=8, title='Category', title_fontsize=9, ncol=1, bbox_to_anchor=(1.02, 1))
save_panel(fig, 'Fig4A_GSEA_bubble', OUT)

# ===========================================================
# 4B — Reactome top pathway lollipop
# ===========================================================
fig, ax = plt.subplots(figsize=(8.5, 7))
r = gsea_r.nsmallest(20, 'pval').sort_values('NES').reset_index(drop=True)
y = np.arange(len(r))
colors = [GOOD if n>0 else BAD for n in r.NES]
# Stems
for i, row in r.iterrows():
    ax.plot([0, row.NES], [i, i], color='#1d3557', lw=1.2, alpha=0.4)
ax.scatter(r.NES, y, s=[300*(-np.log10(p+1e-300))/12 for p in r.pval],
           c=colors, edgecolor='#1d3557', linewidth=1, alpha=0.9, zorder=3)
ax.axvline(0, color='#1d3557', lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([t.replace('REACTOME_','').replace('_',' ').title()[:55] for t in r.pathway], fontsize=8)
ax.set_xlabel('Normalized enrichment score (NES)')
ax.set_title('Top 20 Reactome pathways (good vs bad)\nbubble size ∝ −log10(p)')
add_axis_spines(ax)
save_panel(fig, 'Fig4B_Reactome_lollipop', OUT)

# ===========================================================
# 4C — ssGSEA top 6 pathway raincloud grid
# ===========================================================
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
ssg_m = ssg.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in ssg_m.columns:
    ssg_m = ssg.reset_index(); ssg_m.columns = ['sample_id'] + list(ssg_m.columns[1:])
ssg_m = ssg_m.merge(rna_inv[['sample_id','timepoint','response_bin']], on='sample_id')
top_pre_pw = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(6).pathway.tolist()
for ax, pw in zip(axes.flat, top_pre_pw):
    if pw not in ssg_m.columns: ax.set_visible(False); continue
    pre = ssg_m[ssg_m.timepoint=='pre']
    raincloud(ax, pre, 'response_bin', pw, ['good','bad'], PAL_RESP)
    g = pre[pre.response_bin=='good'][pw].dropna()
    b = pre[pre.response_bin=='bad'][pw].dropna()
    if len(g)>=3 and len(b)>=3:
        p = stats.mannwhitneyu(g, b).pvalue
        stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.005, p, h=0.003)
    ax.set_xlabel(''); ax.set_ylabel('ssGSEA z-score' if ax in [axes[0,0], axes[1,0]] else '')
    title = pw.split(' R-HSA')[0][:42]
    ax.set_title(title, fontsize=10)

fig.suptitle('Top 6 ssGSEA pathways — pre-treatment good vs bad', fontsize=12, fontweight='bold', y=1.0)
fig.tight_layout()
save_panel(fig, 'Fig4C_ssGSEA_grid', OUT)

# ===========================================================
# 4D — ssGSEA pathway heatmap (top 30 × 33 samples)
# ===========================================================
pre = ssg_m[ssg_m.timepoint=='pre']
top30 = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(30).pathway.tolist()
present = [p for p in top30 if p in pre.columns]
mat = pre.set_index('sample_id')[present]
mat_z = mat.sub(mat.mean()).div(mat.std())
order_samples = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
mat_z = mat_z.loc[order_samples]

# multi-tier annotation
resp_a = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_samples], index=order_samples)
col_colors_df = pd.DataFrame({'Response': resp_a.map(PAL_RESP)})

g = sns.clustermap(mat_z.T, row_cluster=True, col_cluster=False, cmap='RdBu_r', center=0,
                    vmin=-2.5, vmax=2.5, col_colors=col_colors_df, figsize=(12, 9),
                    xticklabels=False, yticklabels=[p.split(' R-HSA')[0][:50] for p in present],
                    dendrogram_ratio=(0.07, 0.03), cbar_kws={'label':'z-score','shrink':0.4},
                    linewidths=0.1, linecolor='#fafafa')
g.ax_heatmap.set_xlabel('Pre-treatment samples (sorted: good → bad)')
g.fig.suptitle('Top 30 ssGSEA pathways — response-stratified', y=1.01, fontsize=12, fontweight='bold', color='#1d3557')
g.fig.savefig(OUT/'Fig4D_ssGSEA_heatmap.pdf', bbox_inches='tight')
g.fig.savefig(OUT/'Fig4D_ssGSEA_heatmap.png', dpi=400, bbox_inches='tight', facecolor='white')
plt.close('all')
print('  ✓ Fig4D_ssGSEA_heatmap')

# ===========================================================
# 4E — Pathway category summary radar plot
# ===========================================================
# Aggregate pre-treatment NES by category
h = gsea_h.copy()
h['name'] = h.pathway.str.replace('HALLMARK_','')
h['category'] = h.name.map(HALLMARK_CAT).fillna('Other')
cat_summary = h.groupby('category')['NES'].mean().reindex(['Proliferation','DNA repair','Immune','Stromal/EMT','Metabolism','Signaling','Apoptosis','Stress']).fillna(0)

# Polar plot
fig = plt.figure(figsize=(6.5, 6.5))
ax = fig.add_subplot(111, projection='polar')
N = len(cat_summary)
theta = np.linspace(0, 2*np.pi, N, endpoint=False)
values = cat_summary.values

# Plot
ax.fill_between(theta, np.minimum(values, 0), 0, color=BAD, alpha=0.45)
ax.fill_between(theta, 0, np.maximum(values, 0), color=GOOD, alpha=0.45)
# Outline
ax.plot(np.append(theta, theta[0]), np.append(values, values[0]), color='#1d3557', lw=1.5, marker='o', markersize=8)
# Reference circle at 0
ax.plot(np.linspace(0, 2*np.pi, 100), [0]*100, color='#1d3557', lw=1, ls='--', alpha=0.5)
ax.set_xticks(theta)
ax.set_xticklabels(cat_summary.index, fontsize=10)
ax.set_ylim(min(values)-0.3, max(values)+0.3)
ax.set_yticks([-2,-1,0,1,2])
ax.set_yticklabels([f'{v:+d}' for v in [-2,-1,0,1,2]], fontsize=8)
ax.set_title('Hallmark category enrichment summary (mean NES per category)',
             fontsize=11, fontweight='bold', pad=22, color='#1d3557')
# Legend
fig.text(0.5, 0.02, '↑ outward = UP in good responders ; ↓ inward = UP in bad', ha='center', fontsize=9, color='#1d3557', style='italic')
save_panel(fig, 'Fig4E_radar_categories', OUT)

print('\n=== Fig 4 (5 panels) complete ===')
