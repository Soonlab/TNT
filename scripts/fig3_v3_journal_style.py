"""
Figure 3 v3 — Journal-style RNA-seq immune & stromal signatures.
Motifs:
  3A: Thorsson Immunity 2018 + Bagaev Cancer Cell 2021 — TME radar (signatures per response group)
  3B: Mariathasan Nature 2018 — multi-tier annotated signature heatmap
  3C: Jerby-Arnon Cell 2018 — pathway-categorized volcano with labeled top genes
  3D: Litchfield Cell 2021 — effect-size forest with bootstrap CI
  3E: Tirosh Science 2016 + Sade-Feldman Cell 2018 — bivariate CD8 functional state
  3F: Cabrita Nature 2020 — TLS gene-level heatmap with response panel
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap
from adjustText import adjust_text

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
sig_stats = pd.read_csv(ROOT/'06_rna_immune/sig_response_stats.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
tpm = pd.read_csv(ROOT/'06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')

sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index(); sigs_m.columns = ['sample_id'] + list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

UNMATCHED = [13,15,16,17,18,19,33]

# Signature → category mapping (functional grouping)
SIG_CATEGORY = {
    'CD8_proliferation':'Cytotoxic T-cell',
    'CD8_activation':'Cytotoxic T-cell',
    'CD8_exhaustion':'Cytotoxic T-cell',
    'Cytolytic_activity':'Cytotoxic T-cell',
    'Antigen_presentation':'Antigen presentation',
    'MHC_II':'Antigen presentation',
    'NLRC5_HLA_IFNG':'Antigen presentation',
    'IFNg_Ayers_18':'IFN response',
    'TLS_Cabrita':'B-cell / TLS',
    'B_cell':'B-cell / TLS',
    'NK_cell':'Innate',
    'Mac_M1':'Innate',
    'Mac_M2':'Innate',
    'Treg':'Regulatory',
    'Checkpoint_inhibitory':'Regulatory',
    'TGFb_Mariathasan':'Stromal/EMT',
    'EMT_Mak':'Stromal/EMT',
    'CAF_iCAF':'Stromal/EMT',
    'CAF_myCAF':'Stromal/EMT',
    'Hypoxia_Buffa':'Hypoxia',
    'Stemness_mRNAsi_proxy':'Other',
    'Epithelial':'Other',
}
CAT_COLORS = {
    'Cytotoxic T-cell':'#0f8b78',
    'Antigen presentation':'#118ab2',
    'IFN response':'#06aed5',
    'B-cell / TLS':'#ffd23f',
    'Innate':'#f4a261',
    'Regulatory':'#9467bd',
    'Stromal/EMT':'#e76f51',
    'Hypoxia':'#1d3557',
    'Other':'#8d99ae',
}

# ============================================================
# FIG 3A — TME radar (Thorsson Immunity 2018 / Bagaev 2021 style)
# ============================================================
fig = plt.figure(figsize=(7.5, 7.5))
ax = fig.add_subplot(111, projection='polar')

# Curate 12 high-level signatures for radar
radar_sigs = [
    ('CD8_proliferation','CD8 prolif'),
    ('CD8_activation','CD8 activation'),
    ('Cytolytic_activity','Cytolytic'),
    ('Antigen_presentation','MHC-I'),
    ('MHC_II','MHC-II'),
    ('IFNg_Ayers_18','IFN-γ GEP'),
    ('TLS_Cabrita','TLS'),
    ('B_cell','B cell'),
    ('NK_cell','NK cell'),
    ('Treg','Treg'),
    ('TGFb_Mariathasan','TGF-β'),
    ('EMT_Mak','EMT'),
]
sig_keys = [s[0] for s in radar_sigs]
sig_labels = [s[1] for s in radar_sigs]
N = len(radar_sigs)
theta = np.linspace(0, 2*np.pi, N, endpoint=False)

# Compute medians per response group (pre-treatment only)
pre = sigs_m[sigs_m.timepoint=='pre']
good_med = pre[pre.response_bin=='good'][sig_keys].median().values
bad_med = pre[pre.response_bin=='bad'][sig_keys].median().values

# Closed polygon
theta_c = np.append(theta, theta[0])
good_c = np.append(good_med, good_med[0])
bad_c = np.append(bad_med, bad_med[0])

ax.fill(theta_c, good_c, color=GOOD, alpha=0.32, label='Good responders (median)')
ax.plot(theta_c, good_c, color=GOOD, lw=2.5, marker='o', markersize=7, markeredgecolor='white', markeredgewidth=0.8)
ax.fill(theta_c, bad_c, color=BAD, alpha=0.32, label='Poor responders (median)')
ax.plot(theta_c, bad_c, color=BAD, lw=2.5, marker='o', markersize=7, markeredgecolor='white', markeredgewidth=0.8)

# Reference circle at z=0
ax.plot(np.linspace(0, 2*np.pi, 200), [0]*200, color='#1d3557', lw=1, ls='--', alpha=0.4)

# Significance stars on labels
labels_with_sig = []
for sk, sl in radar_sigs:
    g = pre[pre.response_bin=='good'][sk].dropna()
    b = pre[pre.response_bin=='bad'][sk].dropna()
    p = stats.mannwhitneyu(g,b).pvalue if len(g)>=3 and len(b)>=3 else 1
    star = sig_symbol(p)
    if star and star != 'ns':
        labels_with_sig.append(f'{sl}\n{star}')
    else:
        labels_with_sig.append(sl)

ax.set_xticks(theta)
ax.set_xticklabels(labels_with_sig, fontsize=10)
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticklabels(['', '', '0', '+0.5', '+1'], fontsize=8, color='#6c757d')
ax.set_ylim(min(good_med.min(), bad_med.min())-0.15, max(good_med.max(), bad_med.max())+0.25)
ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
ax.spines['polar'].set_color('#1d3557')
ax.grid(color='#dfe5ec', lw=0.6)
ax.legend(loc='upper right', bbox_to_anchor=(1.32, 1.08), fontsize=9, frameon=False)
ax.set_title('Tumor microenvironment radar — 12 immune & stromal signatures\n(median z-score per response group, pre-treatment)\n* p<0.05, † p<0.1',
             fontsize=11, fontweight='bold', y=1.10, color='#1d3557')
save_panel(fig, 'Fig3A_TME_radar', OUT)

# ============================================================
# FIG 3B — Mariathasan-style multi-tier annotated signature heatmap
# ============================================================
sig_cols = list(SIG_CATEGORY.keys())
sig_cols = [c for c in sig_cols if c in sigs_m.columns]
pre = sigs_m[sigs_m.timepoint=='pre']
mat = pre.set_index('sample_id')[sig_cols]
mat_z = mat.sub(mat.mean()).div(mat.std())

# Order samples: by response then hierarchically within
from scipy.cluster.hierarchy import linkage, leaves_list
order_resp = pre.sort_values('response_bin').sample_id.tolist()
mat_z = mat_z.loc[order_resp]

# Order signatures by category
sig_cols_sorted = sorted(sig_cols, key=lambda s: (SIG_CATEGORY.get(s,'Other'), s))

# Annotation columns
def s2subj(s): return rna_inv[rna_inv.sample_id==s].subject_id.iloc[0]
resp_a = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_resp], index=order_resp)
ct_a = pd.Series([clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in order_resp], index=order_resp)
sex_a = pd.Series([clin[clin.subject_id==s2subj(s)].sex.iloc[0] for s in order_resp], index=order_resp)
age_a = pd.Series([clin[clin.subject_id==s2subj(s)].age.iloc[0] for s in order_resp], index=order_resp)

# TMB
tmb_map = tmb.set_index('sample_id')['TMB_nonsyn_per_Mb']
tmb_a = pd.Series([tmb_map.get(f'{s2subj(s)}-PR', np.nan) for s in order_resp], index=order_resp)

# Build column color DataFrame
col_colors_data = pd.DataFrame({
    'Response': resp_a.map(PAL_RESP),
    'cT': ct_a.map(PAL_STAGE),
    'Sex': sex_a.map({'M':'#264653','F':'#e63946'}),
    'Age': pd.Series([plt.cm.YlOrRd((a-30)/50) for a in age_a], index=order_resp).apply(matplotlib.colors.to_hex),
    'TMB': pd.Series([plt.cm.Purples(min(t,3)/3) if not pd.isna(t) else '#ecf0f1' for t in tmb_a], index=order_resp).apply(matplotlib.colors.to_hex),
})

# Row category colors
row_colors_data = pd.Series({s: CAT_COLORS[SIG_CATEGORY.get(s,'Other')] for s in sig_cols_sorted}, name='Module')

cg = sns.clustermap(mat_z[sig_cols_sorted].T, row_cluster=False, col_cluster=False,
                    cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
                    col_colors=col_colors_data, row_colors=row_colors_data,
                    figsize=(13.5, 8), xticklabels=False, yticklabels=True,
                    dendrogram_ratio=(0.06, 0.04), cbar_kws={'label':'z-score','shrink':0.4},
                    linewidths=0.2, linecolor='#fafafa')
cg.ax_heatmap.set_xlabel('Pre-treatment samples (n = 33; sorted by response)', fontsize=10)
cg.ax_heatmap.set_ylabel('')
cg.ax_heatmap.set_yticklabels([t.get_text().replace('_',' ') for t in cg.ax_heatmap.get_yticklabels()], fontsize=9)

# Module legend
import matplotlib.lines as mlines
mod_handles = [mpatches.Patch(color=c, label=k) for k, c in CAT_COLORS.items()]
cg.fig.legend(handles=mod_handles, loc='center left', bbox_to_anchor=(0.95, 0.5),
              fontsize=8.5, title='Signature module', title_fontsize=9.5, frameon=False)
cg.fig.suptitle('Immune & stromal signature landscape  (signatures grouped by functional module)',
                y=1.01, fontsize=12, fontweight='bold', color='#1d3557')
cg.fig.savefig(OUT/'Fig3B_signature_heatmap.pdf', bbox_inches='tight')
cg.fig.savefig(OUT/'Fig3B_signature_heatmap.png', dpi=400, bbox_inches='tight', facecolor='white')
plt.close('all')
print('  ✓ Fig3B_signature_heatmap')

# ============================================================
# FIG 3C — Annotated volcano (Jerby-Arnon Cell 2018 style)
# ============================================================
# Pathway-categorized gene coloring
CC_GENES = {'MKI67','TOP2A','CCNB1','CCNB2','CDK1','CDC20','CDC25A','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7',
            'BIRC5','CENPF','PLK1','AURKA','AURKB','BUB1','TYMS','UBE2C','CCNE1','PCNA','RRM2'}
DR_GENES = {'BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR','CHEK1','CHEK2','MRE11',
            'RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2','FANCI','FANCL','BLM','BRIP1','EXO1','POLD1','POLE',
            'MSH2','MSH6','MLH1','PMS2','OGG1','MUTYH'}
EMT_GENES = {'VIM','CDH2','SNAI1','SNAI2','TWIST1','TWIST2','FN1','MMP2','MMP3','MMP9','MMP14','ZEB1','ZEB2',
             'TGFB1','TGFB2','TGFBR1','COL1A1','COL3A1','COL4A1','FAP','ACTA2','S100A4','POSTN','SPARC'}
IM_GENES = {'CD8A','CD8B','GZMA','GZMB','GZMK','GZMH','PRF1','GNLY','IFNG','CD3D','CD3E','HLA-A','HLA-B','HLA-C',
            'B2M','NLRC5','TAP1','TAP2','CIITA','HLA-DRA','HLA-DRB1','MS4A1','CD79A','CD79B','BANK1','FOXP3',
            'CTLA4','PDCD1','LAG3','HAVCR2','TIGIT','CD274','CXCL9','CXCL10','CCL19','CCL21','CXCL13'}
def cat(g):
    if g in CC_GENES: return 'Cell cycle'
    if g in DR_GENES: return 'DNA repair'
    if g in IM_GENES: return 'Immune'
    if g in EMT_GENES: return 'EMT/Stromal'
    return 'Other'
deg_p = deg.copy()
deg_p['cat'] = deg_p.gene.apply(cat)
deg_p['-log10p'] = -np.log10(deg_p.pvalue.replace(0, 1e-300))

volcano_colors = {'Cell cycle':'#0f8b78','DNA repair':'#118ab2','Immune':'#ef476f',
                  'EMT/Stromal':'#e76f51','Other':'#cccccc'}

fig, ax = plt.subplots(figsize=(8, 6.5))
# Background
ns = deg_p[deg_p['cat']=='Other']
ax.hexbin(ns.log2FoldChange.clip(-5,5), ns['-log10p'].clip(0,8), gridsize=50, mincnt=2,
          cmap=LinearSegmentedColormap.from_list('h',['#ffffff','#dde2e8']), zorder=1)

# Colored category points
for c in ['Cell cycle','DNA repair','Immune','EMT/Stromal']:
    sub = deg_p[deg_p['cat']==c]
    ax.scatter(sub.log2FoldChange.clip(-5,5), sub['-log10p'].clip(0,8),
               s=22, alpha=0.85, color=volcano_colors[c],
               edgecolor='#1d3557', linewidth=0.3, label=c, zorder=3)
# Other but significant
others_sig = deg_p[(deg_p['cat']=='Other') & (deg_p.pvalue<1e-3)]
ax.scatter(others_sig.log2FoldChange.clip(-5,5), others_sig['-log10p'].clip(0,8),
           s=14, alpha=0.6, color='#6c757d', edgecolor='none', zorder=2)

# Label key categorized genes
to_label = deg_p[(deg_p['cat']!='Other') & (deg_p['-log10p']>2.5)].head(35)
texts = []
for _, r in to_label.iterrows():
    t = ax.text(np.clip(r.log2FoldChange, -5, 5), np.clip(r['-log10p'], 0, 8), r.gene,
                fontsize=8, color=volcano_colors[r['cat']], fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#6c757d', lw=0.3),
                expand_points=(1.3,1.4), force_points=0.4)
except: pass

# Threshold lines
ax.axvline(0, color='#1d3557', lw=0.5, alpha=0.4)
ax.axhline(-np.log10(0.01), color='#6c757d', lw=0.6, ls='--', alpha=0.5)
ax.text(4.7, -np.log10(0.01)+0.1, 'p = 0.01', fontsize=8, color='#6c757d', ha='right')

# Direction arrows
ax.annotate('', xy=(4.5, -0.5), xytext=(0.5, -0.5), xycoords=('data','axes fraction'),
            arrowprops=dict(arrowstyle='->', color=GOOD, lw=2))
ax.text(2.3, -0.62, 'Up in Good', ha='center', fontsize=9.5, color=GOOD, fontweight='bold',
        transform=ax.get_xaxis_transform())
ax.annotate('', xy=(-4.5, -0.5), xytext=(-0.5, -0.5), xycoords=('data','axes fraction'),
            arrowprops=dict(arrowstyle='->', color=BAD, lw=2))
ax.text(-2.3, -0.62, 'Up in Poor', ha='center', fontsize=9.5, color=BAD, fontweight='bold',
        transform=ax.get_xaxis_transform())

ax.set_xlim(-5.2, 5.2); ax.set_ylim(0, 8.2)
ax.set_xlabel('log2 fold change (good vs poor)')
ax.set_ylabel('−log10(p-value)')
ax.set_title('Pre-treatment differential expression (DESeq2)\nlabeled by pathway category', pad=10)
ax.legend(loc='upper left', fontsize=9, title='Category', title_fontsize=10, frameon=False)
add_axis_spines(ax)
save_panel(fig, 'Fig3C_volcano_journal', OUT)

# ============================================================
# FIG 3D — Effect-size forest (Litchfield Cell 2021 style)
# ============================================================
fig, ax = plt.subplots(figsize=(8, 7.5))
pre_st = sig_stats[sig_stats.timepoint=='pre'].sort_values('pvalue').head(20).iloc[::-1].reset_index(drop=True)

# Bootstrap CI (per signature)
np.random.seed(42)
def bootstrap_delta(g_vals, b_vals, n=1000):
    diffs = []
    for _ in range(n):
        gi = np.random.choice(g_vals, len(g_vals), replace=True)
        bi = np.random.choice(b_vals, len(b_vals), replace=True)
        diffs.append(np.mean(gi) - np.mean(bi))
    return np.percentile(diffs, [2.5, 97.5])

cis = []
for _, r in pre_st.iterrows():
    sig_n = r.signature
    g_vals = sigs_m[(sigs_m.timepoint=='pre') & (sigs_m.response_bin=='good')][sig_n].dropna().values
    b_vals = sigs_m[(sigs_m.timepoint=='pre') & (sigs_m.response_bin=='bad')][sig_n].dropna().values
    if len(g_vals)>=3 and len(b_vals)>=3:
        ci = bootstrap_delta(g_vals, b_vals)
    else:
        ci = (0, 0)
    cis.append(ci)
pre_st['ci_low'] = [c[0] for c in cis]
pre_st['ci_high'] = [c[1] for c in cis]

y = np.arange(len(pre_st))
for i, (_, r) in enumerate(pre_st.iterrows()):
    color = GOOD if r.delta_good_minus_bad>0 else BAD
    # CI line
    ax.plot([r.ci_low, r.ci_high], [i, i], color=color, lw=1.6, alpha=0.7)
    # Center dot
    ax.scatter(r.delta_good_minus_bad, i, s=160 if r.pvalue<0.05 else 90,
               color=color, edgecolor='#1d3557', linewidth=1.0, zorder=3)
    # P-value label
    star = sig_symbol(r.pvalue)
    if star == 'ns': star = ''
    label = f'p = {r.pvalue:.3g} {star}'
    ax.text(max(r.ci_high, r.delta_good_minus_bad) + 0.05, i, label,
            va='center', ha='left', fontsize=8, color='#1d3557')

ax.axvline(0, color='#1d3557', lw=0.9)
ax.axvspan(-0.05, 0.05, color='#dee2e6', alpha=0.4)
ax.set_yticks(y)
ax.set_yticklabels([s.replace('_',' ') for s in pre_st.signature], fontsize=9)
ax.set_xlabel('Δ z-score  (good − poor)\n95% bootstrap CI')
ax.set_title('Forest plot of signature effect sizes  (pre-treatment)', pad=10)
add_axis_spines(ax)
ax.set_xlim(min(pre_st.ci_low.min(), -0.5)*1.05, max(pre_st.ci_high.max(), 1.5)*1.5)
save_panel(fig, 'Fig3D_forest_lollipop', OUT)

# ============================================================
# FIG 3E — Bivariate CD8 functional state space (Tirosh/Sade-Feldman)
# ============================================================
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(7.5, 7.5))
gs = gridspec.GridSpec(4, 4, hspace=0.05, wspace=0.05)
ax_main = fig.add_subplot(gs[1:, :3])
ax_top = fig.add_subplot(gs[0, :3], sharex=ax_main)
ax_right = fig.add_subplot(gs[1:, 3], sharey=ax_main)

# X = CD8 cytolytic, Y = CD8 exhaustion  (classic Tirosh axes)
# Use Cytolytic_activity vs CD8_exhaustion if both significant; else proliferation vs exhaustion
x_sig = 'Cytolytic_activity'; y_sig = 'CD8_exhaustion'
pre = sigs_m[sigs_m.timepoint=='pre']

for resp in ['good','bad']:
    sub = pre[pre.response_bin==resp]
    ax_main.scatter(sub[x_sig], sub[y_sig], color=PAL_RESP[resp], s=110, alpha=0.85,
                    edgecolor='white', linewidth=1.2, label=f'{resp} (n={len(sub)})', zorder=4)
    # Density contour
    try:
        sns.kdeplot(x=sub[x_sig], y=sub[y_sig], ax=ax_main, color=PAL_RESP[resp],
                    levels=3, alpha=0.5, lw=1.3)
    except: pass
    # Marginals
    sns.kdeplot(x=sub[x_sig].dropna(), ax=ax_top, color=PAL_RESP[resp], fill=True, alpha=0.4, lw=1.5)
    sns.kdeplot(y=sub[y_sig].dropna(), ax=ax_right, color=PAL_RESP[resp], fill=True, alpha=0.4, lw=1.5)

# Quadrant annotations (Tirosh-style functional zones)
xmid = pre[x_sig].median(); ymid = pre[y_sig].median()
ax_main.axvline(xmid, color='#6c757d', ls=':', lw=0.6, alpha=0.7)
ax_main.axhline(ymid, color='#6c757d', ls=':', lw=0.6, alpha=0.7)
xmax, ymax = pre[x_sig].max(), pre[y_sig].max()
xmin, ymin = pre[x_sig].min(), pre[y_sig].min()
quad_labels = [
    (xmax*0.85, ymax*0.85, 'Active +\nExhausted', '#1d3557'),
    (xmin*0.85, ymax*0.85, 'Cold +\nExhausted', '#1d3557'),
    (xmax*0.85, ymin*0.85, 'Functional\ncytotoxic', GOOD),
    (xmin*0.85, ymin*0.85, 'Cold tumor', BAD),
]
for x, y, lbl, color in quad_labels:
    ax_main.text(x, y, lbl, ha='center', va='center', fontsize=8.5, color=color,
                 fontweight='bold', alpha=0.7,
                 bbox=dict(facecolor='white', edgecolor=color, alpha=0.8, boxstyle='round,pad=0.3'))

ax_main.set_xlabel('CD8 cytolytic activity (z-score)')
ax_main.set_ylabel('CD8 exhaustion (z-score)')
ax_main.legend(loc='lower left', fontsize=9, frameon=False)

ax_top.set_xticks([]); ax_top.set_yticks([])
ax_right.set_xticks([]); ax_right.set_yticks([])
for s in ['top','right','bottom']: ax_top.spines[s].set_visible(False)
for s in ['top','right','left']: ax_right.spines[s].set_visible(False)
add_axis_spines(ax_main)
fig.suptitle('CD8 functional state space (pre-treatment)\n(Tirosh-style cytolytic × exhaustion biaxial)',
             fontsize=11.5, fontweight='bold', y=0.97, color='#1d3557')
save_panel(fig, 'Fig3E_CD8_biaxial', OUT)

# ============================================================
# FIG 3F — TLS deep dive (Cabrita Nature 2020 style)
# ============================================================
TLS_GENES_CABRITA = ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3','CD79B','MS4A1','CCL18','PTGDS','CXCL8']
present_tls = [g for g in TLS_GENES_CABRITA if g in log_tpm.index]

# Pre samples ordered by TLS_Cabrita score
pre = sigs_m[sigs_m.timepoint=='pre'].sort_values(['response_bin','TLS_Cabrita'], ascending=[True, False])
pre_samples = pre.sample_id.tolist()
tls_mat = log_tpm.loc[present_tls, pre_samples]
# z-score per gene
tls_z = tls_mat.sub(tls_mat.mean(axis=1), axis=0).div(tls_mat.std(axis=1), axis=0)

fig = plt.figure(figsize=(13, 6.5))
gs = fig.add_gridspec(4, 1, height_ratios=[0.55, 0.18, 0.18, 4.5], hspace=0.06)

# Top: TLS signature score bar
ax_top = fig.add_subplot(gs[0])
tls_sig_vals = pre.TLS_Cabrita.values
colors = [PAL_RESP[r] for r in pre.response_bin]
ax_top.bar(range(len(pre_samples)), tls_sig_vals, color=colors, edgecolor='white', linewidth=0.5, width=0.92)
ax_top.axhline(0, color='#1d3557', lw=0.6)
ax_top.set_ylabel('TLS signature\nz-score', fontsize=9)
ax_top.set_xticks([]); ax_top.set_xlim(-0.5, len(pre_samples)-0.5)
add_axis_spines(ax_top)
ax_top.tick_params(labelsize=8)

# Annotation strips
ax_a1 = fig.add_subplot(gs[1])
arr = np.array([[matplotlib.colors.to_rgb(PAL_RESP[r]) for r in pre.response_bin]])
ax_a1.imshow(arr, aspect='auto', extent=[-0.5, len(pre_samples)-0.5, -0.5, 0.5])
ax_a1.set_yticks([0]); ax_a1.set_yticklabels(['Response'], fontsize=9)
ax_a1.set_xticks([]); ax_a1.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_a1.spines[s].set_visible(False)

ax_a2 = fig.add_subplot(gs[2])
mhc2_vals = pre.MHC_II.values
v_min, v_max = pre[['MHC_II','TLS_Cabrita']].values.min(), pre[['MHC_II','TLS_Cabrita']].values.max()
norm = (mhc2_vals - v_min)/(v_max-v_min+1e-9)
arr2 = plt.cm.Purples(norm.clip(0,1))[:,:3][np.newaxis,...]
ax_a2.imshow(arr2, aspect='auto', extent=[-0.5, len(pre_samples)-0.5, -0.5, 0.5])
ax_a2.set_yticks([0]); ax_a2.set_yticklabels(['MHC-II'], fontsize=9)
ax_a2.set_xticks([]); ax_a2.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_a2.spines[s].set_visible(False)

# Main heatmap
ax_main = fig.add_subplot(gs[3])
im = ax_main.imshow(tls_z.values, cmap='RdBu_r', aspect='auto', vmin=-2.5, vmax=2.5,
                     interpolation='nearest', extent=[-0.5, len(pre_samples)-0.5, -0.5, len(present_tls)-0.5])
ax_main.set_yticks(range(len(present_tls)))
ax_main.set_yticklabels(present_tls[::-1], fontsize=10, fontstyle='italic')
ax_main.invert_yaxis()
ax_main.set_xticks(range(len(pre_samples)))
ax_main.set_xticklabels(pre_samples, fontsize=6.5, rotation=90)
ax_main.tick_params(length=2)

cbar = fig.colorbar(im, ax=ax_main, shrink=0.5, pad=0.01, fraction=0.025)
cbar.set_label('log2(TPM+1) z-score', fontsize=9)
ax_main.set_xlabel('Pre-treatment samples (sorted by TLS signature within response)')

# Stats annotation
g = pre[pre.response_bin=='good'].TLS_Cabrita
b = pre[pre.response_bin=='bad'].TLS_Cabrita
p_tls = stats.mannwhitneyu(g, b).pvalue
fig.text(0.93, 0.94, f'TLS sig: good vs bad\np = {p_tls:.3f}', ha='right', va='top',
         fontsize=10, color='#1d3557', fontweight='bold',
         bbox=dict(facecolor='white', edgecolor='#1d3557', alpha=0.85, boxstyle='round,pad=0.5'))

fig.suptitle('Tertiary lymphoid structure (TLS) gene panel — Cabrita 12-gene signature',
             fontsize=12, fontweight='bold', y=0.97, color='#1d3557')
save_panel(fig, 'Fig3F_TLS_Cabrita', OUT)

print('\n=== Fig 3 v3 (6 journal-style panels) complete ===')
print(f'Output: {OUT}')
