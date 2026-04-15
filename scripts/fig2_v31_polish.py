"""
Fig 2 v3.1 — deeper/saturated colors + no titles.
Panels preserve motif grounding (Bailey Cell 2018 / Alexandrov Nature 2020 / Cercek NEJM 2022 / TCGA / Knijnenburg Cell 2018).
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# Saturated palette (consistent with Fig 3-6 v3)
GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'; GOLD_DEEP = '#d4a300'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}
PAL_STAGE_DEEP = {'T2':'#7fb0c4', 'T2/T3':'#1c5d7e', 'T3':'#0e2a47', 'T4':'#a01b2b'}
# Saturated mutation type colors (stronger versions)
PAL_MUT_DEEP = {
    'Missense': '#057a64',
    'Nonsense': '#c01a1a',
    'Frameshift': '#0e2a47',
    'In-frame indel': '#7a3aad',
    'Splice site': '#e05a0e',
    'Other protein-altering': '#6a0572',
    'Stop lost': '#a01b2b',
    'Start lost': '#9e240e',
    'Multiple': '#1d3557',
}
# Pathway colors deep
PATH_COLORS_DEEP = {
    'WNT':'#00567d','TP53':'#c96806','RTK/RAS':'#c01a1a','PI3K':'#7a3aad','TGFβ':'#6a4c24',
    'MMR':'#057a64','DDR':'#ad2974','Chromatin':'#4a5565','Cell cycle':'#b8a000','Other':'#0099b8',
}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(ROOT/'02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
variants = pd.read_csv(ROOT/'02_wes_tmb_msi/variant_master.tsv.gz', sep='\t', low_memory=False)
sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_activities_with_meta.tsv', sep='\t')
cnv = pd.read_csv(ROOT/'04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
hrd = pd.read_csv(ROOT/'04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
sbs96 = pd.read_csv(ROOT/'01_wes_signatures/vcf_input/output/SBS/Input_vcffiles.SBS96.all', sep='\t', index_col=0)

UNMATCHED = [13,15,16,17,18,19,33]

# Override raincloud palette
def raincloud_deep(ax, data, x, y, order, palette, width=0.4, alpha=0.75):
    for i, g in enumerate(order):
        vals = pd.to_numeric(data[data[x]==g][y], errors='coerce').dropna().values
        if len(vals) < 2: continue
        color = palette[g]
        parts = ax.violinplot(vals, positions=[i-0.15], widths=width, showmeans=False,
                              showextrema=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(color); pc.set_alpha(0.55); pc.set_edgecolor(color); pc.set_linewidth(0.9)
            m = np.mean(pc.get_paths()[0].vertices[:,0])
            pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], -np.inf, m)
        bp = ax.boxplot([vals], positions=[i+0.08], widths=0.15, patch_artist=True, showfliers=False,
                        medianprops=dict(color='white', linewidth=1.6),
                        boxprops=dict(facecolor=color, alpha=0.95, edgecolor=color, linewidth=0.9),
                        whiskerprops=dict(color=color, linewidth=0.9),
                        capprops=dict(color=color, linewidth=0.9))
        jitter_x = i + 0.28 + np.random.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(jitter_x, vals, s=18, color=color, alpha=alpha, edgecolor='white', linewidth=0.5, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=11, color='#0e2a47')

# ============================================================
# 2A — Oncoprint (Bailey Cell 2018 style, deep colors, no title)
# ============================================================
CRC_DRIVERS = ['APC','TP53','KRAS','PIK3CA','SMAD4','FBXW7','BRAF','NRAS','TCF7L2','KMT2D',
               'KMT2C','ARID1A','SOX9','CTNNB1','AMER1','RNF43','MSH6','MLH1','MSH2','PMS2',
               'POLE','POLD1','ATM','ERBB2','ERBB3','PTEN','CDKN2A','CCND1','MYC','GNAS']
NONSYN_TYPES = {
    'missense_variant': 'Missense', 'stop_gained': 'Nonsense',
    'frameshift_variant': 'Frameshift',
    'inframe_insertion': 'In-frame indel', 'inframe_deletion': 'In-frame indel',
    'splice_acceptor_variant': 'Splice site', 'splice_donor_variant': 'Splice site',
    'protein_altering_variant': 'Other protein-altering',
    'stop_lost': 'Stop lost', 'start_lost': 'Start lost',
}

v_drv = variants[variants.GENE.isin(CRC_DRIVERS) & variants.is_nonsyn].copy()
v_drv['mut_type'] = v_drv.EFFECT_primary.map(NONSYN_TYPES).fillna('Other protein-altering')
pre_samples = wes_inv[wes_inv.timepoint=='pre'].sample_id.tolist()
v_drv_pre = v_drv[v_drv.sample_id.isin(pre_samples)]
gene_carriers = v_drv_pre.groupby('GENE').sample_id.nunique().sort_values(ascending=False)
top_genes = [g for g in gene_carriers.index if gene_carriers[g] >= 1][:25]

def s2subj(s): return int(s.split('-')[0])
sample_meta = []
for s in pre_samples:
    n_drv = v_drv_pre[v_drv_pre.sample_id==s].GENE.nunique()
    resp = wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]
    sample_meta.append({'sample':s,'response':resp,'n_drv':n_drv,'subject':s2subj(s)})
sample_meta = pd.DataFrame(sample_meta).sort_values(['response','n_drv'], ascending=[True, False])
sample_order = sample_meta['sample'].tolist()

tmb_by_type = variants[(variants.is_coding) & variants.sample_id.isin(sample_order)].copy()
tmb_by_type['mt'] = tmb_by_type.EFFECT_primary.map(NONSYN_TYPES).fillna('Synonymous')
tmb_stack = tmb_by_type.groupby(['sample_id','mt']).size().unstack(fill_value=0)
tmb_stack = tmb_stack.reindex(sample_order).fillna(0)

n_samples = len(sample_order); n_genes = len(top_genes)
fig = plt.figure(figsize=(16, 11.5))
gs = fig.add_gridspec(8, 3,
    height_ratios=[1.5, 0.18, 0.18, 0.18, 0.18, n_genes*0.32, 0.5, 0.6],
    width_ratios=[n_samples*0.22, 1.2, 2.6],
    hspace=0.06, wspace=0.04)

# Top TMB stacked bar
ax_tmb = fig.add_subplot(gs[0, 0])
type_order = ['Missense','Nonsense','Frameshift','In-frame indel','Splice site','Other protein-altering']
type_order_present = [t for t in type_order if t in tmb_stack.columns]
bottom = np.zeros(n_samples)
for t in type_order_present:
    vals = tmb_stack[t].values
    ax_tmb.bar(range(n_samples), vals, bottom=bottom, color=PAL_MUT_DEEP[t], width=0.92,
               edgecolor='white', linewidth=0.3, label=t)
    bottom += vals
ax_tmb.set_ylabel('Coding mutations\n(per sample)', fontsize=10, color='#0e2a47', fontweight='bold')
ax_tmb.set_xlim(-0.5, n_samples-0.5); ax_tmb.set_xticks([])
add_axis_spines(ax_tmb)
ax_tmb.tick_params(labelsize=8.5)

def ann_strip(ax, vals, label, palette_map=None, cmap=None, vmin=None, vmax=None):
    if palette_map:
        arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        v_n = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(v_n.replace([np.inf,-np.inf], np.nan).dropna())
        if vmax is None: vmax = np.nanmax(v_n.replace([np.inf,-np.inf], np.nan).dropna())
        norm = ((v_n-vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[-0.5, n_samples-0.5, -0.5, 0.5])
    ax.set_yticks([0]); ax.set_yticklabels([label], fontsize=10, color='#0e2a47', fontweight='bold')
    ax.set_xticks([]); ax.tick_params(length=0)
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(-0.5, n_samples-0.5)

resp_vals = [wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in sample_order]
ct_vals = [clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in sample_order]
matched_vals = [s2subj(s) not in UNMATCHED for s in sample_order]
tmb_vals = [tmb[tmb.sample_id==s].TMB_nonsyn_per_Mb.iloc[0] if (tmb.sample_id==s).any() else np.nan for s in sample_order]
msi_vals = [msi[msi.sample_id==s].MSI_pct.iloc[0] if (msi.sample_id==s).any() else np.nan for s in sample_order]

ax_a1 = fig.add_subplot(gs[1,0]); ann_strip(ax_a1, resp_vals, 'Response', PAL)
ax_a2 = fig.add_subplot(gs[2,0]); ann_strip(ax_a2, ct_vals, 'cT stage', PAL_STAGE_DEEP)
ax_a3 = fig.add_subplot(gs[3,0]); ann_strip(ax_a3, tmb_vals, 'TMB/Mb', cmap=plt.cm.Purples, vmin=0, vmax=3)
ax_a4 = fig.add_subplot(gs[4,0]); ann_strip(ax_a4, msi_vals, 'MSI %', cmap=plt.cm.Greens, vmin=0, vmax=0.3)

# Main oncoprint
ax_onco = fig.add_subplot(gs[5, 0])
mut_lookup = {}
for _, r in v_drv_pre[v_drv_pre.sample_id.isin(sample_order) & v_drv_pre.GENE.isin(top_genes)].iterrows():
    mut_lookup.setdefault((r.sample_id, r.GENE), []).append(r.mut_type)

for i in range(n_genes):
    for j in range(n_samples):
        ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06), 0.92, 0.88,
                                    facecolor='#eef2f5', edgecolor='white', linewidth=0.3))
for (sid, gene), types in mut_lookup.items():
    j = sample_order.index(sid); i = top_genes.index(gene)
    types_unique = list(dict.fromkeys(types))
    if len(types_unique) > 1:
        h_each = 0.88 / len(types_unique)
        for k, t in enumerate(types_unique):
            color = PAL_MUT_DEEP.get(t, '#5a6772')
            ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06+k*h_each), 0.92, h_each,
                                        facecolor=color, edgecolor='white', linewidth=0.3))
    else:
        color = PAL_MUT_DEEP.get(types_unique[0], '#5a6772')
        ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06), 0.92, 0.88,
                                    facecolor=color, edgecolor='white', linewidth=0.3))

ax_onco.set_xlim(-0.05, n_samples); ax_onco.set_ylim(-0.05, n_genes)
ax_onco.set_yticks([n_genes-1-i+0.5 for i in range(n_genes)])
ax_onco.set_yticklabels(top_genes, fontsize=10, fontstyle='italic', color='#0e2a47')
ax_onco.set_xticks([]); ax_onco.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_onco.spines[s].set_visible(False)

# Right marginal: frequency by response
ax_rmg = fig.add_subplot(gs[5, 1])
from scipy.stats import fisher_exact
n_good_total = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
n_bad_total = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
gene_freq = []
for g in top_genes:
    g_mut = sum(1 for s in sample_order if (s,g) in mut_lookup and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
    b_mut = sum(1 for s in sample_order if (s,g) in mut_lookup and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
    fp = fisher_exact([[g_mut, n_good_total-g_mut],[b_mut, n_bad_total-b_mut]]).pvalue
    gene_freq.append({'gene':g,'good_pct':g_mut/n_good_total*100,'bad_pct':b_mut/n_bad_total*100,'p':fp})
gf_df = pd.DataFrame(gene_freq)
y = np.arange(n_genes)[::-1]
bar_h = 0.35
ax_rmg.barh(y-bar_h/2, gf_df.good_pct, height=bar_h, color=GOOD_DEEP, edgecolor='white', linewidth=0.5, label='Good')
ax_rmg.barh(y+bar_h/2, gf_df.bad_pct, height=bar_h, color=BAD_DEEP, edgecolor='white', linewidth=0.5, label='Poor')
for i, r in gf_df.iterrows():
    if r.p < 0.1:
        x_pos = max(r.good_pct, r.bad_pct) + 2
        star = sig_symbol(r.p)
        ax_rmg.text(x_pos, n_genes-1-i, star, va='center', fontsize=11, color='#c01a1a', fontweight='bold')
ax_rmg.set_xlabel('% mutated', fontsize=10, color='#0e2a47', fontweight='bold')
ax_rmg.set_xlim(0, max(gf_df[['good_pct','bad_pct']].max().max()*1.15, 60))
ax_rmg.set_yticks([]); ax_rmg.set_ylim(-0.5, n_genes-0.5)
ax_rmg.legend(fontsize=9, loc='lower right', frameon=False)
add_axis_spines(ax_rmg)
ax_rmg.tick_params(labelsize=8.5)

# Right: pathway annotation
ax_cat = fig.add_subplot(gs[5, 2])
ax_cat.axis('off')
ax_cat.text(0.05, 0.99, 'Pathway / Function', transform=ax_cat.transAxes, fontsize=10.5,
            fontweight='bold', va='top', color='#0e2a47')
GENE_PATHWAY = {
    'APC':'WNT', 'CTNNB1':'WNT', 'AMER1':'WNT', 'RNF43':'WNT', 'TCF7L2':'WNT',
    'TP53':'TP53', 'CDKN2A':'TP53', 'ATM':'DDR',
    'KRAS':'RTK/RAS', 'BRAF':'RTK/RAS', 'NRAS':'RTK/RAS', 'ERBB2':'RTK/RAS', 'ERBB3':'RTK/RAS',
    'PIK3CA':'PI3K', 'PTEN':'PI3K', 'SMAD4':'TGFβ',
    'MLH1':'MMR', 'MSH2':'MMR', 'MSH6':'MMR', 'PMS2':'MMR',
    'POLE':'DDR', 'POLD1':'DDR',
    'KMT2D':'Chromatin', 'KMT2C':'Chromatin', 'ARID1A':'Chromatin',
    'FBXW7':'Cell cycle', 'CCND1':'Cell cycle', 'MYC':'Cell cycle',
    'SOX9':'Other', 'GNAS':'Other'
}
for i, gene in enumerate(top_genes):
    pw = GENE_PATHWAY.get(gene, 'Other')
    color = PATH_COLORS_DEEP[pw]
    ax_cat.add_patch(Rectangle((0.0, (n_genes-1-i+0.1)/n_genes), 0.04, 0.85/n_genes,
                                facecolor=color, edgecolor='white', linewidth=0.3, transform=ax_cat.transAxes))
    ax_cat.text(0.06, (n_genes-1-i+0.5)/n_genes, pw, transform=ax_cat.transAxes,
                fontsize=9, va='center', color='#0e2a47')

# Legend area
ax_leg1 = fig.add_subplot(gs[6, 0])
ax_leg1.axis('off')
mut_legend = [mpatches.Patch(color=PAL_MUT_DEEP[t], label=t) for t in type_order_present]
ax_leg1.legend(handles=mut_legend, loc='upper left', ncol=4, fontsize=9.5,
               title='Mutation type', title_fontsize=10, frameon=False)

ax_leg2 = fig.add_subplot(gs[7, 0])
ax_leg2.axis('off')
ann_items = []
for k, c in PAL.items(): ann_items.append(mpatches.Patch(color=c, label=k))
for k, c in PAL_STAGE_DEEP.items(): ann_items.append(mpatches.Patch(color=c, label=k))
ann_items.append(mpatches.Patch(color=plt.cm.Purples(0.7), label='TMB (gradient)'))
ann_items.append(mpatches.Patch(color=plt.cm.Greens(0.7), label='MSI (gradient)'))
ax_leg2.legend(handles=ann_items, loc='upper left', ncol=5, fontsize=9.5,
               title='Annotation tracks', title_fontsize=10, frameon=False)

save_panel(fig, 'Fig2A_oncoprint_journal', OUT)

# ============================================================
# 2B — TMB raincloud (deep, no title)
# ============================================================
fig, ax = plt.subplots(figsize=(4.8, 4.5))
pre_m = tmb[(tmb.timepoint=='pre') & (~tmb.subject_id.isin(UNMATCHED))]
raincloud_deep(ax, pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ['good','bad'], PAL)
ax.axhline(10, color='#c01a1a', ls='--', lw=1, alpha=0.85)
ax.text(1.5, 10, 'TMB-high (10/Mb)', fontsize=8.5, color='#c01a1a', va='bottom', ha='right', fontweight='bold')
g = pre_m[pre_m.response_bin=='good'].TMB_nonsyn_per_Mb
b = pre_m[pre_m.response_bin=='bad'].TMB_nonsyn_per_Mb
from scipy.stats import mannwhitneyu
p = mannwhitneyu(g, b).pvalue
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.3, p)
ax.set_ylabel('Nonsynonymous TMB (/Mb)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlabel('')
annotate_count(ax, pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ['good','bad'], dy=0.25)
save_panel(fig, 'Fig2B_TMB_raincloud', OUT)

# ============================================================
# 2C — MSI × TMB scatter (deep, no title)
# ============================================================
fig, ax = plt.subplots(figsize=(5.5, 4.8))
mdf = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id', how='left')
mdf = mdf.dropna(subset=['TMB_nonsyn_per_Mb','MSI_pct'])
for resp in ['good','bad']:
    sub = mdf[mdf.response_bin==resp]
    ax.scatter(sub.TMB_nonsyn_per_Mb, sub.MSI_pct, color=PAL[resp], s=85, alpha=0.85,
               edgecolor='white', linewidth=1.0, label=f'{resp} (n={len(sub)})', zorder=3)
ax.axhspan(20, 100, alpha=0.07, color=GOOD_DEEP)
ax.axvspan(10, 100, alpha=0.07, color=GOLD_DEEP)
ax.axhline(20, color=GOOD_DEEP, ls='--', lw=0.9, alpha=0.6)
ax.axvline(10, color=GOLD_DEEP, ls='--', lw=0.9, alpha=0.6)
ax.text(9, 0.01, 'MSS zone\n(TMB-low)', fontsize=9.5, color='#0e2a47', ha='right', va='bottom',
        style='italic', fontweight='bold')
ax.set_xlabel('Nonsynonymous TMB (/Mb)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('MSI percentage (%)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(-0.3, max(tmb.TMB_nonsyn_per_Mb.max()+1, 13))
ax.set_ylim(-0.02, 0.8)
ax.legend(fontsize=9.5, loc='upper right', frameon=False)
save_panel(fig, 'Fig2C_MSI_TMB_scatter', OUT)

# ============================================================
# 2D — MMR/SBS3 reassessment (deep, no title)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
sbs_pre = sbs[sbs.timepoint=='pre'].copy()
sbs_pre['total'] = sbs_pre[[c for c in sbs_pre.columns if c.startswith('SBS')]].sum(axis=1)
sbs_pre['MMR_prop'] = (sbs_pre.get('SBS6',0)+sbs_pre.get('SBS15',0)+sbs_pre.get('SBS20',0)+sbs_pre.get('SBS26',0))/sbs_pre.total.replace(0,np.nan)
sbs_pre['SBS3_prop'] = sbs_pre.get('SBS3',0)/sbs_pre.total.replace(0,np.nan)
CLAIMED_MMR = {1,9,12}; CLAIMED_SBS3 = {5,14}

ax = axes[0]
matched = sbs_pre[~sbs_pre.subject_id.isin(UNMATCHED)]
ax.scatter(matched[~matched.subject_id.isin(CLAIMED_MMR)].subject_id,
           matched[~matched.subject_id.isin(CLAIMED_MMR)].MMR_prop*100,
           color='#5a6772', s=68, alpha=0.82, edgecolor='white', linewidth=0.8)
highlighted = matched[matched.subject_id.isin(CLAIMED_MMR)]
for _, r in highlighted.iterrows():
    ax.scatter(r.subject_id, r.MMR_prop*100, s=230, color=GOLD_DEEP,
               edgecolor='#c01a1a', linewidth=2.2, marker='*', zorder=10)
    ax.annotate(f"S{int(r.subject_id)}\n(prior MMR claim)",
                (r.subject_id, r.MMR_prop*100), xytext=(5, 8),
                textcoords='offset points', fontsize=8.5, color='#c01a1a', fontweight='bold')
ax.axhline(20, color='#c01a1a', ls='--', lw=0.9, alpha=0.6)
ax.text(matched.subject_id.max(), 20.8, 'Clinical MMR-d cutoff', fontsize=8.5, color='#c01a1a',
        ha='right', fontweight='bold')
ax.set_xlabel('Subject ID (matched, pre)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('MMR signature proportion (%)', fontsize=11, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)

ax = axes[1]
ax.scatter(matched[~matched.subject_id.isin(CLAIMED_SBS3)].subject_id,
           matched[~matched.subject_id.isin(CLAIMED_SBS3)].SBS3_prop*100,
           color='#5a6772', s=68, alpha=0.82, edgecolor='white', linewidth=0.8)
highlighted3 = matched[matched.subject_id.isin(CLAIMED_SBS3)]
for _, r in highlighted3.iterrows():
    ax.scatter(r.subject_id, r.SBS3_prop*100, s=230, color=GOLD_DEEP,
               edgecolor='#c01a1a', linewidth=2.2, marker='*', zorder=10)
    ax.annotate(f"S{int(r.subject_id)}\n(prior SBS3 claim)",
                (r.subject_id, r.SBS3_prop*100), xytext=(5, 8),
                textcoords='offset points', fontsize=8.5, color='#c01a1a', fontweight='bold')
ax.set_xlabel('Subject ID (matched, pre)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('SBS3 (HRD) proportion (%)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylim(-0.5, 5)
add_axis_spines(ax)
save_panel(fig, 'Fig2D_SBS_reassessment', OUT)

# ============================================================
# 2E — CNV/HRD raincloud triple (deep, no title)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
cnv_m = cnv[(cnv.timepoint=='pre') & cnv.matched]
hrd_m = hrd[(hrd.timepoint=='pre') & hrd.matched]

ax = axes[0]
raincloud_deep(ax, cnv_m, 'response_bin', 'CIN', ['good','bad'], PAL)
g = cnv_m[cnv_m.response_bin=='good'].CIN; b = cnv_m[cnv_m.response_bin=='bad'].CIN
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.02, mannwhitneyu(g,b).pvalue)
ax.set_ylabel('Fraction genome altered (CIN)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlabel('')

ax = axes[1]
raincloud_deep(ax, hrd_m, 'response_bin', 'LST', ['good','bad'], PAL)
g = hrd_m[hrd_m.response_bin=='good'].LST; b = hrd_m[hrd_m.response_bin=='bad'].LST
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.5, mannwhitneyu(g,b).pvalue)
ax.set_ylabel('Large-scale transitions (LST)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlabel('')

ax = axes[2]
raincloud_deep(ax, hrd_m, 'response_bin', 'HRD_sum', ['good','bad'], PAL)
g = hrd_m[hrd_m.response_bin=='good'].HRD_sum; b = hrd_m[hrd_m.response_bin=='bad'].HRD_sum
stat_bracket(ax, 0, 1, max(g.max(), b.max())+1, mannwhitneyu(g,b).pvalue)
ax.set_ylabel('HRD total score (LST+LOH+TAI)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlabel('')

fig.tight_layout()
save_panel(fig, 'Fig2E_CNV_HRD', OUT)

print('\n=== Fig 2 v3.1 (deep colors, no titles) saved ===')
