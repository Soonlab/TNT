"""
Figure 2 — WES genomic landscape (comprehensive oncoprint)
  2A: Mega oncoprint (49 tumors × 25 drivers + 6 annotation tracks + SBS stacked bar)
  2B: TMB rain cloud (pre matched)
  2C: MSI vs TMB scatter
  2D: SBS signature refit (reassessment of prior MMR/SBS3 claims)
  2E: CNV/CIN + HRD proxy integrated
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import Rectangle

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(ROOT/'02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
variants = pd.read_csv(ROOT/'02_wes_tmb_msi/variant_master.tsv.gz', sep='\t', low_memory=False)
sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_activities_with_meta.tsv', sep='\t')
cnv = pd.read_csv(ROOT/'04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
hrd = pd.read_csv(ROOT/'04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv', sep='\t')

UNMATCHED = [13,15,16,17,18,19,33]

# ===========================================================
# 2A — Mega oncoprint
# ===========================================================
# CRC drivers
CRC_DRIVERS = ['APC','TP53','KRAS','BRAF','PIK3CA','SMAD4','FBXW7','KMT2D','ARID1A','SOX9',
               'AMER1','TCF7L2','NRAS','ATM','PIK3R1','CREBBP','ERBB3','POLE','POLD1',
               'MLH1','MSH2','MSH6','PMS2','CTNNB1','SMAD2']
# Filter non-synonymous
NONSYN = {'missense_variant','stop_gained','stop_lost','start_lost','frameshift_variant',
          'inframe_insertion','inframe_deletion','splice_acceptor_variant','splice_donor_variant',
          'protein_altering_variant'}
v_drv = variants[variants.GENE.isin(CRC_DRIVERS) & variants.is_nonsyn].copy()
# Per sample-gene, take most damaging effect (use primary)
gene_counts = v_drv.GENE.value_counts()
present_drivers = [g for g in CRC_DRIVERS if g in gene_counts.index and gene_counts[g] >= 1]
# Sort by count
present_drivers = sorted(present_drivers, key=lambda g: -gene_counts[g])[:20]

# Build mutation type matrix: rows = genes, cols = samples
pre_samples = wes_inv[wes_inv.timepoint=='pre'].sample_id.tolist()
# Order samples: good first, then bad; within each by subject id
sample_order = sorted(pre_samples, key=lambda s: (
    wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]!='good',
    wes_inv[wes_inv.sample_id==s].subject_id.iloc[0]
))
mut_mat = {}  # sample_id -> gene -> list of effects
for _, r in v_drv[v_drv.sample_id.isin(sample_order)].iterrows():
    mut_mat.setdefault(r.sample_id, {}).setdefault(r.GENE, []).append(r.EFFECT_primary)

# Figure layout
fig = plt.figure(figsize=(15, 11))
n_ann = 5
n_genes = len(present_drivers)
n_samples = len(sample_order)
gs = fig.add_gridspec(n_ann + 3, 2,
                     height_ratios=[0.28]*n_ann + [0.2, n_genes*0.45, 2.0],
                     width_ratios=[n_samples*0.25, 4],
                     hspace=0.1, wspace=0.05)

# Top annotation rows
def ann_bar(ax, vals, label, is_cat, palette_map=None, cmap=None, vmin=None, vmax=None):
    if is_cat:
        arr = np.array([[matplotlib.colors.to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        vals_num = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(vals_num.replace([np.inf,-np.inf], np.nan).dropna())
        if vmax is None: vmax = np.nanmax(vals_num.replace([np.inf,-np.inf], np.nan).dropna())
        norm = ((vals_num - vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest')
    ax.set_yticks([0]); ax.set_yticklabels([label], fontsize=9, color='#1d3557')
    ax.set_xticks([]); ax.tick_params(length=0)
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)

# Build annotation values in sample_order
def s2subj(s): return int(s.split('-')[0])
resp_vals = [wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in sample_order]
ct_vals = [clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in sample_order]
matched_vals = [s2subj(s) not in UNMATCHED for s in sample_order]
tmb_vals = [tmb[tmb.sample_id==s].TMB_nonsyn_per_Mb.iloc[0] if (tmb.sample_id==s).any() else np.nan for s in sample_order]
msi_vals = [msi[msi.sample_id==s].MSI_pct.iloc[0] if (msi.sample_id==s).any() else np.nan for s in sample_order]

ax_ann0 = fig.add_subplot(gs[0,0]); ann_bar(ax_ann0, resp_vals, 'Response', True, PAL_RESP)
ax_ann1 = fig.add_subplot(gs[1,0]); ann_bar(ax_ann1, ct_vals, 'cT stage', True, PAL_STAGE)
ax_ann2 = fig.add_subplot(gs[2,0]); ann_bar(ax_ann2, matched_vals, 'Matched N', True, {True:'#264653', False:'#e9c46a'})
ax_ann3 = fig.add_subplot(gs[3,0]); ann_bar(ax_ann3, tmb_vals, 'TMB/Mb', False, cmap=plt.cm.Purples, vmin=0, vmax=3)
ax_ann4 = fig.add_subplot(gs[4,0]); ann_bar(ax_ann4, msi_vals, 'MSI %', False, cmap=plt.cm.Greens, vmin=0, vmax=0.3)

# Main oncoprint
ax_onco = fig.add_subplot(gs[6,0])
for i, gene in enumerate(present_drivers):
    for j, sid in enumerate(sample_order):
        effects = mut_mat.get(sid, {}).get(gene, [])
        if not effects:
            ax_onco.add_patch(Rectangle((j, n_genes-1-i), 0.92, 0.88,
                              facecolor='#f5f5f5', edgecolor='white', linewidth=0.3))
        else:
            if len(set(effects)) > 1:
                color = PAL_MUT.get('multi')
            else:
                color = PAL_MUT.get(effects[0], '#8d99ae')
            ax_onco.add_patch(Rectangle((j, n_genes-1-i), 0.92, 0.88,
                              facecolor=color, edgecolor='white', linewidth=0.5))

ax_onco.set_xlim(-0.05, n_samples); ax_onco.set_ylim(-0.05, n_genes)
ax_onco.set_yticks([n_genes-1-i+0.44 for i in range(n_genes)])
ax_onco.set_yticklabels(present_drivers, fontsize=9)
ax_onco.set_xticks([])
ax_onco.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_onco.spines[s].set_visible(False)

# Right marginal: gene frequency by response
ax_rmg = fig.add_subplot(gs[6,1])
rows = []
for gene in present_drivers:
    good_carriers = sum(1 for s in sample_order
                        if mut_mat.get(s,{}).get(gene) and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
    bad_carriers = sum(1 for s in sample_order
                       if mut_mat.get(s,{}).get(gene) and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
    n_good = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
    n_bad = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
    rows.append({'gene':gene, 'good':good_carriers/n_good*100, 'bad':bad_carriers/n_bad*100})
freq = pd.DataFrame(rows)
y = np.arange(n_genes)[::-1]
ax_rmg.barh(y-0.2, freq.good, height=0.36, color=GOOD, edgecolor='white', label='Good')
ax_rmg.barh(y+0.2, freq.bad, height=0.36, color=BAD, edgecolor='white', label='Poor')
ax_rmg.set_xlabel('% mutated')
ax_rmg.set_yticks([]); ax_rmg.set_ylim(-0.5, n_genes-0.5)
ax_rmg.legend(fontsize=8, loc='lower right')
ax_rmg.spines['top'].set_visible(False); ax_rmg.spines['right'].set_visible(False)

# Bottom SBS stacked
ax_sbs = fig.add_subplot(gs[7,0])
sig_cols = [c for c in sbs.columns if c.startswith('SBS')]
active = [c for c in sig_cols if (sbs[c]>0).any()]
total_c = sbs[active].sum().sort_values(ascending=False)
keep_sigs = total_c.head(8).index.tolist()
sbs_pre = sbs[sbs.sample_id.isin(sample_order)].set_index('sample_id').reindex(sample_order)
sbs_norm = sbs_pre[keep_sigs].div(sbs_pre[keep_sigs].sum(axis=1).replace(0,np.nan), axis=0).fillna(0)
cmap = plt.cm.tab10(range(len(keep_sigs)))
bottom = np.zeros(n_samples)
for i, sig in enumerate(keep_sigs):
    ax_sbs.bar(range(n_samples), sbs_norm[sig].values, bottom=bottom, label=sig,
               color=cmap[i], width=0.92, edgecolor='white', linewidth=0.2)
    bottom += sbs_norm[sig].values
ax_sbs.set_ylabel('SBS proportion')
ax_sbs.set_xticks(range(n_samples))
ax_sbs.set_xticklabels(sample_order, rotation=90, fontsize=6.5)
ax_sbs.set_xlim(-0.5, n_samples-0.5); ax_sbs.set_ylim(0, 1.02)
ax_sbs.legend(title='SBS', loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=7.5, title_fontsize=8)
for s in ['top','right']: ax_sbs.spines[s].set_visible(False)

# Mutation type legend (bottom right marginal)
ax_leg = fig.add_subplot(gs[7,1])
ax_leg.axis('off')
legend_items = [mpatches.Patch(color=c, label=e.replace('_',' ').title())
                for e, c in PAL_MUT.items() if e != 'multi']
legend_items.append(mpatches.Patch(color=PAL_MUT['multi'], label='Multiple'))
ax_leg.legend(handles=legend_items, loc='upper left', ncol=1, fontsize=7.5,
              title='Mutation type', title_fontsize=8, frameon=False)

fig.suptitle('Genomic landscape: driver mutations, mutational signatures, and clinical/molecular context',
             fontsize=13, fontweight='bold', y=0.93, color='#1d3557')
save_panel(fig, 'Fig2A_oncoprint', OUT)

# ===========================================================
# 2B — TMB rain cloud (pre matched, good vs bad) + density
# ===========================================================
fig, ax = plt.subplots(figsize=(4.5, 4))
pre_m = tmb[(tmb.timepoint=='pre') & (~tmb.subject_id.isin(UNMATCHED))]
raincloud(ax, pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ['good','bad'], PAL_RESP)
# TMB-high reference
ax.axhline(10, color='#d62828', ls='--', lw=0.9, alpha=0.7)
ax.text(1.5, 10, 'TMB-high (10/Mb)', fontsize=8, color='#d62828', va='bottom', ha='right')
# Stat
g = pre_m[pre_m.response_bin=='good'].TMB_nonsyn_per_Mb
b = pre_m[pre_m.response_bin=='bad'].TMB_nonsyn_per_Mb
p = stats.mannwhitneyu(g, b).pvalue
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.3, p)
ax.set_ylabel('Nonsynonymous TMB (/Mb)')
ax.set_xlabel('')
ax.set_title('Tumor mutational burden\n(pre-treatment, matched)')
annotate_count(ax, pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ['good','bad'], dy=0.25)
save_panel(fig, 'Fig2B_TMB_raincloud', OUT)

# ===========================================================
# 2C — MSI × TMB scatter (proves MSS + low TMB)
# ===========================================================
fig, ax = plt.subplots(figsize=(5, 4.5))
mdf = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id', how='left')
mdf = mdf.dropna(subset=['TMB_nonsyn_per_Mb','MSI_pct'])
for resp in ['good','bad']:
    sub = mdf[mdf.response_bin==resp]
    ax.scatter(sub.TMB_nonsyn_per_Mb, sub.MSI_pct, color=PAL_RESP[resp], s=65, alpha=0.75,
               edgecolor='white', linewidth=0.8, label=f'{resp} (n={len(sub)})')
# Shaded zones
ax.axhspan(20, 100, alpha=0.06, color='#2a9d8f', label='MSI-H zone')
ax.axvspan(10, 100, alpha=0.06, color='#f4a261')
ax.axhline(20, color='#2a9d8f', ls='--', lw=0.8, alpha=0.5)
ax.axvline(10, color='#f4a261', ls='--', lw=0.8, alpha=0.5)
ax.text(9, 0.01, 'MSS zone\n(TMB-low)', fontsize=9, color='#6c757d', ha='right', va='bottom', style='italic')
ax.set_xlabel('Nonsynonymous TMB (/Mb)')
ax.set_ylabel('MSI percentage (%)')
ax.set_xlim(-0.3, max(tmb.TMB_nonsyn_per_Mb.max()+1, 13))
ax.set_ylim(-0.02, 0.8)  # zoom since all MSS
ax.legend(fontsize=8.5, loc='upper right')
ax.set_title('MSI × TMB: all samples fall in MSS / TMB-low zone')
save_panel(fig, 'Fig2C_MSI_TMB_scatter', OUT)

# ===========================================================
# 2D — SBS signature refit: prior vs Mutect2 (MMR + SBS3 reassessment)
# ===========================================================
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
sbs_pre = sbs[sbs.timepoint=='pre'].copy()
sbs_pre['total'] = sbs_pre[[c for c in sbs_pre.columns if c.startswith('SBS')]].sum(axis=1)
sbs_pre['MMR_prop'] = (sbs_pre.get('SBS6',0)+sbs_pre.get('SBS15',0)+sbs_pre.get('SBS20',0)+sbs_pre.get('SBS26',0))/sbs_pre.total.replace(0,np.nan)
sbs_pre['SBS3_prop'] = sbs_pre.get('SBS3',0)/sbs_pre.total.replace(0,np.nan)

CLAIMED_MMR = {1,9,12}; CLAIMED_SBS3 = {5,14}

# Panel 2D-left: MMR
ax = axes[0]
matched = sbs_pre[~sbs_pre.subject_id.isin(UNMATCHED)]
ax.scatter(matched[~matched.subject_id.isin(CLAIMED_MMR)].subject_id,
           matched[~matched.subject_id.isin(CLAIMED_MMR)].MMR_prop*100,
           color='#8d99ae', s=55, alpha=0.7, edgecolor='white', label='Other subjects')
highlighted = matched[matched.subject_id.isin(CLAIMED_MMR)]
for _, r in highlighted.iterrows():
    ax.scatter(r.subject_id, r.MMR_prop*100, s=180, color=HIGHLIGHT,
               edgecolor='#d62828', linewidth=2, marker='*', zorder=10)
    ax.annotate(f"S{int(r.subject_id)}\n(prior MMR claim)",
                (r.subject_id, r.MMR_prop*100), xytext=(5, 8),
                textcoords='offset points', fontsize=8, color='#d62828', fontweight='bold')
ax.axhline(20, color='#d62828', ls='--', lw=0.8, alpha=0.5)
ax.text(matched.subject_id.max(), 20.5, 'Clinical MMR-d cutoff', fontsize=7.5, color='#d62828', ha='right')
ax.set_xlabel('Subject ID (matched, pre)')
ax.set_ylabel('MMR signature proportion (%)')
ax.set_title('MMR signatures (SBS6/15/20/26)\nprior claim vs Mutect2 refit')

# Panel 2D-right: SBS3
ax = axes[1]
ax.scatter(matched[~matched.subject_id.isin(CLAIMED_SBS3)].subject_id,
           matched[~matched.subject_id.isin(CLAIMED_SBS3)].SBS3_prop*100,
           color='#8d99ae', s=55, alpha=0.7, edgecolor='white')
highlighted3 = matched[matched.subject_id.isin(CLAIMED_SBS3)]
for _, r in highlighted3.iterrows():
    ax.scatter(r.subject_id, r.SBS3_prop*100, s=180, color=HIGHLIGHT,
               edgecolor='#d62828', linewidth=2, marker='*', zorder=10)
    ax.annotate(f"S{int(r.subject_id)}\n(prior SBS3 claim)",
                (r.subject_id, r.SBS3_prop*100), xytext=(5, 8),
                textcoords='offset points', fontsize=8, color='#d62828', fontweight='bold')
ax.set_xlabel('Subject ID (matched, pre)')
ax.set_ylabel('SBS3 (HRD) proportion (%)')
ax.set_title('SBS3 signature\nall samples = 0 after Mutect2 refit')
ax.set_ylim(-0.5, 5)

save_panel(fig, 'Fig2D_SBS_reassessment', OUT)

# ===========================================================
# 2E — CNV / CIN + HRD proxy integrated view
# ===========================================================
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=False)

# CIN raincloud
ax = axes[0]
cnv_m = cnv[(cnv.timepoint=='pre') & cnv.matched]
raincloud(ax, cnv_m, 'response_bin', 'CIN', ['good','bad'], PAL_RESP)
g = cnv_m[cnv_m.response_bin=='good'].CIN; b = cnv_m[cnv_m.response_bin=='bad'].CIN
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.02, stats.mannwhitneyu(g,b).pvalue)
ax.set_ylabel('Fraction genome altered (CIN)')
ax.set_xlabel('')
ax.set_title('Chromosomal instability')

# LST
ax = axes[1]
hrd_m = hrd[(hrd.timepoint=='pre') & hrd.matched]
raincloud(ax, hrd_m, 'response_bin', 'LST', ['good','bad'], PAL_RESP)
g = hrd_m[hrd_m.response_bin=='good'].LST; b = hrd_m[hrd_m.response_bin=='bad'].LST
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.5, stats.mannwhitneyu(g,b).pvalue)
ax.set_ylabel('Large-scale transitions (LST)')
ax.set_xlabel('')
ax.set_title('HRD proxy: LST')

# HRD sum (LST+LOH+TAI)
ax = axes[2]
raincloud(ax, hrd_m, 'response_bin', 'HRD_sum', ['good','bad'], PAL_RESP)
g = hrd_m[hrd_m.response_bin=='good'].HRD_sum; b = hrd_m[hrd_m.response_bin=='bad'].HRD_sum
stat_bracket(ax, 0, 1, max(g.max(), b.max())+1, stats.mannwhitneyu(g,b).pvalue)
ax.set_ylabel('HRD total score')
ax.set_xlabel('')
ax.set_title('HRD genomic scar (total)')

save_panel(fig, 'Fig2E_CNV_HRD', OUT)

print('\n=== Fig 2 (5 panels) complete ===')
