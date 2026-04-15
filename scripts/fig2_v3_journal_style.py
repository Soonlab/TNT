"""
Figure 2 v3 — Journal-style genomic landscape, motifs from landmark papers.
  2A: Bailey/maftools-style comprehensive oncoprint (top TMB stacked bar + main matrix + right freq + bottom annotation)
  2B: Alexandrov-style SBS96 trinucleotide context profile (representative samples)
  2C: PCAWG signature attribution heatmap (sample × signature matrix with annotations)
  2D: Cercek/NEJM-style TMB/MSI waterfall
  2E: TCGA pan-cancer-style genome-wide CNV heatmap
  2F: Knijnenburg-style HRD scar component breakdown
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap

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

# ============================================================
# FIG 2A — Comprehensive oncoprint (Bailey/maftools/ComplexHeatmap style)
# ============================================================
# Layout strategy (Bailey Cell 2018):
#   1. Top header: title
#   2. Sample-level bar plot: TMB stacked by mutation type
#   3. Main mutation matrix: genes × samples, cells colored by mutation type
#   4. Right marginal: gene-level frequency split good/bad with significance
#   5. Bottom annotation tracks: Response, cT, MSI, HLA LOH
#   6. Legend below

# Curated CRC drivers from literature
CRC_DRIVERS = ['APC','TP53','KRAS','PIK3CA','SMAD4','FBXW7','BRAF','NRAS','TCF7L2','KMT2D',
               'KMT2C','ARID1A','SOX9','CTNNB1','AMER1','RNF43','MSH6','MLH1','MSH2','PMS2',
               'POLE','POLD1','ATM','ERBB2','ERBB3','PTEN','CDKN2A','CCND1','MYC','GNAS']

NONSYN_TYPES = {
    'missense_variant': 'Missense',
    'stop_gained': 'Nonsense',
    'frameshift_variant': 'Frameshift',
    'inframe_insertion': 'In-frame indel',
    'inframe_deletion': 'In-frame indel',
    'splice_acceptor_variant': 'Splice site',
    'splice_donor_variant': 'Splice site',
    'protein_altering_variant': 'Other protein-altering',
    'stop_lost': 'Stop lost',
    'start_lost': 'Start lost',
}
TYPE_COLORS = {
    'Missense': '#2ca02c',
    'Nonsense': '#d62728',
    'Frameshift': '#1f77b4',
    'In-frame indel': '#9467bd',
    'Splice site': '#ff7f0e',
    'Other protein-altering': '#8c564b',
    'Stop lost': '#e377c2',
    'Start lost': '#bcbd22',
    'Multiple': '#1d3557',
}

# Filter to driver mutations
v_drv = variants[variants.GENE.isin(CRC_DRIVERS) & variants.is_nonsyn].copy()
v_drv['mut_type'] = v_drv.EFFECT_primary.map(NONSYN_TYPES).fillna('Other protein-altering')

# Pre samples only for cleaner figure
pre_samples = wes_inv[wes_inv.timepoint=='pre'].sample_id.tolist()

# Genes to show (top by frequency in pre samples)
v_drv_pre = v_drv[v_drv.sample_id.isin(pre_samples)]
gene_carriers = v_drv_pre.groupby('GENE').sample_id.nunique().sort_values(ascending=False)
top_genes = [g for g in gene_carriers.index if gene_carriers[g] >= 1][:25]

# Order samples: by response, then by mutation count (more mutated first within group)
def s2subj(s): return int(s.split('-')[0])
sample_meta = []
for s in pre_samples:
    n_mut_in_drivers = v_drv_pre[v_drv_pre.sample_id==s].GENE.nunique()
    resp = wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]
    sample_meta.append({'sample':s,'response':resp,'n_drv':n_mut_in_drivers,'subject':s2subj(s)})
sample_meta = pd.DataFrame(sample_meta).sort_values(['response','n_drv'], ascending=[True, False])
sample_order = sample_meta['sample'].tolist()

# Per-sample TMB stacked by type
tmb_by_type = variants[(variants.is_coding) & variants.sample_id.isin(sample_order)].copy()
tmb_by_type['mt'] = tmb_by_type.EFFECT_primary.map(NONSYN_TYPES).fillna('Synonymous')
tmb_stack = tmb_by_type.groupby(['sample_id','mt']).size().unstack(fill_value=0)
tmb_stack = tmb_stack.reindex(sample_order).fillna(0)

# Build figure
n_samples = len(sample_order); n_genes = len(top_genes)
fig = plt.figure(figsize=(15.5, 11))
gs = fig.add_gridspec(8, 3,
    height_ratios=[1.5, 0.18, 0.18, 0.18, 0.18, n_genes*0.32, 0.5, 0.6],
    width_ratios=[n_samples*0.22, 1.2, 2.6],
    hspace=0.06, wspace=0.04)

# Row 0: TMB stacked bar (top)
ax_tmb = fig.add_subplot(gs[0, 0])
type_order = ['Missense','Nonsense','Frameshift','In-frame indel','Splice site','Other protein-altering']
type_order_present = [t for t in type_order if t in tmb_stack.columns]
bottom = np.zeros(n_samples)
for t in type_order_present:
    vals = tmb_stack[t].values
    ax_tmb.bar(range(n_samples), vals, bottom=bottom, color=TYPE_COLORS[t], width=0.92,
               edgecolor='white', linewidth=0.3, label=t)
    bottom += vals
ax_tmb.set_ylabel('Coding mutations\n(per sample)', fontsize=9)
ax_tmb.set_xlim(-0.5, n_samples-0.5); ax_tmb.set_xticks([])
add_axis_spines(ax_tmb)
ax_tmb.tick_params(labelsize=8)

# Rows 1-4: annotation strips
def ann_strip(ax, vals, label, palette_map=None, cmap=None, vmin=None, vmax=None):
    if palette_map:
        arr = np.array([[matplotlib.colors.to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        v_n = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(v_n)
        if vmax is None: vmax = np.nanmax(v_n)
        norm = ((v_n-vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[-0.5, n_samples-0.5, -0.5, 0.5])
    ax.set_yticks([0]); ax.set_yticklabels([label], fontsize=9, color='#1d3557')
    ax.set_xticks([]); ax.tick_params(length=0)
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(-0.5, n_samples-0.5)

resp_v = [wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in sample_order]
ct_v = [clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in sample_order]
matched_v = [s2subj(s) not in UNMATCHED for s in sample_order]
msi_v = [msi[msi.sample_id==s].MSI_pct.iloc[0] if (msi.sample_id==s).any() else np.nan for s in sample_order]
loh_v = []
for s in sample_order:
    sub = loh[loh['sample']==s]
    loh_v.append(int(sub.LOH_call.sum()) if len(sub) else 0)

ax_a1 = fig.add_subplot(gs[1,0]); ann_strip(ax_a1, resp_v, 'Response', PAL_RESP)
ax_a2 = fig.add_subplot(gs[2,0]); ann_strip(ax_a2, ct_v, 'cT stage', PAL_STAGE)
ax_a3 = fig.add_subplot(gs[3,0]); ann_strip(ax_a3, msi_v, 'MSI %', cmap=plt.cm.Greens, vmin=0, vmax=0.3)
ax_a4 = fig.add_subplot(gs[4,0]); ann_strip(ax_a4, loh_v, 'HLA LOH', cmap=plt.cm.YlOrRd, vmin=0, vmax=3)

# Main oncoprint
ax_onco = fig.add_subplot(gs[5, 0])
mut_lookup = {}
for _, r in v_drv_pre[v_drv_pre.sample_id.isin(sample_order) & v_drv_pre.GENE.isin(top_genes)].iterrows():
    mut_lookup.setdefault((r.sample_id, r.GENE), []).append(r.mut_type)

# Background grey
for i in range(n_genes):
    for j in range(n_samples):
        ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06), 0.92, 0.88,
                                    facecolor='#f2f2f2', edgecolor='white', linewidth=0.3))

# Mutation cells with type color
for (sid, gene), types in mut_lookup.items():
    j = sample_order.index(sid)
    i = top_genes.index(gene)
    types_unique = list(dict.fromkeys(types))  # preserve order
    if len(types_unique) > 1:
        # Multiple types: split horizontally
        h_each = 0.88 / len(types_unique)
        for k, t in enumerate(types_unique):
            color = TYPE_COLORS.get(t, '#8d99ae')
            ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06+k*h_each), 0.92, h_each,
                                        facecolor=color, edgecolor='white', linewidth=0.3))
    else:
        color = TYPE_COLORS.get(types_unique[0], '#8d99ae')
        ax_onco.add_patch(Rectangle((j+0.04, n_genes-1-i+0.06), 0.92, 0.88,
                                    facecolor=color, edgecolor='white', linewidth=0.3))

ax_onco.set_xlim(0, n_samples); ax_onco.set_ylim(0, n_genes)
ax_onco.set_yticks([n_genes-1-i+0.5 for i in range(n_genes)])
ax_onco.set_yticklabels(top_genes, fontsize=9, fontstyle='italic')
ax_onco.set_xticks([])
ax_onco.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_onco.spines[s].set_visible(False)

# Right marginal: % mutated by response with Fisher
ax_rm = fig.add_subplot(gs[5, 1])
n_good_total = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
n_bad_total = sum(1 for s in sample_order if wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
gene_freq = []
for g in top_genes:
    g_mut = sum(1 for s in sample_order if (s,g) in mut_lookup and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='good')
    b_mut = sum(1 for s in sample_order if (s,g) in mut_lookup and wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]=='bad')
    fisher_p = stats.fisher_exact([[g_mut, n_good_total-g_mut],[b_mut, n_bad_total-b_mut]]).pvalue
    gene_freq.append({'gene':g,'good_pct':g_mut/n_good_total*100,'bad_pct':b_mut/n_bad_total*100,'p':fisher_p})
gf_df = pd.DataFrame(gene_freq)
y = np.arange(n_genes)[::-1]
bar_h = 0.35
ax_rm.barh(y-bar_h/2, gf_df.good_pct, height=bar_h, color=GOOD, edgecolor='white', linewidth=0.5, label='Good')
ax_rm.barh(y+bar_h/2, gf_df.bad_pct, height=bar_h, color=BAD, edgecolor='white', linewidth=0.5, label='Poor')
# Significance star
for i, r in gf_df.iterrows():
    if r.p < 0.1:
        x_pos = max(r.good_pct, r.bad_pct) + 2
        star = sig_symbol(r.p)
        ax_rm.text(x_pos, n_genes-1-i, star, va='center', fontsize=10, color='#d62828', fontweight='bold')
ax_rm.set_xlabel('% mutated', fontsize=9)
ax_rm.set_xlim(0, max(gf_df[['good_pct','bad_pct']].max().max()*1.15, 60))
ax_rm.set_yticks([]); ax_rm.set_ylim(-0.5, n_genes-0.5)
ax_rm.legend(fontsize=8, loc='lower right', frameon=False)
add_axis_spines(ax_rm)
ax_rm.tick_params(labelsize=8)

# Right of right: gene category annotation
ax_cat = fig.add_subplot(gs[5, 2])
ax_cat.axis('off')
ax_cat.text(0.05, 0.97, 'Pathway / Function', transform=ax_cat.transAxes, fontsize=10,
            fontweight='bold', va='top', color='#1d3557')
GENE_PATHWAY = {
    'APC':'WNT', 'CTNNB1':'WNT', 'AMER1':'WNT', 'RNF43':'WNT', 'TCF7L2':'WNT',
    'TP53':'TP53', 'CDKN2A':'TP53', 'ATM':'DDR',
    'KRAS':'RTK/RAS', 'BRAF':'RTK/RAS', 'NRAS':'RTK/RAS', 'ERBB2':'RTK/RAS', 'ERBB3':'RTK/RAS',
    'PIK3CA':'PI3K', 'PTEN':'PI3K',
    'SMAD4':'TGFβ',
    'MLH1':'MMR', 'MSH2':'MMR', 'MSH6':'MMR', 'PMS2':'MMR',
    'POLE':'DDR', 'POLD1':'DDR',
    'KMT2D':'Chromatin', 'KMT2C':'Chromatin', 'ARID1A':'Chromatin',
    'FBXW7':'Cell cycle', 'CCND1':'Cell cycle', 'MYC':'Cell cycle',
    'SOX9':'Other', 'GNAS':'Other'
}
PATH_COLORS = {'WNT':'#1f77b4','TP53':'#ff7f0e','RTK/RAS':'#d62728','PI3K':'#9467bd','TGFβ':'#8c564b',
               'MMR':'#2ca02c','DDR':'#e377c2','Chromatin':'#7f7f7f','Cell cycle':'#bcbd22','Other':'#17becf'}
for i, gene in enumerate(top_genes):
    pw = GENE_PATHWAY.get(gene, 'Other')
    color = PATH_COLORS[pw]
    ax_cat.add_patch(Rectangle((0.0, (n_genes-1-i+0.1)/n_genes), 0.04, 0.85/n_genes,
                                facecolor=color, edgecolor='white', linewidth=0.3, transform=ax_cat.transAxes))
    ax_cat.text(0.06, (n_genes-1-i+0.5)/n_genes, pw, transform=ax_cat.transAxes,
                fontsize=8, va='center', color='#1d3557')

# Bottom legend area (rows 6-7): mutation type + pathway
ax_leg1 = fig.add_subplot(gs[6, 0])
ax_leg1.axis('off')
mut_legend = [mpatches.Patch(color=TYPE_COLORS[t], label=t) for t in type_order_present]
ax_leg1.legend(handles=mut_legend, loc='upper left', ncol=4, fontsize=8.5,
               title='Mutation type', title_fontsize=9, frameon=False)

ax_leg2 = fig.add_subplot(gs[7, 0])
ax_leg2.axis('off')
ann_legend_items = []
for k, c in PAL_RESP.items(): ann_legend_items.append(mpatches.Patch(color=c, label=k))
for k, c in PAL_STAGE.items(): ann_legend_items.append(mpatches.Patch(color=c, label=k))
ann_legend_items.append(mpatches.Patch(color=plt.cm.Greens(0.7), label='MSI (high → green)'))
ann_legend_items.append(mpatches.Patch(color=plt.cm.YlOrRd(0.7), label='HLA LOH (count)'))
ax_leg2.legend(handles=ann_legend_items, loc='upper left', ncol=5, fontsize=8.5,
               title='Annotation tracks', title_fontsize=9, frameon=False)

fig.suptitle('Genomic landscape of MSS LARC (n = 33 pre-treatment tumors)\nstacked TMB · driver mutations × samples · per-gene response stratification · molecular annotations',
             fontsize=12, fontweight='bold', y=0.945, color='#1d3557')
save_panel(fig, 'Fig2A_oncoprint_journal', OUT)

# ============================================================
# FIG 2B — Alexandrov-style SBS96 trinucleotide profile
# ============================================================
# SBS96 matrix: rows = 96 mutation contexts, cols = samples
# Standard ordering: C>A, C>G, C>T, T>A, T>C, T>G with 16 contexts each
# Standard COSMIC colors:
SBS96_COLORS = {
    'C>A': '#03BCEE',  # blue
    'C>G': '#010101',  # black
    'C>T': '#E32925',  # red
    'T>A': '#CAC9C9',  # grey
    'T>C': '#A1CE63',  # green
    'T>G': '#EBC6C5',  # pink
}

# Re-organize SBS96 matrix in COSMIC standard order
def parse_context(ctx):
    # Format: A[C>A]A — first char is 5' flank, [REF>ALT], last char is 3' flank
    sub = ctx[2:5]  # C>A
    flank5 = ctx[0]; flank3 = ctx[-1]
    return sub, flank5, flank3

cosmic_order = []
for sub_type in ['C>A','C>G','C>T','T>A','T>C','T>G']:
    for f5 in ['A','C','G','T']:
        for f3 in ['A','C','G','T']:
            ref = sub_type[0]
            cosmic_order.append(f'{f5}[{sub_type}]{f3}')

# sbs96 has rows like 'A[C>A]A'; ensure we match
sbs96_ordered = sbs96.reindex(cosmic_order)

# Pick representative samples to display:
# - 1-PR: reported MMR signature high
# - 12-PR: claimed MMR but found 0%
# - 5-PR: claimed SBS3 but found 0%
# - 14-PR (good responder, contains some SBS15)
samples_to_show = ['1-PR', '12-PR', '5-PR', '14-PR']
samples_to_show = [s for s in samples_to_show if s in sbs96_ordered.columns]

fig, axes = plt.subplots(len(samples_to_show), 1, figsize=(14, 2.4*len(samples_to_show)),
                         sharex=True)
if len(samples_to_show)==1: axes = [axes]
for ax, sid in zip(axes, samples_to_show):
    counts = sbs96_ordered[sid].values
    if counts.sum() == 0:
        ax.text(48, 0.5, f'{sid}: no mutations', ha='center', fontsize=10); continue
    # Color each bar by sub_type group (every 16 bars)
    bar_colors = []
    for i, ctx in enumerate(cosmic_order):
        sub_type, _, _ = parse_context(ctx)
        bar_colors.append(SBS96_COLORS[sub_type])
    ax.bar(range(96), counts, color=bar_colors, edgecolor='black', linewidth=0.15, width=0.85)
    ax.set_ylabel(f'{sid}\n(n={int(counts.sum())} mutations)', fontsize=10)
    # Top color bar with subtype labels
    for i, sub_type in enumerate(['C>A','C>G','C>T','T>A','T>C','T>G']):
        ax.add_patch(Rectangle((i*16-0.45, ax.get_ylim()[1]*1.02), 16, ax.get_ylim()[1]*0.06,
                                facecolor=SBS96_COLORS[sub_type], clip_on=False))
        ax.text(i*16+7.5, ax.get_ylim()[1]*1.16, sub_type, ha='center', fontsize=10,
                color='white' if sub_type in ('C>G',) else 'black', fontweight='bold')
    ax.set_xlim(-0.6, 95.6)
    add_axis_spines(ax)
    ax.tick_params(labelsize=7)

# X axis labels: trinucleotide context (only on bottom)
xlabels = [c.replace('[','').replace(']','') for c in cosmic_order]
axes[-1].set_xticks(range(96))
axes[-1].set_xticklabels(xlabels, rotation=90, fontsize=5.5, family='monospace')
axes[-1].set_xlabel('Trinucleotide mutation context (5\'-N[REF>ALT]N-3\')', fontsize=10)

fig.suptitle('SBS96 trinucleotide mutational context profiles  (Alexandrov-style)',
             fontsize=12, fontweight='bold', y=0.99, color='#1d3557')
fig.tight_layout()
save_panel(fig, 'Fig2B_SBS96_profile', OUT)

# ============================================================
# FIG 2C — Signature attribution heatmap (PCAWG-style)
# ============================================================
sig_cols = [c for c in sbs.columns if c.startswith('SBS')]
active_sigs = [c for c in sig_cols if (sbs[c]>0).any()]
total_per_sig = sbs[active_sigs].sum(axis=0).sort_values(ascending=False)
top_sigs = total_per_sig.head(15).index.tolist()

# Compute proportion matrix
sbs_mat = sbs.set_index('sample_id')[top_sigs]
sbs_total = sbs_mat.sum(axis=1)
sbs_prop = sbs_mat.div(sbs_total.replace(0,np.nan), axis=0).fillna(0)

# Order samples: pre only, by response then total mutation count desc
pre_samples_sbs = [s for s in sample_order if s in sbs_prop.index]
sbs_prop_pre = sbs_prop.loc[pre_samples_sbs]

# Annotation
resp_a = pd.Series([wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in pre_samples_sbs], index=pre_samples_sbs)
total_a = sbs_total.loc[pre_samples_sbs]

# COSMIC etiology labels (subset)
COSMIC_ETIOLOGY = {
    'SBS1': 'Aging (5-mC deamin.)', 'SBS5': 'Aging (clock-like)', 'SBS6': 'MMR deficiency',
    'SBS10b': 'POLE proofreading', 'SBS15': 'MMR deficiency', 'SBS20': 'MMR deficiency',
    'SBS3': 'HRD (BRCA1/2)', 'SBS4': 'Tobacco smoking', 'SBS7b': 'UV exposure',
    'SBS9': 'Polymerase eta', 'SBS18': 'ROS damage', 'SBS40': 'Aging-like (unknown)',
    'SBS54': 'Possible artefact', 'SBS39': 'Unknown', 'SBS24': 'Aflatoxin (likely artefact)',
    'SBS14': 'MMR + POLE', 'SBS19': 'Unknown', 'SBS29': 'Tobacco chewing',
    'SBS30': 'NTHL1 deficiency', 'SBS31': 'Platinum chemotherapy',
    'SBS48': 'Possible artefact', 'SBS49': 'Possible artefact', 'SBS50': 'Possible artefact',
    'SBS87': 'Thiopurine chemotherapy', 'SBS88': 'Colibactin (E. coli)', 'SBS93': 'Unknown',
}
sig_ylabels = [f'{s}  ({COSMIC_ETIOLOGY.get(s,"unknown")})' for s in top_sigs]

fig = plt.figure(figsize=(13.5, 7.5))
gs = fig.add_gridspec(4, 1, height_ratios=[0.4, 0.18, 0.18, 5.5], hspace=0.05)

# Top: total mutation count bar
ax_top = fig.add_subplot(gs[0])
ax_top.bar(range(len(pre_samples_sbs)), total_a.values, color='#1d3557', edgecolor='white', linewidth=0.4)
ax_top.set_ylabel('Total mutations', fontsize=9)
ax_top.set_xticks([])
ax_top.set_xlim(-0.5, len(pre_samples_sbs)-0.5)
add_axis_spines(ax_top)
ax_top.tick_params(labelsize=8)

# Annotation tracks
ax_a1 = fig.add_subplot(gs[1])
ann_strip(ax_a1, resp_a.values.tolist(), 'Response', PAL_RESP)
ct_a = [clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in pre_samples_sbs]
ax_a2 = fig.add_subplot(gs[2])
ann_strip(ax_a2, ct_a, 'cT stage', PAL_STAGE)

# Main heatmap
ax_main = fig.add_subplot(gs[3])
im = ax_main.imshow(sbs_prop_pre[top_sigs].T.values, cmap='Reds', aspect='auto',
                     vmin=0, vmax=0.5, interpolation='nearest')
ax_main.set_yticks(range(len(top_sigs)))
ax_main.set_yticklabels(sig_ylabels, fontsize=8.5)
ax_main.set_xticks(range(len(pre_samples_sbs)))
ax_main.set_xticklabels(pre_samples_sbs, fontsize=7, rotation=90)
ax_main.tick_params(length=2)

# Colorbar
cbar = fig.colorbar(im, ax=ax_main, shrink=0.5, pad=0.01, fraction=0.025)
cbar.set_label('Signature contribution\n(proportion)', fontsize=9)

# Highlight prior MMR/SBS3 claims (subj 1, 5, 9, 12, 14)
HIGHLIGHTED = {1, 5, 9, 12, 14}
for j, s in enumerate(pre_samples_sbs):
    if s2subj(s) in HIGHLIGHTED:
        ax_main.add_patch(Rectangle((j-0.5, -0.5), 1, len(top_sigs), facecolor='none',
                                     edgecolor='#1d3557', linewidth=1.5, linestyle='-'))

ax_main.set_ylabel('')
fig.suptitle('Mutational signature attribution (COSMIC v3.3)  — pre-treatment samples sorted by response × mutation burden\nblack-bordered samples = prior MMR/HRD signature claims',
             fontsize=11.5, fontweight='bold', y=0.995, color='#1d3557')
save_panel(fig, 'Fig2C_signature_attribution', OUT)

# ============================================================
# FIG 2D — TMB / MSI waterfall (Cercek/Yarchoan style)
# ============================================================
# Per-sample MSI score sorted, with TMB scatter overlay
fig, ax1 = plt.subplots(figsize=(13, 4.5))

# Combine: matched samples only, all timepoints
mer = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id')
mer = mer[~mer.subject_id.isin(UNMATCHED)].copy()
mer = mer.sort_values('MSI_pct', ascending=False).reset_index(drop=True)

x = np.arange(len(mer))
colors = [PAL_RESP[r] for r in mer.response_bin]
bars = ax1.bar(x, mer.MSI_pct, color=colors, edgecolor='white', linewidth=0.5, alpha=0.92)

# MSI-H threshold
ax1.axhline(20, color='#d62828', ls='--', lw=1, alpha=0.7)
ax1.text(len(mer)*0.99, 20.5, 'MSI-H clinical cutoff (20%)', ha='right', va='bottom',
         color='#d62828', fontsize=9, fontweight='bold')

# Annotate the highest sample
top_msi = mer.iloc[0]
ax1.annotate(f'{top_msi.sample_id}\nMSI {top_msi.MSI_pct:.2f}%',
             xy=(0, top_msi.MSI_pct), xytext=(2, 2),
             arrowprops=dict(arrowstyle='->', color='#1d3557', lw=0.8),
             fontsize=9, color='#1d3557')

ax1.set_ylabel('MSI percentage (%)', color='#1d3557', fontsize=10)
ax1.set_ylim(0, 25)
ax1.set_xlim(-0.6, len(mer)-0.4)
ax1.set_xticks(x); ax1.set_xticklabels(mer.sample_id, rotation=90, fontsize=6.5)
ax1.tick_params(labelsize=8)
add_axis_spines(ax1)

# Secondary axis: TMB
ax2 = ax1.twinx()
ax2.scatter(x, mer.TMB_nonsyn_per_Mb, color='#9467bd', s=35, edgecolor='white',
            linewidth=0.5, alpha=0.85, zorder=3, label='TMB')
ax2.axhline(10, color='#9467bd', ls=':', lw=1, alpha=0.6)
ax2.set_ylabel('TMB (nonsynonymous mutations / Mb)', color='#9467bd', fontsize=10)
ax2.tick_params(labelcolor='#9467bd', labelsize=8)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_color('#9467bd')

# Legend
fig.legend(handles=[mpatches.Patch(color=GOOD, label='Good (TRG 0-1)'),
                    mpatches.Patch(color=BAD, label='Poor (TRG 2-3)'),
                    matplotlib.lines.Line2D([0],[0], marker='o', color='w', markerfacecolor='#9467bd',
                                            markersize=8, label='TMB (right axis)')],
           loc='upper right', bbox_to_anchor=(0.98, 0.92), fontsize=9, frameon=False, ncol=3)
ax1.set_title('Sample-level MSI percentage (waterfall) and TMB overlay  (matched tumors)\nAll samples are MSS (MSI < 20%) and TMB-low (< 10 / Mb)',
              fontsize=12, fontweight='bold', color='#1d3557', pad=15)
save_panel(fig, 'Fig2D_TMB_MSI_waterfall', OUT)

# ============================================================
# FIG 2E — Genome-wide CNV heatmap (TCGA pan-cancer style)
# ============================================================
import glob
CNV_DIR = ROOT/'04_wes_cnv_clonal/cnvkit'

# Chr sizes (GRCh38)
CHR_SIZE = {'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,
            '6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,
            '11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,
            '16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,
            '21':46709983,'22':50818468,'X':156040895,'Y':57227415}
chr_list = list(CHR_SIZE.keys())[:22]  # autosomes + X

# Build genome-wide bin matrix (5 Mb bins)
BIN_SIZE = 5_000_000
chr_starts = {}
cum = 0
for c in chr_list:
    chr_starts[c] = cum
    cum += CHR_SIZE[c]
total_len = cum
n_bins = total_len // BIN_SIZE + 1

# Per-sample copy number per bin (focus on pre samples)
cn_mat = []
sample_used = []
for s in pre_samples:
    f = CNV_DIR/f'{s}_DNA.call.cns'
    if not f.exists(): continue
    df = pd.read_csv(f, sep='\t')
    df = df[df.chromosome.astype(str).isin(chr_list)].copy()
    bins = np.full(n_bins, 2.0)  # default = neutral
    for _, r in df.iterrows():
        chrom = str(r.chromosome)
        gs = chr_starts[chrom] + int(r.start)
        ge = chr_starts[chrom] + int(r.end)
        b1 = gs // BIN_SIZE; b2 = ge // BIN_SIZE + 1
        bins[b1:b2] = r.cn
    cn_mat.append(bins)
    sample_used.append(s)
cn_mat = np.array(cn_mat)

# Order samples by response then by overall CN dosage
sample_resp = pd.Series([wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in sample_used], index=range(len(sample_used)))
sort_idx = sample_resp.sort_values().index
cn_mat = cn_mat[sort_idx]
sample_used = [sample_used[i] for i in sort_idx]
sample_resp = sample_resp[sort_idx].reset_index(drop=True)

fig = plt.figure(figsize=(15, 7))
gs = fig.add_gridspec(3, 1, height_ratios=[0.3, 0.18, 6], hspace=0.05)

# Top: chromosome ideogram bar
ax_chr = fig.add_subplot(gs[0])
chr_colors_alt = ['#cccccc','#888888']*15
for i, c in enumerate(chr_list):
    start_b = chr_starts[c]//BIN_SIZE
    width_b = CHR_SIZE[c]//BIN_SIZE
    ax_chr.add_patch(Rectangle((start_b, 0), width_b, 1, facecolor=chr_colors_alt[i],
                                edgecolor='white', linewidth=0.4))
    ax_chr.text(start_b + width_b/2, 0.5, c, ha='center', va='center', fontsize=8.5,
                color='black', fontweight='bold')
ax_chr.set_xlim(0, n_bins); ax_chr.set_ylim(0, 1)
ax_chr.set_xticks([]); ax_chr.set_yticks([])
for s in ['top','right','left','bottom']: ax_chr.spines[s].set_visible(False)

# Annotation track
ax_ann = fig.add_subplot(gs[1])
ann_arr = np.array([[matplotlib.colors.to_rgb(PAL_RESP[r]) for r in sample_resp]])
ax_ann.imshow(ann_arr.transpose(1,0,2), aspect='auto', interpolation='nearest',
              extent=[0, n_bins, 0, len(sample_used)])
# Y label on the side of annotation
ax_ann.set_xticks([]); ax_ann.set_yticks([])
ax_ann.set_xlim(0, n_bins); ax_ann.set_ylim(0, len(sample_used))
for s in ['top','right','left','bottom']: ax_ann.spines[s].set_visible(False)
ax_ann.text(-0.005*n_bins, len(sample_used)/2, 'Resp.', ha='right', va='center', fontsize=8,
            transform=ax_ann.transData)

# Main CNV heatmap
ax_main = fig.add_subplot(gs[2])
# Custom diverging cmap centered on 2 (neutral)
cnv_cmap = LinearSegmentedColormap.from_list('cnv', ['#08306b','#2171b5','#deebf7','white','#fee0d2','#fb6a4a','#67000d'])
im = ax_main.imshow(cn_mat, cmap=cnv_cmap, vmin=0, vmax=4, aspect='auto', interpolation='nearest')
ax_main.set_yticks(range(len(sample_used)))
ax_main.set_yticklabels(sample_used, fontsize=6.5)
ax_main.set_xticks([]); ax_main.set_xlim(0, n_bins)
ax_main.tick_params(length=0)

# Chromosome dividers
for i, c in enumerate(chr_list):
    if i == 0: continue
    ax_main.axvline(chr_starts[c]//BIN_SIZE, color='white', lw=0.6)

cbar = fig.colorbar(im, ax=ax_main, shrink=0.4, pad=0.01, fraction=0.025, ticks=[0,1,2,3,4])
cbar.set_label('Copy number', fontsize=9)
cbar.set_ticklabels(['0\n(homo del)','1\n(loss)','2\n(neutral)','3\n(gain)','4+\n(amp)'])
cbar.ax.tick_params(labelsize=7)

fig.suptitle('Genome-wide copy number landscape  (samples × 5 Mb bins; chromosomes 1–22 + X)',
             fontsize=12, fontweight='bold', y=0.94, color='#1d3557')
save_panel(fig, 'Fig2E_CNV_genome', OUT)

# ============================================================
# FIG 2F — HRD scar component breakdown (Knijnenburg style)
# ============================================================
fig, ax = plt.subplots(figsize=(13, 4.5))

hrd_pre = hrd[hrd.timepoint=='pre'].copy()
hrd_pre = hrd_pre[~hrd_pre.subject_id.isin(UNMATCHED)]
hrd_pre = hrd_pre.sort_values('HRD_sum', ascending=False).reset_index(drop=True)

components = ['LST','LOH','TAI']
comp_colors = {'LST':'#118ab2','LOH':'#06aed5','TAI':'#a8dadc'}
bottom = np.zeros(len(hrd_pre))
x = np.arange(len(hrd_pre))
for c in components:
    ax.bar(x, hrd_pre[c], bottom=bottom, color=comp_colors[c], width=0.85,
           edgecolor='white', linewidth=0.6, label=c)
    bottom += hrd_pre[c]

# Response color bar below
for i, (_, r) in enumerate(hrd_pre.iterrows()):
    ax.add_patch(Rectangle((i-0.5, -1.5), 1, 0.7, color=PAL_RESP[r.response_bin], clip_on=False))
ax.text(-0.5, -1.15, 'Resp.', ha='right', va='center', fontsize=9, fontweight='bold')

# LST significance annotation
g = hrd_pre[hrd_pre.response_bin=='good'].LST
b = hrd_pre[hrd_pre.response_bin=='bad'].LST
p_lst = stats.mannwhitneyu(g, b).pvalue
ax.text(0.98, 0.94, f'LST: good vs bad p = {p_lst:.3f}', transform=ax.transAxes,
        ha='right', va='top', fontsize=10, color='#118ab2', fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#118ab2', alpha=0.8, boxstyle='round,pad=0.4'))

ax.set_xticks(x)
ax.set_xticklabels(hrd_pre.sample_id, rotation=90, fontsize=7)
ax.set_ylabel('HRD scar score (LST + LOH + TAI)', fontsize=10)
ax.set_xlim(-0.6, len(hrd_pre)-0.4)
ax.legend(title='Component', loc='upper right', fontsize=9, title_fontsize=9.5,
          frameon=False, bbox_to_anchor=(0.96, 0.85))
add_axis_spines(ax)
ax.tick_params(labelsize=8)
ax.set_title('HRD genomic scar breakdown (matched pre-treatment tumors)',
             fontsize=12, fontweight='bold', color='#1d3557')
save_panel(fig, 'Fig2F_HRD_breakdown', OUT)

print('\n=== Fig 2 v3 (6 journal-style panels) complete ===')
print(f'Output: {OUT}')
