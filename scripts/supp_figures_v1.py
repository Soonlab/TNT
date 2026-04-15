"""
Supplementary Figures S1–S11 for Genome Medicine submission.
Journal-motif-driven design; deep saturated palette matches main figures.

S1  Cohort & QC (Ellrott Cell Syst 2018 / Lek Nature 2016 / Thorsson Immunity 2018)
S2  SBS signatures (Alexandrov Nature 2020 / Degasperi Science 2022)
S3  CNV + HRD (Knijnenburg Cell Rep 2018 / Davoli Science 2017)
S4  Full oncoprint + VAF (Bailey Cell 2018 / Martincorena Cell 2017)
S5  Full GSEA heatmap + leading-edge (Subramanian PNAS 2005 / Liberzon Cell Syst 2015)
S6  ssGSEA corr + CMS/CRIS overlap (Charoentong Cell Rep 2017 / Guinney Nat Med 2015)
S7  Immune deconv + TRUST4 (Newman Nat Meth / Aran Genome Biol / Thorsson Immunity)
S8  ML model comparison (Chowell Nat Biotech 2021 / Cristescu Science 2018)
S9  GEO cohort QC (Butler Nat Biotech / Conway UpSetR)
S10 HLA + neoantigen detail (McGranahan Cell 2017 / Marty Cell 2017)
S11 PyClone-VI diagnostics (Gillis Bioinformatics 2020 / Roth Nat Methods 2014)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
import warnings; warnings.filterwarnings('ignore')

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3_supp'; OUT.mkdir(parents=True, exist_ok=True)

GOOD_DEEP='#0a7d6e'; BAD_DEEP='#c53e1f'; BLACK_DEEP='#0e2a47'; GOLD_DEEP='#d4a300'
GRAY_DEEP='#5a6772'; BAND='#ecedef'
PAL_RESP_DEEP = {'good': GOOD_DEEP, 'bad': BAD_DEEP, np.nan:'#cfd4da'}

# ==================================================================
# S1 — Cohort & QC
# ==================================================================
def supp1():
    wes = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
    rna = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
    tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
    cnv = pd.read_csv(ROOT/'04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')

    fig = plt.figure(figsize=(14, 7.5))
    gs = gridspec.GridSpec(2, 3, hspace=0.45, wspace=0.35)

    # A: WES coverage proxy via n_total variants (quality correlate) per sample, colored by response
    ax = fig.add_subplot(gs[0,0])
    t = tmb.sort_values('n_total').reset_index(drop=True)
    c = [PAL_RESP_DEEP.get(r, '#cfd4da') for r in t.response_bin]
    ax.bar(range(len(t)), t.n_total, color=c, edgecolor=BLACK_DEEP, lw=0.3, width=0.85)
    ax.set_ylabel('Total variants (WES)', fontsize=9.5, color=BLACK_DEEP)
    ax.set_xlabel('WES samples (sorted)', fontsize=9.5, color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=8.5)

    # B: CIN (fraction genome altered) per sample
    ax = fig.add_subplot(gs[0,1])
    cs = cnv.sort_values('CIN').reset_index(drop=True)
    c = [PAL_RESP_DEEP.get(r, '#cfd4da') for r in cs.response_bin]
    ax.bar(range(len(cs)), cs.CIN, color=c, edgecolor=BLACK_DEEP, lw=0.3, width=0.85)
    ax.set_ylabel('Chromosomal instability  (CIN)', fontsize=9.5, color=BLACK_DEEP)
    ax.set_xlabel('WES samples (sorted)', fontsize=9.5, color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=8.5)

    # C: RNA sample matrix by timepoint × response
    ax = fig.add_subplot(gs[0,2])
    mat = pd.crosstab(rna.timepoint.fillna('NA'), rna.response_bin.fillna('NA'))
    im = ax.imshow(mat.values, cmap='Blues', aspect='auto')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, int(mat.iloc[i,j]), ha='center', va='center',
                    fontsize=11, color=BLACK_DEEP, fontweight='bold')
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, fontsize=10)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=10)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)

    # D: Per-subject data availability tile matrix (WES pre/post, RNA pre/post)
    ax = fig.add_subplot(gs[1,:])
    subj_ids = sorted(set(wes.subject_id) | set(rna.subject_id))
    tracks = ['WES-pre','WES-post','RNA-pre','RNA-post']
    M = np.zeros((len(tracks), len(subj_ids)))
    for j, s in enumerate(subj_ids):
        for i, tk in enumerate(tracks):
            src = wes if tk.startswith('WES') else rna
            tp = 'pre' if 'pre' in tk else 'post'
            M[i,j] = int(((src.subject_id==s) & (src.timepoint==tp)).any())
    for i in range(len(tracks)):
        for j in range(len(subj_ids)):
            col = BLACK_DEEP if M[i,j] else '#f1f3f5'
            ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor=col,
                                   edgecolor='white', lw=0.8))
    # response annotation strip on top
    resp_map = wes.drop_duplicates('subject_id').set_index('subject_id').response_bin.to_dict()
    for j, s in enumerate(subj_ids):
        r = resp_map.get(s, None)
        ax.add_patch(Rectangle((j, -1.1), 0.9, 0.7,
                               facecolor=PAL_RESP_DEEP.get(r,'#cfd4da'),
                               edgecolor='white', lw=0.6))
    ax.text(-0.4, -0.75, 'response', ha='right', va='center', fontsize=9,
            color=BLACK_DEEP, fontweight='bold')
    ax.set_xlim(-0.4, len(subj_ids)+0.1); ax.set_ylim(len(tracks)+0.1, -1.3)
    ax.set_xticks([j+0.45 for j in range(len(subj_ids))])
    ax.set_xticklabels([f'S{int(s)}' for s in subj_ids], fontsize=7.5, rotation=90,
                       color=BLACK_DEEP)
    ax.set_yticks([i+0.45 for i in range(len(tracks))])
    ax.set_yticklabels(tracks, fontsize=9, color=BLACK_DEEP)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)
    save_panel(fig, 'FigS1_cohort_QC', OUT)

# ==================================================================
# S2 — SBS signatures
# ==================================================================
def supp2():
    sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_activities_with_meta.tsv', sep='\t')
    sig_cols = [c for c in sbs.columns if c.startswith('SBS')]
    # top 8 most active signatures across cohort
    totals = sbs[sig_cols].sum().sort_values(ascending=False)
    top = totals.head(8).index.tolist()

    fig = plt.figure(figsize=(14, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.55)

    # A: Stacked relative exposure per sample, ordered by response then dominant sig
    ax = fig.add_subplot(gs[0])
    rel = sbs[top].div(sbs[sig_cols].sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    rel['other'] = 1 - rel.sum(axis=1)
    rel_cols = top + ['other']
    meta = sbs[['sample_id','timepoint','response_bin']].copy()
    meta['dom'] = rel[top].idxmax(axis=1)
    order = meta.sort_values(['response_bin','dom']).index
    rel_s = rel.loc[order, rel_cols]
    import matplotlib.cm as cm
    palette = list(plt.get_cmap('tab10').colors)[:len(rel_cols)]
    palette[-1] = '#c7ccd3'  # 'other' gray
    x = np.arange(len(rel_s))
    bottom = np.zeros(len(rel_s))
    for col, color in zip(rel_cols, palette):
        ax.bar(x, rel_s[col].values, bottom=bottom, color=color, edgecolor='white',
               linewidth=0.3, width=0.92, label=col)
        bottom += rel_s[col].values
    ax.set_ylim(0,1.02); ax.set_xticks(x)
    ax.set_xticklabels(meta.loc[order,'sample_id'].values, rotation=90, fontsize=6.5,
                       color=BLACK_DEEP)
    ax.set_ylabel('Relative SBS exposure', fontsize=10, color=BLACK_DEEP)
    ax.legend(ncol=len(rel_cols), fontsize=7.5, frameon=False,
              bbox_to_anchor=(0.5, 1.03), loc='lower center')
    add_axis_spines(ax); ax.tick_params(labelsize=8)

    # B: Per-signature median exposure by response (good vs bad) horizontal bars
    ax = fig.add_subplot(gs[1])
    top20 = totals.head(20).index
    med_good = sbs[sbs.response_bin=='good'][top20].median()
    med_bad  = sbs[sbs.response_bin=='bad' ][top20].median()
    ypos = np.arange(len(top20))
    ax.barh(ypos-0.2, med_good.values, height=0.38, color=GOOD_DEEP,
            edgecolor=BLACK_DEEP, lw=0.4, label='good (median)', alpha=0.9)
    ax.barh(ypos+0.2, med_bad.values, height=0.38, color=BAD_DEEP,
            edgecolor=BLACK_DEEP, lw=0.4, label='bad (median)', alpha=0.9)
    ax.set_yticks(ypos); ax.set_yticklabels(top20, fontsize=9, color=BLACK_DEEP)
    ax.invert_yaxis()
    ax.set_xlabel('Signature activity (mutations assigned)', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS2_SBS_signatures', OUT)

# ==================================================================
# S3 — CNV + HRD
# ==================================================================
def supp3():
    cnv = pd.read_csv(ROOT/'04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
    hrd = pd.read_csv(ROOT/'04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv', sep='\t')

    fig = plt.figure(figsize=(12.5, 6.8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1,1], width_ratios=[3,2],
                           hspace=0.5, wspace=0.3)

    # A: per-sample stacked HRD subscores
    ax = fig.add_subplot(gs[0,:])
    h = hrd.sort_values(['response_bin','HRD_sum']).reset_index(drop=True)
    x = np.arange(len(h)); bottom = np.zeros(len(h))
    for col, color, lbl in [('LST','#1d3557','LST'),
                             ('LOH','#118ab2','LOH'),
                             ('TAI','#06aed5','TAI')]:
        ax.bar(x, h[col], bottom=bottom, color=color, edgecolor='white',
               lw=0.3, width=0.9, label=lbl)
        bottom += h[col].values
    # response strip
    for i, r in enumerate(h.response_bin):
        ax.add_patch(Rectangle((i-0.45, -h.HRD_sum.max()*0.08),
                               0.9, h.HRD_sum.max()*0.06,
                               facecolor=PAL_RESP_DEEP.get(r,'#cfd4da'),
                               edgecolor='white', lw=0.3, clip_on=False))
    ax.set_xticks(x)
    ax.set_xticklabels(h.sample_id.values, rotation=90, fontsize=6.5, color=BLACK_DEEP)
    ax.set_ylabel('HRD proxy subscore', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, ncol=3, loc='upper left',
              bbox_to_anchor=(0.0, 1.0))
    add_axis_spines(ax); ax.tick_params(labelsize=8)

    # B: HRD_sum by response
    ax = fig.add_subplot(gs[1,0])
    for i, g in enumerate(['good','bad']):
        vals = hrd[hrd.response_bin==g].HRD_sum.values
        c = PAL_RESP_DEEP[g]
        ax.scatter(np.full(len(vals), i) + np.random.uniform(-0.1,0.1,len(vals)),
                   vals, s=50, color=c, edgecolor=BLACK_DEEP, linewidth=0.5, alpha=0.85)
        ax.plot([i-0.25, i+0.25], [np.median(vals)]*2, color=c, lw=2.5)
    p = stats.mannwhitneyu(hrd[hrd.response_bin=='good'].HRD_sum,
                           hrd[hrd.response_bin=='bad'].HRD_sum).pvalue
    ax.text(0.5, 0.98, f'p = {p:.2f}', transform=ax.transAxes, ha='center', va='top',
            fontsize=10, color=BLACK_DEEP, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=BLACK_DEEP, lw=0.6))
    ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'], fontsize=10)
    ax.set_ylabel('HRD sum  (LST + LOH + TAI)', fontsize=10, color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # C: CIN vs HRD_sum scatter
    ax = fig.add_subplot(gs[1,1])
    merged = cnv.merge(hrd[['sample_id','HRD_sum']], on='sample_id')
    for g in ['good','bad']:
        s = merged[merged.response_bin==g]
        ax.scatter(s.CIN, s.HRD_sum, s=55, color=PAL_RESP_DEEP[g], edgecolor=BLACK_DEEP,
                   linewidth=0.5, alpha=0.85, label=g)
    from scipy.stats import spearmanr
    rho, p = spearmanr(merged.CIN, merged.HRD_sum)
    ax.text(0.02, 0.98, f'ρ = {rho:.2f}, p = {p:.2g}', transform=ax.transAxes,
            ha='left', va='top', fontsize=9.5, fontweight='bold', color=BLACK_DEEP,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=BLACK_DEEP, lw=0.6))
    ax.set_xlabel('CIN  (fraction genome altered)', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('HRD sum', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS3_CNV_HRD', OUT)

# ==================================================================
# S4 — Full driver oncoprint + VAF
# ==================================================================
def supp4():
    onc = pd.read_csv(ROOT/'04_wes_cnv_clonal/driver_oncoprint_matrix.tsv', sep='\t')
    drivers = pd.read_csv(ROOT/'04_wes_cnv_clonal/driver_mutations.tsv', sep='\t')
    wes = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')

    genes = onc.GENE.tolist()
    samples = [c for c in onc.columns if c not in ['GENE','total_samples_mutated']]
    M = onc.set_index('GENE')[samples].values
    # order samples by response then mutation burden
    s_info = wes.set_index('sample_id').loc[samples]
    s_info['mut_sum'] = M.sum(axis=0)
    s_info = s_info.sort_values(['response_bin','mut_sum'], ascending=[True, False])
    s_order = s_info.index.tolist()

    fig = plt.figure(figsize=(13, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.1], hspace=0.25)

    # A: Oncoprint
    ax = fig.add_subplot(gs[0])
    onc2 = onc.set_index('GENE')[s_order]
    M2 = onc2.values
    for i in range(M2.shape[0]):
        for j in range(M2.shape[1]):
            c = BLACK_DEEP if M2[i,j]>0 else '#f1f3f5'
            ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor=c,
                                   edgecolor='white', lw=0.6))
    # response track
    for j, s in enumerate(s_order):
        r = s_info.loc[s, 'response_bin']
        ax.add_patch(Rectangle((j, -1.1), 0.9, 0.7,
                               facecolor=PAL_RESP_DEEP.get(r,'#cfd4da'),
                               edgecolor='white', lw=0.4))
    ax.set_xlim(-0.4, len(s_order)+0.1)
    ax.set_ylim(len(genes)+0.2, -1.3)
    ax.set_yticks([i+0.45 for i in range(len(genes))])
    ax.set_yticklabels(genes, fontsize=8, color=BLACK_DEEP)
    ax.set_xticks([j+0.45 for j in range(len(s_order))])
    ax.set_xticklabels(s_order, rotation=90, fontsize=6.5, color=BLACK_DEEP)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    # B: Top mutated driver genes — count per response
    ax = fig.add_subplot(gs[1])
    onc_long = onc2.reset_index().melt(id_vars='GENE', var_name='sample_id',
                                        value_name='mut').query('mut > 0')
    onc_long = onc_long.merge(wes[['sample_id','response_bin']], on='sample_id', how='left')
    top_g = onc2.sum(axis=1).sort_values(ascending=False).head(15).index
    counts = (onc_long[onc_long.GENE.isin(top_g)]
              .groupby(['GENE','response_bin']).size().unstack(fill_value=0))
    counts = counts.reindex(top_g)
    if 'good' not in counts.columns: counts['good']=0
    if 'bad'  not in counts.columns: counts['bad'] =0
    y = np.arange(len(counts))
    ax.barh(y-0.2, counts['good'].values, height=0.38, color=GOOD_DEEP,
            edgecolor=BLACK_DEEP, lw=0.4, label='good', alpha=0.9)
    ax.barh(y+0.2, counts['bad'].values,  height=0.38, color=BAD_DEEP,
            edgecolor=BLACK_DEEP, lw=0.4, label='bad',  alpha=0.9)
    ax.set_yticks(y); ax.set_yticklabels(counts.index, fontsize=8.5, color=BLACK_DEEP)
    ax.invert_yaxis()
    ax.set_xlabel('# mutated samples', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS4_oncoprint_VAF', OUT)

# ==================================================================
# S5 — Full GSEA heatmap + leading edge
# ==================================================================
def supp5():
    gh = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
    gr = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')
    gh['source']='Hallmark'; gr['source']='Reactome'
    g = pd.concat([gh, gr], ignore_index=True)
    g['sig'] = g.padj < 0.25
    # keep top 40 by |NES|
    top = g.reindex(g.NES.abs().sort_values(ascending=False).index).head(40)
    top = top.sort_values('NES').reset_index(drop=True)

    fig = plt.figure(figsize=(11, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.35, 4], wspace=0.02)
    ax_src = fig.add_subplot(gs[0,0]); ax = fig.add_subplot(gs[0,1])

    # source track
    src_col = {'Hallmark':GOLD_DEEP, 'Reactome':'#118ab2'}
    for i, s in enumerate(top.source):
        ax_src.add_patch(Rectangle((0, i), 1, 0.9, facecolor=src_col[s],
                                   edgecolor='white', lw=0.5))
    ax_src.set_xlim(0,1); ax_src.set_ylim(len(top), -0.5)
    for sp in ax_src.spines.values(): sp.set_visible(False)
    ax_src.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # NES bars
    y = np.arange(len(top))
    colors = [GOOD_DEEP if v>0 else BAD_DEEP for v in top.NES]
    ax.barh(y, top.NES, color=colors, edgecolor=BLACK_DEEP, lw=0.3, height=0.78,
            alpha=0.92)
    x_label = top.NES.max()*1.08
    for i, r in top.iterrows():
        star = '**' if r.padj<0.05 else ('*' if r.padj<0.1 else '')
        ax.text(x_label, i,
                f'{r.pathway[:55]}{"…" if len(r.pathway)>55 else ""}   q={r.padj:.2g}{star}',
                va='center', ha='left', fontsize=7.5, color=BLACK_DEEP)
    ax.axvline(0, color=BLACK_DEEP, lw=0.8)
    ax.invert_yaxis()
    ax.set_yticks([]); ax.set_xlabel('NES  (good − bad)', fontsize=10, color=BLACK_DEEP)
    ax.set_xlim(top.NES.min()*1.2, top.NES.max()*4.2)
    # legend
    leg = [mpatches.Patch(color=src_col['Hallmark'], label='Hallmark'),
           mpatches.Patch(color=src_col['Reactome'], label='Reactome')]
    ax.legend(handles=leg, loc='lower right', fontsize=9, frameon=False)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS5_GSEA_full', OUT)

# ==================================================================
# S6 — ssGSEA correlation + CMS distribution
# ==================================================================
def supp6():
    ss = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t')
    cms = pd.read_csv(ROOT/'07_rna_cms/cms_assignments.tsv', sep='\t')

    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1], wspace=0.3)

    # A: correlation heatmap of selected pathway modules
    ax = fig.add_subplot(gs[0])
    # select a curated subset for readability
    cols_keep = [c for c in ss.columns if c != 'sample_id']
    # pick 20 by high variance across samples
    variances = ss[cols_keep].var().sort_values(ascending=False)
    sel = variances.head(20).index
    corr = ss[sel].corr()
    # order by hierarchical cluster
    from scipy.cluster.hierarchy import linkage, leaves_list
    Z = linkage(corr.values, method='average')
    idx = leaves_list(Z)
    corr_o = corr.iloc[idx, idx]
    im = ax.imshow(corr_o.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(len(corr_o)))
    # shorten names
    short = [s.split(' R-HSA')[0][:28] for s in corr_o.columns]
    ax.set_xticklabels(short, rotation=45, ha='right', fontsize=7, color=BLACK_DEEP)
    ax.set_yticks(range(len(corr_o)))
    ax.set_yticklabels(short, fontsize=7, color=BLACK_DEEP)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, label='Spearman ρ')
    # B: CMS distribution by response
    ax = fig.add_subplot(gs[1])
    cms_bin = cms[cms.timepoint=='pre'] if 'timepoint' in cms.columns else cms
    tab = pd.crosstab(cms_bin.prediction, cms_bin.response_bin)
    tab = tab.reindex(index=['CMS1','CMS2','CMS3','CMS4'], fill_value=0)
    pct = tab.div(tab.sum(axis=0), axis=1) * 100
    x = np.arange(len(pct.columns))
    bottom = np.zeros(len(pct.columns))
    cms_col = {'CMS1':'#ee6c4d','CMS2':'#118ab2','CMS3':'#06aed5','CMS4':'#7b2cbf'}
    for cms_k in ['CMS1','CMS2','CMS3','CMS4']:
        vals = pct.loc[cms_k].values
        ax.bar(x, vals, bottom=bottom, color=cms_col[cms_k], edgecolor='white',
               linewidth=1.3, label=cms_k, alpha=0.92, width=0.55)
        for i,(v,tot) in enumerate(zip(vals, tab.loc[cms_k].values)):
            if v>4:
                ax.text(i, bottom[i]+v/2, f'{v:.0f}%\n(n={int(tot)})', ha='center',
                        va='center', fontsize=8.5, color='white', fontweight='bold')
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(pct.columns, fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('% samples', fontsize=10, color=BLACK_DEEP)
    ax.set_ylim(0, 105)
    ax.legend(loc='lower center', ncol=4, fontsize=9, frameon=False,
              bbox_to_anchor=(0.5, -0.2))
    from scipy.stats import chi2_contingency
    try: _,chi_p,_,_ = chi2_contingency(tab.values)
    except: chi_p = np.nan
    ax.text(0.5, 1.02, f'χ² p = {chi_p:.2f}', transform=ax.transAxes, ha='center',
            fontsize=10, color=BLACK_DEEP, fontweight='bold')
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS6_ssGSEA_CMS', OUT)

# ==================================================================
# S7 — Immune signatures per-sample + TRUST4
# ==================================================================
def supp7():
    sig = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t')
    rna = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
    tr = pd.read_csv(ROOT/'06_rna_immune/trust4_summary.tsv', sep='\t')
    sig = sig.merge(rna[['sample_id','response_bin','timepoint']], on='sample_id', how='left')

    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.3, 1], hspace=0.6)

    # A: per-sample signature heatmap (z-scored columns)
    ax = fig.add_subplot(gs[0])
    cols = [c for c in sig.columns if c not in ['sample_id','response_bin','timepoint']]
    pre = sig[sig.timepoint=='pre'].copy().sort_values('response_bin')
    Z = (pre[cols] - pre[cols].mean()) / pre[cols].std()
    im = ax.imshow(Z.T.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
    ax.set_yticks(range(len(cols))); ax.set_yticklabels(cols, fontsize=8, color=BLACK_DEEP)
    ax.set_xticks(range(len(pre))); ax.set_xticklabels(pre.sample_id.values, rotation=90,
                                                       fontsize=6.5, color=BLACK_DEEP)
    # response strip
    for j, r in enumerate(pre.response_bin.values):
        ax.add_patch(Rectangle((j-0.5, len(cols)+0.1), 1, 0.5,
                               facecolor=PAL_RESP_DEEP.get(r,'#cfd4da'),
                               edgecolor='white', lw=0.3, clip_on=False))
    plt.colorbar(im, ax=ax, label='z-score', shrink=0.7)
    # B: TRUST4 clonotype metrics per response
    ax = fig.add_subplot(gs[1])
    tr_pre = tr[tr.timepoint=='pre']
    metrics = [('TRB_n','TRB clonotypes'),
               ('TRB_shannon','TRB Shannon'),
               ('IGH_shannon','IGH Shannon'),
               ('TRB_gini','TRB Gini')]
    positions = np.arange(len(metrics))
    width = 0.35
    for i, (m, lbl) in enumerate(metrics):
        for k, g in enumerate(['good','bad']):
            vals = tr_pre[tr_pre.response_bin==g][m].dropna().values
            if len(vals)==0: continue
            xc = i + (-width/2 if g=='good' else +width/2)
            c = PAL_RESP_DEEP[g]
            ax.scatter(np.full(len(vals), xc) + np.random.uniform(-0.06,0.06,len(vals)),
                       vals / np.nanmax(tr_pre[m]) if m.endswith('_n') else vals,
                       s=40, color=c, edgecolor=BLACK_DEEP, linewidth=0.5, alpha=0.85)
    ax.set_xticks(positions); ax.set_xticklabels([l for _,l in metrics], fontsize=9,
                                                 color=BLACK_DEEP)
    ax.set_ylabel('Value (n-metrics normalized to cohort max)', fontsize=9.5,
                  color=BLACK_DEEP)
    leg = [mpatches.Patch(color=GOOD_DEEP, label='good'),
           mpatches.Patch(color=BAD_DEEP,  label='bad')]
    ax.legend(handles=leg, loc='upper right', fontsize=9, frameon=False)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS7_immune_TRUST4', OUT)

# ==================================================================
# S8 — ML model comparison
# ==================================================================
def supp8():
    ml = pd.read_csv(ROOT/'10_ml_predictor/ml_loocv_results.tsv', sep='\t')
    rf = pd.read_csv(ROOT/'10_ml_predictor/rf_feature_importance.tsv', sep='\t')
    # bootstrap CI for AUC: use Hanley approximate SE given LOOCV AUC, n=35
    n = 35
    def auc_se(auc, n_pos=14, n_neg=21):
        Q1 = auc/(2-auc); Q2 = 2*auc**2/(1+auc)
        se = np.sqrt((auc*(1-auc) + (n_pos-1)*(Q1-auc**2) +
                      (n_neg-1)*(Q2-auc**2)) / (n_pos*n_neg))
        return se

    fig = plt.figure(figsize=(12, 5.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.1, 1.6], wspace=0.85)

    # A: AUC with 95% CI per model
    ax = fig.add_subplot(gs[0])
    ml = ml.sort_values('LOOCV_AUC').reset_index(drop=True)
    y = np.arange(len(ml))
    for i, r in ml.iterrows():
        se = auc_se(r.LOOCV_AUC)
        lo, hi = r.LOOCV_AUC - 1.96*se, r.LOOCV_AUC + 1.96*se
        c = GOOD_DEEP if r.LOOCV_AUC>=0.7 else (GOLD_DEEP if r.LOOCV_AUC>=0.6 else GRAY_DEEP)
        ax.plot([lo, hi], [i,i], color=c, lw=2, alpha=0.85)
        ax.scatter(r.LOOCV_AUC, i, s=160, color=c, edgecolor=BLACK_DEEP, linewidth=0.8,
                   zorder=3)
        ax.text(hi+0.02, i, f'AUC={r.LOOCV_AUC:.3f}  Acc={r.LOOCV_Accuracy:.2f}',
                va='center', fontsize=9, color=BLACK_DEEP, fontweight='bold')
    ax.axvline(0.5, color=GRAY_DEEP, ls='--', lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(ml.model.values, fontsize=10, color=BLACK_DEEP)
    ax.set_xlim(0.3, 1.05)
    ax.set_xlabel('LOOCV AUC  (95% Hanley CI)', fontsize=10, color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # B: RF top feature importance
    ax = fig.add_subplot(gs[1])
    top = rf.sort_values('importance', ascending=True).tail(20).reset_index(drop=True)
    y = np.arange(len(top))
    ax.barh(y, top.importance, color=BLACK_DEEP, edgecolor=BLACK_DEEP, lw=0.3,
            height=0.72, alpha=0.88)
    ax.set_yticks(y)
    ax.set_yticklabels([f[:55] for f in top.feature.values], fontsize=8,
                       color=BLACK_DEEP)
    ax.set_xlabel('RF feature importance (mean decrease Gini)', fontsize=10,
                  color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=8.5)
    save_panel(fig, 'FigS8_ML_model_comparison', OUT)

# ==================================================================
# S9 — GEO cohort QC (summary table + per-cohort n bars)
# ==================================================================
def supp9():
    stat = pd.read_csv(ROOT/'11_external_validation/signature_stats_manual.tsv', sep='\t')
    coh_sum = stat.groupby('gse').agg(n_good=('n_good','first'),
                                       n_bad=('n_bad','first')).reset_index()
    coh_sum['n_total'] = coh_sum.n_good + coh_sum.n_bad
    coh_sum = coh_sum.sort_values('n_total', ascending=True).reset_index(drop=True)

    fig = plt.figure(figsize=(12, 5.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.4)

    # A: stacked good/bad per cohort
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(coh_sum))
    ax.barh(y, coh_sum.n_good, color=GOOD_DEEP, edgecolor=BLACK_DEEP, lw=0.4,
            height=0.7, label='good', alpha=0.92)
    ax.barh(y, coh_sum.n_bad, left=coh_sum.n_good, color=BAD_DEEP,
            edgecolor=BLACK_DEEP, lw=0.4, height=0.7, label='bad', alpha=0.92)
    for i, r in coh_sum.iterrows():
        ax.text(r.n_total+1, i, f'n={int(r.n_total)}  ({int(r.n_good)}/{int(r.n_bad)})',
                va='center', fontsize=9, color=BLACK_DEEP, fontweight='bold')
    ax.set_yticks(y); ax.set_yticklabels(coh_sum.gse, fontsize=9.5, color=BLACK_DEEP)
    ax.set_xlabel('# patients', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # B: per-cohort signature direction tile
    ax = fig.add_subplot(gs[1])
    sig_names = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
    expected = {'DSB_HDR_repair':1,'E2F_MYC_cellcycle':1,'CD8_proliferation':1,'EMT':-1}
    piv = stat.pivot(index='gse', columns='signature', values='delta')[sig_names]
    piv = piv.reindex(coh_sum.gse)
    for i, gse in enumerate(piv.index):
        for j, sig in enumerate(sig_names):
            d = piv.loc[gse, sig]
            match = np.sign(d)==expected[sig]
            c = GOOD_DEEP if match else BAD_DEEP
            alpha = min(1.0, 0.35 + abs(d)*1.2)
            ax.add_patch(Rectangle((j, i), 0.92, 0.92, facecolor=c,
                                   edgecolor='white', lw=1.2, alpha=alpha))
            ax.text(j+0.46, i+0.46, f'{d:+.2f}', ha='center', va='center',
                    fontsize=7.5, color='white', fontweight='bold')
    ax.set_xlim(-0.2, len(sig_names)+0.1)
    ax.set_ylim(len(piv)+0.1, -0.2)
    ax.set_xticks([j+0.46 for j in range(len(sig_names))])
    ax.set_xticklabels([s.replace('_','\n') for s in sig_names], fontsize=8.5,
                       color=BLACK_DEEP)
    ax.set_yticks([i+0.46 for i in range(len(piv))])
    ax.set_yticklabels(piv.index, fontsize=9, color=BLACK_DEEP)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)
    save_panel(fig, 'FigS9_GEO_cohorts', OUT)

# ==================================================================
# S10 — HLA + neoantigen detail
# ==================================================================
def supp10():
    hla = pd.read_csv(ROOT/'03_hla/hla_class_I_typing.tsv', sep='\t')
    neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')

    fig = plt.figure(figsize=(13, 7))
    gs = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.3)

    # A: peptide per site (neoantigens per site) distribution by HLA LOH status
    ax = fig.add_subplot(gs[0,0])
    neo_pre = neo[(neo.timepoint=='pre') & neo.matched]
    for loh_flag, lbl, col in [(False,'HLA retained', GOOD_DEEP),
                                (True, 'HLA LOH',       BAD_DEEP)]:
        vals = neo_pre[neo_pre.HLA_LOH==loh_flag].neoantigens_per_site.dropna().values
        if len(vals)==0: continue
        ax.hist(vals, bins=12, alpha=0.55, color=col, edgecolor=BLACK_DEEP, lw=0.6,
                label=f'{lbl} (n={len(vals)})')
    ax.set_xlabel('Neoantigens per mutation site', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('# samples', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # B: binder yield (n_binders / n_candidates) per sample
    ax = fig.add_subplot(gs[0,1])
    yld = neo_pre.assign(yield_pct=100*neo_pre.n_binders_500nM/neo_pre.n_candidate_peptides.replace(0,np.nan))
    for g in ['good','bad']:
        s = yld[yld.response_bin==g].yield_pct.dropna()
        ax.scatter(np.full(len(s), {'good':0,'bad':1}[g]) + np.random.uniform(-0.12,0.12,len(s)),
                   s, s=55, color=PAL_RESP_DEEP[g], edgecolor=BLACK_DEEP,
                   linewidth=0.5, alpha=0.85)
    ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'], fontsize=10)
    ax.set_ylabel('% candidate peptides that bind (<500 nM)', fontsize=9.5,
                  color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # C: HLA zygosity summary per locus
    ax = fig.add_subplot(gs[1,0])
    hla_g = hla[hla.is_germline==True]
    cats = ['A','B','C']
    bars = [(hla_g[f'homozygous_{c}'].sum(), len(hla_g)-hla_g[f'homozygous_{c}'].sum()) for c in cats]
    x = np.arange(len(cats))
    hom = [b[0] for b in bars]; het = [b[1] for b in bars]
    for i, c in enumerate(cats):
        ax.bar(i, bars[i][1], color={'A':'#118ab2','B':'#1d3557','C':'#06aed5'}[c],
               edgecolor='white', lw=1)
        ax.bar(i, bars[i][0], bottom=bars[i][1], color=GOLD_DEEP, edgecolor='white', lw=1)
        ax.text(i, len(hla_g)+0.5, f'{bars[i][0]} hom', ha='center', fontsize=9,
                color=BLACK_DEEP, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels([f'HLA-{c}' for c in cats], fontsize=10,
                                          color=BLACK_DEEP)
    ax.set_ylabel('# patients', fontsize=10, color=BLACK_DEEP)
    ax.set_ylim(0, len(hla_g)*1.15)
    leg = [mpatches.Patch(color=GOLD_DEEP, label='homozygous'),
           mpatches.Patch(color='#118ab2', label='heterozygous')]
    ax.legend(handles=leg, fontsize=9, frameon=False, loc='lower right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # D: neoantigen count vs HLA homozygosity count
    ax = fig.add_subplot(gs[1,1])
    merged = neo_pre.merge(hla_g[['subject_id','n_homozygous_loci']], on='subject_id',
                           how='left')
    for g in ['good','bad']:
        s = merged[merged.response_bin==g]
        ax.scatter(s.n_homozygous_loci + np.random.uniform(-0.12, 0.12, len(s)),
                   s.n_sites_with_binder, s=60, color=PAL_RESP_DEEP[g],
                   edgecolor=BLACK_DEEP, linewidth=0.5, alpha=0.85, label=g)
    ax.set_xticks([0,1,2,3])
    ax.set_xlabel('# homozygous HLA class I loci', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('# sites with MHC-I binder', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False, loc='upper right')
    add_axis_spines(ax); ax.tick_params(labelsize=9)
    save_panel(fig, 'FigS10_HLA_neoantigen_detail', OUT)

# ==================================================================
# S11 — PyClone-VI diagnostics
# ==================================================================
def supp11():
    pyclone = pd.read_csv(ROOT/'04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\t')
    pyclone = pyclone.sort_values(['response','subject_id']).reset_index(drop=True)
    PY = ROOT/'04_wes_cnv_clonal/pyclone'

    fig = plt.figure(figsize=(13, 8))
    gs = gridspec.GridSpec(3, 4, hspace=0.65, wspace=0.25)

    # A: cluster count vs n_muts
    ax = fig.add_subplot(gs[0,:2])
    for g in ['good','bad']:
        s = pyclone[pyclone.response==g]
        ax.scatter(s.n_muts, s.n_clusters, s=100, color=PAL_RESP_DEEP[g],
                   edgecolor=BLACK_DEEP, linewidth=0.7, alpha=0.88, label=g)
        for _, r in s.iterrows():
            ax.text(r.n_muts+1, r.n_clusters+0.05, f'S{int(r.subject_id)}',
                    fontsize=7.5, color=BLACK_DEEP)
    ax.set_xlabel('# mutations fitted', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('# PyClone clusters', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # B: dominant shrink vs expand
    ax = fig.add_subplot(gs[0,2:])
    for g in ['good','bad']:
        s = pyclone[pyclone.response==g]
        ax.scatter(s.dominant_shrink, s.dominant_expand, s=100, color=PAL_RESP_DEEP[g],
                   edgecolor=BLACK_DEEP, linewidth=0.7, alpha=0.88, label=g)
    ax.axhline(0, color=GRAY_DEEP, ls='--', lw=0.6)
    ax.axvline(0, color=GRAY_DEEP, ls='--', lw=0.6)
    ax.set_xlabel('Dominant Δ CP (shrink)', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('Dominant Δ CP (expand)', fontsize=10, color=BLACK_DEEP)
    ax.legend(fontsize=9, frameon=False)
    add_axis_spines(ax); ax.tick_params(labelsize=9)

    # C: per-subject CCF ridgeplot (Roth motif)
    ax_r = fig.add_subplot(gs[1:, :])
    for i, s in enumerate(pyclone.subject_id.values):
        f = PY/f'results_subj{int(s)}.tsv'
        if not f.exists(): continue
        r = pd.read_csv(f, sep='\t')
        r['tp'] = r.sample_id.str.extract(r'-(PR|PO)$')[0].map({'PR':'pre','PO':'post'})
        for tp, offset in [('pre', -0.25), ('post', +0.25)]:
            vals = r[r.tp==tp].cellular_prevalence.dropna().values
            if len(vals)<3 or np.std(vals)<1e-6: continue
            try:
                kde = stats.gaussian_kde(vals, bw_method=0.4)
                xg = np.linspace(0,1,120)
                d = kde(xg); d = d/d.max()*0.42
                base = i + offset
                c = PAL_RESP_DEEP[pyclone[pyclone.subject_id==s].response.iloc[0]]
                ls_style = '-' if tp=='pre' else '--'
                ax_r.fill_between(xg, base, base-d, color=c, alpha=0.25, linewidth=0)
                ax_r.plot(xg, base-d, color=c, lw=1, ls=ls_style, alpha=0.9)
            except Exception: pass
    ax_r.set_yticks(range(len(pyclone)))
    ax_r.set_yticklabels([f'S{int(s)}' for s in pyclone.subject_id.values], fontsize=8,
                         color=BLACK_DEEP)
    ax_r.invert_yaxis()
    ax_r.set_xlabel('Cellular prevalence', fontsize=10, color=BLACK_DEEP)
    ax_r.set_xlim(0, 1); ax_r.set_ylim(len(pyclone)-0.2, -0.8)
    add_axis_spines(ax_r); ax_r.tick_params(labelsize=9)
    save_panel(fig, 'FigS11_pyclone_diagnostics', OUT)

# ==================================================================
# RUN ALL
# ==================================================================
for fn in [supp1, supp2, supp3, supp4, supp5, supp6, supp7, supp8, supp9, supp10, supp11]:
    try:
        fn()
    except Exception as e:
        print(f'  ✗ {fn.__name__}: {e}')
print('\n=== Supplementary Figures S1–S11 saved to panels_v3_supp/ ===')
