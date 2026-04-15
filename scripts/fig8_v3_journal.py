"""
Figure 8 — HLA class I & neoantigen landscape (35 Korean LARC patients)
Journal-style motifs:
  8A  Per-locus HLA-A/B/C allele frequency small multiples (Chowell Science 2018 Fig 1B-C)
  8B  HLA homozygosity violin+jitter by response (Litchfield Cell 2021 Fig 2)
  8C  Per-subject HLA LOH tile heatmap with response sidebar (McGranahan Cell 2017 Fig 2A/3A)
  8D  Pre-treatment neoantigen burden, 3 metrics (Łuksza Nature 2017 Fig 2A-B)
  8E  Pre→Post paired slope per subject (Anagnostou Cancer Discovery 2017 Fig 2)
  8F  Per-subject neoantigen lollipop with HLA allele annotation (Marty Cell 2017 Fig 1D/4A)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'; GOLD_DEEP = '#d4a300'
GRAY_DEEP = '#5a6772'; BAND = '#ecedef'
PAL_RESP_DEEP = {'good': GOOD_DEEP, 'bad': BAD_DEEP}
LOCUS_COL = {'HLA-A': '#118ab2', 'HLA-B': '#1d3557', 'HLA-C': '#06aed5'}
LOCUS_COL_LIGHT = {'HLA-A': '#9fc8d8', 'HLA-B': '#7d8ba8', 'HLA-C': '#9cd6e2'}

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
hla = pd.read_csv(ROOT/'03_hla/hla_class_I_typing.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')

# keep germline HLA only (one per subject)
hla_g = hla[hla.is_germline==True].copy()

# ============================================================
# 8A — Per-locus HLA-A/B/C allele frequency (Chowell Fig 1B-C)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
for ax, locus in zip(axes, ['A','B','C']):
    freq = pd.concat([hla_g[f'{locus}1'], hla_g[f'{locus}2']]).value_counts().head(10)
    # homozygous count per allele
    hom_mask = hla_g[f'{locus}1']==hla_g[f'{locus}2']
    hom_counts = hla_g[hom_mask][f'{locus}1'].value_counts()
    y = np.arange(len(freq))[::-1]
    total = freq.values
    hom_vals = np.array([hom_counts.get(a,0)*2 for a in freq.index])
    het_vals = total - hom_vals
    base_col = LOCUS_COL[f'HLA-{locus}']
    ax.barh(y, het_vals, color=base_col, edgecolor='white', linewidth=1.0,
            alpha=0.92, label='heterozygous allele count')
    ax.barh(y, hom_vals, left=het_vals, color=GOLD_DEEP, edgecolor='white',
            linewidth=1.0, alpha=0.92, label='from homozygous genotype')
    for i, (yy, v) in enumerate(zip(y, total)):
        ax.text(v+0.25, yy, str(int(v)), va='center', fontsize=8.5,
                color=BLACK_DEEP, fontweight='bold')
    ax.set_yticks(y)
    ax.set_yticklabels([a.replace('HLA-','') for a in freq.index], fontsize=9,
                       color=BLACK_DEEP)
    ax.set_xlabel('Allele count (2 × n_patients)', fontsize=9.5, color=BLACK_DEEP)
    ax.text(0.02, 1.02, f'HLA-{locus}  (top 10)', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=11, fontweight='bold',
            color=LOCUS_COL[f'HLA-{locus}'])
    ax.set_xlim(0, max(total)*1.18)
    add_axis_spines(ax); ax.tick_params(labelsize=8.5)
axes[-1].legend(loc='lower right', fontsize=8, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig8A_HLA_alleles', OUT)

# ============================================================
# 8B — HLA homozygosity discrete dot-stack by response
# (n_homozygous_loci is integer 0–3; stack patients per level/group)
# ============================================================
fig, ax = plt.subplots(figsize=(5.4, 4.3))
order = ['good','bad']
x_offset = {'good': -0.18, 'bad': +0.18}
for g in order:
    vals = hla_g[hla_g.response_bin==g].n_homozygous_loci.values
    col = PAL_RESP_DEEP[g]
    # count per level
    for lvl in [0,1,2,3]:
        n = int((vals==lvl).sum())
        if n==0: continue
        # stacked dots
        for k in range(n):
            ax.scatter(lvl + x_offset[g], k+0.5, s=120, color=col,
                       edgecolor=BLACK_DEEP, linewidth=0.8, alpha=0.92, zorder=3)
        ax.text(lvl + x_offset[g], n+0.45, str(n), ha='center', fontsize=9,
                color=col, fontweight='bold')
    # mean line
    mean_v = np.mean(vals)
    ax.axvline(mean_v, color=col, ls='--', lw=1.1, alpha=0.55, zorder=1)

g_v = hla_g[hla_g.response_bin=='good'].n_homozygous_loci.values
b_v = hla_g[hla_g.response_bin=='bad'].n_homozygous_loci.values
p = stats.mannwhitneyu(g_v, b_v).pvalue
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['0','1','2','3'], fontsize=11, color=BLACK_DEEP)
ax.tick_params(axis='x', pad=8)
ax.set_xlabel('# homozygous HLA class-I loci', fontsize=10, color=BLACK_DEEP)
ax.set_ylabel('Number of patients', fontsize=10, color=BLACK_DEEP)
ax.set_xlim(-0.6, 3.6)
y_top = max([(g_v==l).sum() for l in [0,1,2,3]] +
            [(b_v==l).sum() for l in [0,1,2,3]]) + 2.2
ax.set_ylim(0, y_top)
ax.text(0.98, 0.98, f'Mann–Whitney p = {p:.2f}\n(good mean={g_v.mean():.2f}, bad mean={b_v.mean():.2f})',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        color=BLACK_DEEP, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=BLACK_DEEP, lw=0.6, alpha=0.9))
leg = [mpatches.Patch(color=GOOD_DEEP, label=f'good (n={len(g_v)})'),
       mpatches.Patch(color=BAD_DEEP,  label=f'bad  (n={len(b_v)})')]
ax.legend(handles=leg, loc='upper right', bbox_to_anchor=(0.98, 0.72),
          fontsize=9, frameon=False)
add_axis_spines(ax); ax.tick_params(labelsize=9)
save_panel(fig, 'Fig8B_HLA_homozygosity', OUT)

# ============================================================
# 8C — Per-subject HLA LOH tile heatmap (McGranahan Fig 2A/3A)
#      rows = subjects (pre samples), cols = A/B/C, cell = LOH status
# ============================================================
loh_pre = loh.merge(
    wes_inv[['sample_id','timepoint','response_bin']].rename(columns={'sample_id':'sample'}),
    on='sample', how='left')
loh_pre = loh_pre[loh_pre.timepoint=='pre'].copy()

# pivot
tile = loh_pre.pivot_table(index='subject_id', columns='locus', values='LOH_call',
                           aggfunc='max').astype(float)
# sort by response then total LOH desc
meta_subj = loh_pre.groupby('subject_id')['response_bin'].first()
tile = tile.assign(_resp=meta_subj, _tot=tile.sum(axis=1))
tile = tile.sort_values(['_resp','_tot'], ascending=[True, False])
resp_col = tile._resp.map(PAL_RESP_DEEP)
tile_vals = tile.drop(columns=['_resp','_tot'])

fig = plt.figure(figsize=(11, 3.6))
gs = gridspec.GridSpec(1, 3, width_ratios=[0.25, 10, 2.6], wspace=0.55)
ax_side = fig.add_subplot(gs[0,0])
ax_tile = fig.add_subplot(gs[0,1])
ax_sum  = fig.add_subplot(gs[0,2])

# response sidebar (rows as x-axis swap: rows=loci, cols=subjects for horizontal readability)
tile_T = tile_vals.T   # loci × subjects
subjects = tile_vals.index.tolist()
loci_ord = ['HLA-A','HLA-B','HLA-C']
tile_T = tile_T.loc[loci_ord, subjects]

# draw tiles
for i, loc_name in enumerate(loci_ord):
    for j, s in enumerate(subjects):
        v = tile_T.loc[loc_name, s]
        col = '#f1f3f5' if (v==0 or np.isnan(v)) else BAD_DEEP
        ax_tile.add_patch(Rectangle((j, i), 0.92, 0.92, facecolor=col,
                                    edgecolor='white', lw=1.0))
        if v==1:
            ax_tile.text(j+0.46, i+0.46, '■', ha='center', va='center',
                         fontsize=11, color='white', fontweight='bold')
# response annotation row above
for j, s in enumerate(subjects):
    r = meta_subj.loc[s]
    ax_tile.add_patch(Rectangle((j, -0.95), 0.92, 0.7, facecolor=PAL_RESP_DEEP[r],
                                edgecolor='white', lw=1.0))
ax_tile.text(-0.2, -0.6, 'response', ha='right', va='center', fontsize=9,
             color=BLACK_DEEP, fontweight='bold')
ax_tile.set_xlim(-0.2, len(subjects)+0.1)
ax_tile.set_ylim(len(loci_ord)+0.1, -1.2)
ax_tile.set_xticks([j+0.46 for j in range(len(subjects))])
ax_tile.set_xticklabels([f'S{int(s)}' for s in subjects], fontsize=8.5, rotation=0,
                        color=BLACK_DEEP)
ax_tile.set_yticks([i+0.46 for i in range(len(loci_ord))])
ax_tile.set_yticklabels(loci_ord, fontsize=10, color=BLACK_DEEP)
for spine in ax_tile.spines.values(): spine.set_visible(False)
ax_tile.tick_params(length=0)
ax_side.axis('off')

# Right summary: LOH counts per response
summary = loh_pre.groupby(['response_bin','LOH_call']).size().unstack(fill_value=0)
# test
tbl = pd.crosstab(loh_pre.response_bin, loh_pre.LOH_call.astype(bool))
try:
    fisher_p = stats.fisher_exact(tbl).pvalue
except Exception:
    fisher_p = np.nan
resp_totals = loh_pre.groupby('response_bin').size()
loh_pct = summary.get(True, pd.Series(0, index=summary.index)) / resp_totals * 100
ypos = np.arange(len(loh_pct))
colors = [PAL_RESP_DEEP[r] for r in loh_pct.index]
ax_sum.barh(ypos, loh_pct.values, color=colors, edgecolor=BLACK_DEEP, lw=0.7,
            height=0.55, alpha=0.92)
for i, (r, v) in enumerate(loh_pct.items()):
    n_loh = int(summary.loc[r].get(True,0))
    n_tot = int(resp_totals[r])
    ax_sum.text(v+1, i, f'{n_loh}/{n_tot}  ({v:.0f}%)', va='center',
                fontsize=9, color=BLACK_DEEP, fontweight='bold')
ax_sum.set_yticks(ypos); ax_sum.set_yticklabels(loh_pct.index, fontsize=10,
                                               color=BLACK_DEEP)
ax_sum.invert_yaxis()
ax_sum.set_xlim(0, max(loh_pct.values)*1.6 + 10)
ax_sum.set_xlabel('% loci with LOH', fontsize=9.5, color=BLACK_DEEP)
ax_sum.set_title(f'Fisher p = {fisher_p:.2f}', fontsize=10, color=BLACK_DEEP,
                 fontweight='bold', loc='left')
add_axis_spines(ax_sum); ax_sum.tick_params(labelsize=9)

# Legend for LOH tile
leg = [mpatches.Patch(color=BAD_DEEP, label='LOH'),
       mpatches.Patch(color='#f1f3f5', ec=GRAY_DEEP, label='retained / biallelic'),
       mpatches.Patch(color=GOOD_DEEP, label='good'),
       mpatches.Patch(color=BAD_DEEP, label='bad')]
ax_tile.legend(handles=leg[:2], loc='upper right', bbox_to_anchor=(1.0, -0.08),
               fontsize=8, frameon=False, ncol=2)
fig.tight_layout()
save_panel(fig, 'Fig8C_HLA_LOH', OUT)

# ============================================================
# 8D — Pre-treatment neoantigen burden, 3 metrics (Łuksza Fig 2)
# ============================================================
neo_pre = neo[(neo.timepoint=='pre') & neo.matched].copy()
fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
panels = [
    ('n_sites_with_binder',   '# mutation sites with MHC-I binder (<500 nM)'),
    ('n_strong_binders_50nM', '# strong MHC-I binders (<50 nM)'),
    ('PCN_score',             'PCN score  (HLA-LOH adjusted)')]
for ax, (col, ylab) in zip(axes, panels):
    for i, g in enumerate(['good','bad']):
        vals = pd.to_numeric(neo_pre[neo_pre.response_bin==g][col],
                             errors='coerce').dropna().values
        c = PAL_RESP_DEEP[g]
        # half violin
        if len(vals) >= 2 and np.std(vals) > 1e-6:
            try:
                kde = stats.gaussian_kde(vals, bw_method=0.55)
                xg = np.linspace(vals.min()-1, vals.max()+1, 120)
                d = kde(xg); d = d/d.max()*0.36
                ax.fill_betweenx(xg, i, i-d, color=c, alpha=0.3, linewidth=0)
                ax.plot(i-d, xg, color=c, lw=1.1, alpha=0.85)
            except Exception:
                pass
        q1, med, q3 = np.percentile(vals, [25,50,75])
        ax.add_patch(Rectangle((i+0.05, q1), 0.14, q3-q1, facecolor='white',
                               edgecolor=c, lw=1.3))
        ax.plot([i+0.05, i+0.19], [med, med], color=c, lw=2.2)
        jit = np.random.uniform(0.26, 0.44, len(vals))
        ax.scatter(i+jit, vals, s=38, color=c, edgecolor=BLACK_DEEP,
                   linewidth=0.6, alpha=0.85, zorder=3)
    g = pd.to_numeric(neo_pre[neo_pre.response_bin=='good'][col], errors='coerce').dropna()
    b = pd.to_numeric(neo_pre[neo_pre.response_bin=='bad'][col], errors='coerce').dropna()
    if len(g)>=3 and len(b)>=3:
        p = stats.mannwhitneyu(g,b).pvalue
        ymax = max(g.max(), b.max())*1.10
        ax.plot([0,0,1,1],[ymax, ymax*1.03, ymax*1.03, ymax], color=BLACK_DEEP, lw=1)
        ax.text(0.5, ymax*1.05, f'p = {p:.2f}', ha='center', fontsize=9,
                color=BLACK_DEEP, fontweight='bold')
    ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'], fontsize=11, color=BLACK_DEEP)
    ax.tick_params(axis='x', pad=10)
    ax.set_ylabel(ylab, fontsize=9.5, color=BLACK_DEEP)
    ax.set_xlim(-0.55, 1.55)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
fig.tight_layout()
save_panel(fig, 'Fig8D_neoantigen_pre', OUT)

# ============================================================
# 8E — Pre→Post paired slope per subject (Anagnostou Cancer Discovery 2017 Fig 2)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(9, 4.6), sharey=True)
for ax, metric, ylab in [
    (axes[0], 'n_sites_with_binder', '# sites with MHC-I binder'),
    (axes[1], 'n_strong_binders_50nM','# strong binders (<50 nM)'),
]:
    pre = neo[(neo.timepoint=='pre') & neo.matched].set_index('subject_id')
    post = neo[(neo.timepoint=='post') & neo.matched].set_index('subject_id')
    common = sorted(set(pre.index) & set(post.index))
    for resp in ['good','bad']:
        subs = [s for s in common if clin[clin.subject_id==s].response_bin.iloc[0]==resp]
        c = PAL_RESP_DEEP[resp]
        for s in subs:
            y0, y1 = pre.loc[s, metric], post.loc[s, metric]
            ax.plot([0, 1], [y0, y1], color=c, lw=1.1, alpha=0.40, zorder=2)
            ax.scatter([0,1], [y0,y1], s=38, color=c, edgecolor=BLACK_DEEP,
                       linewidth=0.5, zorder=3, alpha=0.75)
            ax.text(-0.04, y0, f'S{int(s)}', ha='right', va='center',
                    fontsize=7.5, color=c, alpha=0.75)
        # mean line with CI shading + big anchor dots
        if subs:
            y0s = np.array([pre.loc[s, metric] for s in subs])
            y1s = np.array([post.loc[s, metric] for s in subs])
            mean_pre, mean_post = y0s.mean(), y1s.mean()
            se_pre = y0s.std(ddof=1)/np.sqrt(len(y0s)) if len(y0s)>1 else 0
            se_post = y1s.std(ddof=1)/np.sqrt(len(y1s)) if len(y1s)>1 else 0
            # gradient band (1 SE)
            xs = np.linspace(0,1,30)
            mean_line = mean_pre + (mean_post-mean_pre)*xs
            se_line = se_pre + (se_post-se_pre)*xs
            ax.fill_between(xs, mean_line-se_line, mean_line+se_line,
                            color=c, alpha=0.18, zorder=4)
            # bold mean line
            ax.plot([0,1], [mean_pre, mean_post], color=c, lw=3.2, alpha=1.0,
                    zorder=5, solid_capstyle='round')
            ax.scatter([0,1], [mean_pre, mean_post], s=180, color=c,
                       edgecolor=BLACK_DEEP, linewidth=1.3, zorder=6, alpha=1.0)
            # Δ̄ label — stack vertically with small offset so both labels sit
            # near the post endpoints but never overlap each other
            delta = mean_post - mean_pre
            yt_frac = 0.92 if resp == 'good' else 0.82
            ax.annotate(f'{resp}  Δ̄ = {delta:+.1f}',
                        xy=(1.0, mean_post),
                        xytext=(1.08, yt_frac),
                        textcoords=('data', 'axes fraction'),
                        ha='left', va='center',
                        fontsize=9, color=c, fontweight='bold',
                        arrowprops=dict(arrowstyle='-', color=c, lw=0.8, alpha=0.5),
                        bbox=dict(boxstyle='round,pad=0.22', fc='white',
                                  ec=c, lw=0.8, alpha=0.95))
    # Wilcoxon per group — show as legend box at upper right
    p_lines = []
    for resp in ['good','bad']:
        subs = [s for s in common if clin[clin.subject_id==s].response_bin.iloc[0]==resp]
        if len(subs) >= 3:
            y0 = [pre.loc[s, metric] for s in subs]
            y1 = [post.loc[s, metric] for s in subs]
            try:
                p = stats.wilcoxon(y0, y1).pvalue
            except Exception:
                p = np.nan
            p_lines.append((resp, p, len(subs)))
    if p_lines:
        txt = '\n'.join([f'{r}  (n={n}):  p = {p:.2f}' for r,p,n in p_lines])
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, ha='left', va='top',
                fontsize=8.5, color=BLACK_DEEP, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=BLACK_DEEP,
                          lw=0.6, alpha=0.9))
    ax.set_xticks([0,1]); ax.set_xticklabels(['pre','post'], fontsize=11,
                                              color=BLACK_DEEP, fontweight='bold')
    ax.tick_params(axis='x', pad=10)
    ax.set_ylabel(ylab, fontsize=10, color=BLACK_DEEP)
    ax.set_xlim(-0.25, 1.45)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
fig.tight_layout()
save_panel(fig, 'Fig8E_neoantigen_paired', OUT)

# ============================================================
# 8F — Per-subject neoantigen lollipop with HLA allele annotation
#      (Marty Cell 2017 Fig 1D/4A)
# ============================================================
fig = plt.figure(figsize=(10.5, 6.5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 3.4], wspace=0.02)
ax_hla = fig.add_subplot(gs[0,0])
ax_lol = fig.add_subplot(gs[0,1], sharey=ax_hla)

np_pre = neo[(neo.timepoint=='pre') & neo.matched].copy()
np_pre = np_pre.sort_values(['response_bin','n_sites_with_binder'],
                            ascending=[True, True]).reset_index(drop=True)
y = np.arange(len(np_pre))

# lollipop panel (right)
for i, r in np_pre.iterrows():
    c = PAL_RESP_DEEP[r.response_bin]
    loh_flag = bool(r.HLA_LOH) if 'HLA_LOH' in r else False
    ax_lol.plot([0, r.n_sites_with_binder], [i, i], color=c, lw=1.4, alpha=0.6, zorder=2)
    marker = 'X' if loh_flag else 'o'
    ax_lol.scatter(r.n_sites_with_binder, i, s=150, color=c, edgecolor=BLACK_DEEP,
                   linewidth=0.9, alpha=0.95, zorder=3, marker=marker)
    ax_lol.text(r.n_sites_with_binder + 1.5, i, str(int(r.n_sites_with_binder)),
                va='center', fontsize=8, color=BLACK_DEEP, fontweight='bold')
    # strong binder overlay
    ax_lol.scatter(r.n_strong_binders_50nM, i, s=36, color=GOLD_DEEP,
                   edgecolor=BLACK_DEEP, linewidth=0.6, alpha=0.9, zorder=4)
ax_lol.set_yticks(y)
ax_lol.set_yticklabels([f'S{int(r.subject_id)}' for _, r in np_pre.iterrows()],
                       fontsize=8.5, color=BLACK_DEEP)
ax_lol.set_xlabel('# mutation sites with MHC-I binder  (pre)', fontsize=10,
                  color=BLACK_DEEP)
ax_lol.set_xlim(-2, np_pre.n_sites_with_binder.max()*1.15)
ax_lol.invert_yaxis()
add_axis_spines(ax_lol); ax_lol.tick_params(labelsize=9)

# HLA allele track (left) — homozygosity per locus + response
ax_hla.set_xlim(0, 5.2); ax_hla.set_ylim(len(np_pre)-0.5, -0.5)
# columns: response, A, B, C
col_x = {'resp':0.3, 'HLA-A':1.5, 'HLA-B':2.6, 'HLA-C':3.7}
# headers
for k, xx in col_x.items():
    ax_hla.text(xx, -1, k.replace('HLA-',''), ha='center', va='center',
                fontsize=9.5, fontweight='bold', color=BLACK_DEEP)
ax_hla.text(col_x['resp'], -1.7, 'response', ha='center', va='center',
            fontsize=8.5, color=BLACK_DEEP)
ax_hla.text(2.6, -1.7, 'HLA zygosity', ha='center', va='center',
            fontsize=8.5, color=BLACK_DEEP, fontstyle='italic')
for i, r in np_pre.iterrows():
    # response square
    ax_hla.add_patch(Rectangle((col_x['resp']-0.35, i-0.35), 0.7, 0.7,
                               facecolor=PAL_RESP_DEEP[r.response_bin],
                               edgecolor='white', lw=0.8))
    # A/B/C zygosity
    row_g = hla_g[hla_g.subject_id==r.subject_id]
    if len(row_g)==0: continue
    g = row_g.iloc[0]
    for locus in ['A','B','C']:
        hom = bool(g[f'homozygous_{locus}'])
        xx = col_x[f'HLA-{locus}']
        ax_hla.add_patch(Rectangle((xx-0.35, i-0.35), 0.7, 0.7,
                                   facecolor=GOLD_DEEP if hom else LOCUS_COL_LIGHT[f'HLA-{locus}'],
                                   edgecolor='white', lw=0.8))
        ax_hla.text(xx, i, 'hom' if hom else 'het', ha='center', va='center',
                    fontsize=7, color=BLACK_DEEP,
                    fontweight='bold' if hom else 'normal')
for s in ax_hla.spines.values(): s.set_visible(False)
ax_hla.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

# Legend
handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=GRAY_DEEP,
                  markersize=10, markeredgecolor=BLACK_DEEP, label='HLA retained'),
           Line2D([0],[0], marker='X', color='w', markerfacecolor=GRAY_DEEP,
                  markersize=10, markeredgecolor=BLACK_DEEP, label='HLA LOH'),
           Line2D([0],[0], marker='o', color='w', markerfacecolor=GOLD_DEEP,
                  markersize=7, markeredgecolor=BLACK_DEEP, label='strong binders (<50 nM)'),
           mpatches.Patch(color=GOLD_DEEP, label='homozygous locus')]
ax_lol.legend(handles=handles, loc='center left', bbox_to_anchor=(1.01, 0.5),
              fontsize=8, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig8F_neoantigen_lollipop', OUT)

print('\n=== Fig 8 v3 (HLA + neoantigen, journal motifs) complete ===')
