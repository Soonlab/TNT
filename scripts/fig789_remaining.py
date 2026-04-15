"""
Figures 7, 8, 9 — External validation, HLA/Neoantigen, Clonal evolution
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import Rectangle, FancyBboxPatch, Wedge, Circle
import matplotlib.gridspec as gridspec

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
hla = pd.read_csv(ROOT/'03_hla/hla_class_I_typing.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
ext_meta = pd.read_csv(ROOT/'11_external_validation/meta_analysis_manual.tsv', sep='\t')
ext_stats = pd.read_csv(ROOT/'11_external_validation/signature_stats_manual.tsv', sep='\t')
pyclone = pd.read_csv(ROOT/'04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\t')

UNMATCHED = [13,15,16,17,18,19,33]

# =================================================================
# FIGURE 7 — External validation
# =================================================================
# 7A — Forest plot per cohort × signature with 95% CI
fig, axes = plt.subplots(1, 4, figsize=(15, 5), sharey=True)
sig_names = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
sig_titles = ['DSB/HDR repair','E2F/MYC cell cycle','CD8 proliferation','EMT']
expected = {'DSB_HDR_repair':1, 'E2F_MYC_cellcycle':1, 'CD8_proliferation':1, 'EMT':-1}

for ax, sig, title in zip(axes, sig_names, sig_titles):
    sub = ext_stats[ext_stats.signature==sig].sort_values('delta').reset_index(drop=True)
    # Bootstrap CI from delta - approximate via Hedges' g standard error
    sub['n_total'] = sub.n_good + sub.n_bad
    sub['SE'] = np.sqrt((1/sub.n_good + 1/sub.n_bad))  # approx
    sub['CI_low'] = sub.delta - 1.96*sub.SE
    sub['CI_high'] = sub.delta + 1.96*sub.SE
    y = np.arange(len(sub))
    # Color by direction match
    dir_match = np.sign(sub.delta) == expected[sig]
    colors = [GOOD if d>0 else BAD for d in sub.delta]
    # Draw CI lines
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.plot([r.CI_low, r.CI_high], [i, i], color=colors[i], lw=1.5, alpha=0.7)
        ax.scatter(r.delta, i, s=70+200*np.sqrt(r.n_total/100), color=colors[i],
                   edgecolor='#1d3557', linewidth=0.8, zorder=3, alpha=0.92)
    ax.axvline(0, color='#1d3557', lw=0.8)
    # Meta-analysis diamond at bottom
    meta_z = ext_meta[ext_meta.signature==sig].iloc[0]
    meta_delta = sub['delta'].mean()  # approximate
    meta_p = meta_z.p_meta_onesided
    ax.scatter(meta_delta, len(sub)+0.6, marker='D', s=200, color='#1d3557', edgecolor='#1d3557', zorder=4)
    ax.text(meta_delta, len(sub)+1.2, f'meta Z = {meta_z.Z_stouffer:.2f}\np = {meta_p:.2f}',
            ha='center', fontsize=8, color='#1d3557', fontweight='bold')

    ax.set_yticks(list(y) + [len(sub)+0.6])
    ax.set_yticklabels(list(sub.gse) + ['META-ANALYSIS'], fontsize=8.5)
    ax.set_xlabel('Δ z-score (good − bad)')
    ax.set_title(f'{title}\nexpected: {"↑good" if expected[sig]>0 else "↑bad"}', fontsize=10)
    add_axis_spines(ax)
    ax.set_ylim(-0.6, len(sub)+1.8)

fig.suptitle('External validation — 7 GEO LARC/CRC cohorts (n = 290) with meta-analysis',
             fontsize=12, fontweight='bold', y=1.0, color='#1d3557')
fig.tight_layout()
save_panel(fig, 'Fig7A_forest_meta', OUT)

# 7B — Effect size heatmap (cohort × signature)
fig, ax = plt.subplots(figsize=(7, 5))
heatmap_data = ext_stats.pivot_table(index='gse', columns='signature', values='delta')
heatmap_data = heatmap_data[sig_names]
heatmap_data.columns = sig_titles
heatmap_p = ext_stats.pivot_table(index='gse', columns='signature', values='pvalue')
heatmap_p = heatmap_p[sig_names]; heatmap_p.columns = sig_titles
sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, ax=ax, annot=heatmap_p.values,
            fmt='.2f', annot_kws={'fontsize':8.5, 'color':'#1d3557', 'fontweight':'bold'},
            vmin=-0.6, vmax=0.6, linewidths=1, linecolor='white',
            cbar_kws={'label':'Δ (good − bad)','shrink':0.7})
# Annotate which cohorts agree with discovery (for each sig)
for i, gse in enumerate(heatmap_data.index):
    for j, sig_t in enumerate(heatmap_data.columns):
        d = heatmap_data.iloc[i, j]
        sig_orig = sig_names[j]
        match = np.sign(d) == expected[sig_orig]
        if abs(d) > 0.15:
            ax.text(j+0.5, i+0.85, '✓' if match else '✗', ha='center', va='center',
                    fontsize=12, color='#0f8b78' if match else '#d62828', fontweight='bold')
ax.set_title('Per-cohort effect sizes (annotated = p-value)\n✓ matches discovery direction')
save_panel(fig, 'Fig7B_effect_heatmap', OUT)

# 7C — Stouffer Z meta summary
fig, ax = plt.subplots(figsize=(7, 4))
y_pos = np.arange(len(ext_meta))
colors = [GOOD if z>0 else BAD for z in ext_meta.Z_stouffer]
ax.barh(y_pos, ext_meta.Z_stouffer, color=colors, alpha=0.85, edgecolor='#1d3557', linewidth=0.7, height=0.55)
for i, (_, r) in enumerate(ext_meta.iterrows()):
    ax.text(r.Z_stouffer + (0.08 if r.Z_stouffer>0 else -0.08), i,
            f'Z = {r.Z_stouffer:+.2f},  p = {r.p_meta_onesided:.2f}',
            ha='left' if r.Z_stouffer>0 else 'right', va='center', fontsize=9, fontweight='bold')
ax.axvline(0, color='#1d3557', lw=0.9)
ax.axvline(1.96, color='red', ls='--', lw=0.8, alpha=0.5)
ax.axvline(-1.96, color='red', ls='--', lw=0.8, alpha=0.5)
ax.set_yticks(y_pos); ax.set_yticklabels([s.replace('_',' ') for s in ext_meta.signature])
ax.set_xlabel('Stouffer Z (in expected direction)')
ax.set_xlim(-3, 3)
ax.set_title('Meta-analysis Z scores\n(positive = supports discovery direction)')
add_axis_spines(ax)
save_panel(fig, 'Fig7C_meta_Zscore', OUT)

# =================================================================
# FIGURE 8 — HLA & neoantigen
# =================================================================
# 8A — HLA allele frequency, side-by-side per locus
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
locus_colors = {'A':'#118ab2','B':'#1d3557','C':'#06aed5'}
for ax, locus in zip(axes, ['A','B','C']):
    freq = pd.concat([hla[f'{locus}1'], hla[f'{locus}2']]).value_counts().head(10)
    y = np.arange(len(freq))
    bars = ax.barh(y, freq.values, color=locus_colors[locus], edgecolor='white', linewidth=1.2, alpha=0.92)
    for bar, v in zip(bars, freq.values):
        ax.text(v+0.3, bar.get_y()+bar.get_height()/2, str(int(v)), va='center', fontsize=8.5,
                color='#1d3557', fontweight='bold')
    ax.set_yticks(y); ax.set_yticklabels([f.replace('HLA-','') for f in freq.index], fontsize=9)
    ax.set_xlabel('Allele count'); ax.set_title(f'HLA-{locus}\n(top 10 alleles)', fontsize=11)
    add_axis_spines(ax)
fig.suptitle('HLA class I allele frequency in 35 Korean LARC patients',
             fontsize=12, fontweight='bold', y=1.02)
fig.tight_layout()
save_panel(fig, 'Fig8A_HLA_alleles', OUT)

# 8B — HLA homozygosity raincloud
fig, ax = plt.subplots(figsize=(4.5, 4))
raincloud(ax, hla, 'response_bin', 'n_homozygous_loci', ['good','bad'], PAL_RESP)
g = hla[hla.response_bin=='good'].n_homozygous_loci; b = hla[hla.response_bin=='bad'].n_homozygous_loci
p = stats.mannwhitneyu(g, b).pvalue
stat_bracket(ax, 0, 1, 2.5, p, h=0.15)
ax.set_ylabel('# homozygous class I loci (0–3)')
ax.set_xlabel('')
ax.set_title('HLA homozygosity by response')
ax.set_yticks([0,1,2,3])
save_panel(fig, 'Fig8B_HLA_homozygosity', OUT)

# 8C — HLA LOH bar by response with subject-level breakdown
fig, ax = plt.subplots(figsize=(8, 4))
loh_per = loh.groupby(['subject_id','sample'])['LOH_call'].sum().reset_index()
loh_per = loh_per.merge(wes_inv[['sample_id','timepoint','response_bin']].rename(columns={'sample_id':'sample'}), on='sample', how='left')
# Pre LOH only
loh_pre = loh_per[loh_per.timepoint=='pre'].copy()
loh_pre = loh_pre.sort_values(['response_bin','LOH_call'], ascending=[True, False]).reset_index(drop=True)
colors = [PAL_RESP[r] for r in loh_pre.response_bin]
bars = ax.bar(range(len(loh_pre)), loh_pre.LOH_call, color=colors, edgecolor='white', linewidth=1.2)
ax.set_xticks(range(len(loh_pre)))
ax.set_xticklabels([f'S{int(s)}' for s in loh_pre.subject_id], fontsize=8, rotation=90)
ax.set_ylabel('# HLA class I loci with LOH (pre)')
ax.set_yticks([0,1,2,3])
# Annotate notable
for i, (_, r) in enumerate(loh_pre.iterrows()):
    if r.LOH_call >= 2:
        ax.text(i, r.LOH_call+0.05, '★', ha='center', fontsize=12, color='#d62828')
fisher_p = stats.fisher_exact(pd.crosstab(loh_pre.response_bin, loh_pre.LOH_call.astype(bool))).pvalue
ax.set_title(f'HLA class I LOH (pre-treatment, matched samples)  Fisher p = {fisher_p:.2f}')
add_axis_spines(ax)
fig.legend(handles=[mpatches.Patch(color=GOOD, label='Good'), mpatches.Patch(color=BAD, label='Poor')],
           loc='upper right', fontsize=9, bbox_to_anchor=(0.97, 0.92), frameon=False)
save_panel(fig, 'Fig8C_HLA_LOH', OUT)

# 8D — Pre-treatment neoantigen rain cloud
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
neo_pre_m = neo[(neo.timepoint=='pre') & neo.matched]
for ax, col, title, ylabel in [
    (axes[0], 'n_sites_with_binder', 'Mutation sites with MHC-I binder', '# sites (<500 nM)'),
    (axes[1], 'n_strong_binders_50nM', 'Strong MHC-I binders', '# strong binders (<50 nM)'),
    (axes[2], 'PCN_score', 'Presentation-competent neoantigens', 'PCN score (×HLA LOH adj.)'),
]:
    raincloud(ax, neo_pre_m, 'response_bin', col, ['good','bad'], PAL_RESP)
    g = neo_pre_m[neo_pre_m.response_bin=='good'][col]
    b = neo_pre_m[neo_pre_m.response_bin=='bad'][col]
    if len(g)>=3 and len(b)>=3:
        p = stats.mannwhitneyu(g,b).pvalue
        stat_bracket(ax, 0, 1, max(g.max(), b.max())*1.05, p, h=max(g.max(), b.max())*0.04)
    ax.set_ylabel(ylabel); ax.set_xlabel(''); ax.set_title(title, fontsize=10.5)
fig.suptitle('Pre-treatment neoantigen burden by response', fontsize=12, fontweight='bold', y=1.02)
fig.tight_layout()
save_panel(fig, 'Fig8D_neoantigen_pre', OUT)

# 8E — Pre vs Post neoantigen scatter with arrows per subject
fig, ax = plt.subplots(figsize=(6, 5.5))
neo_pre_d = neo[(neo.timepoint=='pre') & neo.matched].set_index('subject_id')['n_sites_with_binder']
neo_post_d = neo[(neo.timepoint=='post') & neo.matched].set_index('subject_id')['n_sites_with_binder']
common = list(set(neo_pre_d.index) & set(neo_post_d.index))
for s in common:
    resp = clin[clin.subject_id==s].response_bin.iloc[0]
    color = PAL_RESP[resp]
    ax.annotate('', xy=(neo_post_d[s], s), xytext=(neo_pre_d[s], s),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.5, mutation_scale=14, alpha=0.85))
    # Subject label
    ax.text(max(neo_pre_d[s], neo_post_d[s])+5, s, f'S{int(s)}', va='center', fontsize=8, color='#1d3557')
ax.set_xlabel('# sites with MHC-I binder')
ax.set_ylabel('Subject ID')
ax.set_title('Pre → Post neoantigen sites per subject\n(arrow = treatment-induced change)')
fig.legend(handles=[mpatches.Patch(color=GOOD, label='Good'), mpatches.Patch(color=BAD, label='Poor')],
           loc='upper right', fontsize=9, bbox_to_anchor=(0.97, 0.92))
add_axis_spines(ax)
save_panel(fig, 'Fig8E_neoantigen_arrows', OUT)

# 8F — Lollipop of subject neoantigen counts
fig, ax = plt.subplots(figsize=(8, 5))
np_pre = neo[(neo.timepoint=='pre') & neo.matched].sort_values(['response_bin','n_sites_with_binder'], ascending=[True, False])
y = np.arange(len(np_pre))
colors = [PAL_RESP[r] for r in np_pre.response_bin]
for i, (_, r) in enumerate(np_pre.iterrows()):
    ax.plot([0, r.n_sites_with_binder], [i, i], color=colors[i], lw=1.3, alpha=0.65)
ax.scatter(np_pre.n_sites_with_binder, y, s=120, c=colors, edgecolor='#1d3557', linewidth=0.8, alpha=0.92, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([f'S{int(r.subject_id)} ({r.timepoint})' for _, r in np_pre.iterrows()], fontsize=8)
ax.set_xlabel('# mutation sites with MHC-I binder (pre, matched)')
ax.set_title('Per-subject neoantigen burden')
add_axis_spines(ax)
save_panel(fig, 'Fig8F_neoantigen_lollipop', OUT)

# =================================================================
# FIGURE 9 — Clonal evolution (PyClone-VI)
# =================================================================
# 9A — Per-subject cluster summary (small multiples grid)
n_subj = len(pyclone)
fig, axes = plt.subplots(3, 4, figsize=(13, 8))
pyclone_sorted = pyclone.sort_values(['response','subject_id']).reset_index(drop=True)
for ax, (_, r) in zip(axes.flat, pyclone_sorted.iterrows()):
    color = PAL_RESP[r.response]
    # Vertical bar showing n_shrinking, n_expanding, n_stable
    n_total = r.n_clusters
    n_shr = r.n_shrinking; n_exp = r.n_expanding; n_stable = n_total - n_shr - n_exp
    cats = ['Shrinking','Stable','Expanding']
    vals = [n_shr, n_stable, n_exp]
    cols = ['#118ab2','#8d99ae','#d62828']
    bars = ax.bar(cats, vals, color=cols, edgecolor='white', linewidth=1.2, alpha=0.9)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.05, str(int(v)), ha='center', fontsize=9, fontweight='bold')
    ax.set_ylim(0, max(pyclone.n_clusters)+0.8)
    ax.set_title(f'S{int(r.subject_id)} ({r.response})\nn={int(r.n_muts)} mut · k={int(r.n_clusters)} clusters\nshrink Δmin={r.dominant_shrink:+.2f}',
                 fontsize=9, color=color, fontweight='bold')
    ax.set_ylabel('Cluster count' if ax in [axes[0,0], axes[1,0], axes[2,0]] else '')
    add_axis_spines(ax)
    ax.tick_params(labelsize=8)
# Hide unused
for k in range(len(pyclone_sorted), 12):
    axes.flat[k].set_visible(False)
fig.suptitle('PyClone-VI clonal cluster dynamics — 12 paired pre+post subjects',
             fontsize=12, fontweight='bold', y=1.0)
fig.tight_layout()
save_panel(fig, 'Fig9A_pyclone_grid', OUT)

# 9B — Dominant shrinkage raincloud
fig, ax = plt.subplots(figsize=(4.5, 4))
raincloud(ax, pyclone, 'response', 'dominant_shrink', ['good','bad'], PAL_RESP)
g = pyclone[pyclone.response=='good'].dominant_shrink; b = pyclone[pyclone.response=='bad'].dominant_shrink
p = stats.mannwhitneyu(g, b).pvalue
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.05, p, h=0.04)
ax.axhline(0, color='#1d3557', ls='--', lw=0.6)
ax.set_ylabel('Min Δ cellular prevalence (post − pre)')
ax.set_xlabel('')
ax.set_title('Dominant clone shrinkage by response\n(more negative = more clearance)')
save_panel(fig, 'Fig9B_dominant_shrink', OUT)

# 9C — Shrinking vs expanding scatter
fig, ax = plt.subplots(figsize=(5.5, 5))
for resp in ['good','bad']:
    sub = pyclone[pyclone.response==resp]
    ax.scatter(sub.n_shrinking, sub.n_expanding, color=PAL_RESP[resp], s=180, alpha=0.85,
               edgecolor='#1d3557', linewidth=1, label=f'{resp} (n={len(sub)})', zorder=3)
    for _, r in sub.iterrows():
        ax.annotate(f'S{int(r.subject_id)}', (r.n_shrinking, r.n_expanding),
                    xytext=(5, 5), textcoords='offset points', fontsize=8.5, color='#1d3557')
ax.plot([0, 4], [0, 4], color='#6c757d', ls='--', lw=0.6, alpha=0.5)
ax.set_xlabel('Shrinking clusters (Δ CP < −0.2)')
ax.set_ylabel('Expanding clusters (Δ CP > +0.2)')
ax.set_xlim(-0.3, 3); ax.set_ylim(-0.3, 2)
ax.set_xticks([0,1,2,3]); ax.set_yticks([0,1,2])
ax.set_title('Per-subject clonal evolution pattern')
ax.legend(fontsize=9, loc='upper right')
add_axis_spines(ax)
save_panel(fig, 'Fig9C_clonal_scatter', OUT)

print('\n=== Figs 7, 8, 9 (all panels) complete ===')
