"""
Figure 7 — External validation (7 GEO LARC/CRC cohorts, n=290)
Journal-style motifs:
  7A  Per-cohort forest + Stouffer pooled diamond (Litchfield Cell 2021 / Chowell Science 2018)
  7B  Cohort × signature effect-size heatmap with ✓/✗ direction agreement (Mariathasan Nature 2018 / Cristescu Science 2018)
  7C  Stouffer Z meta bars, row-aligned to 7B (Thorsson Immunity 2018 marginal bar motif; Bareche Ann Oncol 2022 combined panel)
  7D  Per-cohort concordance score (fraction signatures agreeing with discovery)
  7E  Funnel-style |Δ| vs precision (1/SE) per signature for heterogeneity view
  7F  Discovery vs validation direction dot-summary
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)

# Deep saturated palette to match Fig 1/2/5 polish
GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'; GOLD_DEEP = '#d4a300'
GRAY_DEEP = '#5a6772'; BAND = '#ecedef'
PAL_RESP_DEEP = {'good': GOOD_DEEP, 'bad': BAD_DEEP}

ext_meta = pd.read_csv(ROOT/'11_external_validation/meta_analysis_manual.tsv', sep='\t')
ext_stats = pd.read_csv(ROOT/'11_external_validation/signature_stats_manual.tsv', sep='\t')

sig_names  = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
sig_titles = ['DSB/HDR repair','E2F/MYC cell cycle','CD8 proliferation','EMT']
expected   = {'DSB_HDR_repair':1, 'E2F_MYC_cellcycle':1, 'CD8_proliferation':1, 'EMT':-1}

# Pre-compute CIs
ext_stats = ext_stats.copy()
ext_stats['n_total'] = ext_stats.n_good + ext_stats.n_bad
ext_stats['SE'] = np.sqrt(1/ext_stats.n_good + 1/ext_stats.n_bad)
ext_stats['CI_low']  = ext_stats.delta - 1.96*ext_stats.SE
ext_stats['CI_high'] = ext_stats.delta + 1.96*ext_stats.SE

# =========================================================================
# 7A — Per-cohort forest with Stouffer pooled diamond + heterogeneity band
# =========================================================================
fig, axes = plt.subplots(1, 4, figsize=(15.5, 5.2), sharey=True)
for ax, sig, title in zip(axes, sig_names, sig_titles):
    sub = ext_stats[ext_stats.signature==sig].sort_values('delta').reset_index(drop=True)
    y = np.arange(len(sub))
    exp_dir = expected[sig]
    # Background band
    for i in range(len(sub)):
        if i % 2 == 0:
            ax.axhspan(i-0.5, i+0.5, color=BAND, alpha=0.55, zorder=0)
    # CI + point
    for i, r in sub.iterrows():
        match = np.sign(r.delta) == exp_dir
        col = GOOD_DEEP if match else BAD_DEEP
        ax.plot([r.CI_low, r.CI_high], [i, i], color=col, lw=1.8, alpha=0.85, zorder=2,
                solid_capstyle='round')
        ax.scatter(r.delta, i, s=30 + 4.5*r.n_total, color=col, edgecolor=BLACK_DEEP,
                   linewidth=0.9, zorder=3, alpha=0.95)
    ax.axvline(0, color=BLACK_DEEP, lw=1.0, zorder=1)
    # Meta pooled diamond (Samstein shaded band trick)
    m = ext_meta[ext_meta.signature==sig].iloc[0]
    # Inverse-variance weighted mean Δ
    w = 1/sub.SE**2
    meta_delta = float((sub.delta*w).sum()/w.sum())
    meta_se = float(np.sqrt(1/w.sum()))
    meta_lo, meta_hi = meta_delta - 1.96*meta_se, meta_delta + 1.96*meta_se
    # Heterogeneity Q and I²
    Q = float(((sub.delta - meta_delta)**2 * w).sum()); df = len(sub)-1
    I2 = max(0.0, (Q-df)/Q) if Q>0 else 0.0
    yd = len(sub)+0.9
    ax.axhspan(yd-0.55, yd+0.55, xmin=0, xmax=1, color='#ecedef', alpha=0.5, zorder=0)
    diamond = Polygon([(meta_lo, yd), (meta_delta, yd+0.35), (meta_hi, yd), (meta_delta, yd-0.35)],
                      closed=True, facecolor=BLACK_DEEP, edgecolor=BLACK_DEEP, zorder=4)
    ax.add_patch(diamond)
    # Meta stats text
    agree = '✓' if np.sign(m.Z_stouffer)==exp_dir else '✗'
    col_meta = GOOD_DEEP if agree=='✓' else BAD_DEEP
    ax.text(0.98, 1.18, f'Stouffer Z = {m.Z_stouffer:+.2f}   p = {m.p_meta_onesided:.2f}\nI² = {I2*100:.0f}%   cohorts = {len(sub)}   N = {int(sub.n_total.sum())}   {agree}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5,
            color=col_meta, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=col_meta, lw=0.8, alpha=0.95))
    ax.set_yticks(list(y)+[yd])
    ax.set_yticklabels(list(sub.gse)+['pooled (Stouffer)'], fontsize=8.5, color=BLACK_DEEP)
    ax.set_xlabel('Δ z-score  (good − bad)', fontsize=10, color=BLACK_DEEP)
    exp_arrow = '↑ good' if exp_dir>0 else '↑ bad'
    ax.text(0.02, 1.02, f'{title}   (expected: {exp_arrow})', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=10, fontweight='bold', color=BLACK_DEEP)
    ax.set_ylim(-0.7, yd+0.9)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
# Cohort-size legend on last ax
size_legend = [Line2D([0],[0], marker='o', ls='', markersize=np.sqrt(30+4.5*n)/2.2,
                      color=GRAY_DEEP, markeredgecolor=BLACK_DEEP, label=f'n = {n}')
               for n in [30, 60, 120]]
axes[-1].legend(handles=size_legend, title='cohort size', loc='lower right',
                fontsize=8, title_fontsize=8.5, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.88))
save_panel(fig, 'Fig7A_forest_meta', OUT)

# =========================================================================
# 7B — Cohort × signature effect-size heatmap with direction-agreement borders
# =========================================================================
fig = plt.figure(figsize=(8.5, 5.2))
gs = gridspec.GridSpec(1, 2, width_ratios=[7.0, 0.32], wspace=0.05)
ax_hm = fig.add_subplot(gs[0,0])
ax_cb = fig.add_subplot(gs[0,1])

hm = ext_stats.pivot_table(index='gse', columns='signature', values='delta')[sig_names]
hm_p = ext_stats.pivot_table(index='gse', columns='signature', values='pvalue')[sig_names]
hm_n = ext_stats.pivot_table(index='gse', columns='signature', values='n_total')[sig_names]
hm.columns = sig_titles; hm_p.columns = sig_titles; hm_n.columns = sig_titles

vmax = 0.6
im = ax_hm.imshow(hm.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
# Grid
for i in range(hm.shape[0]+1):
    ax_hm.axhline(i-0.5, color='white', lw=1.5)
for j in range(hm.shape[1]+1):
    ax_hm.axvline(j-0.5, color='white', lw=1.5)
# Cell annotations: value + p + ✓/✗ border (Cristescu tick motif)
for i, gse in enumerate(hm.index):
    for j, sig_t in enumerate(hm.columns):
        d = hm.iloc[i,j]; p = hm_p.iloc[i,j]
        sig_orig = sig_names[j]
        match = np.sign(d) == expected[sig_orig]
        star = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else ''))
        txt_col = 'white' if abs(d)>0.35 else BLACK_DEEP
        ax_hm.text(j, i-0.12, f'{d:+.2f}', ha='center', va='center', fontsize=9,
                   color=txt_col, fontweight='bold')
        ax_hm.text(j, i+0.22, f'p={p:.2f}{star}', ha='center', va='center', fontsize=7.5,
                   color=txt_col)
        if match and abs(d)>0.05:
            ax_hm.add_patch(Rectangle((j-0.48, i-0.48), 0.96, 0.96, fill=False,
                                      edgecolor=GOOD_DEEP, lw=2.2, zorder=5))
ax_hm.set_xticks(range(len(hm.columns)))
ax_hm.set_xticklabels(hm.columns, rotation=18, ha='right', fontsize=10, color=BLACK_DEEP)
ax_hm.set_yticks(range(len(hm.index)))
ax_hm.set_yticklabels([f'{g}  (n={int(hm_n.loc[g].iloc[0])})' for g in hm.index],
                      fontsize=9.5, color=BLACK_DEEP)
ax_hm.tick_params(length=0)
for spine in ax_hm.spines.values(): spine.set_visible(False)
# Colorbar (left narrow column)
cb = plt.colorbar(im, cax=ax_cb)
cb.set_label('Δ z-score (good − bad)', fontsize=9.5, color=BLACK_DEEP)
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor(BLACK_DEEP)
fig.tight_layout()
save_panel(fig, 'Fig7B_heatmap', OUT)

# =========================================================================
# 7C — Stouffer Z meta bars per signature (standalone)
# =========================================================================
fig, ax_bar = plt.subplots(figsize=(7.2, 3.8))
meta = ext_meta.set_index('signature').loc[sig_names].reset_index()
y = np.arange(len(meta))
colors = [GOOD_DEEP if np.sign(z)==e else BAD_DEEP for z,e in zip(meta.Z_stouffer, meta.expected_dir)]
ax_bar.barh(y, meta.Z_stouffer, color=colors, edgecolor=BLACK_DEEP, linewidth=0.8,
            height=0.62, alpha=0.9)
ax_bar.axvline(0, color=BLACK_DEEP, lw=1.0)
ax_bar.axvspan(-1.96, 1.96, color='#f2f3f5', alpha=0.65, zorder=0)
ax_bar.axvline( 1.96, color=GRAY_DEEP, ls='--', lw=0.8)
ax_bar.axvline(-1.96, color=GRAY_DEEP, ls='--', lw=0.8)
LABEL_X = 3.35
for i, r in meta.iterrows():
    bar_end = r.Z_stouffer
    tip_x = bar_end + (0.08 if bar_end>=0 else -0.08)
    ax_bar.plot([tip_x, LABEL_X-0.05], [i, i], color=GRAY_DEEP, lw=0.6, ls=':', alpha=0.7)
    ax_bar.text(LABEL_X, i, f'Z = {r.Z_stouffer:+.2f}   p = {r.p_meta_onesided:.2f}',
                ha='left', va='center', fontsize=9.5, color=BLACK_DEEP, fontweight='bold')
ax_bar.set_yticks(y)
ax_bar.set_yticklabels(sig_titles, fontsize=10, color=BLACK_DEEP)
ax_bar.invert_yaxis()
ax_bar.set_xlim(-3.2, 6.4)
ax_bar.set_xticks([-3,-1.96,0,1.96,3])
ax_bar.set_xticklabels(['-3','-1.96','0','1.96','3'], fontsize=9)
ax_bar.set_xlabel('Stouffer Z  (positive = supports discovery direction)',
                  fontsize=10, color=BLACK_DEEP)
ax_bar.text(1.96, -0.85, '|Z|=1.96', color=GRAY_DEEP, fontsize=7.5, ha='center')
add_axis_spines(ax_bar); ax_bar.tick_params(labelsize=9)
ax_bar.spines['right'].set_visible(False); ax_bar.spines['top'].set_visible(False)
fig.tight_layout()
save_panel(fig, 'Fig7C_meta_Zscore', OUT)

# =========================================================================
# 7D — Per-cohort concordance score  (fraction of signatures agreeing)
# =========================================================================
fig, ax = plt.subplots(figsize=(7, 4.2))
concord = []
for gse, g in ext_stats.groupby('gse'):
    n_match = sum(np.sign(r.delta)==expected[r.signature] for _, r in g.iterrows())
    concord.append({'gse':gse, 'frac':n_match/len(g), 'n_match':n_match,
                    'n_total':len(g), 'patients':int(g.n_total.iloc[0])})
cdf = pd.DataFrame(concord).sort_values('frac', ascending=True).reset_index(drop=True)
y = np.arange(len(cdf))
bar_colors = [GOOD_DEEP if f>=0.5 else BAD_DEEP for f in cdf.frac]
bars = ax.barh(y, cdf.frac*100, color=bar_colors, edgecolor=BLACK_DEEP, linewidth=0.7,
               height=0.65, alpha=0.92)
ax.axvline(50, color=GRAY_DEEP, ls='--', lw=0.8)
for i, r in cdf.iterrows():
    ax.text(r.frac*100 + 1.2, i, f'{r.n_match}/{r.n_total}  (n={r.patients})',
            va='center', fontsize=9, color=BLACK_DEEP, fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(cdf.gse, fontsize=9.5, color=BLACK_DEEP)
ax.set_xlabel('Signatures matching discovery direction (%)', fontsize=10, color=BLACK_DEEP)
ax.set_xlim(0, 115)
ax.set_xticks([0,25,50,75,100])
add_axis_spines(ax); ax.tick_params(labelsize=9)
save_panel(fig, 'Fig7D_cohort_concordance', OUT)

# =========================================================================
# 7E — Funnel plot per signature (2×2 small multiples)
#      x = effect (Δ), y = SE (inverted). 95% CI wedge for null=0.
# =========================================================================
fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
se_max = float(ext_stats.SE.max())*1.15
for ax, sig, title in zip(axes.flat, sig_names, sig_titles):
    sub = ext_stats[ext_stats.signature==sig].copy()
    exp_dir = expected[sig]
    # 95% pseudo-CI wedge around null=0: x = ±1.96 * SE
    se_grid = np.linspace(0, se_max, 50)
    ax.fill_betweenx(se_grid, -1.96*se_grid, 1.96*se_grid,
                     color=GRAY_DEEP, alpha=0.10, zorder=0)
    ax.plot(-1.96*se_grid, se_grid, color=GRAY_DEEP, ls='--', lw=0.9)
    ax.plot( 1.96*se_grid, se_grid, color=GRAY_DEEP, ls='--', lw=0.9)
    ax.axvline(0, color=BLACK_DEEP, lw=1.0, zorder=1)
    # Inverse-variance meta mean
    w = 1/sub.SE**2
    meta_delta = float((sub.delta*w).sum()/w.sum())
    ax.axvline(meta_delta, color=BLACK_DEEP, ls=':', lw=1.1, alpha=0.8)
    ax.text(meta_delta, 0.04, f'pooled Δ = {meta_delta:+.2f}',
            transform=ax.get_xaxis_transform(), ha='center', va='bottom',
            fontsize=8.5, color=BLACK_DEEP, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=BLACK_DEEP, lw=0.6, alpha=0.9))
    # Cohort points
    for _, r in sub.iterrows():
        match = np.sign(r.delta)==exp_dir
        col = GOOD_DEEP if match else BAD_DEEP
        ax.scatter(r.delta, r.SE, s=30+3*r.n_total, color=col,
                   edgecolor=BLACK_DEEP, linewidth=0.9, alpha=0.9, zorder=3)
        ax.text(r.delta, r.SE - se_max*0.03, r.gse.replace('GSE',''),
                ha='center', va='top', fontsize=7.5, color=BLACK_DEEP)
    ax.invert_yaxis()
    ax.set_ylim(se_max, 0)
    ax.set_xlim(-1.2, 1.2)
    ax.set_title(f'{title}', fontsize=10.5, fontweight='bold', color=BLACK_DEEP)
    ax.set_ylabel('Standard error', fontsize=9, color=BLACK_DEEP)
    ax.set_xlabel('Δ z-score (good − bad)', fontsize=9, color=BLACK_DEEP)
    add_axis_spines(ax); ax.tick_params(labelsize=8.5)
# Single legend
handles = [mpatches.Patch(color=GOOD_DEEP, label='Cohort agrees with discovery'),
           mpatches.Patch(color=BAD_DEEP, label='Cohort disagrees'),
           Line2D([0],[0], color=GRAY_DEEP, ls='--', label='95% pseudo-CI (null = 0)'),
           Line2D([0],[0], color=BLACK_DEEP, ls=':', label='Inverse-variance pooled Δ')]
fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=8.5,
           frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.04, 1, 1))
save_panel(fig, 'Fig7E_funnel', OUT)

# =========================================================================
# 7F — Discovery ↔ validation concordance summary
#      Slope-style connector between discovery effect (TNT) and pooled validation Δ,
#      flanked by per-cohort dots; right panel shows concordance counts.
# =========================================================================
# Discovery effect (TNT n=35): approximate from expected direction * fixed mag? Use |Δ|=1 sentinel?
# Better: we don't have numeric discovery Δ here — express discovery as direction only.
# Use layout: left col = per-cohort Δ dots (jittered by signature row), center = pooled Δ bar,
#             right col = concordance pill.
fig = plt.figure(figsize=(11, 4.8))
gs = gridspec.GridSpec(1, 3, width_ratios=[4.5, 3.5, 2.2], wspace=0.35)
ax_L = fig.add_subplot(gs[0,0])
ax_M = fig.add_subplot(gs[0,1], sharey=ax_L)
ax_R = fig.add_subplot(gs[0,2], sharey=ax_L)

y_map = {s: i for i, s in enumerate(sig_titles)}

# --- Left: per-cohort Δ strip, colored by agreement ---
for sig, title in zip(sig_names, sig_titles):
    exp_dir = expected[sig]
    sub = ext_stats[ext_stats.signature==sig]
    yy = y_map[title]
    for _, r in sub.iterrows():
        match = np.sign(r.delta)==exp_dir
        col = GOOD_DEEP if match else BAD_DEEP
        ax_L.scatter(r.delta, yy + np.random.uniform(-0.12, 0.12),
                     s=25+2.2*r.n_total, color=col, edgecolor=BLACK_DEEP,
                     linewidth=0.6, alpha=0.75, zorder=3)
ax_L.axvline(0, color=BLACK_DEEP, lw=1.0)
for i in range(len(sig_titles)):
    if i%2==0: ax_L.axhspan(i-0.5, i+0.5, color=BAND, alpha=0.5, zorder=0)
ax_L.set_xlim(-1.2, 1.2)
ax_L.set_xlabel('Per-cohort Δ z-score', fontsize=10, color=BLACK_DEEP)
ax_L.set_yticks(range(len(sig_titles)))
ax_L.set_yticklabels(sig_titles, fontsize=10, color=BLACK_DEEP)
ax_L.invert_yaxis()
add_axis_spines(ax_L); ax_L.tick_params(labelsize=9)
ax_L.set_title('Cohort-level effect sizes', fontsize=10.5, color=BLACK_DEEP,
               fontweight='bold', loc='left')

# --- Middle: pooled Δ bar per signature ---
for sig, title in zip(sig_names, sig_titles):
    exp_dir = expected[sig]
    sub = ext_stats[ext_stats.signature==sig]
    w = 1/sub.SE**2
    pooled = float((sub.delta*w).sum()/w.sum())
    pooled_se = float(np.sqrt(1/w.sum()))
    yy = y_map[title]
    agree = np.sign(pooled)==exp_dir
    col = GOOD_DEEP if agree else BAD_DEEP
    ax_M.barh(yy, pooled, color=col, edgecolor=BLACK_DEEP, linewidth=0.8,
              alpha=0.9, height=0.55)
    ax_M.plot([pooled-1.96*pooled_se, pooled+1.96*pooled_se], [yy, yy],
              color=BLACK_DEEP, lw=1.3)
    tip = pooled + 1.96*pooled_se
    tail = pooled - 1.96*pooled_se
    label_x = (tip + 0.06) if pooled>=0 else (tail - 0.06)
    ha = 'left' if pooled>=0 else 'right'
    ax_M.text(label_x, yy, f'Δ = {pooled:+.2f}', va='center', ha=ha,
              fontsize=9, color=BLACK_DEEP, fontweight='bold')
ax_M.axvline(0, color=BLACK_DEEP, lw=1.0)
ax_M.set_xlim(-1.1, 1.1)
ax_M.set_xlabel('Pooled Δ (inverse-variance)', fontsize=10, color=BLACK_DEEP)
ax_M.tick_params(labelleft=False, labelsize=9)
add_axis_spines(ax_M)
ax_M.set_title('Meta-analytic pooled effect', fontsize=10.5, color=BLACK_DEEP,
               fontweight='bold', loc='left')

# --- Right: concordance pill + counts ---
for sig, title in zip(sig_names, sig_titles):
    exp_dir = expected[sig]
    sub = ext_stats[ext_stats.signature==sig]
    n_match = int((np.sign(sub.delta)==exp_dir).sum())
    n_tot = len(sub)
    meta_row = ext_meta[ext_meta.signature==sig].iloc[0]
    agree_pool = np.sign(meta_row.Z_stouffer)==exp_dir
    yy = y_map[title]
    # Concordance pie
    frac = n_match/n_tot
    ax_R.barh(yy, 1.0, left=0, height=0.45, color='#e9ecef',
              edgecolor=BLACK_DEEP, linewidth=0.6)
    ax_R.barh(yy, frac, left=0, height=0.45,
              color=GOOD_DEEP if agree_pool else BAD_DEEP,
              edgecolor=BLACK_DEEP, linewidth=0.6, alpha=0.9)
    ax_R.text(0.5, yy, f'{n_match}/{n_tot}', ha='center', va='center',
              fontsize=9, color=BLACK_DEEP, fontweight='bold')
    ax_R.text(1.08, yy, '✓' if agree_pool else '✗', va='center', ha='left',
              fontsize=14, color=GOOD_DEEP if agree_pool else BAD_DEEP,
              fontweight='bold')
    ax_R.text(1.45, yy, f'Z = {meta_row.Z_stouffer:+.2f}\np = {meta_row.p_meta_onesided:.2f}',
              va='center', ha='left', fontsize=8.5, color=BLACK_DEEP)
ax_R.set_xlim(0, 2.8)
ax_R.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
for s in ax_R.spines.values(): s.set_visible(False)
ax_R.set_title('Concordance  (cohorts · pooled)', fontsize=10.5,
               color=BLACK_DEEP, fontweight='bold', loc='left')

fig.tight_layout()
save_panel(fig, 'Fig7F_discovery_validation', OUT)

print('\n=== Fig 7 v3 (external validation, journal motifs) complete ===')
