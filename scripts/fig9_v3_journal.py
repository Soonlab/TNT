"""
Figure 9 — Tumor clonal evolution (PyClone-VI, 12 paired pre+post subjects)
Journal-style motifs:
  9A  Per-subject parallel-coordinate clone CCF trajectories (Roth Nat Methods 2014 Fig 3)
  9B  Pre vs Post CCF scatter per clone, faceted by response (Landau Cell 2013 Fig 3; Morrissy Nature 2016)
  9C  Per-subject stacked shrinking/stable/expanding bars (Jamal-Hanjani NEJM 2017 TRACERx Fig 2-3)
  9D  Dominant clone shrinkage violin by response (Rosenthal Nature 2019 Fig 2-3; Anagnostou CD 2017)
  9E  Shrinking vs expanding scatter with marginal histograms (Riaz Cell 2017 Fig 4)
  9F  Aggregated clone-fate composition by response (Morrissy Nature 2016 fate schema)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)
PY = ROOT/'04_wes_cnv_clonal/pyclone'

GOOD_DEEP='#0a7d6e'; BAD_DEEP='#c53e1f'; BLACK_DEEP='#0e2a47'; GOLD_DEEP='#d4a300'
GRAY_DEEP='#5a6772'; BAND='#ecedef'
PAL_RESP_DEEP = {'good': GOOD_DEEP, 'bad': BAD_DEEP}
# Fate colors
SHRINK = '#118ab2'; STABLE = '#8d99ae'; EXPAND = '#d62828'
FATE_COL = {'shrink': SHRINK, 'stable': STABLE, 'expand': EXPAND}

pyclone = pd.read_csv(PY/'clonal_summary.tsv', sep='\t')
pyclone = pyclone.sort_values(['response','subject_id']).reset_index(drop=True)

# Load per-subject results; compute per-cluster mean CCF per timepoint
def load_clusters(subj):
    f = PY/f'results_subj{int(subj)}.tsv'
    if not f.exists(): return None
    r = pd.read_csv(f, sep='\t')
    # Timepoint from sample_id suffix -PR (pre) / -PO (post)
    r['tp'] = r.sample_id.str.extract(r'-(PR|PO)$')[0].map({'PR':'pre','PO':'post'})
    agg = r.groupby(['cluster_id','tp']).cellular_prevalence.mean().unstack()
    if 'pre' not in agg.columns or 'post' not in agg.columns:
        return None
    # size = n mutations per cluster
    agg['n'] = r.groupby('cluster_id').mutation_id.nunique()
    agg['delta'] = agg.post - agg.pre
    agg['fate'] = np.where(agg.delta < -0.2, 'shrink',
                  np.where(agg.delta > +0.2, 'expand', 'stable'))
    return agg.reset_index()

subj_clusters = {}
for s in pyclone.subject_id:
    c = load_clusters(s)
    if c is not None: subj_clusters[int(s)] = c

# ============================================================
# 9A — Per-subject fishplot-style clonal trajectories
#      (Miller fishplot / Morrissy Nature 2016 fishplot motif)
#      Each clone = smooth ribbon whose height = cellular prevalence,
#      stacked symmetrically around midline (largest clone = outer envelope).
# ============================================================
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

def smooth_ribbon(ax, x0, x1, y0_top, y0_bot, y1_top, y1_bot, color, alpha=0.9):
    """Draw a smooth Bezier ribbon between two vertical segments."""
    xm = (x0+x1)/2
    verts = [
        (x0, y0_top),
        (xm, y0_top), (xm, y1_top), (x1, y1_top),  # top edge cubic
        (x1, y1_bot),
        (xm, y1_bot), (xm, y0_bot), (x0, y0_bot),  # bottom edge cubic
        (x0, y0_top),
    ]
    codes = [MplPath.MOVETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    p = PathPatch(MplPath(verts, codes), facecolor=color, edgecolor='white',
                  lw=0.6, alpha=alpha, zorder=3)
    ax.add_patch(p)

n = len(pyclone); ncol=4; nrow=int(np.ceil(n/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(14, 7.8))
# Cluster palette (ordered distinct colors, fall back to fate-based)
CLUSTER_PAL = ['#1d3557','#118ab2','#06aed5','#0a7d6e','#d4a300','#e07a5f',
               '#c53e1f','#7b2cbf','#6c757d','#8d99ae']

for ax, (_, row) in zip(axes.flat, pyclone.iterrows()):
    subj = int(row.subject_id)
    c = subj_clusters.get(subj)
    if c is None:
        ax.set_visible(False); continue
    resp_c = PAL_RESP_DEEP[row.response]
    # Sort clones largest → smallest by max(pre,post) so big ones form outer shape
    c = c.copy()
    c['size_max'] = c[['pre','post']].max(axis=1)
    c = c.sort_values('size_max', ascending=False).reset_index(drop=True)
    # Stack symmetrically: each clone occupies a band centered on midline,
    # smaller clones offset inward. Use cumulative half-heights.
    # Allocate y-positions: largest gets full envelope; each smaller clone
    # nests inside the remaining space at its own pre/post heights.
    # Simpler scheme: stack from bottom, then reflect to create fish silhouette.
    pre_vals  = c.pre.values.astype(float)
    post_vals = c.post.values.astype(float)
    # cumulative offsets (half-stack)
    pre_cum  = np.concatenate([[0], np.cumsum(pre_vals)])
    post_cum = np.concatenate([[0], np.cumsum(post_vals)])
    pre_total, post_total = pre_cum[-1], post_cum[-1]
    # center the stack vertically around 0.5 envelope based on max total
    for k in range(len(c)):
        col = CLUSTER_PAL[k % len(CLUSTER_PAL)]
        # Symmetric stacking: top/bottom pair
        y0_top = 0.5 + (pre_cum[k+1]  - pre_total/2)
        y0_bot = 0.5 + (pre_cum[k]    - pre_total/2)
        y1_top = 0.5 + (post_cum[k+1] - post_total/2)
        y1_bot = 0.5 + (post_cum[k]   - post_total/2)
        smooth_ribbon(ax, 0, 1, y0_top, y0_bot, y1_top, y1_bot, col, alpha=0.92)
        # Fate marker dot at midpoint of ribbon at post end
        fate = c.iloc[k].fate
        midy_post = (y1_top + y1_bot)/2
        ax.scatter(1.03, midy_post, s=32, color=FATE_COL[fate],
                   edgecolor=BLACK_DEEP, linewidth=0.6, zorder=5, clip_on=False)
    # Timepoint baseline
    ax.axvline(0, color=BLACK_DEEP, lw=0.8); ax.axvline(1, color=BLACK_DEEP, lw=0.8)
    ax.set_xticks([0,1]); ax.set_xticklabels(['pre','post'], fontsize=9, color=BLACK_DEEP)
    ax.set_xlim(-0.05, 1.12); ax.set_ylim(-0.02, 1.05)
    ax.set_yticks([])
    # Header
    ax.text(0.5, 1.03, f'S{subj}  ·  {row.response}  ·  k={int(row.n_clusters)}',
            transform=ax.transAxes, ha='center', va='bottom', fontsize=9.5,
            color=resp_c, fontweight='bold')
    for s_side in ['top','right','left']: ax.spines[s_side].set_visible(False)
    ax.spines['bottom'].set_color(BLACK_DEEP)
    ax.tick_params(left=False, labelsize=8.5)
for k in range(len(pyclone), nrow*ncol):
    axes.flat[k].set_visible(False)
# Legend
leg = [mpatches.Patch(color=CLUSTER_PAL[0], label='Clone ribbon (height ∝ cellular prevalence)'),
       Line2D([0],[0], marker='o', color='w', markerfacecolor=SHRINK,
              markeredgecolor=BLACK_DEEP, markersize=9, label='Shrinking fate'),
       Line2D([0],[0], marker='o', color='w', markerfacecolor=STABLE,
              markeredgecolor=BLACK_DEEP, markersize=9, label='Stable fate'),
       Line2D([0],[0], marker='o', color='w', markerfacecolor=EXPAND,
              markeredgecolor=BLACK_DEEP, markersize=9, label='Expanding fate')]
fig.legend(handles=leg, loc='lower center', ncol=4, fontsize=9.5, frameon=False,
           bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.04, 1, 1))
save_panel(fig, 'Fig9A_clone_trajectories', OUT)

# ============================================================
# 9B — Pre vs Post CCF scatter per clone, faceted by response (Landau CLL Fig 3)
# ============================================================
all_clones = []
for subj, c in subj_clusters.items():
    resp = pyclone.loc[pyclone.subject_id==subj, 'response'].iloc[0]
    cc = c.copy(); cc['subj']=subj; cc['response']=resp
    all_clones.append(cc)
all_clones = pd.concat(all_clones, ignore_index=True)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)
for ax, resp in zip(axes, ['good','bad']):
    sub = all_clones[all_clones.response==resp]
    # Shaded Δ>0.2 bands (expansion zone above diagonal, shrink below)
    xs = np.linspace(0,1,100)
    ax.fill_between(xs, xs+0.2, 1.05, color=EXPAND, alpha=0.07, zorder=0)
    ax.fill_between(xs, -0.05, xs-0.2, color=SHRINK, alpha=0.07, zorder=0)
    ax.plot([0,1],[0,1], color=BLACK_DEEP, ls='--', lw=1, alpha=0.8, zorder=1)
    ax.plot(xs, xs+0.2, color=EXPAND, ls=':', lw=0.8, alpha=0.55, zorder=1)
    ax.plot(xs, xs-0.2, color=SHRINK, ls=':', lw=0.8, alpha=0.55, zorder=1)
    for _, r in sub.iterrows():
        ax.scatter(r.pre, r.post, s=40+2.2*min(r.n,80), color=FATE_COL[r.fate],
                   edgecolor=BLACK_DEEP, linewidth=0.6, alpha=0.85, zorder=3)
    # Fate count annotation
    counts = sub.fate.value_counts().to_dict()
    txt = (f'shrink: {counts.get("shrink",0)}   '
           f'stable: {counts.get("stable",0)}   '
           f'expand: {counts.get("expand",0)}')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, ha='left', va='top',
            fontsize=9, color=BLACK_DEEP, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                      ec=PAL_RESP_DEEP[resp], lw=1.0, alpha=0.95))
    ax.text(0.98, 0.02,
            f'{resp}  ·  n={sub.subj.nunique()} subjects  ·  {len(sub)} clones',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9.5,
            color=PAL_RESP_DEEP[resp], fontweight='bold')
    ax.set_xlabel('Pre-TNT cellular prevalence', fontsize=10, color=BLACK_DEEP)
    ax.set_ylabel('Post-TNT cellular prevalence', fontsize=10, color=BLACK_DEEP)
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
    add_axis_spines(ax); ax.tick_params(labelsize=9)
fig.tight_layout()
save_panel(fig, 'Fig9B_CCF_pre_post', OUT)

# ============================================================
# 9C — Per-subject stacked shrink/stable/expand bars (TRACERx motif)
# ============================================================
py_sorted = pyclone.copy()
py_sorted['n_stable'] = py_sorted.n_clusters - py_sorted.n_shrinking - py_sorted.n_expanding
py_sorted = py_sorted.sort_values(['response','n_shrinking','n_expanding'],
                                   ascending=[True, False, True]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(10, 4.5))
y = np.arange(len(py_sorted))
left = np.zeros(len(py_sorted))
for cat, col, label in [('n_shrinking', SHRINK, 'Shrinking'),
                         ('n_stable',    STABLE, 'Stable'),
                         ('n_expanding', EXPAND, 'Expanding')]:
    vals = py_sorted[cat].values
    ax.barh(y, vals, left=left, height=0.68, color=col,
            edgecolor='white', linewidth=1.2, label=label, alpha=0.92)
    for i,v in enumerate(vals):
        if v>0:
            ax.text(left[i]+v/2, i, str(int(v)), ha='center', va='center',
                    fontsize=8.5, color='white', fontweight='bold')
    left += vals
# Response colored y-tick labels
ax.set_yticks(y)
ax.set_yticklabels([f'S{int(r.subject_id)}' for _, r in py_sorted.iterrows()],
                   fontsize=9.5)
for tick, resp in zip(ax.get_yticklabels(), py_sorted.response):
    tick.set_color(PAL_RESP_DEEP[resp]); tick.set_fontweight('bold')
# Response sidebar
for i, r in py_sorted.iterrows():
    ax.add_patch(Rectangle((-0.8, i-0.34), 0.5, 0.68,
                           facecolor=PAL_RESP_DEEP[r.response],
                           edgecolor='white', lw=0.7, clip_on=False))
ax.set_xlabel('# PyClone clusters', fontsize=10, color=BLACK_DEEP)
ax.set_xlim(0, py_sorted.n_clusters.max()+0.5)
ax.invert_yaxis()
ax.legend(loc='lower right', fontsize=9, frameon=False, ncol=3,
          bbox_to_anchor=(0.98, -0.28))
add_axis_spines(ax); ax.tick_params(labelsize=9)
fig.tight_layout()
save_panel(fig, 'Fig9C_cluster_stacked', OUT)

# ============================================================
# 9D — Dominant clone shrinkage by response (Rosenthal / Anagnostou)
# ============================================================
fig, ax = plt.subplots(figsize=(4.8, 4.3))
for i, g in enumerate(['good','bad']):
    vals = pyclone[pyclone.response==g].dominant_shrink.values
    col = PAL_RESP_DEEP[g]
    if len(vals) >= 2 and np.std(vals) > 1e-6:
        try:
            kde = stats.gaussian_kde(vals, bw_method=0.55)
            xg = np.linspace(vals.min()-0.1, vals.max()+0.1, 120)
            d = kde(xg); d = d/d.max()*0.35
            ax.fill_betweenx(xg, i, i-d, color=col, alpha=0.32, linewidth=0)
            ax.plot(i-d, xg, color=col, lw=1.1, alpha=0.85)
        except Exception: pass
    q1, med, q3 = np.percentile(vals, [25,50,75])
    ax.add_patch(Rectangle((i+0.05, q1), 0.14, q3-q1, facecolor='white',
                           edgecolor=col, lw=1.3))
    ax.plot([i+0.05, i+0.19], [med, med], color=col, lw=2.2)
    jit = np.random.uniform(0.26, 0.44, len(vals))
    ax.scatter(i+jit, vals, s=55, color=col, edgecolor=BLACK_DEEP,
               linewidth=0.6, alpha=0.88, zorder=3)
g = pyclone[pyclone.response=='good'].dominant_shrink.values
b = pyclone[pyclone.response=='bad'].dominant_shrink.values
p = stats.mannwhitneyu(g, b).pvalue
ymax = max(g.max(), b.max()); ymin = min(g.min(), b.min())
y_br = ymax + 0.08
ax.plot([0,0,1,1],[y_br, y_br+0.03, y_br+0.03, y_br], color=BLACK_DEEP, lw=1)
ax.text(0.5, y_br+0.04, f'Mann–Whitney p = {p:.2f}', ha='center', fontsize=9.5,
        color=BLACK_DEEP, fontweight='bold')
ax.axhline(0, color=GRAY_DEEP, ls='--', lw=0.7, alpha=0.7)
ax.text(1.52, 0, '← no change', ha='left', va='center', fontsize=8, color=GRAY_DEEP)
ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'], fontsize=11, color=BLACK_DEEP)
ax.tick_params(axis='x', pad=10)
ax.set_ylabel('Min Δ CP  (post − pre)\n(more negative = more clearance)',
              fontsize=9.5, color=BLACK_DEEP)
ax.set_xlim(-0.55, 1.7)
add_axis_spines(ax)
save_panel(fig, 'Fig9D_dominant_shrink', OUT)

# ============================================================
# 9E — Shrinking vs expanding scatter with marginal counts (Riaz Cell 2017 Fig 4)
# ============================================================
fig = plt.figure(figsize=(6.2, 5.8))
gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                       wspace=0.05, hspace=0.05)
ax = fig.add_subplot(gs[1,0])
ax_tx = fig.add_subplot(gs[0,0], sharex=ax)
ax_ry = fig.add_subplot(gs[1,1], sharey=ax)

for resp in ['good','bad']:
    sub = pyclone[pyclone.response==resp]
    c = PAL_RESP_DEEP[resp]
    # jitter slightly
    jx = sub.n_shrinking + np.random.uniform(-0.08, 0.08, len(sub))
    jy = sub.n_expanding + np.random.uniform(-0.08, 0.08, len(sub))
    ax.scatter(jx, jy, s=200, color=c, edgecolor=BLACK_DEEP, linewidth=1.0,
               alpha=0.88, zorder=3, label=f'{resp} (n={len(sub)})')
    for xi, yi, s in zip(jx, jy, sub.subject_id):
        ax.text(xi+0.08, yi+0.08, f'S{int(s)}', fontsize=8, color=BLACK_DEEP,
                fontweight='bold')
ax.plot([-0.5,4],[-0.5,4], color=GRAY_DEEP, ls='--', lw=0.7, alpha=0.7)
ax.set_xlabel('# Shrinking clusters  (Δ CP < −0.2)', fontsize=10, color=BLACK_DEEP)
ax.set_ylabel('# Expanding clusters  (Δ CP > +0.2)', fontsize=10, color=BLACK_DEEP)
ax.set_xticks([0,1,2,3]); ax.set_yticks([0,1,2])
ax.set_xlim(-0.4, 3.4); ax.set_ylim(-0.4, 2.4)
add_axis_spines(ax); ax.tick_params(labelsize=9)
ax.legend(loc='upper right', fontsize=9, frameon=False)

# Top marginal: shrinking histogram by response
bins = np.arange(-0.25, 3.75, 0.5)
for resp in ['good','bad']:
    vals = pyclone[pyclone.response==resp].n_shrinking.values
    ax_tx.hist(vals, bins=bins, color=PAL_RESP_DEEP[resp], alpha=0.55,
               edgecolor=BLACK_DEEP, linewidth=0.6)
ax_tx.tick_params(bottom=False, labelbottom=False, labelleft=False, left=False,
                  length=0)
for s in ax_tx.spines.values(): s.set_visible(False)

# Right marginal: expanding histogram
binsY = np.arange(-0.25, 2.75, 0.5)
for resp in ['good','bad']:
    vals = pyclone[pyclone.response==resp].n_expanding.values
    ax_ry.hist(vals, bins=binsY, color=PAL_RESP_DEEP[resp], alpha=0.55,
               edgecolor=BLACK_DEEP, linewidth=0.6, orientation='horizontal')
ax_ry.tick_params(left=False, labelleft=False, labelbottom=False, bottom=False,
                  length=0)
for s in ax_ry.spines.values(): s.set_visible(False)

save_panel(fig, 'Fig9E_shrink_expand_scatter', OUT)

# ============================================================
# 9F — Aggregated clone-fate composition by response (Morrissy fate schema)
# ============================================================
fig, ax = plt.subplots(figsize=(6.5, 4.2))
agg = all_clones.groupby(['response','fate']).size().unstack(fill_value=0)
for c in ['shrink','stable','expand']:
    if c not in agg.columns: agg[c] = 0
agg = agg[['shrink','stable','expand']]
totals = agg.sum(axis=1)
pct = agg.div(totals, axis=0) * 100

order = ['good','bad']
y = np.arange(len(order))
left = np.zeros(len(order))
for cat, col, label in [('shrink', SHRINK, 'Shrinking'),
                         ('stable', STABLE, 'Stable'),
                         ('expand', EXPAND, 'Expanding')]:
    vals = pct.loc[order, cat].values
    nvals = agg.loc[order, cat].values
    bars = ax.barh(y, vals, left=left, height=0.55, color=col,
                   edgecolor='white', linewidth=1.5, label=label, alpha=0.92)
    for i,(v,nv) in enumerate(zip(vals, nvals)):
        if v > 4:
            ax.text(left[i]+v/2, i, f'{v:.0f}%\n(n={int(nv)})', ha='center',
                    va='center', fontsize=9, color='white', fontweight='bold')
    left += vals

# Chi-square or Fisher
from scipy.stats import chi2_contingency
try:
    chi2, chi_p, *_ = chi2_contingency(agg.loc[order].values)
except Exception:
    chi_p = np.nan
ax.text(1.02, 0.5, f'χ² p = {chi_p:.2f}\n(fate × response)', transform=ax.transAxes,
        ha='left', va='center', fontsize=9.5, color=BLACK_DEEP, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=BLACK_DEEP, lw=0.7, alpha=0.95))
ax.set_yticks(y); ax.set_yticklabels(
    [f'{r}  (N={int(totals[r])} clones,  {pyclone[pyclone.response==r].shape[0]} subj)' for r in order],
    fontsize=10, color=BLACK_DEEP)
for tick, r in zip(ax.get_yticklabels(), order):
    tick.set_color(PAL_RESP_DEEP[r]); tick.set_fontweight('bold')
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xlabel('% of clones', fontsize=10, color=BLACK_DEEP)
ax.legend(loc='lower center', ncol=3, fontsize=9, frameon=False,
          bbox_to_anchor=(0.5, -0.32))
add_axis_spines(ax); ax.tick_params(labelsize=9)
fig.tight_layout()
save_panel(fig, 'Fig9F_fate_composition', OUT)

print('\n=== Fig 9 v3 (clonal evolution, journal motifs) complete ===')
