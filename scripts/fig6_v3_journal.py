"""
Figure 6 v3 — journal-style treatment-induced response cascade
Motifs:
  6A: Krishna Cancer Cell 2024 + Bi Cancer Cell 2021 — polished swimmer plot with multi-event tracks
  6B: Bassez Nat Med 2021 + Helmink Nature 2020 — paired slope graphs with directional arrows per subject
  6C: Mariathasan Nature 2018 + Litchfield Cell 2021 — Δ effect-size forest plot
  6D: Liu Nat Med 2021 + Krishna Cancer Cell 2024 — per-subject multi-feature delta waterfall
  6E: Oliveira Nature 2022 + Lin Nature 2024 — clonal stream plot (PyClone cluster trajectories)
  6F: illustrative cascade diagram (Mariathasan/Helmink style)

Style: no titles, saturated colors, no overlapping text.
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from scipy.interpolate import interp1d
import warnings; warnings.filterwarnings('ignore')

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'; GOLD = '#d4a300'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
trust = pd.read_csv(ROOT/'06_rna_immune/trust4_summary.tsv', sep='\t')
neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
neo_delta = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')
delta_22 = pd.read_csv(ROOT/'09_integration/paired_delta/delta_22sigs_response.tsv', sep='\t')
delta_trust = pd.read_csv(ROOT/'09_integration/paired_delta/delta_trust4_response.tsv', sep='\t')
sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_summary_key.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
pyclone = pd.read_csv(ROOT/'04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\t')

PAIRED_SUBJ = list(range(1, 15))

# Build subject-level Δ master
def build_delta():
    rows = []
    for s in PAIRED_SUBJ:
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        # WES deltas
        ts = tmb[(tmb.subject_id==s) & tmb.timepoint.isin(['pre','post'])]
        d_tmb = d_miss = np.nan
        if {'pre','post'}.issubset(set(ts.timepoint)):
            pre_t = ts[ts.timepoint=='pre'].iloc[0]; post_t = ts[ts.timepoint=='post'].iloc[0]
            d_tmb = post_t.TMB_nonsyn_per_Mb - pre_t.TMB_nonsyn_per_Mb
            d_miss = post_t.n_nonsyn - pre_t.n_nonsyn
        # SBS5
        ss = sbs[(sbs.subject_id==s) & sbs.timepoint.isin(['pre','post'])]
        d_sbs = np.nan
        if {'pre','post'}.issubset(set(ss.timepoint)):
            d_sbs = ss[ss.timepoint=='post'].iloc[0].SBS5 - ss[ss.timepoint=='pre'].iloc[0].SBS5
        # Neo
        nd = neo_delta[neo_delta.subject_id==s]
        d_neo = nd.delta_binders.iloc[0] if len(nd) else np.nan
        # HLA LOH
        loh_pre = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PR')]
        loh_post = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PO')]
        d_loh = (loh_post.LOH_call.sum() - loh_pre.LOH_call.sum()) if len(loh_pre)+len(loh_post)>0 else np.nan
        # RNA
        rs = rna_inv[(rna_inv.subject_id==s) & rna_inv.timepoint.isin(['pre','post'])]
        d_treg = d_mhc = d_cd8exh = d_igh = np.nan
        if {'pre','post'}.issubset(set(rs.timepoint)):
            pre_id = rs[rs.timepoint=='pre'].sample_id.iloc[0]
            post_id = rs[rs.timepoint=='post'].sample_id.iloc[0]
            try:
                d_treg = sigs.loc[post_id,'Treg'] - sigs.loc[pre_id,'Treg']
                d_mhc = sigs.loc[post_id,'MHC_II'] - sigs.loc[pre_id,'MHC_II']
                d_cd8exh = sigs.loc[post_id,'CD8_exhaustion'] - sigs.loc[pre_id,'CD8_exhaustion']
            except: pass
            try:
                pre_tr = trust[trust.sample_id==pre_id]; post_tr = trust[trust.sample_id==post_id]
                d_igh = post_tr.IGH_n.iloc[0] - pre_tr.IGH_n.iloc[0]
            except: pass
        rows.append({'subject_id':s, 'response':resp,
                     'd_TMB':d_tmb, 'd_missense':d_miss, 'd_SBS5':d_sbs, 'd_neo':d_neo, 'd_HLA_LOH':d_loh,
                     'd_Treg':d_treg, 'd_MHC_II':d_mhc, 'd_CD8_exh':d_cd8exh, 'd_IGH':d_igh})
    return pd.DataFrame(rows)

D = build_delta()

# ============================================================
# 6A — Swimmer plot (Krishna Cancer Cell 2024 inspired)
# ============================================================
D_sw = D.sort_values(['response','subject_id']).reset_index(drop=True)
N_sub = len(D_sw)
fig, ax = plt.subplots(figsize=(13.5, 7))

metrics = [
    ('d_TMB',      'TMB',          'clearance', 0.05),
    ('d_missense', 'Missense',     'clearance', 0.16),
    ('d_SBS5',     'SBS5 muts',    'clearance', 0.27),
    ('d_neo',      'Neoantigens',  'clearance', 0.40),
    ('d_HLA_LOH',  'HLA LOH',      'clearance', 0.53),
    ('d_Treg',     'Treg',         'activation',0.64),
    ('d_MHC_II',   'MHC-II',       'activation',0.73),
    ('d_CD8_exh',  'CD8 exh.',     'activation',0.82),
    ('d_IGH',      'BCR (IGH)',    'activation',0.94),
]

# Row bg strip: response color, light alpha
for i, (_, r) in enumerate(D_sw.iterrows()):
    ax.axhspan(i-0.4, i+0.4, facecolor=PAL[r.response], alpha=0.10, zorder=0)

# Events
for i, (_, r) in enumerate(D_sw.iterrows()):
    for m, lbl, direction, x_pos in metrics:
        v = r[m]
        if pd.isna(v): continue
        if direction=='clearance':
            # negative = good (blue), positive = bad (red)
            color = '#118ab2' if v<0 else '#c11456'
            mag = abs(v)
        else:
            # positive = good (green), negative = bad (orange)
            color = '#0f8b78' if v>0 else '#d96125'
            mag = abs(v)
        # Size scaling per metric family
        if m == 'd_neo': s = min(500, 50 + mag*0.8)
        elif m == 'd_IGH': s = min(500, 50 + mag*0.15)
        elif m == 'd_missense': s = min(500, 50 + mag*3)
        elif m == 'd_SBS5': s = min(500, 50 + mag*3)
        elif m == 'd_HLA_LOH': s = 80 + mag*80
        elif m == 'd_TMB': s = min(500, 50 + mag*40)
        else: s = min(500, 40 + mag*80)
        ax.scatter(x_pos, i, s=s, c=color, alpha=0.88, edgecolor='white', lw=1.2, zorder=5)

# Headers for metrics (top of plot)
for m, lbl, direction, x_pos in metrics:
    ax.text(x_pos, N_sub-0.2, lbl, ha='center', va='bottom', fontsize=9.5,
            color='#0e2a47', fontweight='bold', rotation=35, rotation_mode='anchor')

# Top category separator lines
ax.axvline(0.58, color='#5a6772', lw=0.6, ls=':', alpha=0.5, ymin=0.02, ymax=0.90)
ax.text(0.295, N_sub+0.35, 'Tumor clearance (post − pre < 0 = ↓tumor)', ha='center', fontsize=9.5,
        color='#0e2a47', fontweight='bold', style='italic')
ax.text(0.78, N_sub+0.35, 'Immune activation (post − pre > 0 = ↑immune)', ha='center', fontsize=9.5,
        color='#0e2a47', fontweight='bold', style='italic')

# Y-axis: subject + response
ax.set_yticks(range(N_sub))
ax.set_yticklabels([f"S{int(r.subject_id)}" for _, r in D_sw.iterrows()], fontsize=9.5, color='#0e2a47')
# Response column beside y-labels (response color dot)
for i, (_, r) in enumerate(D_sw.iterrows()):
    ax.scatter(-0.03, i, s=85, c=PAL[r.response], edgecolor='#0e2a47', lw=1.0, clip_on=False, zorder=6)

# Good/poor bracket
good_n = (D_sw.response=='good').sum()
ax.axhline(good_n-0.5, color='#0e2a47', lw=1.2)
ax.text(-0.08, (good_n-1)/2, 'Good', ha='right', va='center', fontsize=11.5,
        color=GOOD_DEEP, fontweight='bold', rotation=90)
ax.text(-0.08, good_n + (N_sub-good_n-1)/2, 'Poor', ha='right', va='center', fontsize=11.5,
        color=BAD_DEEP, fontweight='bold', rotation=90)

ax.set_ylim(-0.7, N_sub+1.2); ax.set_xlim(-0.02, 1.0)
ax.set_xticks([])
for s in ['top','right','bottom','left']: ax.spines[s].set_visible(False)

# Legend (bottom)
leg_elems = [
    plt.scatter([],[], s=160, c='#118ab2', edgecolor='white', lw=1.2, label='↓ Tumor (Δ < 0, favorable)'),
    plt.scatter([],[], s=160, c='#c11456', edgecolor='white', lw=1.2, label='↑ Tumor (Δ > 0, unfavorable)'),
    plt.scatter([],[], s=160, c='#0f8b78', edgecolor='white', lw=1.2, label='↑ Immune (Δ > 0, favorable)'),
    plt.scatter([],[], s=160, c='#d96125', edgecolor='white', lw=1.2, label='↓ Immune (Δ < 0, unfavorable)'),
]
ax.legend(handles=leg_elems, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=4,
          fontsize=9.5, frameon=False)

save_panel(fig, 'Fig6A_swimmer', OUT)

# ============================================================
# 6B — Paired slope graph with directional arrows (Bassez / Helmink)
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharex=False)

def slope_arrow_panel(ax, pre_map, post_map, ylabel, ylabel_units=''):
    xs_pre = 0.0; xs_post = 1.0
    for s in PAIRED_SUBJ:
        if s not in pre_map or s not in post_map: continue
        vp = pre_map[s]; vpo = post_map[s]
        if pd.isna(vp) or pd.isna(vpo): continue
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        color = PAL[resp]
        # Arrow
        arr = FancyArrowPatch((xs_pre, vp), (xs_post, vpo),
                               arrowstyle='-|>', color=color, lw=1.8, alpha=0.8,
                               mutation_scale=13)
        ax.add_patch(arr)
        # Dots
        ax.scatter(xs_pre, vp, s=60, color=color, edgecolor='white', lw=0.8, zorder=3)
        ax.scatter(xs_post, vpo, s=60, color=color, edgecolor='white', lw=0.8, zorder=3)
        # S label (only for extreme)
    ax.set_xticks([xs_pre, xs_post])
    ax.set_xticklabels(['pre', 'post'], fontsize=11, color='#0e2a47', fontweight='bold')
    ax.set_xlim(-0.20, 1.20)
    ax.set_ylabel(ylabel + (f' ({ylabel_units})' if ylabel_units else ''),
                  fontsize=10.5, fontweight='bold', color='#0e2a47')
    add_axis_spines(ax)

# TMB
tmb_pre = {s: tmb[tmb.sample_id==f'{s}-PR'].TMB_nonsyn_per_Mb.iloc[0]
           if (tmb.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ}
tmb_post = {s: tmb[tmb.sample_id==f'{s}-PO'].TMB_nonsyn_per_Mb.iloc[0]
            if (tmb.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ}
slope_arrow_panel(axes[0], tmb_pre, tmb_post, 'TMB', '/Mb')

# Neoantigen binders
neo_pre = {s: neo[neo.sample_id==f'{s}-PR'].n_binders_500nM.iloc[0]
           if (neo.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ}
neo_post = {s: neo[neo.sample_id==f'{s}-PO'].n_binders_500nM.iloc[0]
            if (neo.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ}
slope_arrow_panel(axes[1], neo_pre, neo_post, 'MHC-I binders', '<500 nM')

# IGH clonotype count
rna_pre_ids = {s: rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='pre')].sample_id.iloc[0]
               if len(rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='pre')]) else None for s in PAIRED_SUBJ}
rna_post_ids = {s: rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='post')].sample_id.iloc[0]
                if len(rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='post')]) else None for s in PAIRED_SUBJ}
def trust_get(sid, col):
    if sid is None: return np.nan
    sub = trust[trust.sample_id==sid]
    return sub[col].iloc[0] if len(sub) else np.nan
igh_pre = {s: trust_get(rna_pre_ids[s], 'IGH_n') for s in PAIRED_SUBJ}
igh_post = {s: trust_get(rna_post_ids[s], 'IGH_n') for s in PAIRED_SUBJ}
slope_arrow_panel(axes[2], igh_pre, igh_post, 'BCR (IGH) clonotypes', 'count')

# Treg
def sig_get(sid, col):
    if sid is None or sid not in sigs.index: return np.nan
    return sigs.loc[sid, col]
treg_pre = {s: sig_get(rna_pre_ids[s], 'Treg') for s in PAIRED_SUBJ}
treg_post = {s: sig_get(rna_post_ids[s], 'Treg') for s in PAIRED_SUBJ}
slope_arrow_panel(axes[3], treg_pre, treg_post, 'Treg signature', 'z-score')

# Legend
fig.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good responders'),
                    mpatches.Patch(color=BAD_DEEP, label='Poor responders')],
           loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=2, fontsize=11, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig6B_slope_arrows', OUT)

# ============================================================
# 6C — Δ effect-size forest plot (Mariathasan / Litchfield)
# ============================================================
fig, ax = plt.subplots(figsize=(9, 6.5))
# Top features by paired-delta significance
# Use existing tables
feats_order = [
    ('d_missense',  'Missense mutations Δ',    'tumor'),
    ('d_TMB',       'TMB Δ',                   'tumor'),
    ('d_SBS5',      'SBS5 mutations Δ',        'tumor'),
    ('d_neo',       'MHC-I neoantigens Δ',     'tumor'),
    ('d_HLA_LOH',   'HLA LOH events Δ',        'tumor'),
    ('d_Treg',      'Treg signature Δ',        'immune'),
    ('d_MHC_II',    'MHC-II signature Δ',      'immune'),
    ('d_CD8_exh',   'CD8 exhaustion Δ',        'immune'),
    ('d_IGH',       'IGH (BCR) clonotypes Δ',  'immune'),
]

rows = []
np.random.seed(42)
def boot_median_diff(g_vals, b_vals, n=2000):
    diffs = []
    for _ in range(n):
        gi = np.random.choice(g_vals, len(g_vals), replace=True)
        bi = np.random.choice(b_vals, len(b_vals), replace=True)
        diffs.append(np.median(gi) - np.median(bi))
    return np.median(diffs), np.percentile(diffs, [2.5, 97.5])

for fkey, label, cat in feats_order:
    g_vals = D[D.response=='good'][fkey].dropna().values
    b_vals = D[D.response=='bad'][fkey].dropna().values
    if len(g_vals)<3 or len(b_vals)<3: continue
    from scipy.stats import mannwhitneyu
    p = mannwhitneyu(g_vals, b_vals).pvalue
    diff, (ci_lo, ci_hi) = boot_median_diff(g_vals, b_vals)
    rows.append({'feature':label, 'cat':cat, 'diff':diff, 'ci_lo':ci_lo, 'ci_hi':ci_hi, 'p':p,
                 'g_med':np.median(g_vals), 'b_med':np.median(b_vals)})
forest = pd.DataFrame(rows)

# Normalize effect sizes by feature std for visual comparability
for i, r in forest.iterrows():
    fkey = [f[0] for f in feats_order if f[1]==r.feature][0]
    all_vals = D[fkey].dropna().values
    std = all_vals.std() if all_vals.std()>0 else 1
    forest.loc[i, 'diff_norm'] = r['diff'] / std
    forest.loc[i, 'ci_lo_norm'] = r['ci_lo'] / std
    forest.loc[i, 'ci_hi_norm'] = r['ci_hi'] / std

forest = forest.iloc[::-1].reset_index(drop=True)
y_pos = np.arange(len(forest))
for i, r in forest.iterrows():
    color = GOOD_DEEP if r.diff_norm>0 else BAD_DEEP
    # If category tumor-clearance and diff<0 → good responder clears more → green
    # Actually: for tumor features, MORE NEGATIVE Δ good is favorable → if d_miss good − d_miss bad is NEGATIVE means good cleared more
    # For immune features, MORE POSITIVE Δ good is favorable
    if r['cat'] == 'tumor':
        color = GOOD_DEEP if r.diff_norm<0 else BAD_DEEP
    else:
        color = GOOD_DEEP if r.diff_norm>0 else BAD_DEEP
    ax.plot([r.ci_lo_norm, r.ci_hi_norm], [i, i], color=color, lw=2.4, alpha=0.85, solid_capstyle='round')
    ax.scatter(r.diff_norm, i, s=180 if r.p<0.1 else 100, color=color, edgecolor='white',
               linewidth=1.4, zorder=3)
    star = sig_symbol(r.p); star = star if star != 'ns' else ''
    ax.text(max(r.ci_hi_norm, r.diff_norm) + 0.1, i, f'p = {r.p:.3g} {star}',
            va='center', ha='left', fontsize=9, color='#0e2a47')

ax.axvline(0, color='#0e2a47', lw=1.0)
ax.axvspan(-0.1, 0.1, color='#dee2e6', alpha=0.4)
ax.set_yticks(y_pos)
ax.set_yticklabels(forest.feature, fontsize=10, color='#0e2a47')
ax.set_xlabel('Normalized Δ (good − poor)   95% bootstrap CI',
              fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(forest.ci_lo_norm.min()*1.1, forest.ci_hi_norm.max()*1.8)
add_axis_spines(ax)

# Direction guide
ax.text(0.98, 0.97, 'Green = favorable in good\nCoral = unfavorable',
        transform=ax.transAxes, ha='right', va='top', fontsize=9, color='#0e2a47',
        fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.92, boxstyle='round,pad=0.4'))

save_panel(fig, 'Fig6C_delta_forest', OUT)

# ============================================================
# 6D — Per-subject multi-feature delta waterfall (Liu Nat Med 2021 / Krishna Cancer Cell 2024)
# ============================================================
fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)

# Sort subjects by response then missense delta
D_ord = D.sort_values(['response','d_missense']).reset_index(drop=True)

for ax_i, (col, title, units) in enumerate([
    ('d_missense', 'Missense mutations Δ (post − pre)', 'count'),
    ('d_neo',      'MHC-I neoantigens Δ (post − pre)', 'count'),
    ('d_IGH',      'BCR (IGH) clonotypes Δ (post − pre)', 'count'),
    ('d_Treg',     'Treg signature Δ (post − pre)', 'z-score'),
]):
    ax = axes[ax_i]
    df_sub = D_ord.copy()
    df_sub = df_sub.sort_values(['response', col]).reset_index(drop=True)
    # Color dots by response
    colors = [PAL[r] for r in df_sub.response]
    x = np.arange(len(df_sub))
    ax.bar(x, df_sub[col], color=colors, edgecolor='white', lw=1.0, width=0.82)
    ax.axhline(0, color='#0e2a47', lw=0.9)
    ax.set_ylabel(title, fontsize=10, fontweight='bold', color='#0e2a47')
    add_axis_spines(ax)
    # Value labels on bars
    for i, v in enumerate(df_sub[col]):
        if pd.isna(v): continue
        if abs(v) > df_sub[col].abs().max() * 0.2:
            ax.text(i, v + (v*0.02 if v>0 else v*0.02-0.1), f'{v:+.0f}' if 'count' in units else f'{v:+.1f}',
                    ha='center', va='bottom' if v>0 else 'top',
                    fontsize=7, color='#0e2a47', fontweight='bold')
    ax.set_xlim(-0.6, len(df_sub)-0.4)

axes[-1].set_xticks(range(len(D_ord)))
axes[-1].set_xticklabels([f'S{int(s)}' for s in D_ord.subject_id], fontsize=9, color='#0e2a47')
axes[-1].set_xlabel('Subject', fontsize=11, fontweight='bold', color='#0e2a47')

# Legend
fig.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good responder'),
                    mpatches.Patch(color=BAD_DEEP, label='Poor responder')],
           loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=11, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig6D_waterfall', OUT)

# ============================================================
# 6E — Clonal stream plot (Oliveira Nature 2022 / fishplot-inspired)
# ============================================================
# PyClone cluster trajectories: for each paired subject, show cluster sizes pre→post as stacked "stream"
from pathlib import Path as _Path
pyclone_dir = ROOT/'04_wes_cnv_clonal/pyclone'
paired_subj_pyclone = pyclone.subject_id.tolist()

fig, axes = plt.subplots(3, 4, figsize=(14, 8.5), sharex=True, sharey=True)
axes_flat = axes.flatten()

# Sort subjects by response then subject_id
py_sorted = pyclone.sort_values(['response','subject_id']).reset_index(drop=True)

for i, (_, row) in enumerate(py_sorted.iterrows()):
    if i >= len(axes_flat): break
    ax = axes_flat[i]
    subj = int(row.subject_id)
    res_f = pyclone_dir/f'results_subj{subj}.tsv'
    if not res_f.exists():
        ax.set_visible(False); continue
    df = pd.read_csv(res_f, sep='\t')
    # cluster × sample → mean cellular prevalence
    cp = df.groupby(['cluster_id','sample_id']).cellular_prevalence.mean().unstack(fill_value=0)
    # Separate pre and post columns
    pre_col = [c for c in cp.columns if '-PR' in c]
    post_col = [c for c in cp.columns if '-PO' in c]
    if not pre_col or not post_col:
        ax.set_visible(False); continue
    cp_df = pd.DataFrame({'pre': cp[pre_col[0]], 'post': cp[post_col[0]]})
    # Normalize to fraction of total cells
    cp_df['pre'] = cp_df['pre'] / cp_df['pre'].sum() if cp_df['pre'].sum()>0 else 0
    cp_df['post'] = cp_df['post'] / cp_df['post'].sum() if cp_df['post'].sum()>0 else 0
    # Sort clusters by pre prevalence descending
    cp_df = cp_df.sort_values('pre', ascending=False)
    # Plot as stacked streamgraph (like a fishplot)
    x_pts = np.array([0, 1])
    colors_cluster = plt.cm.tab20(np.linspace(0, 1, len(cp_df)))
    bot_pre = 0; bot_post = 0
    for k, (cl, cp_r) in enumerate(cp_df.iterrows()):
        verts = [(0, bot_pre), (1, bot_post), (1, bot_post + cp_r.post), (0, bot_pre + cp_r.pre)]
        poly = Polygon(verts, facecolor=colors_cluster[k], edgecolor='white', lw=0.6, alpha=0.85)
        ax.add_patch(poly)
        bot_pre += cp_r.pre; bot_post += cp_r.post
    # Subject label with response color
    ax.text(0.5, 1.04, f'S{subj}', transform=ax.transAxes, ha='center', va='bottom',
            fontsize=10.5, color=PAL[row.response], fontweight='bold')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['pre','post'], fontsize=9, color='#0e2a47')
    ax.set_yticks([])
    add_axis_spines(ax, sides=('bottom',))
    ax.spines['left'].set_visible(False)
    # Annotate number of clusters
    ax.text(0.02, 0.98, f'k = {len(cp_df)}', transform=ax.transAxes, ha='left', va='top',
            fontsize=8, color='#0e2a47',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, boxstyle='round,pad=0.2'))

# Hide unused
for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)

fig.tight_layout()
save_panel(fig, 'Fig6E_clonal_stream', OUT)

# ============================================================
# 6F — Cascade schematic (illustrative, Mariathasan/Helmink style)
# ============================================================
fig, ax = plt.subplots(figsize=(13, 7))
ax.axis('off')

def _box(x, y, w, h, label, color, textcolor='white', fontsize=10.5):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.04",
                 facecolor=color, edgecolor='#0e2a47', linewidth=1.6))
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            color=textcolor, fontweight='bold')

def _arrow(x1, y1, x2, y2, color='#0e2a47', lw=2.0, label=''):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=lw, color=color, mutation_scale=22))
    if label:
        ax.text((x1+x2)/2, (y1+y2)/2+0.025, label, ha='center', fontsize=9,
                color=color, style='italic')

# Top: pre-treatment
_box(0.5, 0.92, 0.35, 0.09,
     'Pre-treatment tumor\nProliferative · DNA-repair-proficient · ±HLA-LOH · neoantigen-presenting',
     '#0e2a47', fontsize=10)

# TNT arrow
_arrow(0.5, 0.87, 0.5, 0.78, color=GOLD, lw=3.4)
ax.text(0.52, 0.82, 'TNT (FOLFOX/CAPOX + CRT)', fontsize=10.5, color='#0e2a47',
        fontweight='bold', va='center', ha='left')

# Cascade events in good responder (left branch)
events = [
    (0.09, 0.62, 'Mutation\nclearance\n$\\Delta$missense –67', '#118ab2'),
    (0.26, 0.62, 'Neoantigen\nclone loss\n$\\Delta$binders –312', '#06aed5'),
    (0.43, 0.62, 'HLA-LOH\nclone loss\nS3: 3→1, S4: 2→0', '#0f8b78'),
    (0.60, 0.62, 'Immune\nreprogramming\nTreg+MHC-II ↑', '#c11456'),
    (0.80, 0.62, 'Lymphocyte\ninfiltration\nIGH +1,424', '#7a3aad'),
]
for x, y, lbl, col in events:
    _box(x, y, 0.15, 0.15, lbl, col, fontsize=9)
    _arrow(0.5, 0.72, x, 0.70, color=col, lw=1.6)

# Output — good
_box(0.5, 0.28, 0.5, 0.09,
     'Pathologic complete / near-complete response (TRG 0–1)',
     GOOD_DEEP, fontsize=11)
for x, _, _, col in events:
    ax.annotate('', xy=(0.5, 0.33), xytext=(x, 0.54),
                arrowprops=dict(arrowstyle='->', lw=1.3, color=col, alpha=0.55, mutation_scale=16))

# Poor responder — separate box below
_box(0.5, 0.07, 0.55, 0.09,
     'Poor responder: minimal molecular / immune change → primary treatment insensitivity',
     BAD_DEEP, fontsize=10.5)

ax.set_xlim(0, 1); ax.set_ylim(0, 1)
save_panel(fig, 'Fig6F_cascade', OUT)

print('\n=== Fig 6 v3 (6 journal-style panels) saved ===')
