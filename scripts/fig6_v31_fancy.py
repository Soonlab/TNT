"""
Fig 6 v3.1 — substantial redesign per user feedback:
  6A: fix left labels (outside plot) + staggered metric headers with guide lines
  6B: fancy — slope + violin + data points + bracket stats (RainCloud-like)
  6D: y-labels shortened/repositioned, no overlap
  6E: real smooth-Bezier fishplot with cluster hierarchy (Oliveira Nature 2022 / fishplot R package style)
  6F: complex dual-cascade schematic with embedded mini bar charts + curved connectors
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import (Rectangle, FancyBboxPatch, FancyArrowPatch, Polygon,
                                 PathPatch, Circle, ConnectionPatch)
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap, to_rgb
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
sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_summary_key.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
pyclone = pd.read_csv(ROOT/'04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\t')

PAIRED_SUBJ = list(range(1, 15))

def build_delta():
    rows = []
    for s in PAIRED_SUBJ:
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        ts = tmb[(tmb.subject_id==s) & tmb.timepoint.isin(['pre','post'])]
        d_tmb = d_miss = np.nan
        if {'pre','post'}.issubset(set(ts.timepoint)):
            pre_t = ts[ts.timepoint=='pre'].iloc[0]; post_t = ts[ts.timepoint=='post'].iloc[0]
            d_tmb = post_t.TMB_nonsyn_per_Mb - pre_t.TMB_nonsyn_per_Mb
            d_miss = post_t.n_nonsyn - pre_t.n_nonsyn
        ss = sbs[(sbs.subject_id==s) & sbs.timepoint.isin(['pre','post'])]
        d_sbs = np.nan
        if {'pre','post'}.issubset(set(ss.timepoint)):
            d_sbs = ss[ss.timepoint=='post'].iloc[0].SBS5 - ss[ss.timepoint=='pre'].iloc[0].SBS5
        nd = neo_delta[neo_delta.subject_id==s]
        d_neo = nd.delta_binders.iloc[0] if len(nd) else np.nan
        loh_pre = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PR')]
        loh_post = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PO')]
        d_loh = (loh_post.LOH_call.sum() - loh_pre.LOH_call.sum()) if len(loh_pre)+len(loh_post)>0 else np.nan
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
# 6A — swimmer fixed (labels out of plot, staggered headers)
# ============================================================
D_sw = D.sort_values(['response','subject_id']).reset_index(drop=True)
N_sub = len(D_sw)

fig = plt.figure(figsize=(15, 8))
gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 6, 0.9], width_ratios=[0.10, 1.0],
                      hspace=0.02, wspace=0.02)

# Header axis (row 0, col 1) — metric labels with guide lines, staggered
ax_hdr = fig.add_subplot(gs[0, 1])
ax_hdr.axis('off')
ax_hdr.set_xlim(0, 1); ax_hdr.set_ylim(0, 1)

metrics = [
    ('d_TMB',      'TMB',          'clearance', 0.05, 0.90),
    ('d_missense', 'Missense',     'clearance', 0.15, 0.60),
    ('d_SBS5',     'SBS5 muts',    'clearance', 0.25, 0.90),
    ('d_neo',      'Neoantigens',  'clearance', 0.38, 0.60),
    ('d_HLA_LOH',  'HLA LOH',      'clearance', 0.51, 0.90),
    ('d_Treg',     'Treg',         'activation',0.62, 0.60),
    ('d_MHC_II',   'MHC-II',       'activation',0.72, 0.90),
    ('d_CD8_exh',  'CD8 exh.',     'activation',0.83, 0.60),
    ('d_IGH',      'BCR (IGH)',    'activation',0.94, 0.90),
]
for m, lbl, direction, x_pos, y_lbl in metrics:
    ax_hdr.text(x_pos, y_lbl, lbl, ha='center', va='bottom', fontsize=10,
                color='#0e2a47', fontweight='bold')
    # Guide line from header to top of plot
    ax_hdr.plot([x_pos, x_pos], [y_lbl-0.05, 0.05], color='#aab3bf', lw=0.7, alpha=0.6)

# Category labels above
ax_hdr.text(0.295, 1.05, 'Tumor clearance (post − pre < 0 = ↓tumor)', ha='center', fontsize=10.5,
            color='#0e2a47', fontweight='bold', style='italic',
            bbox=dict(facecolor='#e8f0f8', edgecolor='#118ab2', alpha=0.9, boxstyle='round,pad=0.3'))
ax_hdr.text(0.78, 1.05, 'Immune activation (post − pre > 0 = ↑immune)', ha='center', fontsize=10.5,
            color='#0e2a47', fontweight='bold', style='italic',
            bbox=dict(facecolor='#fff5e8', edgecolor='#d96125', alpha=0.9, boxstyle='round,pad=0.3'))

# Main plot (row 1, col 1)
ax = fig.add_subplot(gs[1, 1])
ax.set_xlim(0, 1); ax.set_ylim(-0.5, N_sub-0.5)
ax.invert_yaxis()

# Row bg strip
for i, (_, r) in enumerate(D_sw.iterrows()):
    ax.axhspan(i-0.4, i+0.4, facecolor=PAL[r.response], alpha=0.10, zorder=0)

# Vertical separator between tumor and immune sections
ax.axvline(0.58, color='#5a6772', lw=0.8, ls=':', alpha=0.55, zorder=0.5)

# Events
for i, (_, r) in enumerate(D_sw.iterrows()):
    for m, lbl, direction, x_pos, _ in metrics:
        v = r[m]
        if pd.isna(v): continue
        if direction=='clearance':
            color = '#118ab2' if v<0 else '#c11456'; mag = abs(v)
        else:
            color = '#0f8b78' if v>0 else '#d96125'; mag = abs(v)
        if m == 'd_neo': s = min(500, 50 + mag*0.8)
        elif m == 'd_IGH': s = min(500, 50 + mag*0.15)
        elif m == 'd_missense': s = min(500, 50 + mag*3)
        elif m == 'd_SBS5': s = min(500, 50 + mag*3)
        elif m == 'd_HLA_LOH': s = 80 + mag*80
        elif m == 'd_TMB': s = min(500, 50 + mag*40)
        else: s = min(500, 40 + mag*80)
        ax.scatter(x_pos, i, s=s, c=color, alpha=0.88, edgecolor='white', lw=1.2, zorder=5)

ax.set_yticks(range(N_sub))
ax.set_yticklabels([f"S{int(r.subject_id)}" for _, r in D_sw.iterrows()], fontsize=10, color='#0e2a47')
ax.set_xticks([])
for s in ['top','right','bottom']: ax.spines[s].set_visible(False)
ax.spines['left'].set_color('#0e2a47')

# Left column (row 1, col 0) — response bracket labels OUTSIDE plot
ax_left = fig.add_subplot(gs[1, 0])
ax_left.axis('off')
ax_left.set_xlim(0, 1); ax_left.set_ylim(0, N_sub)
ax_left.invert_yaxis()
good_n = (D_sw.response=='good').sum()
# Bracket - good
ax_left.add_patch(Rectangle((0.65, 0), 0.15, good_n-0.05, facecolor=GOOD_DEEP, alpha=0.2,
                             edgecolor=GOOD_DEEP, linewidth=1.5))
ax_left.text(0.35, good_n/2 - 0.5, 'Good\nresponders', ha='center', va='center',
             fontsize=11.5, color=GOOD_DEEP, fontweight='bold', rotation=90)
# Bracket - bad
ax_left.add_patch(Rectangle((0.65, good_n+0.05), 0.15, N_sub-good_n-0.1, facecolor=BAD_DEEP, alpha=0.2,
                             edgecolor=BAD_DEEP, linewidth=1.5))
ax_left.text(0.35, good_n + (N_sub-good_n)/2 - 0.5, 'Poor\nresponders',
             ha='center', va='center', fontsize=11.5, color=BAD_DEEP, fontweight='bold', rotation=90)

# Bottom legend
ax_leg = fig.add_subplot(gs[2, 1])
ax_leg.axis('off')
leg_elems = [
    plt.scatter([],[], s=180, c='#118ab2', edgecolor='white', lw=1.2, label='↓ Tumor  (Δ < 0, favorable)'),
    plt.scatter([],[], s=180, c='#c11456', edgecolor='white', lw=1.2, label='↑ Tumor  (Δ > 0, unfavorable)'),
    plt.scatter([],[], s=180, c='#0f8b78', edgecolor='white', lw=1.2, label='↑ Immune  (Δ > 0, favorable)'),
    plt.scatter([],[], s=180, c='#d96125', edgecolor='white', lw=1.2, label='↓ Immune  (Δ < 0, unfavorable)'),
]
ax_leg.legend(handles=leg_elems, loc='upper center', ncol=4, fontsize=10, frameon=False,
              bbox_to_anchor=(0.5, 1.0))

save_panel(fig, 'Fig6A_swimmer', OUT)

# ============================================================
# 6B — fancy slope + violin + stats (RainCloud-like)
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(16, 5.2), sharex=False)

def fancy_slope_panel(ax, pre_map, post_map, ylabel, ylabel_units=''):
    x_pre = 0.0; x_post = 1.0
    v_pre_good = []; v_pre_bad = []; v_post_good = []; v_post_bad = []
    for s in PAIRED_SUBJ:
        if s not in pre_map or s not in post_map: continue
        vp = pre_map[s]; vpo = post_map[s]
        if pd.isna(vp) or pd.isna(vpo): continue
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        if resp=='good': v_pre_good.append(vp); v_post_good.append(vpo)
        else: v_pre_bad.append(vp); v_post_bad.append(vpo)

    # Half violins at pre (left) and post (right)
    for xpos, vals_g, vals_b, side in [(x_pre, v_pre_good, v_pre_bad, 'left'),
                                         (x_post, v_post_good, v_post_bad, 'right')]:
        for vals, color, offset in [(vals_g, GOOD_DEEP, -0.10 if side=='left' else 0.10),
                                      (vals_b, BAD_DEEP, -0.10 if side=='left' else 0.10)]:
            if len(vals) < 2: continue
            parts = ax.violinplot([vals], positions=[xpos + offset],
                                   widths=0.22, showmeans=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(color); pc.set_alpha(0.25); pc.set_edgecolor(color); pc.set_linewidth(1.0)
                # Clip to one side
                m = np.mean(pc.get_paths()[0].vertices[:,0])
                if side=='left':
                    pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], -np.inf, m)
                else:
                    pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], m, np.inf)

    # Connecting slope lines
    for s in PAIRED_SUBJ:
        if s not in pre_map or s not in post_map: continue
        vp = pre_map[s]; vpo = post_map[s]
        if pd.isna(vp) or pd.isna(vpo): continue
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        color = PAL[resp]
        # Slightly curved line using FancyArrowPatch
        arr = FancyArrowPatch((x_pre+0.02, vp), (x_post-0.02, vpo),
                               arrowstyle='-|>', color=color, lw=1.8, alpha=0.75,
                               mutation_scale=12, connectionstyle='arc3,rad=0.08')
        ax.add_patch(arr)
        ax.scatter(x_pre+0.02, vp, s=55, c=color, edgecolor='white', lw=1.0, zorder=4)
        ax.scatter(x_post-0.02, vpo, s=55, c=color, edgecolor='white', lw=1.0, zorder=4)

    # Group means as thick lines (bold emphasis)
    for xpos, vals_g, vals_b in [(x_pre, v_pre_good, v_pre_bad), (x_post, v_post_good, v_post_bad)]:
        if vals_g:
            m = np.median(vals_g)
            ax.plot([xpos-0.16, xpos+0.16], [m, m], color=GOOD_DEEP, lw=3.5, solid_capstyle='round', zorder=5)
        if vals_b:
            m = np.median(vals_b)
            ax.plot([xpos-0.16, xpos+0.16], [m, m], color=BAD_DEEP, lw=3.5, solid_capstyle='round', zorder=5)

    # Within-group paired Wilcoxon
    from scipy.stats import wilcoxon
    try:
        if len(v_pre_good) >= 3:
            p_g = wilcoxon(v_pre_good, v_post_good).pvalue
        else: p_g = np.nan
    except: p_g = np.nan
    try:
        if len(v_pre_bad) >= 3:
            p_b = wilcoxon(v_pre_bad, v_post_bad).pvalue
        else: p_b = np.nan
    except: p_b = np.nan

    ax.set_xticks([x_pre, x_post])
    ax.set_xticklabels(['pre', 'post'], fontsize=11.5, color='#0e2a47', fontweight='bold')
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(ylabel + (f' ({ylabel_units})' if ylabel_units else ''),
                  fontsize=11, fontweight='bold', color='#0e2a47')
    add_axis_spines(ax)
    # Stats annotation — two inline p-values
    p_g_str = f'{p_g:.2g}' if not pd.isna(p_g) else 'NA'
    p_b_str = f'{p_b:.2g}' if not pd.isna(p_b) else 'NA'
    txt = f'Good (paired) p = {p_g_str}\nPoor (paired) p = {p_b_str}'
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha='right', va='top', fontsize=9,
            color='#0e2a47', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.92, boxstyle='round,pad=0.35'))

# Metric maps
tmb_pre = {s: tmb[tmb.sample_id==f'{s}-PR'].TMB_nonsyn_per_Mb.iloc[0]
           if (tmb.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ}
tmb_post = {s: tmb[tmb.sample_id==f'{s}-PO'].TMB_nonsyn_per_Mb.iloc[0]
            if (tmb.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ}
fancy_slope_panel(axes[0], tmb_pre, tmb_post, 'TMB', '/Mb')

neo_pre = {s: neo[neo.sample_id==f'{s}-PR'].n_binders_500nM.iloc[0]
           if (neo.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ}
neo_post = {s: neo[neo.sample_id==f'{s}-PO'].n_binders_500nM.iloc[0]
            if (neo.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ}
fancy_slope_panel(axes[1], neo_pre, neo_post, 'MHC-I binders', '<500 nM')

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
fancy_slope_panel(axes[2], igh_pre, igh_post, 'BCR (IGH) clonotypes', 'count')

def sig_get(sid, col):
    if sid is None or sid not in sigs.index: return np.nan
    return sigs.loc[sid, col]
treg_pre = {s: sig_get(rna_pre_ids[s], 'Treg') for s in PAIRED_SUBJ}
treg_post = {s: sig_get(rna_post_ids[s], 'Treg') for s in PAIRED_SUBJ}
fancy_slope_panel(axes[3], treg_pre, treg_post, 'Treg signature', 'z-score')

fig.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good responders'),
                    mpatches.Patch(color=BAD_DEEP, label='Poor responders'),
                    matplotlib.lines.Line2D([0],[0], color='#5a6772', lw=3.5, label='Group median')],
           loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=3, fontsize=11, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig6B_slope_fancy', OUT)

# ============================================================
# 6D — waterfall with shortened y-labels
# ============================================================
fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
D_ord = D.sort_values(['response','d_missense']).reset_index(drop=True)
for ax_i, (col, short, full, units) in enumerate([
    ('d_missense', 'Missense Δ',  'Missense mutations Δ (post − pre)', 'count'),
    ('d_neo',      'Neoantigen Δ','MHC-I neoantigens Δ (post − pre)', 'count'),
    ('d_IGH',      'IGH Δ',       'BCR (IGH) clonotypes Δ (post − pre)', 'count'),
    ('d_Treg',     'Treg Δ',      'Treg signature Δ (post − pre)', 'z-score'),
]):
    ax = axes[ax_i]
    df_sub = D_ord.copy().sort_values(['response', col]).reset_index(drop=True)
    colors = [PAL[r] for r in df_sub.response]
    x = np.arange(len(df_sub))
    ax.bar(x, df_sub[col], color=colors, edgecolor='white', lw=1.0, width=0.82)
    ax.axhline(0, color='#0e2a47', lw=0.9)
    # Short label left, full label as text above
    ax.set_ylabel(short, fontsize=11, fontweight='bold', color='#0e2a47')
    ax.text(0.02, 0.95, full, transform=ax.transAxes, ha='left', va='top',
            fontsize=10, color='#0e2a47', style='italic',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round,pad=0.3'))
    add_axis_spines(ax)
    for i, v in enumerate(df_sub[col]):
        if pd.isna(v): continue
        if abs(v) > df_sub[col].abs().max() * 0.2:
            fmt = f'{v:+.0f}' if 'count' in units else f'{v:+.1f}'
            ax.text(i, v + (v*0.02 if v>0 else v*0.02-0.1), fmt,
                    ha='center', va='bottom' if v>0 else 'top',
                    fontsize=8, color='#0e2a47', fontweight='bold')
    ax.set_xlim(-0.6, len(df_sub)-0.4)

axes[-1].set_xticks(range(len(D_ord)))
axes[-1].set_xticklabels([f'S{int(s)}' for s in D_ord.subject_id], fontsize=9.5, color='#0e2a47')
axes[-1].set_xlabel('Subject', fontsize=11.5, fontweight='bold', color='#0e2a47')
fig.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good responder'),
                    mpatches.Patch(color=BAD_DEEP, label='Poor responder')],
           loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=11, frameon=False)
fig.tight_layout()
save_panel(fig, 'Fig6D_waterfall', OUT)

# ============================================================
# 6E — FANCY fishplot (Bezier-curve stream with cluster hierarchy)
# ============================================================
pyclone_dir = ROOT/'04_wes_cnv_clonal/pyclone'
py_sorted = pyclone.sort_values(['response','subject_id']).reset_index(drop=True)

fig, axes = plt.subplots(3, 4, figsize=(15.5, 9), sharex=True, sharey=True)
axes_flat = axes.flatten()

def bezier_flow(ax, x0, x1, y0_bot, y0_top, y1_bot, y1_top, color, alpha=0.85):
    """Smooth stream polygon using cubic bezier curves."""
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0_bot),
        (cx, y0_bot), (cx, y1_bot), (x1, y1_bot),  # bottom curve
        (x1, y1_top),
        (cx, y1_top), (cx, y0_top), (x0, y0_top),  # top curve
        (x0, y0_bot)
    ]
    codes = [MPath.MOVETO,
             MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.LINETO,
             MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.CLOSEPOLY]
    path = MPath(verts, codes)
    patch = PathPatch(path, facecolor=color, edgecolor='white', lw=1.0, alpha=alpha)
    ax.add_patch(patch)

for i, (_, row) in enumerate(py_sorted.iterrows()):
    if i >= len(axes_flat): break
    ax = axes_flat[i]
    subj = int(row.subject_id)
    res_f = pyclone_dir/f'results_subj{subj}.tsv'
    if not res_f.exists():
        ax.set_visible(False); continue
    df = pd.read_csv(res_f, sep='\t')
    cp = df.groupby(['cluster_id','sample_id']).cellular_prevalence.mean().unstack(fill_value=0)
    pre_col = [c for c in cp.columns if '-PR' in c]
    post_col = [c for c in cp.columns if '-PO' in c]
    if not pre_col or not post_col:
        ax.set_visible(False); continue
    cp_df = pd.DataFrame({'pre': cp[pre_col[0]], 'post': cp[post_col[0]]})
    # Normalize
    cp_df['pre'] = cp_df['pre'] / max(cp_df['pre'].sum(), 1e-9)
    cp_df['post'] = cp_df['post'] / max(cp_df['post'].sum(), 1e-9)
    cp_df = cp_df.sort_values('pre', ascending=False)
    # Cluster colors — varied palette
    n_cl = len(cp_df)
    palette_cl = plt.cm.Spectral(np.linspace(0.1, 0.9, n_cl))

    x0, x1 = 0.0, 1.0
    bot_pre = 0; bot_post = 0
    for k, (cl_id, cp_r) in enumerate(cp_df.iterrows()):
        y0b = bot_pre; y0t = bot_pre + cp_r.pre
        y1b = bot_post; y1t = bot_post + cp_r.post
        bezier_flow(ax, x0, x1, y0b, y0t, y1b, y1t, palette_cl[k], alpha=0.92)
        bot_pre += cp_r.pre; bot_post += cp_r.post

    # Subject title
    ax.text(0.5, 1.06, f'S{subj}', transform=ax.transAxes, ha='center', va='bottom',
            fontsize=11.5, color=PAL[row.response], fontweight='bold')
    # Key stats: dominant shrink + expand
    ax.text(0.02, 0.98, f'k = {int(row.n_clusters)}\n↓{int(row.n_shrinking)}  ↑{int(row.n_expanding)}',
            transform=ax.transAxes, ha='left', va='top', fontsize=8, color='#0e2a47',
            fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.25'))
    # Axis
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['pre','post'], fontsize=9.5, color='#0e2a47', fontweight='bold')
    ax.set_yticks([])
    for s in ['top','right','left']: ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color('#0e2a47')

for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)

fig.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good responder subjects'),
                    mpatches.Patch(color=BAD_DEEP, label='Poor responder subjects')],
           loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=11, frameon=False)
# Common y-axis label using figure text
fig.text(0.005, 0.5, 'Clonal composition (fraction of tumor cells)', rotation=90,
         va='center', fontsize=11, fontweight='bold', color='#0e2a47')
fig.tight_layout()
save_panel(fig, 'Fig6E_fishplot', OUT)

# ============================================================
# 6F — complex dual-pathway cascade with mini bar charts
# ============================================================
fig = plt.figure(figsize=(16, 9))
gs_main = fig.add_gridspec(1, 1)
ax = fig.add_subplot(gs_main[0])
ax.axis('off')
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def fancy_box(x, y, w, h, label, color, textcolor='white', fontsize=10, style='round'):
    if style=='round':
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                     boxstyle="round,pad=0.005,rounding_size=0.025",
                     facecolor=color, edgecolor='#0e2a47', linewidth=1.6))
    else:
        ax.add_patch(Rectangle((x-w/2, y-h/2), w, h, facecolor=color,
                                edgecolor='#0e2a47', linewidth=1.6))
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            color=textcolor, fontweight='bold')

def mini_bar(cx, cy, w, h, d_good, d_bad, ylabel, units='', favorable='neg'):
    """Mini inset bar chart at (cx, cy) center, width w, height h, showing good vs bad median Δ."""
    # Box background
    ax.add_patch(Rectangle((cx-w/2, cy-h/2), w, h, facecolor='white',
                            edgecolor='#5a6772', linewidth=0.8, alpha=0.95))
    # Two bars at x positions within the mini axes
    bar_w = w * 0.25
    bar_max = max(abs(d_good), abs(d_bad)) * 1.3 if max(abs(d_good),abs(d_bad))>0 else 1
    y_base = cy - h*0.12
    # Good bar
    hg = (d_good / bar_max) * (h * 0.35)
    ax.add_patch(Rectangle((cx - w*0.28 - bar_w/2, y_base), bar_w, hg,
                            facecolor=GOOD_DEEP, alpha=0.9, edgecolor='white'))
    # Bad bar
    hb = (d_bad / bar_max) * (h * 0.35)
    ax.add_patch(Rectangle((cx + w*0.08 - bar_w/2, y_base), bar_w, hb,
                            facecolor=BAD_DEEP, alpha=0.9, edgecolor='white'))
    # Zero line
    ax.plot([cx - w*0.4, cx + w*0.4], [y_base, y_base], color='#0e2a47', lw=0.8)
    # Labels (good/bad values)
    ax.text(cx - w*0.28, cy + h*0.32, f'{d_good:+.0f}' if 'count' in units else f'{d_good:+.1f}',
            ha='center', fontsize=7.5, color=GOOD_DEEP, fontweight='bold')
    ax.text(cx + w*0.08, cy + h*0.32, f'{d_bad:+.0f}' if 'count' in units else f'{d_bad:+.1f}',
            ha='center', fontsize=7.5, color=BAD_DEEP, fontweight='bold')
    # Y label at top of mini
    ax.text(cx, cy + h*0.48, ylabel, ha='center', fontsize=8, fontweight='bold', color='#0e2a47')
    # G / B labels below bars
    ax.text(cx - w*0.28, y_base - h*0.06, 'G', ha='center', fontsize=7, color=GOOD_DEEP, fontweight='bold')
    ax.text(cx + w*0.08, y_base - h*0.06, 'B', ha='center', fontsize=7, color=BAD_DEEP, fontweight='bold')

def curve_arrow(x1, y1, x2, y2, color='#0e2a47', lw=1.8, style='-|>', mutation=18, rad=0.0):
    connectionstyle = f'arc3,rad={rad}'
    arr = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, color=color,
                           lw=lw, alpha=0.95, mutation_scale=mutation,
                           connectionstyle=connectionstyle)
    ax.add_patch(arr)

# ===== TITLE REMOVED (per rule) =====

# Top: Pre-treatment tumor (shared)
fancy_box(0.5, 0.91, 0.35, 0.10,
          'Pre-treatment tumor\n(proliferative · DNA-repair-proficient · HLA-LOH± · neoantigen-presenting)',
          '#0e2a47', fontsize=10)

# TNT arrow
curve_arrow(0.5, 0.85, 0.5, 0.77, color=GOLD, lw=4.0, mutation=22)
ax.text(0.53, 0.81, 'TNT\n(FOLFOX/CAPOX + CRT)', fontsize=10, color='#0e2a47',
        fontweight='bold', va='center', ha='left')

# Branch label
ax.text(0.24, 0.69, 'Good responders', fontsize=12, color=GOOD_DEEP, fontweight='bold', ha='center',
        bbox=dict(facecolor='white', edgecolor=GOOD_DEEP, linewidth=1.5, boxstyle='round,pad=0.3'))
ax.text(0.76, 0.69, 'Poor responders', fontsize=12, color=BAD_DEEP, fontweight='bold', ha='center',
        bbox=dict(facecolor='white', edgecolor=BAD_DEEP, linewidth=1.5, boxstyle='round,pad=0.3'))

# Branching arrows
curve_arrow(0.45, 0.73, 0.25, 0.68, color=GOOD_DEEP, lw=2.5, rad=-0.15)
curve_arrow(0.55, 0.73, 0.75, 0.68, color=BAD_DEEP, lw=2.5, rad=0.15)

# Good responder cascade events (stacked down the left) with mini bar charts
d_miss_g = D[D.response=='good'].d_missense.median()
d_miss_b = D[D.response=='bad'].d_missense.median()
d_neo_g = D[D.response=='good'].d_neo.median()
d_neo_b = D[D.response=='bad'].d_neo.median()
d_sbs_g = D[D.response=='good'].d_SBS5.median()
d_sbs_b = D[D.response=='bad'].d_SBS5.median()
d_igh_g = D[D.response=='good'].d_IGH.median()
d_igh_b = D[D.response=='bad'].d_IGH.median()
d_treg_g = D[D.response=='good'].d_Treg.median()
d_treg_b = D[D.response=='bad'].d_Treg.median()
d_mhc_g = D[D.response=='good'].d_MHC_II.median()
d_mhc_b = D[D.response=='bad'].d_MHC_II.median()

# Mini panels layout (good path, curving down-left)
events_good = [
    (0.12, 0.58, 'Mutation\nclearance', '#118ab2', d_miss_g, d_miss_b, 'Δ missense', 'count'),
    (0.12, 0.45, 'Neoantigen\nclone loss', '#06aed5', d_neo_g, d_neo_b, 'Δ binders', 'count'),
    (0.12, 0.32, 'SBS5 mutation\nclearance', '#0f8b78', d_sbs_g, d_sbs_b, 'Δ SBS5', 'count'),
    (0.12, 0.19, 'Treg / MHC-II\nreprogramming', '#c11456', d_treg_g, d_treg_b, 'Δ Treg z', 'z'),
    (0.30, 0.13, 'BCR (IGH)\ninfiltration', '#7a3aad', d_igh_g, d_igh_b, 'Δ IGH', 'count'),
]
# Draw good cascade with boxes + mini bars
prev_x, prev_y = 0.25, 0.67
for cx, cy, lbl, col, dg, db, ylabel, units in events_good:
    # Event box
    fancy_box(cx, cy, 0.14, 0.08, lbl, col, textcolor='white', fontsize=9)
    # Mini bar to its right
    mini_bar(cx+0.10, cy, 0.11, 0.075, dg, db, ylabel, units)
    # Connecting curved arrow from previous
    curve_arrow(prev_x, prev_y-0.035, cx, cy+0.04, color=col, lw=1.8, rad=-0.3)
    prev_x, prev_y = cx, cy

# Good → outcome
fancy_box(0.35, 0.04, 0.42, 0.07,
          'Complete / near-complete response (TRG 0–1)', GOOD_DEEP, fontsize=11)
curve_arrow(prev_x+0.08, prev_y-0.035, 0.33, 0.08, color=GOOD_DEEP, lw=2.3, rad=0.15)

# Poor responder path (right side, static)
fancy_box(0.76, 0.50, 0.22, 0.20,
          'Minimal molecular change\n\n• Missense Δ ≈ −9\n• Neoantigen Δ ≈ −100\n• Treg Δ ≈ 0\n• No B-cell influx\n\n→ Primary resistance',
          BAD_DEEP, fontsize=9.5)
curve_arrow(0.76, 0.63, 0.76, 0.61, color=BAD_DEEP, lw=2.0, rad=0)
fancy_box(0.76, 0.26, 0.30, 0.07,
          'Partial / minimal response (TRG 2–3)', BAD_DEEP, fontsize=11)
curve_arrow(0.76, 0.38, 0.76, 0.30, color=BAD_DEEP, lw=2.0, rad=0)

# Good vs Bad summary divider
ax.plot([0.55, 0.55], [0.02, 0.68], color='#aab3bf', lw=0.8, ls='--', alpha=0.6)

save_panel(fig, 'Fig6F_cascade', OUT)

print('\n=== Fig 6 v3.1 fancy revisions saved ===')
