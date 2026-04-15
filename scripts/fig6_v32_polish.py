"""
Fig 6 v3.2 — targeted fixes:
  6A: category boxes (tumor/immune) pushed higher, Good/Poor brackets moved further left
  6B: p-value text moved OUTSIDE top of data (ylim expansion)
  6F: clean dual-pathway layout — wider figure, proper spacing, no overlap
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import to_rgb
import warnings; warnings.filterwarnings('ignore')
from scipy.stats import wilcoxon

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
# 6A — push category boxes higher + widen left bracket area
# ============================================================
D_sw = D.sort_values(['response','subject_id']).reset_index(drop=True)
N_sub = len(D_sw)

fig = plt.figure(figsize=(16, 8.5))
# Taller header row, wider left col
gs = fig.add_gridspec(3, 2, height_ratios=[1.7, 6, 0.9], width_ratios=[0.16, 1.0],
                      hspace=0.02, wspace=0.02)

# ----- Header axis -----
ax_hdr = fig.add_subplot(gs[0, 1])
ax_hdr.axis('off')
ax_hdr.set_xlim(0, 1); ax_hdr.set_ylim(0, 1)

# Category boxes moved to TOP of header (y ~ 0.92)
ax_hdr.text(0.295, 0.92, 'Tumor clearance (post − pre < 0 = ↓tumor)', ha='center', fontsize=11,
            color='#0e2a47', fontweight='bold', style='italic',
            bbox=dict(facecolor='#e8f0f8', edgecolor='#118ab2', alpha=0.95, boxstyle='round,pad=0.35'))
ax_hdr.text(0.78, 0.92, 'Immune activation (post − pre > 0 = ↑immune)', ha='center', fontsize=11,
            color='#0e2a47', fontweight='bold', style='italic',
            bbox=dict(facecolor='#fff5e8', edgecolor='#d96125', alpha=0.95, boxstyle='round,pad=0.35'))

# Metric labels middle of header (y ~ 0.45, staggered 0.30 / 0.55)
metrics = [
    ('d_TMB',      'TMB',          'clearance', 0.05, 0.55),
    ('d_missense', 'Missense',     'clearance', 0.15, 0.30),
    ('d_SBS5',     'SBS5 muts',    'clearance', 0.25, 0.55),
    ('d_neo',      'Neoantigens',  'clearance', 0.38, 0.30),
    ('d_HLA_LOH',  'HLA LOH',      'clearance', 0.51, 0.55),
    ('d_Treg',     'Treg',         'activation',0.62, 0.30),
    ('d_MHC_II',   'MHC-II',       'activation',0.72, 0.55),
    ('d_CD8_exh',  'CD8 exh.',     'activation',0.83, 0.30),
    ('d_IGH',      'BCR (IGH)',    'activation',0.94, 0.55),
]
for m, lbl, direction, x_pos, y_lbl in metrics:
    ax_hdr.text(x_pos, y_lbl, lbl, ha='center', va='bottom', fontsize=10,
                color='#0e2a47', fontweight='bold')
    # Guide line from label DOWN to plot top
    ax_hdr.plot([x_pos, x_pos], [y_lbl-0.05, -0.02], color='#aab3bf', lw=0.7, alpha=0.6)

# ----- Main plot -----
ax = fig.add_subplot(gs[1, 1])
ax.set_xlim(0, 1); ax.set_ylim(-0.5, N_sub-0.5)
ax.invert_yaxis()
for i, (_, r) in enumerate(D_sw.iterrows()):
    ax.axhspan(i-0.4, i+0.4, facecolor=PAL[r.response], alpha=0.10, zorder=0)
ax.axvline(0.58, color='#5a6772', lw=0.8, ls=':', alpha=0.55, zorder=0.5)

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

# ----- Left bracket (wider + further left) -----
ax_left = fig.add_subplot(gs[1, 0])
ax_left.axis('off')
ax_left.set_xlim(0, 1); ax_left.set_ylim(0, N_sub)
ax_left.invert_yaxis()
good_n = (D_sw.response=='good').sum()

# Good bracket — moved to LEFT edge (x=0.05 — 0.22)
ax_left.add_patch(Rectangle((0.55, 0), 0.25, good_n-0.05, facecolor=GOOD_DEEP, alpha=0.22,
                             edgecolor=GOOD_DEEP, linewidth=1.8))
ax_left.text(0.22, good_n/2 - 0.5, 'Good\nresponders', ha='center', va='center',
             fontsize=11.5, color=GOOD_DEEP, fontweight='bold', rotation=90)
# Bad bracket
ax_left.add_patch(Rectangle((0.55, good_n+0.05), 0.25, N_sub-good_n-0.1, facecolor=BAD_DEEP, alpha=0.22,
                             edgecolor=BAD_DEEP, linewidth=1.8))
ax_left.text(0.22, good_n + (N_sub-good_n)/2 - 0.5, 'Poor\nresponders',
             ha='center', va='center', fontsize=11.5, color=BAD_DEEP, fontweight='bold', rotation=90)

# ----- Bottom legend -----
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
# 6B — fancy slope with p-value OUTSIDE top (ylim expanded)
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(16, 5.5), sharex=False)

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

    # Determine y-range for proper ylim (need space ABOVE data for p-value)
    all_vals = v_pre_good + v_pre_bad + v_post_good + v_post_bad
    if not all_vals:
        return
    y_max = max(all_vals); y_min = min(all_vals)
    y_range = y_max - y_min if y_max != y_min else 1
    # Expand ylim: top +30% for p-value; bottom +5%
    ax.set_ylim(y_min - y_range*0.05, y_max + y_range*0.40)

    # Half violins
    for xpos, vals_g, vals_b, side in [(x_pre, v_pre_good, v_pre_bad, 'left'),
                                         (x_post, v_post_good, v_post_bad, 'right')]:
        for vals, color in [(vals_g, GOOD_DEEP), (vals_b, BAD_DEEP)]:
            if len(vals) < 2: continue
            offset = -0.10 if side=='left' else 0.10
            parts = ax.violinplot([vals], positions=[xpos + offset],
                                   widths=0.22, showmeans=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(color); pc.set_alpha(0.25); pc.set_edgecolor(color); pc.set_linewidth(1.0)
                m = np.mean(pc.get_paths()[0].vertices[:,0])
                if side=='left':
                    pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], -np.inf, m)
                else:
                    pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], m, np.inf)

    # Slope arrows
    for s in PAIRED_SUBJ:
        if s not in pre_map or s not in post_map: continue
        vp = pre_map[s]; vpo = post_map[s]
        if pd.isna(vp) or pd.isna(vpo): continue
        resp = clin[clin.subject_id==s].response_bin.iloc[0]
        color = PAL[resp]
        arr = FancyArrowPatch((x_pre+0.02, vp), (x_post-0.02, vpo),
                               arrowstyle='-|>', color=color, lw=1.8, alpha=0.75,
                               mutation_scale=12, connectionstyle='arc3,rad=0.08')
        ax.add_patch(arr)
        ax.scatter(x_pre+0.02, vp, s=55, c=color, edgecolor='white', lw=1.0, zorder=4)
        ax.scatter(x_post-0.02, vpo, s=55, c=color, edgecolor='white', lw=1.0, zorder=4)

    # Median bars
    for xpos, vals_g, vals_b in [(x_pre, v_pre_good, v_pre_bad), (x_post, v_post_good, v_post_bad)]:
        if vals_g:
            m = np.median(vals_g)
            ax.plot([xpos-0.16, xpos+0.16], [m, m], color=GOOD_DEEP, lw=3.5, solid_capstyle='round', zorder=5)
        if vals_b:
            m = np.median(vals_b)
            ax.plot([xpos-0.16, xpos+0.16], [m, m], color=BAD_DEEP, lw=3.5, solid_capstyle='round', zorder=5)

    # Paired Wilcoxon p-values
    try:
        p_g = wilcoxon(v_pre_good, v_post_good).pvalue if len(v_pre_good)>=3 else np.nan
    except: p_g = np.nan
    try:
        p_b = wilcoxon(v_pre_bad, v_post_bad).pvalue if len(v_pre_bad)>=3 else np.nan
    except: p_b = np.nan

    ax.set_xticks([x_pre, x_post])
    ax.set_xticklabels(['pre', 'post'], fontsize=11.5, color='#0e2a47', fontweight='bold')
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(ylabel + (f' ({ylabel_units})' if ylabel_units else ''),
                  fontsize=11, fontweight='bold', color='#0e2a47')
    add_axis_spines(ax)
    # p-value text at TOP — now has room (ylim expanded)
    p_g_str = f'{p_g:.2g}' if not pd.isna(p_g) else 'NA'
    p_b_str = f'{p_b:.2g}' if not pd.isna(p_b) else 'NA'
    txt = f'Good (paired) p = {p_g_str}\nPoor (paired) p = {p_b_str}'
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha='right', va='top', fontsize=9,
            color='#0e2a47', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.95, boxstyle='round,pad=0.35'))

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
# 6F — clean dual-pathway cascade, no overlaps
# ============================================================
# New layout principle:
#   Top: pre-treatment shared
#   TNT arrow
#   Branch: GOOD (left column) vs POOR (right column)
#   Each cascade event in GOOD column: box + mini bar integrated into a row-based grid
#   Outcome boxes at bottom
#   No backward arrows

fig = plt.figure(figsize=(18, 11))
gs = fig.add_gridspec(1, 1)
ax = fig.add_subplot(gs[0])
ax.axis('off')
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def fancy_box(x, y, w, h, label, color, textcolor='white', fontsize=10):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle="round,pad=0.005,rounding_size=0.025",
                 facecolor=color, edgecolor='#0e2a47', linewidth=1.6))
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            color=textcolor, fontweight='bold')

def curve_arrow(x1, y1, x2, y2, color='#0e2a47', lw=1.8, style='-|>', mutation=18, rad=0.0, alpha=0.95):
    arr = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, color=color,
                           lw=lw, alpha=alpha, mutation_scale=mutation,
                           connectionstyle=f'arc3,rad={rad}')
    ax.add_patch(arr)

def mini_bar(cx, cy, w, h, d_good, d_bad, title, units=''):
    """Clean mini bar chart."""
    ax.add_patch(Rectangle((cx-w/2, cy-h/2), w, h, facecolor='white',
                            edgecolor='#5a6772', linewidth=0.9, alpha=0.95))
    ax.text(cx, cy + h*0.38, title, ha='center', fontsize=8.5, fontweight='bold', color='#0e2a47')
    mag = max(abs(d_good), abs(d_bad)) * 1.3 if max(abs(d_good), abs(d_bad))>0 else 1
    bar_w = w * 0.22
    y_base = cy - h*0.05
    hg = (d_good / mag) * (h * 0.30)
    hb = (d_bad / mag) * (h * 0.30)
    ax.add_patch(Rectangle((cx - w*0.22 - bar_w/2, y_base), bar_w, hg,
                            facecolor=GOOD_DEEP, alpha=0.9, edgecolor='white'))
    ax.add_patch(Rectangle((cx + w*0.12 - bar_w/2, y_base), bar_w, hb,
                            facecolor=BAD_DEEP, alpha=0.9, edgecolor='white'))
    ax.plot([cx - w*0.42, cx + w*0.42], [y_base, y_base], color='#0e2a47', lw=0.8)
    # Value labels — above bars
    fmt = (lambda v: f'{v:+.0f}') if 'count' in units else (lambda v: f'{v:+.1f}')
    y_lbl_g = y_base + hg + h*0.02 if hg>=0 else y_base + hg - h*0.04
    y_lbl_b = y_base + hb + h*0.02 if hb>=0 else y_base + hb - h*0.04
    ax.text(cx - w*0.22, y_lbl_g, fmt(d_good), ha='center',
            va='bottom' if hg>=0 else 'top',
            fontsize=7.5, color=GOOD_DEEP, fontweight='bold')
    ax.text(cx + w*0.12, y_lbl_b, fmt(d_bad), ha='center',
            va='bottom' if hb>=0 else 'top',
            fontsize=7.5, color=BAD_DEEP, fontweight='bold')
    # G/B labels
    ax.text(cx - w*0.22, y_base - h*0.10, 'Good', ha='center', fontsize=7,
            color=GOOD_DEEP, fontweight='bold')
    ax.text(cx + w*0.12, y_base - h*0.10, 'Poor', ha='center', fontsize=7,
            color=BAD_DEEP, fontweight='bold')

# Compute medians
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

# ===== TOP: pre-treatment shared =====
fancy_box(0.5, 0.93, 0.45, 0.08,
          'Pre-treatment LARC tumor\n(proliferative · DNA-repair-proficient · HLA-LOH± · neoantigen-presenting)',
          '#0e2a47', fontsize=10.5)

# TNT arrow down to branching point
curve_arrow(0.5, 0.88, 0.5, 0.82, color=GOLD, lw=4.0, mutation=24)
ax.text(0.535, 0.85, 'TNT\n(FOLFOX / CAPOX + CRT)',
        fontsize=10.5, color='#0e2a47', fontweight='bold', va='center', ha='left')

# ===== Branch labels =====
ax.text(0.26, 0.76, 'Good responders', fontsize=13, color='white', fontweight='bold', ha='center',
        bbox=dict(facecolor=GOOD_DEEP, edgecolor='#0e2a47', linewidth=1.6, boxstyle='round,pad=0.35'))
ax.text(0.74, 0.76, 'Poor responders', fontsize=13, color='white', fontweight='bold', ha='center',
        bbox=dict(facecolor=BAD_DEEP, edgecolor='#0e2a47', linewidth=1.6, boxstyle='round,pad=0.35'))

# Branching arrows from TNT point
curve_arrow(0.47, 0.81, 0.27, 0.77, color=GOOD_DEEP, lw=2.6, rad=-0.18)
curve_arrow(0.53, 0.81, 0.73, 0.77, color=BAD_DEEP, lw=2.6, rad=0.18)

# ===== GOOD pathway (left column) — 5 cascade events stacked vertically =====
# Each event row: event box at x=0.13, mini bar at x=0.34
# Spaced at y = 0.68, 0.55, 0.42, 0.29, 0.16
good_events = [
    (0.68, 'Mutation\nclearance',       '#118ab2', d_miss_g, d_miss_b, 'Δ missense', 'count'),
    (0.55, 'Neoantigen\nclone loss',    '#06aed5', d_neo_g,  d_neo_b,  'Δ MHC-I binders','count'),
    (0.42, 'SBS5 mutation\nclearance',  '#0f8b78', d_sbs_g,  d_sbs_b,  'Δ SBS5', 'count'),
    (0.29, 'Treg / MHC-II\nreprogramming','#c11456',d_treg_g, d_treg_b, 'Δ Treg z', 'z'),
    (0.16, 'BCR (IGH)\ninfiltration',   '#7a3aad', d_igh_g,  d_igh_b,  'Δ IGH count', 'count'),
]
for y_pos, lbl, col, dg, db, mlbl, munit in good_events:
    fancy_box(0.13, y_pos, 0.15, 0.08, lbl, col, textcolor='white', fontsize=9.5)
    mini_bar(0.34, y_pos, 0.14, 0.09, dg, db, mlbl, munit)

# Connecting arrows DOWN through good cascade (from event to event)
for i in range(len(good_events)-1):
    y1 = good_events[i][0] - 0.04
    y2 = good_events[i+1][0] + 0.04
    curve_arrow(0.13, y1, 0.13, y2, color=good_events[i+1][2], lw=2.0, rad=0, mutation=16)

# First arrow: branch → first event
curve_arrow(0.23, 0.74, 0.13, good_events[0][0]+0.04, color=good_events[0][2], lw=2.2, rad=-0.15)

# Good → outcome (bottom)
fancy_box(0.24, 0.05, 0.32, 0.07,
          'Complete / near-complete response (TRG 0–1)',
          GOOD_DEEP, fontsize=11)
curve_arrow(0.13, good_events[-1][0]-0.04, 0.22, 0.09, color=GOOD_DEEP, lw=2.2, rad=0.10)

# ===== POOR pathway (right column) — simplified, one big box =====
fancy_box(0.78, 0.50, 0.25, 0.42,
          'Minimal molecular change\n\n'
          f'• Δ missense  ≈  {d_miss_b:+.0f}\n'
          f'• Δ neoantigens  ≈  {d_neo_b:+.0f}\n'
          f'• Δ SBS5  ≈  {d_sbs_b:+.0f}\n'
          f'• Δ Treg z  ≈  {d_treg_b:+.1f}\n'
          f'• Δ IGH  ≈  {d_igh_b:+.0f}\n\n'
          '→ Primary treatment\n   insensitivity',
          BAD_DEEP, fontsize=10.5)
curve_arrow(0.78, 0.76, 0.78, 0.72, color=BAD_DEEP, lw=2.3, rad=0)

fancy_box(0.78, 0.14, 0.30, 0.07, 'Partial / minimal response (TRG 2–3)', BAD_DEEP, fontsize=11)
curve_arrow(0.78, 0.28, 0.78, 0.18, color=BAD_DEEP, lw=2.3, rad=0)

# Central divider
ax.plot([0.55, 0.55], [0.02, 0.78], color='#aab3bf', lw=0.8, ls='--', alpha=0.6)

save_panel(fig, 'Fig6F_cascade', OUT)

print('=== Fig 6 v3.2 targeted fixes saved ===')
