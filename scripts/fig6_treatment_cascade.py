"""
Figure 6 — Treatment-induced response cascade (CORE FIGURE)
  6A: Swimmer plot - per-subject timeline of all key changes
  6B: Slope graph - paired pre→post for key features
  6C: Waterfall - delta missense/neoantigen/IGH per subject sorted
  6D: Multi-axis paired delta heatmap (subjects × features)
  6E: Cascade flow diagram (good responder)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
trust = pd.read_csv(ROOT/'06_rna_immune/trust4_summary.tsv', sep='\t')
neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
neo_delta = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')
delta_22 = pd.read_csv(ROOT/'09_integration/paired_delta/delta_22sigs_response.tsv', sep='\t')
ssg = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce')
sbs = pd.read_csv(ROOT/'01_wes_signatures/sbs_summary_key.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')

UNMATCHED = [13,15,16,17,18,19,33]
PAIRED_SUBJ = list(range(1, 15))  # subj 1-14 with WES paired (mostly)

# Build per-subject delta master
def build_delta(subjects):
    rows = []
    for s in subjects:
        # WES TMB
        ts = tmb[(tmb.subject_id==s) & tmb.timepoint.isin(['pre','post'])]
        wes_delta_tmb = wes_delta_miss = np.nan
        if {'pre','post'}.issubset(set(ts.timepoint)):
            pre_t = ts[ts.timepoint=='pre'].iloc[0]; post_t = ts[ts.timepoint=='post'].iloc[0]
            wes_delta_tmb = post_t.TMB_nonsyn_per_Mb - pre_t.TMB_nonsyn_per_Mb
            wes_delta_miss = post_t.n_nonsyn - pre_t.n_nonsyn
        # SBS5
        ss = sbs[(sbs.subject_id==s) & sbs.timepoint.isin(['pre','post'])]
        sbs_delta = np.nan
        if {'pre','post'}.issubset(set(ss.timepoint)):
            sbs_delta = ss[ss.timepoint=='post'].iloc[0].SBS5 - ss[ss.timepoint=='pre'].iloc[0].SBS5
        # Neoantigen delta
        nd_row = neo_delta[neo_delta.subject_id==s]
        neo_dlt_b = nd_row.delta_binders.iloc[0] if len(nd_row) else np.nan
        # HLA LOH delta
        loh_pre = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PR')]
        loh_post = loh[(loh.subject_id==s) & loh['sample'].str.endswith('-PO')]
        loh_delta = (loh_post.LOH_call.sum() - loh_pre.LOH_call.sum()) if len(loh_pre)+len(loh_post)>0 else np.nan
        # RNA signatures (Treg, MHC_II, IGH from TRUST4)
        rs = rna_inv[(rna_inv.subject_id==s) & rna_inv.timepoint.isin(['pre','post'])]
        sig_treg = sig_mhc2 = sig_cd8ex = np.nan
        igh_n = np.nan
        if {'pre','post'}.issubset(set(rs.timepoint)):
            pre_id = rs[rs.timepoint=='pre'].sample_id.iloc[0]
            post_id = rs[rs.timepoint=='post'].sample_id.iloc[0]
            try:
                sig_treg = sigs.loc[post_id,'Treg'] - sigs.loc[pre_id,'Treg']
                sig_mhc2 = sigs.loc[post_id,'MHC_II'] - sigs.loc[pre_id,'MHC_II']
                sig_cd8ex = sigs.loc[post_id,'CD8_exhaustion'] - sigs.loc[pre_id,'CD8_exhaustion']
            except: pass
            try:
                tr_pre = trust[trust.sample_id==pre_id]
                tr_post = trust[trust.sample_id==post_id]
                igh_n = tr_post.IGH_n.iloc[0] - tr_pre.IGH_n.iloc[0]
            except: pass
        rows.append({'subject_id':s, 'response': clin[clin.subject_id==s].response_bin.iloc[0],
                     'd_TMB':wes_delta_tmb, 'd_missense':wes_delta_miss, 'd_SBS5':sbs_delta,
                     'd_neo_binders':neo_dlt_b, 'd_HLA_LOH':loh_delta,
                     'd_Treg':sig_treg, 'd_MHC_II':sig_mhc2, 'd_CD8_exh':sig_cd8ex,
                     'd_IGH_n':igh_n})
    return pd.DataFrame(rows)

D = build_delta(PAIRED_SUBJ)

# ===========================================================
# 6A — Swimmer plot
# ===========================================================
fig, ax = plt.subplots(figsize=(13, 6))
# Sort: good first, by subject
D_sw = D.sort_values(['response','subject_id']).reset_index(drop=True)
N = len(D_sw)
y_pos = np.arange(N)

# Background bars by response (long horizontal track)
for i, (_, r) in enumerate(D_sw.iterrows()):
    color = PAL_RESP[r.response]
    ax.barh(i, 1.0, height=0.65, color=color, alpha=0.18, edgecolor=color, linewidth=0.5)
# Dots for events: x positions for each metric
metric_x = {'d_TMB':0.10, 'd_missense':0.20, 'd_SBS5':0.30, 'd_neo_binders':0.42,
            'd_HLA_LOH':0.54, 'd_Treg':0.66, 'd_MHC_II':0.74, 'd_CD8_exh':0.82, 'd_IGH_n':0.92}
metric_label = {'d_TMB':'TMB','d_missense':'Missense','d_SBS5':'SBS5','d_neo_binders':'Neoantigen',
               'd_HLA_LOH':'HLA LOH','d_Treg':'Treg','d_MHC_II':'MHC II','d_CD8_exh':'CD8 exh','d_IGH_n':'IGH (BCR)'}

for i, (_, r) in enumerate(D_sw.iterrows()):
    for m, x in metric_x.items():
        v = r[m]
        if pd.isna(v): continue
        # Direction: for d_TMB/missense/SBS5/neo/HLA_LOH negative = good (clearance)
        # For d_Treg/MHC_II/IGH/CD8_exh positive = good (immune activation)
        is_clearance = m in ['d_TMB','d_missense','d_SBS5','d_neo_binders','d_HLA_LOH']
        if is_clearance:
            magnitude = abs(v); color = '#118ab2' if v<0 else '#d62828'  # blue=cleared, red=gained
        else:
            magnitude = abs(v); color = '#0f8b78' if v>0 else '#c1272d'  # green=activated, red=lost
        size = min(800, 50 + magnitude*8 if m=='d_neo_binders' else 50 + magnitude*1.5 if m=='d_IGH_n' else 50 + abs(v)*100)
        ax.scatter(x, i, s=size, c=color, alpha=0.85, edgecolor='white', linewidth=1, zorder=5)

# Top header: metric labels
for m, x in metric_x.items():
    ax.text(x, N-0.2, metric_label[m], ha='center', va='bottom', fontsize=8.5,
            color='#1d3557', fontweight='bold', rotation=45)

# Y-axis: subject labels
ax.set_yticks(y_pos)
ax.set_yticklabels([f"S{int(r.subject_id)}  ({r.response})" for _, r in D_sw.iterrows()], fontsize=8.5)
ax.set_ylim(-0.5, N+1.5)
ax.set_xlim(0, 1)
ax.set_xticks([])

# Grouping bracket on left
good_n = (D_sw.response=='good').sum()
ax.axhline(good_n-0.5, color='#1d3557', lw=1.2)
ax.text(-0.05, (good_n-1)/2, 'Good\nresponders', ha='right', va='center', fontsize=10,
        fontweight='bold', color=GOOD)
ax.text(-0.05, good_n + (N-good_n-1)/2, 'Poor\nresponders', ha='right', va='center', fontsize=10,
        fontweight='bold', color=BAD)

# Legend
leg = [
    plt.scatter([],[], s=120, c='#118ab2', alpha=0.85, edgecolor='white', lw=1, label='Tumor clearance (Δ < 0)'),
    plt.scatter([],[], s=120, c='#d62828', alpha=0.85, edgecolor='white', lw=1, label='Tumor gain (Δ > 0)'),
    plt.scatter([],[], s=120, c='#0f8b78', alpha=0.85, edgecolor='white', lw=1, label='Immune activation (Δ > 0)'),
    plt.scatter([],[], s=120, c='#c1272d', alpha=0.85, edgecolor='white', lw=1, label='Immune loss (Δ < 0)'),
]
ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=4, fontsize=8, frameon=False)

ax.set_title('Subject-level treatment-induced changes (post − pre)\nbubble size ∝ |Δ|',
             fontsize=12, fontweight='bold', color='#1d3557', pad=20)
for s in ['top','right','bottom','left']: ax.spines[s].set_visible(False)
save_panel(fig, 'Fig6A_swimmer', OUT)

# ===========================================================
# 6B — Slope graph (paired pre → post)
# ===========================================================
fig, axes = plt.subplots(1, 4, figsize=(15, 4.2))

def slope_panel(ax, df_pre, df_post, ylabel, title, transform=lambda x: x):
    pairs = []
    for s in PAIRED_SUBJ:
        if s in df_pre.index and s in df_post.index:
            v_pre = transform(df_pre.loc[s]); v_post = transform(df_post.loc[s])
            resp = clin[clin.subject_id==s].response_bin.iloc[0]
            pairs.append((s, resp, v_pre, v_post))
    for s, resp, vp, vpo in pairs:
        color = PAL_RESP[resp]
        ax.plot([0, 1], [vp, vpo], color=color, lw=1.8, alpha=0.7, marker='o', markersize=8,
                markeredgecolor='white', markeredgewidth=0.8, zorder=3)
        # Label end point if extreme
        if abs(vpo - vp) > 0:
            ax.annotate(f'S{s}', xy=(1.04, vpo), va='center', fontsize=7, color=color)
    ax.set_xticks([0,1]); ax.set_xticklabels(['pre','post'])
    ax.set_xlim(-0.15, 1.4); ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10.5)
    add_axis_spines(ax)

# Build per-subject pre/post metrics
def get_metric(df, sample_ids, col):
    return pd.Series({s: df.loc[df.sample_id==sid, col].iloc[0] if (df.sample_id==sid).any() else np.nan for s, sid in sample_ids.items()})

# WES sample IDs
wes_pre_ids = {s: f'{s}-PR' for s in PAIRED_SUBJ}
wes_post_ids = {s: f'{s}-PO' for s in PAIRED_SUBJ}
rna_pre_ids = {s: rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='pre')].sample_id.iloc[0]
               if len(rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='pre')]) else None for s in PAIRED_SUBJ}
rna_post_ids = {s: rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='post')].sample_id.iloc[0]
                if len(rna_inv[(rna_inv.subject_id==s) & (rna_inv.timepoint=='post')]) else None for s in PAIRED_SUBJ}

# Panel 1: TMB
tmb_pre = pd.Series({s: tmb[tmb.sample_id==f'{s}-PR'].TMB_nonsyn_per_Mb.iloc[0]
                     if (tmb.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ})
tmb_post = pd.Series({s: tmb[tmb.sample_id==f'{s}-PO'].TMB_nonsyn_per_Mb.iloc[0]
                      if (tmb.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ})
slope_panel(axes[0], tmb_pre, tmb_post, 'TMB (/Mb)', 'TMB pre → post')

# Panel 2: Neoantigen binders
neo_pre = pd.Series({s: neo[neo.sample_id==f'{s}-PR'].n_binders_500nM.iloc[0]
                     if (neo.sample_id==f'{s}-PR').any() else np.nan for s in PAIRED_SUBJ})
neo_post = pd.Series({s: neo[neo.sample_id==f'{s}-PO'].n_binders_500nM.iloc[0]
                      if (neo.sample_id==f'{s}-PO').any() else np.nan for s in PAIRED_SUBJ})
slope_panel(axes[1], neo_pre, neo_post, 'MHC-I binders (<500 nM)', 'Neoantigens pre → post')

# Panel 3: IGH n
def trust_get(s, sid, col):
    if sid is None: return np.nan
    sub = trust[trust.sample_id==sid]
    return sub[col].iloc[0] if len(sub) else np.nan
igh_pre = pd.Series({s: trust_get(s, rna_pre_ids[s], 'IGH_n') for s in PAIRED_SUBJ})
igh_post = pd.Series({s: trust_get(s, rna_post_ids[s], 'IGH_n') for s in PAIRED_SUBJ})
slope_panel(axes[2], igh_pre, igh_post, 'IGH clonotype count', 'BCR (IGH) pre → post')

# Panel 4: Treg
def sig_get(sid, col):
    if sid is None or sid not in sigs.index: return np.nan
    return sigs.loc[sid, col]
treg_pre = pd.Series({s: sig_get(rna_pre_ids[s], 'Treg') for s in PAIRED_SUBJ})
treg_post = pd.Series({s: sig_get(rna_post_ids[s], 'Treg') for s in PAIRED_SUBJ})
slope_panel(axes[3], treg_pre, treg_post, 'Treg signature (z)', 'Treg pre → post')

# Legend on outer side
fig.legend(handles=[mpatches.Patch(color=GOOD, label='Good responder'),
                    mpatches.Patch(color=BAD, label='Poor responder')],
           loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=2, fontsize=10, frameon=False)
fig.suptitle('Paired pre → post trajectories per subject (each line = one patient)',
             fontsize=12, fontweight='bold', y=1.07, color='#1d3557')
fig.tight_layout()
save_panel(fig, 'Fig6B_slope', OUT)

# ===========================================================
# 6C — Waterfall of delta missense, delta neoantigen, delta IGH per subject
# ===========================================================
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

# Compute deltas
df_wf = pd.DataFrame({'subject_id': PAIRED_SUBJ})
df_wf['response'] = df_wf.subject_id.map(lambda s: clin[clin.subject_id==s].response_bin.iloc[0])
df_wf['d_missense'] = df_wf.subject_id.map(lambda s:
    (tmb[tmb.sample_id==f'{s}-PO'].n_nonsyn.iloc[0] - tmb[tmb.sample_id==f'{s}-PR'].n_nonsyn.iloc[0])
    if (tmb.sample_id==f'{s}-PR').any() and (tmb.sample_id==f'{s}-PO').any() else np.nan)
df_wf['d_neo'] = df_wf.subject_id.map(lambda s:
    neo_delta[neo_delta.subject_id==s].delta_binders.iloc[0] if (neo_delta.subject_id==s).any() else np.nan)
df_wf['d_IGH'] = df_wf.subject_id.map(lambda s:
    trust_get(s, rna_post_ids[s], 'IGH_n') - trust_get(s, rna_pre_ids[s], 'IGH_n')
    if rna_pre_ids[s] and rna_post_ids[s] else np.nan)

# Sort each by value
for ax, col, title, ylabel in [
    (axes[0], 'd_missense', 'Δ missense mutations (post − pre)', 'Δ missense'),
    (axes[1], 'd_neo', 'Δ MHC-I neoantigen binders', 'Δ binders (<500 nM)'),
    (axes[2], 'd_IGH', 'Δ B-cell receptor (IGH) clonotype count', 'Δ IGH clonotypes'),
]:
    df_sorted = df_wf.dropna(subset=[col]).sort_values(col).reset_index(drop=True)
    colors = [PAL_RESP[r] for r in df_sorted.response]
    bars = ax.bar(range(len(df_sorted)), df_sorted[col], color=colors, edgecolor='white', linewidth=1.2)
    ax.axhline(0, color='#1d3557', lw=0.9)
    ax.set_xticks(range(len(df_sorted)))
    ax.set_xticklabels([f'S{int(s)}' for s in df_sorted.subject_id], fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    add_axis_spines(ax)
    # Annotate extremes
    for i, v in enumerate(df_sorted[col]):
        ax.text(i, v + (10 if v>0 else -25), f'{int(v)}', ha='center', va='bottom' if v>0 else 'top', fontsize=7, color='#1d3557')

axes[-1].set_xlabel('Subject')
fig.suptitle('Waterfall: per-subject treatment-induced change\n(green = good responder, coral = poor responder)',
             fontsize=12, fontweight='bold', y=1.0)
fig.tight_layout()
save_panel(fig, 'Fig6C_waterfall', OUT)

# ===========================================================
# 6D — Multi-axis paired delta heatmap
# ===========================================================
fig, ax = plt.subplots(figsize=(11, 6))
metrics_to_show = ['d_TMB','d_missense','d_SBS5','d_neo_binders','d_HLA_LOH',
                   'd_Treg','d_MHC_II','d_CD8_exh','d_IGH_n']
metric_labels = ['Δ TMB','Δ missense','Δ SBS5','Δ neoantigens','Δ HLA LOH',
                 'Δ Treg','Δ MHC II','Δ CD8 exhaust','Δ IGH (BCR)']
D_sorted = D.sort_values(['response','subject_id']).reset_index(drop=True)
mat = D_sorted[metrics_to_show].copy()
# Z-score per column for visual scaling
for c in metrics_to_show:
    m, s = mat[c].mean(), mat[c].std()
    mat[c] = (mat[c] - m) / (s if s else 1)

im = ax.imshow(mat.values.T, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5, interpolation='nearest')
ax.set_yticks(range(len(metric_labels))); ax.set_yticklabels(metric_labels, fontsize=9.5)
ax.set_xticks(range(len(D_sorted)))
ax.set_xticklabels([f"S{int(r.subject_id)}\n{r.response[0].upper()}" for _, r in D_sorted.iterrows()], fontsize=8)

# Response color bar above
for j, (_, r) in enumerate(D_sorted.iterrows()):
    ax.add_patch(Rectangle((j-0.5, -1.5), 1, 0.6, color=PAL_RESP[r.response], clip_on=False))
ax.text(-0.7, -1.2, 'Response →', ha='right', va='center', fontsize=9)

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label('z-score Δ (post − pre)', fontsize=9)
ax.set_xlabel('Subject (good first, then poor)')
ax.set_title('Treatment-induced multi-axis change (z-scored Δ per subject)\nBlue → tumor clearance / immune activation; Red → opposite',
             fontsize=11, fontweight='bold')
for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
save_panel(fig, 'Fig6D_delta_heatmap', OUT)

# ===========================================================
# 6E — Cascade flow diagram (good responder schematic)
# ===========================================================
fig, ax = plt.subplots(figsize=(13, 6))
ax.axis('off')

def cascade_box(ax, x, y, w, h, label, color, textcolor='white'):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                                facecolor=color, edgecolor='#1d3557', linewidth=1.5))
    ax.text(x, y, label, ha='center', va='center', fontsize=10.5, color=textcolor, fontweight='bold')

def cascade_arrow(ax, x1, y1, x2, y2, label='', color='#1d3557'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=2.2, color=color, mutation_scale=20))
    if label:
        ax.text((x1+x2)/2, (y1+y2)/2 + 0.025, label, ha='center', fontsize=9, color=color, style='italic')

# Top row: pre-treatment state
cascade_box(ax, 0.5, 0.92, 0.4, 0.10,
            'Pre-treatment good responder:\nproliferative · DNA-repair-proficient · ±HLA-LOH · neoantigen-presenting',
            '#264653')

# TNT arrow down
ax.annotate('', xy=(0.5, 0.78), xytext=(0.5, 0.85),
            arrowprops=dict(arrowstyle='->', lw=3, color=HIGHLIGHT, mutation_scale=25))
ax.text(0.55, 0.815, 'TNT (FOLFOX/CAPOX + CRT)', fontsize=10, color='#1d3557',
        fontweight='bold', va='center')

# Middle row: cascade events
events = [
    (0.10, 0.62, 'Mutation\nclearance\nΔmissense −67', '#118ab2'),
    (0.30, 0.62, 'Neoantigen\nclone elimination\nΔbinders −312', '#06aed5'),
    (0.50, 0.62, 'HLA-LOH\nclone contraction\n(S3, S4)', '#0f8b78'),
    (0.70, 0.62, 'Immune\nreprogramming\nTreg+MHC-II↑', '#ef476f'),
    (0.90, 0.62, 'Lymphocyte\ninfiltration\nIGH +1,424', '#ad2831'),
]
for x, y, lbl, col in events:
    cascade_box(ax, x, y, 0.16, 0.16, lbl, col)
    # connect to TNT arrow above
    cascade_arrow(ax, 0.5, 0.78, x, 0.70)

# Bottom: outcome
cascade_box(ax, 0.5, 0.30, 0.5, 0.12,
            'Pathologic complete / near-complete response (TRG 0–1)',
            GOOD, 'white')
for x, _, _, col in events:
    ax.annotate('', xy=(0.5, 0.36), xytext=(x, 0.54),
                arrowprops=dict(arrowstyle='->', lw=1.4, color=col, alpha=0.6, mutation_scale=15))

# Side comparison: bad responder
cascade_box(ax, 0.5, 0.08, 0.6, 0.10,
            'Poor responder: minimal molecular or immune change → primary treatment insensitivity',
            BAD, 'white')

ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title('Treatment-induced response cascade in good responders',
             fontsize=13, fontweight='bold', color='#1d3557', pad=8)
save_panel(fig, 'Fig6E_cascade_diagram', OUT)

print('\n=== Fig 6 (5 panels) complete ===')
