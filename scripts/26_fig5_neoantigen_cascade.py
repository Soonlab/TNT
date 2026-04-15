#!/usr/bin/env python3
"""
Fig5 — Neoantigen cascade (new main figure for Genome Medicine revision).

Four panels:
  A. Pre-treatment missense (n_nonsyn) good vs bad — box + strip
  B. Pre-treatment pVACseq strong-binder sites good vs bad — box + strip
  C. Paired Δ pVACseq binders good vs bad — slopegraph + box
  D. HLA LOH tumor allelic ratio subj 3 & 4 (paired pre -> post, HLA-B/C)

Style inspired by Cristescu Science 2018, Mariathasan Nature 2018,
Chen Nature 2021 — clean, thin spines, Arial-like sans, high contrast.
"""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ANALYSIS = Path('/mnt/sda1/data/TNT/analysis')
OUTDIR = ANALYSIS / 'figures' / 'panels'
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------- Style helper (Nature/Cell/Science aesthetic) -------------
GOOD = '#2E86AB'   # slate blue
BAD  = '#E63946'   # coral red
GREY = '#6c757d'

def setup_style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',   # Arial/Helvetica alternative present
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 10,
        'axes.titleweight': 'bold',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.6,
        'axes.edgecolor': '#222222',
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'legend.fontsize': 7,
        'legend.frameon': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'savefig.dpi': 600,
        'savefig.bbox': 'tight',
        'figure.facecolor': 'white',
    })

def panel_letter(ax, letter):
    ax.text(-0.18, 1.08, letter, transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='top', ha='left')

def box_strip(ax, df, val_col, group_col='response_bin',
              order=('good','bad'), colors=(GOOD, BAD), ylabel='',
              title=''):
    data = []
    for g in order:
        v = pd.to_numeric(df.loc[df[group_col]==g, val_col], errors='coerce').dropna().values
        data.append(v)
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.4,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.1),
                    whiskerprops=dict(color='#222', linewidth=0.6),
                    capprops=dict(color='#222', linewidth=0.6))
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c); patch.set_alpha(0.25)
        patch.set_edgecolor(c); patch.set_linewidth(0.9)
    rng = np.random.default_rng(7)
    for i, (v, c) in enumerate(zip(data, colors)):
        jx = i + rng.uniform(-0.09, 0.09, size=len(v))
        ax.scatter(jx, v, s=18, color=c, edgecolor='white',
                   linewidth=0.6, zorder=3, alpha=0.95)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f'{g}\n(n={len(v)})' for g,v in zip(order, data)])
    ax.set_ylabel(ylabel)
    if len(data[0])>=2 and len(data[1])>=2:
        p = stats.mannwhitneyu(data[0], data[1], alternative='two-sided').pvalue
    else:
        p = np.nan
    # p value top-right
    ax.text(0.98, 0.97, f'P = {p:.3f}' if not np.isnan(p) else 'P = NA',
            transform=ax.transAxes, ha='right', va='top', fontsize=7.5,
            color='#222222')
    if title:
        ax.set_title(title, loc='left')
    return p

def slope_with_box(ax, df_delta, val_col='delta_binders', group_col='response',
                   order=('good','bad'), colors=(GOOD, BAD),
                   ylabel='Δ pVACseq binders (post − pre)', title=''):
    # Slopegraph: show zero vs delta for each subject; we instead show subject-level bars?
    # User requested slopegraph + box for paired delta — use subject-level jitter + horizontal box
    rng = np.random.default_rng(11)
    xs, ys, cs = [], [], []
    data_by_g = {}
    for i, g in enumerate(order):
        v = pd.to_numeric(df_delta.loc[df_delta[group_col]==g, val_col], errors='coerce').dropna().values
        data_by_g[g] = v
    # Draw connecting "zero-to-delta" slopes on the left (illustrates direction)
    for i, g in enumerate(order):
        vals = data_by_g[g]; c = colors[i]
        for v in vals:
            ax.plot([i-0.18, i+0.18], [0, v], color=c, alpha=0.45, lw=0.7, zorder=1)
        # endpoints
        ax.scatter([i-0.18]*len(vals), [0]*len(vals), color=GREY, s=14,
                   edgecolor='white', linewidth=0.5, zorder=2)
        ax.scatter([i+0.18]*len(vals), vals, color=c, s=22,
                   edgecolor='white', linewidth=0.6, zorder=3)
    # Overlay narrow boxplot of deltas
    box_data = [data_by_g[g] for g in order]
    bp = ax.boxplot(box_data, positions=[i+0.55 for i in range(len(order))],
                    widths=0.18, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.0),
                    whiskerprops=dict(color='#222', linewidth=0.5),
                    capprops=dict(color='#222', linewidth=0.5))
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c); patch.set_alpha(0.25); patch.set_edgecolor(c)
    ax.axhline(0, color='#888', lw=0.5, ls='--', zorder=0)
    ax.set_xticks([i+0.15 for i in range(len(order))])
    ax.set_xticklabels([f'{g}\n(n={len(data_by_g[g])})' for g in order])
    ax.set_ylabel(ylabel)
    # stats
    p = stats.mannwhitneyu(data_by_g['good'], data_by_g['bad'],
                           alternative='two-sided').pvalue
    ax.text(0.98, 0.97, f'P = {p:.3f}', transform=ax.transAxes,
            ha='right', va='top', fontsize=7.5)
    if title: ax.set_title(title, loc='left')
    return p

def loh_panel(ax, loh):
    """Panel D: paired HLA LOH tumor_ratio pre->post subj 3 & 4 per locus."""
    subjects = [3, 4]
    loci = ['HLA-A','HLA-B','HLA-C']
    sub_x = {3: 0, 4: 1}
    offsets = {'HLA-A': -0.22, 'HLA-B': 0.0, 'HLA-C': 0.22}
    marker = {'HLA-A': 'o', 'HLA-B': 's', 'HLA-C': '^'}
    for s in subjects:
        for L in loci:
            pre  = loh[(loh.subject_id==s)&(loh.locus==L)&
                       (loh['sample'].str.endswith('-PR'))]
            post = loh[(loh.subject_id==s)&(loh.locus==L)&
                       (loh['sample'].str.endswith('-PO'))]
            if len(pre)==0 or len(post)==0: continue
            # use allelic imbalance |ratio - 0.5| (higher = more LOH-like)
            y_pre  = abs(float(pre['tumor_ratio'].iloc[0]) - 0.5)
            y_post = abs(float(post['tumor_ratio'].iloc[0]) - 0.5)
            x0 = sub_x[s] + offsets[L]
            ax.plot([x0-0.06, x0+0.06], [y_pre, y_post],
                    color='#888', lw=0.8, alpha=0.6, zorder=1)
            ax.scatter([x0-0.06], [y_pre],  marker=marker[L], s=34,
                       color='#1d3557', edgecolor='white', linewidth=0.6, zorder=3)
            ax.scatter([x0+0.06], [y_post], marker=marker[L], s=34,
                       color='#2a9d8f', edgecolor='white', linewidth=0.6, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Subject 3\n(good)', 'Subject 4\n(good)'])
    ax.set_ylabel('Allelic imbalance  |tumor ratio − 0.5|')
    ax.set_title('HLA LOH clone clearance (pre → post)', loc='left')
    # legend
    from matplotlib.lines import Line2D
    leg = [
        Line2D([0],[0], marker='o', color='w', label='HLA-A',
               markerfacecolor='#1d3557', markersize=6),
        Line2D([0],[0], marker='s', color='w', label='HLA-B',
               markerfacecolor='#1d3557', markersize=6),
        Line2D([0],[0], marker='^', color='w', label='HLA-C',
               markerfacecolor='#1d3557', markersize=6),
        Line2D([0],[0], marker='o', color='w', label='pre',
               markerfacecolor='#1d3557', markersize=6),
        Line2D([0],[0], marker='o', color='w', label='post',
               markerfacecolor='#2a9d8f', markersize=6),
    ]
    ax.legend(handles=leg, loc='upper right', fontsize=6.5, ncol=2)
    ax.set_ylim(-0.01, ax.get_ylim()[1]*1.05)

# ------------- Load data -------------
neo  = pd.read_csv(ANALYSIS/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
dneo = pd.read_csv(ANALYSIS/'03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')
tmb  = pd.read_csv(ANALYSIS/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
loh  = pd.read_csv(ANALYSIS/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')

# Pre-treatment subsets
pre_tmb = tmb[tmb['timepoint']=='pre'].copy()
pre_neo = neo[neo['timepoint']=='pre'].copy()

# ------------- Figure -------------
setup_style()
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
for ax in axes.ravel():
    ax.tick_params(length=2.5, width=0.6)

# A
p_A = box_strip(axes[0,0], pre_tmb, 'n_nonsyn',
                ylabel='Pre-treatment missense mutations',
                title='Pre missense burden')
panel_letter(axes[0,0], 'A')

# B
p_B = box_strip(axes[0,1], pre_neo, 'n_sites_with_strong',
                ylabel='Pre pVACseq strong-binder sites (< 50 nM)',
                title='Pre strong neoantigen sites')
panel_letter(axes[0,1], 'B')

# C
p_C = slope_with_box(axes[1,0], dneo, val_col='delta_binders',
                     ylabel='Δ binders (post − pre)',
                     title='Paired Δ neoantigen binders')
panel_letter(axes[1,0], 'C')

# D
loh_panel(axes[1,1], loh)
panel_letter(axes[1,1], 'D')

# subtitle / n summary
fig.suptitle('Figure 5  |  Neoantigen cascade in TNT response',
             fontsize=10.5, fontweight='bold', y=1.01, x=0.06, ha='left')

fig.tight_layout(rect=[0,0,1,0.98])
fig.savefig(OUTDIR/'Fig5_neoantigen_cascade.pdf', bbox_inches='tight')
fig.savefig(OUTDIR/'Fig5_neoantigen_cascade.png', dpi=600, bbox_inches='tight')
plt.close(fig)

# Summary stats file
with open(OUTDIR/'Fig5_neoantigen_cascade_stats.txt','w') as f:
    f.write('Fig5 panel statistics (Mann-Whitney U two-sided)\n')
    f.write(f'A pre n_nonsyn good vs bad:              P = {p_A:.4f}\n')
    f.write(f'B pre n_sites_with_strong good vs bad:   P = {p_B:.4f}\n')
    f.write(f'C delta_binders good vs bad:             P = {p_C:.4f}\n')
    f.write('D HLA LOH clone clearance subj 3 & 4: descriptive (paired PR->PO)\n')
print('Fig5 written to', OUTDIR/'Fig5_neoantigen_cascade.pdf')
print(f'P values: A={p_A:.4f}  B={p_B:.4f}  C={p_C:.4f}')
