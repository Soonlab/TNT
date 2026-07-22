#!/usr/bin/env python3
"""
Supplementary external validation figures (Task 3 of GM revision).
Recreate Fig7A-C in Nature/Cell style as SuppFig_external_forest / _meta_zscore / _GSE150082_DSB.
"""
import pandas as pd, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ANALYSIS = Path('/mnt/sda1/data/TNT/analysis')
EXT = ANALYSIS/'11_external_validation'
OUT = ANALYSIS/'figures'/'supp'
OUT.mkdir(parents=True, exist_ok=True)

GOOD = '#2E86AB'
BAD  = '#E63946'
SIG_COLORS = {
    'DSB_HDR_repair':    '#118ab2',
    'E2F_MYC_cellcycle': '#2a9d8f',
    'CD8_proliferation': '#ef476f',
    'EMT':               '#e76f51',
}

def setup():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans','font.size':8,
        'axes.labelsize':8,'axes.titlesize':10,'axes.titleweight':'bold',
        'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':0.6,'axes.edgecolor':'#222',
        'xtick.labelsize':7,'ytick.labelsize':7,
        'xtick.major.width':0.6,'ytick.major.width':0.6,
        'pdf.fonttype':42,'ps.fonttype':42,
        'savefig.dpi':600,'savefig.bbox':'tight',
        'figure.facecolor':'white'
    })

setup()

sig = pd.read_csv(EXT/'signature_stats_manual.tsv', sep='\t')
meta = pd.read_csv(EXT/'meta_analysis_manual.tsv', sep='\t')
cohort = pd.read_csv(EXT/'cohort_summary_manual.tsv', sep='\t')

# ---------- Supp Fig: Forest plot (per signature) ----------
# We plot delta (good - bad) with a 1.96/sqrt(n) pseudo-SE as a visual aid.
sig = sig.copy()
sig['se'] = 1/np.sqrt(sig['n_good']) + 1/np.sqrt(sig['n_bad'])
sig['ci_lo'] = sig['delta'] - 1.96*sig['se']
sig['ci_hi'] = sig['delta'] + 1.96*sig['se']

signatures = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
fig, axes = plt.subplots(1, len(signatures), figsize=(9.5, 3.2), sharey=True)
gse_order = sorted(sig['gse'].unique())
for ax, s in zip(axes, signatures):
    sub = sig[sig.signature==s].set_index('gse').loc[gse_order].reset_index()
    ys = np.arange(len(sub))
    c = SIG_COLORS[s]
    ax.hlines(ys, sub['ci_lo'], sub['ci_hi'], color=c, lw=1.1, alpha=0.9)
    ax.scatter(sub['delta'], ys, s=24, color=c, edgecolor='white', linewidth=0.6, zorder=3)
    ax.axvline(0, color='#888', lw=0.5, ls='--')
    ax.set_yticks(ys); ax.set_yticklabels(sub['gse'])
    ax.set_title(s.replace('_',' '), loc='left')
    ax.set_xlabel('Δ score (good − bad)')
    ax.tick_params(length=2.5, width=0.6)
axes[0].invert_yaxis()
fig.suptitle('Supplementary Figure  |  External validation forest plots (7 GEO cohorts, N=290)',
             fontsize=10, fontweight='bold', y=1.03, x=0.02, ha='left')
fig.tight_layout()
fig.savefig(OUT/'SuppFig_external_forest.pdf', bbox_inches='tight')
fig.savefig(OUT/'SuppFig_external_forest.png', dpi=600, bbox_inches='tight')
plt.close(fig)

# ---------- Supp Fig: Meta Z-score barplot ----------
fig, ax = plt.subplots(figsize=(5.0, 3.0))
meta_sorted = meta.sort_values('Z_stouffer')
colors = [SIG_COLORS.get(s, '#888') for s in meta_sorted['signature']]
bars = ax.barh(meta_sorted['signature'], meta_sorted['Z_stouffer'],
               color=colors, edgecolor='white', linewidth=0.6, alpha=0.85)
ax.axvline(0, color='#888', lw=0.5, ls='--')
ax.axvline(1.96, color='#aaa', lw=0.4, ls=':')
ax.axvline(-1.96, color='#aaa', lw=0.4, ls=':')
ax.set_xlabel("Stouffer's Z (meta, 7 cohorts)")
ax.set_title('Cross-cohort meta-analysis', loc='left', fontsize=9)
ax.tick_params(length=2.5, width=0.6)
# Put P values consistently at the positive x-axis end (outside bar) to avoid
# overlap with y-tick labels for negative bars.
xmax = max(abs(meta_sorted['Z_stouffer'].min()), abs(meta_sorted['Z_stouffer'].max())) + 0.5
ax.set_xlim(-xmax, xmax)
for bar, z, p in zip(bars, meta_sorted['Z_stouffer'], meta_sorted['p_meta_onesided']):
    # Always annotate at right side of axis to avoid collision with bar body
    # and y-tick labels (bars can be either direction from 0).
    ax.text(xmax - 0.08, bar.get_y()+bar.get_height()/2,
            f'P={p:.2f}', va='center', ha='right', fontsize=6.5, color='#222')
fig.tight_layout()
fig.savefig(OUT/'SuppFig_meta_zscore.pdf', bbox_inches='tight')
fig.savefig(OUT/'SuppFig_meta_zscore.png', dpi=600, bbox_inches='tight')
plt.close(fig)

# ---------- Supp Fig: GSE150082 DSB signature good vs bad ----------
scores = pd.read_csv(EXT/'GSE150082_signature_scores.tsv', sep='\t')
# find DSB column / response column
resp_col = [c for c in scores.columns if c.lower() in ('response','response_bin','resp')]
# more permissive
if not resp_col:
    for c in scores.columns:
        if 'response' in c.lower(): resp_col=[c]; break
resp_col = resp_col[0] if resp_col else None
dsb_col = None
for c in scores.columns:
    if 'DSB' in c or 'HDR' in c or 'DNA_repair' in c or 'dsb' in c.lower():
        dsb_col = c; break
print('GSE150082 resp_col=', resp_col, ' dsb_col=', dsb_col, ' cols=', scores.columns.tolist()[:15])

fig, ax = plt.subplots(figsize=(3.2, 3.2))
if resp_col and dsb_col:
    g = pd.to_numeric(scores.loc[scores[resp_col].astype(str).str.lower().str.startswith('good'), dsb_col], errors='coerce').dropna().values
    b = pd.to_numeric(scores.loc[scores[resp_col].astype(str).str.lower().str.startswith('bad'), dsb_col], errors='coerce').dropna().values
    data=[g,b]; labels=['good','bad']; colors=[GOOD,BAD]
    bp = ax.boxplot(data, widths=0.4, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black',linewidth=1.1))
    for p,c in zip(bp['boxes'],colors): p.set_facecolor(c); p.set_alpha(0.25); p.set_edgecolor(c)
    rng=np.random.default_rng(3)
    for i,(v,c) in enumerate(zip(data,colors)):
        ax.scatter(i+1+rng.uniform(-0.08,0.08,size=len(v)), v, s=16, color=c,
                   edgecolor='white', linewidth=0.5, zorder=3)
    ax.set_xticks([1,2]); ax.set_xticklabels([f'{l}\n(n={len(v)})' for l,v in zip(labels,data)])
    ax.set_ylabel(dsb_col)
    if len(g)>=2 and len(b)>=2:
        p=stats.mannwhitneyu(g,b,alternative='two-sided').pvalue
        ax.text(0.98,0.97,f'P = {p:.3f}',transform=ax.transAxes,ha='right',va='top',fontsize=7.5)
else:
    ax.text(0.5,0.5,'GSE150082 response/DSB columns not found',ha='center',va='center',transform=ax.transAxes)
ax.set_title('GSE150082  |  DSB/HDR repair signature', loc='left')
ax.tick_params(length=2.5, width=0.6)
fig.tight_layout()
fig.savefig(OUT/'SuppFig_GSE150082_DSB.pdf', bbox_inches='tight')
fig.savefig(OUT/'SuppFig_GSE150082_DSB.png', dpi=600, bbox_inches='tight')
plt.close(fig)

print('Supp figs saved to', OUT)
