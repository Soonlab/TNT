"""Forest plot for corrected external meta-analysis (CD8-cytotoxic axis)."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
FIG = Path('/mnt/sda1/data/TNT/analysis/figures/supp'); FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(OUT/'v3_signature_response_stats.tsv', sep='\t')
meta = pd.read_csv(OUT/'v3_meta_overall.tsv', sep='\t')

SIGS_ORDER = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration',
              'Tumor_cellcycle','DSB_HDR_repair','E2F_MYC_cellcycle','EMT']

fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), gridspec_kw={'width_ratios':[1.0,0.9,0.9]})

# Panel 1: CD8_cytotoxic forest (per-cohort)
ax = axes[0]
sub = df[df.signature=='CD8_cytotoxic'].copy().sort_values('delta')
y = np.arange(len(sub))
# approximate 95% CI from two-sided p and delta
ci_halfwidth = []
for _,r in sub.iterrows():
    z = stats.norm.isf(max(r.pvalue,1e-300)/2)
    se = abs(r.delta)/max(z,0.01)
    ci_halfwidth.append(1.96*se)
ci_halfwidth = np.array(ci_halfwidth)
ax.errorbar(sub.delta, y, xerr=ci_halfwidth, fmt='s', color='#2E86AB', capsize=3, ms=7, lw=1.5)
for yi, (_,r) in zip(y, sub.iterrows()):
    ax.text(r.delta, yi+0.25, f'n={r.n_good+r.n_bad}', ha='center', va='bottom', fontsize=8, color='gray')
ax.axvline(0, color='k', lw=0.5)
mrow = meta[meta.signature=='CD8_cytotoxic'].iloc[0]
# diamond for meta
meta_delta = sub.delta.mean()
ax.plot([meta_delta-0.05, meta_delta, meta_delta+0.05, meta_delta], [-0.7, -0.5, -0.7, -0.9], 'D-',
        color='#E63946', ms=0, lw=1.5)
ax.fill([meta_delta-0.05, meta_delta, meta_delta+0.05, meta_delta],
        [-0.7, -0.5, -0.7, -0.9], color='#E63946', alpha=0.8)
ax.text(meta_delta, -1.3, f'Meta Z={mrow.Z:+.2f}, p={mrow.p_meta:.3f}',
        ha='center', fontsize=10, color='#E63946', fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(sub.gse)
ax.set_ylim(-1.8, len(sub)+0.5)
ax.set_xlabel('Δ signature score (good − bad)')
ax.set_title('CD8_cytotoxic (pure CD8 effector markers)\n9 cohorts, N=721', fontsize=11)
ax.spines[['top','right']].set_visible(False)

# Panel 2: meta Z across all signatures (bar)
ax = axes[1]
order_meta = meta.set_index('signature').loc[[s for s in SIGS_ORDER if s in meta.signature.values]]
colors = ['#2E86AB' if z>0 else '#E63946' for z in order_meta.Z]
ax.barh(range(len(order_meta)), order_meta.Z, color=colors, edgecolor='k', lw=0.5)
ax.axvline(0, color='k', lw=0.5)
ax.axvline(1.96, color='gray', ls='--', lw=0.5)
ax.axvline(-1.96, color='gray', ls='--', lw=0.5)
ax.set_yticks(range(len(order_meta)))
ax.set_yticklabels(order_meta.index)
ax.set_xlabel('Meta Z-score (Stouffer, √N weighted)')
ax.set_title('Signature-level meta-analysis\n(9 nCRT cohorts)', fontsize=11)
for i,(z,p) in enumerate(zip(order_meta.Z, order_meta.p_meta)):
    ax.text(z + (0.1 if z>0 else -0.1), i, f'p={p:.3f}',
            va='center', ha='left' if z>0 else 'right', fontsize=8)
ax.invert_yaxis()
ax.spines[['top','right']].set_visible(False)

# Panel 3: CD8 vs Tumor_cellcycle per-cohort scatter — showing decoupling
ax = axes[2]
cd8 = df[df.signature=='CD8_cytotoxic'].set_index('gse')['delta']
prol = df[df.signature=='Tumor_cellcycle'].set_index('gse')['delta']
shared = cd8.index.intersection(prol.index)
ax.scatter(prol.loc[shared], cd8.loc[shared], s=70, color='#2E86AB', edgecolor='k', lw=0.7, zorder=3)
for g in shared:
    ax.annotate(g.replace('GSE',''), (prol.loc[g], cd8.loc[g]),
                xytext=(5,5), textcoords='offset points', fontsize=8)
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel('Δ Tumor_cellcycle (good − bad)')
ax.set_ylabel('Δ CD8_cytotoxic (good − bad)')
ax.set_title('CD8 effector vs tumor proliferation\ndecouple across cohorts', fontsize=11)
ax.spines[['top','right']].set_visible(False)

plt.tight_layout()
plt.savefig(FIG/'SuppFig_v3_CD8_meta_forest.png', dpi=600, bbox_inches='tight')
plt.savefig(FIG/'SuppFig_v3_CD8_meta_forest.pdf', bbox_inches='tight')
print('saved', FIG/'SuppFig_v3_CD8_meta_forest.png')
print('\nFinal meta-analysis:')
print(meta.to_string(index=False))
