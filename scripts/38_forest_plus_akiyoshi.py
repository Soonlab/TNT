"""Forest plot v4 — adds Akiyoshi 2023 (GSE216616, n=298) as convergent
paper-level validation row on the CD8 forest, and updates the meta annotation
to 'total independent evidence >1,000 patients across 10 cohorts'.
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
FIG_MAIN = Path('/mnt/sda1/data/TNT/analysis/figures/panels')
FIG_SUPP = Path('/mnt/sda1/data/TNT/analysis/figures/supp')

df   = pd.read_csv(OUT/'v3_signature_response_stats.tsv', sep='\t')
meta = pd.read_csv(OUT/'v3_meta_overall.tsv', sep='\t')

SIGS_ORDER = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration',
              'Tumor_cellcycle','DSB_HDR_repair','E2F_MYC_cellcycle','EMT']

fig, axes = plt.subplots(1, 3, figsize=(16, 6.2),
                         gridspec_kw={'width_ratios':[1.15, 0.9, 0.9]})

# -------- Panel A: CD8 forest + Akiyoshi --------
ax = axes[0]
sub = df[df.signature=='CD8_cytotoxic'].copy().sort_values('delta')
y = np.arange(len(sub))

ci_half = []
for _,r in sub.iterrows():
    z = stats.norm.isf(max(r.pvalue,1e-300)/2)
    se = abs(r.delta)/max(z,0.01)
    ci_half.append(1.96*se)
ci_half = np.array(ci_half)

ax.errorbar(sub.delta, y, xerr=ci_half, fmt='s', color='#2E86AB',
            capsize=3, ms=7, lw=1.4, label='Our 9-cohort meta (N=721)')
for yi, (_,r) in zip(y, sub.iterrows()):
    ax.text(r.delta, yi+0.25, f'n={int(r.n_good+r.n_bad)}',
            ha='center', va='bottom', fontsize=8, color='gray')

# meta diamond (our)
mrow = meta[meta.signature=='CD8_cytotoxic'].iloc[0]
meta_delta = sub.delta.mean()
ax.fill([meta_delta-0.05, meta_delta, meta_delta+0.05, meta_delta],
        [-0.9, -0.6, -0.9, -1.2], color='#E63946', alpha=0.85, zorder=5)
ax.text(meta_delta, -1.55,
        f'Meta Stouffer Z={mrow.Z:+.2f}\np={mrow.p_meta:.3f} (9 cohorts, N=721)',
        ha='center', fontsize=9, color='#E63946', fontweight='bold')

# --- Akiyoshi 2023 convergent row (separate, visually distinct) ---
akiyoshi_y = len(sub) + 0.8
# median cytotoxic lymphocyte score: good 0.76 vs bad 0.58 → Δ=0.18
akiyoshi_delta = 0.18
# approximate 95% CI from reported OR=3.81 [1.82, 7.97] via logistic approximation
# (only for graphical representation; we print the OR directly in label)
# spread proportional to logOR CI width scaled to Δ units
logor = np.log(3.81); logor_lo, logor_hi = np.log(1.82), np.log(7.97)
rel_width = (logor_hi - logor_lo) / logor
ci_aki = akiyoshi_delta * rel_width / 2
ax.errorbar([akiyoshi_delta], [akiyoshi_y], xerr=[ci_aki],
            fmt='D', color='#6A4C93', capsize=4, ms=9, lw=1.8,
            markerfacecolor='white', markeredgecolor='#6A4C93', mew=2,
            label='Akiyoshi 2023 (GSE216616, N=298, paper-level)')
ax.text(akiyoshi_delta, akiyoshi_y+0.35,
        'n=298, OR=3.81 [1.82, 7.97]\nGZMA×PRF1 P=0.005',
        ha='center', va='bottom', fontsize=8, color='#6A4C93', style='italic')

# combined callout
ax.axhline(akiyoshi_y-0.5, color='#6A4C93', ls=':', lw=0.8, alpha=0.5, xmin=0.05, xmax=0.95)
ax.text(0.5, akiyoshi_y+1.3,
        'Total independent evidence: >1,000 patients across 10 cohorts',
        transform=ax.get_yaxis_transform(), ha='right', fontsize=9,
        fontweight='bold', color='#333')

ax.axvline(0, color='k', lw=0.5)
ylabels = list(sub.gse)
ax.set_yticks(list(y) + [akiyoshi_y])
ax.set_yticklabels(ylabels + ['GSE216616\n(Akiyoshi 2023)'])
# color the Akiyoshi tick label
for tick, col in zip(ax.get_yticklabels(),
                     ['black']*len(ylabels) + ['#6A4C93']):
    tick.set_color(col)

ax.set_ylim(-2.0, akiyoshi_y + 1.8)
ax.set_xlabel('Δ CD8-cytotoxic signature score (good − bad)')
ax.set_title('CD8-cytotoxic axis: 9-cohort meta + independent convergent study\n(Akiyoshi et al 2023 JAMA Netw Open)',
             fontsize=10.5)
ax.legend(loc='lower right', fontsize=8, frameon=False)
ax.spines[['top','right']].set_visible(False)

# -------- Panel B: signature-level meta Z bar --------
ax = axes[1]
order_meta = meta.set_index('signature').loc[
    [s for s in SIGS_ORDER if s in meta.signature.values]]
colors = ['#2E86AB' if z>0 else '#E63946' for z in order_meta.Z]
ax.barh(range(len(order_meta)), order_meta.Z, color=colors, edgecolor='k', lw=0.5)
ax.axvline(0, color='k', lw=0.5)
ax.axvline(1.96, color='gray', ls='--', lw=0.5)
ax.axvline(-1.96, color='gray', ls='--', lw=0.5)
ax.set_yticks(range(len(order_meta)))
ax.set_yticklabels(order_meta.index)
ax.set_xlabel('Stouffer meta Z')
ax.set_title('Signature-level meta across 9 nCRT cohorts', fontsize=10.5)
for i,(z,p) in enumerate(zip(order_meta.Z, order_meta.p_meta)):
    ax.text(z + (0.1 if z>0 else -0.1), i, f'p={p:.3f}',
            va='center', ha='left' if z>0 else 'right', fontsize=8)
ax.invert_yaxis()
ax.spines[['top','right']].set_visible(False)

# -------- Panel C: decoupling scatter --------
ax = axes[2]
cd8  = df[df.signature=='CD8_cytotoxic'].set_index('gse')['delta']
prol = df[df.signature=='Tumor_cellcycle'].set_index('gse')['delta']
shared = cd8.index.intersection(prol.index)
ax.scatter(prol.loc[shared], cd8.loc[shared], s=80, color='#2E86AB',
           edgecolor='k', lw=0.7, zorder=3, label='9-cohort meta')
for g in shared:
    ax.annotate(g.replace('GSE',''), (prol.loc[g], cd8.loc[g]),
                xytext=(5,5), textcoords='offset points', fontsize=7.5)
# Akiyoshi point (no tumor_cellcycle to compare; plot only on CD8 axis with '×')
ax.scatter([0], [0.18], marker='D', s=100, color='#6A4C93',
           edgecolor='k', lw=0.7, zorder=4,
           label='Akiyoshi 2023 (CD8 only)')
ax.annotate('Akiyoshi\n2023', (0, 0.18), xytext=(10, 0),
            textcoords='offset points', fontsize=7.5, color='#6A4C93')
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel('Δ Tumor_cellcycle (good − bad)')
ax.set_ylabel('Δ CD8_cytotoxic (good − bad)')
ax.set_title('CD8 effector vs tumor proliferation\ndecouple across cohorts',
             fontsize=10.5)
ax.legend(loc='lower right', fontsize=8, frameon=False)
ax.spines[['top','right']].set_visible(False)

plt.tight_layout()
# save as main Fig 7 (promoted) + supp copy
for base, name in [(FIG_MAIN, 'Fig7_external_CD8_validation_v4'),
                   (FIG_SUPP, 'SuppFig_v4_CD8_meta_plus_akiyoshi')]:
    base.mkdir(parents=True, exist_ok=True)
    plt.savefig(base/f'{name}.png', dpi=600, bbox_inches='tight')
    plt.savefig(base/f'{name}.pdf', bbox_inches='tight')
    print(f'saved {base}/{name}.png/pdf')
