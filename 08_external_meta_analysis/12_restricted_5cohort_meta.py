"""
Restricted meta-analysis on the 5 cohorts that were >=3/4 concordant for Thread 1
(GSE56699, GSE45404, GSE133057, GSE87211, GSE35452; total N = 556).

Plus a paper-level Akiyoshi 2023 (GSE216616, n=298) row for the CD8 axis (Thread 2),
since GSE216616 GEO metadata does NOT contain per-sample TRG labels (paper Supplement
only, JAMA paywalled) so Stouffer combination is not possible at the per-sample level.

Output:
  restricted5_meta_thread1.tsv   — meta Z + p for the 4 Thread 1 features
  restricted5_meta_thread2.tsv   — meta Z + p for the 3 immune (Thread 2) features
  Fig_restricted5_meta_forest.{pdf,png}
"""
import os
import numpy as np, pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/data/data/TNT/analysis/260418_add'
SRC = '/data/data/TNT/analysis/11_external_validation/v3_signature_response_stats.tsv'

KEEP = ['GSE56699','GSE45404','GSE133057','GSE87211','GSE35452']
THREAD1 = ['DSB_HDR_repair','E2F_MYC_cellcycle','Tumor_cellcycle','EMT']
THREAD2 = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration']
EXPECTED = {'DSB_HDR_repair':+1,'E2F_MYC_cellcycle':+1,'Tumor_cellcycle':+1,'EMT':-1,
            'CD8_cytotoxic':+1,'Tcell_infiltration':+1,'Bcell_infiltration':+1}

df = pd.read_csv(SRC, sep='\t')
sub = df[df.gse.isin(KEEP)].copy()

def stouffer(rows):
    """Two-sided p, signed by concordance with expected direction."""
    w = np.sqrt(rows.n_good.values + rows.n_bad.values)
    z_per = np.array([
        norm.isf(p/2) * np.sign(EXPECTED[s] * d)
        for s, p, d in zip(rows.signature, rows.pvalue, rows.delta)
    ])
    Z = np.sum(w * z_per) / np.sqrt(np.sum(w**2))
    p = 2 * (1 - norm.cdf(abs(Z)))
    return float(Z), float(p)

# ---- Thread 1 restricted meta ----
rows = []
for sig in THREAD1:
    s = sub[sub.signature == sig]
    Z, p = stouffer(s)
    rows.append({'thread': 'Thread1_tumor_intrinsic', 'signature': sig,
                 'n_cohorts': len(s), 'total_n': int((s.n_good + s.n_bad).sum()),
                 'Z': round(Z, 2), 'p_meta': round(p, 4),
                 'cohorts': ','.join(s.gse.tolist())})

# ---- Thread 2 restricted meta ----
for sig in THREAD2:
    s = sub[sub.signature == sig]
    Z, p = stouffer(s)
    rows.append({'thread': 'Thread2_immune', 'signature': sig,
                 'n_cohorts': len(s), 'total_n': int((s.n_good + s.n_bad).sum()),
                 'Z': round(Z, 2), 'p_meta': round(p, 4),
                 'cohorts': ','.join(s.gse.tolist())})

R = pd.DataFrame(rows)
R.to_csv(f'{OUT}/restricted5_meta_combined.tsv', sep='\t', index=False)
print('=== Restricted 5-cohort meta (>=3/4 concordant on Thread 1) ===')
print(f'Cohorts: {KEEP}')
print(f'Total N (good+bad across all features): {int((sub[sub.signature=="CD8_cytotoxic"].n_good + sub[sub.signature=="CD8_cytotoxic"].n_bad).sum())}')
print()
print(R.to_string(index=False))

# ---- Akiyoshi 2023 paper-level row (Thread 2 only, CD8 axis) ----
print('\n--- Akiyoshi 2023 (GSE216616, n=298) — paper-level ---')
print('Cytotoxic lymphocyte score: good 0.76 vs bad 0.58 (Δ=0.18)')
print('OR=3.81 [1.82, 7.97], GZMA×PRF1 P=0.005')
print('NOT Stouffer-combined (per-sample TRG labels not in GEO; JAMA Supp only)')

# Save the per-cohort detail too
sub.to_csv(f'{OUT}/restricted5_per_cohort_detail.tsv', sep='\t', index=False)

# ---------- Forest plot ----------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                          gridspec_kw={'width_ratios':[1.0, 1.0]})

palette_t1 = {'DSB_HDR_repair':'#1f77b4', 'E2F_MYC_cellcycle':'#9467bd',
              'Tumor_cellcycle':'#2ca02c', 'EMT':'#d62728'}
palette_t2 = {'CD8_cytotoxic':'#E63946', 'Tcell_infiltration':'#FCA311',
              'Bcell_infiltration':'#6A4C93'}

def plot_forest(ax, sigs, palette, title):
    cohorts = KEEP
    n_sig = len(sigs)
    for i, sig in enumerate(sigs):
        rs = sub[sub.signature == sig]
        for j, gse in enumerate(cohorts):
            row = rs[rs.gse == gse]
            if len(row) == 0: continue
            d = row.delta.values[0]
            n = row.n_good.values[0] + row.n_bad.values[0]
            se = 1.96 / np.sqrt(n)  # rough 95% CI proxy
            y = j + (i - (n_sig-1)/2) * 0.16
            color = palette[sig]
            ax.errorbar([d], [y], xerr=[se], fmt='o', color=color, ms=6, lw=1.0,
                        capsize=3, markerfacecolor=color, markeredgecolor='black', mew=0.4)
        # meta diamond
        meta = R[R.signature == sig].iloc[0]
        meta_y = len(cohorts) + 0.5 + (i - (n_sig-1)/2) * 0.16
        # diamond from delta-weighted mean for graphical
        weights = np.sqrt(rs.n_good.values + rs.n_bad.values)
        meta_d = np.sum(rs.delta.values * weights) / np.sum(weights)
        ax.plot([meta_d - 0.05, meta_d, meta_d + 0.05, meta_d, meta_d - 0.05],
                [meta_y, meta_y + 0.06, meta_y, meta_y - 0.06, meta_y],
                color=color, lw=1.5)
        ax.fill([meta_d - 0.05, meta_d, meta_d + 0.05, meta_d],
                [meta_y, meta_y + 0.06, meta_y, meta_y - 0.06],
                color=color, alpha=0.5)
        ax.text(meta_d + 0.08, meta_y, f'{sig}: Z={meta.Z:+.2f}, P={meta.p_meta:.3f}',
                fontsize=8, va='center', color=color)
    ax.axvline(0, color='black', lw=0.6)
    ax.set_yticks(list(range(len(cohorts))) + [len(cohorts) + 0.5])
    ax.set_yticklabels(cohorts + ['Stouffer meta'], fontsize=8.5)
    ax.set_xlabel('Δ (good − bad), z-score units')
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.9, 1.1)
    for s in ['top','right']: ax.spines[s].set_visible(False)

plot_forest(axes[0], THREAD1, palette_t1, 'Thread 1 — tumor-intrinsic\n(restricted to 5 concordant cohorts)')
plot_forest(axes[1], THREAD2, palette_t2, 'Thread 2 — immune (CD8 axis)\n(same 5 cohorts)')

# Akiyoshi row on Thread 2 panel
ax = axes[1]
ak_y = len(KEEP) + 1.4
ax.errorbar([0.18], [ak_y], xerr=[0.07], fmt='D', color='#6A4C93', ms=8, lw=1.6,
            capsize=4, markerfacecolor='white', markeredgecolor='#6A4C93', mew=2)
ax.text(0.18 + 0.10, ak_y,
        'Akiyoshi 2023 (GSE216616, n=298)\nOR=3.81 [1.82, 7.97], paper-level',
        fontsize=8, va='center', color='#6A4C93', style='italic')
# Add Akiyoshi tick by extending y-tick set
existing_ticks = list(range(len(KEEP))) + [len(KEEP) + 0.5]
existing_labels = KEEP + ['Stouffer meta']
ax.set_yticks(existing_ticks + [ak_y])
ax.set_yticklabels(existing_labels + ['Akiyoshi paper-lvl'], fontsize=8.5)
ax.set_ylim(-0.5, ak_y + 0.6)

fig.suptitle(f'Restricted 5-cohort external meta (N=556) + Akiyoshi 2023 (paper-level)',
             fontsize=11)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_restricted5_meta_forest.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'\nWrote restricted5_meta_combined.tsv, restricted5_per_cohort_detail.tsv, Fig_restricted5_meta_forest')
