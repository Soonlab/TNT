"""
Targeted convergence test (the user's actual question, not a fishing expedition).

The hypothesis we wanted to falsify is specific:
  "Baseline tumor-intrinsic features (the LASSO winners: DSB Repair, HDR, E2F,
   G2M, Myc V2, frac_amp) predict the magnitude of the post-CRT CD8/immune cascade."

If true → discovery predictor (Thread 1) → triggers external CD8 axis (Thread 2)
         → manifested as paired cascade (Thread 3) — convergent causal chain.

If false → three threads are three independent observations.

To avoid the multiple-testing dilution of script 08 (312 pairs), here we test only
the **9 baseline × 4 cascade pairs that the convergence hypothesis specifies**:

  Baseline (the four LASSO winners + the strongest GSEA hits):
    - DSB Repair (Reactome) ← LASSO winner
    - HDR
    - frac_amp ← LASSO winner
    - MHC_II ← LASSO winner
    - MSI_pct ← LASSO winner
    - E2F Targets
    - G2-M Checkpoint
    - Myc Targets V2
    - DNA Repair (Reactome)

  Cascade Δ (the externally validated + paired-significant axes):
    - CD8_cytotoxic_delta ← externally validated axis (key question)
    - Treg_delta ← only between-group robust paired finding
    - MHC_II_delta
    - IGH_n_delta ← BCR infiltration

Output:
  targeted_convergence_test.tsv
  Fig_targeted_convergence_heatmap.{pdf,png}
  Fig_targeted_DSB_vs_CD8cyt_scatter.{pdf,png}  ← THE plot the narrative depends on
"""
import os, warnings
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = '/data/data/TNT/analysis/260418_add'
M = pd.read_csv(f'{OUT}/integrated_subject_master_v2.tsv', sep='\t')
M['subject_id'] = M['subject_id'].astype(str)

# Cascade
LONG = pd.read_csv('/data/data/TNT/analysis/09_integration/paired_delta/paired_feature_long.tsv', sep='\t')
LONG['subject_id'] = LONG['subject_id'].astype(str)
LONG['delta'] = LONG['post'] - LONG['pre']
delta_legacy = LONG.pivot(index='subject_id', columns='feature', values='delta')
delta_legacy.columns = [f'{c}_delta' for c in delta_legacy.columns]
DNEW = pd.read_csv(f'{OUT}/paired_immune_delta_per_subject.tsv', sep='\t')
DNEW['subject_id'] = DNEW['subject_id'].astype(str)
delta_new = DNEW.set_index('subject_id')[[c for c in DNEW.columns if c.endswith('_delta')]]
delta = delta_legacy.join(delta_new, how='outer')

# Targeted feature sets
BASE = ['DNA Double-Strand Break Repair R-HSA-5693532',
        'HDR Thru Homologous Recombination (HRR) R-HSA-5685942',
        'DNA Repair R-HSA-73894',
        'E2F Targets', 'G2-M Checkpoint', 'Myc Targets V2',
        'frac_amp', 'MHC_II', 'MSI_pct']
CASC = ['CD8_cytotoxic_delta', 'Treg_delta', 'MHC_II_delta', 'IGH_n_delta']

paired_subjects = sorted(delta.index, key=int)
B = M.set_index('subject_id').loc[paired_subjects, BASE + ['response_bin']].copy()
B['y_good'] = (B['response_bin']=='good').astype(int)

def partial_corr_spearman(x, y, z):
    df = pd.DataFrame({'x':x,'y':y,'z':z}).dropna()
    if len(df) < 5: return np.nan, np.nan
    rx = stats.rankdata(df.x); ry = stats.rankdata(df.y)
    z_ = df.z.values - df.z.values.mean()
    rx_res = rx - np.polyval(np.polyfit(z_, rx, 1), z_)
    ry_res = ry - np.polyval(np.polyfit(z_, ry, 1), z_)
    if rx_res.std()==0 or ry_res.std()==0: return np.nan, np.nan
    return stats.pearsonr(rx_res, ry_res)

rows = []
for bf in BASE:
    for cf in CASC:
        x = B[bf].values
        y = delta[cf].reindex(B.index).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 6: continue
        rs, ps = stats.spearmanr(x[valid], y[valid])
        pr, pp = partial_corr_spearman(x[valid], y[valid], B['y_good'].values[valid])
        rows.append({'baseline': bf, 'cascade': cf, 'n_paired': int(valid.sum()),
                     'spearman_r': round(rs,3), 'spearman_p': round(ps,4),
                     'partial_r': round(pr,3) if not np.isnan(pr) else np.nan,
                     'partial_p': round(pp,4) if not np.isnan(pp) else np.nan})
R = pd.DataFrame(rows)
qvals = multipletests(R['partial_p'].fillna(1.0), method='fdr_bh')[1]
R['partial_q_BH'] = qvals.round(3)
R = R.sort_values('partial_p')
R.to_csv(f'{OUT}/targeted_convergence_test.tsv', sep='\t', index=False)
print(f'targeted pairs tested: {len(R)} (= {len(BASE)} baseline x {len(CASC)} cascade)')
print(R.to_string(index=False))

# ---------- Heatmap ----------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
Hmat = R.pivot(index='baseline', columns='cascade', values='partial_r').reindex(index=BASE, columns=CASC)
Pmat = R.pivot(index='baseline', columns='cascade', values='partial_p').reindex(index=BASE, columns=CASC)
fig, ax = plt.subplots(figsize=(5.5, 5.0))
v = max(0.6, np.nanmax(np.abs(Hmat.values)))
im = ax.imshow(Hmat.values, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto')
ax.set_xticks(range(len(CASC))); ax.set_xticklabels(CASC, rotation=30, ha='right', fontsize=8)
ax.set_yticks(range(len(BASE))); ax.set_yticklabels([c[:35] for c in BASE], fontsize=8)
for i in range(Hmat.shape[0]):
    for j in range(Hmat.shape[1]):
        v_ = Hmat.values[i,j]; p_ = Pmat.values[i,j]
        mark = '**' if (not np.isnan(p_) and p_<0.01) else ('*' if (not np.isnan(p_) and p_<0.05) else '')
        ax.text(j, i, f'{v_:+.2f}{mark}', ha='center', va='center',
                color='white' if abs(v_)>0.45 else 'black', fontsize=8)
plt.colorbar(im, ax=ax, label='partial-r (response-adj)')
ax.set_title(f'Targeted convergence test ({len(R)} pre-specified pairs, n=12-14 paired)\n* p<0.05, ** p<0.01, BH-q reported in TSV', fontsize=9)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_targeted_convergence_heatmap.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)

# ---------- THE scatter: DSB Repair baseline vs CD8_cytotoxic Δ ----------
fig, axes = plt.subplots(2, 2, figsize=(7, 6))
key_pairs = [
    ('DNA Double-Strand Break Repair R-HSA-5693532', 'CD8_cytotoxic_delta', 'DSB Repair (baseline) → CD8 cytotoxic Δ'),
    ('DNA Double-Strand Break Repair R-HSA-5693532', 'Treg_delta',          'DSB Repair (baseline) → Treg Δ'),
    ('HDR Thru Homologous Recombination (HRR) R-HSA-5685942', 'CD8_cytotoxic_delta', 'HDR (baseline) → CD8 cytotoxic Δ'),
    ('E2F Targets',                                  'CD8_cytotoxic_delta', 'E2F Targets (baseline) → CD8 cytotoxic Δ'),
]
for ax, (bf, cf, title) in zip(axes.flatten(), key_pairs):
    x = B[bf].values
    y = delta[cf].reindex(B.index).values
    valid = ~(np.isnan(x) | np.isnan(y))
    x, y = x[valid], y[valid]
    cols = ['#2E86AB' if g else '#E63946' for g in B['y_good'].values[valid]]
    ax.scatter(x, y, c=cols, s=55, edgecolor='black', lw=0.4)
    if len(x) >= 3:
        m, b_ = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, m*xs+b_, color='gray', lw=0.8, ls='--')
    rs, ps = stats.spearmanr(x, y)
    ax.set_xlabel(bf[:35], fontsize=8.5)
    ax.set_ylabel(cf, fontsize=8.5)
    ax.set_title(f'{title}\nSpearman r={rs:+.2f}, p={ps:.3f} (n={len(x)})', fontsize=9)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle('Convergence test: do baseline tumor-intrinsic features drive the immune cascade?\n(blue = good responder, red = bad responder)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_targeted_DSB_vs_CD8cyt_scatter.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('\nWrote heatmap + 4-panel key scatter')

# ---------- Verdict ----------
print(f'\n=== TARGETED VERDICT ===')
n_sig05 = (R['partial_p']<0.05).sum()
n_sig01 = (R['partial_p']<0.01).sum()
n_sigq  = (R['partial_q_BH']<0.10).sum()
print(f'  pairs tested:            {len(R)}')
print(f'  partial p<0.05:          {n_sig05}')
print(f'  partial p<0.01:          {n_sig01}')
print(f'  partial BH-q<0.10:       {n_sigq}')
print(f'  expected by chance @0.05: {round(0.05*len(R),1)}')
print(f'\n  -> If 0 hits at FDR<0.10 AND <= chance hits at p<0.05, convergence hypothesis FAILS.')
print(f'     The discovery predictor and the paired cascade are observationally independent.')
