"""
Convergence test: do baseline tumor-intrinsic features predict the magnitude of the
post-CRT immune cascade?

Question (from user 2026-04-18): if DSB Repair / HDR / E2F / G2M / Myc V2 are high at
baseline, do those subjects show a stronger post-treatment immune cascade (CD8 cytotoxic,
Treg, MHC II, etc.)?

If yes: discovery LASSO (tumor-intrinsic) → triggers external-validated cascade (CD8 axis)
        → three threads converge on a single mechanism.
If no:  discovery and external validation are observationally independent — paper should be
        reframed as two parallel observations rather than one unified causal model.

Inputs:
  /data/data/TNT/analysis/260418_add/integrated_subject_master_v2.tsv  (35 subjects, 40 features)
  /data/data/TNT/analysis/09_integration/paired_delta/paired_feature_long.tsv  (long-format
      per-subject paired Δ for 8 legacy features, 14 paired subjects)
  /data/data/TNT/analysis/260418_add/paired_immune_delta_per_subject.tsv  (4 new immune sigs Δ
      for 12 paired subjects)

Outputs:
  baseline_vs_delta_corr_long.tsv   — every (baseline, Δ) Spearman + Pearson + partial-r
                                       controlled for response_bin
  baseline_vs_delta_corr_top.tsv    — sorted by |partial-r|, FDR
  fig_baseline_vs_cascade_heatmap.{pdf,png}
  fig_baseline_vs_cascade_top_scatter.{pdf,png}  — top 6 (baseline, Δ) pairs
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

# ---------- Load baseline (pre-treatment, per-subject) ----------
M = pd.read_csv(f'{OUT}/integrated_subject_master_v2.tsv', sep='\t')
M['subject_id'] = M['subject_id'].astype(str)

# Tumor-intrinsic candidates from baseline (the LASSO predictor's universe)
TUMOR_INTRINSIC = [
    'TMB_nonsyn_per_Mb', 'CIN', 'frac_amp', 'frac_del', 'MSI_pct',
    'CD8_proliferation',  # legacy contaminated, included for sanity check
    'MHC_II', 'Antigen_presentation', 'NLRC5_HLA_IFNG',
    'EMT_Mak', 'Hypoxia_Buffa', 'Stemness_mRNAsi_proxy',
    'DNA Double-Strand Break Repair R-HSA-5693532',
    'HDR Thru Homologous Recombination (HRR) R-HSA-5685942',
    'DNA Repair R-HSA-73894',
    'E2F Targets', 'G2-M Checkpoint', 'Myc Targets V2',
    'Epithelial Mesenchymal Transition', 'TGF-beta Signaling', 'Hypoxia',
    'CD8_cytotoxic', 'Tcell_infiltration', 'Bcell_infiltration',
]
TUMOR_INTRINSIC = [c for c in TUMOR_INTRINSIC if c in M.columns]
print(f'baseline features: {len(TUMOR_INTRINSIC)}')

# ---------- Load paired Δ — long format ----------
LONG = pd.read_csv('/data/data/TNT/analysis/09_integration/paired_delta/paired_feature_long.tsv', sep='\t')
LONG['subject_id'] = LONG['subject_id'].astype(str)
LONG['delta'] = LONG['post'] - LONG['pre']
delta_legacy = LONG.pivot(index='subject_id', columns='feature', values='delta')
delta_legacy.columns = [f'{c}_delta' for c in delta_legacy.columns]
print(f'legacy paired Δ: {delta_legacy.shape}')

# Add new immune Δ
DNEW = pd.read_csv(f'{OUT}/paired_immune_delta_per_subject.tsv', sep='\t')
DNEW['subject_id'] = DNEW['subject_id'].astype(str)
new_delta_cols = [c for c in DNEW.columns if c.endswith('_delta')]
delta_new = DNEW.set_index('subject_id')[new_delta_cols]
print(f'new immune Δ: {delta_new.shape}')

# Combine all Δ
delta = delta_legacy.join(delta_new, how='outer')
print(f'combined Δ: {delta.shape}, subjects: {sorted(delta.index, key=int)}')

# Merge with baseline (only paired subjects)
paired_subjects = sorted(delta.index, key=int)
B = M.set_index('subject_id').loc[paired_subjects, TUMOR_INTRINSIC + ['response_bin']].copy()
B['y_good'] = (B['response_bin']=='good').astype(int)
print(f'paired baseline: {B.shape}')

# ---------- Compute correlations ----------
def partial_corr_spearman(x, y, z):
    """Spearman partial correlation of x and y controlling for z (binary)."""
    df = pd.DataFrame({'x':x,'y':y,'z':z}).dropna()
    if len(df) < 5: return np.nan, np.nan
    # Rank, then linear regress out z
    rx = stats.rankdata(df.x)
    ry = stats.rankdata(df.y)
    z_ = df.z.values - df.z.values.mean()
    rx_res = rx - np.polyval(np.polyfit(z_, rx, 1), z_)
    ry_res = ry - np.polyval(np.polyfit(z_, ry, 1), z_)
    if rx_res.std()==0 or ry_res.std()==0: return np.nan, np.nan
    r, p = stats.pearsonr(rx_res, ry_res)
    return r, p

rows = []
delta_cols = list(delta.columns)
for bf in TUMOR_INTRINSIC:
    for df_col in delta_cols:
        x = B[bf].values
        y = delta[df_col].reindex(B.index).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 6: continue
        x_v, y_v = x[valid], y[valid]
        z_v = B['y_good'].values[valid]
        rs, ps = stats.spearmanr(x_v, y_v)
        rp, pp_ = stats.pearsonr(x_v, y_v)
        partial_r, partial_p = partial_corr_spearman(x_v, y_v, z_v)
        rows.append({'baseline': bf, 'delta': df_col, 'n_paired': int(valid.sum()),
                     'spearman_r': round(rs,3), 'spearman_p': round(ps,4),
                     'pearson_r':  round(rp,3),  'pearson_p':  round(pp_,4),
                     'partial_r_resp_adj': round(partial_r,3) if not np.isnan(partial_r) else np.nan,
                     'partial_p_resp_adj': round(partial_p,4) if not np.isnan(partial_p) else np.nan})

R = pd.DataFrame(rows)
# FDR within partial_p (response-adjusted)
mask = R['partial_p_resp_adj'].notna()
qvals = np.full(len(R), np.nan)
qvals[mask] = multipletests(R.loc[mask,'partial_p_resp_adj'], method='fdr_bh')[1]
R['partial_q_BH'] = qvals.round(4)
R = R.sort_values('partial_p_resp_adj')
R.to_csv(f'{OUT}/baseline_vs_delta_corr_long.tsv', sep='\t', index=False)

print(f'\n=== top 20 by response-adjusted partial Spearman P ===')
print(R.head(20).to_string(index=False))
print(f'\n=== count of |partial_r| > 0.5 with p<0.05 ===')
hit = R[(R.partial_r_resp_adj.abs()>0.5) & (R.partial_p_resp_adj<0.05)]
print(f'  {len(hit)} pairs')
hit.to_csv(f'{OUT}/baseline_vs_delta_corr_top.tsv', sep='\t', index=False)

# ---------- Heatmap (baseline × Δ partial-r) ----------
Hmat = R.pivot(index='baseline', columns='delta', values='partial_r_resp_adj')
Hmat = Hmat.reindex(index=TUMOR_INTRINSIC)
Hmat = Hmat.dropna(how='all', axis=1).dropna(how='all', axis=0)

fig, ax = plt.subplots(figsize=(max(7, 0.6*Hmat.shape[1]+3), max(5, 0.32*Hmat.shape[0]+1.5)))
v = max(0.6, np.nanmax(np.abs(Hmat.values)))
im = ax.imshow(Hmat.values, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto')
ax.set_xticks(range(Hmat.shape[1])); ax.set_xticklabels(Hmat.columns, rotation=45, ha='right', fontsize=8)
ax.set_yticks(range(Hmat.shape[0])); ax.set_yticklabels([c[:42] for c in Hmat.index], fontsize=8)
# Annotate significant cells
Pmat = R.pivot(index='baseline', columns='delta', values='partial_p_resp_adj').reindex_like(Hmat)
for i in range(Hmat.shape[0]):
    for j in range(Hmat.shape[1]):
        v_ = Hmat.values[i,j]; p_ = Pmat.values[i,j]
        if not np.isnan(v_):
            mark = '**' if (not np.isnan(p_) and p_<0.01) else ('*' if (not np.isnan(p_) and p_<0.05) else '')
            ax.text(j, i, f'{v_:+.2f}{mark}', ha='center', va='center',
                    color='white' if abs(v_)>0.45 else 'black', fontsize=6.5)
plt.colorbar(im, ax=ax, label='Spearman partial-r (response-adjusted)')
ax.set_title(f'Baseline tumor-intrinsic vs paired Δ — n={len(B)} paired subjects\n(* p<0.05, ** p<0.01, response-adjusted partial Spearman)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_baseline_vs_cascade_heatmap.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote Fig_baseline_vs_cascade_heatmap')

# ---------- Top scatter panels ----------
top = R.dropna(subset=['partial_p_resp_adj']).head(9)
ncol = 3; nrow = int(np.ceil(len(top)/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.0*ncol, 2.7*nrow))
axes = axes.flatten() if nrow*ncol > 1 else [axes]
for ax, (_, row) in zip(axes, top.iterrows()):
    bf, df_col = row['baseline'], row['delta']
    x = B[bf].values
    y = delta[df_col].reindex(B.index).values
    valid = ~(np.isnan(x) | np.isnan(y))
    x, y = x[valid], y[valid]
    cols = ['#2E86AB' if g else '#E63946' for g in B['y_good'].values[valid]]
    ax.scatter(x, y, c=cols, s=40, edgecolor='black', lw=0.4)
    if len(x) >= 3:
        m, b_ = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, m*xs+b_, color='gray', lw=0.8, ls='--')
    ax.set_xlabel(bf[:32], fontsize=7.5)
    ax.set_ylabel(df_col, fontsize=7.5)
    ax.set_title(f"r={row['spearman_r']:+.2f} (p={row['spearman_p']:.3f})\npartial r={row['partial_r_resp_adj']:+.2f} (p={row['partial_p_resp_adj']:.3f})",
                 fontsize=7.5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
for ax in axes[len(top):]: ax.set_visible(False)
fig.suptitle('Top 9 baseline × cascade Δ pairs by response-adjusted partial Spearman P\n(blue = good, red = bad)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_baseline_vs_cascade_top_scatter.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote Fig_baseline_vs_cascade_top_scatter')

print(f'\n=== quick verdict ===')
n_sig05 = (R['partial_p_resp_adj']<0.05).sum()
n_sig01 = (R['partial_p_resp_adj']<0.01).sum()
n_sigq  = (R['partial_q_BH']<0.10).sum()
print(f'  pairs tested:                    {len(R)}')
print(f'  partial p<0.05:                  {n_sig05}')
print(f'  partial p<0.01:                  {n_sig01}')
print(f'  partial BH-q<0.10:               {n_sigq}')
print(f'  expected by chance at p<0.05:    {round(0.05*len(R),1)}')
