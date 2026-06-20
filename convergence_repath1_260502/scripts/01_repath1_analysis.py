"""
Repath-1 convergence test (2026-05-02).

Tests two suspicions about why the original 36-pair convergence test (manuscript §3.9 / Supp Fig S20)
failed to detect baseline → cascade dependence:
  Suspicion 1: response-group adjustment kills between-group dependence.
  Suspicion 2: 36-pair multiple testing burden dilutes mechanistic signal.

Four configurations, all on the same 12 RNA-paired subjects (intersection):
  Original   — partial Spearman, response-adjusted, BH/36 (replicate manuscript)
  Config A   — pooled Spearman + Pearson, no group adjustment, BH/36
  Config B   — 5 pre-specified mechanistic pairs, plain + partial side by side, BH/5
  Config C   — 5 pairs, between-group direction concordance + sign test
"""
import os, json, sys
from datetime import datetime
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/data/data/TNT/analysis'
OUT  = f'{ROOT}/convergence_repath1_260502'
LOG  = open(f'{OUT}/logs/01_run.log', 'w')

def log(*a):
    msg = ' '.join(str(x) for x in a)
    print(msg); LOG.write(msg+'\n'); LOG.flush()

log(f'[{datetime.now().isoformat(timespec="seconds")}] Repath-1 convergence test starting')
log(f'OUT = {OUT}')

# ---------- Load data ----------
M = pd.read_csv(f'{ROOT}/260418_add/integrated_subject_master_v2.tsv', sep='\t')
M['subject_id'] = M['subject_id'].astype(str)

LONG = pd.read_csv(f'{ROOT}/09_integration/paired_delta/paired_feature_long.tsv', sep='\t')
LONG['subject_id'] = LONG['subject_id'].astype(str)
LONG['delta'] = LONG['post'] - LONG['pre']
delta_legacy = LONG.pivot(index='subject_id', columns='feature', values='delta')
delta_legacy.columns = [f'{c}_delta' for c in delta_legacy.columns]

DNEW = pd.read_csv(f'{ROOT}/260418_add/paired_immune_delta_per_subject.tsv', sep='\t')
DNEW['subject_id'] = DNEW['subject_id'].astype(str)
new_cols = [c for c in DNEW.columns if c.endswith('_delta')]
delta_new = DNEW.set_index('subject_id')[new_cols]

# ---------- Define paired subject set: intersection (RNA-paired AND legacy paired) ----------
rna_paired = set(delta_new.index)
legacy_paired = set(delta_legacy.index)
intersect = sorted(rna_paired & legacy_paired, key=int)
log(f'RNA-paired subjects (n={len(rna_paired)}): {sorted(rna_paired, key=int)}')
log(f'Legacy paired subjects (n={len(legacy_paired)}): {sorted(legacy_paired, key=int)}')
log(f'Intersection used for ALL configs (n={len(intersect)}): {intersect}')

# Combine all Δ columns; subset to intersection rows
delta = delta_legacy.join(delta_new, how='outer').loc[intersect]
log(f'Delta matrix shape (intersection): {delta.shape}; Δ columns: {sorted(delta.columns)}')

# Baseline matrix (intersection rows)
B_full = M.set_index('subject_id').loc[intersect]
B_full['y_good'] = (B_full['response_bin']=='good').astype(int)
log(f'Baseline matrix (intersection): {B_full.shape}')
log(f'Response: good={int(B_full["y_good"].sum())} / bad={int((1-B_full["y_good"]).sum())}')

# ---------- Helper: partial Spearman controlling for binary z ----------
def partial_corr_spearman(x, y, z):
    df = pd.DataFrame({'x':x,'y':y,'z':z}).dropna()
    if len(df) < 5: return np.nan, np.nan, 0
    rx = stats.rankdata(df.x); ry = stats.rankdata(df.y)
    z_ = df.z.values - df.z.values.mean()
    rx_res = rx - np.polyval(np.polyfit(z_, rx, 1), z_)
    ry_res = ry - np.polyval(np.polyfit(z_, ry, 1), z_)
    if rx_res.std()==0 or ry_res.std()==0: return np.nan, np.nan, len(df)
    r, p = stats.pearsonr(rx_res, ry_res)
    return r, p, len(df)

# ============================================================
# Configuration ORIGINAL — replicate manuscript §3.9 / Supp Fig S20
# 9 baseline × 4 cascade Δ = 36 pairs, partial Spearman group-adjusted, BH/36
# ============================================================
log('\n' + '='*60)
log('Configuration ORIGINAL (replicate of script 09)')
log('='*60)
BASE_36 = ['DNA Double-Strand Break Repair R-HSA-5693532',
           'HDR Thru Homologous Recombination (HRR) R-HSA-5685942',
           'DNA Repair R-HSA-73894',
           'E2F Targets', 'G2-M Checkpoint', 'Myc Targets V2',
           'frac_amp', 'MHC_II', 'MSI_pct']
CASC_36 = ['CD8_cytotoxic_delta', 'Treg_delta', 'MHC_II_delta', 'IGH_n_delta']

rows = []
for bf in BASE_36:
    for cf in CASC_36:
        x = B_full[bf].values
        y = delta[cf].reindex(B_full.index).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 6: continue
        x_v, y_v = x[valid], y[valid]
        z_v = B_full['y_good'].values[valid]
        rs, ps = stats.spearmanr(x_v, y_v)
        pr, pp, n_pp = partial_corr_spearman(x_v, y_v, z_v)
        rows.append({'baseline_feature': bf, 'cascade_feature': cf, 'n': int(valid.sum()),
                     'spearman_r': round(rs,3), 'spearman_p': round(ps,4),
                     'partial_r': round(pr,3) if not np.isnan(pr) else np.nan,
                     'partial_P': round(pp,4) if not np.isnan(pp) else np.nan})
R_orig = pd.DataFrame(rows)
qvals = multipletests(R_orig['partial_P'].fillna(1.0), method='fdr_bh')[1]
R_orig['BH_q'] = qvals.round(4)
R_orig = R_orig.sort_values('partial_P').reset_index(drop=True)
R_orig.to_csv(f'{OUT}/tables/config_original_36pair.tsv', sep='\t', index=False)
log(f'Wrote: tables/config_original_36pair.tsv ({len(R_orig)} rows)')

# Sanity check: headline pair
hl = R_orig[(R_orig.baseline_feature=='DNA Double-Strand Break Repair R-HSA-5693532') &
            (R_orig.cascade_feature=='CD8_cytotoxic_delta')].iloc[0]
log(f'  Headline pair (DSB Repair × Δ CD8-cyt): partial_r={hl["partial_r"]:+.3f}, partial_P={hl["partial_P"]:.3f}, n={hl["n"]}')
log(f'  Manuscript expectation: partial_r ≈ −0.07, partial_P ≈ 0.83')
log(f'  ORIGINAL config: pairs<0.05={(R_orig.partial_P<0.05).sum()} / BH q<0.05={(R_orig.BH_q<0.05).sum()} / BH q<0.10={(R_orig.BH_q<0.10).sum()}')

# ============================================================
# Configuration A — pooled Spearman + Pearson, no group adjustment, BH/36
# ============================================================
log('\n' + '='*60)
log('Configuration A (pooled, no group adjustment)')
log('='*60)
rows = []
for bf in BASE_36:
    for cf in CASC_36:
        x = B_full[bf].values
        y = delta[cf].reindex(B_full.index).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 6: continue
        x_v, y_v = x[valid], y[valid]
        rs, ps = stats.spearmanr(x_v, y_v)
        rp, pp_pe = stats.pearsonr(x_v, y_v)
        rows.append({'baseline_feature': bf, 'cascade_feature': cf, 'n': int(valid.sum()),
                     'spearman_r': round(rs,3), 'spearman_P': round(ps,4),
                     'pearson_r': round(rp,3),  'pearson_P': round(pp_pe,4)})
R_A = pd.DataFrame(rows)
R_A['BH_q_spearman'] = multipletests(R_A['spearman_P'].fillna(1.0), method='fdr_bh')[1].round(4)
R_A['BH_q_pearson']  = multipletests(R_A['pearson_P'].fillna(1.0), method='fdr_bh')[1].round(4)
R_A = R_A.sort_values('spearman_P').reset_index(drop=True)
R_A.to_csv(f'{OUT}/tables/config_A_pooled_36pair.tsv', sep='\t', index=False)
log(f'Wrote: tables/config_A_pooled_36pair.tsv ({len(R_A)} rows)')

n_sp_05 = (R_A.spearman_P<0.05).sum(); n_sp_q05 = (R_A.BH_q_spearman<0.05).sum()
n_pe_05 = (R_A.pearson_P<0.05).sum();  n_pe_q05 = (R_A.BH_q_pearson<0.05).sum()
log(f'  Pooled Spearman:  {n_sp_05}/36 P<0.05, {n_sp_q05}/36 BH q<0.05')
log(f'  Pooled Pearson:   {n_pe_05}/36 P<0.05, {n_pe_q05}/36 BH q<0.05')
top5 = R_A.reindex(R_A.spearman_r.abs().sort_values(ascending=False).index).head(5)
log(f'  Top 5 |spearman_r|:')
for _, r in top5.iterrows():
    log(f'    {r.baseline_feature[:38]:38s} × {r.cascade_feature:24s}  r={r.spearman_r:+.3f} P={r.spearman_P:.4f} q={r.BH_q_spearman:.4f}')

# ---- Heatmap (Config A pooled Spearman) ----
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
Hmat = R_A.pivot(index='baseline_feature', columns='cascade_feature', values='spearman_r').reindex(index=BASE_36, columns=CASC_36)
Pmat = R_A.pivot(index='baseline_feature', columns='cascade_feature', values='spearman_P').reindex(index=BASE_36, columns=CASC_36)
fig, ax = plt.subplots(figsize=(5.5, 5.5))
v = max(0.6, np.nanmax(np.abs(Hmat.values)))
im = ax.imshow(Hmat.values, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto')
ax.set_xticks(range(len(CASC_36))); ax.set_xticklabels(CASC_36, rotation=30, ha='right', fontsize=8)
ax.set_yticks(range(len(BASE_36))); ax.set_yticklabels([c[:35] for c in BASE_36], fontsize=8)
for i in range(Hmat.shape[0]):
    for j in range(Hmat.shape[1]):
        v_ = Hmat.values[i,j]; p_ = Pmat.values[i,j]
        mark = '**' if (not np.isnan(p_) and p_<0.01) else ('*' if (not np.isnan(p_) and p_<0.05) else '')
        ax.text(j, i, f'{v_:+.2f}{mark}', ha='center', va='center',
                color='white' if abs(v_)>0.45 else 'black', fontsize=8)
plt.colorbar(im, ax=ax, label='Spearman r (pooled, no group adj.)')
ax.set_title(f'Configuration A: pooled Spearman, no group adjustment\n9 baseline × 4 cascade Δ (n=12), * P<0.05, ** P<0.01', fontsize=9)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/figures/FigA_pooled_heatmap.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
log('Wrote: figures/FigA_pooled_heatmap.{pdf,png}')

# ============================================================
# Configuration B — 5 mechanistic pairs, plain + partial side by side, BH/5
# ============================================================
log('\n' + '='*60)
log('Configuration B (5 mechanistic pairs, BH/5)')
log('='*60)

# Mechanistic 5 pairs (column-name mapping noted in README)
MECHANISTIC = [
    ('DNA Double-Strand Break Repair R-HSA-5693532', 'SBS5_delta',
     'DSB repair → radiation-induced SBS5 mutation clearance'),
    ('DNA Double-Strand Break Repair R-HSA-5693532', 'neo_binders_delta',
     'DSB repair → cytotoxic clearance → neoantigen depletion'),
    ('E2F Targets', 'missense_delta',
     'proliferation → mutation accumulation rate (TMB ≈ missense count)'),
    ('CD8_cytotoxic', 'CD8_exhaustion_delta',
     'cytotoxic engagement → chronic antigen → exhaustion'),
    ('E2F Targets', 'Treg_delta',
     'proliferation-driven antigen release → Treg infiltration'),
]

rows = []
for bf, cf, rationale in MECHANISTIC:
    x = B_full[bf].values
    y = delta[cf].reindex(B_full.index).values
    valid = ~(np.isnan(x) | np.isnan(y))
    n = int(valid.sum())
    x_v, y_v = x[valid], y[valid]
    z_v = B_full['y_good'].values[valid]
    rs, ps = stats.spearmanr(x_v, y_v)
    pr, pp, _ = partial_corr_spearman(x_v, y_v, z_v)
    rows.append({'baseline_feature': bf, 'cascade_feature': cf, 'n': n,
                 'mechanistic_rationale': rationale,
                 'spearman_r_plain':   round(rs,3), 'spearman_P_plain':   round(ps,4),
                 'partial_r_groupadj': round(pr,3) if not np.isnan(pr) else np.nan,
                 'partial_P_groupadj': round(pp,4) if not np.isnan(pp) else np.nan})
R_B = pd.DataFrame(rows)
R_B['BH_q_plain']    = multipletests(R_B['spearman_P_plain'].fillna(1.0), method='fdr_bh')[1].round(4)
R_B['BH_q_groupadj'] = multipletests(R_B['partial_P_groupadj'].fillna(1.0), method='fdr_bh')[1].round(4)
R_B.to_csv(f'{OUT}/tables/config_B_mechanistic_5pair.tsv', sep='\t', index=False)
log(f'Wrote: tables/config_B_mechanistic_5pair.tsv ({len(R_B)} rows)')
for _, r in R_B.iterrows():
    log(f'  {r.baseline_feature[:38]:38s} × {r.cascade_feature:24s}  '
        f'plain r={r.spearman_r_plain:+.3f} P={r.spearman_P_plain:.4f} q={r.BH_q_plain:.4f}  | '
        f'partial r={r.partial_r_groupadj:+.3f} P={r.partial_P_groupadj:.4f} q={r.BH_q_groupadj:.4f}  n={r.n}')

# ---- 5-panel scatter (Config B) ----
fig, axes = plt.subplots(2, 3, figsize=(11, 7))
axes_list = axes.flatten().tolist()
for ax, (bf, cf, rationale), (_, row) in zip(axes_list[:5], MECHANISTIC, R_B.iterrows()):
    x = B_full[bf].values
    y = delta[cf].reindex(B_full.index).values
    valid = ~(np.isnan(x) | np.isnan(y))
    x, y = x[valid], y[valid]
    cols = ['#0a7d6e' if g else '#c53e1f' for g in B_full['y_good'].values[valid]]
    ax.scatter(x, y, c=cols, s=55, edgecolor='black', lw=0.5)
    if len(x) >= 3:
        m, b_ = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, m*xs+b_, color='gray', lw=0.8, ls='--')
    ax.set_xlabel(bf[:42], fontsize=8.5)
    ax.set_ylabel(cf, fontsize=8.5)
    ax.set_title(f"plain r={row.spearman_r_plain:+.2f} P={row.spearman_P_plain:.3f}\n"
                 f"partial r={row.partial_r_groupadj:+.2f} P={row.partial_P_groupadj:.3f}  (n={row.n})",
                 fontsize=8.5)
    ax.axhline(0, color='lightgray', lw=0.6, ls=':', zorder=0)
    for s in ['top','right']: ax.spines[s].set_visible(False)
axes_list[5].set_visible(False)
fig.suptitle('Configuration B: 5 pre-specified mechanistic pairs (good=teal, bad=coral)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/figures/FigB_mechanistic_scatter.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
log('Wrote: figures/FigB_mechanistic_scatter.{pdf,png}')

# ============================================================
# Configuration C — between-group direction concordance, sign test
# ============================================================
log('\n' + '='*60)
log('Configuration C (between-group direction concordance + sign test)')
log('='*60)

good_idx = B_full[B_full.y_good==1].index
bad_idx  = B_full[B_full.y_good==0].index
log(f'  good n={len(good_idx)}: {sorted(good_idx, key=int)}')
log(f'  bad  n={len(bad_idx)}:  {sorted(bad_idx, key=int)}')

rows = []
for bf, cf, rationale in MECHANISTIC:
    bx_good = B_full.loc[good_idx, bf].mean()
    bx_bad  = B_full.loc[bad_idx,  bf].mean()
    cy_good = delta[cf].reindex(good_idx).mean()
    cy_bad  = delta[cf].reindex(bad_idx).mean()
    base_dir = 'good>bad' if bx_good > bx_bad else 'bad>good'
    casc_dir = 'good>bad' if cy_good > cy_bad else 'bad>good'
    concordant = (base_dir == casc_dir)
    rows.append({'baseline_feature': bf, 'cascade_feature': cf,
                 'mechanistic_rationale': rationale,
                 'baseline_good_mean': round(bx_good, 4),
                 'baseline_bad_mean':  round(bx_bad,  4),
                 'baseline_direction': base_dir,
                 'cascade_good_mean':  round(cy_good, 4),
                 'cascade_bad_mean':   round(cy_bad,  4),
                 'cascade_direction':  casc_dir,
                 'concordant': bool(concordant)})
R_C = pd.DataFrame(rows)

# Sign test
n_pairs = len(R_C)
n_concord = int(R_C['concordant'].sum())
# One-sided binomial: H0 p=0.5, H1 p>0.5
sign_p = stats.binomtest(n_concord, n_pairs, 0.5, alternative='greater').pvalue
sign_p = float(sign_p)
log(f'  Direction concordant: {n_concord}/{n_pairs}; one-sided sign test (binomial p=0.5) P={sign_p:.4f}')

R_C.loc[len(R_C)] = {'baseline_feature': '__SIGN_TEST__', 'cascade_feature': '',
                     'mechanistic_rationale': f'one-sided binomial(p=0.5), H1 p>0.5',
                     'baseline_good_mean': np.nan, 'baseline_bad_mean': np.nan,
                     'baseline_direction': '', 'cascade_good_mean': np.nan, 'cascade_bad_mean': np.nan,
                     'cascade_direction': f'{n_concord}/{n_pairs} concordant',
                     'concordant': sign_p}
R_C.to_csv(f'{OUT}/tables/config_C_directionconc_5pair.tsv', sep='\t', index=False)
log(f'Wrote: tables/config_C_directionconc_5pair.tsv (5 + 1 summary row)')
for _, r in R_C.iloc[:-1].iterrows():
    log(f'  {r.baseline_feature[:38]:38s} × {r.cascade_feature:24s}  '
        f'baseline {r.baseline_direction} ({r.baseline_good_mean:+.3f} vs {r.baseline_bad_mean:+.3f})  '
        f'cascade {r.cascade_direction} ({r.cascade_good_mean:+.3f} vs {r.cascade_bad_mean:+.3f})  '
        f'concord={r.concordant}')

# ---------- Summary JSON for README ----------
summary = {
    'date': datetime.now().isoformat(timespec='seconds'),
    'n_subjects_intersection': len(intersect),
    'subjects_intersection': intersect,
    'good_subjects': sorted(good_idx.tolist(), key=int),
    'bad_subjects':  sorted(bad_idx.tolist(),  key=int),
    'original': {
        'headline_pair_partial_r': float(hl['partial_r']),
        'headline_pair_partial_P': float(hl['partial_P']),
        'headline_pair_n': int(hl['n']),
        'n_pairs_P_lt_0_05':  int((R_orig.partial_P<0.05).sum()),
        'n_pairs_BH_q_lt_0_05': int((R_orig.BH_q<0.05).sum()),
        'n_pairs_BH_q_lt_0_10': int((R_orig.BH_q<0.10).sum()),
    },
    'config_A': {
        'spearman_n_P_lt_0_05': int((R_A.spearman_P<0.05).sum()),
        'spearman_n_BH_q_lt_0_05': int((R_A.BH_q_spearman<0.05).sum()),
        'pearson_n_P_lt_0_05':  int((R_A.pearson_P<0.05).sum()),
        'pearson_n_BH_q_lt_0_05':  int((R_A.BH_q_pearson<0.05).sum()),
        'top5_strongest_abs_r': top5[['baseline_feature','cascade_feature','spearman_r','spearman_P','BH_q_spearman']].to_dict('records'),
    },
    'config_B': {
        'pairs': R_B.to_dict('records'),
    },
    'config_C': {
        'n_concordant':   n_concord,
        'n_total':        n_pairs,
        'sign_test_P':    sign_p,
        'pair_directions': R_C.iloc[:-1][['baseline_feature','cascade_feature','baseline_direction','cascade_direction','concordant']].to_dict('records'),
    },
}
with open(f'{OUT}/tables/_summary.json', 'w') as f:
    json.dump(summary, f, indent=2, default=str)
log('Wrote: tables/_summary.json')

log(f'\n[{datetime.now().isoformat(timespec="seconds")}] Repath-1 analysis complete')
LOG.close()
