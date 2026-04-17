"""
Per-cohort breakdown of Thread 1 (tumor-intrinsic) feature reproducibility across the 9
external validation cohorts (N = 721 total).

Discovery directions:
  DSB_HDR_repair      UP in good (LASSO β = +0.89, GSEA NES strongly +)
  E2F_MYC_cellcycle   UP in good (NES = 2.78 / 2.36)
  Tumor_cellcycle     UP in good (NES = 2.46 G2M)
  EMT                 DOWN in good (NES = -2.16)

Concordance = sign(external delta) matches sign expected from discovery.
'Significant concordant' = p<0.05 AND concordant direction.
'Significant discordant' = p<0.05 AND opposite direction (counts AGAINST discovery).
"""
import pandas as pd
import numpy as np

OUT = '/data/data/TNT/analysis/260418_add'
SRC = '/data/data/TNT/analysis/11_external_validation/v3_signature_response_stats.tsv'
M   = '/data/data/TNT/analysis/11_external_validation/v3_cohort_summary.tsv'

df = pd.read_csv(SRC, sep='\t')
meta = pd.read_csv(M, sep='\t')[['gse','n_samples','n_good','n_bad','regimen']]

# Thread 1 features and their expected directions (sign of expected delta good - bad)
THREAD1 = {'DSB_HDR_repair': +1, 'E2F_MYC_cellcycle': +1,
           'Tumor_cellcycle': +1, 'EMT': -1}

t1 = df[df.signature.isin(THREAD1)].copy()
t1['expected_sign'] = t1['signature'].map(THREAD1)
t1['observed_sign'] = np.sign(t1['delta'])
t1['concordant'] = (t1['expected_sign'] == t1['observed_sign']).astype(int)
t1['sig_concordant'] = ((t1['pvalue'] < 0.05) & (t1['concordant'] == 1)).astype(int)
t1['sig_discordant'] = ((t1['pvalue'] < 0.05) & (t1['concordant'] == 0)).astype(int)

# Per-cohort summary
per_cohort = t1.groupby('gse').agg(
    n_concordant=('concordant', 'sum'),
    n_sig_concordant=('sig_concordant', 'sum'),
    n_sig_discordant=('sig_discordant', 'sum'),
    deltas=('delta', lambda v: ', '.join(f'{x:+.2f}' for x in v)),
    pvalues=('pvalue', lambda v: ', '.join(f'{x:.3f}' for x in v)),
).reset_index()
per_cohort = per_cohort.merge(meta, on='gse', how='left')
per_cohort = per_cohort.sort_values(['n_sig_concordant','n_concordant'], ascending=[False, False])

# Tag verdict
def tag(row):
    if row.n_sig_concordant >= 1 and row.n_concordant >= 3:
        return 'STRONG concordant'
    if row.n_concordant == 4 and row.n_sig_discordant == 0:
        return 'all-trend concordant'
    if row.n_concordant == 3 and row.n_sig_discordant == 0:
        return 'mostly concordant'
    if row.n_sig_discordant >= 1:
        return 'SIGNIFICANT DISCORDANT'
    return 'mostly discordant'
per_cohort['verdict'] = per_cohort.apply(tag, axis=1)

per_cohort.to_csv(f'{OUT}/thread1_per_cohort_summary.tsv', sep='\t', index=False)

# Pretty print to console
print('Thread 1 per-cohort breakdown (DSB / E2F_MYC / Tumor_cellcycle / EMT)')
print('Discovery direction: DSB+ E2F+ Cellcycle+ EMT-')
print()
cols = ['gse','n_samples','n_good','n_bad','n_concordant','n_sig_concordant','n_sig_discordant','verdict','deltas','pvalues']
print(per_cohort[cols].to_string(index=False))

# Wide table per signature
wide_d = t1.pivot(index='gse', columns='signature', values='delta')
wide_p = t1.pivot(index='gse', columns='signature', values='pvalue')
wide = pd.concat({'delta': wide_d, 'pvalue': wide_p}, axis=1)
wide = wide.swaplevel(axis=1).sort_index(axis=1).reset_index()
wide.to_csv(f'{OUT}/thread1_per_cohort_wide.tsv', sep='\t', index=False)

# Concordant-cohort meta (for context, NOT to be used as cherry-picked main result)
concordant_cohorts = per_cohort[per_cohort.n_concordant >= 3]['gse'].tolist()
print(f'\nCohorts with >=3/4 concordant directions: {concordant_cohorts}')
print('  total N (good+bad):', per_cohort[per_cohort.gse.isin(concordant_cohorts)][['n_good','n_bad']].sum().sum())
print('\n=== If meta restricted to concordant cohorts (informative, NOT for cherry-picking) ===')
for sig in THREAD1:
    sub = t1[(t1.signature==sig) & (t1.gse.isin(concordant_cohorts))]
    n = sub.shape[0]
    if n == 0: continue
    weights = np.sqrt(sub.n_good.values + sub.n_bad.values)
    z_per = np.array([np.sign(s.expected_sign * s.delta) * abs(np.sqrt(2)*np.abs(np.percentile if False else np.log(max(s.pvalue,1e-300))/(-2)))
                     for _, s in sub.iterrows()])  # placeholder simplified
    # Use simple Stouffer with two-sided p but signed by expected direction
    from scipy.stats import norm
    z_per = np.array([norm.isf(s.pvalue/2) * np.sign(s.expected_sign * s.delta) for _, s in sub.iterrows()])
    Z = np.sum(weights * z_per) / np.sqrt(np.sum(weights**2))
    p = 2 * (1 - norm.cdf(abs(Z)))
    print(f'  {sig:22s}  Stouffer Z = {Z:+.2f}  P = {p:.3f}  ({n} cohorts)')
print('\nNB: the unrestricted 9-cohort meta gives Z ~0.7-1.3 (p>0.19) for all four Thread 1 features.')
