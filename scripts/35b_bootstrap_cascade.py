"""BCa bootstrap 95% CIs for cascade paired-delta claims.
Also computes between-group bootstrap on (good_delta - bad_delta).
Addresses reviewer Major #1.
"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats as st
import warnings; warnings.filterwarnings('ignore')

BASE = Path('/mnt/sda1/data/TNT/analysis')
TBL  = BASE/'tables'
RNG  = np.random.RandomState(42)
N_BOOT = 5000

long = pd.read_csv(BASE/'09_integration/paired_delta/paired_feature_long.tsv', sep='\t')
# compute per-subject delta
long['delta'] = long['post'] - long['pre']
# Keep subjects with both pre and post (non-NaN delta)
long = long.dropna(subset=['delta'])

def bca_ci(samples, stat_fn, n_boot=N_BOOT, alpha=0.05):
    samples = np.asarray(samples, float)
    samples = samples[~np.isnan(samples)]
    if len(samples) < 3: return (np.nan,)*3
    theta = stat_fn(samples); n = len(samples)
    boots = np.array([stat_fn(samples[RNG.randint(0, n, n)]) for _ in range(n_boot)])
    boots = boots[~np.isnan(boots)]
    if len(boots) < 100: return theta, np.nan, np.nan
    z0 = st.norm.ppf(np.clip((boots < theta).mean(), 1/len(boots), 1 - 1/len(boots)))
    jk = np.array([stat_fn(np.delete(samples, i)) for i in range(n)])
    num = ((jk.mean() - jk)**3).sum()
    den = 6 * (((jk.mean() - jk)**2).sum())**1.5
    a = num/den if den > 0 else 0
    za = st.norm.ppf(alpha/2); zb = st.norm.ppf(1 - alpha/2)
    a1 = np.clip(st.norm.cdf(z0 + (z0 + za)/(1 - a*(z0+za))), 0, 1)
    a2 = np.clip(st.norm.cdf(z0 + (z0 + zb)/(1 - a*(z0+zb))), 0, 1)
    return theta, np.quantile(boots, a1), np.quantile(boots, a2)

rows = []
for feat in sorted(long.feature.unique()):
    sub = long[long.feature == feat]
    g = sub[sub.response=='good']['delta'].values
    b = sub[sub.response=='bad']['delta'].values
    # median deltas per group + BCa
    med_g, lo_g, hi_g = bca_ci(g, np.median)
    med_b, lo_b, hi_b = bca_ci(b, np.median)
    # between-group: bootstrap median(good) - median(bad)
    diff_samples = []
    if len(g) >= 3 and len(b) >= 3:
        for _ in range(N_BOOT):
            gs = g[RNG.randint(0, len(g), len(g))]
            bs = b[RNG.randint(0, len(b), len(b))]
            diff_samples.append(np.median(gs) - np.median(bs))
        diff_samples = np.array(diff_samples)
        d_med = np.median(diff_samples)
        d_lo, d_hi = np.quantile(diff_samples, [0.025, 0.975])
        # MW test (for reference)
        mw_p = st.mannwhitneyu(g, b, alternative='two-sided').pvalue
    else:
        d_med = d_lo = d_hi = mw_p = np.nan
    rows.append({
        'feature': feat,
        'n_good': len(g), 'n_bad': len(b),
        'good_delta_median': round(med_g, 3) if not np.isnan(med_g) else np.nan,
        'good_95CI': f'[{lo_g:+.2f}, {hi_g:+.2f}]' if not np.isnan(lo_g) else '',
        'bad_delta_median':  round(med_b, 3) if not np.isnan(med_b) else np.nan,
        'bad_95CI':  f'[{lo_b:+.2f}, {hi_b:+.2f}]' if not np.isnan(lo_b) else '',
        'diff_median_good_minus_bad': round(d_med, 3) if not np.isnan(d_med) else np.nan,
        'diff_95CI': f'[{d_lo:+.2f}, {d_hi:+.2f}]' if not np.isnan(d_lo) else '',
        'diff_CI_excludes_zero': int(not np.isnan(d_lo) and ((d_lo > 0) or (d_hi < 0))),
        'MW_p': round(mw_p, 3) if not np.isnan(mw_p) else np.nan,
        'evidence_level': 'robust (CI excl. 0)' if (not np.isnan(d_lo) and ((d_lo>0) or (d_hi<0)))
                          else ('exploratory (CI spans 0)' if len(g)>=3 and len(b)>=3 else 'insufficient_n'),
    })

df = pd.DataFrame(rows)
df.to_csv(TBL/'TableS8_cascade_BCa_bootstrap.tsv', sep='\t', index=False)

print('=== Cascade BCa 95% bootstrap CIs (n=14 paired, 7 good vs 7 bad) ===\n')
cols = ['feature','n_good','n_bad','good_delta_median','good_95CI','bad_delta_median','bad_95CI',
        'diff_median_good_minus_bad','diff_95CI','diff_CI_excludes_zero','MW_p','evidence_level']
print(df[cols].to_string(index=False))
print(f'\nSaved: {TBL}/TableS8_cascade_BCa_bootstrap.tsv')
