"""BCa bootstrap 95% CIs for small-n cascade claims.
Addresses reviewer Major #1: n=6 vs 6 paired-delta claims need uncertainty quantification.
"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats as st
import warnings; warnings.filterwarnings('ignore')

BASE = Path('/mnt/sda1/data/TNT/analysis')
PD   = BASE/'09_integration/paired_delta'
TBL  = BASE/'tables'; TBL.mkdir(exist_ok=True)
RNG  = np.random.RandomState(42)
N_BOOT = 5000

def bca_ci(samples, stat_fn, n_boot=N_BOOT, alpha=0.05):
    """BCa bootstrap CI for a statistic computed from an array-like of samples."""
    samples = np.asarray(samples)
    theta_hat = stat_fn(samples)
    # bootstrap resamples
    n = len(samples)
    boots = np.array([stat_fn(samples[RNG.randint(0, n, n)]) for _ in range(n_boot)])
    # bias correction
    z0 = st.norm.ppf((boots < theta_hat).mean())
    # acceleration via jackknife
    jk = np.array([stat_fn(np.delete(samples, i)) for i in range(n)])
    jk_mean = jk.mean()
    num = ((jk_mean - jk)**3).sum()
    den = 6 * (((jk_mean - jk)**2).sum())**1.5
    a = num / den if den > 0 else 0
    za = st.norm.ppf(alpha/2); zb = st.norm.ppf(1 - alpha/2)
    a1 = st.norm.cdf(z0 + (z0 + za)/(1 - a*(z0+za)))
    a2 = st.norm.cdf(z0 + (z0 + zb)/(1 - a*(z0+zb)))
    return theta_hat, np.quantile(boots, a1), np.quantile(boots, a2)

# ---- Paired delta (22 signatures + ssGSEA) ----
rows = []
for fname, label in [('delta_22sigs_response.tsv','22sigs'),
                     ('delta_ssgsea_response.tsv','ssGSEA'),
                     ('delta_sbs_response.tsv','SBS'),
                     ('delta_tmb_response.tsv','TMB'),
                     ('delta_trust4_response.tsv','TRUST4')]:
    f = PD/fname
    if not f.exists(): continue
    df = pd.read_csv(f, sep='\t')
    # recover per-subject deltas from paired_feature_long
    pass

# Use paired_feature_long for per-subject deltas
long = pd.read_csv(PD/'paired_feature_long.tsv', sep='\t')
print('Long columns:', list(long.columns)[:15])
print('Shape:', long.shape)
print('Sample:')
print(long.head(3))

# Compute BCa CIs for key cascade claims
key_features = ['Treg','MHC_II','CD8_exhaustion','IGH_n_prod','IGK_n_prod','IGL_n_prod',
                'SBS5','n_missense','binder_sites','binders_total']

results = []
for feat in key_features:
    sub = long[long.feature == feat].copy() if 'feature' in long.columns else pd.DataFrame()
    if len(sub) == 0:
        # try different structure
        continue
    # stratify good vs bad and compute deltas
    if 'response_bin' not in sub.columns or 'delta' not in sub.columns: continue
    for grp, d in sub.groupby('response_bin'):
        vals = d['delta'].dropna().values
        if len(vals) < 3: continue
        med, lo, hi = bca_ci(vals, np.median)
        results.append({'feature':feat,'group':grp,'n':len(vals),
                        'median_delta':round(med,3),
                        'BCa_lo':round(lo,3),'BCa_hi':round(hi,3),
                        'CI_excludes_zero': int((lo > 0) or (hi < 0))})

# Fallback: reload directly from per-delta TSVs with subject-level deltas if paired_feature_long
# does not have per-subject deltas
if not results:
    print('paired_feature_long not in long-format; reconstructing from 22sigs + ssgsea direct stats')
    # The delta_{22sigs,ssgsea}_response.tsv files already contain median deltas per group.
    # We'll bootstrap from the long-format extraction below.
    for fname, tag in [('delta_22sigs_response.tsv','22sig'),('delta_ssgsea_response.tsv','ssGSEA')]:
        df = pd.read_csv(PD/fname, sep='\t')
        print(f'{fname}:', df.columns.tolist()[:8])
        print(df.head(2))

pd.DataFrame(results).to_csv(TBL/'cascade_bootstrap_BCa_CIs.tsv', sep='\t', index=False)
print('\n=== Cascade BCa 95% CIs ===')
print(pd.DataFrame(results).to_string(index=False))
