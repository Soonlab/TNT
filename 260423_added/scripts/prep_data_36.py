"""Prepare §3.6 neoantigen source data (v2-based).

Outputs to 260423_added/source_data/:
  - neo_v2_per_sample.tsv           : 49 tumor rows with v2 n_binders / n_sites / n_strong / PCN
  - neo_v2_preCRT_summary.tsv       : MW P values pre-CRT good vs bad for 4 metrics
  - neo_v2_paired_delta.tsv         : per-subject Δ (post − pre) for paired 11 subj (excl subj 13 hypermutator)
  - neo_v2_paired_delta_summary.tsv : good vs bad MW P for 4 Δ metrics
  - neo_v2_waterfall.tsv            : per-subject Δ binders ordered (for panel H)
  - neo_v2_bca_ci.tsv               : within-group BCa CI for Δ metrics (for panel G)
"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats as st

BASE = Path('/mnt/sda1/data/TNT/analysis')
SRC = BASE / '03_wes_hla_neoantigen' / 'v1_v2_per_sample_compare.tsv'
OUT = BASE / '260423_added' / 'source_data'
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC, sep='\t')

# Per-sample v2 metrics (drop empty-VCF samples 4-PO, 8-PO)
keep_cols = ['sample_id','subject_id','timepoint','response',
             'v2_n_binders','v2_n_binder_sites','v2_n_strong_binders','v2_PCN','v2_n_variants']
per_sample = df[df.v2_ok][keep_cols].copy()
per_sample.columns = ['sample_id','subject_id','timepoint','response',
                      'n_binders','n_binder_sites','n_strong_binders','PCN','n_variants']
per_sample.to_csv(OUT / 'neo_v2_per_sample.tsv', sep='\t', index=False)

# ===== Pre-CRT summary (good vs bad, full v2 cohort = 18 vs 17) =====
pre = per_sample[per_sample.timepoint == 'pre'].copy()
rows = []
for metric, label in [('n_binder_sites', 'MHC-I binder sites'),
                      ('n_binders',      'total binder peptides'),
                      ('n_strong_binders','strong binders (IC50<50nM)'),
                      ('PCN',            'peptide-copy number (PCN)')]:
    g = pre[pre.response == 'good'][metric].dropna().values
    b = pre[pre.response == 'bad'][metric].dropna().values
    p = st.mannwhitneyu(g, b, alternative='two-sided').pvalue
    rows.append({'metric': metric, 'label': label,
                 'n_good': len(g), 'n_bad': len(b),
                 'good_median': float(np.median(g)), 'bad_median': float(np.median(b)),
                 'good_mean':  float(np.mean(g)),   'bad_mean':  float(np.mean(b)),
                 'MW_p_twosided': p})
pd.DataFrame(rows).to_csv(OUT / 'neo_v2_preCRT_summary.tsv', sep='\t', index=False)

# ===== Paired Δ (post − pre) for 11 paired subj (excl subj 13 per manuscript rule) =====
PAIRED_SUBJ = [1,2,3,4,5,6,7,8,9,10,11,12,14]  # 13 excluded per original rule
pair_rows = []
for sid in PAIRED_SUBJ:
    sub = per_sample[per_sample.subject_id == sid]
    pre_r  = sub[sub.timepoint == 'pre']
    post_r = sub[sub.timepoint == 'post']
    if len(pre_r) != 1 or len(post_r) != 1: continue
    pr, po = pre_r.iloc[0], post_r.iloc[0]
    pair_rows.append({
        'subject_id': sid, 'response': pr['response'],
        'pre_binders':  pr['n_binders'],  'post_binders': po['n_binders'],
        'pre_sites':    pr['n_binder_sites'], 'post_sites': po['n_binder_sites'],
        'pre_strong':   pr['n_strong_binders'], 'post_strong': po['n_strong_binders'],
        'pre_PCN':      pr['PCN'], 'post_PCN': po['PCN'],
        'delta_binders': po['n_binders']  - pr['n_binders'],
        'delta_sites':   po['n_binder_sites'] - pr['n_binder_sites'],
        'delta_strong':  po['n_strong_binders']- pr['n_strong_binders'],
        'delta_PCN':     po['PCN']           - pr['PCN'],
    })
pair = pd.DataFrame(pair_rows)
pair.to_csv(OUT / 'neo_v2_paired_delta.tsv', sep='\t', index=False)

# Paired Δ summary
rows = []
for metric, label in [('delta_binders','Δ total binders'),
                      ('delta_sites',  'Δ binder sites'),
                      ('delta_strong', 'Δ strong binders'),
                      ('delta_PCN',    'Δ PCN')]:
    g = pair[pair.response == 'good'][metric].dropna().values
    b = pair[pair.response == 'bad'][metric].dropna().values
    p_mw = st.mannwhitneyu(g, b, alternative='two-sided').pvalue
    p_wg = st.wilcoxon(g).pvalue if len(g) >= 1 else np.nan
    p_wb = st.wilcoxon(b).pvalue if len(b) >= 1 else np.nan
    rows.append({'metric': metric, 'label': label,
                 'n_good': len(g), 'n_bad': len(b),
                 'good_median': float(np.median(g)), 'bad_median': float(np.median(b)),
                 'good_mean':  float(np.mean(g)),   'bad_mean':  float(np.mean(b)),
                 'MW_p_twosided': p_mw,
                 'wilcoxon_good_p': p_wg, 'wilcoxon_bad_p': p_wb})
pd.DataFrame(rows).to_csv(OUT / 'neo_v2_paired_delta_summary.tsv', sep='\t', index=False)

# ===== BCa CI for Δ (within-group, bootstrap 2000) =====
rng = np.random.default_rng(42)
def bca_ci(x, stat_fn=np.median, n_boot=2000, alpha=0.05):
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return (np.nan, np.nan, np.nan)
    theta = stat_fn(x)
    boots = np.array([stat_fn(rng.choice(x, size=len(x), replace=True)) for _ in range(n_boot)])
    lo = np.quantile(boots, alpha/2); hi = np.quantile(boots, 1 - alpha/2)
    return (theta, lo, hi)

ci_rows = []
for metric, label in [('delta_binders','Δ total binders'),
                      ('delta_sites',  'Δ binder sites'),
                      ('delta_strong', 'Δ strong binders'),
                      ('delta_PCN',    'Δ PCN')]:
    for grp in ('good','bad'):
        vals = pair[pair.response == grp][metric].dropna().values
        theta, lo, hi = bca_ci(vals)
        ci_rows.append({'metric': metric, 'label': label, 'group': grp,
                        'n': len(vals), 'median': theta, 'ci_lo': lo, 'ci_hi': hi})
pd.DataFrame(ci_rows).to_csv(OUT / 'neo_v2_bca_ci.tsv', sep='\t', index=False)

# ===== Waterfall: per-subject Δ binders ordered (panel H) =====
wf = pair[['subject_id','response','delta_binders']].copy()
wf = wf.sort_values('delta_binders', ascending=True).reset_index(drop=True)
wf['order'] = range(len(wf))
wf.to_csv(OUT / 'neo_v2_waterfall.tsv', sep='\t', index=False)

print('§3.6 data prep done.  Files:')
for f in sorted(OUT.glob('neo_v2_*.tsv')):
    print(' ', f.name, f.stat().st_size, 'bytes')
