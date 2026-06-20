"""Compare pVACseq v1 vs v2 neoantigen outputs — table only (no figures).

v1 output dir: 03_wes_hla_neoantigen/pvacseq/{sample}/MHC_Class_I/{sample}_DNA.MHC_I.all_epitopes.tsv
v2 output dir: 03_wes_hla_neoantigen/pvacseq_v2/{sample}/MHC_Class_I/{sample}_DNA.MHC_I.all_epitopes.tsv

Per-sample metrics (IC50 binder: Best MT IC50 Score < 500 nM):
  - n_binders           : total peptide-HLA predictions with IC50 < 500 nM
  - n_binder_sites      : unique (Chrom, Start, Stop, Variant) mutation sites with >=1 binder
  - n_strong_binders    : peptide-HLA predictions with IC50 < 50 nM
  - PCN                 : peptide-copy number = sum over mutation sites of (VAF * binder_count_at_site) approx — use tumor DNA VAF mean per site * n_binders_at_site, aggregated
Note: v1 code used 'Tumor DNA VAF' column for PCN-style weighting; follow same.
"""
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path('/mnt/sda1/data/TNT/analysis')
V1 = BASE/'03_wes_hla_neoantigen/pvacseq'
V2 = BASE/'03_wes_hla_neoantigen/pvacseq_v2'
OUT_DIR = BASE/'03_wes_hla_neoantigen'
OUT_DIR.mkdir(exist_ok=True)

INV = pd.read_csv(BASE/'00_cohort/wes_inventory.tsv', sep='\t')
CLIN = pd.read_csv(BASE/'00_cohort/clinical_master.tsv', sep='\t')

def summarize(all_ep_tsv):
    if not all_ep_tsv.exists():
        return None
    df = pd.read_csv(all_ep_tsv, sep='\t', low_memory=False)
    if len(df) == 0:
        return dict(n_binders=0, n_binder_sites=0, n_strong_binders=0, PCN=0.0, n_variants=0)
    ic50 = pd.to_numeric(df['Best MT IC50 Score'], errors='coerce')
    vaf  = pd.to_numeric(df['Tumor DNA VAF'], errors='coerce')
    binder = ic50 < 500
    strong = ic50 < 50
    site_key = df[['Chromosome','Start','Stop','Variant']].astype(str).agg('|'.join, axis=1)
    pcn_per_row = np.where(binder, vaf.fillna(0.0)*2.0, 0.0)  # 2*VAF as copy-number proxy
    return dict(
        n_binders       = int(binder.sum()),
        n_binder_sites  = int(site_key[binder].nunique()),
        n_strong_binders= int(strong.sum()),
        PCN             = float(np.nansum(pcn_per_row)),
        n_variants      = int(site_key.nunique()),
    )

rows = []
tumor_ids = INV[INV.timepoint!='normal'].sort_values('sample_id')['sample_id'].tolist()

for sid in tumor_ids:
    v1_path = V1/sid/'MHC_Class_I'/f'{sid}_DNA.MHC_I.all_epitopes.tsv'
    v2_path = V2/sid/'MHC_Class_I'/f'{sid}_DNA.MHC_I.all_epitopes.tsv'
    v1_sum = summarize(v1_path)
    v2_sum = summarize(v2_path)
    row = {'sample_id': sid, 'v1_ok': v1_sum is not None, 'v2_ok': v2_sum is not None}
    for k, d in [('v1', v1_sum), ('v2', v2_sum)]:
        if d is None:
            for m in ['n_binders','n_binder_sites','n_strong_binders','PCN','n_variants']:
                row[f'{k}_{m}'] = np.nan
        else:
            for m, v in d.items():
                row[f'{k}_{m}'] = v
    rows.append(row)

df = pd.DataFrame(rows)
# Add subject/response/timepoint
sid2sub = dict(zip(INV.sample_id, INV.subject_id))
sid2tp  = dict(zip(INV.sample_id, INV.timepoint))
sub2resp = dict(zip(CLIN.subject_id, CLIN.response_bin))
df['subject_id'] = df.sample_id.map(sid2sub)
df['timepoint']  = df.sample_id.map(sid2tp)
df['response']   = df.subject_id.map(sub2resp)

# Per-sample Δ (v2 − v1) and relative change
for m in ['n_binders','n_binder_sites','n_strong_binders','PCN']:
    df[f'delta_{m}']     = df[f'v2_{m}'] - df[f'v1_{m}']
    df[f'pctchg_{m}']    = 100 * df[f'delta_{m}'] / df[f'v1_{m}'].replace(0, np.nan)

df.to_csv(OUT_DIR/'v1_v2_per_sample_compare.tsv', sep='\t', index=False)

# ===== Paired subject-level Δ (post − pre) recomputation =====
PAIRED_SUBJ = [1,2,3,4,5,6,7,8,9,10,11,12,14]

pair_rows = []
for subj in PAIRED_SUBJ:
    sub_df = df[df.subject_id == subj]
    pre = sub_df[sub_df.timepoint == 'pre']
    post = sub_df[sub_df.timepoint == 'post']
    if len(pre) != 1 or len(post) != 1:
        continue
    pre_r = pre.iloc[0]; post_r = post.iloc[0]
    for ver in ['v1','v2']:
        pair_rows.append({
            'subject_id': subj,
            'response':   pre_r['response'],
            'version':    ver,
            'pre_binders':  pre_r[f'{ver}_n_binders'],
            'post_binders': post_r[f'{ver}_n_binders'],
            'delta_binders':  post_r[f'{ver}_n_binders']  - pre_r[f'{ver}_n_binders'],
            'delta_sites':    post_r[f'{ver}_n_binder_sites'] - pre_r[f'{ver}_n_binder_sites'],
            'delta_strong':   post_r[f'{ver}_n_strong_binders']- pre_r[f'{ver}_n_strong_binders'],
            'delta_PCN':      post_r[f'{ver}_PCN']           - pre_r[f'{ver}_PCN'],
        })
pair = pd.DataFrame(pair_rows)
pair.to_csv(OUT_DIR/'v1_v2_paired_delta_compare.tsv', sep='\t', index=False)

# ===== Summary table for manuscript section §3.6 =====
# Recompute: n_binder_sites_preCRT MW (good vs bad, all pre tumors incl unmatched)
# and within-good / within-bad paired Δ medians with BCa-like CI (n-based CI skip here — just medians)
from scipy import stats as st

def preCRT_mwu(df, col_key):
    pre = df[df.timepoint == 'pre'].dropna(subset=[col_key])
    g = pre[pre.response == 'good'][col_key].values
    b = pre[pre.response == 'bad'][col_key].values
    if len(g) < 2 or len(b) < 2: return (np.nan,)*4
    p = st.mannwhitneyu(g, b, alternative='two-sided').pvalue
    return (np.median(g), np.median(b), p, f'{len(g)}v{len(b)}')

rows = []
for metric_key, label in [('n_binder_sites','sites with ≥1 MHC-I binder'),
                          ('n_binders','total binder peptides'),
                          ('PCN','peptide-copy number (PCN proxy)')]:
    for ver in ['v1','v2']:
        col = f'{ver}_{metric_key}'
        g_med, b_med, p, n = preCRT_mwu(df, col)
        rows.append({
            'metric': label, 'version': ver,
            'good_median_pre': g_med, 'bad_median_pre': b_med,
            'MW_p_pre': p, 'n_pre': n,
        })
pre_tbl = pd.DataFrame(rows)
pre_tbl.to_csv(OUT_DIR/'v1_v2_preCRT_summary.tsv', sep='\t', index=False)

# Paired within-group medians (n=11 paired, excluding subj 13 since it's tumor-only)
pair_tbl_rows = []
for ver in ['v1','v2']:
    sub = pair[pair.version == ver]
    sub = sub[sub.subject_id != 13]  # neoantigen-style exclusion
    g = sub[sub.response == 'good']
    b = sub[sub.response == 'bad']
    for metric in ['delta_binders','delta_sites','delta_strong','delta_PCN']:
        g_vals = g[metric].dropna().values
        b_vals = b[metric].dropna().values
        if len(g_vals) < 2 or len(b_vals) < 2: continue
        p_mw = st.mannwhitneyu(g_vals, b_vals, alternative='two-sided').pvalue
        pair_tbl_rows.append({
            'version': ver, 'metric': metric,
            'n_good': len(g_vals), 'n_bad': len(b_vals),
            'good_median': np.median(g_vals),
            'bad_median':  np.median(b_vals),
            'good_mean':   np.mean(g_vals),
            'bad_mean':    np.mean(b_vals),
            'MW_p': p_mw,
        })
pair_tbl = pd.DataFrame(pair_tbl_rows)
pair_tbl.to_csv(OUT_DIR/'v1_v2_paired_summary.tsv', sep='\t', index=False)

# ===== Print tables =====
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

print('\n============== Per-sample v1 vs v2 completion status ==============')
print(df[['sample_id','response','timepoint','v1_ok','v2_ok',
         'v1_n_binders','v2_n_binders','delta_n_binders','pctchg_n_binders']].to_string(index=False))

print('\n============== §3.6 Pre-CRT summary (good vs bad, all pre tumors) ==============')
print(pre_tbl.to_string(index=False))

print('\n============== §3.6 Paired Δ summary (n=11 excl. subj 13) ==============')
print(pair_tbl.to_string(index=False))

print('\n============== Per-subject paired Δ binders (v1 vs v2) ==============')
pv = pair.pivot_table(index=['subject_id','response'], columns='version',
                      values=['delta_binders','delta_sites'], aggfunc='first')
print(pv.to_string())
print('\nSaved:')
print('  - v1_v2_per_sample_compare.tsv')
print('  - v1_v2_paired_delta_compare.tsv')
print('  - v1_v2_preCRT_summary.tsv')
print('  - v1_v2_paired_summary.tsv')
