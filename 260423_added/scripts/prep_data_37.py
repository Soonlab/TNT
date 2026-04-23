"""Prepare §3.7 HLA class I typing + HLA-LOH source data.

Outputs to 260423_added/source_data/:
  - hla_homozygosity.tsv       : per-subject homozygosity flags × response
  - hla_loh_subject_summary.tsv: per-subject any-LOH (strict + lite) × response
  - hla_loh_per_locus.tsv      : per-locus LOH frequency good vs bad (A/B/C)
  - hla_loh_clearance.tsv      : paired pre→post LOH resolution trajectory (strict)
  - hla_loh_concordance.tsv    : per-locus strict vs lite concordance
  - hla_loh_tests.tsv          : Fisher contingency & P for strict/lite × any-LOH
"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats as st

BASE = Path('/mnt/sda1/data/TNT/analysis')
OUT  = BASE / '260423_added' / 'source_data'
OUT.mkdir(parents=True, exist_ok=True)

hla  = pd.read_csv(BASE / '03_hla' / 'hla_class_I_typing.tsv', sep='\t')
loh  = pd.read_csv(BASE / '03_hla' / 'loh_stricter' / 'hla_loh_per_locus_strict.tsv', sep='\t')
clin = pd.read_csv(BASE / '00_cohort' / 'clinical_master.tsv', sep='\t')

# ===== 1. Homozygosity =====
hom = hla[['subject_id','homozygous_A','homozygous_B','homozygous_C',
           'n_homozygous_loci','response_bin']].copy()
hom['any_homozygous'] = hom[['homozygous_A','homozygous_B','homozygous_C']].any(axis=1)
hom.to_csv(OUT / 'hla_homozygosity.tsv', sep='\t', index=False)

# Homozygosity Fisher per-locus good vs bad
rows = []
for locus in ['A','B','C']:
    col = f'homozygous_{locus}'
    g_pos = int(hom[(hom.response_bin=='good') & hom[col]].shape[0])
    g_neg = int(hom[(hom.response_bin=='good') & ~hom[col]].shape[0])
    b_pos = int(hom[(hom.response_bin=='bad') & hom[col]].shape[0])
    b_neg = int(hom[(hom.response_bin=='bad') & ~hom[col]].shape[0])
    odds, p = st.fisher_exact([[g_pos, g_neg], [b_pos, b_neg]])
    rows.append({'metric': f'homozygous_{locus}', 'good_pos': g_pos, 'good_neg': g_neg,
                 'bad_pos': b_pos, 'bad_neg': b_neg,
                 'good_freq': g_pos/(g_pos+g_neg), 'bad_freq': b_pos/(b_pos+b_neg),
                 'OR': odds, 'fisher_p': p})
# any
g_pos = int(hom[(hom.response_bin=='good') & hom.any_homozygous].shape[0])
g_neg = int(hom[(hom.response_bin=='good') & ~hom.any_homozygous].shape[0])
b_pos = int(hom[(hom.response_bin=='bad') & hom.any_homozygous].shape[0])
b_neg = int(hom[(hom.response_bin=='bad') & ~hom.any_homozygous].shape[0])
odds, p = st.fisher_exact([[g_pos, g_neg], [b_pos, b_neg]])
rows.append({'metric': 'any_homozygous', 'good_pos': g_pos, 'good_neg': g_neg,
             'bad_pos': b_pos, 'bad_neg': b_neg,
             'good_freq': g_pos/(g_pos+g_neg), 'bad_freq': b_pos/(b_pos+b_neg),
             'OR': odds, 'fisher_p': p})
pd.DataFrame(rows).to_csv(OUT / 'hla_homozygosity_tests.tsv', sep='\t', index=False)

# ===== 2. Per-subject LOH summary (pre-CRT tumors) =====
# Identify pre samples by sample-id convention
loh['timepoint'] = loh['sample'].apply(
    lambda s: 'pre' if (s.endswith('-PR') or s.endswith('-P')) else
              ('post' if s.endswith('-PO') else 'other'))
pre_loh = loh[loh.timepoint == 'pre'].copy()

subj_summary = (pre_loh.groupby('subject_id')
                .agg(n_het_loci=('is_het_normal','sum'),
                     n_loh_lite=('loh_lite','sum'),
                     n_loh_strict=('loh_strict','sum'))
                .reset_index())
subj_summary = subj_summary.merge(clin[['subject_id','response_bin']], on='subject_id', how='left')
subj_summary['any_loh_lite'] = subj_summary.n_loh_lite > 0
subj_summary['any_loh_strict'] = subj_summary.n_loh_strict > 0
subj_summary.to_csv(OUT / 'hla_loh_subject_summary.tsv', sep='\t', index=False)

# Fisher tests
tests = []
for call_col, label in [('any_loh_strict','strict any-LOH'),
                        ('any_loh_lite','lite any-LOH')]:
    g_pos = int(subj_summary[(subj_summary.response_bin=='good') & subj_summary[call_col]].shape[0])
    g_neg = int(subj_summary[(subj_summary.response_bin=='good') & ~subj_summary[call_col]].shape[0])
    b_pos = int(subj_summary[(subj_summary.response_bin=='bad') & subj_summary[call_col]].shape[0])
    b_neg = int(subj_summary[(subj_summary.response_bin=='bad') & ~subj_summary[call_col]].shape[0])
    odds, p = st.fisher_exact([[g_pos, g_neg], [b_pos, b_neg]])
    tests.append({'call': label, 'good_pos': g_pos, 'good_neg': g_neg,
                  'bad_pos': b_pos, 'bad_neg': b_neg,
                  'good_freq': g_pos/(g_pos+g_neg) if (g_pos+g_neg)>0 else np.nan,
                  'bad_freq':  b_pos/(b_pos+b_neg) if (b_pos+b_neg)>0 else np.nan,
                  'OR': odds, 'fisher_p': p})
pd.DataFrame(tests).to_csv(OUT / 'hla_loh_tests.tsv', sep='\t', index=False)

# ===== 3. Per-locus LOH frequency (pre-CRT, by response) =====
loc_rows = []
for locus in ['HLA-A','HLA-B','HLA-C']:
    loc_df = pre_loh[(pre_loh.locus == locus) & (pre_loh.is_het_normal)]
    loc_df = loc_df.merge(clin[['subject_id','response_bin']], on='subject_id', how='left')
    for call_col, label in [('loh_strict','strict'), ('loh_lite','lite')]:
        g = loc_df[loc_df.response_bin == 'good']
        b = loc_df[loc_df.response_bin == 'bad']
        g_pos = int(g[call_col].sum()); g_n = len(g)
        b_pos = int(b[call_col].sum()); b_n = len(b)
        loc_rows.append({'locus': locus, 'call': label,
                         'good_pos': g_pos, 'good_n_het': g_n,
                         'bad_pos': b_pos, 'bad_n_het': b_n,
                         'good_freq': g_pos/g_n if g_n>0 else 0,
                         'bad_freq': b_pos/b_n if b_n>0 else 0})
pd.DataFrame(loc_rows).to_csv(OUT / 'hla_loh_per_locus.tsv', sep='\t', index=False)

# ===== 4. Clearance trajectory: paired pre→post subjects (strict LOH) =====
# For each paired subject (1-14), per-locus pre/post strict LOH call
clearance_rows = []
for sid in [1,2,3,4,5,6,7,8,9,10,11,12,13,14]:
    for locus in ['HLA-A','HLA-B','HLA-C']:
        sub = loh[(loh.subject_id == sid) & (loh.locus == locus)]
        pre = sub[sub.timepoint == 'pre']
        post = sub[sub.timepoint == 'post']
        if len(pre) != 1 or len(post) != 1: continue
        pr, po = pre.iloc[0], post.iloc[0]
        clearance_rows.append({
            'subject_id': sid, 'locus': locus,
            'pre_ratio': pr.tumor_ratio, 'post_ratio': po.tumor_ratio,
            'pre_imbalance': abs(pr.tumor_ratio - 0.5),
            'post_imbalance': abs(po.tumor_ratio - 0.5),
            'pre_strict': bool(pr.loh_strict), 'post_strict': bool(po.loh_strict),
            'pre_lite': bool(pr.loh_lite), 'post_lite': bool(po.loh_lite),
        })
clearance = pd.DataFrame(clearance_rows)
clearance = clearance.merge(clin[['subject_id','response_bin']], on='subject_id', how='left')
clearance.to_csv(OUT / 'hla_loh_clearance.tsv', sep='\t', index=False)

# ===== 5. Strict vs lite per-locus concordance =====
conc = pre_loh[pre_loh.is_het_normal][
    ['subject_id','sample','locus','delta_ratio','fisher_p','fisher_p_bonf',
     'loh_lite','loh_strict']].copy()
conc = conc.merge(clin[['subject_id','response_bin']], on='subject_id', how='left')
conc['concordance'] = conc.apply(
    lambda r: ('both' if r.loh_strict and r.loh_lite
               else ('lite_only' if r.loh_lite else 'neither')), axis=1)
conc.to_csv(OUT / 'hla_loh_concordance.tsv', sep='\t', index=False)

print('§3.7 data prep done.  Files:')
for f in sorted(OUT.glob('hla_*.tsv')):
    print(' ', f.name, f.stat().st_size, 'bytes')
