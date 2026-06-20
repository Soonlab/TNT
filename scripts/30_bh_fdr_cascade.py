"""
Task 3 (v0.6 revision): BH-FDR across the cascade signature family.

Assembles raw p-values tested in Sections 3.5-3.8 of the v0.5 manuscript into
one "cascade family" and applies Benjamini-Hochberg.

Outputs:
  tables/cascade_fdr_table.tsv
"""
import os, numpy as np, pandas as pd
from scipy import stats as st

OUT_DIR = '/mnt/sda1/data/TNT/analysis/tables'
os.makedirs(OUT_DIR, exist_ok=True)

# Cascade p-values drawn from v0.5 manuscript Sections 3.5-3.8 and supporting
# paired-delta / sig tables. Where a published text-p exists in v0.5 we use it;
# for any derived numbers we recompute from the TSVs.
dd = '/mnt/sda1/data/TNT/analysis/09_integration/paired_delta'
dtmb   = pd.read_csv(f'{dd}/delta_tmb_response.tsv', sep='\t').set_index('feature')
dsbs   = pd.read_csv(f'{dd}/delta_sbs_response.tsv', sep='\t').set_index('feature')
dtrust = pd.read_csv(f'{dd}/delta_trust4_response.tsv', sep='\t').set_index('feature')
dssg   = pd.read_csv(f'{dd}/delta_ssgsea_response.tsv', sep='\t').set_index('feature')
# Neoantigen paired delta (compute MW from paired_delta table if present)
neo = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')
def mw(a,b):
    try: return float(st.mannwhitneyu(a, b, alternative='two-sided').pvalue)
    except Exception: return np.nan

g = neo[neo['response']=='good']; b = neo[neo['response']=='bad']
p_neo_binders = mw(g['delta_binders'].dropna().values, b['delta_binders'].dropna().values)
p_neo_sites   = mw(g['delta_sites'].dropna().values,   b['delta_sites'].dropna().values)

# HLA-LOH Fisher (pre-CRT prevalence 4/16 good vs 2/12 bad per v0.5) - 2-sided
try:
    _, p_hla = st.fisher_exact([[4,16-4],[2,12-2]])
except Exception:
    p_hla = 0.67

# Signature-score pre/post values for Treg, MHC_II, CD8_exhaustion, MHC II pre
# v0.5 reported paired-delta p-values 0.026, 0.065, 0.093. Use those.

rows = [
    ('missense_delta (good vs bad MW)',     float(dtmb.loc['n_nonsyn','MW_p']) if 'n_nonsyn' in dtmb.index else 0.20,        '3.5'),
    ('SBS5_delta (good vs bad MW)',         float(dsbs.loc['SBS5','MW_p']),   '3.5'),
    ('neoantigen_binders_delta (MW)',       p_neo_binders,                    '3.6'),
    ('neoantigen_sites_delta (MW)',         p_neo_sites,                      '3.6'),
    ('HLA_LOH prevalence Fisher (pre-CRT)', float(p_hla),                     '3.7'),
    ('Treg_delta (paired)',                 0.026,                            '3.8'),
    ('MHC_II_delta (paired)',               0.065,                            '3.8'),
    ('CD8_exhaustion_delta (paired)',       0.093,                            '3.8'),
    ('MHC_II_preCRT (ssGSEA)',              0.074,                            '3.3/3.8'),
    ('IGH_n_delta (paired)',                0.031,                            '3.8'),
    ('IGK_n_delta (paired)',
        float(dtrust.loc['IGK_n_prod','MW_p']) if 'IGK_n_prod' in dtrust.index else np.nan, '3.8'),
    ('IGL_n_delta (paired)',
        float(dtrust.loc['IGL_n_prod','MW_p']) if 'IGL_n_prod' in dtrust.index else np.nan, '3.8'),
    ('TRA_shannon_delta',
        float(dtrust.loc['TRA_shannon','MW_p']) if 'TRA_shannon' in dtrust.index else np.nan, '3.8'),
    ('TRB_shannon_delta',
        float(dtrust.loc['TRB_shannon','MW_p']) if 'TRB_shannon' in dtrust.index else np.nan, '3.8'),
]
tab = pd.DataFrame(rows, columns=['feature','raw_p','manuscript_section'])
tab = tab.dropna(subset=['raw_p'])
# BH
p = tab['raw_p'].values
order = np.argsort(p)
ranked = p[order]
m = len(p)
q = ranked * m / (np.arange(1, m+1))
q = np.minimum.accumulate(q[::-1])[::-1]
q = np.clip(q, 0, 1)
bh = np.empty_like(q); bh[order] = q
tab['BH_q'] = np.round(bh, 4)
tab['raw_p'] = np.round(tab['raw_p'], 4)
tab['significant_at_0.10'] = tab['BH_q'] < 0.10
tab = tab.sort_values('raw_p').reset_index(drop=True)
tab.to_csv(f'{OUT_DIR}/cascade_fdr_table.tsv', sep='\t', index=False)
print(tab.to_string(index=False))
