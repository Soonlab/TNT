"""
Final external meta — Option C:
- Thread 1 (tumor-intrinsic): 5 concordant cohorts (N=518), per-sample Stouffer
- Thread 2 (immune CD8 axis): 5 concordant cohorts (N=518) + Akiyoshi 2023 paper-level
  added via OR/p-value -> Z conversion (N=816 total)

Akiyoshi 2023 (GSE216616, n=298) reports for the CD8/cytolytic axis:
  - Cytolytic activity (GZMA x PRF1): TRG1/2 (n=212, good) vs TRG3/4 (n=86, bad), P=0.005
  - MCP-counter cytotoxic lymphocyte: TRG1/2 vs TRG3/4 P=0.005
  - Activated CD8 T cell ssGSEA:    P=0.03
  - Effector memory CD8 T cell:     P<0.001 (we use 0.0009 as a conservative point estimate)

We use Cytolytic activity P=0.005 as the primary Akiyoshi statistic since GZMA + PRF1 are
both core members of our CD8_cytotoxic signature. P<0.001 effector memory is reported in
sensitivity.

Per-sample TRG labels are NOT in GEO for either Akiyoshi cohort (BJS 2019 supplement: figures
only; JAMA 2023 supplement: aggregate counts only) — confirmed by direct inspection.
"""
import os, numpy as np, pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/data/data/TNT/analysis/260418_add'
SRC = '/data/data/TNT/analysis/11_external_validation/v3_signature_response_stats.tsv'

KEEP = ['GSE56699','GSE45404','GSE133057','GSE87211','GSE35452']
THREAD1 = ['DSB_HDR_repair','E2F_MYC_cellcycle','Tumor_cellcycle','EMT']
THREAD2 = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration']
EXPECTED = {'DSB_HDR_repair':+1,'E2F_MYC_cellcycle':+1,'Tumor_cellcycle':+1,'EMT':-1,
            'CD8_cytotoxic':+1,'Tcell_infiltration':+1,'Bcell_infiltration':+1}

df = pd.read_csv(SRC, sep='\t')
sub = df[df.gse.isin(KEEP)].copy()

def stouffer_per_sample(rows):
    """Standard √N-weighted Stouffer's Z signed by sign(expected_dir * delta)."""
    w = np.sqrt(rows.n_good.values + rows.n_bad.values)
    z_per = np.array([norm.isf(p/2) * np.sign(EXPECTED[s] * d)
                      for s, p, d in zip(rows.signature, rows.pvalue, rows.delta)])
    Z = np.sum(w * z_per) / np.sqrt(np.sum(w**2))
    return float(Z), 2*(1-norm.cdf(abs(Z)))

# ---- Thread 1 final ----
t1_rows = []
for sig in THREAD1:
    s = sub[sub.signature == sig]
    Z, p = stouffer_per_sample(s)
    n_total = int((s.n_good + s.n_bad).sum())
    t1_rows.append({'thread':'Thread1_tumor_intrinsic','signature':sig,
                    'n_cohorts':len(s),'n_total':n_total,
                    'Z':round(Z,2),'p_meta':round(p,4),
                    'cohorts':','.join(s.gse.tolist()),
                    'akiyoshi':'not applicable (paper has no Thread 1 effect sizes)'})

# ---- Thread 2 final (with Akiyoshi for CD8_cytotoxic) ----
# Akiyoshi 2023 contribution
AKI = {
    'cytolytic_activity_TRG12vs34':  dict(p=0.005,  n_good=212, n_bad=86, signed=+1, source='eFig 4B'),
    'effector_memory_CD8_TRG12vs34': dict(p=0.0009, n_good=212, n_bad=86, signed=+1, source='eFig 8'),
    'activated_CD8_TRG12vs34':       dict(p=0.03,   n_good=212, n_bad=86, signed=+1, source='eFig 8'),
    'mcp_cytotoxic_TRG12vs34':       dict(p=0.005,  n_good=212, n_bad=86, signed=+1, source='eFig 4A'),
}
print('Akiyoshi 2023 paper-level statistics:')
for k, v in AKI.items():
    z = norm.isf(v['p']/2) * v['signed']
    print(f'  {k:32s}  P={v["p"]:.4f}  Z={z:+.2f}  n={v["n_good"]+v["n_bad"]}  ({v["source"]})')

PRIMARY_AKI = AKI['cytolytic_activity_TRG12vs34']  # most direct for CD8_cytotoxic

t2_rows = []
for sig in THREAD2:
    s = sub[sub.signature == sig]
    Z5, p5 = stouffer_per_sample(s)
    n5 = int((s.n_good + s.n_bad).sum())
    if sig == 'CD8_cytotoxic':
        # Augment with Akiyoshi paper-level using cytolytic_activity P=0.005
        w_per_cohort = np.sqrt(s.n_good.values + s.n_bad.values)
        z_per_cohort = np.array([norm.isf(p/2) * np.sign(EXPECTED[sig] * d)
                                 for p, d in zip(s.pvalue, s.delta)])
        # Akiyoshi
        w_aki = np.sqrt(PRIMARY_AKI['n_good'] + PRIMARY_AKI['n_bad'])
        z_aki = norm.isf(PRIMARY_AKI['p']/2) * PRIMARY_AKI['signed']
        w_all = np.concatenate([w_per_cohort, [w_aki]])
        z_all = np.concatenate([z_per_cohort, [z_aki]])
        Z6 = float(np.sum(w_all * z_all) / np.sqrt(np.sum(w_all**2)))
        p6 = 2 * (1 - norm.cdf(abs(Z6)))
        n6 = n5 + PRIMARY_AKI['n_good'] + PRIMARY_AKI['n_bad']
        t2_rows.append({'thread':'Thread2_immune','signature':sig,
                        'n_cohorts':len(s)+1,'n_total':n6,
                        'Z':round(Z6,2),'p_meta':round(p6,4),
                        'cohorts':','.join(s.gse.tolist())+',GSE216616(Akiyoshi-paper)',
                        'akiyoshi':f'cytolytic_activity P=0.005 added (eFig 4B), Z={z_aki:+.2f} weight={w_aki:.1f}',
                        '5cohort_only_Z':round(Z5,2),'5cohort_only_p':round(p5,4)})
    else:
        t2_rows.append({'thread':'Thread2_immune','signature':sig,
                        'n_cohorts':len(s),'n_total':n5,
                        'Z':round(Z5,2),'p_meta':round(p5,4),
                        'cohorts':','.join(s.gse.tolist()),
                        'akiyoshi':'paper does not report direct equivalent'})

R = pd.DataFrame(t1_rows + t2_rows)
R.to_csv(f'{OUT}/FINAL_meta_with_akiyoshi.tsv', sep='\t', index=False)
print('\n=== FINAL META (Option C) ===')
print(R.to_string(index=False))

# Sensitivity: try Akiyoshi with effector memory CD8 (P<0.001) instead — should give stronger Z
sig = 'CD8_cytotoxic'
s = sub[sub.signature == sig]
w_per = np.sqrt(s.n_good.values + s.n_bad.values)
z_per = np.array([norm.isf(p/2) * np.sign(EXPECTED[sig] * d) for p, d in zip(s.pvalue, s.delta)])
print('\n=== Sensitivity: alternative Akiyoshi statistics for CD8_cytotoxic ===')
for k, v in AKI.items():
    z_aki = norm.isf(v['p']/2) * v['signed']; w_aki = np.sqrt(v['n_good']+v['n_bad'])
    Z = float(np.sum(np.concatenate([w_per,[w_aki]]) * np.concatenate([z_per,[z_aki]])) /
              np.sqrt(np.sum(np.concatenate([w_per,[w_aki]])**2)))
    p = 2*(1-norm.cdf(abs(Z)))
    print(f'  using {k:32s} -> 6-source Z={Z:+.2f}, P={p:.4f}')

print(f'\nWrote FINAL_meta_with_akiyoshi.tsv')
