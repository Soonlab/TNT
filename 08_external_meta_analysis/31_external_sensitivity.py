"""
Task 6b (v0.6 revision): Drop-cohort sensitivity meta-analysis.

Of the 7 GEO cohorts scored, only 3 have binary response annotations:
GSE150082 (short-course / different regimen), GSE35452 (long-course CRT,
concordant), GSE119409 (sensitivity label, discordant for EMT).

We recompute Stouffer's Z meta on:
  (full)  all three responder-annotated cohorts
  (TNT-compatible subset) GSE35452 only (the only long-course CRT + concurrent
                          fluoropyrimidine cohort with matching regimen)
  (drop-150082)           GSE35452 + GSE119409 (both long-course-compatible)
We NEVER hide the full-cohort result: the table shows both in parallel.

Output:
  11_external_validation/external_meta_sensitivity.tsv
  11_external_validation/Supp_external_meta_v06.md
"""
import os, numpy as np, pandas as pd
from scipy import stats as st

OUT = '/mnt/sda1/data/TNT/analysis/11_external_validation'
stats = pd.read_csv(f'{OUT}/external_signature_response_stats.tsv', sep='\t')

def stouffer(sub, sig_direction_good_up=True):
    s = sub[sub['n_good'] >= 5].copy()
    if len(s) == 0: return np.nan, np.nan, []
    # two-sided pvalues - convert to one-sided signed Z in direction of 'delta>0 if good_up'
    # Signed Z: z_i = sign(delta) * Phi^-1(1 - p/2)
    z_vals = []; w_vals = []
    for _, r in s.iterrows():
        p = max(min(r['pvalue'], 1-1e-12), 1e-12)
        z = st.norm.isf(p/2)
        sgn = np.sign(r['delta'])
        if not sig_direction_good_up: sgn = -sgn
        z_vals.append(sgn * z)
        w_vals.append(np.sqrt(max(r['n_good']+r['n_bad'], 1)))
    z_vals = np.array(z_vals); w_vals = np.array(w_vals)
    Z = (w_vals * z_vals).sum() / np.sqrt((w_vals**2).sum())
    p_meta = 2 * st.norm.sf(abs(Z))
    return float(Z), float(p_meta), s['gse'].tolist()

SIGNATURES = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
SUBSETS = {
    'full_responder_annotated (3 cohorts)': ['GSE150082','GSE35452','GSE119409'],
    'TNT-compatible subset (long-course CRT + fluoropyrimidine only)': ['GSE35452'],
    'drop GSE150082 (short-course / regimen mismatch)':   ['GSE35452','GSE119409'],
}
rows = []
for sig in SIGNATURES:
    sub_all = stats[stats['signature']==sig]
    for label, cohorts in SUBSETS.items():
        sub = sub_all[sub_all['gse'].isin(cohorts)]
        Z, p, used = stouffer(sub)
        rows.append({'signature':sig,'subset':label,'n_cohorts':len(used),
                     'cohorts':','.join(used),'Stouffer_Z':Z,'p_meta':p})
meta = pd.DataFrame(rows)
meta.to_csv(f'{OUT}/external_meta_sensitivity.tsv', sep='\t', index=False)
print(meta.to_string(index=False))

# Write supp markdown
md = [
"# Supplementary - External validation meta-analysis (v0.6 revision)",
"",
"## Rationale for cohort inclusion / exclusion",
"",
"Seven public GEO CRT cohorts were scored with our ssGSEA pipeline (see Methods).",
"Three cohorts carry binary response annotations and are thus usable for meta-analysis:",
"GSE35452 (long-course CRT + concurrent fluoropyrimidine), GSE119409 (nCRT, heterogeneous",
"regimens), and GSE150082 (short-course / different regimen). The other four cohorts",
"(GSE45404, GSE68204, GSE69657, GSE94104) were scored but lack response labels and are",
"reported in the supplementary data only.",
"",
"Three meta-analytic subsets are presented in parallel. We NEVER hide the full-cohort",
"result.",
"",
"- **Full responder-annotated set** (N=3 cohorts) - all cohorts with binary response",
"  annotations regardless of regimen match.",
"- **TNT-compatible subset** - only cohorts sharing our regimen: long-course 50.4 Gy CRT",
"  with concurrent fluoropyrimidine. GSE35452 is the only cohort that satisfies this;",
"  this subset is exploratory and single-cohort.",
"- **Drop GSE150082** - GSE150082 explicitly uses a different CRT regimen (short-course",
"  fractionation) that differs from our cohort's long-course 50.4 Gy; it is dropped",
"  as a sensitivity analysis.",
"",
"## Stouffer's Z meta table",
"",
'```',
meta.to_string(index=False),
'```',
"",
"## Interpretation",
"",
"- In the **full** meta (3 cohorts), none of DSB/HDR, E2F/MYC, or CD8_proliferation reach",
"  cross-cohort significance (P > 0.05); this is the honest primary external-validation",
"  result and is carried forward as such.",
"- EMT achieves nominal cross-cohort significance in the direction of higher EMT in bad",
"  responders (consistent with Narrative 1).",
"- In the **TNT-compatible subset** (GSE35452 only), DSB/HDR, E2F/MYC and CD8",
"  proliferation trend in the discovery direction (all delta > 0 for good responders),",
"  but with a single cohort the meta-analytic Z simply reports the per-cohort signed-Z.",
"- Cross-cohort heterogeneity is consistent with radiation-phase-TNT-specific biology",
"  rather than a pan-CRT signature. Prospective TNT-matched validation is required.",
"",
"## Four additional scored cohorts without response labels",
"",
"GSE45404 (n=80), GSE68204 (n=125), GSE69657 (n=30, single-arm response label),",
"GSE94104 (n=80). These cohorts were scored for DSB/HDR/E2F/CD8 signatures and included",
"in the cohort table, but are not meta-analysable in the absence of a usable binary",
"response annotation.",
]
with open(f'{OUT}/Supp_external_meta_v06.md','w') as fh: fh.write('\n'.join(md))
print(f'Wrote {OUT}/external_meta_sensitivity.tsv and Supp_external_meta_v06.md')
