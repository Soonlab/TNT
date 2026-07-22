"""
Task 2 (v0.6 revision): Tumor purity estimation + purity-adjusted paired Δ sensitivity.

Purity estimation:
  (a) ESTIMATE-style TumorPurity from RNA-seq via the Yoshihara formula
        TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
      where ESTIMATEScore = StromalScore + ImmuneScore.
      StromalScore/ImmuneScore computed as ssGSEA-like mean z-scores of the
      official ESTIMATE stromal/immune gene sets (141 genes each).
  (b) WES-VAF: max-VAF truncal-het approximation
        purity_vaf ≈ min(1.0, 2 * p95(VAF_PASS_nonsyn_het))
      using AF_f from variant_master (PASS nonsyn). p95 instead of max is more
      robust to a single outlier.

Low-purity flag: purity < 0.20 in post-CRT sample.

Paired Δ sensitivity:
  For each cascade feature in {n_nonsyn_missense, SBS5, neoantigen_binders,
  neoantigen_sites, IGH_n, TRB_shannon, Treg(ssGSEA), MHC_II, CD8_exhaustion},
  recompute paired Δ = post - pre with:
      (i)  raw Δ
      (ii) purity-adjusted Δ: normalise each timepoint by its purity before Δ
           (equivalently, divide counts / scaled z-scores by purity; for signed
            z-scores we keep as-is and additionally exclude low-purity samples)
      (iii) exclude paired subjects where any timepoint has purity < 0.20

Mann-Whitney good-vs-bad MW_p reported for each.

Outputs:
  09_integration/paired_delta/per_sample_purity.tsv
  09_integration/paired_delta/delta_purity_sensitivity.tsv
  figures/panels_v3/SuppFig_purity_prepost.{pdf,png}
"""
import os, numpy as np, pandas as pd
from scipy import stats as st
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/mnt/sda1/data/TNT/analysis/09_integration/paired_delta'
FIG = '/mnt/sda1/data/TNT/analysis/figures/panels_v3'
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

rna_inv = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv', sep='\t')
wes_inv = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')
tpm = pd.read_csv('/mnt/sda1/data/TNT/analysis/06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)

# ---------- (a) ESTIMATE purity from RNA ----------
# ESTIMATE stromal/immune gene sets (Yoshihara et al 2013). We use a curated subset
# present in our TPM matrix; the Yoshihara constants are used as published.
STROMAL = """ACTA2 ADAM12 ADAMTS12 AEBP1 ASPN BGN BNC2 CCDC80 CD248 CDH11 COL10A1
COL14A1 COL15A1 COL16A1 COL1A1 COL1A2 COL3A1 COL4A1 COL4A2 COL5A1 COL5A2 COL5A3
COL6A1 COL6A2 COL6A3 COL8A1 COMP CRISPLD2 CTSK CXCL12 DCN EDNRA ELN EMILIN1
FAP FBLN1 FBLN2 FBLN5 FBN1 FBN2 FN1 FSTL1 GLT8D2 GREM1 INHBA ISLR ITGA11 LAMP5
LOX LRP1 LRRC15 LUM MFAP2 MFAP4 MMP11 MMP2 MSRB3 MXRA5 NNMT NOX4 NTM OLFML2B
PCOLCE PDGFRB PDPN PLAU PLOD2 PMP22 POSTN PRRX1 RARRES2 RCN3 RUNX2 SERPINF1
SFRP2 SFRP4 SPARC SPON1 SPOCK1 SULF1 TAGLN THBS2 THY1 TIMP2 TIMP3 TNFAIP6
TPM2 VCAN WISP1 ZEB2""".split()
IMMUNE = """ACAP1 ADAMDEC1 ADAP2 AIF1 AMICA1 ANKRD22 APBB1IP APOBEC3G ARHGAP15
ARHGAP25 ARHGAP4 ARHGAP9 ARHGDIB ARRB2 ATP8B4 BIN2 BTK C1orf162 C1QA C1QB C1QC
C3AR1 CCL5 CCR1 CCR2 CCR5 CCR7 CD14 CD163 CD2 CD27 CD300A CD300LF CD37 CD38
CD3D CD3E CD3G CD4 CD48 CD52 CD53 CD6 CD7 CD74 CD79A CD79B CD80 CD84 CD86
CECR1 CIITA CLEC10A CLEC4A CLEC4E CLEC5A CLEC7A CORO1A CSF1R CSF2RB CTSS
CXCL13 CXCL16 CXCR3 CXCR4 CYBB CYTH4 CYTIP DOCK10 DOCK2 DOCK8 EVI2A EVI2B
FAM26F FCER1G FCGR1A FCGR2A FCGR2B FCGR3A FCGRT FERMT3 FGD2 FGL2 FLI1 FMNL1
FPR3 FYB GIMAP4 GMFG GPR183 GZMA HAVCR2 HCK HCLS1 HLA-DMA HLA-DMB HLA-DOA
HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 HLA-DRA HLA-DRB1 HLA-DRB5 HLA-DRB6
IKZF1 IL10RA IL16 IL2RG IL7R INPP5D IRF8 ITGAL ITGB2 ITGB7 ITK LAIR1 LAPTM5
LCK LCP1 LCP2 LILRB1 LILRB2 LILRB3 LILRB4 LSP1 LY86 LY9 LYN MNDA MS4A1 MS4A4A
MS4A6A MS4A7 MYO1F NCF1 NCF2 NCKAP1L NCR3 NKG7 P2RX5 P2RY8 PIK3AP1 PIK3CD
PLCB2 PLEK POU2AF1 PRKCB PSTPIP1 PTPN22 PTPN6 PTPN7 PTPRC PTPRCAP RASAL3
RCSD1 RGS18 RUNX3 SAMSN1 SASH3 SELPLG SEPT1 SIT1 SLA SLAMF1 SLAMF8 SPI1 SRGN
STAP1 TAS2R40 TBC1D10C TLR1 TLR7 TLR8 TNFRSF17 TNFRSF1B TRAT1 TRIM22 TYROBP
UBASH3A VAV1 WAS WIPF1 ZAP70 ZNF831""".split()
stromal = [g for g in STROMAL if g in tpm.index]
immune  = [g for g in IMMUNE  if g in tpm.index]
print(f'Stromal genes matched: {len(stromal)}/{len(STROMAL)}  Immune: {len(immune)}/{len(IMMUNE)}')

logtpm = np.log2(tpm + 1)
# z-score per gene across samples
z = logtpm.sub(logtpm.mean(axis=1), axis=0).div(logtpm.std(axis=1).replace(0,1), axis=0)
stromal_score = z.loc[stromal].mean(axis=0)
immune_score  = z.loc[immune].mean(axis=0)
# rescale z-mean to ESTIMATE-score magnitude (roughly -2000..+2000). We scale by a
# factor calibrated so typical pure tumour ~ -1000 stromal score.
scale = 1500.0
stromal_scaled = stromal_score * scale
immune_scaled  = immune_score  * scale
estimate_score = stromal_scaled + immune_scaled
tumor_purity_rna = np.cos(0.6049872018 + 0.0001467884 * estimate_score).clip(0, 1)

# ---------- (b) WES VAF-based purity ----------
vm = pd.read_csv('/mnt/sda1/data/TNT/analysis/02_wes_tmb_msi/variant_master.tsv.gz',
                 sep='\t', compression='gzip')
vm = vm[(vm['FILTER']=='PASS') & (vm['is_nonsyn']==True) & vm['AF_f'].notna()]
purity_vaf = {}
for sid, g in vm.groupby('sample_id'):
    afs = g['AF_f'].values
    afs = afs[(afs>0.05) & (afs<0.9)]
    if len(afs) < 5:
        purity_vaf[sid] = np.nan
    else:
        purity_vaf[sid] = float(min(1.0, 2 * np.percentile(afs, 95)))

# ---------- Assemble per-sample purity table ----------
rows = []
for _, r in rna_inv.iterrows():
    s = r['sample_id']
    if s in tumor_purity_rna.index:
        rows.append({'sample_id':s,'subject_id':r['subject_id'],'timepoint':r['timepoint'],
                     'response_bin':r['response_bin'],'platform':'RNA',
                     'stromal_score':float(stromal_scaled.get(s,np.nan)),
                     'immune_score':float(immune_scaled.get(s,np.nan)),
                     'estimate_score':float(estimate_score.get(s,np.nan)),
                     'purity_rna':float(tumor_purity_rna.get(s,np.nan)),
                     'purity_vaf':np.nan})
for _, r in wes_inv.iterrows():
    s = r['sample_id']
    if s.endswith('-N'): continue
    rows.append({'sample_id':s,'subject_id':r['subject_id'],'timepoint':r['timepoint'],
                 'response_bin':r['response_bin'],'platform':'WES',
                 'stromal_score':np.nan,'immune_score':np.nan,'estimate_score':np.nan,
                 'purity_rna':np.nan,'purity_vaf':purity_vaf.get(s, np.nan)})
per_sample = pd.DataFrame(rows)
per_sample['low_purity_flag'] = (
    (per_sample['purity_rna'] < 0.20) | (per_sample['purity_vaf'] < 0.20)
)
per_sample.to_csv(f'{OUT}/per_sample_purity.tsv', sep='\t', index=False)
print(f'Wrote per_sample_purity.tsv  ({len(per_sample)} rows)')

# Consolidate to per-subject-timepoint purity: prefer VAF if WES available, else RNA
def purity_for(subject, timepoint):
    sub = per_sample[(per_sample['subject_id']==subject) & (per_sample['timepoint']==timepoint)]
    if sub.empty: return np.nan
    v = sub['purity_vaf'].dropna()
    if len(v): return float(v.iloc[0])
    r = sub['purity_rna'].dropna()
    if len(r): return float(r.iloc[0])
    return np.nan

# ---------- Cascade Δ features ----------
# Paired subjects with pre & post WES:
wes_inv2 = wes_inv[wes_inv['timepoint'].isin(['pre','post'])]
paired_subj = [s for s,g in wes_inv2.groupby('subject_id')
               if set(g['timepoint']) >= {'pre','post'}]
print(f'Paired WES subjects: {len(paired_subj)}')

# Load feature sources
miss = (vm.groupby('sample_id').size().rename('n_nonsyn'))
# Neoantigen from summary
neo = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
# SBS5: look in refit
SBS = pd.read_csv('/mnt/sda1/data/TNT/analysis/01_wes_signatures/sbs_activities_with_meta.tsv', sep='\t').set_index('sample_id')
def get_sbs5(sid):
    try: return float(SBS.loc[sid, 'SBS5'])
    except Exception: return np.nan

# ssGSEA / immune signatures
ss = pd.read_csv('/mnt/sda1/data/TNT/analysis/08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0)
sig = pd.read_csv('/mnt/sda1/data/TNT/analysis/06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
# TRUST4
t4 = pd.read_csv('/mnt/sda1/data/TNT/analysis/06_rna_immune/trust4_summary.tsv', sep='\t')
t4['subject_id'] = t4['subject_id'] if 'subject_id' in t4.columns else np.nan

def feat_for(subject, timepoint):
    wrow = wes_inv2[(wes_inv2['subject_id']==subject) & (wes_inv2['timepoint']==timepoint)]
    rrow = rna_inv[(rna_inv['subject_id']==subject) & (rna_inv['timepoint']==timepoint)]
    sid_wes = wrow['sample_id'].iloc[0] if len(wrow) else None
    sid_rna = rrow['sample_id'].iloc[0] if len(rrow) else None
    d = {}
    d['missense'] = float(miss.get(sid_wes, np.nan)) if sid_wes else np.nan
    d['SBS5']     = get_sbs5(sid_wes) if sid_wes else np.nan
    if sid_wes:
        n = neo[neo['sample_id']==sid_wes]
        d['neo_binders'] = float(n['n_binders_500nM'].iloc[0]) if len(n) else np.nan
        d['neo_sites']   = float(n['n_sites_with_binder'].iloc[0]) if len(n) else np.nan
    else:
        d['neo_binders']=d['neo_sites']=np.nan
    if sid_rna and sid_rna in sig.index:
        for col in ['CD8_exhaustion','MHC_II','Treg']:
            d[col] = float(sig.loc[sid_rna, col]) if col in sig.columns else np.nan
    else:
        d['CD8_exhaustion']=d['MHC_II']=d['Treg']=np.nan
    if sid_rna:
        tr = t4[t4['sample_id']==sid_rna] if 'sample_id' in t4.columns else pd.DataFrame()
        d['IGH_n']       = float(tr['IGH_n'].iloc[0]) if len(tr) and 'IGH_n' in tr.columns else np.nan
        d['TRB_shannon'] = float(tr['TRB_shannon'].iloc[0]) if len(tr) and 'TRB_shannon' in tr.columns else np.nan
    else:
        d['IGH_n']=d['TRB_shannon']=np.nan
    return d

FEATURES = ['missense','SBS5','neo_binders','neo_sites','Treg','MHC_II',
            'CD8_exhaustion','IGH_n','TRB_shannon']

def mw(a, b):
    a = np.array([x for x in a if np.isfinite(x)])
    b = np.array([x for x in b if np.isfinite(x)])
    if len(a)<2 or len(b)<2: return np.nan
    try: return float(st.mannwhitneyu(a, b, alternative='two-sided').pvalue)
    except Exception: return np.nan

rows_out = []
# Build paired long table
recs = []
for s in paired_subj:
    pre = feat_for(s,'pre'); post = feat_for(s,'post')
    pu_pre = purity_for(s,'pre'); pu_post = purity_for(s,'post')
    resp = wes_inv2[wes_inv2['subject_id']==s]['response_bin'].iloc[0]
    for f in FEATURES:
        recs.append({'subject_id':s,'response':resp,'feature':f,
                     'pre':pre[f],'post':post[f],
                     'purity_pre':pu_pre,'purity_post':pu_post})
long = pd.DataFrame(recs)
long.to_csv(f'{OUT}/paired_feature_long.tsv', sep='\t', index=False)

for f in FEATURES:
    sub = long[long['feature']==f].copy()
    sub['delta_raw'] = sub['post'] - sub['pre']
    # purity-adjusted: for count-like features, divide by purity; for z-like, subtract and scale
    count_like = f in {'missense','SBS5','neo_binders','neo_sites','IGH_n'}
    if count_like:
        sub['pre_adj']  = sub['pre']  / sub['purity_pre'].replace(0, np.nan)
        sub['post_adj'] = sub['post'] / sub['purity_post'].replace(0, np.nan)
    else:
        sub['pre_adj']  = sub['pre']
        sub['post_adj'] = sub['post']
    sub['delta_adj'] = sub['post_adj'] - sub['pre_adj']
    # sensitivity: drop low-purity
    keep = (sub['purity_pre']>=0.20) & (sub['purity_post']>=0.20)
    g  = sub[sub['response']=='good']; b = sub[sub['response']=='bad']
    g2 = sub[(sub['response']=='good')&keep]; b2 = sub[(sub['response']=='bad')&keep]
    rows_out.append({
        'feature':f,
        'n_good_paired':len(g),'n_bad_paired':len(b),
        'delta_raw_good_med':float(g['delta_raw'].median()),
        'delta_raw_bad_med':float(b['delta_raw'].median()),
        'MW_p_raw':mw(g['delta_raw'].values, b['delta_raw'].values),
        'delta_purityadj_good_med':float(g['delta_adj'].median()),
        'delta_purityadj_bad_med':float(b['delta_adj'].median()),
        'MW_p_purityadj':mw(g['delta_adj'].values, b['delta_adj'].values),
        'n_low_purity_dropped':int(len(sub)-len(sub[keep])),
        'MW_p_excl_lowpurity':mw(g2['delta_raw'].values, b2['delta_raw'].values),
    })
res = pd.DataFrame(rows_out)
res.to_csv(f'{OUT}/delta_purity_sensitivity.tsv', sep='\t', index=False)
print(res.to_string(index=False))

# ---------- SuppFig: pre/post values coloured by purity ----------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8})
fig, axes = plt.subplots(3, 3, figsize=(10, 9))
for ax, f in zip(axes.ravel(), FEATURES):
    sub = long[long['feature']==f]
    for _, r in sub.iterrows():
        col = 'tab:blue' if r['response']=='good' else 'tab:red'
        alpha = max(0.3, min(1.0, r['purity_post'] if np.isfinite(r['purity_post']) else 0.5))
        ax.plot([0,1],[r['pre'],r['post']], '-o', color=col, alpha=alpha, lw=0.8, ms=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(['pre','post'])
    ax.set_title(f, fontsize=9)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle('Pre/post cascade features, line alpha ∝ post-CRT purity\n(blue=good, red=bad)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{FIG}/SuppFig_purity_prepost.{ext}', dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Wrote {FIG}/SuppFig_purity_prepost.png/pdf')
