"""
SC-RT external validation on GSE254249 (Gao, Ling et al. Cancer Cell 2025,
PMID 41202810) — the only public SC-RT-TNT bulk RNA-seq dataset identified by
the 2026-04-20 search (see SCRT_external_search.md).

Cohort (bulk RNA-seq slice, 11 samples / 8 patients, all in nRCT arm = TNT =
short-course RT 5x5 Gy + 6-8 cycles FOLFOXIRI):
  - 3 PRE-treatment tumor samples  (CRC15, CRC24, CRC25; SampleTimePoint="BL")
  - 8 POST-treatment tumor samples (CRC15-T..CRC26-T; SampleTimePoint="nRCT")
  - Paired pre+post: 3 subjects (CRC15 CR, CRC24 non-CR, CRC25 non-CR)
  - Efficacy label: CR / non-CR, propagated from -T to paired pre sample

Three analyses (mirroring our discovery + paired framework):
  A. Pre-treatment baseline × response     (n=3, directional only)
  B. Post-treatment × response             (n=8, 5 CR vs 3 non-CR)
  C. Paired Δ(post-pre) target engagement  (n=3 paired, discovery direction
     = DSB down / cellcycle down / EMT up after TNT)

Signatures re-used verbatim from scripts/32, 260418_add/13_gse109057_score.py
— 4 Thread 1 (tumor-intrinsic) + 3 Thread 2 (immune).

Output:
  gse254249_bulk_pheno.tsv
  gse254249_scores.tsv                (per-sample z-score per signature)
  gse254249_pre_response_stats.tsv    (A)
  gse254249_post_response_stats.tsv   (B)
  gse254249_paired_delta_stats.tsv    (C)
  Fig_GSE254249_pre_boxplot.{pdf,png}
  Fig_GSE254249_post_boxplot.{pdf,png}
  Fig_GSE254249_paired_delta.{pdf,png}
"""
import os, gzip, warnings
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT    = '/data/data/TNT/analysis/260418_add'
INDIR  = f'{OUT}/gse254249'
EXPR_F = f'{INDIR}/GSE254249_bulkRNA_logTPM.tsv.gz'
META_F = f'{INDIR}/GSE254249_bulkRNA_metadata.tsv.gz'

GOOD = '#0a7d6e'
BAD  = '#c53e1f'

# ---- Signature panels (match scripts/32 + 13_gse109057_score.py) ----
SIGS = {
    'DSB_HDR_repair':    ['BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR',
                          'CHEK1','CHEK2','MRE11','RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2',
                          'FANCI','FANCL','BLM','BRIP1','EXO1','DNA2','POLD1'],
    'E2F_MYC_cellcycle': ['E2F1','E2F2','E2F3','MYC','MYCN','MCM3','MCM4','MCM6','MCM7',
                          'CCNE1','CCNE2','CDC20','CDC25A','CDC45','CDK2','CDK4','CDK6'],
    'Tumor_cellcycle':   ['MKI67','TOP2A','STMN1','TYMS','UBE2C','BIRC5','CCNB1','CCNB2','CDK1',
                          'MCM2','MCM5','PCNA','CENPF','KIF20A','AURKA','AURKB','PLK1','BUB1'],
    'EMT':               ['VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','FN1','ITGB6','MMP2','MMP3','MMP9',
                          'SOX10','ZEB1','ZEB2','TWIST2','TGFB1','TGFB2','COL1A1','COL1A2','COL3A1','FAP',
                          'ACTA2','S100A4'],
    'CD8_cytotoxic':     ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1','IFNG','NKG7','GNLY',
                          'CXCL9','CXCL10','CXCL11','TBX21','EOMES','KLRK1','KLRD1'],
    'Tcell_infiltration':['CD3D','CD3E','CD3G','CD2','CD4','CD8A','CD8B','LCK','ZAP70','ITK'],
    'Bcell_infiltration':['CD19','CD20','MS4A1','CD79A','CD79B','CD22','TCL1A','FCRL5','BLK','FCER2'],
}
# Discovery/external direction convention:
#   baseline (pre-treatment) good responders: Thread 1 UP, EMT DOWN, Thread 2 UP
#   target engagement Δ(post - pre) in good responders: cellcycle/DSB DOWN, EMT UP
EXPECTED_PRE   = {'DSB_HDR_repair':+1,'E2F_MYC_cellcycle':+1,'Tumor_cellcycle':+1,'EMT':-1,
                  'CD8_cytotoxic':+1,'Tcell_infiltration':+1,'Bcell_infiltration':+1}
EXPECTED_DELTA = {'DSB_HDR_repair':-1,'E2F_MYC_cellcycle':-1,'Tumor_cellcycle':-1,'EMT':+1,
                  'CD8_cytotoxic':+1,'Tcell_infiltration':+1,'Bcell_infiltration':+1}
THREAD = {'DSB_HDR_repair':1,'E2F_MYC_cellcycle':1,'Tumor_cellcycle':1,'EMT':1,
          'CD8_cytotoxic':2,'Tcell_infiltration':2,'Bcell_infiltration':2}

# ---- 1. Load metadata + expression ----
pheno = pd.read_csv(META_F, sep='\t', index_col=0)
pheno.columns = [c.strip() for c in pheno.columns]
pheno.index.name = 'sample_id'
print('\n=== pheno ===')
print(pheno.to_string())

# map SampleTimePoint -> timepoint
tp_map = {'BL':'pre', 'nRCT':'post'}
pheno['timepoint'] = pheno['SampleTimePoint'].map(tp_map)
pheno['subject']   = pheno.index.to_series().str.replace(r'-T$', '', regex=True)

# propagate -T efficacy back to paired pre sample
resp_by_subj = (pheno[pheno.timepoint=='post']
                  .reset_index().set_index('subject')['Efficacy']
                  .to_dict())
pheno['response'] = pheno.apply(
    lambda r: r['Efficacy'] if r['timepoint']=='post' else resp_by_subj.get(r['subject']),
    axis=1)
# bin: CR -> good, non-CR -> bad (align with discovery convention)
pheno['response_bin'] = pheno['response'].map({'CR':'good','non-CR':'bad'})
print('\n=== pheno enriched ===')
print(pheno.to_string())
pheno.to_csv(f'{OUT}/gse254249_bulk_pheno.tsv', sep='\t')

expr = pd.read_csv(EXPR_F, sep='\t', index_col=0)
expr.index.name = 'gene'
print(f'\nexpr: {expr.shape}  (log1p(TPM), per paper Methods)')

# ---- 2. Per-sample z-score of mapped genes (ssGSEA-like mean z) ----
zexpr = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1), axis=0).dropna(how='any')
print(f'zexpr: {zexpr.shape}')

scores, cov = {}, []
for sig, genes in SIGS.items():
    found = [g for g in genes if g in zexpr.index]
    cov.append({'signature':sig,'n_total':len(genes),'n_found':len(found),
                'pct':round(100*len(found)/len(genes),1)})
    if found:
        scores[sig] = zexpr.loc[found].mean(axis=0)
score_df = pd.DataFrame(scores)
print('\n=== gene coverage ===')
print(pd.DataFrame(cov).to_string(index=False))

score_df = score_df.merge(pheno[['subject','timepoint','response_bin','Efficacy']],
                           left_index=True, right_index=True)
score_df.to_csv(f'{OUT}/gse254249_scores.tsv', sep='\t')
print('\n=== scores (first rows) ===')
print(score_df.to_string())

def fmt(x, n=3):
    return 'NA' if pd.isna(x) else round(float(x), n)

def mw_test(g, b, expect):
    if len(g)<2 or len(b)<2:
        return {'n_good':len(g),'n_bad':len(b),'mean_good':fmt(np.mean(g) if len(g) else np.nan),
                'mean_bad':fmt(np.mean(b) if len(b) else np.nan),
                'delta':fmt((np.mean(g)-np.mean(b)) if (len(g) and len(b)) else np.nan),
                'mw_p':'NA','expected_dir':expect,'concordant':'NA'}
    u = stats.mannwhitneyu(g, b, alternative='two-sided')
    delta = np.mean(g) - np.mean(b)
    conc = int(np.sign(delta) == expect) if delta != 0 else 0
    return {'n_good':len(g),'n_bad':len(b),'mean_good':fmt(np.mean(g)),
            'mean_bad':fmt(np.mean(b)),'delta':fmt(delta),
            'mw_p':fmt(u.pvalue,4),'expected_dir':expect,'concordant':conc}

# ---- 3A. Pre-treatment × response ----
pre = score_df[score_df.timepoint=='pre']
print(f'\n=== PRE samples (n={len(pre)}): {pre.index.tolist()} ===')
print(f'  response split: {pre.response_bin.value_counts(dropna=False).to_dict()}')
rows = []
for sig in SIGS:
    g = pre[pre.response_bin=='good'][sig].dropna().values
    b = pre[pre.response_bin=='bad'][sig].dropna().values
    r = mw_test(g, b, EXPECTED_PRE[sig])
    rows.append({'signature':sig,'thread':THREAD[sig], **r})
pre_stats = pd.DataFrame(rows)
pre_stats.to_csv(f'{OUT}/gse254249_pre_response_stats.tsv', sep='\t', index=False)
print('\n=== A. PRE × response stats ===')
print(pre_stats.to_string(index=False))

# ---- 3B. Post-treatment × response ----
post = score_df[score_df.timepoint=='post']
print(f'\n=== POST samples (n={len(post)}) ===')
print(f'  response split: {post.response_bin.value_counts().to_dict()}')
rows = []
for sig in SIGS:
    g = post[post.response_bin=='good'][sig].dropna().values
    b = post[post.response_bin=='bad'][sig].dropna().values
    r = mw_test(g, b, EXPECTED_PRE[sig])  # same "good-up" direction as pre
    rows.append({'signature':sig,'thread':THREAD[sig], **r})
post_stats = pd.DataFrame(rows)
post_stats.to_csv(f'{OUT}/gse254249_post_response_stats.tsv', sep='\t', index=False)
print('\n=== B. POST × response stats ===')
print(post_stats.to_string(index=False))

# ---- 3C. Paired Δ(post - pre) target engagement ----
paired_subjects = sorted(set(pre.subject) & set(post.subject))
print(f'\n=== PAIRED subjects (n={len(paired_subjects)}): {paired_subjects} ===')
delta_rows = []
for sig in SIGS:
    deltas = []
    for s in paired_subjects:
        pre_v  = pre[pre.subject==s][sig].values
        post_v = post[post.subject==s][sig].values
        if len(pre_v)==1 and len(post_v)==1:
            deltas.append({'subject':s, 'delta':post_v[0]-pre_v[0],
                           'response':post[post.subject==s].response_bin.iloc[0]})
    if not deltas:
        continue
    dd = pd.DataFrame(deltas)
    # discovery predicted direction in target engagement
    exp = EXPECTED_DELTA[sig]
    n_conc = int((np.sign(dd['delta']) == exp).sum())
    # sign test (binomial)
    if len(dd) >= 2:
        binom_p = stats.binomtest(n_conc, len(dd), 0.5, alternative='two-sided').pvalue
    else:
        binom_p = np.nan
    delta_rows.append({'signature':sig,'thread':THREAD[sig],'n_paired':len(dd),
                       'mean_delta':fmt(dd['delta'].mean()),
                       'expected_dir':exp,
                       'n_concordant':n_conc,'binom_p':fmt(binom_p,4),
                       'deltas':';'.join([f'{r.subject}:{r.delta:+.2f}' for r in dd.itertuples()])})
paired_stats = pd.DataFrame(delta_rows)
paired_stats.to_csv(f'{OUT}/gse254249_paired_delta_stats.tsv', sep='\t', index=False)
print('\n=== C. PAIRED Δ(post-pre) target engagement ===')
print(paired_stats.to_string(index=False))

# ---- 4. Figures ----
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})

def panel_box(ax, g, b, sig, expect, title_prefix):
    if len(g) and len(b):
        parts = ax.boxplot([g, b], positions=[0,1], widths=0.55, patch_artist=True,
                           medianprops=dict(color='black', lw=1.2),
                           boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7),
                           capprops=dict(lw=0.7), flierprops=dict(marker='o', ms=3))
        for patch, c in zip(parts['boxes'], [GOOD, BAD]):
            patch.set_facecolor(c); patch.set_alpha(0.55)
        rng = np.random.RandomState(0)
        ax.scatter(rng.normal(0, 0.06, len(g)), g, color=GOOD, s=22, edgecolor='black', lw=0.3, zorder=3)
        ax.scatter(rng.normal(1, 0.06, len(b)), b, color=BAD, s=22, edgecolor='black', lw=0.3, zorder=3)
    try:
        u = stats.mannwhitneyu(g, b)
        p = u.pvalue
    except Exception:
        p = np.nan
    delta = (np.mean(g)-np.mean(b)) if (len(g) and len(b)) else np.nan
    marker = '✓' if (not np.isnan(delta)) and (np.sign(delta)==expect) else '✗'
    expsym = '+' if expect>0 else '-'
    ax.set_xticks([0,1]); ax.set_xticklabels([f'good\n(n={len(g)})',f'bad\n(n={len(b)})'])
    ax.set_title(f'{sig}\nΔ={delta:+.2f} (exp {expsym}) {marker}\n'
                 f'MW P={p:.3f}' if not np.isnan(p) else f'{sig}\nn insufficient',
                 fontsize=8)
    ax.set_ylabel('z-score')
    for s in ['top','right']: ax.spines[s].set_visible(False)

# --- Fig: pre × response ---
fig, axes = plt.subplots(1, 7, figsize=(16, 3.3), sharey=False)
for ax, sig in zip(axes, SIGS):
    g = pre[pre.response_bin=='good'][sig].dropna().values
    b = pre[pre.response_bin=='bad'][sig].dropna().values
    panel_box(ax, g, b, sig, EXPECTED_PRE[sig], 'pre')
fig.suptitle(f'A. Pre-treatment baseline × response — GSE254249 nRCT arm (SC-RT+FOLFOXIRI), n={len(pre)}',
             fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_GSE254249_pre_boxplot.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)

# --- Fig: post × response ---
fig, axes = plt.subplots(1, 7, figsize=(16, 3.3), sharey=False)
for ax, sig in zip(axes, SIGS):
    g = post[post.response_bin=='good'][sig].dropna().values
    b = post[post.response_bin=='bad'][sig].dropna().values
    panel_box(ax, g, b, sig, EXPECTED_PRE[sig], 'post')
fig.suptitle(f'B. Post-TNT × response — GSE254249 nRCT arm, n={len(post)} (5 CR vs 3 non-CR)',
             fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_GSE254249_post_boxplot.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)

# --- Fig: paired Δ slopegraph ---
fig, axes = plt.subplots(1, 7, figsize=(16, 3.5), sharey=False)
for ax, sig in zip(axes, SIGS):
    for s in paired_subjects:
        pre_v  = pre[pre.subject==s][sig].values
        post_v = post[post.subject==s][sig].values
        if len(pre_v)==1 and len(post_v)==1:
            resp = post[post.subject==s].response_bin.iloc[0]
            c = GOOD if resp=='good' else BAD
            ax.plot([0,1], [pre_v[0], post_v[0]], '-o', color=c, lw=1.3,
                    markersize=6, markeredgecolor='black', markeredgewidth=0.4, alpha=0.85)
            ax.text(1.05, post_v[0], s, fontsize=7, va='center')
    exp = EXPECTED_DELTA[sig]
    expsym = '+' if exp>0 else '-'
    stat_row = paired_stats[paired_stats.signature==sig].iloc[0]
    ax.set_xticks([0,1]); ax.set_xticklabels(['pre','post'])
    ax.set_title(f'{sig}\nΔ mean={stat_row.mean_delta} (exp {expsym})\n'
                 f'{stat_row.n_concordant}/{stat_row.n_paired} concord; sign P={stat_row.binom_p}',
                 fontsize=8)
    ax.set_ylabel('z-score')
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle(f'C. Paired Δ(post-pre) target engagement — GSE254249 nRCT (SC-RT+FOLFOXIRI), n={len(paired_subjects)}',
             fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_GSE254249_paired_delta.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)

print('\n✓ Outputs:')
for f in ['gse254249_bulk_pheno.tsv','gse254249_scores.tsv',
          'gse254249_pre_response_stats.tsv','gse254249_post_response_stats.tsv',
          'gse254249_paired_delta_stats.tsv',
          'Fig_GSE254249_pre_boxplot.{png,pdf}',
          'Fig_GSE254249_post_boxplot.{png,pdf}',
          'Fig_GSE254249_paired_delta.{png,pdf}']:
    print(f'  {OUT}/{f}')
