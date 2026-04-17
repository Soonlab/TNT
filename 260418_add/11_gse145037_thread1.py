"""
Download GSE145037 (Zhang et al, Front Cell Dev Biol 2021; PMID 34485268) and test whether
the discovery's Thread 1 four-feature panel (DSB_HDR_repair, E2F_MYC_cellcycle,
Tumor_cellcycle, EMT) separates the original paper's binary response label.

GSE145037 — n=31 LARC pre-CRT biopsies, standard long-course nCRT, AJCC TRG -> binary
"response to the crt: response / non-response" stored in GEO metadata.
Platform: Affymetrix PrimeView GPL15207 / GPL28140.

Output:
  gse145037_pheno.tsv          — per-sample phenotype + response classification
  gse145037_thread1_scores.tsv — per-sample 4-signature ssGSEA-like z-score
  gse145037_thread1_stats.tsv  — per-signature good vs bad MW + delta + concordance
  Fig_GSE145037_thread1.{pdf,png}  — 4-panel boxplot good vs bad
"""
import os, sys, re, warnings
import numpy as np, pandas as pd
import GEOparse
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = '/data/data/TNT/analysis/260418_add'
CACHE = '/data/data/TNT/analysis/11_external_validation/geo_cache'
os.makedirs(CACHE, exist_ok=True)

GID = 'GSE145037'

# Thread 1 signatures (same as scripts/32)
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
}
EXPECTED_DIR = {'DSB_HDR_repair': +1, 'E2F_MYC_cellcycle': +1,
                'Tumor_cellcycle': +1, 'EMT': -1}

# ---- 1. Download ----
print(f'Fetching {GID} ...')
gse = GEOparse.get_GEO(geo=GID, destdir=CACHE, silent=False)
print(f'  platforms: {list(gse.gpls.keys())}')
print(f'  n samples: {len(gse.gsms)}')

# ---- 2. Phenotype ----
pheno_rows = []
for gsm_id, gsm in gse.gsms.items():
    row = {'sample_id': gsm_id}
    for c in gsm.metadata.get('characteristics_ch1', []):
        if ':' in c:
            k, v = c.split(':', 1)
            row[k.strip().lower()] = v.strip()
    pheno_rows.append(row)
pheno = pd.DataFrame(pheno_rows).set_index('sample_id')
print(f'\npheno columns: {list(pheno.columns)}')
print(pheno.head().to_string())

# Find response column
resp_col_candidates = [c for c in pheno.columns if 'response' in c.lower() or 'crt' in c.lower()]
print(f'\nresponse column candidates: {resp_col_candidates}')

# Use the first match
resp_col = resp_col_candidates[0] if resp_col_candidates else None
if resp_col is None:
    raise SystemExit('no response col found')
print(f'\nUsing response col: "{resp_col}"')
print('value counts:'); print(pheno[resp_col].value_counts())

def classify(v):
    s = str(v).strip().lower()
    if any(k in s for k in ['non-response','nonresponse','non response','no response',
                            'non-responder','poor','resistant','bad']):
        return 'bad'
    if any(k in s for k in ['response','responder','complete','sensitive','good']):
        return 'good'
    return None
pheno['response_bin'] = pheno[resp_col].apply(classify)
print(f'\nclassified:\n{pheno.response_bin.value_counts(dropna=False)}')
pheno.to_csv(f'{OUT}/gse145037_pheno.tsv', sep='\t')

# ---- 3. Expression ----
# GPL28140 ID column already contains gene symbols (this dataset was deposited as gene-level)
gpl = list(gse.gpls.values())[0]
print(f'\nplatform = {gpl.name}; annotation shape = {gpl.table.shape}')
print(f'platform cols: {list(gpl.table.columns)} (ID already = gene symbol)')

exprs = []
for gsm_id, gsm in gse.gsms.items():
    s = gsm.table.set_index('ID_REF')['VALUE']
    s.name = gsm_id
    exprs.append(s)
expr = pd.concat(exprs, axis=1)  # genes x samples (ID == symbol)
print(f'  expression matrix: {expr.shape}')
expr = expr.dropna(how='any')
print(f'  after dropna: {expr.shape}')

median_max = expr.max().median()
if median_max > 50:
    print(f'  median sample max = {median_max:.0f} -> log2 transform')
    expr = np.log2(expr + 1)
else:
    print(f'  median sample max = {median_max:.2f} -> already log-scale')

# ID == symbol directly; collapse duplicates if any
gene_expr = expr.groupby(expr.index).mean()
print(f'  gene-level matrix: {gene_expr.shape}')

# ---- 4. Score signatures (z-score of mean across mapped genes per sample) ----
# Match scripts/32 approach: per-sample = mean(z-scores of mapped genes)
zexpr = gene_expr.sub(gene_expr.mean(axis=1), axis=0).div(gene_expr.std(axis=1), axis=0)
zexpr = zexpr.dropna(how='any')

scores = {}
coverage = []
for sig, genes in SIGS.items():
    found = [g for g in genes if g in zexpr.index]
    coverage.append({'signature': sig, 'n_total': len(genes), 'n_found': len(found),
                     'pct': round(100*len(found)/len(genes),1)})
    if found:
        scores[sig] = zexpr.loc[found].mean(axis=0)
score_df = pd.DataFrame(scores)
print(f'\n=== gene coverage ===\n{pd.DataFrame(coverage).to_string(index=False)}')

score_df.index.name = 'sample_id'
score_df = score_df.merge(pheno[['response_bin']], left_index=True, right_index=True)
score_df.to_csv(f'{OUT}/gse145037_thread1_scores.tsv', sep='\t')

# ---- 5. Per-signature MW test ----
rows = []
for sig in SIGS:
    g = score_df[score_df.response_bin == 'good'][sig].dropna().values
    b = score_df[score_df.response_bin == 'bad'][sig].dropna().values
    if len(g) < 3 or len(b) < 3:
        rows.append({'signature': sig, 'n_good': len(g), 'n_bad': len(b),
                     'mean_good': np.nan, 'mean_bad': np.nan, 'delta': np.nan,
                     'mw_p': np.nan, 'expected_dir': EXPECTED_DIR[sig], 'concordant': np.nan})
        continue
    u = stats.mannwhitneyu(g, b, alternative='two-sided')
    delta = g.mean() - b.mean()
    obs_dir = int(np.sign(delta))
    rows.append({'signature': sig, 'n_good': len(g), 'n_bad': len(b),
                 'mean_good': round(g.mean(),3), 'mean_bad': round(b.mean(),3),
                 'delta': round(delta,3), 'mw_p': round(u.pvalue, 4),
                 'expected_dir': EXPECTED_DIR[sig],
                 'concordant': int(obs_dir == EXPECTED_DIR[sig])})
stats_df = pd.DataFrame(rows)
stats_df.to_csv(f'{OUT}/gse145037_thread1_stats.tsv', sep='\t', index=False)
print(f'\n=== Thread 1 vs response in GSE145037 ===')
print(stats_df.to_string(index=False))
n_conc = stats_df['concordant'].sum()
print(f'\nConcordance with discovery direction: {int(n_conc)}/{len(SIGS)}')
n_sig = ((stats_df['mw_p'] < 0.05) & (stats_df['concordant']==1)).sum()
n_sig_disc = ((stats_df['mw_p'] < 0.05) & (stats_df['concordant']==0)).sum()
print(f'  significant concordant:    {int(n_sig)}')
print(f'  significant DISCORDANT:    {int(n_sig_disc)}')

# ---- 6. Figure: 4-panel boxplot ----
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})
fig, axes = plt.subplots(1, 4, figsize=(11, 3.4))
for ax, sig in zip(axes, SIGS):
    g = score_df[score_df.response_bin=='good'][sig].dropna().values
    b = score_df[score_df.response_bin=='bad'][sig].dropna().values
    if len(g) >= 3 and len(b) >= 3:
        parts = ax.boxplot([g, b], positions=[0, 1], widths=0.55, patch_artist=True,
                            medianprops=dict(color='black', lw=1.2),
                            boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7),
                            capprops=dict(lw=0.7), flierprops=dict(marker='o', ms=3))
        for patch, c in zip(parts['boxes'], ['#2E86AB','#E63946']):
            patch.set_facecolor(c); patch.set_alpha(0.55)
        rng = np.random.RandomState(0)
        ax.scatter(rng.normal(0, 0.06, len(g)), g, color='#2E86AB', s=18, edgecolor='black', lw=0.3)
        ax.scatter(rng.normal(1, 0.06, len(b)), b, color='#E63946', s=18, edgecolor='black', lw=0.3)
        u = stats.mannwhitneyu(g, b)
        delta = g.mean() - b.mean()
        exp = '+' if EXPECTED_DIR[sig]==+1 else '-'
        obs = '+' if delta>0 else '-'
        marker = '✓' if (np.sign(delta)==EXPECTED_DIR[sig]) else '✗'
    ax.set_xticks([0,1]); ax.set_xticklabels([f'good\n(n={len(g)})', f'bad\n(n={len(b)})'])
    ax.set_title(f'{sig}\nΔ={delta:+.3f} (exp {exp}) {marker}\nMW P={u.pvalue:.3f}', fontsize=8.5)
    ax.set_ylabel('z-score')
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle(f'Thread 1 (tumor-intrinsic) signatures vs response — GSE145037 (Zhang 2021, n={len(score_df)})',
             fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_GSE145037_thread1.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'\nWrote Fig_GSE145037_thread1, gse145037_thread1_stats.tsv, gse145037_thread1_scores.tsv')
