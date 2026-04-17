"""
Download GSE109057 (Akiyoshi T et al, Br J Surg 2019;106(10):1381-92, PMID 31268561)
and score 4 Thread 1 + 3 Thread 2 signatures.

GSE109057 — n=90 LARC pre-CRT biopsies, Affymetrix PrimeView (GPL15207). Original paper
reports MCP-counter cytotoxic lymphocyte score TRG1 vs TRG3/4 P=0.01 (eFig 5 of Akiyoshi
2023 JAMA Netw Open). 81/90 patients overlap with GSE216616 but the platform is different
(microarray here vs RNA-seq there) and we already use GSE216616 paper-level only.

Output:
  gse109057_pheno.tsv
  gse109057_thread12_scores.tsv
  gse109057_thread12_stats.tsv
  Fig_GSE109057_thread12.{pdf,png}  (4 Thread1 + 3 Thread2 panels)
"""
import os, re, warnings
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
GID = 'GSE109057'

# Thread 1 + Thread 2 signatures
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
EXPECTED = {'DSB_HDR_repair':+1,'E2F_MYC_cellcycle':+1,'Tumor_cellcycle':+1,'EMT':-1,
            'CD8_cytotoxic':+1,'Tcell_infiltration':+1,'Bcell_infiltration':+1}
THREAD = {'DSB_HDR_repair':1,'E2F_MYC_cellcycle':1,'Tumor_cellcycle':1,'EMT':1,
          'CD8_cytotoxic':2,'Tcell_infiltration':2,'Bcell_infiltration':2}

# ---- 1. Download ----
print(f'Fetching {GID} ...')
gse = GEOparse.get_GEO(geo=GID, destdir=CACHE, silent=True)
print(f'  platforms: {list(gse.gpls.keys())}, n samples: {len(gse.gsms)}')

# ---- 2. Phenotype ----
pheno_rows = []
for gsm_id, gsm in gse.gsms.items():
    row = {'sample_id': gsm_id}
    for c in gsm.metadata.get('characteristics_ch1', []):
        if ':' in c:
            k, v = c.split(':', 1)
            row[k.strip().lower()] = v.strip()
    row['title'] = gsm.metadata.get('title', [''])[0]
    pheno_rows.append(row)
pheno = pd.DataFrame(pheno_rows).set_index('sample_id')
print(f'\npheno columns: {list(pheno.columns)}')
print(pheno.head().to_string())

# Find TRG / response column
trg_col = None
for c in pheno.columns:
    cl = c.lower()
    if 'trg' in cl or 'regression' in cl or 'response' in cl or 'group' in cl:
        trg_col = c; break
print(f'\ntrg/response col candidate: "{trg_col}"')
if trg_col:
    print(f'value counts:\n{pheno[trg_col].value_counts(dropna=False)}')

# Akiyoshi 2023 eFig 5 says GSE109057 has TRG1 (n=54) vs TRG3/4 (n=36) split. Map to good/bad
def classify(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return None
    s = str(v).strip().lower()
    # explicit TRG numeric
    if re.match(r'^(trg)?\s*1$', s) or 'trg1' in s and 'trg3' not in s and 'trg4' not in s:
        return 'good'
    if re.match(r'^(trg)?\s*[34]$', s) or 'trg3' in s or 'trg4' in s or 'trg 3' in s or 'trg 4' in s:
        return 'bad'
    if re.match(r'^(trg)?\s*2$', s) or 'trg2' in s:
        return 'good'  # Akiyoshi groups TRG1/2 as good
    if any(k in s for k in ['responder','complete','sensitive','good']) and \
       not any(k in s for k in ['non','poor','resistant','bad']):
        return 'good'
    if any(k in s for k in ['non-responder','nonresponder','non response','poor','resistant','bad']):
        return 'bad'
    return None

pheno['response_bin'] = pheno[trg_col].apply(classify) if trg_col else None
print(f'\nclassified:\n{pheno.response_bin.value_counts(dropna=False)}')
pheno.to_csv(f'{OUT}/gse109057_pheno.tsv', sep='\t')

# ---- 3. Expression ----
gpl = list(gse.gpls.values())[0]
print(f'\nplatform = {gpl.name}, annotation shape = {gpl.table.shape}')
print(f'annotation cols: {list(gpl.table.columns)[:15]}')

# PrimeView: find symbol col
sym_col = None
for c in gpl.table.columns:
    if c.lower() in ('gene symbol','gene_symbol','symbol','genesymbol'):
        sym_col = c; break
if sym_col is None:
    for c in gpl.table.columns:
        if 'symbol' in c.lower(): sym_col = c; break
print(f'  symbol column: "{sym_col}"')
probe2sym = gpl.table[['ID', sym_col]].dropna().rename(columns={sym_col:'symbol'})
probe2sym['symbol'] = probe2sym['symbol'].astype(str).str.split(' /// ').str[0].str.strip()
probe2sym = probe2sym[probe2sym.symbol != ''].set_index('ID')['symbol']
print(f'  probe -> symbol: {len(probe2sym)}')

# Build expression
exprs = []
for gsm_id, gsm in gse.gsms.items():
    s = gsm.table.set_index('ID_REF')['VALUE']; s.name = gsm_id
    exprs.append(s)
expr = pd.concat(exprs, axis=1).dropna(how='any')
print(f'  expression matrix: {expr.shape}')
median_max = expr.max().median()
if median_max > 50:
    print(f'  median max = {median_max:.0f} -> log2'); expr = np.log2(expr + 1)
else:
    print(f'  median max = {median_max:.2f} -> already log')

# Probe -> gene aggregation
expr2 = expr.copy(); expr2['symbol'] = expr2.index.map(probe2sym)
expr2 = expr2.dropna(subset=['symbol'])
gene_expr = expr2.groupby('symbol').mean()
print(f'  gene-level matrix: {gene_expr.shape}')

# z-score
zexpr = gene_expr.sub(gene_expr.mean(axis=1), axis=0).div(gene_expr.std(axis=1), axis=0).dropna(how='any')

# Score
scores = {}; coverage = []
for sig, genes in SIGS.items():
    found = [g for g in genes if g in zexpr.index]
    coverage.append({'signature': sig, 'n_total': len(genes), 'n_found': len(found),
                     'pct': round(100*len(found)/len(genes),1)})
    scores[sig] = zexpr.loc[found].mean(axis=0) if found else np.nan
score_df = pd.DataFrame(scores); score_df.index.name='sample_id'
score_df = score_df.merge(pheno[['response_bin']], left_index=True, right_index=True)
score_df.to_csv(f'{OUT}/gse109057_thread12_scores.tsv', sep='\t')
print(f'\n=== gene coverage ===\n{pd.DataFrame(coverage).to_string(index=False)}')

# ---- 5. MW test ----
rows = []
for sig in SIGS:
    g = score_df[score_df.response_bin == 'good'][sig].dropna().values
    b = score_df[score_df.response_bin == 'bad'][sig].dropna().values
    if len(g) < 3 or len(b) < 3:
        rows.append({'thread': THREAD[sig], 'signature': sig, 'n_good': len(g), 'n_bad': len(b),
                     'mean_good': np.nan, 'mean_bad': np.nan, 'delta': np.nan,
                     'mw_p': np.nan, 'expected_dir': EXPECTED[sig], 'concordant': np.nan})
        continue
    u = stats.mannwhitneyu(g, b, alternative='two-sided')
    delta = g.mean() - b.mean()
    obs_dir = int(np.sign(delta))
    rows.append({'thread': THREAD[sig], 'signature': sig, 'n_good': len(g), 'n_bad': len(b),
                 'mean_good': round(g.mean(),3), 'mean_bad': round(b.mean(),3),
                 'delta': round(delta,3), 'mw_p': round(u.pvalue, 4),
                 'expected_dir': EXPECTED[sig],
                 'concordant': int(obs_dir == EXPECTED[sig])})
stats_df = pd.DataFrame(rows)
stats_df.to_csv(f'{OUT}/gse109057_thread12_stats.tsv', sep='\t', index=False)
print(f'\n=== Thread 1+2 vs response in GSE109057 ===')
print(stats_df.to_string(index=False))
n_t1_conc = stats_df[stats_df.thread==1]['concordant'].sum()
n_t2_conc = stats_df[stats_df.thread==2]['concordant'].sum()
print(f'\nThread 1 concordance: {int(n_t1_conc)}/4')
print(f'Thread 2 concordance: {int(n_t2_conc)}/3')

# ---- 6. Figure ----
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})
fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.5))
all_sigs = list(SIGS.keys())
for ax, sig in zip(axes.flatten(), all_sigs + [None]):
    if sig is None:
        ax.set_visible(False); continue
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
        ax.scatter(rng.normal(0, 0.06, len(g)), g, color='#2E86AB', s=14, edgecolor='black', lw=0.3)
        ax.scatter(rng.normal(1, 0.06, len(b)), b, color='#E63946', s=14, edgecolor='black', lw=0.3)
        u = stats.mannwhitneyu(g, b)
        delta = g.mean() - b.mean()
        marker = '✓' if (np.sign(delta)==EXPECTED[sig]) else '✗'
        thread_tag = 'T1' if THREAD[sig]==1 else 'T2'
    ax.set_xticks([0,1]); ax.set_xticklabels([f'good\nn={len(g)}', f'bad\nn={len(b)}'])
    ax.set_title(f'[{thread_tag}] {sig}\nΔ={delta:+.3f} {marker}  P={u.pvalue:.3f}', fontsize=8.5)
    ax.set_ylabel('z-score', fontsize=8)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle(f'GSE109057 (Akiyoshi 2019, Br J Surg) — Thread 1 + Thread 2 vs response (n={len(score_df)})',
             fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/Fig_GSE109057_thread12.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'\nWrote pheno + scores + stats + Fig_GSE109057_thread12')
