"""
External validation v2: better response classification + more cohorts
"""
import GEOparse, pandas as pd, numpy as np, os, re
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

OUT = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT/'geo_cache'; CACHE.mkdir(exist_ok=True)

# Expanded list (confirmed accessible + with response annotation)
COHORTS = [
  'GSE150082','GSE94104','GSE87211','GSE15781','GSE35452','GSE45404',
  'GSE3493','GSE68204','GSE119409','GSE69657','GSE104645','GSE143985',
  'GSE41258','GSE46862','GSE56699','GSE133057','GSE190826','GSE119174',
  'GSE93375','GSE5851','GSE15960','GSE22058','GSE37892','GSE59857',
  'GSE30363','GSE45582','GSE73255','GSE76473','GSE39582','GSE39396',
]

SIGS = {
  'DSB_HDR_repair': ['BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR','CHEK1','CHEK2','MRE11','RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2','FANCI','FANCL','BLM','BRIP1','EXO1','DNA2','POLD1'],
  'E2F_MYC_cellcycle': ['E2F1','E2F2','E2F3','MYC','MYCN','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','CCNB1','CCNB2','CDK1','CDK2','CDK4','CDK6','CCNE1','CCNE2','CDC20','CDC25A','CDC45','PCNA','TOP2A','MKI67','PLK1','AURKA','AURKB','BUB1','BUB1B','TYMS'],
  'CD8_proliferation': ['MKI67','TOP2A','STMN1','TYMS','TUBB','UBE2C','BIRC5','CCNB1','CCNB2','CDK1','MCM2','MCM5','PCNA','CENPF','KIF20A'],
  'EMT': ['VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','FN1','ITGB6','MMP2','MMP3','MMP9','SOX10','ZEB1','ZEB2','TWIST2','TGFB1','TGFB2','COL1A1','COL1A2','COL3A1','FAP','ACTA2','S100A4'],
}

def score_sig(expr, genes):
  genes = [g for g in genes if g in expr.index]
  if len(genes)<3: return None
  sub = expr.loc[genes]
  z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0,np.nan), axis=0)
  return z.mean(axis=0)

def detect_response_col(meta_df):
  # priority: explicit response labels
  prefs = ['trg','tumor regression','response','responder','pcr','sensitivity',
           'chemoresponse','responsiveness','pathologic','ypt','recurrence',
           'relapse','survival','outcome','treatment response']
  cols = [str(c).lower() for c in meta_df.columns]
  for p in prefs:
    for i, c in enumerate(cols):
      if p in c: return meta_df.columns[i]
  return None

def classify_response(val):
  s = str(val).strip().lower()
  # Explicit non-responder check FIRST (before 'responder' partial match)
  if any(k in s for k in ['non-responder','nonresponder','non responder','no response',
                          'non-response','no-response','poor','resistant','trg2','trg3','trg4',
                          'trg g2','trg g3','relapse:yes','recurrence:yes','bad']):
    return 'bad'
  if any(k in s for k in ['pr ','partial response','pd ','progressive']): return 'bad'
  # Then responder / complete
  if any(k in s for k in ['responder','response:yes','complete response','cr ','pcr','near-cr',
                          'trg0','trg1','trg g0','trg g1','sensitive','good','no recurrence','no relapse']):
    return 'good'
  # Numeric TRG
  if re.match(r'^\s*(trg)?\s*0\s*$', s): return 'good'
  if re.match(r'^\s*(trg)?\s*1\s*$', s): return 'good'
  if re.match(r'^\s*(trg)?\s*[234]\s*$', s): return 'bad'
  return None

def process_gse(gse_id):
  try:
    print(f'\n===== {gse_id} =====')
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(CACHE), silent=True)
  except Exception as e:
    print(f' FAIL get_GEO: {e}')
    return None
  samples = gse.gsms
  pheno_rows = []
  for gsm_id, gsm in samples.items():
    meta = gsm.metadata
    chars = meta.get('characteristics_ch1', [])
    row = {'sample_id': gsm_id}
    for c in chars:
      if ':' in c:
        k, v = c.split(':', 1)
        row[k.strip().lower()] = v.strip()
    pheno_rows.append(row)
  pheno = pd.DataFrame(pheno_rows).set_index('sample_id')
  # Expression
  plats = gse.gpls
  expr = None
  for gpl_id, gpl in plats.items():
    tbl = gpl.table
    sym_col = None
    for c in ['Gene Symbol','gene_symbol','Symbol','GeneSymbol','Gene_Symbol','ILMN_Gene','SYMBOL','GENE_SYMBOL','gene_assignment']:
      if c in tbl.columns: sym_col=c; break
    if sym_col is None: continue
    probe_sym = dict(zip(tbl['ID'].astype(str), tbl[sym_col].astype(str).str.split('///').str[0].str.split(' // ').str[1 if 'gene_assignment' in str(sym_col) else 0].str.strip() if sym_col=='gene_assignment' else tbl[sym_col].astype(str).str.split('///').str[0].str.strip()))
    mat_rows=[]
    for gsm_id, gsm in samples.items():
      t = gsm.table
      if t is None or len(t)==0: continue
      id_col = 'ID_REF' if 'ID_REF' in t.columns else t.columns[0]
      val_col = 'VALUE' if 'VALUE' in t.columns else t.columns[1]
      s = pd.Series(t[val_col].values, index=t[id_col].astype(str), name=gsm_id)
      mat_rows.append(s)
    if not mat_rows: continue
    mat = pd.concat(mat_rows, axis=1)
    mat.index = mat.index.map(lambda x: probe_sym.get(str(x), None))
    mat = mat[mat.index.notna() & (mat.index != '') & (mat.index != 'nan')]
    mat = mat.groupby(level=0).max()
    mat = mat.apply(pd.to_numeric, errors='coerce').dropna(how='all')
    if mat.shape[0] < 2000: continue
    expr = mat; break
  if expr is None:
    print(f' no usable expression')
    return None
  if expr.max().max() > 50:
    expr = np.log2(expr + 1)
  resp_col = detect_response_col(pheno)
  classified = None
  if resp_col is not None:
    classified = pheno[resp_col].apply(classify_response)
  # Compute sigs
  scores = {name: score_sig(expr, genes) for name, genes in SIGS.items()}
  scores = {k:v for k,v in scores.items() if v is not None}
  sc_df = pd.DataFrame(scores)
  sc_df = sc_df.join(pheno)
  sc_df['response_bin'] = classified
  sc_df.to_csv(OUT/f'{gse_id}_signature_scores.tsv', sep='\t')
  # Response stats
  stats_rows = []
  if classified is not None and (classified=='good').sum()>=3 and (classified=='bad').sum()>=3:
    for sig in scores:
      g = sc_df[sc_df.response_bin=='good'][sig].dropna()
      b = sc_df[sc_df.response_bin=='bad'][sig].dropna()
      if len(g)>=3 and len(b)>=3:
        u = stats.mannwhitneyu(g,b)
        stats_rows.append({'gse':gse_id,'signature':sig,'n_good':len(g),'n_bad':len(b),
                           'mean_good':float(g.mean()),'mean_bad':float(b.mean()),
                           'delta':float(g.mean()-b.mean()),'pvalue':float(u.pvalue)})
  n_good = (classified=='good').sum() if classified is not None else 0
  n_bad = (classified=='bad').sum() if classified is not None else 0
  print(f' {gse_id}: probes={expr.shape[0]} samples={expr.shape[1]} resp_col={resp_col} good={n_good} bad={n_bad}')
  return stats_rows, expr.shape[1], resp_col, n_good, n_bad

all_stats = []
summary = []
for gid in COHORTS:
  r = process_gse(gid)
  if r is None:
    summary.append({'gse':gid,'status':'FAILED'}); continue
  st, n, rc, ng, nb = r
  all_stats.extend(st)
  summary.append({'gse':gid,'status':'OK','n_samples':n,'resp_col':rc,'n_good':ng,'n_bad':nb,'n_stats':len(st)})

pd.DataFrame(summary).to_csv(OUT/'external_cohort_summary.tsv', sep='\t', index=False)
if all_stats:
  S = pd.DataFrame(all_stats)
  # Meta combine p-values per signature (Stouffer/Fisher)
  S_wide = S.pivot_table(index='gse', columns='signature', values='delta')
  S_p = S.pivot_table(index='gse', columns='signature', values='pvalue')
  print('\n=== All stats ===')
  print(S.to_string(index=False))
  # Meta by signature
  meta_rows=[]
  for sig in S['signature'].unique():
    sub = S[S['signature']==sig]
    # Stouffer with sample size weighting, one-sided alternative = expect delta>0 for DSB/E2F/CD8, <0 for EMT
    if sig == 'EMT':
      z = [stats.norm.isf(min(max(p/2,1e-10),1-1e-10)) * (-1 if d>0 else 1) for d,p in zip(sub['delta'],sub['pvalue'])]
    else:
      z = [stats.norm.isf(min(max(p/2,1e-10),1-1e-10)) * (1 if d>0 else -1) for d,p in zip(sub['delta'],sub['pvalue'])]
    z = np.array(z, dtype=float)
    w = np.sqrt(sub['n_good']+sub['n_bad'])
    Z = np.sum(z*w)/np.sqrt(np.sum(w**2))
    p_meta = 2*(1-stats.norm.cdf(abs(Z)))
    meta_rows.append({'signature':sig,'n_cohorts':len(sub),'Z':Z,'p_meta':p_meta,
                     'deltas':','.join(f'{x:.2f}' for x in sub['delta'])})
  pd.DataFrame(meta_rows).to_csv(OUT/'external_meta_analysis.tsv', sep='\t', index=False)
  print('\n=== Meta-analysis ===')
  print(pd.DataFrame(meta_rows).to_string(index=False))
  S.to_csv(OUT/'external_signature_response_stats.tsv', sep='\t', index=False)
print('\nDone. Cohort summary:')
print(pd.DataFrame(summary).to_string(index=False))
