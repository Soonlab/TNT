"""
External validation of TNT response signature in public GEO LARC cohorts.
Tests the key TNT-response signature (DSB repair + HDR + E2F/Myc + CD8 proliferation
+ EMT-down) in multiple independent datasets.
"""
import GEOparse, pandas as pd, numpy as np, os, re
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

OUT = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT/'geo_cache'; CACHE.mkdir(exist_ok=True)

# Candidate LARC neoadjuvant CRT/TNT response cohorts
# key: GSE_id, response column hint, tumor type, context
COHORTS = [
  ('GSE150082', 'LARC, neoadjuvant chemoradiotherapy, Response'),
  ('GSE94104',  'LARC nCRT response'),
  ('GSE87211',  'LARC pre/post nCRT'),
  ('GSE15781',  'LARC pre/post nCRT'),
  ('GSE35452',  'LARC nCRT response'),
  ('GSE45404',  'LARC response'),
  ('GSE3493',   'LARC CRT early'),
  ('GSE68204',  'LARC nCRT'),
  ('GSE46862',  'LARC/CRC CRT'),
  ('GSE56699',  'LARC nCRT'),
  ('GSE133057', 'LARC nCRT'),
  ('GSE190826', 'LARC TNT'),
  ('GSE119409', 'LARC TNT'),
  ('GSE119174', 'rectal nCRT'),
  ('GSE93375',  'LARC nCRT'),
  ('GSE69657',  'LARC nCRT'),
  ('GSE5851',   'colorectal cetuximab — negative control'),
]

# Our key signature (from session findings)
# Use TPM space gene symbols
SIGS = {
  'DSB_HDR_repair': ['BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR','CHEK1','CHEK2','MRE11','RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2','FANCI','FANCL','BLM','BRIP1','EXO1','DNA2','POLD1'],
  'E2F_MYC_cellcycle': ['E2F1','E2F2','E2F3','MYC','MYCN','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','CCNB1','CCNB2','CDK1','CDK2','CDK4','CDK6','CCNE1','CCNE2','CDC20','CDC25A','CDC45','PCNA','TOP2A','MKI67','PLK1','AURKA','AURKB','BUB1','BUB1B','TYMS'],
  'CD8_proliferation': ['MKI67','TOP2A','STMN1','TYMS','TUBB','UBE2C','BIRC5','CCNB1','CCNB2','CDK1','MCM2','MCM5','PCNA','CENPF','KIF20A'],
  'EMT': ['VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','FN1','ITGB6','MMP2','MMP3','MMP9','SOX10','ZEB1','ZEB2','TWIST2','TGFB1','TGFB2','COL1A1','COL1A2','COL3A1','FAP','ACTA2','S100A4'],
}

def score_sig_sample(expr, genes):
  genes = [g for g in genes if g in expr.index]
  if len(genes)<3: return None
  sub = expr.loc[genes]
  z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0,np.nan), axis=0)
  return z.mean(axis=0)

def detect_response_col(meta_df):
  for c in meta_df.columns:
    lc = str(c).lower()
    if any(k in lc for k in ['response','responder','trg','tumor regression','tumor_regression','pcr','pathologic','ypt','ypn','responsiveness','resistan']):
      return c
  return None

def classify_response(val):
  s = str(val).lower()
  if any(k in s for k in ['cr','complete','trg0','trg1','responder','good','sensitive']): return 'good'
  if any(k in s for k in ['pr','partial','trg2','trg3','trg4','non-responder','nonresponder','poor','bad','resist']): return 'bad'
  if re.match(r'^[01]$', s.strip()): return 'good' if s.strip()=='0' else 'bad'
  return None

def process_gse(gse_id, note):
  try:
    print(f'\n===== {gse_id} {note} =====')
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(CACHE), silent=True)
    # Assemble expression + phenotype
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
    try:
      expr = None
      # Try fetching series_matrix style
      import re as _re
      # platform-based
      plats = gse.gpls
      for gpl_id, gpl in plats.items():
        tbl = gpl.table
        # Map probe → gene symbol
        sym_col = None
        for c in ['Gene Symbol','gene_symbol','Symbol','GeneSymbol','Gene_Symbol','ILMN_Gene','SYMBOL','GENE_SYMBOL']:
          if c in tbl.columns: sym_col=c; break
        if sym_col is None: continue
        probe_sym = dict(zip(tbl['ID'].astype(str), tbl[sym_col].astype(str).str.split('///').str[0].str.strip()))
        # Build expression matrix
        mat_rows = []
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
        mat = mat[mat.index.notna() & (mat.index != '')]
        mat = mat.groupby(level=0).max()
        # Convert to numeric
        mat = mat.apply(pd.to_numeric, errors='coerce')
        mat = mat.dropna(how='all')
        expr = mat
        break
      if expr is None or expr.shape[0]<1000:
        print(f' {gse_id}: no usable expression ({None if expr is None else expr.shape})')
        return None
      # Log2 if not already
      if (expr.max().max() > 50):
        expr = np.log2(expr + 1)
      # Detect response column
      resp_col = detect_response_col(pheno)
      if resp_col is None:
        print(f' {gse_id}: no response column; chars cols = {list(pheno.columns)[:10]}')
      classified = None
      if resp_col is not None:
        classified = pheno[resp_col].apply(classify_response)
        ok = classified.notna().sum()
        print(f' {gse_id}: resp_col={resp_col}, classified n={ok} (good={sum(classified=="good")}, bad={sum(classified=="bad")})')
      # Compute sig scores
      scores = {}
      for name, genes in SIGS.items():
        s = score_sig_sample(expr, genes)
        if s is not None: scores[name] = s
      sc_df = pd.DataFrame(scores)
      sc_df = sc_df.join(pheno)
      sc_df['response_bin'] = classified
      sc_df.to_csv(OUT/f'{gse_id}_signature_scores.tsv', sep='\t')
      # Response association
      stats_rows = []
      if classified is not None and classified.notna().sum() >= 10:
        for sig in scores:
          g = sc_df[sc_df.response_bin=='good'][sig].dropna()
          b = sc_df[sc_df.response_bin=='bad'][sig].dropna()
          if len(g)>=3 and len(b)>=3:
            u = stats.mannwhitneyu(g,b)
            stats_rows.append({'gse':gse_id,'signature':sig,'n_good':len(g),'n_bad':len(b),
                               'mean_good':g.mean(),'mean_bad':b.mean(),'delta':g.mean()-b.mean(),
                               'pvalue':float(u.pvalue)})
      return stats_rows, sc_df, expr.shape, resp_col
    except Exception as e:
      print(f' {gse_id}: expr error {e}')
      return None
  except Exception as e:
    print(f' {gse_id}: FAILED - {e}')
    return None

all_stats = []
summary_rows = []
for gid, note in COHORTS:
  res = process_gse(gid, note)
  if res is None:
    summary_rows.append({'gse':gid,'status':'FAILED','n_sig':None,'n_resp':None})
    continue
  stats_rows, sc_df, shape, resp_col = res
  all_stats.extend(stats_rows)
  summary_rows.append({'gse':gid,'status':'OK','n_probes':shape[0],'n_samples':shape[1],'resp_col':resp_col,
                       'n_good':(sc_df.get('response_bin')=='good').sum(),
                       'n_bad':(sc_df.get('response_bin')=='bad').sum()})

pd.DataFrame(summary_rows).to_csv(OUT/'external_cohort_summary.tsv', sep='\t', index=False)
pd.DataFrame(all_stats).to_csv(OUT/'external_signature_response_stats.tsv', sep='\t', index=False)

print('\n=== External validation summary ===')
print(pd.DataFrame(summary_rows).to_string(index=False))
print('\n=== Signature response stats across cohorts ===')
if all_stats:
  df = pd.DataFrame(all_stats).sort_values(['signature','pvalue'])
  print(df.to_string(index=False))
