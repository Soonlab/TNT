"""
External validation v3 — CD8 + RT focused re-analysis.

Fixes in this version (2026-04-15):
1. classify_response: 'Non-responder' bug fixed (non-* check before 'responder' substring).
2. Signature redesign:
   - CD8_cytotoxic  = pure CD8 T-cell effector markers (CD8A/B, GZMs, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11)
   - Tumor_cellcycle = proliferation markers (MKI67, TOP2A, MCMs, CCNBs, CDK1) — this is what the
     old 'CD8_proliferation' signature was really measuring (tumor-intrinsic, not lymphoid).
   - Keep DSB_HDR, E2F_MYC, EMT, and add B_cell_infiltration (per GSE150082 title emphasis).
3. Per-cohort manual TRG-scale override (Dworak vs Mandard vs CAP/Ryan).
4. Per-cohort regimen tag (nCRT-short / nCRT-long / TNT-matched) for stratified meta.
5. Add new cohorts: GSE87211, GSE133057, GSE46862, GSE56699, GSE233517 (if accessible).
"""
import GEOparse, pandas as pd, numpy as np, re, warnings
from pathlib import Path
from scipy import stats
warnings.filterwarnings('ignore')

OUT = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT/'geo_cache'; CACHE.mkdir(exist_ok=True)

# ---- Cohort metadata (manually curated regimen/scale) ----
# regimen: nCRT-long (long-course 45-50 Gy + concurrent chemo), nCRT-short (25 Gy/5fx),
#          TNT-like (CRT + induction or consolidation chemo), chemo-only (control)
COHORT_META = {
    'GSE150082': dict(regimen='nCRT-long', scale='explicit_good_poor',  tissue='rectum'),
    'GSE35452':  dict(regimen='nCRT-long', scale='responder_nonresp',   tissue='rectum'),
    'GSE119409': dict(regimen='nCRT-long', scale='sensitive_resistant', tissue='rectum'),
    'GSE45404':  dict(regimen='nCRT-long', scale='mandard',             tissue='rectum'),
    'GSE68204':  dict(regimen='nCRT-long', scale='mandard',             tissue='rectum'),
    'GSE94104':  dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE3493':   dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE56699':  dict(regimen='nCRT-long', scale='rcrg_3class',         tissue='rectum'),
    'GSE46862':  dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE133057': dict(regimen='nCRT-long', scale='survival',            tissue='rectum'),
    'GSE87211':  dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE119174': dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE69657':  dict(regimen='chemo-only',scale='responder_nonresp',   tissue='crc'),
    'GSE15781':  dict(regimen='nCRT-long', scale='paired_prepost',      tissue='rectum'),
    # New additions
    'GSE233517': dict(regimen='nCRT-long', scale='unknown',             tissue='rectum'),
    'GSE190826': dict(regimen='TNT-like',  scale='unknown',             tissue='rectum'),
}

# ---- Redesigned signatures ----
SIGS = {
  # Genuine CD8 effector / cytotoxic T cell markers (NOT mixed with cell cycle)
  'CD8_cytotoxic': ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1','IFNG','NKG7','GNLY',
                    'CXCL9','CXCL10','CXCL11','TBX21','EOMES','KLRK1','KLRD1'],
  # T cell infiltration (broader)
  'Tcell_infiltration': ['CD3D','CD3E','CD3G','CD2','CD4','CD8A','CD8B','LCK','ZAP70','ITK'],
  # B cell infiltration (GSE150082 emphasis)
  'Bcell_infiltration': ['CD19','CD20','MS4A1','CD79A','CD79B','CD22','TCL1A','FCRL5','BLK','FCER2'],
  # Tumor-intrinsic proliferation (what old CD8_proliferation was actually measuring)
  'Tumor_cellcycle': ['MKI67','TOP2A','STMN1','TYMS','UBE2C','BIRC5','CCNB1','CCNB2','CDK1',
                      'MCM2','MCM5','PCNA','CENPF','KIF20A','AURKA','AURKB','PLK1','BUB1'],
  # DNA repair / DSB-HDR
  'DSB_HDR_repair': ['BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR',
                     'CHEK1','CHEK2','MRE11','RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2',
                     'FANCI','FANCL','BLM','BRIP1','EXO1','DNA2','POLD1'],
  # E2F/MYC transcription programs
  'E2F_MYC_cellcycle': ['E2F1','E2F2','E2F3','MYC','MYCN','MCM3','MCM4','MCM6','MCM7',
                        'CCNE1','CCNE2','CDC20','CDC25A','CDC45','CDK2','CDK4','CDK6'],
  # EMT
  'EMT': ['VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','FN1','ITGB6','MMP2','MMP3','MMP9',
          'SOX10','ZEB1','ZEB2','TWIST2','TGFB1','TGFB2','COL1A1','COL1A2','COL3A1','FAP',
          'ACTA2','S100A4'],
}

# ---- Fixed response classifier (order: non- first, then positive) ----
def classify_response(val, scale=None):
    s = str(val).strip().lower()
    # scale-specific numeric mapping
    if scale == 'mandard':
        if re.match(r'^(trg)?\s*[12]$', s): return 'good'
        if re.match(r'^(trg)?\s*[345]$', s): return 'bad'
    if scale == 'dworak':
        # Dworak: 0=no regression=bad, 4=complete=good
        if re.match(r'^(trg)?\s*[34]$', s): return 'good'
        if re.match(r'^(trg)?\s*[012]$', s): return 'bad'
    if scale == 'cap_ryan':
        # CAP/Ryan: 0=CR=good, 3=no regression=bad
        if re.match(r'^(trg)?\s*[01]$', s): return 'good'
        if re.match(r'^(trg)?\s*[23]$', s): return 'bad'
    # Negative / resistant checked FIRST
    if any(k in s for k in ['non-responder','nonresponder','non responder','no response',
                            'non-response','no-response','poor','resistant','bad',
                            'trg2','trg3','trg4','trg5','trg g2','trg g3','trg g4','trg g5',
                            'relapse:yes','recurrence:yes','progressive']):
        return 'bad'
    if re.match(r'^pr\b', s) or re.match(r'^pd\b', s): return 'bad'
    # Positive / responder
    if any(k in s for k in ['responder','response:yes','complete response','pcr','near-cr',
                            'trg0','trg1','trg g0','trg g1','sensitive','good',
                            'no recurrence','no relapse','complete_clinical_response']):
        return 'good'
    if re.match(r'^cr\b', s): return 'good'
    return None

def score_sig(expr, genes):
    genes_hit = [g for g in genes if g in expr.index]
    if len(genes_hit) < 3: return None, 0
    sub = expr.loc[genes_hit]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), len(genes_hit)

def detect_response_col(pheno):
    prefs = ['response','responder','trg','tumor regression','pcr','sensitivity',
             'chemoresponse','pathologic','treatment response','ptrg','outcome',
             'class','regression','recurr','ajcc']
    cols = list(pheno.columns)
    low  = [str(c).lower() for c in cols]
    for p in prefs:
        for i, c in enumerate(low):
            if p in c: return cols[i]
    return None

# Per-cohort manual response mapping overrides (col_name, value->label dict)
MANUAL_RESP = {
    'GSE45404':  ('class', {'Responder':'good','Non Responder':'bad'}),
    'GSE46862':  ('chemoradiation therapy response',
                   # MO=Moderate, TO=Total (complete), MI=Minor, NT=No response
                   {'TO':'good','MO':'good','MI':'bad','NT':'bad'}),
    'GSE87211':  ('cancer recurrance after surgery',
                   # 0 = no recurrence = good, 1 = recurrence = bad (surrogate outcome)
                   {'0':'good','1':'bad'}),
    'GSE133057': ('ajcc score',
                   # AJCC TRG: 0=CR=good, 1=near-CR=good, 2=partial=bad, 3=no regression=bad
                   {'0':'good','1':'good','2':'bad','3':'bad'}),
    'GSE94104':  ('tumour regression grade',
                   # Rodel/Dworak in French LARC papers often: 1=poor,2=partial,3=good/complete — BUT
                   # the study (Rimini) uses Dworak where 3=very good+ complete response (good), 1=poor (bad)
                   # Rödel TRG: 0/1=bad, 2=intermediate, 3/4=good. With only 1/2/3 here, assume 3=good, 1=bad, 2=bad
                   {'1':'bad','2':'bad','3':'good'}),
    'GSE68204':  ('disease status', {}),  # no direct response field, skip
}

def load_expression(gse):
    for gpl_id, gpl in gse.gpls.items():
        tbl = gpl.table
        sym_col = None
        for c in ['Gene Symbol','gene_symbol','Symbol','GeneSymbol','Gene_Symbol',
                  'ILMN_Gene','SYMBOL','GENE_SYMBOL']:
            if c in tbl.columns: sym_col = c; break
        if sym_col is None and 'gene_assignment' in tbl.columns:
            # Affymetrix HuGene / HTA style
            probe_sym = dict(zip(tbl['ID'].astype(str),
                tbl['gene_assignment'].astype(str).str.split(' // ').str[1].str.strip()))
        elif sym_col is None:
            continue
        else:
            probe_sym = dict(zip(tbl['ID'].astype(str),
                tbl[sym_col].astype(str).str.split('///').str[0].str.strip()))
        mat_rows = []
        for gsm_id, gsm in gse.gsms.items():
            t = gsm.table
            if t is None or len(t)==0: continue
            id_col  = 'ID_REF' if 'ID_REF' in t.columns else t.columns[0]
            val_col = 'VALUE'  if 'VALUE'  in t.columns else t.columns[1]
            s = pd.Series(t[val_col].values, index=t[id_col].astype(str), name=gsm_id)
            mat_rows.append(s)
        if not mat_rows: continue
        mat = pd.concat(mat_rows, axis=1)
        mat.index = mat.index.map(lambda x: probe_sym.get(str(x), None))
        mat = mat[mat.index.notna() & (mat.index != '') & (mat.index != 'nan')]
        mat = mat.groupby(level=0).max()
        mat = mat.apply(pd.to_numeric, errors='coerce').dropna(how='all')
        if mat.shape[0] < 2000: continue
        if mat.max().max() > 50: mat = np.log2(mat + 1)
        return mat
    return None

def process(gid):
    meta = COHORT_META.get(gid, {})
    try:
        gse = GEOparse.get_GEO(geo=gid, destdir=str(CACHE), silent=True)
    except Exception as e:
        print(f'{gid}: FAIL download {e}'); return None
    pheno_rows = []
    for gsm_id, gsm in gse.gsms.items():
        row = {'sample_id': gsm_id}
        for c in gsm.metadata.get('characteristics_ch1', []):
            if ':' in c:
                k, v = c.split(':', 1)
                row[k.strip().lower()] = v.strip()
        pheno_rows.append(row)
    pheno = pd.DataFrame(pheno_rows).set_index('sample_id')
    expr = load_expression(gse)
    if expr is None:
        print(f'{gid}: no expression matrix')
        return None
    if gid in MANUAL_RESP and MANUAL_RESP[gid][1]:
        mcol, mmap = MANUAL_RESP[gid]
        if mcol in pheno.columns:
            classified = pheno[mcol].astype(str).map(lambda v: mmap.get(v.strip()))
            rcol = f'{mcol} [manual]'
            scale = 'manual'
        else:
            rcol = detect_response_col(pheno); scale = meta.get('scale')
            classified = pheno[rcol].apply(lambda v: classify_response(v, scale=scale)) if rcol else None
    else:
        rcol = detect_response_col(pheno)
        scale = meta.get('scale')
        classified = pheno[rcol].apply(lambda v: classify_response(v, scale=scale)) if rcol else None
    ng = int((classified=='good').sum()) if classified is not None else 0
    nb = int((classified=='bad').sum())  if classified is not None else 0
    print(f'{gid}: resp_col={rcol} scale={scale} good={ng} bad={nb}')
    scores = {}
    n_genes = {}
    for name, genes in SIGS.items():
        s, n = score_sig(expr, genes)
        if s is not None:
            scores[name] = s
            n_genes[name] = n
    sc = pd.DataFrame(scores).join(pheno)
    if classified is not None:
        sc['response_bin'] = classified
    sc['gse'] = gid
    sc['regimen'] = meta.get('regimen','?')
    sc.to_csv(OUT/f'{gid}_v3_scores.tsv', sep='\t')
    stats_rows = []
    if classified is not None and (ng>=3 and nb>=3):
        for sig in scores:
            g = sc[sc.response_bin=='good'][sig].dropna()
            b = sc[sc.response_bin=='bad'][sig].dropna()
            if len(g)>=3 and len(b)>=3:
                u = stats.mannwhitneyu(g,b)
                stats_rows.append({'gse':gid,'regimen':meta.get('regimen','?'),
                    'signature':sig,'n_good':len(g),'n_bad':len(b),
                    'mean_good':float(g.mean()),'mean_bad':float(b.mean()),
                    'delta':float(g.mean()-b.mean()),
                    'pvalue':float(u.pvalue),
                    'n_genes_in_sig':n_genes.get(sig,0)})
    return stats_rows, {'gse':gid,'n_samples':expr.shape[1],'n_probes':expr.shape[0],
                        'resp_col':rcol,'scale':scale,'regimen':meta.get('regimen','?'),
                        'n_good':ng,'n_bad':nb}

all_stats, summaries = [], []
for gid in COHORT_META:
    r = process(gid)
    if r is None:
        summaries.append({'gse':gid,'status':'FAILED'}); continue
    srows, summ = r
    summ['status']='OK'; summaries.append(summ); all_stats.extend(srows)

pd.DataFrame(summaries).to_csv(OUT/'v3_cohort_summary.tsv', sep='\t', index=False)
pd.DataFrame(all_stats).to_csv(OUT/'v3_signature_response_stats.tsv', sep='\t', index=False)

# ---- Meta-analysis (overall + stratified by regimen) ----
def stouffer(df, sig):
    sub = df[df.signature==sig].copy()
    if len(sub)==0: return None
    # one-sided Z in direction of discovery (good > bad)
    # Convert two-sided p to signed z using sign(delta)
    zs = []; ws = []
    for _,r in sub.iterrows():
        # two-sided p -> one-sided p matching sign of delta
        p_two = max(r.pvalue, 1e-300)
        z_mag = stats.norm.isf(p_two/2)  # magnitude
        z = z_mag * np.sign(r.delta)
        w = np.sqrt(r.n_good + r.n_bad)
        zs.append(z); ws.append(w)
    zs=np.array(zs); ws=np.array(ws)
    Z = float((ws*zs).sum() / np.sqrt((ws**2).sum()))
    p_meta = float(stats.norm.sf(abs(Z))*2)
    return {'signature':sig,'n_cohorts':len(sub),'Z':Z,'p_meta':p_meta,
            'deltas':','.join(f'{d:+.2f}' for d in sub.delta),
            'cohorts':','.join(sub.gse)}

sdf = pd.DataFrame(all_stats)
if len(sdf):
    meta_overall = [stouffer(sdf, s) for s in SIGS if stouffer(sdf, s)]
    pd.DataFrame(meta_overall).to_csv(OUT/'v3_meta_overall.tsv', sep='\t', index=False)
    strat_rows = []
    for reg, sub in sdf.groupby('regimen'):
        for s in SIGS:
            r = stouffer(sub, s)
            if r: r['regimen']=reg; strat_rows.append(r)
    pd.DataFrame(strat_rows).to_csv(OUT/'v3_meta_stratified.tsv', sep='\t', index=False)
    print('\n=== META OVERALL ===')
    print(pd.DataFrame(meta_overall).to_string(index=False))
    print('\n=== META STRATIFIED ===')
    print(pd.DataFrame(strat_rows).to_string(index=False))

print('\n=== Per-cohort summary ===')
print(pd.DataFrame(summaries).to_string(index=False))
