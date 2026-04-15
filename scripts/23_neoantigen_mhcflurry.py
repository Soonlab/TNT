"""
Direct neoantigen prediction via MHCflurry on VEP-annotated Mutect2 PASS VCFs.
For each missense variant, generate 8-11mer peptides spanning the altered residue,
predict MHC-I binding for subject-specific alleles, count high-binding neoantigens.
"""
import pandas as pd, numpy as np, subprocess as sp, os, re, gzip
from pathlib import Path
from scipy import stats

VEP_DIR = Path('/mnt/sda1/data/TNT/analysis/03_hla/neoantigen/vep')
OUT = Path('/mnt/sda1/data/TNT/analysis/03_hla/neoantigen/mhcflurry_direct'); OUT.mkdir(parents=True, exist_ok=True)
HLA = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv', sep='\t')
INV = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')

def parse_vep_vcf(path):
    """Extract protein changes from VEP-annotated VCF. Returns list of dicts with gene, protein context."""
    out=[]
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        ann_fields = None
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                # Extract field format
                m = re.search(r'Format: ([^"]+)', line)
                if m: ann_fields = m.group(1).split('|')
            if line.startswith('#'): continue
            parts = line.rstrip().split('\t')
            chrom,pos,_,ref,alt = parts[:5]
            info = parts[7]
            # Find CSQ
            csq = None
            for ent in info.split(';'):
                if ent.startswith('CSQ='):
                    csq = ent[4:]; break
            if csq is None or ann_fields is None: continue
            for trans_str in csq.split(','):
                v = trans_str.split('|')
                if len(v) != len(ann_fields): continue
                d = dict(zip(ann_fields, v))
                csq_terms = d.get('Consequence','')
                if 'missense' not in csq_terms: continue
                # Amino acids + position
                aa = d.get('Amino_acids','')  # e.g. R/Q
                pp = d.get('Protein_position','')
                if '/' not in aa: continue
                wt, mt = aa.split('/')
                if len(wt)!=1 or len(mt)!=1: continue
                sym = d.get('SYMBOL','')
                # Protein sequence needed for peptide — we'll only have local context via Amino_acids for now
                # Use ENSP and HGVSp if available (not full protein)
                out.append({'chrom':chrom,'pos':int(pos),'ref':ref,'alt':alt,
                            'gene':sym,'wt_aa':wt,'mt_aa':mt,'prot_pos':pp,
                            'ensp':d.get('ENSP',''),
                            'hgvsp':d.get('HGVSp',''),
                            'transcript':d.get('Feature','')})
    return out

def fetch_protein(ensp):
    """Fetch protein sequence from local VEP cache fasta or Ensembl REST."""
    import urllib.request, json
    if not ensp or not ensp.startswith('ENSP'): return None
    try:
        url = f'https://rest.ensembl.org/sequence/id/{ensp}?type=protein'
        req = urllib.request.Request(url, headers={'Content-Type':'application/json'})
        data = json.loads(urllib.request.urlopen(req, timeout=10).read())
        return data.get('seq','')
    except Exception as e:
        return None

PROT_CACHE = OUT/'protein_seq_cache.tsv'
protein_cache = {}
if PROT_CACHE.exists():
    c = pd.read_csv(PROT_CACHE, sep='\t')
    protein_cache = dict(zip(c.ensp, c.seq))

def get_peptides(mut):
    """Generate 8-11mer peptides covering the mutated residue."""
    if not mut['ensp'] or not mut['prot_pos']: return []
    seq = protein_cache.get(mut['ensp'])
    if seq is None:
        seq = fetch_protein(mut['ensp'])
        protein_cache[mut['ensp']] = seq or ''
    if not seq: return []
    try:
        pos = int(str(mut['prot_pos']).split('-')[0]) - 1
    except: return []
    if pos<0 or pos>=len(seq): return []
    if seq[pos] != mut['wt_aa']: return []  # mismatch
    mt_seq = seq[:pos] + mut['mt_aa'] + seq[pos+1:]
    peps=[]
    for L in [8,9,10,11]:
        for start in range(max(0,pos-L+1), min(len(mt_seq)-L+1, pos+1)):
            pep = mt_seq[start:start+L]
            wt_pep = seq[start:start+L]
            if 'X' in pep or '*' in pep: continue
            if pep == wt_pep: continue  # not neoantigen
            peps.append({'pep':pep,'wt_pep':wt_pep,'length':L,
                         'gene':mut['gene'],'chrom':mut['chrom'],'pos':mut['pos'],
                         'ensp':mut['ensp'],'prot_pos':mut['prot_pos']})
    return peps

# Step 1: enumerate peptides per sample
tumors = INV[(INV.timepoint!='normal') & (~INV.subject_id.isin([13,15,16,17,18,19,33]))]
print(f'{len(tumors)} matched tumors')

for _, t in tumors.iterrows():
    sid = t['sample_id']
    vcf = VEP_DIR/f'{sid}.vep.vcf.gz'
    if not vcf.exists(): continue
    pep_file = OUT/f'{sid}.peptides.tsv'
    if pep_file.exists(): continue
    muts = parse_vep_vcf(vcf)
    if not muts:
        pd.DataFrame(columns=['pep','wt_pep','length','gene','chrom','pos','ensp','prot_pos']).to_csv(pep_file, sep='\t', index=False)
        continue
    peps = []
    for m in muts:
        peps.extend(get_peptides(m))
    pd.DataFrame(peps).to_csv(pep_file, sep='\t', index=False)
    print(f'{sid}: {len(muts)} missense → {len(peps)} peptides')

# Save cache
pd.DataFrame({'ensp':list(protein_cache.keys()),'seq':list(protein_cache.values())}).to_csv(PROT_CACHE, sep='\t', index=False)

# Step 2: MHCflurry prediction per sample × subject alleles
import subprocess
MHCF = '/home/soon/miniconda3/envs/pvactools/bin/mhcflurry-predict'

results_rows = []
for _, t in tumors.iterrows():
    sid = t['sample_id']; subj = t['subject_id']
    pep_file = OUT/f'{sid}.peptides.tsv'
    if not pep_file.exists(): continue
    peps = pd.read_csv(pep_file, sep='\t')
    if len(peps)==0:
        results_rows.append({'sample_id':sid,'subject_id':subj,'n_peptides':0,'n_strong_binders':0,'n_weak_binders':0}); continue
    # Subject HLA alleles
    hla_sub = HLA[HLA.subject_id==subj]
    if len(hla_sub)==0: continue
    hs = hla_sub.iloc[0]
    alleles = []
    for c in ['A1','A2','B1','B2','C1','C2']:
        a = hs[c]
        if pd.isna(a): continue
        parts = str(a).split(':')
        if len(parts)>=2: alleles.append(':'.join(parts[:2]))
    alleles = list(dict.fromkeys(alleles))
    if not alleles: continue
    # Write peptide list
    plist = OUT/f'{sid}.pep_list.txt'
    peps['pep'].to_csv(plist, index=False, header=False)
    # Run mhcflurry-predict for each allele
    pred_out = OUT/f'{sid}.predictions.csv'
    if not pred_out.exists():
        cmd = [MHCF,'--alleles']+alleles+['--peptides']+peps['pep'].tolist()[:50000]+['--out',str(pred_out)]
        # too many args — use file input instead
        with open(plist,'w') as f:
            for p in peps['pep'].tolist(): f.write(p+'\n')
        cmd = [MHCF,'--alleles']+alleles+['--peptides-file',str(plist),'--out',str(pred_out)]
        r = sp.run(cmd, capture_output=True, text=True, timeout=900)
        if r.returncode != 0:
            print(f'{sid} mhcflurry FAIL: {r.stderr[-300:]}')
            continue
    pred = pd.read_csv(pred_out)
    # Binders
    strong = pred[pred['mhcflurry_affinity']<50]
    weak = pred[(pred['mhcflurry_affinity']>=50) & (pred['mhcflurry_affinity']<500)]
    results_rows.append({'sample_id':sid,'subject_id':subj,
                         'n_variants':peps['gene'].nunique(),
                         'n_peptides':len(peps),
                         'n_binders_500nM':(pred['mhcflurry_affinity']<500).sum(),
                         'n_strong_50nM':len(strong),
                         'n_weak_50_500nM':len(weak)})
    print(f'{sid}: {len(peps)} peps, {(pred["mhcflurry_affinity"]<500).sum()} binders <500nM')

res = pd.DataFrame(results_rows)
res = res.merge(INV[['sample_id','timepoint','response_bin']], on='sample_id')
res.to_csv(OUT/'neoantigen_summary.tsv', sep='\t', index=False)

print('\n=== Response association (pre, matched) ===')
pre = res[res.timepoint=='pre']
for col in ['n_variants','n_peptides','n_binders_500nM','n_strong_50nM']:
    g = pre[pre.response_bin=='good'][col]
    b = pre[pre.response_bin=='bad'][col]
    if len(g)>=3 and len(b)>=3:
        u = stats.mannwhitneyu(g, b)
        print(f'  {col}: good med={g.median():.0f} vs bad med={b.median():.0f}  p={u.pvalue:.3f}')
