#!/bin/bash
# pVACseq rerun (v2) — IDENTICAL config to v1 (2026-04-14) for reproducibility check
# v1 output: 03_wes_hla_neoantigen/pvacseq/
# v2 output: 03_wes_hla_neoantigen/pvacseq_v2/
set -u

source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate pvactools

BASE=/mnt/sda1/data/TNT/analysis
VEP_DIR=$BASE/03_wes_hla_neoantigen/vep
OUT=$BASE/03_wes_hla_neoantigen/pvacseq_v2
LOG_DIR=$BASE/03_wes_hla_neoantigen/logs_v2
mkdir -p $OUT $LOG_DIR

HLA_TSV=$BASE/03_hla/hla_class_I_typing.tsv
WES_INV=$BASE/00_cohort/wes_inventory.tsv
PIPE_LOG=$LOG_DIR/pipeline_v2.log
PVACSEQ=/home/soon/miniconda3/envs/pvactools/bin/pvacseq
BCFTOOLS=/home/soon/miniconda3/envs/wes_somatic/bin/bcftools
IEDB_DIR=/home/soon/miniconda3/envs/pvactools

echo "===== $(date '+%F %T') v2 pVACseq rerun start =====" > $PIPE_LOG
echo "v1 config replicated: MHCflurry, -e1 8,9,10,11, --pass-only, HLA 2-field" >> $PIPE_LOG

/home/soon/miniconda3/envs/pvactools/bin/python3 << PYEOF 2>&1 | tee -a $PIPE_LOG
import pandas as pd, subprocess as sp
from pathlib import Path
import time

BASE = Path('$BASE')
VEP_DIR = Path('$VEP_DIR')
OUT = Path('$OUT')
LOG_DIR = Path('$LOG_DIR')

hla = pd.read_csv('$HLA_TSV', sep='\t')
inv = pd.read_csv('$WES_INV', sep='\t')

# v1 rule: all tumors (including unmatched/hypermutator), no exclusion
tumors = inv[inv.timepoint != 'normal'].sort_values('sample_id')
print(f'[{time.strftime("%T")}] tumors to process: {len(tumors)}', flush=True)

for _, r in tumors.iterrows():
    sid = r['sample_id']; subj = int(r['subject_id'])
    vep_vcf = VEP_DIR/f'{sid}.vep.vcf.gz'
    if not vep_vcf.exists():
        print(f'[skip] {sid}: VEP VCF missing', flush=True); continue

    outd = OUT/sid
    done_marker = outd/'MHC_Class_I'/f'{sid}_DNA.MHC_I.all_epitopes.aggregated.tsv'
    if done_marker.exists():
        print(f'[skip] {sid}: already complete', flush=True); continue
    outd.mkdir(parents=True, exist_ok=True)

    hla_sub = hla[hla.subject_id == subj]
    if len(hla_sub) == 0:
        print(f'[skip] {sid}: no HLA typing', flush=True); continue
    hs = hla_sub.iloc[0]
    alleles = []
    for c in ['A1','A2','B1','B2','C1','C2']:
        a = hs[c]
        if pd.isna(a): continue
        parts = str(a).split(':')
        if len(parts) >= 2:
            alleles.append(':'.join(parts[:2]))
    alleles = list(dict.fromkeys(alleles))
    if not alleles:
        print(f'[skip] {sid}: no usable HLA alleles', flush=True); continue
    allele_str = ','.join(alleles)

    proc = sp.run(['$BCFTOOLS','query','-l',str(vep_vcf)], capture_output=True, text=True)
    samples = [s for s in proc.stdout.strip().split('\n') if s]
    tumor_sm = next((s for s in samples
                     if not s.endswith('-N') and not s.endswith('-N_DNA')), samples[0])

    cmd = ['$PVACSEQ','run', str(vep_vcf), tumor_sm, allele_str,
           'MHCflurry', str(outd),
           '-e1','8,9,10,11',
           '-t','2',
           '--pass-only']

    t0 = time.time()
    print(f'[{time.strftime("%T")}] [run] {sid} subj={subj} tumor={tumor_sm} alleles={allele_str}', flush=True)
    with open(LOG_DIR/f'pvacseq_{sid}.log','w') as lf:
        r2 = sp.run(cmd, stdout=lf, stderr=sp.STDOUT, timeout=3600)
    dt = time.time() - t0
    print(f'[{time.strftime("%T")}] [done] {sid} exit={r2.returncode} ({dt:.0f}s)', flush=True)

print(f'[{time.strftime("%T")}] ===== v2 pipeline finished =====', flush=True)
PYEOF

echo "===== $(date '+%F %T') v2 pVACseq rerun end =====" >> $PIPE_LOG
