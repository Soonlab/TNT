"""
SBS96 mutational signature refit via SigProfilerAssignment (COSMIC v3.3).
Input: /mnt/sda1/data/TNT/analysis/01_wes_signatures/vcf_input/*.vcf (49 Mutect2 PASS VCFs)
"""
import os, sys

VCF_DIR = '/mnt/sda1/data/TNT/analysis/01_wes_signatures/vcf_input'
OUT_DIR = '/mnt/sda1/data/TNT/analysis/01_wes_signatures/sbs_refit'
os.makedirs(OUT_DIR, exist_ok=True)

# Install GRCh38 genome for SigProfilerMatrixGenerator (first time only)
print('Installing GRCh38 for SigProfiler if needed...', flush=True)
try:
    from SigProfilerMatrixGenerator import install as genInstall
    genInstall.install('GRCh38', rsync=False, bash=True)
except Exception as e:
    print(f'genInstall skipped: {e}', flush=True)

print('Running SigProfilerAssignment refit...', flush=True)
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples=VCF_DIR,
                   output=OUT_DIR,
                   input_type='vcf',
                   context_type='96',
                   genome_build='GRCh38',
                   cosmic_version=3.3,
                   collapse_to_SBS96=True,
                   make_plots=True,
                   verbose=True)
print('DONE. Output:', OUT_DIR)
