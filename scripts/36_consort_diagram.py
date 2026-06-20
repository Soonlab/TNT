"""CONSORT-style sample-flow diagram reconciling all analysis n's.
Addresses reviewer Major #6 (Methods internal inconsistency).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

FIG = Path('/mnt/sda1/data/TNT/analysis/figures/supp'); FIG.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(11, 13))
ax.set_xlim(0, 12); ax.set_ylim(0, 18)
ax.axis('off')

def box(x, y, w, h, text, color='#e8f0fe', edge='#1f77b4', lw=1.2, fs=9):
    p = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle='round,pad=0.05',
                       fc=color, ec=edge, lw=lw)
    ax.add_patch(p)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, wrap=True)

def arrow(x1, y1, x2, y2, label='', style='-|>', lw=1.2, color='k'):
    a = FancyArrowPatch((x1,y1),(x2,y2), arrowstyle=style, mutation_scale=12, lw=lw, color=color)
    ax.add_patch(a)
    if label:
        ax.text((x1+x2)/2+0.15, (y1+y2)/2, label, fontsize=8, color='gray', style='italic')

# ---- Title ----
ax.text(6, 17.3, 'Sample flow — TNT discovery cohort (N=35 patients)',
        ha='center', fontsize=13, fontweight='bold')
ax.text(6, 16.8, 'Reconciles per-analysis n across WES, RNA-seq, paired, and integrated tables',
        ha='center', fontsize=9, color='gray', style='italic')

# ---- Enrollment ----
box(6, 15.8, 5.5, 0.7,
    '35 MSS LARC patients enrolled\nAge 34–78 (median 60); M 25 / F 10; cT2–T4',
    color='#fff3cd', edge='#b58a00', fs=9)

# ---- Response label ----
box(6, 14.7, 5.5, 0.7,
    'Final TNT response (Dworak TRG, post-consolidation) — binary good vs bad\n'
    'good (TRG 0–1) n=18 | bad (TRG 2–3) n=17',
    color='#fff3cd', edge='#b58a00', fs=9)
arrow(6, 15.4, 6, 15.05)

# ---- Paired-set structure ----
box(2.8, 13.4, 4.8, 0.8,
    'prepost_set = Y (paired)\nn=14 patients\n'
    '3 samples each (normal/pre-CRT/post-CRT) = 42',
    color='#d4edda', edge='#2e7d32', fs=8.5)
box(9.2, 13.4, 4.8, 0.8,
    'prepost_set = N (single timepoint)\nn=21 patients\n'
    '30 with normal+pre + 6 pre-only = 36',
    color='#f8d7da', edge='#b71c1c', fs=8.5)
arrow(6, 14.35, 3.3, 13.8); arrow(6, 14.35, 8.7, 13.8)

# ---- WES panel ----
box(3, 11.8, 4.4, 0.7,
    'WES sequenced samples: 77\n(14×3 paired + 36 non-paired; subj 13-N missing)',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(2.8, 13.0, 3.0, 12.15)
arrow(9.2, 13.0, 3.2, 12.15, color='gray')

box(3, 10.6, 4.4, 0.7,
    'Mutect2 somatic calling\n49 tumors PASS (41 T-N matched + 8 tumor-only/PoN)',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(3, 11.45, 3, 10.95)

box(3, 9.4, 4.4, 0.7,
    'MSI (msisensor-pro): 41 matched — all MSS\n'
    'SBS refit: 49 samples\nCNVkit CNV/CIN: 49 samples',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(3, 10.25, 3, 9.75)

box(3, 8.1, 4.4, 0.8,
    'Driver oncoprint, TMB: 49 tumors\n'
    '(pre-CRT TMB analysis uses 35 pre-CRT tumors)\n'
    'HLA class I OptiType: 35 subjects (from normals)',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(3, 9.05, 3, 8.5)

box(3, 6.7, 4.4, 0.8,
    'Paired pre→post (WES-derived):\n'
    '14 subjects — SBS5 Δ, HLA-LOH (28 tumors)\n'
    'PyClone-VI converged on 12/14 subjects',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(3, 7.7, 3, 7.1)

box(3, 5.3, 4.4, 0.8,
    'pVACseq MHCflurry neoantigens: 49 tumors\n'
    'Paired Δ binders — 11 paired subjects\n'
    '(subj 3, 11 excluded: post-CRT VCF incomplete)',
    color='#e3f2fd', edge='#1976d2', fs=8.5)
arrow(3, 6.3, 3, 5.7)

# ---- RNA-seq panel ----
box(9, 11.8, 4.4, 0.7,
    'RNA-seq sequenced samples: 56\n(27 paired subj 1–14 + 29 non-paired subj 15–35)',
    color='#fce4ec', edge='#c2185b', fs=8.5)
arrow(2.8, 13.0, 9.0, 12.15, color='gray')
arrow(9.2, 13.0, 9.0, 12.15)

box(9, 10.6, 4.4, 0.8,
    'DESeq2 DEG / fgsea Hallmark: pre-CRT n=33\n(14 pre-CRT in Y + 19 pre-only in N)',
    color='#fce4ec', edge='#c2185b', fs=8.5)
arrow(9, 11.45, 9, 11.0)

box(9, 9.3, 4.4, 0.8,
    'ssGSEA (95 pathways): 56 samples\n'
    'TRUST4 TCR/BCR: 56 samples\n'
    'CMScaller / immune 22 sigs: 56 samples',
    color='#fce4ec', edge='#c2185b', fs=8.5)
arrow(9, 10.2, 9, 9.7)

box(9, 8.0, 4.4, 0.8,
    'Paired pre→post (RNA-derived):\n'
    '12 paired subjects with both timepoints\n'
    '(subj 3 pre-only, 11 post-only in RNA)',
    color='#fce4ec', edge='#c2185b', fs=8.5)
arrow(9, 8.9, 9, 8.4)

box(9, 6.6, 4.4, 0.8,
    'Paired 22-sig Δ (Treg/MHC-II/etc):\n'
    'Both timepoints + matched response = 12\n'
    'For response-stratified Δ: 6 good + 6 bad\n'
    '(subj 5 post-only, 14 pre-only in 22-sig set)',
    color='#fce4ec', edge='#c2185b', fs=8.5)
arrow(9, 7.6, 9, 7.0)

# ---- Integration ----
box(6, 5.0, 5.8, 0.9,
    'INTEGRATED PER-SUBJECT MASTER TABLE\n35 subjects × 37 features\n'
    '(WES-derived 17 + RNA-derived 20); imputation only for missing neoantigen on tumor-only subjects',
    color='#fff8e1', edge='#f57c00', fs=9)
arrow(3, 4.9, 5.3, 5.3)
arrow(9, 6.2, 6.7, 5.3)

# ---- Predictor ----
box(3.5, 3.3, 5.5, 0.9,
    'LASSO logistic regression, features = 37\n'
    'Nested LOOCV (outer 35) + inner 5-fold CV\n'
    'Outer AUC 0.755 (95% CI via bootstrap, perm P reported)',
    color='#e8eaf6', edge='#303f9f', fs=8.5)
arrow(6, 4.55, 4.2, 3.8)

box(8.5, 3.3, 5.5, 0.9,
    'External validation (nCRT meta)\n9 GEO cohorts × 7 signatures\n'
    'N = 721 patients (long-course CRT, rectal)',
    color='#e8eaf6', edge='#303f9f', fs=8.5)
arrow(6, 4.55, 7.8, 3.8)

# ---- Survival note ----
box(6, 1.7, 8.5, 0.8,
    'DFS / OS analysis DEFERRED\n'
    'Clinical outcome data not yet mature in this recently-accrued cohort.\n'
    'TRG-based final TNT response used as endpoint throughout this report.',
    color='#fce4ec', edge='#b71c1c', fs=8.5)

# ---- Caption ----
ax.text(6, 0.5,
    'Fig S1. CONSORT-style sample flow diagram reconciling all per-analysis sample counts. '
    'WES branch (blue) includes 77 sequenced samples collapsing to 49 PASS tumor somatic VCFs; '
    'RNA-seq branch (pink) includes 56 samples collapsing to 33 pre-CRT for DEG/GSEA and '
    '12 paired pre/post for Δ analysis; integration produces a 35×37 subject-level master.',
    ha='center', fontsize=8, color='dimgray', style='italic', wrap=True)

plt.tight_layout()
for ext in ('png','pdf'):
    plt.savefig(FIG/f'SuppFig_consort_sample_flow.{ext}', dpi=300, bbox_inches='tight')
print(f'Wrote: {FIG}/SuppFig_consort_sample_flow.pdf/png')
