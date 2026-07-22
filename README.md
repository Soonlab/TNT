# TNT — Multi-omics analysis code

Analysis code for a multi-omics study of **total neoadjuvant therapy (TNT) response in
microsatellite-stable rectal cancer**. Patients received short-course radiotherapy
(SC-RT, 25 Gy / 5 fractions) followed by consolidation FOLFOX/CAPOX (N = 35; WES + RNA-seq,
with paired pre-/post-radiotherapy biopsies in a subset).

Each top-level folder corresponds to a subsection of the **Methods** in the manuscript, and
contains only the final version of the analysis code for that subsection.

Sequencing data are deposited under controlled access (SRA/GEO); accessions
are listed in the manuscript.

## Folders (mapped to Methods)

| Folder | Methods subsection |
|--------|--------------------|
| — *(no code; clinical metadata only)* | Patients, samples, and response grading |
| `01_wes_somatic_variant_calling/` | Whole-exome sequencing and somatic variant calling |
| `02_signatures_cnv_hrd_hla_neoantigen/` | Mutational signatures, copy number, HRD proxies, and HLA/neoantigen analysis |
| `03_rnaseq_deg/` | RNA-seq alignment, quantification, and differential expression |
| `04_pathway_immune_repertoire/` | Pathway, immune-signature, and repertoire analysis |
| `05_nanostring_immune/` | NanoString nCounter immune profiling |
| `06_integrated_classifier/` | Integrated classifier |
| `07_paired_radiation_phase/` | Paired radiation-phase analyses |
| `08_external_meta_analysis/` | External meta-analysis across nCRT cohorts |
| `09_statistical_considerations/` | Statistical considerations (sensitivity analyses) |

## Contents

- **01 — WES somatic calling:** pseudo-somatic extraction, reference prep, GATK Mutect2,
  snpEff annotation, variant master / TMB / drivers, MSI (msisensor-pro).
- **02 — Signatures / CNV / HRD / HLA / neoantigen:** SBS refit (SigProfilerAssignment),
  CNVkit, scarHRD, OptiType HLA typing, HLA-LOH (lite + stricter), pVACseq neoantigens
  (with VEP), neoantigen summary.
- **03 — RNA-seq DEG:** TPM + signature build, DESeq2 differential expression + fgsea,
  CMS classification (CMScaller).
- **04 — Pathway / immune / repertoire:** ssGSEA (GSVA), TRUST4 BCR/TCR reconstruction,
  purified immune-signature scoring, IGH V-gene directional coherence.
- **05 — NanoString:** pre-specified primary/secondary/tertiary tests and exploratory analysis
  of the extreme-phenotype 6-patient subset.
- **06 — Integrated classifier:** multi-omic master-table integration, extended nested-LOOCV
  with permutation, drop-vs-swap ablation, and per-fold feature-stability sweep
  (28-feature elastic-net reference).
- **07 — Paired radiation-phase analyses:** paired pre→post Δ, purity-adjusted Δ, BH-FDR across
  the cascade family, BCa bootstrap CIs, paired immune Δ, baseline-vs-cascade convergence tests,
  and baseline-factor pharmacodynamics.
- **08 — External meta-analysis:** per-cohort reproducibility, restricted 5-cohort nCRT meta,
  final CD8-cytotoxic meta incorporating Akiyoshi et al. (GSE216616, n=298), the v3 CD8-axis
  re-analysis, drop-cohort sensitivity, and the GSE254249 SC-RT validation.
- **09 — Statistical considerations:** sensitivity analysis excluding the 8 tumor-only WES samples.

## Environment

Python via miniconda; R with Bioconductor (R 4.3.x / Bioconductor 3.18). Key external tools:
BWA, HISAT2, StringTie, GATK/Mutect2, snpEff, CNVkit, scarHRD, msisensor-pro,
SigProfilerAssignment, DESeq2, fgsea, GSVA, CMScaller, TRUST4, OptiType, pVACseq, MHCflurry.
Exact tool versions are recorded in the manuscript's supplementary methods.
