# Supplementary Text S1 — Extended Methods

This document supplements the main Methods section with technical detail required for full reproducibility.

## S1.1 Sample collection and processing

Endoscopic biopsies were collected at two timepoints per paired patient (n = 14): pre-treatment (within 14 days before initiation of chemoradiotherapy [CRT]) and post-CRT (after completion of 50.4 Gy radiation but before initiation of consolidation chemotherapy, typically 4–6 weeks after the final radiation fraction). Single-timepoint pre-treatment biopsies were collected for the remaining 21 patients. Whole-blood normal samples were drawn at the pre-CRT visit. Tissue was preserved in RNAlater (Thermo) and snap-frozen within 4 hours of collection.

## S1.2 Whole-exome sequencing

DNA was extracted using AllPrep DNA/RNA Mini (Qiagen). Libraries were prepared with Agilent SureSelect Human All Exon V5 (50 Mb capture). Sequencing on Illumina NovaSeq 6000 generated paired-end 101 bp reads. Median target coverage: tumor 150× (interquartile range 122–176), normal 90× (75–108). Per-sample QC: Picard CollectHsMetrics, FASTQC, Qualimap; samples with mean target coverage <50× were excluded (none).

**Alignment & preprocessing.** BWA-MEM 0.7.17 to GRCh38.p13 (no-alt). MarkDuplicates (Picard 2.27), BaseRecalibrator (GATK 4.6.2) with known sites from gnomAD v3.1, dbSNP build 154, and Mills indels.

**Somatic variant calling.** GATK Mutect2 4.6.2 with a cohort panel-of-normals (PoN) constructed from 28 TNT cohort blood normals (per-normal Mutect2 calls with `--max-mnp-distance 0` and `--germline-resource af-only-gnomad`, joined via GenomicsDBImport + CreateSomaticPanelOfNormals). Tumor-normal pairs (n = 41) were called with `--germline-resource af-only-gnomad`, `--panel-of-normals tnt_cohort_pon.vcf.gz`, `-L exome_targets.bed --interval-padding 100`, `--f1r2-tar-gz` (for read-orientation modelling). Tumor-only calls (n = 8 unmatched tumors) used identical Mutect2 invocation, simply omitting `-I normal -normal`; **all numeric thresholds were left at GATK defaults** — `--max-population-af 0.01`, `--af-of-alleles-not-in-resource 1e-6`, `--tumor-lod-to-emit 3.0`. Per-sample post-call filtering: LearnReadOrientationModel → GetPileupSummaries / CalculateContamination (small_exac_common_3 sites) → FilterMutectCalls (all default thresholds, including `--tumor-lod 2.0` log10, `--threshold-strategy OPTIMAL_F_SCORE`, `--unique-alt-read-count 0`) → `bcftools view -f PASS`. Tumor-only calls were not subject to additional stricter post-filters in this pipeline; their interpretation in driver and TMB tabulations was conservatively flagged in downstream tables. snpEff 5.1d (GRCh38.99) annotation. Final PASS somatic call set: 18,580 variants across 49 tumors.

**Mutational signatures.** SigProfilerAssignment v0.1.6 with COSMIC v3.3 SBS reference, refit mode (no de novo extraction; cohort too small for stable extraction).

**Microsatellite instability.** msisensor-pro v1.2.0 paired mode on 41 matched pairs; baseline computed from gnomAD-derived microsatellite catalog.

**Copy-number variation.** CNVkit v0.9.10 batch mode, pooled reference from blood normals. Segmentation: CBS. HRD proxy scores (LST = large-scale state transitions; TAI = telomeric allelic imbalance; LOH = total LOH segment count) computed from CNV segments using a custom script that mirrors the scarHRD algorithm (scarHRD itself failed installation under the conda environment; results validated by independent re-implementation against scarHRD published examples).

## S1.3 RNA-sequencing

RNA was extracted with AllPrep DNA/RNA Mini (Qiagen) and assessed on Bioanalyzer 2100 (Agilent); RIN ≥ 7.0 required. Stranded total RNA-seq libraries were prepared with TruSeq Stranded Total RNA with Ribo-Zero (Illumina). Sequencing: NovaSeq 6000 PE101, target 40M paired reads per sample.

**Alignment & quantification.** HISAT2 2.2.1 to GRCh38 / GENCODE v39. StringTie 2.2.1 for transcript assembly and quantification with the GENCODE v39 reference annotation. Gene-level TPM matrix: 46,425 protein-coding and lncRNA features × 56 samples.

**Differential expression.** DESeq2 1.42 with the design formula `~ sex + cT_simple + response_bin`. Independent filtering enabled. Wald test, Benjamini–Hochberg adjustment.

**Pathway analysis.** fgsea 1.28 (Hallmark, Reactome MSigDB v2024.1.Hs); gseapy 1.1 ssGSEA on a curated panel of 95 sets (immune, DNA repair, cell cycle, EMT, hypoxia, metabolism). P-values from fgsea adaptive multilevel test that fall below 10⁻¹⁰ are reported as `P < 10⁻¹⁰` to avoid spurious precision.

**Immune deconvolution.** MCP-counter, EPIC, and quanTIseq (immunedeconv 2.1) for cross-method validation; primary results from MCP-counter and ssGSEA-derived signatures (CD8 cytotoxic, activation, proliferation, exhaustion; MHC I/II; NLRC5; TLS Cabrita; TGF-β Mariathasan; EMT Mak; hypoxia Buffa).

**Consensus molecular subtyping.** CMScaller 2.0 with default parameters.

**TCR/BCR repertoire.** TRUST4 1.1.5 with the IMGT human reference. Primary metrics: number of productive clonotypes, Shannon entropy, Gini coefficient, top-clone proportion (per IGH, IGK, IGL, TRA, TRB).

## S1.4 HLA class I typing

OptiType 1.3.5 on bait-extracted HLA reads (BWA-MEM realignment to a 6th-class HLA reference). 4-digit resolution required. Discordant calls between tumor and normal (5/35 subjects) were resolved using the normal-derived call.

## S1.5 Statistical reporting standards

- Continuous comparisons: Mann–Whitney U (two-sided) for between-group; Wilcoxon signed-rank (two-sided) for within-paired-group; BCa bootstrap for confidence intervals on small-n medians.
- Categorical: Fisher exact (two-sided).
- Multiple-testing: Benjamini–Hochberg across feature panels (Table 2; Supp Text S4 for cascade BH FDR).
- Effect sizes reported alongside P-values throughout; for cascade claims (n = 14 paired) the BCa 95% CI is the primary inferential statistic and Mann–Whitney P is an auxiliary check.
