# Tumor DNA Repair Proficiency and Cell-Cycle Activity Predict Response to Total Neoadjuvant Therapy in Locally Advanced Rectal Cancer

## Authors
(박지원, 서울대학교병원, et al.)

## Abstract (draft ~250 words)

**Background.** Total neoadjuvant therapy (TNT) is standard for locally advanced rectal cancer (LARC), but response is heterogeneous. Molecular predictors of response in MSS-dominant LARC remain ill-defined.

**Methods.** We performed whole-exome sequencing (N=77, 28 matched normals) and RNA-seq (N=56) on tumors from 35 LARC patients pre- and post-TNT. Patients were categorized by TRG as good responders (TRG0–1; n=18) or poor (TRG2–3; n=17). Somatic variants were called with GATK Mutect2 using a cohort panel-of-normals and matched tumor-normal pairs (unmatched n=8 used tumor-only+PoN). MSI was assessed with msisensor-pro; mutational signatures by SigProfilerAssignment against COSMIC v3.3; CNV by CNVkit; HLA class I by OptiType. Transcriptomic analyses included DESeq2 (adjusting for sex and cT), fgsea Hallmark/Reactome GSEA, and ssGSEA. Integration used Spearman correlation and Mann-Whitney association testing on 37 integrated per-subject features.

**Results.** All 41 matched tumors were microsatellite stable (MSI <0.2%) and carried low TMB (median 1.6/Mb). Mutational signatures were dominated by ageing (SBS5/1), with sporadic MMR-related SBS15 contributions in both responder strata; SBS3 was absent. CRC driver mutations (APC 30/49, TP53 20/49, KRAS 14/49) showed classic MSS-CRC patterns. **Good responders exhibited markedly higher pre-treatment expression of DNA double-strand break repair (p=0.007), homology-directed repair, general DNA repair, E2F/Myc/G2-M cell cycle programs, and a CD8 proliferation signature (all p<0.05), with concomitant lower epithelial-to-mesenchymal transition (EMT) pathway activity (GSEA NES=−2.16).** Hallmark GSEA identified proliferation/DNA-repair upregulation (E2F_TARGETS NES=2.78, p=8×10⁻²⁶) and mesenchymal/stromal downregulation as the most pathway-robust signal. CNV burden, MSI, and CMS subtype did not discriminate responders.

**Conclusion.** In this MSS LARC cohort, response to TNT is primarily governed by tumor-intrinsic DNA-repair proficiency and proliferative capacity, not by MSI/TMB/neoantigenicity or classical CMS stratification. EMT activation marks a chemoradioresistant state. These transcriptomic features provide a practical response predictor and nominate candidates for response-adapted TNT de-escalation or alternative neoadjuvant strategies.

---

## 1. Introduction (~800 words)

- LARC + TNT standard of care (Dutch TRIGGER, OPRA, PRODIGE 23)
- Heterogeneous response — ~15-20% pathologic complete response
- Prior molecular predictors: CEA, inflammation, CMS subtype, immune signatures (mostly underpowered)
- MSI-H LARC responds dramatically to ICB (Cercek 2022 NEJM) — but MSS majority lacks clear biomarker
- Knowledge gap: **MSS LARC response biology is unresolved**
- Our study: matched WES+RNA-seq of 35 patients, proper somatic calling (Mutect2 T-N), integrated multi-omics

## 2. Methods (~1500 words)

### 2.1 Patients and cohort
- 35 LARC patients, Seoul National University Hospital
- TNT regimen: neoadjuvant chemo (FOLFOX or CAPOX) + long-course chemoradiation
- TRG grading 0-3 (CR/near-CR vs PR/poor): good=18, bad=17
- Sample structure: 14 paired pre+post+normal, 15 paired normal+pre, 6 pre-only

### 2.2 WES
- Agilent SureSelect V5 (50Mb), Illumina paired-end 101bp, BWA-MEM → BQSR
- Somatic calling: GATK Mutect2 + cohort PoN (28 normals) + gnomAD germline resource
- T-N for 41 matched tumors; tumor-only+PoN for 8 unmatched
- FilterMutectCalls + LearnReadOrientationModel + contamination correction
- Annotation: snpEff GRCh38.99
- Mutational signatures: SigProfilerAssignment, COSMIC v3.3, SBS96 refit
- MSI: msisensor-pro T-N
- CNV: CNVkit batch mode with cohort reference
- HLA class I typing: OptiType on MHC-region-extracted reads

### 2.3 RNA-seq
- HiSeq/NovaSeq paired-end, HISAT2 → StringTie → Gene-level TPM
- 46,425 unique symbols × 56 samples
- DEG: DESeq2, design ~ sex + cT_simple + response_bin, pre-samples only
- GSEA: fgsea Hallmark + Reactome (msigdbr)
- ssGSEA: gseapy, Hallmark + curated Reactome immune/repair/cycle sets
- Signature scoring: mean z-score across member genes (CD8 prolif/activation/exhaustion, MHC, NLRC5-HLA-IFNG, TLS Cabrita, TGFβ Mariathasan, EMT Mak, Hypoxia Buffa)
- CMS subtype: CMScaller (Entrez input, log2TPM+1)

### 2.4 Statistics & Integration
- Mann-Whitney for continuous features, Fisher exact for categorical
- BH FDR for multi-feature correction
- Spearman correlation for feature interaction
- Per-subject integrated master table (37 features) for response association

---

## 3. Results

### 3.1 Cohort characteristics
- 35 patients, 18 good / 17 bad
- cT4 enriched in bad (41% vs 11%, p=0.086)
- No significant age/sex difference

### 3.2 The cohort is microsatellite stable and low-TMB
- MSI-H: 0/41; all MSI <0.2%
- TMB (matched): median 1.6/Mb, trend higher in good (1.85 vs 1.40, p=0.186)
- **Key message: Classical ICB-response biomarkers (MSI/TMB) do not apply**

### 3.3 Mutational signatures are dominated by ageing
- SBS5 + SBS1 account for >60% of mutations
- Sporadic MMR signature (SBS15/SBS6) in both good and bad; MMR proportion not associated with response
- SBS3 (HRD) absent from cohort
- CRC driver mutations follow classical pattern (APC 61%, TP53 41%, KRAS 29%); FBXW7 mutation trend ↑ in good (OR=3.7, p=0.36)

### 3.4 Copy number burden does not discriminate response
- CIN median: good 0.20 vs bad 0.23 (p=0.659, pre)
- No significant amp/del fraction difference

### 3.5 **Pre-treatment DNA repair and cell-cycle activity predict response** (main finding)
- Hallmark GSEA: E2F_TARGETS NES=2.78 (p=8×10⁻²⁶), G2M_CHECKPOINT NES=2.46, MYC_TARGETS_V1 NES=2.36 → UP in good
- EMT NES=−2.16 → DOWN in good
- Reactome: Cell_Cycle_Checkpoints, M_Phase, **HDR (p=10⁻⁹)**, DNA_Double_Strand_Break_Repair UP in good
- ssGSEA: DSB Repair (p=0.007), HDR (p=0.020), Myc Targets V2 (p=0.018), E2F (p=0.035)
- **Biological interpretation**: Proliferating, DNA-repair-proficient tumors respond to chemoradiation; EMT-high mesenchymal tumors resist.

### 3.6 CD8 proliferation signature confirms immune correlate
- CD8 proliferation (MKI67/TOP2A/CCNB etc.) ↑ good (Δz = +0.50, p=0.035)
- MHC_II paradoxically ↓ good (p=0.074) — suggesting non-classical immune interaction
- Post-treatment: all immune axes (IFNγ, Treg, checkpoint, exhaustion) elevated in good responders → consistent with treatment-induced immune activation of responders

### 3.7 Integrated predictor
- 37-feature integrated table → DSB repair, HDR, E2F, Myc, G2M, CD8 prolif all point same direction
- Spearman correlation heatmap shows tight cluster of DNA-repair + cell-cycle features; separate cluster for EMT/TGFβ (low in good)

### 3.8 CMS and HLA
- CMS distribution not associated with response (Fisher p=1, pre)
- HLA class I homozygosity not associated (p=0.31)

---

## 4. Discussion (~1500 words)

### Key findings
- MSS LARC response is **not** driven by classical ICB biomarkers
- It is driven by tumor-intrinsic DNA-repair/cell-cycle state
- Proliferating + DNA-repair-proficient + non-EMT = good responder phenotype
- This explains why chemoradiation works: it exploits active cell cycle and functioning apoptotic/repair machinery

### Contrast with prior dogma
- ICB dogma: high TMB + MSI-H + inflamed → good response
- TNT in MSS: opposite — intact DNA repair helps responding to induced damage
- Highlights that neoadjuvant chemoradiation and ICB operate through distinct molecular logics

### Clinical implications
- Pre-treatment RNA-seq transcriptomic classifier (DSB/HDR/E2F/Myc/CD8_prolif) could guide response prediction
- EMT-high tumors candidates for intensified/alternative neoadjuvant (e.g., add taxane, anti-TGFβ)
- MSI-H subset (not present in this cohort but relevant) should still receive dostarlimab per Cercek 2022

### Limitations
- N=35, single-center
- Unmatched WES samples (8) bias somatic calling; interpret with care
- No HLA LOH (LOHHLA pending) or neoantigen burden (pVACtools pending) at this draft
- Need external validation in LARC TNT cohorts (GSE150082, others)

### Future directions
- Validate DSB/HDR/E2F signature in external LARC TNT cohort
- Prospective trial stratifying by the signature
- Single-cell RNA-seq subtype-specific TME analysis

---

## 5. Figures/Tables plan

### Figures
- **Fig 1**: Cohort overview (flowchart + Table 1 graphic + response distribution)
- **Fig 2**: WES landscape — oncoprint (drivers) + TMB distribution + MSI status + SBS signatures bar
- **Fig 3**: RNA immune axes — signature boxplots + heatmap (pre/post × good/bad)
- **Fig 4**: **Main result** — GSEA Hallmark barplot + ssGSEA top pathways + DSB/HDR boxplots by response
- **Fig 5**: Integrated response association barplot (−log10p ranked) + correlation heatmap
- **Fig 6**: Model — "proliferative/DNA-repair-proficient vs EMT-resistant" schematic

### Tables
- **Table 1**: Clinical characteristics (good vs bad)
- **Table 2**: Integrated feature vs response — top 20 with p/q values
- **Supp Table S1**: All 37 features per subject
- **Supp Table S2**: SBS signature activities per sample
- **Supp Table S3**: Full GSEA Hallmark + Reactome results
- **Supp Table S4**: CRC driver mutations per sample
- **Supp Table S5**: HLA class I types per subject

---

## 6. Outstanding analyses to complete before submission

1. **HLA LOH (LOHHLA)** — waiting on pipeline setup
2. **Neoantigen prediction (pVACtools + NetMHCpan)** — queued
3. **HRD genomic scar score (scarHRD)** — queued
4. **TCR repertoire (TRUST4)** — queued
5. **External validation**: GSE150082, GSE190826 (neoadjuvant CRT cohorts)
6. **ML predictor**: elastic net / RF with LOOCV using integrated features
7. **Mediation analysis**: HDR → CD8 prolif → response causal path test
8. **Clonal evolution (PyClone-VI)** for 13 paired pre+post subjects — treatment-sensitive clones

---

## Data/code availability
- GitHub repo: TBD (per user GitHub-sync preference; create TNT_manuscript repo)
- Raw data: Macrogen HN00249207 (WES) + HN00249209 (RNA-seq), matched metadata
- Analysis workspace: `/mnt/sda1/data/TNT/analysis/` on local server

---

*Draft v0.1 — 2026-04-14. Subject to revision after completing outstanding analyses.*
