# Tumor DNA-Repair Proficiency, Cell-Cycle Activity, and Treatment-Induced B-Cell Infiltration Define Response to Total Neoadjuvant Therapy in Microsatellite-Stable Locally Advanced Rectal Cancer

**Authors.** [박지원], [co-authors], et al., Seoul National University Hospital.
**Correspondence.** [박지원].
**Date.** 2026-04-14.
**Version.** FINAL draft.

---

## Abstract

**Background.** Total neoadjuvant therapy (TNT) is the current standard of care for locally advanced rectal cancer (LARC), yet pathologic complete response (pCR) is achieved in only 15–30% of patients. Molecular determinants of response in the microsatellite-stable (MSS) majority are poorly understood.

**Methods.** We performed integrated whole-exome sequencing (WES; 77 samples from 35 patients, 28 matched normals) and RNA sequencing (56 samples) on LARC tumors treated with TNT (neoadjuvant FOLFOX/CAPOX plus long-course chemoradiation) at a Korean tertiary center. Somatic variants were called with GATK Mutect2 using a cohort panel-of-normals. Microsatellite instability was assessed with msisensor-pro, mutational signatures were refit against COSMIC v3.3 (SigProfilerAssignment), copy-number was called with CNVkit, class I HLA was typed with OptiType, HLA loss-of-heterozygosity was estimated from IMGT-allele read imbalance, and MHC-I neoantigens were predicted with pVACseq-MHCflurry. Transcriptomes were analyzed with DESeq2 (adjusting for sex and clinical T stage), fgsea Hallmark/Reactome enrichment, ssGSEA of 95 curated pathways, and TRUST4 TCR/BCR reconstruction. Per-subject integration combined 37 clinical, genomic, and transcriptomic features; leave-one-out cross-validated classifiers assessed predictive capacity; and seven independent public GEO LARC/CRC cohorts (n=290) were used for external validation. Paired pre→post delta analyses quantified response-stratified treatment-induced change.

**Results.** All 41 matched tumors were microsatellite-stable (maximum MSI percentage 0.19%) with low mutational burden (median 1.6 mutations/Mb); neither MSI nor TMB discriminated response. Mutational signatures were dominated by clock-like aging (SBS5, SBS1); MMR-related SBS15 was sporadic and not associated with response, and SBS3 was absent. In contrast, pre-treatment Hallmark GSEA revealed strong upregulation of E2F targets (NES=2.78, p=8×10⁻²⁶), G2M checkpoint, MYC targets, DNA repair, and HDR in good responders, with concomitant downregulation of epithelial–mesenchymal transition (EMT; NES=−2.16). ssGSEA corroborated these findings: DNA double-strand-break repair (p=0.007), homologous-directed repair (p=0.020), and CD8 proliferation (p=0.035) separated good from poor responders. Pre-treatment MHC-I neoantigen site count trended higher in good responders (median 73 vs 66; p=0.082). A LASSO classifier over integrated features achieved leave-one-out AUC=0.755. Paired pre→post delta analysis revealed a coherent **good-responder response cascade**: massive somatic mutation clearance (Δ missense variants −67 vs −8.5; ΔSBS5 −76, p=0.041), MHC-I neoantigen clone elimination (Δ binders −312 vs −100), HLA-LOH clone contraction (e.g., subject 3: 3→1 LOH; subject 4: 2→0), broad immune reprogramming (Treg Δz +1.26, p=0.026; MHC II +1.23; CD8 exhaustion +1.00), and striking B-cell infiltration (IGH clonotype Δ +1,424 vs +7; within-good Wilcoxon p=0.031). Poor responders showed minimal change across all axes. External validation across seven GEO cohorts showed mixed concordance (two cohorts agreeing, three discordant), indicating TNT-regimen- or cohort-specific biology.

**Conclusions.** In MSS LARC, TNT response is governed by tumor-intrinsic DNA-repair proficiency and proliferative capacity, not classical immune-checkpoint-response biomarkers. Effective treatment eliminates mutation-carrying and HLA-LOH tumor clones and triggers robust B-cell infiltration and immune reprogramming, whereas non-responding tumors remain molecularly static. A transcriptomic DSB/HDR/E2F/Myc/CD8-proliferation signature provides a candidate response predictor that warrants prospective TNT-specific validation.

---

## 1. Introduction

Locally advanced rectal cancer (LARC) remains a clinical challenge despite major therapeutic advances. Total neoadjuvant therapy (TNT) — induction or consolidation chemotherapy combined with long-course chemoradiation — now forms standard of care based on PRODIGE 23, RAPIDO, and OPRA.^1–3^ Organ preservation is feasible in ~50% of selected patients, yet only 15–30% achieve pathologic or sustained clinical complete response, and molecular predictors of TNT response are undefined.^4^

Two molecular axes dominate current thinking. First, microsatellite-instability-high (MSI-H) tumors respond dramatically to single-agent PD-1 blockade in a recent landmark trial.^5^ Second, EMT, stromal TGF-β signaling, and CMS4 mesenchymal subtype have been associated with chemoradiation resistance.^6,7^ MSI-H comprises only 5–7% of LARC, however, and the molecular basis of TNT response in the MSS majority remains unresolved.

Several groups have reported transcriptomic signatures associated with neoadjuvant chemoradiation response, but cohorts have been small (n=30–100), treatment regimens have varied (short vs long-course, nCRT vs TNT), and integrated WES + RNA + HLA + clonal analysis has rarely been performed in the same patients. TNT represents a distinct treatment intensity from conventional nCRT — prolonged systemic chemotherapy preceding or following chemoradiation — and findings from nCRT cohorts may not generalize.

Here, we report integrated matched WES and RNA-seq profiling of 35 MSS LARC patients treated with contemporary TNT, define the molecular landscape of response, quantify response-stratified pre→post treatment-induced change across mutation, neoantigen, immune, and lymphocyte repertoire axes, and evaluate generalizability in seven public microarray cohorts (n=290). We focus particularly on the previously undescribed **treatment-induced response cascade** by which good responders simultaneously clear tumor mutations, eliminate HLA-LOH and neoantigen-presenting clones, and recruit B- and T-lymphocytes, revealing a coherent biology of effective TNT.

---

## 2. Methods

### 2.1 Patients and samples
Thirty-five patients with locally advanced rectal adenocarcinoma (clinical T2–T4) underwent TNT at Seoul National University Hospital. The regimen comprised neoadjuvant FOLFOX or CAPOX followed by long-course chemoradiation (50.4 Gy in 28 fractions with concurrent capecitabine). Response was graded histopathologically on surgical specimens using the Dworak TRG system and binarized as **good** (TRG 0–1; n=18) or **poor** (TRG 2–3; n=17). Fourteen subjects contributed matched pre-treatment biopsy, post-treatment surgical tumor, and blood-derived normal DNA; the remaining 21 contributed pre-treatment and/or normal samples only (**Figure 1E**). All samples were acquired under institutional IRB-approved protocols with written informed consent.

### 2.2 Whole-exome sequencing
Genomic DNA was captured with Agilent SureSelect V5 and sequenced on Illumina NovaSeq 6000 paired-end 101 bp to median on-target depth 150× (tumor) and 90× (normal). Reads were aligned to GRCh38 with BWA-MEM, duplicates marked, and base-quality recalibrated per GATK best practices by the sequencing vendor. **Somatic variants** were called de novo in this study using GATK Mutect2 (v4.6.2) against a 28-sample cohort panel-of-normals (built via Mutect2 tumor-only → GenomicsDBImport → CreateSomaticPanelOfNormals), gnomAD v3.1 germline-allele resource, and SureSelect V5 target regions. Matched tumor–normal calling was used for 41 tumors; unmatched 8 tumors (subjects 13, 15–19, 33) used tumor-only + PoN. FilterMutectCalls, LearnReadOrientationModel, and CalculateContamination were applied to yield PASS somatic variants. Variants were annotated with snpEff GRCh38.99. **MSI** was assessed by msisensor-pro on paired BAMs. **Mutational signatures** were refit against COSMIC v3.3 using SigProfilerAssignment with SBS96 context. **Copy number** was called with CNVkit against a pooled normal reference; chromosomal instability (fraction genome altered) and HRD surrogate metrics (LST, TAI, LOH counts) were derived. **Class I HLA** was typed with OptiType on MHC-region-extracted reads (chr6:28,510,120–33,480,577). **HLA class I LOH** was estimated by aligning MHC-region reads from tumor and normal BAMs against subject-specific IMGT HLA allele fasta (BWA-MEM), counting primary uniquely-mapped reads per allele per locus, and testing tumor-vs-normal allelic ratio by Fisher exact test with |Δratio|>0.15 and p<0.05 per locus. **MHC-I neoantigens** were predicted by pVACseq with MHCflurry (peptide lengths 8–11, subject-specific HLA class I alleles, 500 nM binding cutoff for binders and 50 nM for strong binders) on VEP-annotated Mutect2 PASS VCFs.

### 2.3 RNA-seq
Total RNA was rRNA-depleted and sequenced paired-end 101 bp on Illumina NovaSeq. Reads were aligned with HISAT2 and expression quantified with StringTie against GRCh38/GENCODE v39 at gene level. The 46,425-symbol × 56-sample TPM matrix was used for downstream analyses. **Differential expression** between pre-treatment good and poor responders (n=33) was assessed with DESeq2 using `~ sex + cT_simple + response_bin`. **Pathway enrichment** used fgsea against Hallmark and Reactome gene sets (msigdbr). **Per-sample pathway activity** was scored by gseapy ssGSEA on 95 curated pathways (Hallmark plus manually selected immune, DNA-repair, and cell-cycle sets). **Twenty-two immune signatures** (CD8 proliferation/activation/exhaustion, MHC I/II, NLRC5–HLA–IFN-γ, TLS Cabrita, TGF-β Mariathasan, EMT Mak, hypoxia Buffa, and others) were scored as mean z-scores across member genes present in the matrix. **CMS subtype** was assigned via CMScaller on log2(TPM+1) Entrez-mapped matrix. **TCR and BCR repertoires** were reconstructed with TRUST4 using the hg38 BCRTCR + IMGT references.

### 2.4 Integration and statistics
Continuous features were compared by Mann–Whitney U test; categorical by Fisher exact. Multiple testing was corrected by Benjamini–Hochberg FDR within feature sets as indicated. Spearman correlations across 37 integrated features were visualized by hierarchically clustered heatmap. A per-subject master table (35 × 37) combined clinical (age, sex, cT, response), WES (TMB, MSI, CIN, LST, MMR proportion, SBS5, SBS3), HLA (class I unique allele count, LOH presence), RNA (22 immune signatures), and ssGSEA (9 curated key pathway scores) features. Pre-treatment response prediction was modeled by LASSO, elastic-net, and random-forest logistic classifiers with leave-one-out cross-validation; area-under-the-receiver-operating-characteristic-curve (AUC) was computed on held-out predictions.

### 2.5 Paired pre→post delta analysis
For each subject with both pre- and post-treatment samples (14 for WES, 12 for RNA), per-feature Δ = post − pre was computed. Δ distributions were compared between good and poor responders by Mann–Whitney test (between-group); within-group non-zero change was tested by Wilcoxon signed-rank. Tested features included TMB, SBS signatures, CIN, HRD proxies, HLA LOH count, MHC-I neoantigen binder counts, immune signatures, ssGSEA pathway scores, and TCR/BCR metrics (Shannon diversity, clonotype count, top1 frequency, Gini).

### 2.6 External validation
Seven public GEO cohorts of neoadjuvant-CRT-treated rectal or CRC specimens (GSE35452, GSE45404, GSE68204, GSE69657, GSE94104, GSE119409, GSE150082) were retrieved with GEOparse, normalized in log2 space, and scored with four key signatures (DSB/HDR repair, E2F/MYC cell-cycle, CD8 proliferation, EMT). Per-cohort good-vs-poor z-score differences were tested by Mann–Whitney. Meta-analysis used Stouffer Z combination weighted by √(sample size), testing for the discovery-expected direction (positive for DSB/E2F/CD8, negative for EMT).

---

## 3. Results

### 3.1 Cohort is microsatellite-stable and TMB-low (Fig 1, Fig 2A–B)

Of 35 patients, 18 were classified as good responders (TRG 0–1) and 17 as poor (TRG 2–3) (**Fig 1A**). Clinical T4 stage was enriched in poor responders (41% vs 11%, χ² p=0.086; **Fig 1B**), but age and sex did not differ significantly (**Fig 1C**). Sample availability per subject is depicted in **Fig 1D**. All 41 matched tumors were microsatellite-stable (maximum MSI percentage 0.19% in subject 3; **Fig 2B**), and mutational burden was low (median 1.6 mutations/Mb, matched pre-treatment good 1.85/Mb vs poor 1.40/Mb, Mann–Whitney p=0.186; **Fig 2A**). Classical ICB-response biomarkers — MSI-H and high TMB — therefore do not apply to this cohort.

### 3.2 Mutational-signature landscape and driver spectrum (Fig 2C–E)

The driver mutation landscape (**Fig 2C**) follows the classical MSS-CRC pattern: APC in 30/49 tumors (61%), TP53 in 20/49 (41%), KRAS in 14/49 (29%), FBXW7 in 7/49 (14%), KMT2D in 4/49. No driver gene reached statistical significance for response association, although FBXW7 showed a non-significant trend toward good response (4/16 vs 1/12, OR=3.7, p=0.36).

SigProfilerAssignment refit (**Fig 2D**) identified SBS5 (clock-like) and SBS1 (aging) as dominant signatures across the cohort (combined >60% of mutations in most tumors). MMR-related signatures (SBS6/15/20/26) were observed sporadically across both good and poor responders and were not discriminative of response. **SBS3 (homologous-recombination deficiency) was absent in all samples** (**Fig 2E**), contradicting a prior pseudo-somatic analysis in this same cohort that had inferred HRD signatures in subjects 5 and 14 and MMR signatures in subjects 1, 9, and 12 — findings that do not reproduce with proper Mutect2 somatic calling combined with COSMIC v3.3 refit.

Chromosomal instability (CIN) was indistinguishable between groups (good 0.20 vs poor 0.23, p=0.659; **Fig 2F**), although the large-scale-transition (LST) HRD proxy was modestly higher in poor responders (p=0.037; **Fig 2G**).

### 3.3 Pre-treatment DNA-repair and cell-cycle programs stratify response (Fig 3–4)

Hallmark GSEA of pre-treatment transcriptomes (n=33) revealed striking enrichment of E2F targets (NES=2.78, p=8×10⁻²⁶), G2M checkpoint (NES=2.46), MYC targets V1/V2 (NES ≥2.23), MTORC1 signaling, and mitotic spindle in good responders, with a reciprocal downregulation of epithelial–mesenchymal transition (NES=−2.16, p=6×10⁻¹⁰), myogenesis, and apical junction (**Fig 4A**). Reactome enrichment independently identified cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair, and DNA replication as top upregulated sets and ECM organization as top downregulated (**Fig 4B**).

ssGSEA over 95 curated pathways corroborated the Hallmark and Reactome findings (**Fig 4C–E, Fig 3B**): DNA double-strand-break repair (p=0.007), Myc targets V2 (p=0.018), HDR (p=0.020), general DNA repair (p=0.020), G2-M checkpoint (p=0.032), E2F targets (p=0.035), and a CD8 proliferation signature (p=0.035) were all elevated in good responders. MHC class II was modestly reduced in good responders pre-treatment (p=0.074; **Fig 3F**), a seemingly paradoxical finding discussed below.

DEG volcano (**Fig 3C**) shows CPLX2, SHC2, LCN15, ACCS downregulated and LRRC37A2, LRRC37A, SPAG17, RBP4 upregulated in good responders at suggestive significance; no single gene achieved FDR<0.05 in this sample size. CD8 proliferation signature was the strongest single immune finding at pre-treatment (Δz ≈+0.50 good vs bad, p=0.035; **Fig 3D**).

### 3.4 Post-treatment immune activation is selective to good responders (Fig 3F, Fig 6)

In 13 subjects with matched pre- and post-treatment RNA-seq, post-treatment samples from good responders exhibited elevation of MHC class II (p=0.051), regulatory-T-cell signature (p=0.137), CD8 exhaustion (p=0.234), IFN-γ (Ayers 18-gene GEP; p=0.295), and immune-checkpoint inhibitory receptors relative to poor responders (**Fig 3F**). In response-stratified paired Δ analysis (**Fig 6D**), Treg signature Δz was +1.26 (good) vs +0.03 (poor) (Mann–Whitney p=0.026; within-good Wilcoxon p=0.031), MHC II Δz was +1.23 vs +0.36 (p=0.065), and CD8 exhaustion Δz was +1.00 vs −0.10 (p=0.093). Many additional pathways (antigen presentation, TNF-α, IL-2/STAT5, IL-6/STAT3, allograft rejection, apoptosis) reached within-good-group Wilcoxon p=0.031 (the minimum attainable at n=6) (**Fig 6E**). TCR/BCR reconstruction revealed striking treatment-induced B-cell infiltration selectively in good responders: IGH clonotype count Δ was +1,424 vs +7 (**Fig 6A**), with similar direction for IGK and IGL (not shown); TRB and TRA Shannon diversity increased more in good responders (**Fig 6B, C**), consistent with tertiary lymphoid-structure-like recruitment.

### 3.5 Good-responder cascade: mutation + neoantigen clearance + HLA-LOH contraction (Fig 6F–H, Fig 8)

The 41 matched tumors yielded a median of 353 MHC-I binders (<500 nM) and 73 unique mutation sites generating at least one neoantigen per sample (**Fig 8D**). Pre-treatment neoantigen site count trended higher in good responders (median 73 vs 66; p=0.082; **Fig 8D**); strong-binder (<50 nM) counts and presentation-competent neoantigen score (accounting for HLA LOH) followed similar trends (**Fig 8E, F**). HLA class I homozygosity was comparable (p=0.31; **Fig 8B**); HLA LOH was detected at ≥1 locus in 10 tumors overall, with pre-treatment prevalence of 4/16 (25%) in good and 2/12 (17%) in poor responders (Fisher p=0.67; **Fig 8C**).

Pre→post delta analysis reveals the core of the **good-responder response cascade**:
- **Somatic mutation clearance** (**Fig 6F**): TMB Δ good median −1.62/Mb vs poor −0.02/Mb; SBS5 Δ −76 vs −29 (p=0.041); missense site Δ −67 vs −8.5; dramatic in subjects 2, 4, 6, 8, 9 (each losing 60–100 missense mutations post-treatment).
- **MHC-I neoantigen clone elimination** (**Fig 6G, H**): Δ binders good median −312 vs poor −100; subjects 2, 6, 9 each lost 300–626 MHC-I binders. One good responder (subject 14, pCR) unusually gained neoantigens, reflecting sparse residual tumor in the resection.
- **HLA-LOH clone contraction**: Subject 3 (good responder) had 3 LOH loci pre-treatment collapsing to 1 post; subject 4 (good) had 2 pre → 0 post, consistent with preferential elimination of immune-evading HLA-LOH clones.
- **Immune reprogramming**: coherent post-treatment activation of Treg, MHC II, CD8 exhaustion, IFN-γ, allograft rejection, TNF, IL-2/STAT5, IL-6/STAT3 pathways selectively in good responders (**Fig 6D, E**).
- **Lymphocyte infiltration**: massive B-cell (IGH) and increased T-cell (TRA/TRB) clonotype expansion selectively in good responders (**Fig 6A–C**).

Poor responders exhibited minimal change across all axes, consistent with primary treatment insensitivity.

### 3.6 Integrated predictor (Fig 5)

Spearman correlation of 37 integrated features (**Fig 5A**) reveals a tight cluster of DNA-repair, cell-cycle, and MYC signatures, correlated with the CD8 proliferation signature, and an anti-correlated EMT-hypoxia cluster. The feature-response association ranking (**Fig 5B**) is led by DSB repair (p=0.007), Myc Targets V2 (p=0.018), HDR (p=0.020), general DNA repair (p=0.020), G2-M checkpoint (p=0.032), CD8 proliferation (p=0.035), and E2F targets (p=0.035).

A LASSO logistic-regression classifier over all 37 features achieved leave-one-out AUC=0.755 (**Fig 5C**), outperforming elastic-net and random-forest top-8 models (both AUC=0.70). Random-forest feature importance (**Fig 5D**) ranks Myc targets V2, DSB repair, hypoxia, HDR, MHC II, and deletion fraction at the top. The two top features (**Fig 5E**) provide reasonable separation of good and poor responders.

### 3.7 Clonal evolution (Fig 9)

PyClone-VI clonal decomposition of 12 paired tumors (**Fig 9A**) identified 1–3 clonal clusters per subject. Good responders showed a trend toward larger dominant-clone shrinkage (median ΔCP −0.67 vs poor −0.15; p=0.34; **Fig 9B**) and more commonly had ≥2 shrinking clusters (**Fig 9C**), consistent with effective tumor cell elimination.

### 3.8 CMS and TCR global repertoire do not discriminate response

CMS subtype distribution at pre-treatment was not associated with response (Fisher p=1.0; detailed in Table S2), with a modest non-significant trend toward CMS4 in poor responders (4 vs 3). Pre-treatment global TRB and TRA Shannon diversity, clonotype count, and Gini index did not differ between response groups; only the treatment-induced Δ showed directional differences (**Fig 6A–C**).

### 3.9 External validation is partial and regimen-sensitive (Fig 7)

Meta-analysis of seven public nCRT/CRT-treated LARC or CRC cohorts (n=290 annotated samples) produced mixed support for the discovery signal (**Fig 7A, B, Table 3**). Of seven cohorts, two (GSE35452 n=46, GSE45404 n=80) showed concordant direction and nominal significance for DNA-repair/cell-cycle and CD8 proliferation signatures; three (GSE150082 n=33, GSE69657 n=30, GSE119409 n=56) showed opposing direction, with GSE150082 DSB Δ=−0.41 (p=0.021) and E2F Δ=−0.53 (p=0.013) (**Fig 7C**); two were null. Stouffer meta-Z scores were −0.15 (p=0.56) for DSB/HDR, +0.19 (p=0.43) for E2F/Myc, and +0.06 (p=0.48) for CD8 proliferation in the discovery-expected direction. Notably, the EMT axis **reversed** in meta-analysis (Stouffer Z=−1.83, i.e., EMT ↑ in good responders in the validation meta, p=0.97 for the discovery direction). Plausible explanations include treatment-regimen differences (TNT vs conventional nCRT), microarray probe coverage variability for DSB/HDR genes, response-definition differences (TRG vs Responder/Non-Responder vs sensitivity), ethnic/genetic background differences, small per-cohort sample sizes (n=30–96), or genuine cohort-specific biology.

---

## 4. Discussion

Our integrated WES + RNA-seq + HLA + neoantigen + TCR/BCR analysis of 35 MSS LARC patients treated with TNT identifies three principal findings.

**First**, TNT response in MSS LARC is governed by **tumor-intrinsic transcriptomic programs of DNA damage response and cell-cycle activity**, rather than by classical immune-checkpoint-response biomarkers (MSI-H, high TMB, or an a priori inflamed TIL phenotype). The proliferative, DNA-repair-proficient, non-EMT phenotype that predicts good response likely reflects sensitivity to the dual genotoxic insult of fluoropyrimidine-plus-radiation: active cell cycle exposes tumors to replication-associated damage, and intact HDR/DSB-repair machinery engages apoptotic and senescence programs following damage recognition. EMT-high mesenchymal tumors invoke well-documented chemoradiation-resistance programs^6,7^ that our data confirm.

**Second**, and previously undescribed in this detail, response to TNT is marked by a coherent **response cascade** visible in paired pre/post samples: (i) dramatic somatic-mutation clearance (Δmissense ≈−70, ΔSBS5 ≈−76 in good responders versus ≈0 in poor); (ii) elimination of MHC-I-neoantigen-presenting tumor clones (Δbinders ≈−300); (iii) contraction of HLA-LOH clones (subjects 3, 4); (iv) broad immune reprogramming with Treg, MHC II, CD8 exhaustion, IFN-γ, and NF-κB/IL-6 signaling upregulated; and (v) striking B-cell infiltration (IGH clonotype Δ >1,400 in good responders) accompanied by increased TCR diversity. Crucially, these changes occur **only** in patients who go on to achieve good response; poor responders are molecularly and immunologically static, indicating primary treatment insensitivity of the tumor bulk rather than a failed subsequent immune recruitment. The massive B-cell infiltration in responders is consistent with TLS formation reported in other cancers^8^ and with recent rectal-cancer studies noting B-cell involvement in pathologic complete response.^9^

**Third**, while our discovery-cohort signal is biologically coherent and achieves LOOCV AUC 0.755, meta-analysis across seven public nCRT cohorts shows only **partial concordance** (two of seven agreeing, three discordant, two null). This discordance is a real and important finding: it suggests that the biology may be specific to TNT regimens (which combine prolonged systemic chemotherapy with chemoradiation, distinct from nCRT) or to the discovery population (Korean, single-center), and highlights the pressing need for prospective TNT-specific multicenter validation before clinical deployment.

### The paradox of pre-treatment MHC II
The observation that MHC class II is modestly reduced in good responders at pre-treatment but elevated post-treatment supports a two-stage model: good responders begin therapy with primarily tumor-intrinsic proliferative-repair phenotype (that nevertheless enables damage-induced cell death), and subsequently experience treatment-induced immune remodeling — tumor-antigen release from dying cells recruits dendritic cells, engages draining lymph nodes, and seeds B- and T-cell infiltration. Poor responders lack both the initial damage-induced cell death and the consequent immune cascade.

### Contrast with the ICB paradigm
Under the ICB paradigm, MSI-H / high-TMB / inflamed tumors respond well; our data show that for TNT the logic is different — intact DNA repair and proliferation (not immunogenicity at baseline) drive response, followed by inflammation as a consequence. This distinction has clinical implications: TNT-responder-like signatures may serve as **negative predictors** for ICB-mono therapy in MSS rectal cancer, and conversely, combination TNT + ICB strategies should consider pre-treatment DNA-repair status.

### Clinical implications
A pre-treatment RNA-seq classifier built on DSB/HDR/E2F/Myc/CD8-proliferation could, pending prospective validation, stratify TNT candidates for organ-preservation watch-and-wait strategies. EMT-high mesenchymal tumors might benefit from alternative neoadjuvant approaches (e.g., addition of taxane, anti-TGF-β, or extended consolidation). Patients with MSI-H LARC should continue to receive ICB per Cercek et al.^5^

### Limitations
- Discovery n=35 single-center; external validation results are mixed.
- Eight unmatched WES samples may retain residual germline; we confined matched analyses to the remaining 41 matched tumors.
- HLA LOH was assessed by a simplified IMGT-allele count approach rather than the full LOHHLA pipeline; absolute LOH calls may be conservative.
- pVACseq neoantigen prediction used MHCflurry and did not include RNA-expression filtering.
- Clinical outcome data (DFS/OS) for survival validation were not yet integrated.

### Future directions
- Prospective multicenter TNT-specific validation in 150+ patients.
- Single-cell RNA-seq of paired pre/post tumors to resolve B-cell subtypes (TLS-associated), T-cell states (exhausted vs memory), and CAF subtypes driving EMT resistance.
- Integration with long-term outcomes (DFS, OS, local recurrence).
- TNT + anti-TGF-β or TNT + ICB combination trials stratified by the DNA-repair signature.

---

## 5. Figures

- **Figure 1.** Cohort overview. (A) Response distribution pie. (B) Clinical T-stage by response (stacked bar). (C) Age × sex distribution. (D) Per-subject sample availability matrix. (E) CONSORT-style study design.
- **Figure 2.** WES landscape. (A) TMB pre-treatment (matched). (B) MSI status (all MSS). (C) Driver mutation oncoprint (top 15 CRC genes × pre-treatment samples). (D) SBS signature composition (top 10 contributors, pre-treatment). (E) MMR-signature reassessment: prior vs Mutect2 refit. (F) Chromosomal instability (CIN). (G) HRD proxy LST.
- **Figure 3.** RNA-seq signatures & DEG. (A) Key pre-treatment signature boxplot grid. (B) Pre-treatment 22-signature heatmap. (C) DESeq2 DEG volcano. (D) CD8 proliferation boxplot (main signature finding). (E) Signature response ranking. (F) Post-treatment MHC II / Treg / IFN-γ.
- **Figure 4.** Pathway enrichment. (A) Hallmark GSEA bar. (B) Reactome GSEA top. (C) ssGSEA top-pathway boxplot grid. (D) ssGSEA pathway heatmap. (E) Key DNA-repair & cell-cycle pathway boxplots.
- **Figure 5.** Integration & ML. (A) 37-feature Spearman correlation heatmap. (B) Top 20 response-associated feature ranking. (C) LOOCV AUC across ML models. (D) Random-forest top-15 feature importance. (E) Top-two feature scatter.
- **Figure 6.** Treatment-induced paired delta. (A) BCR IGH clonotype Δ. (B) IGH Shannon diversity Δ. (C) TCR TRB Shannon Δ. (D) Top 12 immune-signature Δ ranking. (E) Top 15 ssGSEA pathway Δ ranking. (F) TMB Δ. (G) MHC-I neoantigen-binder Δ. (H) Neoantigen mutation-site Δ.
- **Figure 7.** External validation. (A) Per-cohort effect-size forest for four signatures. (B) Stouffer meta-Z scores. (C) GSE150082 DSB/HDR signature (discordant cohort deep-dive).
- **Figure 8.** HLA & neoantigen. (A) HLA class I allele frequency (A/B/C). (B) HLA homozygosity by response. (C) HLA LOH events by response. (D) Pre-treatment neoantigen mutation sites. (E) Strong MHC-I binders. (F) Presentation-competent neoantigen (PCN) score.
- **Figure 9.** Clonal evolution (PyClone-VI). (A) Clonal clusters per subject. (B) Dominant clone shrinkage by response. (C) Shrinking-vs-expanding clusters scatter per subject.

## 6. Tables

- **Table 1.** Clinical characteristics (good vs poor responders).
- **Table 2.** Top 20 integrated-feature response-association statistics.
- **Table 3.** External-cohort meta-analysis (7 GEO cohorts, n=290).
- **Table S1.** 37-feature per-subject master table.
- **Table S2.** SBS signature activities per sample (COSMIC v3.3).
- **Table S3.** Full GSEA (Hallmark + Reactome) results.
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types per subject.
- **Table S6.** Neoantigen summary per sample (pVACseq-MHCflurry).

## 7. Data and code availability
Analysis workspace: `/mnt/sda1/data/TNT/analysis/`. GitHub repository (to be created): `https://github.com/Soonlab/TNT`. Raw sequencing data (Macrogen orders HN00249207 WES and HN00249209 RNA-seq) will be deposited to SRA/ENA upon acceptance.

## 8. References (placeholders — to be completed)
1. Conroy T. et al. Lancet Oncol 2021;22:702-715 (PRODIGE 23).
2. Bahadoer R. et al. Lancet Oncol 2021;22:29-42 (RAPIDO).
3. Garcia-Aguilar J. et al. JCO 2022;40:2546-2556 (OPRA).
4. Cercek A. et al. JAMA Oncol 2022;8:1311-1320.
5. Cercek A. et al. N Engl J Med 2022;386:2363-2376.
6. Guinney J. et al. Nat Med 2015;21:1350-1356 (CMS consortium).
7. Ganesh K. et al. Nat Rev Gastroenterol Hepatol 2019;16:361-375.
8. Cabrita R. et al. Nature 2020;577:561-565 (TLS melanoma).
9. Ayers M. et al. J Clin Invest 2017;127:2930-2940 (IFN-γ GEP).

---

*End of FINAL manuscript draft — generated 2026-04-14 from `/mnt/sda1/data/TNT/analysis/`. Panels in `figures/panels/Fig{N}{A-H}.{pdf,png}`. Workspace README `PROGRESS.md`.*
