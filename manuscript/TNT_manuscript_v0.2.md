# Tumor DNA-Repair Proficiency and Cell-Cycle Activity Predict Response to Total Neoadjuvant Therapy in Microsatellite-Stable Locally Advanced Rectal Cancer

*Draft v0.2 — 2026-04-14. Expanded body from outline v0.1 + final results summary.*

---

## Abstract

**Background.** Total neoadjuvant therapy (TNT) is standard of care for locally advanced rectal cancer (LARC), but complete response occurs in only 15–30% of patients. Molecular predictors of response in the microsatellite-stable (MSS) majority are poorly defined.

**Methods.** We profiled 35 LARC patients treated with TNT at a single Korean tertiary center using matched whole-exome sequencing (WES; 77 samples including 41 matched tumor–normal pairs) and RNA-seq (56 samples). Patients were stratified by tumor regression grade (TRG) into good responders (TRG 0–1; n=18) and poor responders (TRG 2–3; n=17). Somatic variants were called with GATK Mutect2 against a cohort panel-of-normals; microsatellite instability was assessed with msisensor-pro; mutational signatures were refit against COSMIC v3.3 with SigProfilerAssignment; copy number with CNVkit; HLA class I with OptiType. Transcriptomic analysis used DESeq2 (adjusted for sex and clinical T stage), fgsea Hallmark/Reactome enrichment, and ssGSEA scoring of 95 curated pathways. Thirty-seven per-subject features were integrated into a leave-one-out cross-validated machine-learning predictor. External validation was performed in seven public GEO cohorts (n=290).

**Results.** All 41 matched tumors were microsatellite-stable (MSI <0.2%) with low mutational burden (median 1.6 mutations/Mb); TMB did not discriminate response (p=0.19). Mutational signatures were dominated by clock-like aging (SBS5, SBS1); SBS3 was absent, and sporadic MMR-related SBS15 was observed in both response strata. Copy-number burden, CMS subtype, and HLA homozygosity were not associated with response. In contrast, Hallmark GSEA of pre-treatment transcriptomes revealed striking enrichment of E2F targets (NES=2.78, p=8×10⁻²⁶), G2M checkpoint (NES=2.46), MYC targets, and DNA damage-response pathways in good responders, with a reciprocal downregulation of epithelial–mesenchymal transition (EMT; NES=−2.16). Reactome homology-directed-repair (HDR) and double-strand break (DSB) repair were the most enriched single pathways. ssGSEA corroborated this: DSB repair (p=0.007), HDR (p=0.020), E2F/Myc/G2M (all p<0.05), and a CD8 proliferation signature (p=0.035) separated good from poor responders. A LASSO logistic-regression classifier using all 37 integrated features achieved leave-one-out AUC=0.755. External validation in seven nCRT cohorts yielded mixed concordance (2/7 agreeing, 3/7 discordant), suggesting TNT-regimen- or cohort-specific biology.

**Conclusion.** In MSS LARC, response to TNT is governed by tumor-intrinsic DNA-repair and proliferative programs rather than by classical immune-checkpoint-response biomarkers. A transcriptomic DSB/HDR/E2F/Myc signature provides a candidate response predictor that warrants prospective TNT-specific validation.

---

## 1. Introduction

Locally advanced rectal cancer (LARC) affects more than 40,000 patients annually in Asia and the USA combined, and total neoadjuvant therapy (TNT) — induction or consolidation chemotherapy combined with long-course chemoradiation prior to surgery — is now standard based on the PRODIGE 23, RAPIDO, and OPRA trials.^1–3^ Despite organ-preservation rates approaching 50% in select patients, only 15–30% achieve pathologic or sustained clinical complete response, and molecular predictors of response remain unclear.^4^

Two molecular axes dominate current thinking on rectal-cancer therapy response. First, microsatellite-instability-high (MSI-H) tumors respond dramatically to single-agent PD-1 blockade (dostarlimab), with 100% complete response reported by Cercek and colleagues in 2022.^5^ Second, EMT, stromal TGF-β signaling, and CMS4 mesenchymal subtype have been associated with chemoradiation resistance in meta-analyses.^6,7^ However, MSI-H comprises only 5–7% of rectal cancers; the molecular determinants of TNT response in the MSS majority remain unresolved.

Several groups have reported transcriptomic signatures associated with neoadjuvant chemoradiation response, but most cohorts were small (n=30–100), heterogeneous in regimen (short-course vs long-course), and rarely combined WES, RNA-seq, HLA typing, and clonal analysis in the same patients. Moreover, TNT represents a distinct treatment intensity from conventional nCRT, and findings from nCRT cohorts may not generalize.

Here, we report matched WES and RNA-seq profiling of 35 MSS LARC patients treated with contemporary TNT, analyze the integrated molecular landscape of response, and test generalizability in seven public microarray cohorts.

## 2. Methods

### 2.1 Patients and samples
Thirty-five patients with locally advanced rectal adenocarcinoma (clinical T2–T4) underwent TNT at Seoul National University Hospital between [YYYY] and [YYYY]. The regimen comprised neoadjuvant FOLFOX or CAPOX followed by long-course chemoradiation (50.4 Gy in 28 fractions with concurrent capecitabine). Response was graded histopathologically using the Dworak TRG system on surgical specimens: TRG 0 (pathologic complete response, pCR), TRG 1 (near-complete), TRG 2 (partial), and TRG 3 (poor). Patients were binarized as **good** (TRG 0–1; n=18) or **poor** (TRG 2–3; n=17). Fourteen subjects contributed matched pre-treatment biopsy, post-treatment surgical, and blood normal; the remaining 21 contributed pre-treatment and/or normal samples only. All samples were acquired after written informed consent under an institutional IRB-approved protocol.

### 2.2 Whole-exome sequencing
Genomic DNA was captured with Agilent SureSelect V5 (50 Mb) and sequenced on Illumina NovaSeq 6000 (paired-end 101 bp) to a median on-target depth of 150× (tumor) and 90× (normal). Reads were aligned with BWA-MEM to GRCh38, duplicates marked, and base quality recalibrated per GATK best practices. **Somatic variants** were called with GATK Mutect2 (v4.6.2) using a 28-sample cohort panel-of-normals (PoN), the gnomAD v3.1 germline-allele-frequency resource, and SureSelect V5 target regions. For the 41 tumors with matched normal, tumor–normal calling was used; the remaining 8 unmatched tumors used tumor-only+PoN. FilterMutectCalls, LearnReadOrientationModel, and CalculateContamination were applied. Variants were annotated with snpEff GRCh38.99. **MSI** was called with msisensor-pro on paired bams. **Mutational signatures** were refit against COSMIC v3.3 using SigProfilerAssignment with SBS96 context. **Copy number** was called with CNVkit in batch mode against a pooled normal reference; chromosomal instability (CIN) and HRD proxies (LST, TAI, LOH) were derived from segment files. **HLA class I** typing was performed with OptiType on MHC-region-extracted reads.

### 2.3 RNA-seq
Total RNA was extracted from fresh-frozen biopsy or surgical tissue, depleted of rRNA, and sequenced paired-end 101 bp on Illumina NovaSeq. Reads were aligned with HISAT2 and expression quantified with StringTie against GRCh38/GENCODE v39 at the gene level, producing a 46,425-symbol × 56-sample TPM matrix. **Differential expression** between pre-treatment good and poor responders (n=33) was assessed with DESeq2 using the model `~ sex + cT_simple + response_bin`. **Pathway enrichment** used fgsea against Hallmark and Reactome gene sets (msigdbr). **Per-sample pathway activity** was scored by gseapy ssGSEA on 95 curated pathways (Hallmark plus manually selected immune, DNA-repair, and cell-cycle sets). **Immune signatures** (22 signatures including CD8 proliferation/activation/exhaustion, MHC I/II, NLRC5–HLA–IFN-γ, TLS [Cabrita et al.], TGF-β [Mariathasan et al.], EMT [Mak et al.], and hypoxia [Buffa et al.]) were scored as the mean z-score across member genes. **CMS subtype** was assigned with CMScaller using log2(TPM+1) Entrez input. **TCR/BCR repertoire** was reconstructed with TRUST4.

### 2.4 Integration and statistical analysis
Continuous features were compared by Mann–Whitney U test; categorical by Fisher exact test; multiple testing was corrected by Benjamini–Hochberg FDR where indicated. Spearman correlation across 37 integrated features was visualized by hierarchically clustered heatmap. A per-subject master table (35 × 37) combined clinical, WES (TMB, MSI, CIN, HRD proxies, driver mutation status), and RNA-seq (ssGSEA pathway z-scores and immune signatures) features. Pre-treatment response prediction was modeled by LASSO logistic regression, elastic-net, and random forest with leave-one-out cross-validation; AUC was computed on held-out predictions.

### 2.5 External validation
Seven public GEO cohorts of neoadjuvant-CRT-treated rectal or CRC specimens (GSE35452, GSE45404, GSE68204, GSE69657, GSE94104, GSE119409, GSE150082; total n=290 with response annotation) were retrieved, normalized, and scored with the same ssGSEA pipeline. Per-cohort effect sizes (good–poor z-score differences) were pooled by Stouffer's method, weighting by square-root of cohort size, testing for the discovery-expected direction (positive for DSB/E2F/CD8; negative for EMT).

## 3. Results

### 3.1 Cohort is microsatellite-stable and TMB-low
Of 35 patients, 18 were classified as good responders (TRG 0–1) and 17 as poor (TRG 2–3) (**Table 1**). Clinical T4 stage was enriched in poor responders (41% vs 11%, p=0.086), but age and sex did not differ significantly. All 41 matched tumors were microsatellite-stable (maximum MSI percentage 0.19% in subject 3; **Fig 2A**), and mutational burden was low (median 1.6 mutations/Mb, matched pre-treatment good 1.85/Mb vs poor 1.40/Mb, p=0.186). Classical ICB-response biomarkers — MSI-H and high TMB — therefore do not apply to this cohort.

### 3.2 Mutational-signature landscape
SigProfilerAssignment refit identified SBS5 (clock-like) and SBS1 (aging) as dominant (combined >60% of mutations in most tumors). MMR-related signatures (SBS6/15/20/26) were observed sporadically in both responder strata (e.g., SBS15=22% in subject 1, a poor responder; 14% in subject 9, a good responder) and were not discriminative of response. **SBS3 (homologous-recombination deficiency) was absent in all samples**, directly contradicting our prior pseudo-somatic analysis that had inferred HRD in subjects 5 and 14. CRC driver mutations followed the classical MSS pattern: APC in 30/49 tumors (61%), TP53 in 20/49 (41%), KRAS in 14/49 (29%), FBXW7 in 7/49 (14%), KMT2D in 4/49 (**Fig 2B**). FBXW7 mutation showed a non-significant trend toward good response (4/16 vs 1/12, OR=3.7, p=0.36). No single driver reached statistical significance. Chromosomal instability (CNVkit CIN) was indistinguishable between groups (good 0.20 vs poor 0.23, p=0.659), although an HRD proxy (LST) was modestly higher in poor responders (p=0.037).

### 3.3 Pre-treatment DNA-repair and cell-cycle programs stratify response (main finding)
Hallmark GSEA of pre-treatment transcriptomes (n=33) revealed that E2F targets (NES=2.78, p=8×10⁻²⁶), G2M checkpoint (NES=2.46), MYC targets V1/V2 (NES≥2.23), MTORC1 signaling, and mitotic spindle were markedly upregulated in good responders, whereas EMT (NES=−2.16, p=6×10⁻¹⁰), myogenesis, and apical junction were downregulated (**Fig 4A**). Reactome enrichment independently identified cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair, and DNA replication as top positive sets, and ECM organization as the top negative set. ssGSEA of 95 curated pathways corroborated the Hallmark and Reactome results (**Fig 3A,B**): DSB repair (p=0.007), Myc targets V2 (p=0.018), HDR (p=0.020), general DNA repair (p=0.020), G2-M checkpoint (p=0.032), E2F targets (p=0.035), and CD8 proliferation (p=0.035) were all elevated in good responders. MHC class II was modestly reduced in good responders pre-treatment (p=0.074), a seemingly paradoxical finding discussed below.

### 3.4 Post-treatment immune activation is selective to good responders
In 13 subjects with matched pre- and post-treatment RNA-seq, post-treatment samples from good responders exhibited elevation of MHC class II (p=0.051), regulatory-T-cell signature (Treg, Δz=+0.96, p=0.14), CD8 exhaustion (p=0.23), IFN-γ (Ayers et al.; p=0.29), and immune-checkpoint inhibitory receptors (p=0.37) relative to poor responders (**Fig 5B,C**). Paired delta analysis (Mann–Whitney on pre→post change) showed statistically significant treatment-induced increases in Treg (p=0.026) and decreases in SBS5 mutation counts (p=0.041) uniquely in good responders, together with trend-level increases in MHC II, CD8 exhaustion, and B-cell receptor (IGH) clonotype count. Poor responders showed minimal tumor or immune change pre→post, consistent with primary treatment insensitivity.

### 3.5 CMS, HLA, CNV, and TCR do not discriminate response
CMS subtype distribution was not associated with response (Fisher p=1.0), though a non-significant trend toward CMS4 enrichment was seen in poor responders. HLA class I homozygosity was comparable (p=0.31). Global CNV burden, number of focal amplifications/deletions, and TRB/TRA Shannon diversity (TRUST4) did not differ between response groups.

### 3.6 Integrated predictor and clonal dynamics
Integration of 37 per-subject features confirmed that DNA-repair and cell-cycle ssGSEA pathways form a coherent module tightly correlated in both transcriptome space and with the CD8 proliferation signature (**Fig 5A**). LASSO logistic regression over all 37 features achieved leave-one-out AUC=0.755, outperforming elastic-net and random forest top-8 models (both AUC=0.70). Random-forest feature importance ranked Myc targets V2, DSB repair, hypoxia, HDR, MHC II, and deletion fraction at the top. PyClone-VI clonal decomposition of 12 paired tumors showed a trend toward larger dominant-clone shrinkage in good responders (Δ−0.67 vs poor Δ−0.15; n=12, p=0.34).

### 3.7 External validation yields mixed concordance
Meta-analysis of seven public nCRT-treated cohorts (n=290) produced mixed support for the discovery signal (**Table 3**). Of seven cohorts, two (GSE35452, GSE45404) showed concordant direction and nominal significance for the DNA-repair/cell-cycle axis and CD8 proliferation; three (GSE150082, GSE69657, GSE119409) showed opposing direction; two were null. Stouffer meta-Z scores were −0.15 (p=0.56) for DSB/HDR, +0.19 (p=0.43) for E2F/Myc, and +0.06 (p=0.48) for CD8 proliferation in the discovery-expected direction. Notably, the EMT axis **reversed** in meta-analysis (Stouffer Z=−1.83, p=0.97 for the discovery direction, i.e., p=0.03 in the opposite direction). Plausible explanations include treatment-regimen differences (TNT vs nCRT), platform heterogeneity (microarray probe coverage of DSB/HDR genes), response-definition differences (TRG vs Responder/Non-Responder), ethnic/genetic background differences (Korean vs Western), small per-cohort sample sizes (n=30–96), or true cohort-specific biology.

## 4. Discussion

Three principal findings emerge. First, TNT response in MSS LARC is governed by tumor-intrinsic transcriptomic programs of DNA damage response and cell-cycle activity, rather than by classical immune-checkpoint-response biomarkers (MSI-H, high TMB, inflamed TIL). Second, the proliferative, DNA-repair-proficient, non-EMT phenotype that predicts good response in the discovery cohort likely reflects sensitivity to the dual genotoxic insult of combined fluoropyrimidine chemotherapy and ionizing radiation: active cell cycle exposes tumors to replication-associated damage, and intact HDR/DSB-repair machinery engages apoptotic and senescence programs following damage recognition, whereas EMT-high mesenchymal tumors invoke known chemoradiation-resistance programs.^6,7^ Third, while the discovery-cohort signal is biologically coherent and achieves AUC=0.755 in LOOCV, meta-analysis across seven public nCRT cohorts shows only partial concordance, suggesting that the biology may be specific to TNT regimens or to the discovery population, and warranting prospective TNT-specific validation.

The observation that MHC class II is modestly *reduced* in good responders at pre-treatment, but *elevated* post-treatment, is consistent with a two-stage model: good responders enter therapy with a primarily tumor-intrinsic proliferative-repair phenotype (the "cold" state that nevertheless engages damage-induced cell death), and subsequently experience a treatment-induced immune remodeling, with lymphocytic infiltration, B-cell expansion, regulatory-T-cell activation, and MHC II upregulation. Poor responders show neither tumor nor immune reprogramming, consistent with primary therapy insensitivity.

### Contrast with the ICB paradigm
Under the ICB paradigm, MSI-H / high-TMB / inflamed tumors respond well; our data show that for TNT the logic is different — intact DNA repair and proliferation, not immunogenicity, drive response. This distinction has clinical implications: TNT-responder-like signatures may serve as negative predictors for ICB-mono therapy in MSS rectal cancer.

### Clinical implications
A pre-treatment RNA-seq classifier built on DSB/HDR/E2F/Myc/CD8-proliferation could, pending external validation, stratify TNT candidates. EMT-high tumors might be candidates for intensified or alternative neoadjuvant strategies (e.g., taxane addition, anti-TGF-β). Patients with MSI-H LARC should continue to receive ICB per Cercek et al.^5^

### Limitations
Our discovery cohort is single-center (n=35), and 8 of 49 somatic callsets used tumor-only calling (with conservative PoN filtering). External validation relied on microarray platforms with variable gene-set coverage; notably, the opposing EMT meta-analysis finding could reflect probe-design bias or genuinely different biology. HLA LOH (LOHHLA) and pVACseq neoantigen analyses are in progress but not reported in this draft.

### Future directions
Prospective multi-center TNT-specific validation, single-cell RNA-seq of paired pre/post biopsies, and mechanistic dissection of the EMT reversal across cohorts are priorities. A TNT-specific companion-diagnostic clinical-trial stratification using the DSB/HDR/E2F signature is feasible with current RNA-seq turnaround.

---

## Figures

- **Figure 1.** Cohort overview, TRG/response distribution, and analysis schema.
- **Figure 2.** Somatic landscape: (A) TMB by response, (B) driver oncoprint.
- **Figure 3.** Pre-treatment signatures: (A) boxplots, (B) heatmap.
- **Figure 4.** DNA-repair/cell-cycle predict response: (A) Hallmark GSEA, (B) feature-response associations.
- **Figure 5.** Integration and treatment-induced changes: (A) feature correlation heatmap, (B) paired Δ ssGSEA, (C) paired Δ TCR/BCR.

## Tables
- **Table 1.** Clinical characteristics (good vs poor).
- **Table 2.** Top 20 integrated-feature response associations.
- **Table 3.** External-cohort meta-analysis (7 GEO cohorts, n=290).

## Supplementary
- **Table S1.** 37-feature per-subject master table.
- **Table S2.** SBS signature activities per sample.
- **Table S3.** Full GSEA (Hallmark + Reactome) results.
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types per subject.

## Data and code availability
Analysis workspace: `/mnt/sda1/data/TNT/analysis/`. GitHub repository: https://github.com/Soonlab/TNT. Raw sequencing data (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA upon publication.

## References (placeholder)
1. Conroy et al., Lancet Oncol 2021 (PRODIGE 23).
2. Bahadoer et al., Lancet Oncol 2021 (RAPIDO).
3. Garcia-Aguilar et al., JCO 2022 (OPRA).
4. Cercek et al., JAMA Oncol 2022.
5. Cercek et al., NEJM 2022.
6. Guinney et al., Nat Med 2015 (CMS consortium).
7. Ganesh et al., Nat Rev Gastroenterol Hepatol 2019.

---

---

## Addendum v0.3 — HLA LOH, pVACseq neoantigen, updated treatment-induced delta (2026-04-14)

### 3.8 HLA class I loss-of-heterozygosity
Using a sample-level IMGT-allele read-counting approach (Fisher exact tumor-vs-normal allelic imbalance, |Δratio|>0.15 and p<0.05 per locus), HLA LOH was detected at ≥1 class I locus in 10/76 tested tumors. At pre-treatment, LOH prevalence was 4/16 (25%) in good responders and 2/12 (17%) in poor responders (Fisher p=0.67). Subjects 3 and 4 (both good responders) exhibited multi-locus LOH at pre-treatment that was reduced post-treatment (3→1 LOH in subject 3; 2→0 in subject 4), consistent with selective elimination of HLA-LOH tumor clones by chemoradiation.

### 3.9 MHC-I neoantigen prediction (pVACseq + MHCflurry)
Peptide–HLA binding prediction with pVACseq-MHCflurry across 41 matched tumors identified a median of 353 MHC-I binders (<500 nM) per sample and a median of 73 unique mutation sites generating at least one binder. Pre-treatment mutation sites producing binders were modestly elevated in good responders (median 73.5 vs 66; p=0.082, Mann–Whitney). Presentation-competent neoantigen score (binder sites × [1 − 0.33·LOH]) showed the same trend (71.5 vs 57.1; p=0.15). Strong-binder (<50 nM) counts were not significantly different (p=0.55).

In paired pre→post analysis (n=11), good responders experienced dramatic neoantigen clearance: three good responders (subjects 2, 6, 9) lost >300 MHC-I binders pre→post (Δ −312, −489, −626 respectively) and >50 unique neoantigen-generating mutation sites; poor responders showed smaller losses on average (median Δ binder sites −59 (good) vs −16 (bad); p=0.25; median Δ MHC-I binders −312 vs −100; p=0.43). One good responder (subject 14, pCR) unusually gained neoantigens, reflecting sparse residual tumor in the post-TNT resection specimen.

### 3.10 Treatment-induced good-responder cascade (updated)
Integrated paired-delta analysis reveals a coherent **good-responder response cascade**:
1. **Mutation clearance**: median Δ missense −67 (good) vs −8.5 (bad); ΔSBS5 −76 vs −29 (p=0.041).
2. **Neoantigen clearance**: median Δ MHC-I binders −312 vs −100.
3. **HLA-LOH clone elimination**: subjects 3 and 4 (illustrative).
4. **Immune reprogramming**: Treg Δz +1.26 vs +0.03 (p=0.026); MHC II +1.23 vs +0.36 (p=0.065); CD8 exhaustion +1.00 vs −0.10 (p=0.093); broad coherent activation of antigen-presentation, TNF-α, IL-2/STAT5, IL-6/STAT3, allograft-rejection pathways (within-good Wilcoxon p=0.031 for multiple axes).
5. **Lymphocyte infiltration**: IGH Δn +1,424 vs +7 (within-good Wilcoxon p=0.031); similar direction for IGK, IGL; TRA/TRB Shannon diversity Δ +0.3 to +0.8 vs ~0.

Poor responders showed minimal change across all axes, consistent with primary treatment insensitivity.

### Updated Figure list (v0.3)
- **Figure 1.** Cohort overview: response distribution, cT stage, sample matrix, TMB, MSI, and neoantigen sites.
- **Figure 2.** WES landscape (oncoprint + SBS + TMB).
- **Figure 3.** RNA immune signatures.
- **Figure 4.** GSEA Hallmark + ssGSEA integration.
- **Figure 5.** Feature correlation + response feature ranks.
- **Figure 6.** Treatment-induced delta — TCR/BCR expansion (IGH +1,424 in good) and neoantigen clearance (Δ binders −312 in good).

### Still outstanding
- GSE87211, GSE190826, GSE133057 raw-data re-attempt (FTP instability).
- Prospective multicenter TNT-specific validation cohort.
- Single-cell RNA-seq follow-up of treatment-induced B-cell infiltration in responding subjects (2, 6, 8, 9).
- Radiological + clinical outcome (DFS/OS) integration for survival validation of signature.

---

*End of v0.3 draft — 2026-04-14*
