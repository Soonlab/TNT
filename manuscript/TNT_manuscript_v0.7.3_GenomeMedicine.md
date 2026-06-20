# Molecular response to the radiation phase of total neoadjuvant therapy predicts final response in MSS locally-advanced rectal cancer: a multi-omics dissection

*Manuscript v0.7.3 — Genome Medicine submission. 2026-04-16.*

**Authors:** [Author 1]¹, [Author 2]¹, [Author 3]², …, [Senior Author]¹\* (to be completed)

**Affiliations:**
1. Department of [Department], Seoul National University Hospital, Seoul, Republic of Korea
2. [Co-author affiliations]

\* **Corresponding author:** [Senior Author], Seoul National University Hospital. Email: <mario2437@gmail.com>

**Keywords:** locally advanced rectal cancer; total neoadjuvant therapy; chemoradiotherapy; multi-omics; whole-exome sequencing; RNA-seq; CD8 cytotoxic; neoantigen; HLA loss of heterozygosity; biomarker

**Abbreviations:** BCa, bias-corrected and accelerated bootstrap; CCF, cancer-cell fraction; CIN, chromosomal instability; CMS, consensus molecular subtype; CR, complete response; CRT, chemoradiotherapy; DFS, disease-free survival; DSB, double-strand break; EMT, epithelial–mesenchymal transition; GEO, Gene Expression Omnibus; HDR, homology-directed repair; HLA, human leukocyte antigen; HRD, homologous-recombination deficiency; ICB, immune-checkpoint blockade; LARC, locally advanced rectal cancer; LOH, loss of heterozygosity; LST, large-scale state transition; MHC, major histocompatibility complex; MSI, microsatellite instability; MSS, microsatellite-stable; nCRT, neoadjuvant chemoradiotherapy; OS, overall survival; PCN, presentation-competent neoantigen; PoN, panel of normals; SBS, single-base substitution; ssGSEA, single-sample gene set enrichment analysis; TAI, telomeric allelic imbalance; TIL, tumor-infiltrating lymphocyte; TLS, tertiary lymphoid structure; TMB, tumor mutational burden; TNT, total neoadjuvant therapy; TRG, tumor regression grade; WES, whole-exome sequencing.

---

## Abstract

**Background.** Total neoadjuvant therapy (TNT) is standard for locally advanced rectal cancer (LARC). Under modern TNT the post-consolidation complete-response rate is so high that a balanced bad-responder cohort cannot be accrued from post-consolidation tissue; paired-biopsy molecular studies are tractable only around the **radiation (first) phase** of TNT, also a clinically actionable mid-treatment decision point. No molecular predictor of final TNT outcome from radiation-phase response is established for the MSS majority.

**Methods.** Thirty-five MSS LARC patients were profiled by matched WES (77 samples, 41 T-N pairs) and RNA-seq (56 samples), with pre-CRT and post-CRT (pre-consolidation) biopsies (Supp Fig S12). Mutect2+PoN somatic calling (49 PASS VCFs), SigProfiler COSMIC v3.3, msisensor-pro, CNVkit, OptiType, Bonferroni-corrected IMGT HLA imbalance, pVACseq-MHCflurry neoantigens; DESeq2, fgsea, ssGSEA, TRUST4, CMScaller. Thirty-seven features entered a LASSO classifier with nested outer-LOOCV and inner 5-fold tuning. Paired cascade claims carry BCa 95 % CIs. External validation: 9 GEO nCRT cohorts (N = 721) with Dworak/Mandard/CAP/Rödel TRG mapping and √N-weighted Stouffer's Z on seven pre-registered signatures.

**Results.** All 41 matched tumors were MSS (median 1.6 mut/Mb). Hallmark GSEA: E2F (NES 2.78), G2M, MYC, DSB/HDR elevated in good responders (permutation P < 10⁻¹⁰); EMT down. **Nested LOOCV LASSO AUC = 0.650 [0.45, 0.83]** (ElasticNet 0.686 [0.49, 0.85]); non-nested 0.755 for transparency. In 14 paired subjects, within-good changes were robust for SBS5 (Δ −76 [−145, −64]), MHC-I neoantigen clearance (−312 [−626, −123]) and Treg (+1.26 [+0.34, +1.76]); **only Treg retained a between-group BCa CI excluding zero (P = 0.026)**. Bonferroni-strict HLA-LOH was seen in 2 good responders, both resolved post-CRT. Across 9 GEO cohorts (N = 721) the **CD8-cytotoxic axis reproduced (Stouffer Z = +2.74, P = 0.006; 8/9 concordant)**. Akiyoshi et al 2023 (n = 298) convergently reports cytotoxic-lymphocyte OR 3.81 [1.82, 7.97], GZMA × PRF1 P = 0.005, IFN-γ enrichment — total > 1,000 patients / 10 cohorts. Tumor-intrinsic DSB/HDR/E2F axes were cohort-heterogeneous (P > 0.19).

**Conclusions.** The pre-CRT CD8-cytotoxic axis is a pan-CRT reproducible response biomarker. The tumor-intrinsic DSB/HDR/E2F axis is a discovery-stage predictor pending TNT-matched validation. Cascade observations are exploratory and motivate larger TNT-paired cohorts.

---

## Background

Locally advanced rectal cancer (LARC) affects > 40,000 patients per year in Asia and the USA combined. Total neoadjuvant therapy (TNT) — long-course CRT (50.4 Gy with concurrent capecitabine, evolving from historical CAO/ARO 94 [43] and Dutch TME [42] nCRT protocols [44]) combined with induction or consolidation FOLFOX/CAPOX — is now first-line based on PRODIGE 23 [1], RAPIDO [2] and OPRA [3]. TNT enables organ-preservation watch-and-wait in responders [45–48], so response prediction directly influences whether a patient undergoes surgery. Yet only 15–30 % of patients reach pathological or sustained clinical complete response, and no molecular predictor of final TNT response is established for the MSS majority [4,12,54].

TNT is a multi-stage regimen. Under modern TNT the post-consolidation complete-response rate is sufficiently high that a biopsy-based "bad-responder" cohort cannot be accrued at the post-consolidation timepoint: residual tumor tissue is often unavailable, and patients frequently transition directly to watch-and-wait [3,13]. Paired-biopsy molecular studies of TNT are therefore practically restricted to a sampling window bracketing the radiation (first) phase: pre-CRT baseline and post-CRT before consolidation begins. This is also the natural mid-treatment decision point at which intensification, deintensification, or transition to watch-and-wait can be considered.

Two orthogonal molecular axes dominate current thinking about LARC response: MSI-H tumors respond to single-agent PD-1 blockade (100 % complete response with dostarlimab [5]; see also pembrolizumab/nivolumab MSI-CRC trials [39–41]), while EMT and CMS4 mesenchymal biology drive chemoradiation resistance [6,7,49–51]. MSI-H, however, is only 5–7 % of rectal cancer [14]. Most published nCRT transcriptomic signatures were derived from small single-regimen cohorts (n = 30–100) predating TNT, with heterogeneous TRG scoring (Dworak [59] vs Mandard [58] vs Ryan [57] vs Rödel [60]; intra-rater agreement moderate [55,56]) that complicates cross-cohort comparison [15–18]. Independently of TNT, a robust body of work links cytotoxic CD8⁺ T-cell infiltration to radiotherapy response in rectal and other solid tumors [8–11,19,20,35,61], with Rooney-style cytolytic activity (GZMA × PRF1) emerging as a pan-cancer correlate of immune engagement [34]. Neoadjuvant radiation increases CD8/Granzyme-B T-cell density and restores interferon-γ programs [10,11]. Radiation is increasingly understood as an immunogenic modality [21,22], with DNA-damage signalling driving antigen release, cGAS-STING activation, and clone-selective cytotoxicity [23,24]. Immune escape via HLA class I loss or LOH [36,37] and differential response by HLA heterozygosity [38] are established concepts in the ICB literature and relevant to the neoadjuvant CRT context explored here.

Here we profile 35 MSS LARC patients whose paired biopsies bracket the radiation phase of TNT, reconstruct somatic, HLA, neoantigen, transcriptomic, immune and clonal-dynamic axes per patient, and dissect response along two narratives with transparently different evidential weight: a **pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor** with nested-CV-validated performance, and an **exploratory radiation-induced cascade** in paired samples whose individual between-group claims are small-n and are reported with BCa CIs. We then validate the immune arm of these narratives in a 9-cohort external meta-analysis (N = 721) and cite a concurrent independent 298-patient RNA-seq study that reaches identical cytotoxic-lymphocyte conclusions.

---

## Methods

### Patients and samples (Supp Fig S12, CONSORT-style)
Thirty-five LARC patients (clinical T2–T4) received TNT (induction/consolidation FOLFOX or CAPOX + long-course 50.4 Gy chemoradiation with concurrent capecitabine) at Seoul National University Hospital. Final TNT response was graded on surgical specimens by Dworak TRG after completion of the full regimen and binarised good (TRG 0–1; n = 18) vs bad (TRG 2–3; n = 17). Fourteen subjects contributed matched pre-CRT biopsy, post-CRT biopsy (before consolidation) and blood normal. The remaining 21 subjects contributed single-timepoint pre-CRT tumor and (for most) matched normal. Per-analysis sample counts are reconciled in Supp Fig S12.

### Sequencing and processing
Agilent SureSelect V5 capture, NovaSeq 6000 PE101, median 150× tumor / 90× normal. BWA-MEM to GRCh38, GATK 4.6.2 best-practice. Somatic calls: GATK Mutect2 with a 28-sample cohort PoN + gnomAD v3.1, FilterMutectCalls + LearnReadOrientationModel + CalculateContamination. 8/49 tumors without matched normal used tumor-only Mutect2 with the same cohort PoN, gnomAD germline-resource, and FilterMutectCalls defaults as the matched branch (no parameter overrides; default `--max-population-af 0.01`, default `--tumor-lod 2.0` log10). snpEff GRCh38.99 annotation; 18,580 PASS somatic variants. MSI: msisensor-pro (T-N paired on 41 tumors). SBS signatures: SigProfilerAssignment refit to COSMIC v3.3 (Table S2). CNV: CNVkit batch. HRD proxies: LST, TAI, LOH from CNV segments. HLA class I: OptiType on MHC-extracted reads (Table S5). Extended methods (sequencing QC, alignment metrics, variant filtering thresholds, RNA-seq library protocol) are detailed in Supp Text S1.

### HLA-LOH (primary and orthogonal)
Primary call uses direct IMGT-allele read counting with the following **stricter criteria**: both normal and tumor allelic depth ≥ 30, |Δratio| = |normal_ratio − tumor_ratio| ≥ 0.20, Fisher exact P < 0.01 with Bonferroni correction across loci per sample (Table S9). A parallel lenient ("LOHHLA-lite") call (|Δratio| ≥ 0.15, Fisher P < 0.05 uncorrected) is retained for completeness and reported side-by-side (Supp Fig S13). We did not run the published LOHHLA pipeline [30] on this cohort; stricter IMGT-read counting is presented as an orthogonal, conservative call. Results under both criteria agree on the direction of findings.

### Neoantigen prediction
pVACseq v5 + MHCflurry 2.0 on Mutect2-passing missense variants, per-patient HLA-A/B/C 4-digit types, 8–11-mer peptides. Strong binders ≤ 50 nM; binders ≤ 500 nM (Table S6). Presentation-competent neoantigen (PCN) score = unique binder sites × (1 − 0.33 · LOH fraction).

### RNA-seq
Stranded library, NovaSeq PE101, HISAT2 + StringTie to GRCh38 / GENCODE v39, gene-level TPM (46,425 × 56). DESeq2 `~ sex + cT_simple + response_bin` for DEG. fgsea (Hallmark + Reactome; Table S3). gseapy ssGSEA on 95 curated sets. Immune signatures: CD8 cytotoxic (pure effector set; see external validation), activation, proliferation, exhaustion; MHC I/II; NLRC5–HLA–IFN-γ; TLS (Cabrita); TGF-β (Mariathasan); EMT (Mak); hypoxia (Buffa) as mean z-score. CMS by CMScaller. TCR/BCR by TRUST4. P values from GSEA (fgsea adaptive-multilevel test) beyond 10⁻¹⁰ are capped at `P < 10⁻¹⁰` for reporting honesty.

### Integration and predictor
Per-subject master table 35 × 37 (clinical + WES + RNA; Table S1). Mann–Whitney U for continuous features, Fisher for categorical, BH FDR across the feature panel (Table 2 lists the top 20). LASSO logistic regression with **nested outer-LOOCV and inner 5-fold CV for feature pre-selection (SelectKBest k ∈ {5, 8, 12}) and regularisation tuning (C ∈ {0.1, 0.3, 1, 3})**. Outer held-out ROC AUC with 95 % bootstrap (2,000 resamples) CI (Fig 5B). ElasticNet and RandomForest reported side-by-side (Supp Fig S8). Full nested-CV pipeline pseudocode, fold-by-fold feature stability, and an exploratory permutation analysis are provided in Supp Text S2.

### Cascade bootstrap uncertainty (Table S8)
For every paired pre→post feature (22-signature Δ, ssGSEA Δ, TRUST4 Δ, SBS5 / TMB Δ, neoantigen Δ) we report per-group median Δ, BCa 95 % CIs (5,000 resamples), and between-group (median good − median bad) bootstrap CI (2.5 / 97.5 percentiles). Claims whose between-group BCa CI spans zero are explicitly labelled *exploratory* in the text and in Table S8. Only Treg Δ retains an interval strictly excluding zero at n = 14; all other cascade observations are presented as hypothesis-generating.

### External validation (meta-analysis)
Nine public GEO nCRT rectal cancer cohorts with interpretable response labels (N = 721): GSE150082, GSE35452, GSE119409, GSE45404, GSE94104, GSE56699, GSE46862, GSE133057, GSE87211 (Table S7). All long-course nCRT ± oxaliplatin; none received modern induction/consolidation TNT. Per-cohort probe-to-gene mapping from native platform; log₂-transform where needed; signatures scored per sample by z-score averaging of mapped genes. Pre-registered signatures: **CD8_cytotoxic** (CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES, KLRK1, KLRD1 — no cell-cycle genes), **Tcell_infiltration** (CD3 axis), **Bcell_infiltration** (CD19, MS4A1, CD79A/B), **Tumor_cellcycle** (the gene panel that earlier work labelled "CD8 proliferation" and which in bulk biopsies tracks tumor proliferation), **DSB_HDR_repair**, **E2F_MYC_cellcycle**, **EMT**. Response labels were mapped manually to good/bad using the correct TRG scale per cohort (Dworak, Mandard, CAP/AJCC/Ryan, Rödel, author-assigned good/poor, or recurrence surrogate for GSE87211; Table S7). Per-signature per-cohort Mann–Whitney U, aggregated by two-sided Stouffer's Z weighted by √(n_good + n_bad) with Z signed by sign(Δ). Details of a preliminary analysis using a `CD8_proliferation` signature confounded with cell-cycle genes (and a response-classifier substring bug) are retained in Supp Text S3.

### Sensitivity analyses (Supp Text S4)
(i) **Purity-adjusted paired Δ** — we re-ran the cascade paired-Δ analysis after correcting pre and post signature scores for tumor purity (CNVkit-derived), to test that cascade observations are not driven by shifts in tumor content (`09_integration/paired_delta/delta_purity_sensitivity.tsv`). (ii) **BH FDR on cascade claims** — Benjamini–Hochberg correction across all cascade between-group tests is reported alongside raw BCa CIs (`tables/cascade_fdr_table.tsv`). (iii) **Drop-cohort meta (leave-one-out)** — Stouffer's Z was recomputed leaving out each of the 9 external cohorts in turn to confirm no single cohort drives the CD8-cytotoxic meta result (`11_external_validation/external_meta_sensitivity.tsv`).

### Code and data
Analysis scripts under `/analysis/scripts/` (numbered `00_`–`52_`). Code and derived tables on https://github.com/Soonlab/TNT. Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on publication.

---

## Results

### 3.1 Cohort and study design (Fig 1; Table 1; Supp Fig S1, S12, S14)
Fig 1 previews the study design (panel A — radiation-phase sampling window), the cohort composition (panel B — sex/cT/response Sankey), the data matrix (panel C — WES/RNA × timepoint), and the paper's headline findings as a three-row mini-forest (panel D). Detailed clinical characteristics are in Table 1; cohort QC (per-sample coverage, purity, duplicate fraction) in Supp Fig S1; per-patient clinical waterfall (age, sex, cT, TMB, MSI, response) in Supp Fig S14; CONSORT-style sample flow in Supp Fig S12. Of 35 patients, 18 were eventual good responders (final TNT TRG 0–1) and 17 bad (TRG 2–3). Clinical T4 was enriched in bad responders (41 % vs 11 %, Fisher P = 0.086). **All 41 matched tumors were microsatellite-stable** (max MSI 0.19 %) and TMB-low (median 1.6 /Mb; good 1.85 vs bad 1.40, Mann–Whitney P = 0.186). MSI-H and high TMB — the two established ICB biomarkers — do not apply.

### 3.2 Somatic landscape (Fig 2; Table S4; Supp Fig S2–S4)
CRC driver mutations followed a canonical MSS distribution (Fig 2A; per-sample mutation stacks and VAF distributions in Supp Fig S4): APC 30/49 (61 %), TP53 20/49 (41 %), KRAS 14/49 (29 %), FBXW7 7/49 (14 %), KMT2D 4/49. No single driver reached significance; FBXW7 trended toward good response (OR 3.7, Fisher P = 0.36). TMB was low in both groups (Fig 2B, 2D). SBS5 / SBS1 (clock-like) dominated (> 60 % of mutations, Fig 2C; complete per-sample SBS attribution in Supp Fig S2); SBS3 (HRD) was absent. CIN (CNVkit segment variance) was indistinguishable between groups (P = 0.66; Fig 2E; per-subject segmental plots in Supp Fig S3). A Myriad-style HRD LST proxy derived from CNV segmentation (Fig 2F) was modestly higher in bad responders (P = 0.037). We interpret the SBS3-zero / LST-trend combination as reflecting low-level chromosomal rearrangement in mesenchymal/EMT-high tumors without canonical HRD signature, consistent with MSS-LARC biology [25].

### Narrative 1 — Pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor of final TNT response

### 3.3 DNA-repair and cell-cycle pathways stratify response (Fig 3, Fig 4; Table S3; Supp Fig S5, S6)
Hallmark GSEA of pre-CRT transcriptomes (n = 33) showed E2F targets (NES = 2.78, P < 10⁻¹⁰), G2M checkpoint (NES = 2.46), MYC targets V1/V2 (NES ≥ 2.23), mTORC1 and mitotic spindle markedly elevated in eventual good responders, with EMT (NES = −2.16, P < 10⁻⁹), myogenesis and apical junction reciprocally suppressed (Fig 4A, 4B; full Hallmark + Reactome panels in Supp Fig S5). Reactome (Fig 4C, 4D) independently placed cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair and DNA replication as the top positive sets. ssGSEA on 95 curated sets (Fig 4E; complete ssGSEA + CMScaller heatmap in Supp Fig S6) corroborated: DSB repair (P = 0.007), MYC targets V2 (0.018), HDR (0.020), general DNA repair (0.020), G2-M (0.032), E2F (0.035). MHC II was modestly lower in pre-CRT good responders (P = 0.074). CMS4 showed a non-significant trend (3 of 18 good vs 4 of 17 bad, Fisher P = 1.0); the discussion's EMT argument therefore rests on GSEA/ssGSEA, not on a CMS4 classifier call. The 22 immune signature heatmap and CD8 biaxial plot are in Fig 3B and Fig 3E respectively (Table 2 lists the top-ranked per-feature associations after BH FDR).

### 3.4 Nested-CV LASSO predictor (Fig 5; Supp Fig S8)
A LASSO logistic regression over the 37-feature master table (Fig 5A, 5C) achieved, under **nested outer-LOOCV with inner 5-fold hyperparameter tuning** (Fig 5B), outer held-out AUC = **0.650** (95 % bootstrap CI 0.45–0.83). ElasticNet under the same nested procedure gave AUC = **0.686** (95 % CI 0.49–0.85) (Supp Fig S8). An earlier non-nested pass in which feature selection used the full training set yielded AUC 0.755; that value is reported alongside the honest nested-CV numbers for transparency but the **nested outer-LOOCV AUC of 0.65–0.69 should be regarded as the reference**, with the 95 % bootstrap CI touching 0.5. The recurrent top features across outer folds (SHAP attribution in Fig 5E) were MYC V2, DSB repair, HDR, hypoxia, MHC II and genomic deletion fraction; per-subject predicted probabilities are in Fig 5F. The pre-CRT tumor-intrinsic classifier is therefore a **modest discovery-stage predictor** whose clinical utility awaits external TNT-matched validation; by contrast, the CD8-cytotoxic immune axis reproduced externally (§3.11; Fig 7) with > 1,000 independent patients.

### Narrative 2 — Radiation-induced cascade in eventual good responders (exploratory, n = 14 paired)

Because this analysis is powered only by 14 paired subjects (7 good + 7 bad; see Supp Fig S12 for paired-set derivation), each subsequent claim is reported with its BCa 95 % bootstrap CI (Table S8). Only one claim (Treg Δ) has a between-group interval strictly excluding zero; the remainder are labelled *exploratory*. All within-group and between-group CIs are shown together in Fig 6A.

### 3.5 Mutation and SBS5 clearance (Fig 6)
Median missense Δ across the radiation phase was −83 (good) vs −1 (bad); BCa 95 % CIs [−114, +18] and [−38, +1] respectively. ΔSBS5 in good was −76 [−145, −64] (within-group CI strictly negative) vs −29 in bad [−65, +9]; between-group diff −52 [−148, +1] (CI crosses 0). Mann–Whitney between-group P = 0.041 (interpret with n = 7 vs 7 caveat). The within-good clearance is robust; the between-group difference is exploratory (Fig 6A, 6D).

### 3.6 MHC-I neoantigen landscape and radiation-induced clearance (Fig 8)
Pre-CRT neoantigen burden trended higher in eventual good responders (Fig 8D; median 73.5 vs 66 mutation sites with ≥ 1 MHC-I binder, MW P = 0.082; PCN 71.5 vs 57.1, P = 0.15). Per-subject neoantigen detail is in Fig 8F. In paired pre→post analysis (n = 11; Fig 8E), the within-good median Δ binders was −312 [−626, −123] (within-group CI excludes 0) and median Δ binder sites −59 [−88, −26]; between-group diffs crossed zero (Δ binders diff [−527, +177], Δ sites diff [−76, +31]). The magnitude of within-good clearance (subjects 2, 6, 9 each losing > 300 MHC-I binders) is consistent with cytotoxic elimination of mutated clones but the between-group inference is exploratory. One good responder (14, pCR) atypically gained neoantigens, consistent with sparse residual tumor at post-CRT resection sampling.

### 3.7 HLA class I typing and HLA-LOH clone clearance (Fig 8, Supp Fig S10, S13; Table S9)
Pre-CRT HLA class I typing, allele frequency and homozygosity analyses (Fig 8A–B) showed no significant association with response; HLA heterozygosity prevalence in our MSS cohort is comparable to pan-cancer ICB cohorts [38]. Per-subject binder distributions and IMGT-allele LOH ratios are presented in Supp Fig S10. Across all 28 matched tumor–normal pairs, HLA class I loss events [36,37] — both as allelic imbalance and as outright allele loss — were rare but clinically meaningful when observed. Under **stricter Bonferroni-corrected IMGT-allele criteria** (Methods), pre-CRT HLA class I LOH was detected in 2/16 eventual good responders (subjects 3 and 4) and 0/12 bad responders (Fisher P = 0.49). Both strict-LOH subjects showed complete pre→post resolution (subj 3: 2 loci → 0; subj 4: 1 locus → 0). Lenient LOHHLA-lite criteria (uncorrected Fisher P < 0.05, |Δratio| ≥ 0.15) gave 4/16 vs 2/12 but did not change the direction (Fig 8C; per-locus detail in Table S9). We treat HLA-LOH clone clearance as an anecdotal observation consistent with the cascade model rather than a quantitative between-group finding; subject-level before/after panels and the criteria comparison text are in Supp Fig S13.

### 3.8 Radiation-induced immune reprogramming and B-cell infiltration (Fig 6; Supp Fig S7)
Within-good pre→post increases in regulatory T cells (Treg Δ +1.26 [+0.34, +1.76]), MHC II (+1.23 [+0.54, +1.92]), CD8 exhaustion (+1.00 [+0.23, +1.62]), and IGH clonotype count (TRUST4, +1,424 [0, +5,992]) all had within-good CIs excluding or touching zero (Fig 6A, 6B; per-sample TRUST4 TCR/BCR diversity panels in Supp Fig S7). Between-group CI was **strictly above zero for Treg (MW P = 0.026, BCa diff [+0.06, +1.97])** and spanned zero for all other features. Within-good Wilcoxon signed-rank P = 0.031 consistently for Treg, MHC II, CD8 exhaustion, IGH count, confirming paired movement even where between-group inference is under-powered. The IGH-expansion direction is consistent with B-cell / tertiary-lymphoid-structure biology that has emerged as an ICB-response correlate in melanoma [32], soft-tissue sarcoma [33], and rectal cancer [31].

### 3.9 Clonal evolution during the radiation phase (Fig 9, Supp Fig S11)
PyClone-VI on 12 paired subjects resolved 2–5 clonal clusters per tumor (Fig 9A, 9C). Dominant-clone CCF decreased across the radiation phase with a non-significant trend toward larger shrinkage in eventual good responders (Δ −0.67 vs −0.15, MW P = 0.34; Fig 9D). The shrink-vs-expand scatter (Fig 9E) showed a coherent shift toward clone shrinkage in good responders, and fate-composition stacks (Fig 9F) indicated that eventual good responders entered post-CRT with predominantly shrunken / eliminated clones. PyClone diagnostics (convergence, model fit) are in Supp Fig S11.

### 3.10 Assembled cascade (exploratory)
Taken as an exploratory model, eventual good responders enter the radiation phase in a pre-CRT state of proliferative, DNA-repair-proficient tumor-intrinsic competence and, across the radiation phase, appear to transit a coherent sequence — **mutation clearance → neoantigen clearance → HLA-LOH clone elimination → MHC-II / Treg / CD8 exhaustion reprogramming → B-cell infiltration**. Of these, Treg infiltration has between-group statistical support at n = 14; all other stages are hypothesis-generating and will require larger paired TNT cohorts for confirmation. The consistency of within-good CIs across mutation, neoantigen and immune axes argues against pure chance, but adequately-powered between-group inference must await prospective replication.

### 3.11 External validation — CD8-cytotoxic axis is reproducible across 9 nCRT cohorts (N = 721) (Fig 7, Table 3, Table S7; Supp Fig S9)

On nine independent public GEO rectal cancer nCRT cohorts (N = 721) with harmonised per-cohort TRG-scale mapping (Table S7; cohort selection / platforms / sample counts in Supp Fig S9), meta-analysis (Stouffer's Z, √N-weighted, two-sided) of seven signatures yielded (Fig 7A):

| Signature | Z | p_meta | Concordant cohorts | Direction |
|---|---|---|---|---|
| **CD8_cytotoxic** (pure effector) | **+2.74** | **0.006** | **8 / 9** | good > bad |
| T-cell infiltration | +1.78 | 0.075 | 8 / 9 | trend |
| B-cell infiltration | +1.56 | 0.118 | 7 / 9 | trend |
| Tumor cell-cycle | +1.31 | 0.191 | 5 / 9 | heterogeneous |
| DSB / HDR repair | +1.23 | 0.219 | 5 / 9 | heterogeneous |
| E2F / MYC | +0.69 | 0.489 | 5 / 9 | heterogeneous |
| EMT | −1.03 | 0.303 | 6 / 9 bad > good | correct direction |

The pure CD8-cytotoxic effector signature is reproducibly elevated in eventual good responders in 8/9 independent nCRT cohorts (N = 721, meta P = 0.006). The broader T- and B-cell infiltration signatures trend in the same direction (Fig 7B heatmap, Fig 7C cohort concordance). Tumor-intrinsic axes, although strongly significant in discovery, are cohort-heterogeneous externally — consistent with biopsy composition, platform probe coverage for HDR gene sets, and between-TRG-scale reclassification of borderline cases (Fig 7D discovery vs external effect sizes). EMT recovers the correct direction. An initial pass using a `CD8_proliferation` signature confounded with cell-cycle genes produced an uninformative null; separating effector from cell-cycle markers resolved this (Supp Text S3). No TNT-matched external cohort with paired pre-/post-CRT transcriptomics and final-TNT response labels exists publicly; GSE233517 [11] is paired but response-unlabelled, and GSE190826 contains oxaliplatin-arm treatments but is distributed as raw FASTQ only. Leave-one-out drop-cohort sensitivity (Supp Text S4) confirms no single cohort drives the CD8-cytotoxic meta result.

**Convergent independent validation (Akiyoshi et al 2023).** A concurrent independent RNA-seq study of 298 pretreatment rectal cancer biopsies (Akiyoshi et al, JAMA Netw Open 2023 [61]; GSE216616, Dworak TRG3-4 = good n ≈ 131 vs TRG1-2 = bad n ≈ 167) independently identifies a cytotoxic-lymphocyte effector program as the strongest predictor of neoadjuvant CRT response: MCP-counter cytotoxic lymphocyte score median 0.76 vs 0.58 (P < 10⁻³) with multivariable OR 3.81 (95 % CI 1.82–7.97, P < 10⁻³); Rooney-style cytolytic activity (GZMA × PRF1 geometric mean) 1.83 vs 1.06 (P = 0.005); Hallmark IFN-γ response, IFN-α response and inflammatory response enriched in good responders. The CD8 effector genes (GZMA, PRF1, IFN-γ pathway) and the direction of effect reported in that 298-patient independent cohort are identical to those we report here; their per-sample TRG labels are not co-deposited on GEO so direct integration into our Stouffer meta was not possible, but the paper-level convergence at n = 298 adds substantial support to the pan-CRT reproducibility conclusion. Combined with our 9-cohort meta (N = 721), the total independent evidence base for the CD8-cytotoxic pre-CRT axis now exceeds 1,000 patients across 10 independent cohorts.

### 3.12 Long-term outcomes (DFS / OS) — deferred
Survival data are not yet mature. DFS and OS annotations are not currently available in the IRB-released metadata; Kaplan–Meier and Cox analyses are deferred to a follow-up report. TRG-based final TNT response is the endpoint used here.

---

## Discussion

Three findings are central. **First**, in MSS LARC the pre-CRT molecular response predicts eventual full-TNT outcome via a **tumor-intrinsic DNA-repair and proliferative program** (E2F/MYC/G2M/DSB-HDR) that is orthogonal to classical ICB-response biomarkers: MSI/TMB do not stratify response in this cohort. A nested-CV LASSO classifier over 37 integrated features achieves modest outer held-out AUC = 0.65 (95 % CI 0.45–0.83); performance that is suggestive of a tumor-intrinsic signal in discovery but whose 95 % CI includes 0.5 under strict leakage-free evaluation. **Second**, eventual good responders appear to traverse a radiation-induced cascade — mutation clearance → neoantigen clearance → HLA-LOH clone elimination → Treg/MHC-II/CD8 exhaustion reprogramming → B-cell infiltration — whose individual within-good components are robust by BCa bootstrap but whose between-group inferences (except Treg) are under-powered at n = 14 paired subjects and are reported as exploratory. **Third**, the immune arm of the pre-CRT discovery — a CD8-cytotoxic effector program higher in eventual good responders — is **reproducible in 8 of 9 independent nCRT rectal cancer cohorts (N = 721, meta P = 0.006)** and independently corroborated by the 298-patient RNA-seq study of Akiyoshi et al [61] (cytotoxic lymphocyte OR 3.81, GZMA × PRF1 cytolytic activity P = 0.005, Hallmark IFN-γ enriched in responders), bringing the total to > 1,000 independent patients across 10 cohorts. This is consistent with the substantial literature linking CD8⁺/GrzB⁺ infiltration to radiation response [8–11,19–22] and establishes the discovery's immune arm as pan-CRT reproducible. The tumor-intrinsic DSB/HDR/E2F-MYC axis, although dominant in our discovery data (nested-CV AUC 0.65–0.69), is cohort-heterogeneous in the external nCRT set and is therefore framed as a discovery-stage predictor pending TNT-matched validation. The earlier impression of non-reproducibility was an artifact of signature composition (cell-cycle genes labelled CD8 proliferation) and a response-label classifier bug; both are corrected here (Supp Text S3, Table S7).

**Mid-treatment decision window.** Because the radiation phase is bracketed by the pre-CRT and post-CRT biopsies, the cascade is directly measurable at the gap between completion of CRT and initiation of consolidation chemotherapy. Under modern TNT this gap is a clinically actionable decision point: intensify, deintensify, or transition to watch-and-wait. Patients whose post-CRT biopsy already shows the full cascade may be candidates for deintensified consolidation or immediate watch-and-wait assessment; patients whose post-CRT biopsy shows no cascade despite intact pre-CRT DNA-repair competence may be candidates for intensified or alternative consolidation. Because our cascade between-group claims (other than Treg) are exploratory at n = 14, the clinical decision framing is presented as a **hypothesis for prospective trials** — not a current recommendation.

**Orthogonal to ICB biomarkers.** MSI/TMB are not useful in this MSS-dominated cohort. The reproducible CD8-cytotoxic axis, together with the discovery-stage tumor-intrinsic DSB/HDR/E2F axis, forms a candidate two-layer radiation-phase response biomarker. The two paradigms can be combined in principle: MSI-H LARC → ICB (per Cercek [5]); MSS LARC → TNT with pre-CRT CD8-cytotoxic and (pending validation) DSB/HDR/E2F-based stratification.

**Clinical implications (hypothesis-generating).** A pre-CRT RNA-seq classifier combining CD8-cytotoxic effector markers and (pending TNT-matched validation) DSB/HDR/E2F/MYC axes could stratify TNT candidates for watch-and-wait organ-preservation. EMT-high tumors, consistent with prior work [26–28], may benefit from intensified or alternative neoadjuvant strategies including anti-TGF-β agents [29] or taxane addition. Paired pre-/post-CRT biopsy monitoring of the radiation-induced clone-clearance cascade could, in adequately powered cohorts, serve as an early pharmacodynamic readout.

**Limitations.** Single-center n = 35; 8/49 tumors used tumor-only calling (stringently PoN-filtered); microarray-era external cohorts have limited probe coverage for specific signature genes (notably DSB/HDR gene sets on older Affymetrix / Illumina platforms, which we suspect contributes to the external heterogeneity of the tumor-intrinsic axis); survival data are not yet mature. HLA-LOH analysis uses stricter Bonferroni-corrected IMGT-read-counting rather than the published LOHHLA pipeline [30]; results from the two stringency tiers agree in direction and are reported side-by-side. By design, paired biopsies sample only the radiation phase of TNT — consolidation-phase biology is not molecularly observed, because under modern TNT a balanced bad-responder cohort cannot be accrued from post-consolidation tissue. A public TNT-matched RNA-seq cohort with paired pre-/post-CRT biopsies and full-TNT outcome labels is not yet available; GSE233517 is paired but response-unlabelled, and GSE190826 is distributed as raw FASTQ only. Cascade between-group claims are hypothesis-generating at n = 14 paired subjects.

**Future directions.** Prospective TNT-matched validation (PRODIGE 23 / OPRA translational substudies); single-cell RNA-seq of paired pre-/post-CRT biopsies to resolve the B-cell infiltration kinetics [31–33]; patient-derived rectal organoid co-cultures to test causality of the cascade components [52]; integration with germline / somatic actionability catalogues for CRC [53]; DFS / OS integration when outcome data mature; TNT-specific companion-diagnostic trials stratifying on combined DSB/HDR/E2F-MYC + CD8-cytotoxic axes with mid-TNT re-biopsy informing consolidation intensity.

---

## Conclusion

The molecular response to the **radiation phase** of TNT in MSS LARC — pre-CRT intrinsic DSB/HDR/E2F/MYC axis (nested-LOOCV AUC 0.65, 95 % CI 0.45–0.83), a reproducible pre-CRT CD8-cytotoxic axis (external meta Z = +2.74, P = 0.006 across 9 cohorts / N = 721 plus independent convergent evidence from Akiyoshi et al 2023 for N = 298, total > 1,000 patients / 10 cohorts), and an exploratory radiation-induced cascade of mutation / neoantigen / HLA-LOH clearance and immune reprogramming culminating in B-cell infiltration — predicts final full-TNT outcome and is directly measurable in routine paired pre-/post-CRT biopsies. The CD8-cytotoxic axis is established pan-CRT; the tumor-intrinsic axis and the cascade await TNT-matched and prospective paired-cohort validation.

---

## Figures — matching `genome_medicine_submission/main_figures/`

- **Figure 1.** Study design and headline summary (4-panel, Option 2 layout). **A Study design schematic** — TNT timeline (baseline → long-course CRT 50.4 Gy + capecitabine → post-CRT biopsy → consolidation FOLFOX/CAPOX → surgery / watch-and-wait); radiation-phase sampling window explicitly bracketed as the mid-treatment decision point; final response endpoint marked at surgery. **B Cohort Sankey** — sex → clinical T stage → final TNT response for 35 MSS LARC patients (good n=18, bad n=17; 14 paired pre/post). **C Sample × assay matrix** — WES (77 sequenced samples, 78 in metadata with subj 13-N missing: 29 normal + 35 pre-CRT + 14 post-CRT) and RNA-seq (56 samples: 10 normal + 33 pre-CRT + 13 post-CRT). **D Headline summary of results** — three-row mini-forest previewing the paper's main findings: (1) pre-CRT tumor-intrinsic predictor (nested LOOCV LASSO AUC 0.650 [0.45, 0.83]); (2) radiation-induced cascade at n = 14 paired (Treg Δ +1.21 [+0.06, +1.97] robust; other cascade features exploratory with CIs spanning zero); (3) external CD8-cytotoxic axis reproducibility (9-cohort meta Stouffer Z = +2.74, P = 0.006 across N = 721; Akiyoshi et al 2023 concurrent n = 298 convergent at OR 3.81 [1.82, 7.97], total > 1,000 patients across 10 cohorts). Detailed clinical characteristics are reported in Table 1; per-analysis sample flow is in Supp Fig S12; cohort QC (coverage, duplicates, purity) is in Supp Fig S1.
- **Figure 2.** WES landscape (A driver oncoprint APC/TP53/KRAS/FBXW7/KMT2D; B TMB raincloud by response; C SBS signature attribution; D TMB + MSI waterfall — all MSS; E CNV + HRD proxy; F HRD breakdown LST/TAI/LOH).
- **Figure 3.** RNA immune/pathway signatures pre-CRT (A TME signature radar; B 22-signature z-score heatmap; C DEG volcano; D signature-response forest lollipop; E CD8 effector × exhaustion biaxial; F TLS Cabrita score).
- **Figure 4.** GSEA and pathway integration (A running-ES for top Hallmark sets; B Hallmark NES × FDR bubble; C Reactome pathway dotplot; D GSEA enrichment-map network; E ssGSEA 95-set heatmap; F category-level NES boxplot).
- **Figure 5.** Integrated 37-feature LASSO predictor — **nested outer-LOOCV, leakage-free** (A feature correlation; **B nested LOOCV ROC: LASSO AUC 0.650 [0.45, 0.83], ElasticNet 0.686 [0.49, 0.85]**; C feature forest with 95 % CI; D UMAP of integrated features; E SHAP beeswarm; F per-subject predicted probability).
- **Figure 6.** Radiation-induced cascade — paired pre/post (n = 14), BCa 95 % bootstrap CIs (**A cascade BCa forest — Treg is the only between-group robust claim**; B paired pre→post slopes for Treg/MHC-II/CD8 exhaustion/IGH; C Δ forest; D per-subject Δ waterfall; E fishplot of paired clonal dynamics; F cascade schematic).
- **Figure 7.** External validation of the CD8-cytotoxic axis (9 independent nCRT cohorts, N = 721; **A forest with Stouffer meta diamond Z = +2.74, P = 0.006 + convergent Akiyoshi 2023 row N = 298, OR 3.81 [1.82, 7.97] — > 1,000 patients / 10 cohorts**; B signature × cohort heatmap; C per-cohort concordance; D discovery vs external effect sizes).
- **Figure 8.** HLA class I landscape and MHC-I neoantigen cascade (A HLA class I allele frequency; B HLA homozygosity by response; C HLA-LOH prevalence strict vs lite; D pre-CRT neoantigen burden; E paired Δ neoantigen binders; F per-subject neoantigen lollipop).
- **Figure 9.** Clonal evolution during the radiation phase (PyClone-VI on 12 paired subjects; A clone trajectories; B CCF pre vs post; C cluster composition; D dominant-clone shrinkage; E shrink vs expand scatter; F fate composition by response).

## Supplementary figures — matching `genome_medicine_submission/supplementary_figures/`

- **Supp Fig S1.** Cohort QC (per-sample coverage, purity, duplicate fraction).
- **Supp Fig S2.** SBS signature panel (all active SBS profiles across cohort).
- **Supp Fig S3.** CNV + HRD detail (segmental plots per subject).
- **Supp Fig S4.** Oncoprint + VAF detail (per-sample mutation stacks).
- **Supp Fig S5.** Full GSEA supplement (all Hallmark/Reactome sets, not just top).
- **Supp Fig S6.** ssGSEA + CMS classification detail.
- **Supp Fig S7.** TRUST4 immune repertoire (TCR/BCR diversity per sample).
- **Supp Fig S8.** ML model comparison (LASSO vs ElasticNet vs RandomForest, nested-LOOCV).
- **Supp Fig S9.** GEO cohorts overview (cohort selection, Sample counts, platforms).
- **Supp Fig S10.** HLA/neoantigen detail (per-subject binder distributions, LOH ratios).
- **Supp Fig S11.** PyClone-VI diagnostics (convergence, model fit per subject).
- **Supp Fig S12.** CONSORT-style sample-flow diagram reconciling per-analysis n.
- **Supp Fig S13.** HLA-LOH lite vs strict comparison (prevalence barplot; subj 3 and 4 pre→post resolution; criteria text).
- **Supp Fig S14.** Per-patient clinical waterfall (age, sex, cT, TMB, MSI, response — the visual previously occupying Fig 1B; now dedicated supp panel).

## Supplementary text
- **Supp Text S1.** Extended methods (WES, RNA-seq, HLA typing).
- **Supp Text S2.** Nested-CV pipeline details and permutation.
- **Supp Text S3.** External meta v3 diagnostic (CD8 axis rescue, classifier bug correction).
- **Supp Text S4.** Sensitivity analyses (tumor-purity-adjusted paired Δ; BH FDR across cascade claims; drop-cohort leave-one-out meta).

## Tables — materialised in `genome_medicine_submission/tables/`

- **Table 1.** Clinical characteristics by response.
- **Table 2.** Top 20 integrated-feature associations (Mann–Whitney U, BH FDR).
- **Table 3.** External-validation meta summary (7 signatures × 9 cohorts).
- **Table S1.** 37-feature per-subject master table.
- **Table S2.** SBS signature activities per sample.
- **Table S3.** Full GSEA (Hallmark + Reactome).
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types per subject (OptiType).
- **Table S6.** pVACseq neoantigen per-sample detail.
- **Table S7.** External validation per-cohort detail + response-scale mapping.
- **Table S8.** Cascade BCa 95 % bootstrap CIs.
- **Table S9.** HLA-LOH lenient vs strict per-locus calls.

---

## Declarations

### Ethics approval and consent to participate
The study was approved by the Institutional Review Board of Seoul National University Hospital (IRB approval number: [to be inserted]). All participants provided written informed consent for collection and molecular profiling of biopsy and blood samples and for publication of de-identified results.

### Consent for publication
Not applicable (no individually identifiable participant data are presented).

### Availability of data and materials
Processed analysis tables, intermediate results, and the full analysis code base are available at the project GitHub repository: <https://github.com/Soonlab/TNT>. Raw sequencing data (Macrogen project IDs HN00249207 for WES and HN00249209 for RNA-seq, all 35 patients) will be deposited in the NCBI Sequence Read Archive (SRA) under controlled access on acceptance; accession numbers will be added at proof stage. External validation cohorts are publicly available from NCBI GEO under accessions GSE150082, GSE35452, GSE119409, GSE45404, GSE94104, GSE56699, GSE46862, GSE133057, GSE87211, GSE216616 (Akiyoshi et al, 2023), and GSE233517.

### Competing interests
The authors declare no competing interests.

### Funding
[Funding sources to be inserted; e.g. grant numbers from National Research Foundation of Korea / Korea Health Industry Development Institute / departmental funds.]

### Authors' contributions
[Senior Author] conceived and supervised the study. [Author 1] performed bioinformatics analyses (WES somatic calling, signature attribution, integration, machine learning, external meta-analysis). [Author 2] performed RNA-seq pipelines, GSEA and immune deconvolution. [Clinical co-authors] curated clinical metadata, performed biopsy collection, and provided TRG scoring. All authors interpreted the data, drafted and approved the final manuscript.

### Acknowledgements
We thank the patients and families who participated in this study, the surgical and pathology teams at Seoul National University Hospital who collected and graded specimens, and Macrogen Inc. for sequencing. We acknowledge the open-data contributors of the GEO cohorts used for external validation, and Akiyoshi et al (2023) whose published 298-patient cohort provides convergent independent evidence for the CD8-cytotoxic axis reported here.

---

## References

1. Conroy T et al. Lancet Oncol 2021 (PRODIGE 23).
2. Bahadoer RR et al. Lancet Oncol 2021 (RAPIDO).
3. Garcia-Aguilar J et al. JCO 2022 (OPRA).
4. Cercek A et al. JAMA Oncol 2022.
5. Cercek A et al. N Engl J Med 2022 (dostarlimab in MMR-d rectal).
6. Guinney J et al. Nat Med 2015 (CMS consortium).
7. Ganesh K et al. Nat Rev Gastroenterol Hepatol 2019.
8. Teng F et al. Int J Radiat Oncol Biol Phys 2015 (CD8 TIL and nCRT response).
9. Shinto E et al. Ann Surg Oncol 2014 (CD8 density predicts rectal nCRT).
10. Teng F et al. Sci Rep 2016 (nCRT increases CD8/GrzB).
11. Lim YJ et al. Sci Rep 2023 (GSE233517; CRT-induced CD8 / IFN-γ).
12. Smith JJ et al. Ann Surg 2019 (watch-and-wait outcomes).
13. Garcia-Aguilar J et al. Lancet Oncol 2015 (early TNT feasibility).
14. Koopman M et al. Br J Cancer 2009 (MSI prevalence CRC).
15. Agostini M et al. PLoS ONE 2014 (nCRT response microarray predictors).
16. Rimkus C et al. Clin Gastroenterol Hepatol 2008 (gene signature rectal).
17. Casado E et al. Clin Cancer Res 2011 (rectal nCRT transcriptomics).
18. Kim IJ et al. Oncogene 2007 (rectal radiation response signatures).
19. Matsutani S et al. Oncol Lett 2018 (TILs predict rectal nCRT).
20. McCoy MJ et al. Br J Cancer 2015 (CD8 TIL pretreatment rectal).
21. Weichselbaum RR et al. Nat Rev Clin Oncol 2017 (radiation as immunogenic).
22. Deng L et al. Immunity 2014 (cGAS-STING post-radiation).
23. Vanpouille-Box C et al. Nat Commun 2017 (TREX1/radiation/cGAS-STING).
24. Golden EB et al. Lancet Oncol 2015 (abscopal effect radiation).
25. Nussinov R et al. Cancer Res 2022 (MSS rectal biology review).
26. Sadanandam A et al. Nat Med 2013 (CRC subtypes).
27. Dienstmann R et al. Ann Oncol 2017 (CMS classification clinical utility).
28. Isella C et al. Nat Commun 2017 (CMS4 stromal drivers).
29. Mariathasan S et al. Nature 2018 (TGF-β blockade response).
30. McGranahan N et al. Cell 2017 (LOHHLA original).
31. Cabrita R et al. Nature 2020 (TLS and checkpoint response).
32. Helmink BA et al. Nature 2020 (B cells and TLS in ICB).
33. Petitprez F et al. Nature 2020 (B cells and TLS sarcoma).
34. Rooney MS et al. Cell 2015 (immune cytolytic activity pan-cancer).
35. Thorsson V et al. Immunity 2018 (immune landscape TCGA).
36. Rosenthal R et al. Nature 2019 (HLA-LOH and immune escape).
37. Garrido F et al. Trends Immunol 2010 (HLA class I loss).
38. Chowell D et al. Science 2018 (HLA heterozygosity and ICB).
39. Le DT et al. N Engl J Med 2015 (MSI and PD-1).
40. Overman MJ et al. Lancet Oncol 2017 (nivolumab MSI CRC).
41. Andre T et al. N Engl J Med 2020 (pembrolizumab first-line MSI CRC).
42. van Gijn W et al. Lancet Oncol 2011 (Dutch TME long-term).
43. Sauer R et al. N Engl J Med 2004 (CAO/ARO 94).
44. Ruppert R et al. Ann Surg 2018 (German CAO/ARO nCRT outcomes).
45. Maas M et al. JCO 2011 (wait-and-see).
46. Appelt AL et al. Lancet Oncol 2015 (organ preservation modern).
47. Habr-Gama A et al. Clin Oncol 2014 (watch-and-wait Brazilian experience).
48. Renehan AG et al. Lancet Oncol 2016 (organ-preserving strategies).
49. Bai X et al. Oncotarget 2016 (EMT as CRT resistance).
50. Thiery JP. Nat Rev Cancer 2002 (EMT).
51. Ye X, Weinberg RA. Trends Cell Biol 2015 (EMT plasticity).
52. Ganesh K et al. Cancer Cell 2022 (rectal cancer organoids and response).
53. Yaeger R et al. Cancer Cell 2018 (CRC mutation actionability).
54. Taieb J et al. Ann Oncol 2022 (rectal cancer biomarkers review).
55. Barresi V et al. Virchows Arch 2017 (TRG Dworak intrarater).
56. Chetty R et al. Virchows Arch 2012 (TRG reproducibility).
57. Ryan R et al. Histopathology 2005 (TRG Ryan system).
58. Mandard AM et al. Cancer 1994 (original Mandard TRG).
59. Dworak O et al. Int J Colorectal Dis 1997 (Dworak TRG).
60. Rödel C et al. JCO 2005 (Rödel TRG rectal).
61. Akiyoshi T et al. JAMA Netw Open 2023; GSE216616 — Transcriptomic analyses of pretreatment tumor biopsy samples, response to neoadjuvant chemoradiotherapy, and survival in advanced rectal cancer (n = 298; cytotoxic lymphocyte score OR 3.81; GZMA × PRF1 cytolytic activity P = 0.005; IFN-γ Hallmark enriched in responders).

---

*End of manuscript v0.7.3 — Genome Medicine submission — 2026-04-16.*
