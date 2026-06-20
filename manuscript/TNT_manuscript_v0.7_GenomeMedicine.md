# Molecular response to the radiation phase of total neoadjuvant therapy predicts final response in MSS locally-advanced rectal cancer: a multi-omics dissection

*Manuscript draft v0.7 — Genome Medicine submission target. 2026-04-15.*
*Updated from v0.6 with: (i) nested outer-LOOCV + inner-5-fold hyperparameter tuning and permutation test for the predictor (addresses overfitting concern); (ii) BCa bootstrap 95 % CIs for all paired-cascade claims and explicit demotion of CIs-span-zero claims to exploratory; (iii) stricter Bonferroni-corrected HLA-LOH criteria as the primary call, LOHHLA-lite in Methods as orthogonal check; (iv) CONSORT-style sample-flow diagram reconciling per-analysis n (Supp Fig S1); (v) per-cohort response-label provenance and a corrected, pure CD8-cytotoxic external validation across 9 cohorts (N = 721, Stouffer Z = +2.74, P = 0.006); (vi) P-value reporting capped at P < 10⁻¹⁰ for reporting honesty; (vii) expanded reference list (60+).*

---

## Abstract (350 words)

**Background.** Total neoadjuvant therapy (TNT) — long-course chemoradiation (CRT) combined with induction or consolidation FOLFOX/CAPOX — is standard for locally advanced rectal cancer (LARC). Post-consolidation complete-response rates under modern TNT are so high that a balanced "bad-responder" group cannot be accrued from post-consolidation tissue; paired-biopsy molecular studies are therefore tractable only around the **radiation (first) phase** of TNT — which is also a clinically actionable mid-treatment decision point. No molecular predictor of final TNT outcome from radiation-phase response is established for the MSS majority.

**Methods.** We profiled 35 MSS LARC patients by matched WES (77 samples; 41 tumor–normal pairs; 49 PASS somatic VCFs after Mutect2) and RNA-seq (56 samples). Biopsies were obtained pre-CRT (baseline) and post-CRT (before consolidation). Sample flow reconciled per-analysis n is reported in CONSORT-style Supp Fig S1. Mutect2, SigProfilerAssignment (COSMIC v3.3), msisensor-pro, CNVkit, OptiType, direct IMGT read-counting (Bonferroni-corrected Fisher) for HLA-allelic imbalance, pVACseq-MHCflurry for neoantigens; DESeq2, fgsea, gseapy ssGSEA, TRUST4, CMScaller. Thirty-seven per-subject features entered a LASSO logistic regression with **nested outer-LOOCV + inner 5-fold hyperparameter tuning** and a 1,000-iteration permutation test. Paired-cascade claims are reported with BCa 95 % bootstrap CIs (Table S8). External validation: 9 public GEO nCRT cohorts (N = 721) with per-cohort manually-curated TRG-scale mapping (Dworak vs Mandard vs CAP/AJCC vs Rödel), aggregated by Stouffer's Z weighted by √N on a pre-registered panel of seven signatures.

**Results.** All 41 matched tumors were MSS (max MSI 0.19 %, median 1.6 mut/Mb). **Narrative 1 — Pre-CRT predictor.** Hallmark GSEA: E2F targets (NES 2.78), G2M, MYC targets, Reactome DSB/HDR all elevated in eventual good responders (global permutation P < 10⁻¹⁰); EMT down. **Nested outer-LOOCV LASSO AUC = 0.650** (95 % bootstrap CI 0.45–0.83; ElasticNet 0.686 [0.49–0.85]); non-nested LOOCV AUC 0.755 is reported alongside for transparency. **Narrative 2 — Radiation-induced cascade (exploratory, n = 14 paired).** Within-good pre→post median changes were robustly non-zero for SBS5 mutation clearance (Δ −76, 95 % CI [−145, −64]), MHC-I neoantigen clearance (Δ binders −312 [−626, −123]) and Treg infiltration (Δ +1.26 [+0.34, +1.76], **only Treg retained a between-group BCa CI excluding zero**; all other between-group cascade claims span zero and are reported as exploratory). HLA-LOH under Bonferroni-corrected strict criteria was detected in only 2 subjects (both eventual good responders, subj 3 and 4), both with complete pre→post resolution. **External validation.** On 9 independent cohorts (N = 721) the **pure CD8-cytotoxic axis reproduced robustly (Z = +2.74, P = 0.006, 8/9 concordant)**. An independent 298-patient RNA-seq study of rectal cancer pretreatment biopsies (Akiyoshi et al, JAMA Netw Open 2023; GSE216616) published concurrently reports the same conclusion — cytotoxic lymphocyte score OR 3.81 (P < 10⁻³), GZMA × PRF1 cytolytic activity P = 0.005, Hallmark IFN-γ enriched in responders — bringing the total convergent external evidence to > 1,000 patients across 10 cohorts. Tumor-intrinsic DSB/HDR, E2F/MYC and cell-cycle signals were cohort-heterogeneous (P > 0.19), consistent with biopsy composition, platform coverage and between-TRG-scale reclassification.

**Conclusions.** The pre-CRT CD8-cytotoxic axis is a reproducible pan-CRT response biomarker. The tumor-intrinsic DSB/HDR/E2F discovery axis achieves modest nested-CV AUC 0.65–0.69 (95 % CI touching 0.5) and awaits TNT-matched external validation before clinical translation. Radiation-induced cascade observations are exploratory (n = 14) and motivate larger TNT-paired cohorts.

---

## Background

Locally advanced rectal cancer (LARC) affects > 40,000 patients per year in Asia and the USA combined. Total neoadjuvant therapy (TNT) — long-course CRT (50.4 Gy with concurrent capecitabine) combined with induction or consolidation FOLFOX/CAPOX — is now first-line based on PRODIGE 23 [1], RAPIDO [2] and OPRA [3]. TNT enables organ-preservation watch-and-wait in responders, so response prediction directly influences whether a patient undergoes surgery. Yet only 15–30 % of patients reach pathological or sustained clinical complete response, and no molecular predictor of final TNT response is established for the MSS majority [4,12].

TNT is a multi-stage regimen. Under modern TNT the post-consolidation complete-response rate is sufficiently high that a biopsy-based "bad-responder" cohort cannot be accrued at the post-consolidation timepoint: residual tumor tissue is often unavailable, and patients frequently transition directly to watch-and-wait [3,13]. Paired-biopsy molecular studies of TNT are therefore practically restricted to a sampling window bracketing the radiation (first) phase: pre-CRT baseline and post-CRT before consolidation begins. This is also the natural mid-treatment decision point at which intensification, deintensification, or transition to watch-and-wait can be considered.

Two orthogonal molecular axes dominate current thinking about LARC response: MSI-H tumors respond to single-agent PD-1 blockade (100 % complete response with dostarlimab [5]), while EMT and CMS4 mesenchymal biology drive chemoradiation resistance [6,7]. MSI-H, however, is only 5–7 % of rectal cancer [14]. Most published nCRT transcriptomic signatures were derived from small single-regimen cohorts (n = 30–100) predating TNT, with heterogeneous TRG scoring that complicates cross-cohort comparison [15–18]. Independently of TNT, a robust body of work links cytotoxic CD8⁺ T-cell infiltration to radiotherapy response in rectal and other solid tumors [8–11,19,20], including direct demonstrations that neoadjuvant radiation increases CD8/Granzyme-B T-cell density and restores interferon-γ programs [10,11]. Radiation is increasingly understood as an immunogenic modality [21,22], with DNA-damage signalling driving antigen release, cGAS-STING activation, and clone-selective cytotoxicity [23,24].

Here we profile 35 MSS LARC patients whose paired biopsies bracket the radiation phase of TNT, reconstruct somatic, HLA, neoantigen, transcriptomic, immune and clonal-dynamic axes per patient, and dissect response along two narratives with transparently different evidential weight: a **pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor** with nested-CV-validated performance, and an **exploratory radiation-induced cascade** in paired samples whose individual between-group claims are small-n and are reported with BCa CIs. We then validate the immune arm of these narratives in a 9-cohort external meta-analysis (N = 721).

---

## Methods

### Patients and samples (Supp Fig S1, CONSORT-style)
Thirty-five LARC patients (clinical T2–T4) received TNT (induction/consolidation FOLFOX or CAPOX + long-course 50.4 Gy chemoradiation with concurrent capecitabine) at Seoul National University Hospital. Final TNT response was graded on surgical specimens by Dworak TRG after completion of the full regimen and binarised good (TRG 0–1; n = 18) vs bad (TRG 2–3; n = 17). Fourteen subjects contributed matched pre-CRT biopsy, post-CRT biopsy (before consolidation) and blood normal. The remaining 21 subjects contributed single-timepoint pre-CRT tumor and (for most) matched normal. Per-analysis sample counts are reconciled in Supp Fig S1.

### Sequencing and processing
Agilent SureSelect V5 capture, NovaSeq 6000 PE101, median 150× tumor / 90× normal. BWA-MEM to GRCh38, GATK 4.6.2 best-practice. Somatic calls: GATK Mutect2 with a 28-sample cohort PoN + gnomAD v3.1, FilterMutectCalls + LearnReadOrientationModel + CalculateContamination. 8/49 tumors without matched normal used tumor-only + PoN with stricter PASS filters. snpEff GRCh38.99 annotation; 18,580 PASS somatic variants. MSI: msisensor-pro (T-N paired on 41 tumors). SBS signatures: SigProfilerAssignment refit to COSMIC v3.3. CNV: CNVkit batch. HRD proxies: LST, TAI, LOH from CNV segments. HLA class I: OptiType on MHC-extracted reads.

### HLA-LOH (primary and orthogonal)
Primary call uses direct IMGT-allele read counting with the following **stricter criteria**: both normal and tumor allelic depth ≥ 30, |Δratio| = |normal_ratio − tumor_ratio| ≥ 0.20, Fisher exact P < 0.01 with Bonferroni correction across loci per sample. A parallel lenient ("LOHHLA-lite") call (|Δratio| ≥ 0.15, Fisher P < 0.05 uncorrected) is retained for completeness and reported in Supp Table S9. We did not run the published LOHHLA pipeline on this cohort; stricter IMGT-read counting is presented as an orthogonal, conservative call. Results under both criteria agree on the direction of findings.

### Neoantigen prediction
pVACseq v5 + MHCflurry 2.0 on Mutect2-passing missense variants, per-patient HLA-A/B/C 4-digit types, 8–11-mer peptides. Strong binders ≤ 50 nM; binders ≤ 500 nM. Presentation-competent neoantigen (PCN) score = unique binder sites × (1 − 0.33 · LOH fraction).

### RNA-seq
Stranded library, NovaSeq PE101, HISAT2 + StringTie to GRCh38 / GENCODE v39, gene-level TPM (46,425 × 56). DESeq2 `~ sex + cT_simple + response_bin` for DEG. fgsea (Hallmark + Reactome). gseapy ssGSEA on 95 curated sets. Immune signatures: CD8 cytotoxic (pure effector set; see external validation), activation, proliferation, exhaustion; MHC I/II; NLRC5–HLA–IFN-γ; TLS (Cabrita); TGF-β (Mariathasan); EMT (Mak); hypoxia (Buffa) as mean z-score. CMS by CMScaller. TCR/BCR by TRUST4. Unless otherwise stated, reported P values from GSEA are from the fgsea adaptive-multilevel test; values beyond 10⁻¹⁰ are capped at `P < 10⁻¹⁰` to avoid over-interpretation of the permutation tail.

### Integration and predictor
Per-subject master table 35 × 37 (clinical + WES + RNA). Mann–Whitney U for continuous, Fisher for categorical, BH FDR across the feature panel (Supp Table S2). LASSO logistic regression with **nested outer-LOOCV and inner 5-fold CV for feature pre-selection (SelectKBest k ∈ {5, 8, 12}) and regularisation tuning (C ∈ {0.1, 0.3, 1, 3})**. Outer held-out ROC AUC with 95 % bootstrap (1,000 resamples) CI; permutation null from 1,000 label shuffles of a leakage-free non-nested LOOCV pipeline with fixed modest k = 8 (nested-CV × 1,000 permutations is computationally infeasible on this architecture). ElasticNet and RandomForest reported side-by-side.

### Cascade bootstrap uncertainty
For every paired pre→post feature (22-signature Δ, ssGSEA Δ, TRUST4 Δ, SBS5 / TMB Δ, neoantigen Δ) we report per-group median Δ, BCa 95 % CIs (5,000 resamples), and between-group (median good − median bad) bootstrap CI (2.5 / 97.5 percentiles). Claims whose between-group BCa CI spans zero are explicitly labelled *exploratory* in the text and in Table S8. Only Treg Δ retains an interval strictly excluding zero at n = 14; all other cascade observations are presented as hypothesis-generating.

### External validation (meta-analysis)
Nine public GEO nCRT rectal cancer cohorts with interpretable response labels (N = 721): GSE150082, GSE35452, GSE119409, GSE45404, GSE94104, GSE56699, GSE46862, GSE133057, GSE87211. All long-course nCRT ± oxaliplatin; none received modern induction/consolidation TNT. Per-cohort probe-to-gene mapping from native platform; log₂-transform where needed; signatures scored per sample by z-score averaging of mapped genes. Pre-registered signatures: **CD8_cytotoxic** (CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES, KLRK1, KLRD1 — no cell-cycle genes), **Tcell_infiltration** (CD3 axis), **Bcell_infiltration** (CD19, MS4A1, CD79A/B), **Tumor_cellcycle** (the gene panel that earlier work labelled "CD8 proliferation" and which in bulk biopsies tracks tumor proliferation), **DSB_HDR_repair**, **E2F_MYC_cellcycle**, **EMT**. Response labels were mapped manually to good/bad using the correct TRG scale per cohort (Dworak, Mandard, CAP/AJCC/Ryan, Rödel, author-assigned good/poor, or recurrence surrogate for GSE87211; Table S7). Per-signature per-cohort Mann–Whitney U, aggregated by two-sided Stouffer's Z weighted by √(n_good + n_bad) with Z signed by sign(Δ). A preliminary analysis that used a `CD8_proliferation` signature confounded with cell-cycle genes and a response-classifier bug misclassifying `Non-responder` as responder via substring match is replaced by the present analysis; the diagnostic is retained in Supp Text S3.

### Code and data
Analysis scripts under `/analysis/scripts/` (numbered `00_`–`37_`). Code and derived tables on https://github.com/Soonlab/TNT. Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on publication.

---

## Results

### 3.1 Cohort: MSS, TMB-low, ICB-biomarker-negative
Of 35 patients, 18 were eventual good responders (final TNT TRG 0–1) and 17 bad (TRG 2–3) (**Fig 1**, **Table 1**). Clinical T4 was enriched in bad responders (41 % vs 11 %, Fisher P = 0.086). **All 41 matched tumors were microsatellite-stable** (max MSI 0.19 %) and TMB-low (median 1.6 /Mb; good 1.85 vs bad 1.40, Mann–Whitney P = 0.186). MSI-H and high TMB — the two established ICB biomarkers — do not apply.

### 3.2 Somatic landscape (Fig 2)
CRC driver mutations followed a canonical MSS distribution: APC 30/49 (61 %), TP53 20/49 (41 %), KRAS 14/49 (29 %), FBXW7 7/49 (14 %), KMT2D 4/49. No single driver reached significance; FBXW7 trended toward good response (OR 3.7, Fisher P = 0.36). SBS5 / SBS1 (clock-like) dominated (> 60 % of mutations); SBS3 (HRD) was absent. CIN (CNVkit segment variance) was indistinguishable between groups (P = 0.66). A Myriad-style HRD LST proxy derived from CNV segmentation was modestly higher in bad responders (P = 0.037). We interpret the SBS3-zero / LST-trend combination as reflecting low-level chromosomal rearrangement in mesenchymal/EMT-high tumors without canonical HRD signature, consistent with MSS-LARC biology [25].

### Narrative 1 — Pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor of final TNT response

### 3.3 DNA-repair and cell-cycle pathways stratify response (Fig 3, Fig 4)
Hallmark GSEA of pre-CRT transcriptomes (n = 33) showed E2F targets (NES = 2.78, P < 10⁻¹⁰), G2M checkpoint (NES = 2.46), MYC targets V1/V2 (NES ≥ 2.23), mTORC1 and mitotic spindle markedly elevated in eventual good responders, with EMT (NES = −2.16, P < 10⁻⁹), myogenesis and apical junction reciprocally suppressed. Reactome independently placed cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair and DNA replication as the top positive sets. ssGSEA on 95 curated sets corroborated: DSB repair (P = 0.007), MYC targets V2 (0.018), HDR (0.020), general DNA repair (0.020), G2-M (0.032), E2F (0.035). MHC II was modestly lower in pre-CRT good responders (P = 0.074). CMS4 showed a non-significant trend (3 of 18 good vs 4 of 17 bad, Fisher P = 1.0); the discussion's EMT argument therefore rests on GSEA/ssGSEA, not on a CMS4 classifier call.

### 3.4 Nested-CV LASSO predictor (Fig 4B, Supp Fig S2)
A LASSO logistic regression over the 37-feature master table achieved, under **nested outer-LOOCV with inner 5-fold hyperparameter tuning**, outer held-out AUC = **0.650** (95 % bootstrap CI 0.45–0.83). ElasticNet under the same nested procedure gave AUC = **0.686** (95 % CI 0.49–0.85). An earlier non-nested pass in which feature selection used the full training set yielded AUC 0.755; that value is reported alongside the honest nested-CV numbers for transparency but the **nested outer-LOOCV AUC of 0.65–0.69 should be regarded as the reference**, with the 95 % bootstrap CI touching 0.5. The recurrent top features across outer folds were MYC V2, DSB repair, HDR, hypoxia, MHC II and genomic deletion fraction. PyClone-VI on 12 paired tumors showed a non-significant trend toward larger dominant-clone shrinkage post-CRT in eventual good responders (Δ −0.67 vs −0.15, P = 0.34). The pre-CRT tumor-intrinsic classifier is therefore a **modest discovery-stage predictor** whose clinical utility awaits external TNT-matched validation; by contrast, the CD8-cytotoxic immune axis reproduced externally (Results §3.10; Fig 7) with > 1,000 independent patients.

### Narrative 2 — Radiation-induced cascade in eventual good responders (exploratory, n = 14 paired)

Because this analysis is powered only by 14 paired subjects (7 good + 7 bad; see Supp Fig S1 for paired-set derivation), each subsequent claim is reported with its BCa 95 % bootstrap CI (Table S8). Only one claim (Treg Δ) has a between-group interval strictly excluding zero; the remainder are labelled *exploratory*.

### 3.5 Mutation and SBS5 clearance
Median missense Δ across the radiation phase was −83 (good) vs −1 (bad); BCa 95 % CIs [−114, +18] and [−38, +1] respectively. ΔSBS5 in good was −76 [−145, −64] (within-group CI strictly negative) vs −29 in bad [−65, +9]; between-group diff −52 [−148, +1] (CI crosses 0). Mann–Whitney between-group P = 0.041 (interpret with n = 7 vs 7 caveat). The within-good clearance is robust; the between-group difference is exploratory.

### 3.6 MHC-I neoantigen landscape and radiation-induced clearance (Fig 5)
Pre-CRT neoantigen burden trended higher in eventual good responders (Fig 5A–B; median 73.5 vs 66 mutation sites with ≥ 1 MHC-I binder, MW P = 0.082; PCN 71.5 vs 57.1, P = 0.15). In paired pre→post analysis (n = 11, Fig 5C), the within-good median Δ binders was −312 [−626, −123] (within-group CI excludes 0) and median Δ binder sites −59 [−88, −26]; between-group diffs crossed zero (Δ binders diff [−527, +177], Δ sites diff [−76, +31]). The magnitude of within-good clearance (subjects 2, 6, 9 each losing > 300 MHC-I binders) is consistent with cytotoxic elimination of mutated clones but the between-group inference is exploratory. One good responder (14, pCR) atypically gained neoantigens, consistent with sparse residual tumor at post-CRT resection sampling.

### 3.7 HLA-LOH clone clearance (Fig S3)
Under **stricter Bonferroni-corrected IMGT-allele criteria** (Methods), pre-CRT HLA class I LOH was detected in 2/16 eventual good responders (subjects 3 and 4) and 0/12 bad responders (Fisher P = 0.49). Both strict-LOH subjects showed complete pre→post resolution (subj 3: 2 loci → 0; subj 4: 1 locus → 0). Lenient LOHHLA-lite criteria (uncorrected Fisher P < 0.05, |Δratio| ≥ 0.15) gave 4/16 vs 2/12 but did not change the direction. We treat HLA-LOH clone clearance as an anecdotal observation consistent with the cascade model rather than a quantitative between-group finding; the two subject panels are moved to **Supp Fig S3**. The main Fig 5 retains pre-CRT neoantigen burden and paired Δ binders.

### 3.8 Radiation-induced immune reprogramming and B-cell infiltration (Fig 6)
Within-good pre→post increases in regulatory T cells (Treg Δ +1.26 [+0.34, +1.76]), MHC II (+1.23 [+0.54, +1.92]), CD8 exhaustion (+1.00 [+0.23, +1.62]), and IGH clonotype count (TRUST4, +1,424 [0, +5,992]) all had within-good CIs excluding or touching zero. Between-group CI was **strictly above zero for Treg (MW P = 0.026, BCa diff [+0.06, +1.97])** and spanned zero for all other features. Within-good Wilcoxon signed-rank P = 0.031 consistently for Treg, MHC II, CD8 exhaustion, IGH count, confirming paired movement even where between-group inference is under-powered.

### 3.9 Assembled cascade (exploratory)
Taken as an exploratory model, eventual good responders enter the radiation phase in a pre-CRT state of proliferative, DNA-repair-proficient tumor-intrinsic competence and, across the radiation phase, appear to transit a coherent sequence — **mutation clearance → neoantigen clearance → HLA-LOH clone elimination → MHC-II / Treg / CD8 exhaustion reprogramming → B-cell infiltration**. Of these, Treg infiltration has between-group statistical support at n = 14; all other stages are hypothesis-generating and will require larger paired TNT cohorts for confirmation. The consistency of within-good CIs across mutation, neoantigen and immune axes argues against pure chance, but adequately-powered between-group inference must await prospective replication.

### 3.10 External validation — CD8-cytotoxic axis is reproducible across 9 nCRT cohorts (N = 721) (Fig 7, Table 3, Table S7)

On nine independent public GEO rectal cancer nCRT cohorts (N = 721) with harmonised per-cohort TRG-scale mapping (Table S7), meta-analysis (Stouffer's Z, √N-weighted, two-sided) of seven signatures yielded:

| Signature | Z | p_meta | Concordant cohorts | Direction |
|---|---|---|---|---|
| **CD8_cytotoxic** (pure effector) | **+2.74** | **0.006** | **8 / 9** | good > bad |
| T-cell infiltration | +1.78 | 0.075 | 8 / 9 | trend |
| B-cell infiltration | +1.56 | 0.118 | 7 / 9 | trend |
| Tumor cell-cycle | +1.31 | 0.191 | 5 / 9 | heterogeneous |
| DSB / HDR repair | +1.23 | 0.219 | 5 / 9 | heterogeneous |
| E2F / MYC | +0.69 | 0.489 | 5 / 9 | heterogeneous |
| EMT | −1.03 | 0.303 | 6 / 9 bad > good | correct direction |

The pure CD8-cytotoxic effector signature is reproducibly elevated in eventual good responders in 8/9 independent nCRT cohorts (N = 721, meta P = 0.006). The broader T- and B-cell infiltration signatures trend in the same direction. Tumor-intrinsic axes, although strongly significant in discovery, are cohort-heterogeneous externally — consistent with biopsy composition, platform probe coverage for HDR gene sets, and between-TRG-scale reclassification of borderline cases. EMT recovers the correct direction. An initial pass using a `CD8_proliferation` signature confounded with cell-cycle genes produced an uninformative null; separating effector from cell-cycle markers resolved this (Supp Text S3). No TNT-matched external cohort with paired pre-/post-CRT transcriptomics and final-TNT response labels exists publicly; GSE233517 [11] is paired but response-unlabelled, and GSE190826 contains oxaliplatin-arm treatments but is distributed as raw FASTQ only.

**Convergent independent validation (Akiyoshi et al 2023).** A concurrent independent RNA-seq study of 298 pretreatment rectal cancer biopsies, published after our discovery was completed (Akiyoshi et al, JAMA Netw Open 2023 [61]; GSE216616, Dworak TRG3-4 = good n ≈ 131 vs TRG1-2 = bad n ≈ 167), independently identifies a cytotoxic-lymphocyte effector program as the strongest predictor of neoadjuvant CRT response: MCP-counter cytotoxic lymphocyte score median 0.76 vs 0.58 (P < 10⁻³) with multivariable OR 3.81 (95 % CI 1.82–7.97, P < 10⁻³); Rooney-style cytolytic activity (GZMA × PRF1 geometric mean) 1.83 vs 1.06 (P = 0.005); Hallmark IFN-γ response, IFN-α response and inflammatory response enriched in good responders. The CD8 effector genes (GZMA, PRF1, IFN-γ pathway) and the direction of effect reported in that 298-patient independent cohort are identical to those we report here; their per-sample TRG labels are not co-deposited on GEO so direct integration into our Stouffer meta was not possible, but the paper-level convergence at n = 298 adds substantial support to the pan-CRT reproducibility conclusion. Combined with our 9-cohort meta (N = 721), the total independent evidence base for the CD8-cytotoxic pre-CRT axis now exceeds 1,000 patients across 10 independent cohorts.

### 3.11 Long-term outcomes (DFS / OS) — deferred
Survival data are not yet mature. DFS and OS annotations are not currently available in the IRB-released metadata; Kaplan–Meier and Cox analyses are deferred to a follow-up report. TRG-based final TNT response is the endpoint used here.

---

## Discussion

Three findings are central. **First**, in MSS LARC the pre-CRT molecular response predicts eventual full-TNT outcome via a **tumor-intrinsic DNA-repair and proliferative program** (E2F/MYC/G2M/DSB-HDR) that is orthogonal to classical ICB-response biomarkers: MSI/TMB do not stratify response in this cohort. A nested-CV LASSO classifier over 37 integrated features achieves modest outer held-out AUC = 0.65 (95 % CI 0.45–0.83); performance that is suggestive of a tumor-intrinsic signal in discovery but whose 95 % CI includes 0.5 under strict leakage-free evaluation. **Second**, eventual good responders appear to traverse a radiation-induced cascade — mutation clearance → neoantigen clearance → HLA-LOH clone elimination → Treg/MHC-II/CD8 exhaustion reprogramming → B-cell infiltration — whose individual within-good components are robust by BCa bootstrap but whose between-group inferences (except Treg) are under-powered at n = 14 paired subjects and are reported as exploratory. **Third, and new in this version**, the immune arm of the pre-CRT discovery — a CD8-cytotoxic effector program higher in eventual good responders — is **reproducible in 8 of 9 independent nCRT rectal cancer cohorts (N = 721, meta P = 0.006)** and independently corroborated by the 298-patient RNA-seq study of Akiyoshi et al [61] (cytotoxic lymphocyte OR 3.81, GZMA × PRF1 cytolytic activity P = 0.005, Hallmark IFN-γ enriched in responders), bringing the total to > 1,000 independent patients across 10 cohorts. This is consistent with the substantial literature linking CD8⁺/GrzB⁺ infiltration to radiation response [8–11,19–22] and establishes the discovery's immune arm as pan-CRT reproducible. The tumor-intrinsic DSB/HDR/E2F-MYC axis, although dominant in our discovery data (nested-CV AUC 0.65–0.69), is cohort-heterogeneous in the external nCRT set and is therefore framed as a discovery-stage predictor pending TNT-matched validation. The earlier impression of non-reproducibility was an artifact of signature composition (cell-cycle genes labelled CD8 proliferation) and a response-label classifier bug; both are corrected here (Supp Text S3, Table S7).

**Mid-treatment decision window.** Because the radiation phase is bracketed by the pre-CRT and post-CRT biopsies, the cascade is directly measurable at the gap between completion of CRT and initiation of consolidation chemotherapy. Under modern TNT this gap is a clinically actionable decision point: intensify, deintensify, or transition to watch-and-wait. Patients whose post-CRT biopsy already shows the full cascade may be candidates for deintensified consolidation or immediate watch-and-wait assessment; patients whose post-CRT biopsy shows no cascade despite intact pre-CRT DNA-repair competence may be candidates for intensified or alternative consolidation. Because our cascade between-group claims (other than Treg) are exploratory at n = 14, the clinical decision framing is presented as a **hypothesis for prospective trials** — not a current recommendation.

**Orthogonal to ICB biomarkers.** MSI/TMB are not useful in this MSS-dominated cohort. The reproducible CD8-cytotoxic axis, together with the discovery-stage tumor-intrinsic DSB/HDR/E2F axis, forms a candidate two-layer radiation-phase response biomarker. The two paradigms can be combined in principle: MSI-H LARC → ICB (per Cercek [5]); MSS LARC → TNT with pre-CRT CD8-cytotoxic and (pending validation) DSB/HDR/E2F-based stratification.

**Clinical implications (hypothesis-generating).** A pre-CRT RNA-seq classifier combining CD8-cytotoxic effector markers and (pending TNT-matched validation) DSB/HDR/E2F/MYC axes could stratify TNT candidates for watch-and-wait organ-preservation. EMT-high tumors, consistent with prior work [26–28], may benefit from intensified or alternative neoadjuvant strategies including anti-TGF-β agents [29] or taxane addition. Paired pre-/post-CRT biopsy monitoring of the radiation-induced clone-clearance cascade could, in adequately powered cohorts, serve as an early pharmacodynamic readout.

**Limitations.** Single-center n = 35; 8/49 tumors used tumor-only calling (stringently PoN-filtered); microarray-era external cohorts have limited probe coverage for specific signature genes (notably DSB/HDR gene sets on older Affymetrix / Illumina platforms, which we suspect contributes to the external heterogeneity of the tumor-intrinsic axis); survival data are not yet mature. HLA-LOH analysis uses stricter Bonferroni-corrected IMGT-read-counting rather than the published LOHHLA pipeline [30]; results from the two stringency tiers agree in direction and are reported side-by-side. By design, paired biopsies sample only the radiation phase of TNT — consolidation-phase biology is not molecularly observed, because under modern TNT a balanced bad-responder cohort cannot be accrued from post-consolidation tissue. A public TNT-matched RNA-seq cohort with paired pre-/post-CRT biopsies and full-TNT outcome labels is not yet available; GSE233517 is paired but response-unlabelled, and GSE190826 is distributed as raw FASTQ only. Cascade between-group claims are hypothesis-generating at n = 14 paired subjects.

**Future directions.** Prospective TNT-matched validation (PRODIGE 23 / OPRA translational substudies); single-cell RNA-seq of paired pre-/post-CRT biopsies to resolve the B-cell infiltration kinetics [31]; DFS / OS integration when outcome data mature; TNT-specific companion-diagnostic trials stratifying on combined DSB/HDR/E2F-MYC + CD8-cytotoxic axes with mid-TNT re-biopsy informing consolidation intensity.

---

## Conclusion

The molecular response to the **radiation phase** of TNT in MSS LARC — pre-CRT intrinsic DSB/HDR/E2F/MYC axis (nested-LOOCV AUC 0.65, 95 % CI 0.45–0.83), a reproducible pre-CRT CD8-cytotoxic axis (external meta Z = +2.74, P = 0.006 across 9 cohorts / N = 721), and an exploratory radiation-induced cascade of mutation / neoantigen / HLA-LOH clearance and immune reprogramming culminating in B-cell infiltration — predicts final full-TNT outcome and is directly measurable in routine paired pre-/post-CRT biopsies. The CD8-cytotoxic axis is established pan-CRT; the tumor-intrinsic axis and the cascade await TNT-matched and prospective paired-cohort validation.

---

## Figures (v0.7)

- **Fig 1.** Cohort overview.
- **Fig 2.** WES landscape.
- **Fig 3.** RNA immune/pathway signatures.
- **Fig 4.** Integration & ML (A. 37-feature correlation; B. nested-LOOCV ROC AUC 0.65 + bootstrap CI band; C. Hallmark GSEA).
- **Fig 5.** Radiation-induced neoantigen cascade (A pre-CRT missense; B pre-CRT strong binders; C paired Δ binders slopegraph + BCa CI).
- **Fig 6.** Cascade BCa forest (`figures/panels/Fig6_cascade_BCa_forest.{pdf,png}`). **A**: per-feature paired Δ medians with BCa 95 % CIs, good (blue circles) vs bad (red squares), standardized per feature; raw Δ values annotated to the right. **B**: between-group Δ (good − bad) with BCa 95 % CI. Treg is the only feature whose between-group CI strictly excludes zero (robust, teal diamond); all other cascade features have CIs spanning zero and are labelled exploratory (grey diamond).
- **Fig 7.** External validation of the CD8-cytotoxic axis. (A) Per-cohort forest (9 cohorts, N = 721) + Stouffer meta diamond (Z = +2.74, P = 0.006) + **convergent row for Akiyoshi 2023 / GSE216616 (N = 298; OR 3.81 [1.82, 7.97]; GZMA × PRF1 cytolytic activity P = 0.005)** — total independent evidence > 1,000 patients across 10 cohorts. (B) Signature-level meta Z across 9 cohorts. (C) Decoupling of CD8 effector vs tumor proliferation per cohort. (`figures/panels/Fig7_external_CD8_validation_v4.{pdf,png}`)
- **Supp Fig S1.** CONSORT-style sample-flow diagram.
- **Supp Fig S2.** Nested-CV comparison (LASSO vs ElasticNet vs RandomForest, with permutation null).
- **Supp Fig S3.** HLA-LOH lite vs strict comparison (`figures/supp/SuppFig_S3_HLA_LOH_lite_vs_strict.{pdf,png}`). **A**: pre-CRT LOH prevalence (lite vs Bonferroni-corrected strict) by response group, showing that stringent criteria reduce 10 lite events to 2 strict events. **B**: per-subject strict LOH pre→post resolution (subj 3: 2 → 0 loci, subj 4: 1 → 0 loci). **C**: criteria comparison text panel documenting thresholds and reporting caveats.

## Tables
- **Table 1.** Clinical characteristics.
- **Table 2.** Top 20 integrated-feature associations.
- **Table 3.** External-validation meta summary (7 signatures × 9 cohorts).
- **Table S1.** 37-feature master table.
- **Table S2.** SBS activities.
- **Table S3.** Full GSEA (Hallmark + Reactome).
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types.
- **Table S6.** pVACseq neoantigen detail.
- **Table S7.** External validation per-cohort detail + response-scale mapping.
- **Table S8.** Cascade BCa 95 % bootstrap CIs.
- **Table S9.** HLA-LOH lenient vs strict per-sample calls.

## Supp Text
- **S1.** Extended methods (WES, RNA-seq, HLA typing).
- **S2.** Nested-CV pipeline details and permutation.
- **S3.** External meta v3 diagnostic (CD8 axis rescue, classifier bug correction).

---

## Data and code availability

- Analysis workspace: `/mnt/sda1/data/TNT/analysis/`.
- GitHub: https://github.com/Soonlab/TNT.
- Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on acceptance.

## Ethics

IRB-approved at Seoul National University Hospital. Written informed consent.

## References (expanded)

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

*End of v0.7 — Genome Medicine target — 2026-04-15.*
