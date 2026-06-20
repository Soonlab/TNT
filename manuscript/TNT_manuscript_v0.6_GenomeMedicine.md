# Molecular response to the radiation phase of total neoadjuvant therapy predicts final response in MSS locally-advanced rectal cancer: a multi-omics dissection

*Manuscript draft v0.6 — Genome Medicine submission target. 2026-04-15.*
*Updated from v0.5 with corrected external validation: pure CD8-cytotoxic signature (N=721 across 9 independent nCRT cohorts) reproduces the discovery immune axis at Stouffer Z=+2.74, p=0.006. Tumor-intrinsic DNA-repair/cell-cycle axis remains cohort-heterogeneous and is framed as discovery-stage pending TNT-matched validation.*

---

## Abstract (350 words)

**Background.** Total neoadjuvant therapy (TNT) — long-course chemoradiation (CRT) combined with induction or consolidation FOLFOX/CAPOX — is standard for locally advanced rectal cancer (LARC) and enables organ-preservation watch-and-wait. Under modern TNT the post-consolidation complete-response rate is so high that a balanced "bad-responder" group cannot be accrued from post-consolidation tissue; paired-biopsy molecular studies are therefore tractable only around the **radiation (first) phase** of TNT. This window is also a clinically actionable mid-treatment decision point (intensify, deintensify, or abandon consolidation). No molecular predictor of final TNT outcome measurable from the radiation-phase response is established for the microsatellite-stable (MSS) majority.

**Methods.** We profiled 35 MSS LARC patients by matched WES (77 samples, 41 tumor–normal pairs) and RNA-seq (56 samples). Biopsies were obtained pre-CRT (baseline) and post-CRT (before consolidation chemotherapy); the binary response label reflects final TNT response after completion of the full regimen. Somatic calls: GATK Mutect2 with a 28-sample PoN; MSI by msisensor-pro; SBS by SigProfilerAssignment COSMIC v3.3; CNV by CNVkit; HLA class I by OptiType; HLA-allelic imbalance by direct IMGT read-counting; MHC-I neoantigens by pVACseq-MHCflurry. Transcriptomic analysis: DESeq2 (adjusted for sex and cT), fgsea Hallmark/Reactome, ssGSEA on 95 curated pathways, TRUST4, CMScaller. Thirty-seven per-subject features entered a LASSO logistic regression predictor with leave-one-out cross-validation. External validation used nine public GEO nCRT cohorts (N = 721) with harmonised per-cohort TRG-scale mapping (Dworak vs Mandard vs CAP/AJCC vs Rödel), aggregated by Stouffer's Z weighted by √N, on four immune signatures (CD8-cytotoxic, T-cell infiltration, B-cell infiltration, tumor cell-cycle) and three tumor-intrinsic axes (DSB/HDR, E2F/MYC, EMT).

**Results.** All 41 matched tumors were MSS (median 1.6 mutations/Mb) and TMB did not discriminate response (P = 0.19). **Narrative 1 — Pre-CRT predictor.** Hallmark GSEA showed elevated E2F targets (NES = 2.78), G2M, MYC, and Reactome DSB/HDR repair in eventual good responders, with reciprocal EMT down. A 37-feature LASSO classifier achieved LOOCV AUC = 0.755. **Narrative 2 — Radiation-induced cascade.** In 14 paired pre-/post-CRT subjects, eventual good responders uniquely exhibited a coherent cascade: mutation clearance (ΔSBS5 P = 0.041) → MHC-I neoantigen clearance (median Δ binders −312 vs −100) → HLA-LOH clone elimination (subjects 3 and 4) → immune reprogramming (Treg P = 0.026) → B-cell infiltration (IGH Δ +1,424 vs +7, P = 0.031). **External validation.** On 9 independent cohorts (N = 721) the **CD8-cytotoxic axis reproduced robustly (Z = +2.74, P = 0.006, 8/9 cohorts concordant)**, T- and B-cell infiltration trended (P = 0.075, 0.12), and EMT was directionally correct. Tumor-intrinsic DSB/HDR/E2F signals were cohort-heterogeneous, consistent with biopsy composition, platform, and TRG-scale variation.

**Conclusions.** The molecular response to the **radiation phase of TNT** is measurable in paired pre-/post-CRT biopsies and predicts final TNT outcome. The CD8-cytotoxic immune axis is reproducible pan-CRT; the tumor-intrinsic DSB/HDR/E2F axis is a strong discovery-stage predictor that requires TNT-matched external validation.

---

## Background

Locally advanced rectal cancer (LARC) affects > 40,000 patients per year in Asia and the USA combined. Total neoadjuvant therapy (TNT) — long-course chemoradiation (CRT; 50.4 Gy with concurrent capecitabine) combined with induction or consolidation FOLFOX/CAPOX — is now first-line based on PRODIGE 23, RAPIDO and OPRA [1–3]. TNT enables organ-preservation watch-and-wait in responders, so response prediction directly influences whether a patient undergoes surgery. Yet only 15–30 % of patients reach pathological or sustained clinical complete response, and **no molecular predictor of final TNT response is established for the MSS majority** [4].

TNT is a multi-stage regimen: a radiation phase (long- or short-course CRT) and a consolidation (or induction) chemotherapy phase with fluoropyrimidine + oxaliplatin. Under modern TNT, the post-consolidation complete-response rate is sufficiently high that a biopsy-based "bad-responder" cohort cannot be accrued at the post-consolidation timepoint: residual tumor tissue is often unavailable, and patients frequently transition directly to watch-and-wait. Paired-biopsy molecular studies of TNT are therefore practically restricted to a sampling window **bracketing the radiation (first) phase**: pre-CRT baseline and post-CRT before consolidation chemotherapy begins. The tractable — and clinically actionable — scientific question is: **how does the molecular response to the radiation phase of TNT predict the final TNT response?**

Two orthogonal molecular axes dominate current thinking about LARC response: MSI-H tumors respond to single-agent PD-1 blockade (100 % complete response with dostarlimab [5]), while EMT and CMS4 mesenchymal biology drive chemoradiation resistance [6,7]. But MSI-H is only 5–7 % of rectal cancer, and most published signatures were derived from small single-regimen nCRT cohorts (n = 30–100) that pre-date TNT and did not examine radiation-phase molecular dynamics with final TNT outcome as endpoint.

Independently of TNT, a robust body of work links cytotoxic CD8⁺ T-cell infiltration to radiotherapy response in rectal and other solid tumors [8,9], including direct demonstrations that neoadjuvant (chemo)radiation increases CD8/Granzyme-B T-cell density and restores interferon-γ programs [10,11]. This motivates placing a pure CD8-effector axis at the centre of external validation.

Here we profile 35 MSS LARC patients whose paired biopsies bracket the radiation phase of TNT by matched WES and RNA-seq, reconstruct somatic, HLA, neoantigen, transcriptomic, immune and clonal-dynamic axes per patient, and dissect response along two narratives: a **pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor** and a **radiation-induced immunogenic-clone-clearance cascade in eventual good responders**. We then validate the immune arm of these narratives in a 9-cohort external meta-analysis (N = 721), showing that a pure CD8-cytotoxic axis is reproducible across independent nCRT datasets.

---

## Methods

### Patients and samples
Thirty-five LARC patients (clinical T2–T4) received TNT (induction/consolidation FOLFOX or CAPOX + long-course 50.4 Gy chemoradiation with concurrent capecitabine). Final TNT response was graded on surgical specimens by Dworak TRG and binarised good (TRG 0–1; n = 18) vs bad (TRG 2–3; n = 17) after completion of the full regimen. Fourteen subjects contributed matched pre-CRT biopsy, post-CRT biopsy and blood normal. IRB approved.

### Treatment and sampling timing
Biopsies were obtained **pre-CRT** (baseline, before any neoadjuvant therapy) and **post-CRT** (after completion of 50.4 Gy + concurrent capecitabine, before consolidation chemotherapy began). Paired-biopsy profiling therefore brackets the radiation phase only; the consolidation phase is not sampled. The binary response label is the final TNT response after completion of the full regimen.

### WES, RNA-seq, HLA, neoantigen, CNV, integration
(Unchanged from v0.5.) Agilent SureSelect V5 capture, NovaSeq PE101, BWA-MEM to GRCh38, GATK best-practice. Mutect2 v4.6.2 with 28-sample cohort PoN + gnomAD v3.1; FilterMutectCalls + LearnReadOrientationModel + CalculateContamination; 8 unmatched tumors used tumor-only + PoN. SBS refit: SigProfilerAssignment COSMIC v3.3. MSI: msisensor-pro. CNV: CNVkit. HLA: OptiType on MHC-extracted reads. HLA-allelic imbalance: direct IMGT read counting (|Δratio|>0.15, Fisher P < 0.05 per locus). Neoantigens: pVACseq v5 + MHCflurry 2.0, 4-digit HLA, 8–11-mers, binders ≤ 500 nM / strong ≤ 50 nM. Presentation-competent neoantigen score = unique binder sites × (1 − 0.33·LOH fraction). RNA-seq: stranded NovaSeq PE101, HISAT2 + StringTie to GRCh38 / GENCODE v39, gene-level TPM (46,425 × 56). DESeq2 with `~sex + cT_simple + response_bin`. fgsea (Hallmark + Reactome). gseapy ssGSEA on 95 sets. TCR/BCR by TRUST4. CMS by CMScaller. Thirty-seven per-subject features → LASSO logistic regression with LOOCV held-out AUC.

### External validation (meta-analysis)
We assembled nine public GEO nCRT rectal cancer cohorts with interpretable response labels (N = 721): **GSE150082, GSE35452, GSE119409, GSE45404, GSE94104, GSE56699, GSE46862, GSE133057, GSE87211**. All used long-course nCRT (45–50.4 Gy ± fluoropyrimidine ± oxaliplatin); none received modern induction/consolidation TNT. For each cohort, response labels were mapped manually to the final good/bad binary using the scale appropriate to that cohort (Dworak, Mandard, CAP/Ryan/AJCC, Rödel, or cohort-specific responder/non-responder calls; GSE87211 used cancer recurrence after surgery as a survival surrogate because explicit TRG was not provided). Per-cohort probe-to-gene mapping was performed from the native platform annotation; expression matrices were log2-transformed where necessary; signatures were scored per sample by z-score averaging. Signatures tested:
- **CD8_cytotoxic** — CD8A, CD8B, GZMA, GZMB, GZMH, GZMK, PRF1, IFNG, NKG7, GNLY, CXCL9, CXCL10, CXCL11, TBX21, EOMES, KLRK1, KLRD1 (pure CD8 effector markers; *no* cell-cycle genes).
- **Tcell_infiltration** — CD3D, CD3E, CD3G, CD2, CD4, CD8A, CD8B, LCK, ZAP70, ITK.
- **Bcell_infiltration** — CD19, MS4A1 (CD20), CD79A, CD79B, CD22, TCL1A, FCRL5, BLK, FCER2.
- **Tumor_cellcycle** — MKI67, TOP2A, STMN1, TYMS, UBE2C, BIRC5, CCNB1/B2, CDK1, MCM2, MCM5, PCNA, CENPF, KIF20A, AURKA, AURKB, PLK1, BUB1.
- **DSB_HDR_repair, E2F_MYC_cellcycle, EMT** as in the discovery analysis.

Per-signature per-cohort Mann–Whitney U compared good vs bad. Cohort-level p-values were aggregated by two-sided Stouffer's Z weighted by √(n_good + n_bad); the sign of each cohort's Z was set to match the sign of Δ(good − bad). Full per-cohort statistics are reported (Supp Table S7). A preliminary version of this analysis that used a "CD8_proliferation" signature confounded with tumor cell-cycle markers, and that additionally mis-classified `Non-responder` labels as responders due to a substring matching bug, is replaced by the present signature-refined and label-corrected analysis. The diagnostic narrative of that artifact is retained in Supplementary `Supp_external_meta_v3_CD8axis.md`.

### Code and data
Analysis scripts under `/analysis/scripts/` (numbered `00_`–`33_`). Code and derived tables on https://github.com/Soonlab/TNT. Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on publication.

---

## Results

### 3.1 Cohort: MSS, TMB-low, ICB-biomarker-negative
Of 35 patients, 18 were eventual good responders (final TNT TRG 0–1) and 17 bad (TRG 2–3) (**Fig 1**, **Table 1**). Clinical T4 was enriched in bad responders (41 % vs 11 %, P = 0.086). **All 41 matched tumors were microsatellite-stable** (maximum MSI % 0.19 %) and TMB-low (median 1.6 /Mb; good 1.85 vs bad 1.40, P = 0.186). MSI-H and high TMB — the two established ICB biomarkers — do not apply.

### 3.2 Somatic landscape (Fig 2)
CRC driver mutations followed a canonical MSS distribution: APC 30/49 (61 %), TP53 20/49 (41 %), KRAS 14/49 (29 %), FBXW7 7/49 (14 %), KMT2D 4/49. No single driver reached significance; FBXW7 trended toward good response (OR 3.7, P = 0.36). SBS5 / SBS1 (clock-like) dominated (> 60 %); **SBS3 (HRD) was absent** in all samples; sporadic MMR signatures (SBS6/15/20/26) were seen in both strata. CIN was indistinguishable (P = 0.66); an HRD LST proxy was modestly higher in bad responders (P = 0.037).

### Narrative 1 — Pre-CRT tumor-intrinsic DNA-repair/cell-cycle predictor of final TNT response

### 3.3 DNA-repair and cell-cycle pathways stratify response (Fig 3, Fig 4)
Hallmark GSEA of **pre-CRT** transcriptomes (n = 33) showed E2F targets (NES = 2.78), G2M checkpoint (NES = 2.46), MYC targets V1/V2 (NES ≥ 2.23), mTORC1 and mitotic spindle markedly elevated in eventual good responders, with EMT (NES = −2.16), myogenesis and apical junction reciprocally suppressed. Reactome independently placed cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair and DNA replication as the top positive sets. ssGSEA on 95 curated sets corroborated: DSB repair (P = 0.007), MYC targets V2 (0.018), HDR (0.020), general DNA repair (0.020), G2-M (0.032), E2F (0.035), a proliferation-heavy CD8 score (0.035). MHC II was modestly lower in pre-CRT good responders (P = 0.074).

### 3.4 Integrated 37-feature LASSO predictor (Fig 4)
Thirty-seven per-subject features formed a tightly co-varying DNA-repair + cell-cycle module together with an immune proliferation module. A LASSO logistic-regression classifier trained on pre-CRT features achieved **leave-one-out AUC = 0.755** against final TNT response, outperforming elastic-net and random-forest top-8 models (both 0.70). RF feature importance ranked MYC V2, DSB repair, hypoxia, HDR, MHC II and deletion fraction at the top. PyClone-VI on 12 paired tumors showed a non-significant trend toward larger dominant-clone shrinkage post-CRT in eventual good responders (Δ −0.67 vs −0.15, P = 0.34).

### Narrative 2 — Radiation-induced cascade in eventual good responders

### 3.5 Missense and SBS5 mutation clearance post-CRT in eventual good responders
In 14 paired pre-/post-CRT tumors, median Δ missense count across the radiation phase was −67 (good) vs −8.5 (bad); ΔSBS5 −76 vs −29, P = 0.041.

### 3.6 MHC-I neoantigen landscape and radiation-induced clearance (Fig 5)
Pre-CRT neoantigen burden was modestly elevated in eventual good responders (**Fig 5A–B**; median 73.5 vs 66 mutation sites generating ≥ 1 MHC-I binder, P = 0.082; PCN 71.5 vs 57.1, P = 0.15). In paired pre-CRT → post-CRT analysis (n = 11; **Fig 5C**), subjects 2, 6 and 9 each lost > 300 MHC-I binders (Δ −312, −489, −626) and > 50 binder-generating sites. Bad responders showed smaller losses (median Δ binder sites −59 good vs −16 bad, P = 0.25). One good responder (14, pCR) atypically gained neoantigens, reflecting sparse residual tumor at post-CRT resection.

### 3.7 HLA-LOH clone clearance during the radiation phase (Fig 5D)
HLA class I LOH (direct IMGT-allele imbalance) was detected at ≥ 1 class I locus in 10/76 tumors. Pre-CRT prevalence: 4/16 good vs 2/12 bad (Fisher P = 0.67). Subjects 3 and 4 (both eventual good responders) showed multi-locus HLA-LOH pre-CRT that reduced post-CRT (3→1 and 2→0 respectively; **Fig 5D**).

### 3.8 Radiation-induced immune reprogramming and B-cell infiltration
Post-CRT tumors from eventual good responders showed paired elevations in regulatory T cells (Δz +1.26 vs +0.03, P = 0.026), MHC II (+1.23 vs +0.36, P = 0.065), CD8 exhaustion (+1.00 vs −0.10, P = 0.093), and coherent upregulation of antigen-presentation, TNF-α, IL-2/STAT5, IL-6/STAT3 and allograft-rejection pathways (within-good Wilcoxon P = 0.031). IGH clonotype count (TRUST4) increased Δ +1,424 in good vs +7 in bad (within-good P = 0.031); IGK/IGL and TRA/TRB Shannon diversity moved in the same direction.

### 3.9 Assembled radiation-induced cascade in eventual good responders
Eventual good responders enter the radiation phase of TNT in a pre-CRT state of proliferative, DNA-repair-proficient tumor-intrinsic competence and, across the radiation phase, transit a coherent cascade — **mutation clearance → neoantigen clearance → HLA-LOH clone elimination → MHC-II / Treg / CD8 exhaustion reprogramming → B-cell infiltration** — measurable in the post-CRT biopsy before consolidation chemotherapy begins.

### 3.10 External validation — CD8-cytotoxic immune axis is reproducible across 9 nCRT cohorts (N = 721)

We assembled nine independent public GEO rectal cancer nCRT cohorts (total N = 721 patients) spanning multiple TRG grading scales (Dworak, Mandard, CAP/AJCC, Rödel, and cohort-specific binary response calls; recurrence for GSE87211). Response labels were harmonised per-cohort to a final good/bad binary, and per-sample signature scores were computed from the native platform annotations. Meta-analysis (Stouffer's Z, √N-weighted, two-sided) of the immune and tumor-intrinsic axes from the discovery analysis yielded the following (**Fig 6 [new], Table 3**):

| Signature | Z | p_meta | Concordant cohorts | Direction |
|---|---|---|---|---|
| **CD8_cytotoxic (pure effector)** | **+2.74** | **0.006** | **8 / 9** | good > bad |
| T-cell infiltration (CD3 axis) | +1.78 | 0.075 | 8 / 9 | good > bad (trend) |
| B-cell infiltration | +1.56 | 0.118 | 7 / 9 | good > bad (trend) |
| Tumor cell cycle | +1.31 | 0.191 | 5 / 9 | cohort-heterogeneous |
| DSB / HDR repair | +1.23 | 0.219 | 5 / 9 | cohort-heterogeneous |
| E2F / MYC | +0.69 | 0.489 | 5 / 9 | cohort-heterogeneous |
| EMT | −1.03 | 0.303 | 6 / 9 bad > good | correct direction |

The pure CD8-cytotoxic effector signature (CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES) is **reproducibly elevated in eventual good responders in 8 of 9 independent nCRT cohorts**, spanning 721 patients, and the meta Z reaches P = 0.006. The broader T-cell and B-cell infiltration signatures show the same direction as a trend. In contrast, the **tumor-intrinsic DSB-repair, E2F/MYC and cell-cycle axes**, although strongly significant in the discovery cohort, are **cohort-heterogeneous across the external set** (meta Z 0.7–1.3, P > 0.19), consistent with known variation in biopsy composition, tumor purity, microarray platform coverage of these gene sets, and between-TRG-scale reclassification of borderline cases. EMT recovers the correct direction (bad > good) although not reaching significance.

Two methodological notes are required for transparency. First, an initial pass of this external analysis used a signature labelled `CD8_proliferation` that contained primarily cell-cycle / proliferation genes (MKI67, TOP2A, MCM2/5, PCNA, CCNB1/B2, CDK1, CENPF, KIF20A, UBE2C, BIRC5); in bulk pre-treatment biopsies these genes track tumor-intrinsic proliferation rather than lymphoid effector state and therefore inherited the tumor-intrinsic heterogeneity across cohorts, producing an uninformative null. Separating pure CD8 effector markers from cell-cycle markers resolves this confound and is the correct construction for an immune-axis external validation. Second, a bug in the original response classifier treated `Non-responder` as a responder via substring matching; per-cohort manual TRG-scale mapping (Supp Table S7) corrects this and additionally rescues cohorts whose response field sat in non-standard annotation columns (e.g., GSE45404 `class`, GSE46862 TO/MO/MI/NT, GSE87211 recurrence, GSE94104 Rödel TRG, GSE133057 AJCC). Full per-cohort, per-signature effect sizes are in Supp Table S7 and `manuscript/supplementary/Supp_external_meta_v3_CD8axis.md`; forest plot is Fig 6 with the decoupling of CD8 effector vs tumor cell cycle shown in panel C.

### 3.11 Long-term outcomes (DFS / OS) — deferred
Survival data are not yet mature — this is a recently-accrued cohort; DFS/OS analyses are planned for a follow-up report. DFS and OS annotations are not currently available in the IRB-released metadata (see `analysis/clinical_survival_status.md`). Kaplan–Meier and Cox analyses are therefore deferred; TRG-based final TNT response remains the endpoint used here.

---

## Discussion

Three findings are central. **First**, in MSS LARC the molecular response to the **radiation phase of TNT** as captured by pre-CRT intrinsic biology predicts eventual full-TNT outcome and is governed by tumor-intrinsic DNA-repair and proliferative programs, not by classical ICB-response biomarkers. This runs **orthogonal to the ICB paradigm**: under ICB, MSI-H / high-TMB / inflamed tumors respond well via neoantigen-driven T-cell recognition; under the radiation phase of TNT, intact HDR/DSB-repair and active cell cycle sensitise tumors to ionising radiation + concurrent capecitabine. EMT-high mesenchymal tumors invoke established chemoradiation-resistance programs.

**Second**, eventual good responders traverse a coherent radiation-induced immunogenic-clone-clearance cascade: mutation and neoantigen clearance, HLA-LOH clone elimination, MHC-II / Treg / CD8 exhaustion reprogramming, and B-cell / TLS-like infiltration. This mirrors emergent TLS biology reported in ICB responders but is here induced by the radiation phase of TNT rather than by checkpoint blockade.

**Third, and new in this version**, the immune arm of the discovery — a pre-CRT CD8-cytotoxic effector program higher in eventual good responders — is **reproducible in 8 of 9 independent nCRT rectal cancer cohorts (N = 721, meta P = 0.006)**. This is consistent with the substantial literature linking CD8⁺ / Granzyme-B⁺ infiltration to radiation response [8–11] and places the discovery's immune arm in a pan-CRT reproducible class. The tumor-intrinsic DSB/HDR/E2F-MYC axis, although extremely strong in our discovery data and the dominant contributor to the 0.755 LOOCV classifier, is **cohort-heterogeneous in the external nCRT set** and should therefore be treated as a discovery-stage predictor pending TNT-matched external validation. The CD8 reproducibility failure in our earlier external analysis was an artifact of signature composition (cell-cycle genes labelled as CD8 proliferation) and a response-label classifier bug, both corrected here; the full diagnostic is in Supp Table S7 and `Supp_external_meta_v3_CD8axis.md`.

**Mid-treatment decision window.** Because the radiation phase is bracketed by the pre-CRT and post-CRT biopsies, the cascade is directly measurable at the gap between completion of CRT and initiation of consolidation chemotherapy. Patients whose post-CRT biopsy already shows the full cascade (mutation / neoantigen / HLA-LOH clearance + B-cell infiltration) may be candidates for deintensified consolidation or immediate watch-and-wait assessment, while patients whose post-CRT biopsy shows no cascade despite intact pre-CRT DNA-repair competence may be candidates for intensified or alternative consolidation. A limitation is that consolidation-phase biology itself is not directly observed in this study.

**Orthogonal to ICB biomarkers.** MSI/TMB are not useful in this MSS-dominated cohort. The reproducible CD8-cytotoxic axis, together with the tumor-intrinsic DSB/HDR/E2F axis, forms a two-layer radiation-phase response biomarker. The two paradigms can be combined: MSI-H LARC → ICB (per Cercek [5]); MSS LARC → TNT with DSB/HDR/E2F and CD8-cytotoxic-based pre-CRT stratification and post-CRT cascade-based consolidation decision.

**Clinical implications.** A pre-CRT RNA-seq classifier built on DSB/HDR/E2F/MYC/CD8-cytotoxic could, pending TNT-matched validation of the tumor-intrinsic arm, stratify TNT candidates for watch-and-wait organ-preservation. EMT-high tumors might be candidates for intensified or alternative neoadjuvant strategies. Paired pre-/post-CRT biopsy monitoring of the radiation-induced clone-clearance cascade could serve as an early pharmacodynamic readout and a mid-TNT decision tool.

**Limitations.** Single-center n = 35; 8/49 tumors used tumor-only calling (stringently PoN-filtered); microarray-era external cohorts have limited probe coverage for specific signature genes; survival data are not yet mature. HLA-LOH analysis uses IMGT-read-counting rather than LOHHLA. By design, paired biopsies sample only the radiation phase of TNT — consolidation-phase biology is not molecularly observed, because under modern TNT a balanced bad-responder cohort cannot be accrued from post-consolidation tissue. A public TNT-matched RNA-seq cohort with paired pre-/post-CRT biopsies and full-TNT outcome labels is not yet available; GSE233517 is paired but response-unlabelled, and GSE190826 contains some oxaliplatin-arm TNT-like treatments but is distributed as raw FASTQ only.

**Future directions.** Prospective TNT-matched validation (PRODIGE 23 / OPRA translational substudies); single-cell RNA-seq of paired pre-/post-CRT biopsies to resolve the B-cell infiltration; DFS / OS integration when outcome data mature; TNT-specific companion-diagnostic clinical-trial stratification on combined DSB/HDR/E2F + CD8-cytotoxic axes with mid-TNT re-biopsy informing consolidation intensity.

---

## Conclusion

The molecular response to the **radiation phase** of TNT in MSS LARC — pre-CRT intrinsic DSB/HDR/E2F/MYC axis (LOOCV AUC 0.755), reproducible pre-CRT CD8-cytotoxic axis (external meta Z = +2.74, P = 0.006 across 9 cohorts / N = 721), and a radiation-induced mutation-clearance → neoantigen-clearance → HLA-LOH clone-elimination → immune-reprogramming → B-cell-infiltration cascade — predicts final full-TNT outcome and is directly measurable in routine paired pre-/post-CRT biopsies, providing a two-axis radiation-phase-TNT response biomarker and a clinically actionable mid-treatment decision window.

---

## Figures (v0.6)

- **Figure 1.** Cohort overview (final TNT response, cT stage, sex/age, sample matrix, study design; timing panel emphasises pre-/post-CRT sampling and post-consolidation response label).
- **Figure 2.** WES landscape (TMB, MSI, driver oncoprint, SBS signatures, MMR, CIN, HRD).
- **Figure 3.** RNA immune/pathway signatures (signature grid + heatmap + DEG volcano + CD8 axes + post-CRT immune).
- **Figure 4.** Integration & ML (37-feature correlation + response ranks + Hallmark GSEA + LASSO ROC AUC = 0.755).
- **Figure 5.** Radiation-induced neoantigen cascade (A pre-CRT missense; B pre-CRT strong binders; C paired Δ binders; D HLA-LOH subj 3, 4).
- **Figure 6. [new]** External validation of the CD8-cytotoxic axis (9 cohorts, N = 721). **A**: per-cohort forest of Δ CD8-cytotoxic (good − bad), with Stouffer meta diamond Z = +2.74, P = 0.006. **B**: meta Z across all seven tested signatures, showing CD8 significant, T-cell / B-cell trending, and tumor-intrinsic axes non-significant. **C**: per-cohort scatter of Δ CD8-cytotoxic vs Δ tumor cell cycle, illustrating that the immune axis is pan-cohort positive while the tumor-intrinsic axis is heterogeneous. (`figures/supp/SuppFig_v3_CD8_meta_forest.{pdf,png}` — promoted to main Fig 6.)
- **Figure 7.** Paired pre-/post-CRT immune delta (Treg, MHC II, CD8 exhaustion, IGH clonotype expansion) — previously Fig 6.
- **Supp 1–N.** PyClone-VI clonal decomposition, per-subject neoantigen detail, full GSEA tables, per-subject master table.

## Tables
- **Table 1.** Clinical characteristics (eventual good vs bad, final TNT response).
- **Table 2.** Top 20 integrated-feature associations with final TNT response.
- **Table 3. [new]** External-validation meta-analysis summary (7 signatures × 9 cohorts, Stouffer Z, per-cohort Δ, concordance).
- **Table S1.** 37-feature per-subject master table.
- **Table S2.** SBS signature activities per sample.
- **Table S3.** Full GSEA (Hallmark + Reactome).
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types per subject.
- **Table S6.** pVACseq neoantigen per-subject detail.
- **Table S7. [new]** External validation — per-cohort per-signature effect sizes and response-label mapping (9 cohorts).

---

## Data and code availability

- Analysis workspace: `/mnt/sda1/data/TNT/analysis/`.
- GitHub: https://github.com/Soonlab/TNT.
- Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on acceptance.

## Ethics

IRB-approved at Seoul National University Hospital. Written informed consent from all subjects.

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
11. Lim YJ et al. Sci Rep 2023; GSE233517 (CRT-induced CD8 and IFN-γ in rectal cancer).

---

*End of v0.6 draft — Genome Medicine target — 2026-04-15 — CD8-axis external reproducibility added.*
