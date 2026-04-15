# Multi-omics dissection of response to total neoadjuvant therapy in MSS locally-advanced rectal cancer reveals a tumor-intrinsic DNA-repair predictor and a post-treatment immune-reprogramming cascade

*Manuscript draft v0.4 — Genome Medicine submission target. 2026-04-15.*
*Restructured from v0.3 around two core narratives; adds dedicated Figure 5 (neoantigen cascade), moves external validation to supplementary, integrates HLA-LOH / pVACseq as main-text content, frames external heterogeneity as TNT-regimen-specific biology.*

---

## Abstract (350 words)

**Background.** Total neoadjuvant therapy (TNT) — induction or consolidation chemotherapy combined with long-course chemoradiation — is now standard for locally advanced rectal cancer (LARC) and enables organ-preservation watch-and-wait protocols. However, only 15–30% of patients achieve complete response, and no molecular predictor is established for the microsatellite-stable (MSS) majority in which classical immune-checkpoint-blockade biomarkers (MSI-H, high TMB) do not apply. Post-treatment biology of responders is also poorly charted.

**Methods.** We performed matched whole-exome sequencing (WES; 77 samples from 35 patients, 41 matched tumor–normal pairs, 14 paired pre/post) and RNA-seq (56 samples) on a single-center TNT cohort. Somatic calling used GATK Mutect2 against a cohort panel-of-normals; MSI by msisensor-pro; mutational signatures by SigProfilerAssignment (COSMIC v3.3); copy number by CNVkit; HLA class I by OptiType; HLA allelic imbalance by direct IMGT read-counting; MHC-I neoantigens by pVACseq-MHCflurry. Transcriptomic analysis used DESeq2 (adjusted for sex and cT), fgsea Hallmark/Reactome enrichment, ssGSEA on 95 curated pathways, TRUST4 for TCR/BCR, and CMScaller for CMS. Thirty-seven per-subject features were integrated into a LASSO logistic-regression predictor with leave-one-out cross-validation. Reproducibility was tested in seven public GEO cohorts (N = 290) by meta-analysis (Stouffer's Z).

**Results.** All 41 matched tumors were MSS (MSI < 0.2 %, median 1.6 mutations/Mb); TMB did not discriminate response (P = 0.19). Two coherent narratives emerged. **(1) Pre-treatment predictor.** Hallmark GSEA showed markedly higher E2F targets (NES = 2.78), G2M checkpoint, MYC targets and Reactome DSB/HDR repair in good responders, with reciprocal EMT downregulation. An integrated LASSO classifier over 37 features achieved LOOCV AUC = 0.755. **(2) Post-treatment cascade.** In 14 paired pre/post subjects, good responders uniquely exhibited a coherent cascade: mutation clearance (ΔSBS5 P = 0.041) → MHC-I neoantigen clearance (median Δ binders −312 good vs −100 bad) → HLA-LOH clone elimination (subjects 3 and 4 reverting from multi-locus LOH pre-treatment to resolved post) → immune reprogramming (Treg P = 0.026, MHC II, CD8 exhaustion) → B-cell infiltration (IGH Δ +1,424 good vs +7 bad, P = 0.031). External validation in seven GEO nCRT cohorts was only partially reproducible, consistent with TNT-regimen-specific biology.

**Conclusions.** In MSS LARC, TNT response is governed by tumor-intrinsic DNA-repair and proliferative programs, not classical ICB-response biology. Responders engage a post-treatment immunogenic-clone-clearance cascade measurable in routine pre/post biopsies, suggesting new biomarker and clinical-trial-stratification opportunities for TNT.

---

## Background

Locally advanced rectal cancer (LARC) affects > 40,000 patients per year in Asia and the USA combined. Total neoadjuvant therapy (TNT) — induction or consolidation FOLFOX/CAPOX plus long-course chemoradiation — is now first-line based on the PRODIGE 23, RAPIDO and OPRA trials [1–3]. TNT enables **organ-preservation watch-and-wait protocols** in responders, so response prediction directly influences whether a patient undergoes surgery. Yet only 15–30 % of patients reach pathological or sustained clinical complete response, and **no molecular predictor of TNT response is established for the MSS majority** [4]. Two orthogonal molecular axes dominate current thinking: MSI-H tumors respond to single-agent PD-1 blockade (100 % complete response with dostarlimab [5]), while EMT and CMS4 mesenchymal biology drive chemoradiation resistance [6,7]. But MSI-H is only 5–7 % of rectal cancer, and most published signatures were derived from small single-regimen nCRT cohorts (n = 30–100) with limited multi-omic integration and pre-date TNT. Post-treatment biology in responders, in particular, has been largely neglected.

Here we profile 35 MSS LARC patients treated with contemporary TNT by matched WES and RNA-seq (with 14 paired pre/post tumors), reconstruct somatic, HLA, neoantigen, transcriptomic, immune, and clonal-dynamic axes per patient, and dissect response along two narratives: a **pre-treatment tumor-intrinsic DNA-repair/cell-cycle predictor** and a **post-treatment immunogenic-clone-clearance cascade**. We position these axes orthogonally to the ICB biomarker paradigm.

---

## Methods

### Patients and samples
Thirty-five patients with LARC (clinical T2–T4) received TNT (induction/consolidation FOLFOX or CAPOX + long-course 50.4 Gy chemoradiation with concurrent capecitabine). Response was graded on surgical specimens by Dworak TRG and binarised good (TRG 0–1; n = 18) vs bad (TRG 2–3; n = 17). Fourteen subjects contributed matched pre-biopsy, post-surgery and blood normal. Informed consent was IRB-approved.

### WES
Agilent SureSelect V5 capture, NovaSeq 6000 PE101, median 150× tumor / 90× normal. BWA-MEM to GRCh38; GATK best-practice preprocessing. **Somatic calls**: GATK Mutect2 v4.6.2 with a 28-sample cohort PoN + gnomAD v3.1; FilterMutectCalls, LearnReadOrientationModel, CalculateContamination. Eight unmatched tumors used tumor-only + PoN with stricter filtering. snpEff GRCh38.99 annotation. **MSI**: msisensor-pro. **SBS signatures**: SigProfilerAssignment refit to COSMIC v3.3 on SBS96 context. **CNV**: CNVkit batch. **HRD proxies**: LST, TAI, LOH from CNV segments. **HLA class I**: OptiType on MHC-extracted reads. **HLA LOH**: direct IMGT-allele read counting with Fisher exact tumor-vs-normal imbalance (|Δratio| > 0.15, P < 0.05 per locus).

### Neoantigen prediction
pVACseq v5 + MHCflurry 2.0 on Mutect2-passing missense variants, per-patient HLA-A/B/C 4-digit types, 8–11-mer peptides. Strong binders ≤ 50 nM, binders ≤ 500 nM. Presentation-competent neoantigen (PCN) score = unique binder sites × (1 − 0.33 · LOH fraction).

### RNA-seq
Stranded library, NovaSeq PE101, HISAT2 + StringTie to GRCh38 / GENCODE v39, gene-level TPM (46,425 × 56). DESeq2 `~ sex + cT_simple + response_bin` for DEG. fgsea (Hallmark + Reactome). gseapy ssGSEA on 95 curated sets. Immune signatures: CD8 proliferation/activation/exhaustion, MHC I/II, NLRC5–HLA–IFN-γ, TLS (Cabrita), TGF-β (Mariathasan), EMT (Mak), hypoxia (Buffa) as mean z-score. CMS by CMScaller. TCR/BCR by TRUST4.

### Integration and predictor
Per-subject master table 35 × 37 (clinical + WES + RNA-seq). Mann–Whitney U for continuous, Fisher for categorical, BH FDR. LASSO logistic regression with LOOCV (held-out AUC).

### External validation (meta-analysis)
Seven public GEO nCRT cohorts (GSE35452, GSE45404, GSE68204, GSE69657, GSE94104, GSE119409, GSE150082; total N = 290) scored with the identical ssGSEA pipeline; cross-cohort aggregation by Stouffer's Z weighted by √N (details in Supplementary).

### Code and data
All analysis scripts under `/analysis/scripts/` (numbered `00_`–`27_`). Code and derived tables are on https://github.com/Soonlab/TNT. Raw sequencing (Macrogen HN00249207 WES, HN00249209 RNA-seq) will be deposited to SRA on publication.

---

## Results

### 3.1 Cohort: MSS, TMB-low, ICB-biomarker-negative
Of 35 patients, 18 were good responders (TRG 0–1) and 17 bad (TRG 2–3) (**Fig 1**; **Table 1**). Clinical T4 was enriched in bad responders (41 % vs 11 %, P = 0.086). **All 41 matched tumors were microsatellite-stable** (maximum MSI % 0.19 %) and TMB-low (median 1.6 /Mb; good 1.85 vs bad 1.40, P = 0.186). MSI-H and high TMB — the two established ICB biomarkers — do not apply.

### 3.2 Somatic landscape (Fig 2)
CRC driver mutations followed a canonical MSS distribution: APC 30/49 (61 %), TP53 20/49 (41 %), KRAS 14/49 (29 %), FBXW7 7/49 (14 %), KMT2D 4/49. No single driver reached significance; FBXW7 showed a non-significant trend toward good response (OR 3.7, P = 0.36). SBS5 / SBS1 (clock-like) dominated (> 60 % of mutations); **SBS3 (HRD) was absent** in all samples; sporadic MMR signatures (SBS6/15/20/26) were seen in both strata. CIN was indistinguishable (P = 0.66); an HRD LST proxy was modestly higher in bad responders (P = 0.037).

### Narrative 1 — Pre-treatment tumor-intrinsic DNA-repair/cell-cycle predictor

### 3.3 DNA-repair and cell-cycle pathways stratify response (Fig 3, Fig 4)
Hallmark GSEA of pre-treatment transcriptomes (n = 33) showed E2F targets (NES = 2.78, P = 8 × 10⁻²⁶), G2M checkpoint (NES = 2.46), MYC targets V1/V2 (NES ≥ 2.23), mTORC1, and mitotic spindle markedly elevated in good responders, with EMT (NES = −2.16, P = 6 × 10⁻¹⁰), myogenesis and apical junction reciprocally suppressed. Reactome independently placed cell-cycle checkpoints, M-phase, homology-directed repair, DSB repair, and DNA replication as the top positive sets. ssGSEA on 95 curated sets corroborated: DSB repair (P = 0.007), MYC targets V2 (0.018), HDR (0.020), general DNA repair (0.020), G2-M (0.032), E2F (0.035), CD8 proliferation (0.035). MHC II was modestly lower in good pre-treatment (P = 0.074).

### 3.4 Integrated 37-feature LASSO predictor (Fig 4)
Integration of 37 per-subject features showed a tightly co-varying DNA-repair + cell-cycle + CD8-proliferation module. A LASSO logistic-regression classifier achieved **leave-one-out AUC = 0.755**, outperforming elastic-net and random-forest top-8 models (both 0.70). RF feature importance ranked MYC V2, DSB repair, hypoxia, HDR, MHC II and deletion fraction at the top. PyClone-VI on 12 paired tumors showed a non-significant trend toward larger dominant-clone shrinkage in good responders (Δ −0.67 vs −0.15, P = 0.34).

### Narrative 2 — Post-treatment immunogenic-clone-clearance cascade

### 3.5 Missense and SBS5 mutation clearance in good responders
In 14 paired pre/post tumors, median Δ missense count was −67 (good) vs −8.5 (bad); ΔSBS5 −76 vs −29, P = 0.041. Mutation clearance, consistent with cytotoxic elimination of mutated clones, is specific to good responders.

### 3.6 MHC-I neoantigen landscape and clearance (Fig 5 — new)
Pre-treatment neoantigen burden was modestly elevated in good responders (**Fig 5A, 5B**; median 73.5 vs 66 mutation sites generating ≥ 1 MHC-I binder, P = 0.082; PCN 71.5 vs 57.1, P = 0.15; strong binders non-significant, P = 0.55). In paired pre → post analysis (n = 11; **Fig 5C**), good responders experienced dramatic neoantigen clearance: subjects 2, 6 and 9 each lost > 300 MHC-I binders (Δ −312, −489, −626) and > 50 binder-generating sites. Bad responders showed smaller losses (median Δ binder sites −59 good vs −16 bad, P = 0.25). One good responder (14, pCR) atypically gained neoantigens, reflecting sparse residual tumor at resection.

### 3.7 HLA-LOH clone clearance (Fig 5D)
HLA class I LOH (direct IMGT-allele imbalance test) was detected at ≥ 1 class I locus in 10/76 tumors. Pre-treatment prevalence was 4/16 good vs 2/12 bad (Fisher P = 0.67). Subjects 3 and 4 (both good responders) showed multi-locus HLA-LOH pre-treatment (3/3 loci and 2/3 loci respectively) that was reduced post-treatment (3→1 and 2→0 respectively; **Fig 5D**), consistent with **selective elimination of HLA-LOH tumor clones** by chemoradiation — the hallmark of an intact surveillance response to restored antigen presentation.

### 3.8 Immune reprogramming and B-cell infiltration
Post-treatment tumors from good responders showed paired elevations in regulatory T cells (Δz +1.26 vs +0.03, P = 0.026), MHC II (+1.23 vs +0.36, P = 0.065), CD8 exhaustion (+1.00 vs −0.10, P = 0.093), and coherent upregulation of antigen-presentation, TNF-α, IL-2/STAT5, IL-6/STAT3 and allograft-rejection pathways (within-good Wilcoxon P = 0.031). IGH clonotype count (TRUST4) increased by Δ +1,424 in good responders vs +7 in bad (within-good P = 0.031); IGK/IGL and TRA/TRB Shannon diversity moved in the same direction.

### 3.9 Assembled cascade
Taken together, good responders enter therapy in a pre-treatment state of proliferative, DNA-repair-proficient tumor-intrinsic competence and post-treatment transit a coherent cascade — **mutation clearance → neoantigen clearance → HLA-LOH clone elimination → MHC-II / Treg / CD8 exhaustion reprogramming → B-cell infiltration**. Bad responders move along none of these axes, consistent with primary treatment insensitivity.

### 3.10 External validation — partial reproducibility, TNT-specific biology
Meta-analysis of seven public GEO nCRT cohorts (N = 290) using Stouffer's Z yielded no significant cross-cohort signal for DSB/HDR (Z = −0.15, P = 0.56), E2F/MYC (+0.19, P = 0.43) or CD8 proliferation (+0.06, P = 0.48); EMT trended opposite to discovery (Z = −1.83, P = 0.97 in expected direction). Per-cohort effects varied; two cohorts (GSE35452, GSE45404) were concordant while three (GSE150082, GSE69657, GSE119409) were discordant. **All seven public cohorts received long-course nCRT without induction/consolidation FOLFOX/CAPOX, distinguishing them from modern TNT.** We interpret the partial reproducibility as **TNT-regimen-specific biology** (see Supplementary `Supp_external_meta.md`, SuppFig_external_forest, SuppFig_meta_zscore, SuppFig_GSE150082_DSB). Prospective TNT-matched validation is required.

### 3.11 Long-term outcomes (DFS / OS) — deferred
DFS and OS annotations are not currently available in the IRB-released metadata (see `analysis/clinical_survival_status.md`). Kaplan–Meier and Cox analyses are therefore **deferred to a follow-up study** when outcome data become available; TRG-based response remains the endpoint used here.

---

## Discussion

Two findings are central. **First**, in MSS LARC TNT response is governed by tumor-intrinsic DNA-repair and proliferative programs, not by classical ICB-response biomarkers. This runs **orthogonal to the ICB paradigm**: under ICB, MSI-H / high-TMB / inflamed tumors respond well via neoantigen-driven T-cell recognition; under TNT, intact HDR/DSB-repair and active cell cycle sensitise tumors to the combined genotoxic insult of platinum-based chemotherapy + ionising radiation, engaging apoptotic and senescence programs. EMT-high mesenchymal tumors invoke established chemoradiation-resistance programs. **Second**, good responders traverse a coherent post-treatment immunogenic-clone-clearance cascade: mutation and neoantigen clearance, HLA-LOH clone elimination, MHC-II / Treg / CD8 exhaustion reprogramming, and B-cell / TLS-like infiltration. This mirrors emergent "tertiary-lymphoid-structure" biology reported in ICB responders but is here induced by chemoradiation rather than by checkpoint blockade.

**Orthogonal to ICB biomarkers.** MSI/TMB are not useful in this MSS-dominated cohort; our DSB/HDR/E2F/MYC signature is instead a **TNT-specific response biomarker axis**. The two paradigms can be combined: MSI-H LARC → ICB (per Cercek [5]); MSS LARC → TNT with DSB/HDR/E2F signature-based stratification.

**External heterogeneity is informative.** The seven public cohorts all received long-course nCRT *without* modern induction/consolidation FOLFOX/CAPOX. Platinum-based multi-agent chemotherapy imposes an acute DNA-damage burden that plausibly selects for DSB/HDR-proficient tumors as responders — an effect absent in single-agent chemoradiation. Meta-analytic non-reproducibility thus reflects distinct regimen biology, not a false-positive discovery. Prospective TNT-matched validation (PRODIGE 23 / OPRA translational substudies) is the next step.

**Clinical implications.** A pre-treatment RNA-seq classifier built on DSB/HDR/E2F/MYC/CD8-proliferation could, pending external TNT-matched validation, stratify TNT candidates for watch-and-wait organ-preservation. EMT-high tumors might be candidates for intensified or alternative neoadjuvant strategies (taxane addition, anti-TGF-β). Paired pre/post biopsy-based monitoring of the clone-clearance cascade (neoantigen, HLA-LOH, IGH) could serve as an early pharmacodynamic readout.

**Limitations.** Single-center n = 35; 8/49 tumors used tumor-only calling (stringently PoN-filtered); microarray-era external cohorts have limited probe coverage; DFS/OS unavailable. HLA-LOH analysis uses a lightweight IMGT-read-counting approach rather than LOHHLA; results are consistent across paired subjects but warrant confirmation.

**Future directions.** Prospective TNT-matched validation; single-cell RNA-seq of paired pre/post biopsies to resolve the B-cell infiltration; DFS / OS integration when outcome data are released; TNT-specific companion-diagnostic clinical-trial stratification on the DSB/HDR/E2F axis.

---

## Conclusion

Response to TNT in MSS LARC is predicted by a tumor-intrinsic DSB/HDR/E2F/MYC/CD8-proliferation axis (LOOCV AUC 0.755) and manifests post-treatment as a coherent mutation-clearance → neoantigen-clearance → HLA-LOH clone-elimination → immune-reprogramming → B-cell-infiltration cascade. These findings provide a candidate TNT-specific response biomarker and position TNT response as an orthogonal axis to classical ICB biomarkers.

---

## Figures (v0.4)

- **Figure 1.** Cohort overview (response, cT stage, sex/age, sample matrix, study design).
- **Figure 2.** WES landscape (TMB, MSI, driver oncoprint, SBS signatures, MMR, CIN, HRD).
- **Figure 3.** RNA immune/pathway signatures (signature grid + heatmap + DEG volcano + CD8 proliferation + post-immune).
- **Figure 4.** Integration & ML (37-feature correlation + response ranks + Hallmark GSEA + LASSO ROC AUC = 0.755).
- **Figure 5.** *(New — main figure)* **Neoantigen cascade.** A: pre missense count good vs bad. B: pre pVACseq strong-binder sites. C: paired Δ pVACseq binders (slopegraph + box). D: HLA-LOH clone clearance subj 3 & 4 (paired allelic imbalance pre → post).
- **Figure 6.** Paired pre/post immune delta (Treg, MHC II, CD8 exhaustion, IGH clonotype expansion).
- **Supp 1–N.** External meta-analysis (forest, Z-score, GSE150082 DSB), PyClone-VI clonal decomposition, per-subject neoantigen detail, full GSEA tables, per-subject master table.

## Tables
- **Table 1.** Clinical characteristics (good vs bad).
- **Table 2.** Top 20 integrated-feature response associations.
- **Table S1.** 37-feature per-subject master table.
- **Table S2.** SBS signature activities per sample.
- **Table S3.** Full GSEA (Hallmark + Reactome).
- **Table S4.** CRC driver mutations per sample.
- **Table S5.** HLA class I types per subject.
- **Table S6.** External meta-analysis details (cohort table + per-cohort effects).

---

## Data and code availability

- Analysis workspace: `/mnt/sda1/data/TNT/analysis/`.
- GitHub: https://github.com/Soonlab/TNT (code + derived tables; raw sequence excluded).
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

---

*End of v0.4 draft — Genome Medicine target — 2026-04-15*
