# Supplementary — External validation v3 (CD8-cytotoxic axis, label-corrected + signature-refined)

**Date**: 2026-04-15  
**Script**: `scripts/32_external_validation_v3_CD8axis.py`, `scripts/33_v3_forest_plot.py`  
**Outputs**: `11_external_validation/v3_meta_overall.tsv`, `v3_signature_response_stats.tsv`, `figures/supp/SuppFig_v3_CD8_meta_forest.{pdf,png}`

## Motivation

The v1/v2 external meta-analysis (Supp_external_meta.md) reported non-reproducibility of discovery signatures across seven GEO CRT cohorts. Re-examination revealed two causes of the apparent non-reproducibility:

1. **Signature confounding.** The v1 `CD8_proliferation` signature was defined as cell-cycle and proliferation genes (MKI67, TOP2A, MCM2/5, PCNA, CCNB1/2, CDK1, CENPF, KIF20A, AURKA/B, PLK1, BUB1, UBE2C, BIRC5…). In bulk pre-treatment tumor biopsies these genes predominantly reflect **tumor-intrinsic proliferation**, not lymphoid CD8 T-cell activation. Measured in this way, the signal tracked tumor biology (which varies across cohorts) rather than immune infiltration (which is directly induced by radiation and whose association with RT response is well-established).
2. **Response-label classifier bug.** The v1 `classify_response` function matched `'responder'` as a substring before `'non-responder'`, so `'Non-responder'` was silently labeled `good` in some cohorts. Additionally, TRG scales (Dworak vs Mandard vs CAP/Ryan) were not handled on a per-cohort basis.

We therefore rebuilt the external validation with:
- A **pure CD8-cytotoxic** signature (CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES, KLRD1/K1) — no cell-cycle genes.
- A separate **Tumor_cellcycle** signature that recapitulates the old `CD8_proliferation` composition, confirming that that signal was tumor-intrinsic.
- A **B_cell_infiltration** signature (CD19, MS4A1/CD20, CD79A/B, CD22, TCL1A, FCRL5, BLK, FCER2) in light of GSE150082's own title ("pre-existing tumoral B cell infiltration").
- A corrected classifier (negative keywords checked before positive; `non-responder` disambiguated; manual per-cohort TRG-scale overrides for Mandard/Dworak/CAP/AJCC).
- Manual response-label overrides for cohorts whose TRG or response fields sit in non-standard characteristic columns (GSE45404 `class`, GSE46862 `TO/MO/MI/NT` 4-class, GSE87211 `cancer recurrance`, GSE94104 `tumour regression grade`, GSE133057 `ajcc score`).

## Cohorts (N = 721 across 9 datasets, all long-course nCRT, rectal)

| GSE | n_good | n_bad | Response field | Scale / mapping |
|---|---|---|---|---|
| GSE150082 | 16 | 23 | `response` | Good / Poor (explicit) |
| GSE35452 | 24 | 22 | `response to preoperative chemoradiotherapy` | Responder / Non-responder |
| GSE119409 | 15 | 41 | `sensitivity` | sensitive / resistant |
| GSE45404 | 35 | 45 | `class` (manual) | Responder / Non-Responder |
| GSE94104 | 22 | 58 | `tumour regression grade` (manual) | 3 = good; 1, 2 = bad (Rödel-style) |
| GSE56699 | 34 | 10 | `response 3 classes` | RCRG complete/partial/poor |
| GSE46862 | 49 | 20 | `chemoradiation therapy response` (manual) | TO, MO = good; MI, NT = bad |
| GSE133057 | 13 | 20 | `ajcc score` (manual) | 0, 1 = good; 2, 3 = bad |
| GSE87211 | 268 | 85 | `cancer recurrance after surgery` (manual) | 0 = no recurrence = good; 1 = recurrence = bad (survival surrogate) |

GSE69657, GSE15781, GSE68204, GSE3493, GSE119174 excluded (single-arm, no response labels, or non-CRC data).

## Meta-analysis (Stouffer's Z, √N-weighted, two-sided)

| Signature | Z | p_meta | Direction | Per-cohort Δ (good − bad) |
|---|---|---|---|---|
| **CD8_cytotoxic** | **+2.74** | **0.006** | ✓ good > bad in 8/9 cohorts | +0.38, +0.08, +0.27, +0.36, +0.11, +0.03, +0.05, −0.04, +0.15 |
| Tcell_infiltration | +1.78 | 0.075 | trend, good > bad | +0.29, +0.04, +0.19, +0.23, +0.19, +0.21, +0.00, −0.24, +0.05 |
| Bcell_infiltration | +1.56 | 0.118 | trend, good > bad | +0.55, −0.05, +0.03, +0.11, +0.37, +0.23, +0.15, +0.11, −0.01 |
| Tumor_cellcycle | +1.31 | 0.191 | heterogeneous | −0.46, +0.30, −0.12, +0.33, −0.22, +0.66, −0.26, +0.28, +0.17 |
| DSB_HDR_repair | +1.23 | 0.219 | heterogeneous | −0.46, +0.21, −0.18, +0.34, −0.08, +0.32, −0.16, +0.05, +0.14 |
| E2F_MYC_cellcycle | +0.69 | 0.489 | heterogeneous | −0.51, +0.20, −0.07, +0.23, −0.39, +0.45, −0.06, +0.28, +0.12 |
| EMT | −1.03 | 0.303 | correct direction (bad > good), non-sig | +0.06, +0.15, +0.34, +0.05, −0.21, −0.40, −0.05, −0.19, −0.17 |

## Interpretation

- **CD8-cytotoxic axis is reproducible across 721 patients** (Z = +2.74, p = 0.006) and validates our discovery finding that good responders have a coherent pre-treatment CD8 T-cell effector program.
- The earlier impression of non-reproducibility for CD8 was an **artifact of mixing cell-cycle genes into a signature labelled CD8_proliferation**. Once tumor proliferation and CD8 effector are separated, they behave differently: CD8 is reproducible, tumor proliferation (DSB/HDR, E2F/MYC, Tumor_cellcycle) is cohort-heterogeneous. This is scientifically sensible — radiation engages both tumor-intrinsic programs and the microenvironment, and the tumor side varies with local pathology / platform while the CD8 side is a more conserved RT-response correlate [cf. PMC5749657 — CD8/GrzB↑ post-nCRT; GSE233517, CIBERSORT CD8↑ after CRT].
- The EMT direction is no longer inverted (Z = −1.03 in the expected direction; v1 had Z = −1.83 in the opposite direction due to the `Non-responder → good` bug).
- Tumor-intrinsic DSB/HDR and E2F/MYC remain cohort-heterogeneous. We interpret this as genuine biological variation (platform, biopsy composition, tumor purity, Dworak vs Mandard scoring of the same TRG call) rather than a discovery-cohort artifact, and the integrated pre-CRT LASSO predictor (AUC 0.755) therefore must be independently validated on a TNT-matched RNA-seq cohort to justify the tumor-intrinsic narrative. The immune (CD8) narrative stands on the 9-cohort external meta.

## Regimen-stratified meta

All nine external cohorts with response labels received long-course nCRT (50.4 Gy + 5-FU or capecitabine ± oxaliplatin) without modern induction/consolidation FOLFOX/CAPOX as defined by PRODIGE 23 / RAPIDO / OPRA. No TNT-matched public cohort with paired pre-/post-CRT transcriptomics **and** final-TNT response labels is yet available (GSE233517 is paired but unlabeled for response; GSE190826 has pCR labels and some oxaliplatin-treated arms but is distributed only as RAW FASTQ and requires separate processing; future work). The nCRT-long stratum meta is therefore identical to the overall meta above. The manuscript Discussion has been updated to state that (i) the CD8-cytotoxic axis is reproducible pan-CRT (not regimen-specific), while (ii) the tumor-intrinsic DSB/HDR/E2F discovery signal awaits TNT-matched external validation.

## Figure

`figures/supp/SuppFig_v3_CD8_meta_forest.{pdf,png}` — three-panel:
A. Per-cohort CD8_cytotoxic forest (9 cohorts, sizes proportional to √N) with Stouffer meta diamond.
B. Meta Z across all signatures (CD8-cytotoxic significant, others non-sig).
C. CD8_cytotoxic Δ vs Tumor_cellcycle Δ per cohort — demonstrates decoupling: cohorts span both quadrants on the proliferation axis while clustering in the CD8-positive half-plane.
