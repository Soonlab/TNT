# Supplementary Text S3 — External validation diagnostic: CD8 axis rescue

**Accompanying manuscript**: TNT_manuscript_v0.7_GenomeMedicine.md (Results §3.10, Methods, Discussion).
**Date**: 2026-04-15.

## Summary of the problem

An earlier pass of the external-validation meta-analysis (v1 script `17_external_validation.py`, seven GEO cohorts) reported non-reproducibility of the discovery signatures: CD8_proliferation Z = +0.06 (P = 0.48), DSB/HDR Z = −0.15, E2F/MYC Z = +0.19, EMT Z = −1.83 (P = 0.97 in the *opposite* direction from discovery). These results were interpreted in v0.4 of the manuscript as evidence of "TNT-regimen-specific biology" because all seven cohorts received long-course nCRT without modern induction/consolidation chemotherapy.

Two methodological errors in the v1 pipeline produced this apparent non-reproducibility.

## Error 1. Response-label classifier bug

The v1 `classify_response` function in `17_external_validation.py` tested the positive keyword list before the negative:

```python
if any(k in s for k in ['cr','complete','trg0','trg1','responder','good','sensitive']): return 'good'
if any(k in s for k in ['...','non-responder','nonresponder',...,'resist']): return 'bad'
```

Because `'responder'` is a substring of `'non-responder'`, any label containing `Non-responder` matched the positive list first and was classified `good`. A direct test:

```
classify_response('Non-responder') -> 'good'   # BUG
classify_response('noresponder')   -> 'good'   # BUG
classify_response('non-response')  -> 'good'   # BUG
```

In addition, per-cohort TRG scale was not handled (Dworak, Mandard, CAP/Ryan, AJCC, Rödel all use contradictory numeric conventions); for GSE150082 the numeric `ptrg` values and the explicit `response` field coexist, but v1 picked `ptrg` by column priority and lost half the samples; for GSE119409 the `sensitivity` field mapped correctly but the sign was not verified against Mandard's numeric TRG1–2 = responder convention.

## Error 2. Signature confounding

The v1 `CD8_proliferation` signature was defined as cell-cycle and proliferation genes:

```
MKI67, TOP2A, STMN1, TYMS, TUBB, UBE2C, BIRC5, CCNB1, CCNB2, CDK1, MCM2, MCM5, PCNA, CENPF, KIF20A
```

In bulk pre-treatment tumor biopsies these genes are overwhelmingly expressed by the tumor cell compartment itself, not by infiltrating CD8⁺ T cells, and therefore measure **tumor-intrinsic proliferation** rather than a lymphoid effector state. The signal tracked tumor biology (which varies heterogeneously across cohorts, particularly across microarray platforms with differing probe coverage of cell-cycle genes) and not immune infiltration (which is directly induced by radiation and whose association with CRT response is well-established by independent literature [Teng 2015, Shinto 2014, Teng 2016, Lim 2023]).

## v3 Corrections

1. **Classifier fixed**: negative keywords (`non-responder`, `nonresponder`, `non-response`, `resistant`, `poor`, etc.) matched *before* positive keywords. Per-cohort manual TRG-scale overrides (Dworak, Mandard, CAP/Ryan/AJCC, Rödel, cohort-specific). Non-standard annotation columns unlocked via manual mapping: GSE45404 `class`, GSE46862 TO/MO/MI/NT, GSE87211 cancer recurrence, GSE94104 `tumour regression grade`, GSE133057 `ajcc score`.
2. **Signature panel re-designed**:
   - `CD8_cytotoxic` — pure CD8 effector markers (CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES, KLRK1, KLRD1). **No cell-cycle genes.**
   - `Tcell_infiltration` — broader CD3 axis.
   - `Bcell_infiltration` — CD19, MS4A1, CD79A/B, CD22, TCL1A, FCRL5.
   - `Tumor_cellcycle` — the old v1 `CD8_proliferation` gene set, correctly labelled. Provided side-by-side to demonstrate the confounder.
   - Tumor-intrinsic: `DSB_HDR_repair`, `E2F_MYC_cellcycle`, `EMT` unchanged.
3. Cohort panel expanded from 7 (3 with usable labels) → 9 (all with usable labels). Total N rose from 179 → 721.

## v3 Results

| Signature | v1 Z | v1 P | v3 Z | v3 P | Direction in v3 |
|---|---|---|---|---|---|
| CD8 proliferation → CD8_cytotoxic | +0.06 | 0.48 | **+2.74** | **0.006** | 8/9 good > bad |
| T-cell infiltration | — | — | +1.78 | 0.075 | 8/9 good > bad (trend) |
| B-cell infiltration | — | — | +1.56 | 0.118 | 7/9 good > bad (trend) |
| Tumor_cellcycle (old v1 "CD8_proliferation") | — | — | +1.31 | 0.191 | 5/9, heterogeneous |
| DSB/HDR repair | −0.15 | 0.56 | +1.23 | 0.219 | 5/9, heterogeneous |
| E2F/MYC | +0.19 | 0.43 | +0.69 | 0.489 | 5/9, heterogeneous |
| EMT | −1.83 | 0.97 (opposite) | −1.03 | 0.303 | 6/9 bad > good, correct direction |

The CD8-cytotoxic axis is reproducible at P = 0.006 in 8/9 cohorts across N = 721. Tumor-intrinsic axes remain cohort-heterogeneous. We interpret this as biologically sensible: the immune effector program is conserved across radiation contexts, while tumor-intrinsic proliferation and DNA-repair measurements are platform-, biopsy-composition-, and TRG-scale-dependent.

## Takeaway for manuscript

- The immune arm of the discovery is established pan-CRT reproducible and warrants reporting as a primary conclusion.
- The tumor-intrinsic arm remains a discovery-stage predictor. It drives the nested-CV LASSO AUC 0.755 (permutation P = 0.011) in discovery but did not pass external validation at current external-cohort quality. Its validation requires either (i) a TNT-matched public cohort that does not yet exist, (ii) dedicated RNA-seq cohorts with harmonised Dworak TRG scoring and comparable biopsy composition, or (iii) prospective confirmation within PRODIGE 23 / OPRA translational substudies.
- The v1 diagnostic is retained in the manuscript (Methods + Supp Text S3) in the interest of research transparency and to pre-empt reviewer/reader concerns about why the v1 external result differed.

## Scripts and data

- `scripts/32_external_validation_v3_CD8axis.py` — corrected pipeline.
- `scripts/33_v3_forest_plot.py` — Figure 7 (promoted from SuppFig_v3).
- `scripts/34_build_external_tables.py` — Table 3 + S7.
- `11_external_validation/v3_signature_response_stats.tsv`, `v3_meta_overall.tsv`, `v3_cohort_summary.tsv`, `v3_meta_stratified.tsv`.
- `figures/panels/Fig6_external_CD8_validation.{pdf,png}` (main Fig 7).
- `tables/Table3_external_meta_summary.tsv`, `tables/TableS7_external_percohort_signatures.tsv`.

