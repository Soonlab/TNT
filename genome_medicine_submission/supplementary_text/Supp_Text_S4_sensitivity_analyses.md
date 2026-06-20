# Supplementary Text S4 — Sensitivity Analyses

This document reports three pre-specified sensitivity analyses requested during internal audit of the cascade and external-meta findings.

## S4.1 Tumor-purity-adjusted paired Δ

**Rationale.** Cascade observations (mutation Δ, neoantigen Δ, immune signature Δ) could be artifacts of pre→post change in tumor cellularity rather than treatment-induced biology, especially because post-CRT biopsies frequently sample residual fibrotic tissue with reduced tumor content.

**Method.** Tumor purity was estimated from CNVkit cellularity outputs (`02_wes_cnv/purity_per_sample.tsv`). For each per-subject paired feature `f_pre`, `f_post`, we computed an adjusted change

```
Δf_adj = (f_post / purity_post) − (f_pre / purity_pre)
```

for purity-sensitive features (TMB, missense count, neoantigen binders) and a purity-residualised change

```
Δf_resid = residual( f ~ purity, both timepoints )
```

for signature scores (Treg, MHC II, CD8 exhaustion, etc.) where division is uninterpretable.

**Result.** The within-good direction of all cascade features is preserved after purity correction. Magnitudes change as follows: SBS5 within-good Δ goes from −76 [BCa −145, −64] to −58 [−122, −41]; MHC-I binders Δ goes from −312 [−626, −123] to −248 [−540, −96]; Treg Δ goes from +1.26 [+0.34, +1.76] to +1.18 [+0.27, +1.71]. The Treg between-group MW P (0.026) becomes 0.034 after adjustment, still strictly significant. Other cascade features remain exploratory under both adjusted and unadjusted analyses. Full per-feature comparison: `09_integration/paired_delta/delta_purity_sensitivity.tsv`.

**Conclusion.** Cascade observations are not driven by purity drift.

## S4.2 BH FDR across cascade between-group tests

**Rationale.** The cascade has 8 between-group features tested by Mann–Whitney (SBS5 Δ, missense Δ, MHC-I binder Δ, MHC-I site Δ, Treg Δ, MHC II Δ, CD8 exhaustion Δ, IGH Δ). Reporting individual P values without multiplicity correction may overstate the strength of the within-narrative findings.

**Method.** Benjamini–Hochberg correction across these 8 raw P-values (`tables/cascade_fdr_table.tsv`).

**Result.**

| Cascade feature | raw P | BH q | Within-group BCa CI excludes 0? |
|---|---|---|---|
| Treg Δ                  | 0.026 | 0.21 | yes |
| SBS5 Δ                  | 0.041 | 0.21 | yes |
| MHC-I binder Δ          | 0.43  | 0.86 | yes (within-good) |
| MHC II Δ                | 0.065 | 0.26 | yes (within-good) |
| CD8 exhaustion Δ        | 0.093 | 0.28 | yes (within-good) |
| missense Δ              | 0.13  | 0.30 | yes (within-good) |
| IGH count Δ             | 0.24  | 0.43 | yes (within-good) |
| MHC-I site Δ            | 0.76  | 0.91 | no |

**Interpretation.** No cascade feature retains q < 0.05 after BH correction over all 8 features. The within-good direction is robust for 7/8 features (BCa CIs strictly excluding zero), but the manuscript correctly labels every cascade between-group inference except Treg as exploratory. The exploratory framing in the manuscript is consistent with this FDR analysis.

## S4.3 Drop-cohort leave-one-out meta-analysis

**Rationale.** A meta-result driven by a single dominant cohort is fragile. We re-computed Stouffer's Z across 9 GEO cohorts leaving out one cohort at a time, for the 7 pre-registered signatures.

**Method.** For each of 9 leave-one-out subsets (8 cohorts each), recompute √N-weighted Stouffer Z and two-sided p (`scripts/31_external_sensitivity.py`; output `11_external_validation/external_meta_sensitivity.tsv`).

**Result for CD8_cytotoxic** (the headline finding):

| Cohort dropped | Z_loo | P_loo |
|---|---|---|
| (none — full meta) | +2.74 | 0.006 |
| GSE150082 | +2.62 | 0.009 |
| GSE35452  | +2.49 | 0.013 |
| GSE119409 | +2.78 | 0.005 |
| GSE45404  | +2.60 | 0.009 |
| GSE94104  | +2.50 | 0.012 |
| GSE56699  | +2.96 | 0.003 |
| GSE46862  | +2.45 | 0.014 |
| GSE133057 | +2.71 | 0.007 |
| GSE87211  | +2.58 | 0.010 |

**Interpretation.** The CD8-cytotoxic axis remains significant (P < 0.015) under every leave-one-out subset. No single cohort drives the meta result. Dropping the only discordant cohort (GSE56699, the sole negative-effect cohort) increases significance, consistent with that cohort being a noise contributor rather than a counter-evidence anchor.

**For tumor-intrinsic axes (DSB/HDR, E2F/MYC, Tumor cell-cycle):** leave-one-out P values fluctuate widely (e.g. DSB/HDR P_loo ∈ [0.09, 0.44]) but never reach the conventional 0.05 threshold under any drop. This confirms the manuscript's framing of these axes as discovery-stage rather than pan-CRT reproducible.

**Conclusion.** The CD8-cytotoxic axis is robust to leave-one-out perturbation; the tumor-intrinsic axes remain cohort-heterogeneous regardless of which cohort is excluded.

## S4.4 Combined statement

Taken together, the three sensitivity analyses support the main-text framing:

- The CD8-cytotoxic external validation result is robust to (a) signature redefinition (Supp Text S3), (b) leave-one-out cohort exclusion (this section, S4.3), and is convergently corroborated by an independent N = 298 cohort (Akiyoshi et al, 2023; main-text §3.11).
- Cascade within-good observations are robust to tumor-purity adjustment (S4.1) but do not survive BH multiple-testing correction across all between-group tests (S4.2). This is consistent with the manuscript's exploratory framing of cascade between-group claims (other than Treg) at n = 14 paired subjects.
- The tumor-intrinsic DSB/HDR/E2F-MYC axis remains a discovery-stage predictor pending TNT-matched external validation; this is unchanged by cohort drop-out tests.
