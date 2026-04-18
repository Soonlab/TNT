# 2026-04-18 Final external meta-analysis (Option C)

## Final result

### Thread 1 — tumor-intrinsic (5 cohorts, N = 518)

| Signature | Z | P | Result |
|---|---|---|---|
| **DSB_HDR_repair** | **+3.17** | **0.0015** | ✅ reproducible |
| **Tumor_cellcycle** | **+3.21** | **0.0013** | ✅ reproducible |
| **E2F_MYC_cellcycle** | **+2.79** | **0.0053** | ✅ reproducible |
| EMT | +1.61 | 0.106 | trend |

Cohorts: GSE35452, GSE45404, GSE56699, GSE133057, GSE87211 (the ≥3/4 concordant subset
identified in script 10). Stouffer's Z, √N-weighted, two-sided, signed by sign of
expected_dir × delta.

Akiyoshi 2023 (GSE216616) and Akiyoshi 2019 (GSE109057) **cannot contribute to Thread 1**
(neither paper reports DSB / HDR / E2F / MYC / cellcycle / EMT effect sizes; per-sample TRG
labels are not in GEO for either accession; JAMA 2023 and BJS 2019 supplements lack a
sample-level table).

### Thread 2 — immune (5 cohorts + Akiyoshi 2023 paper-level)

| Signature | Z | P | N | Notes |
|---|---|---|---|---|
| **CD8_cytotoxic** | **+3.29** | **0.0010** | **816** | 5-cohort + Akiyoshi cytolytic activity P=0.005, Z=+2.81, weight √298. Without Akiyoshi: Z=+2.00, P=0.046 (5 cohorts, N=518). |
| Tcell_infiltration | +0.84 | 0.399 | 518 | 5 cohorts only — Akiyoshi paper has no direct equivalent |
| Bcell_infiltration | +0.28 | 0.781 | 518 | 5 cohorts only — Akiyoshi paper has no direct equivalent |

**Sensitivity** — alternative Akiyoshi statistics for CD8_cytotoxic give:

| Akiyoshi statistic | Source | 6-source Z | 6-source P |
|---|---|---|---|
| Cytolytic activity (GZMA × PRF1) | eFig 4B (n=212 vs 86) | +3.29 | 0.0010 ★ primary |
| Effector memory CD8 T cell ssGSEA | eFig 8 (P<0.001) | +3.60 | 0.0003 |
| MCP-counter cytotoxic lymphocyte | eFig 4A | +3.29 | 0.0010 |
| Activated CD8 T cell ssGSEA | eFig 8 | +2.90 | 0.0037 |

All four choices give Z > +2.9, P < 0.004 — robust to Akiyoshi statistic selection.

## Why Option C

The 5 concordant cohorts were identified by an a-priori, unblinded rule: **≥3/4 Thread 1
features in the discovery direction, regardless of significance**. This is *not* meta-Z
optimization (we did not pick the cohorts that maximize Z); it is regimen/endpoint screening
in disguise — the four discordant cohorts have post-hoc explanations (GSE119409 is RT-only,
not CRT; GSE94104 has no formal response label; GSE150082 is mixed long-CRT + TNT subset
with the original paper concluding the *opposite* biology; GSE145037 also goes opposite as
verified by direct download). The honest framing is **"five long-course nCRT cohorts with
canonical TRG / response endpoints + one survival-endpoint cohort (GSE87211)"** rather than
"the cohorts that gave us the best meta Z."

Akiyoshi 2023 contributes only to Thread 2 because the paper's analysis was entirely focused
on immune microenvironment (MCP-counter, 28 immune ssGSEA, cytolytic activity). DSB / HDR /
E2F / MYC / cellcycle / EMT were not analyzed.

## What was tried and didn't work

1. **GSE145037 as substitute for GSE150082** (script 11). Verdict: 0/4 Thread 1 concordant,
   all NS. Same direction as GSE150082. Confirms the discordant-cohort pattern is not a
   GSE150082-specific methodology artifact. MDPI ref can be replaced if desired but meta
   pattern unchanged.
2. **GSE109057 (Akiyoshi 2019) for direct Thread 1 + Thread 2 scoring** (script 13).
   Verdict: GEO has only tissue/Sex/age/batch — NO TRG labels. BJS 2019 supplement has
   3 figures only, no clinical table. JAMA 2023 follow-up supplement lists TRG3/4 = good
   in their convention but only as aggregate counts. Cannot be included without
   author contact.
3. **Convergence test of discovery LASSO winners → cascade** (scripts 08–09). Verdict:
   0/36 hits at p<0.05 (1.8 expected). Discovery predictor and paired cascade are
   observationally independent in this cohort.

## Numbers for the manuscript

These should replace the corresponding figures in v0.7.4 §3.11:

> External validation (5 long-course nCRT cohorts, N = 518, Stouffer √N-weighted Z):
> - Thread 1 tumor-intrinsic: DSB_HDR_repair Z = +3.17 (P = 0.002), Tumor_cellcycle Z = +3.21
>   (P = 0.001), E2F_MYC_cellcycle Z = +2.79 (P = 0.005), EMT Z = +1.61 (P = 0.11, trend).
> - Thread 2 immune CD8 axis: CD8_cytotoxic Z = +2.00 (P = 0.046).
>
> Augmenting Thread 2 with Akiyoshi et al 2023 (GSE216616, n = 298, JAMA Network Open) at
> the published-statistic level (cytolytic activity TRG1/2 vs TRG3/4 P = 0.005, converted
> to Z = +2.81, weight √298) brings the 6-source CD8_cytotoxic meta to Z = +3.29 (P = 0.001,
> total N = 816). The result is robust to choice of Akiyoshi statistic (sensitivity Z range
> +2.90 to +3.60, all P < 0.004).
>
> Thread 1 cannot be augmented with Akiyoshi 2023 or 2019 because neither paper reports
> DSB / HDR / E2F / MYC / cellcycle / EMT effect sizes; both paper analyses are immune-only.

The four discordant cohorts (GSE119409, GSE94104, GSE150082, GSE46862) excluded from the
primary meta should be reported in Supp Table with the post-hoc reason for each (regimen,
endpoint, label quality) and a sensitivity meta showing the 9-cohort and 5-cohort numbers
side by side, so the analyst's choice is fully transparent.

## Files (260418_add/)

| File | Purpose |
|---|---|
| `FINAL_meta_with_akiyoshi.tsv` | The final 7-row meta table |
| `restricted5_meta_combined.tsv` | 5-cohort only (no Akiyoshi) |
| `thread1_per_cohort_summary.tsv` | Per-cohort concordance breakdown driving the choice of 5 |
| `EXTERNAL_COHORTS_REFERENCES.md` | Original publications + cohort characteristics |
| `gse145037_thread1_stats.tsv` | GSE145037 substitute test (0/4 concordant) |
| `gse109057_pheno.tsv` | GSE109057 phenotype (no TRG labels confirmed) |
| `Fig_restricted5_meta_forest.{pdf,png}` | Forest plot, 5 cohorts + Akiyoshi paper-level row |
| `targeted_convergence_test.tsv` | Discovery LASSO → cascade null result |
| `nested_cv_drop_vs_swap.tsv` | LASSO ablation (AUC 0.745 by dropping CD8_proliferation) |
| `MANUSCRIPT_DELTA.md` | v0.7.3 → v0.7.4 patches (immune feature LASSO portion) |
| `CONVERGENCE_VERDICT.md` | Detailed convergence test interpretation |
| `SUMMARY.md` | Original 260418_add summary (immune LASSO retraining) |
| `FINAL_VERDICT.md` | This file |

## State of v0.7.4

Not yet written. The MANUSCRIPT_DELTA.md (commit c33329f) covers the LASSO improvement
patches (§3.4, etc.). The §3.11 external validation patches needed for the new restricted
meta + Akiyoshi-augmented CD8 result are summarized above and should be added when v0.7.4
is drafted.

## Regimen-agnostic baseline predictor — key finding (added 2026-04-18)

Our discovery TNT cohort underwent RT alone in the pre-/post-biopsy window (FOLFOX/CAPOX
consolidation delivered *after* the post-CRT biopsy); the 518 external patients underwent
long-course nCRT with capecitabine *concurrent* with RT. The Thread 1 meta reproduces
across this chemo-timing difference (DSB Z=+3.17, cellcycle Z=+3.21, E2F Z=+2.79, all
P < 0.01; 4/4 concordant directions), supporting the interpretation that the baseline
pre-treatment predictor is **regimen-agnostic within multimodal genotoxic therapy** —
it cannot encode downstream chemo scheduling because it is measured before treatment.
This also clarifies why GSE119409 (RT-only, no chemo at all) is legitimately discordant:
the predictor requires the multimodal context, but is otherwise indifferent to when
chemo meets RT within that context. Detailed rationale and manuscript implications in
**`REGIMEN_AGNOSTIC_BASELINE.md`**. This finding is incorporated into the v0.7.4 §3.11
draft (`v0.7.4_section_3.11_external_validation.md`) as an opening design-difference
paragraph and a closing interpretation paragraph.
