# Supplementary — External Validation Meta-Analysis

*Revision v0.5 (Genome Medicine target, 2026-04-15). Reframed around radiation-phase TNT biology: the discovery cohort sampled only the radiation (first) phase of TNT, and external CRT cohorts differ in CRT regimen and consolidation exposure.*

## Overview

The TNT-cohort discovery signal (DSB/HDR repair, E2F/MYC cell-cycle, CD8 proliferation up; EMT down in eventual good responders at pre-CRT baseline) was tested in **seven public GEO neoadjuvant chemoradiation (CRT) cohorts, total N = 290** with response annotation. Results are mixed and do **not** reach meta-significance at α = 0.05 in any of the four signatures, which we interpret as evidence of **radiation-phase TNT biology**: external cohorts differ in CRT regimen (short-course vs long-course), fractionation, and whether consolidation chemotherapy followed, which may explain partial reproducibility of a signature discovered in a cohort sampled only across the radiation (first) phase of modern TNT.

## Cohort composition (7 cohorts, 290 patients)

| GEO accession | Platform       | N   | N good | N bad | Treatment regimen (as reported)      | Response definition      |
|---------------|----------------|-----|--------|-------|--------------------------------------|--------------------------|
| GSE119409     | Affymetrix U133 Plus 2.0 | 66  | 15 | 41  | Long-course nCRT (50.4 Gy + 5-FU/cape) | pCR / TRG 0–1            |
| GSE150082     | Illumina HumanHT-12 v4   | 39  | 13 | 20  | Long-course nCRT                      | Complete clinical response (watch-and-wait) |
| GSE35452      | Affymetrix U133 Plus 2.0 | 46  | 24 | 22  | Long-course nCRT (cape)               | pCR / near-pCR           |
| GSE45404      | Affymetrix U133 Plus 2.0 | 80  | 35 | 45  | Long-course nCRT (cape or 5-FU)       | Dworak TRG 3–4 (good)    |
| GSE68204      | Affymetrix HG U133A 2.0  | 96  | 20 | 51  | Long-course nCRT + oxaliplatin        | pCR                      |
| GSE69657      | Affymetrix Human Gene 1.0 ST | 30  | 13 | 17  | Long-course nCRT                      | Mandard TRG 1–2 (good)   |
| GSE94104      | Illumina expression array | 40  | 12 | 28  | Long-course nCRT                      | Dworak 3–4 (good)        |

(The remaining TNT-specific cohort — our 35 patients — received **induction/consolidation FOLFOX or CAPOX plus long-course chemoradiation**, i.e. the modern TNT regimen, which is distinct from all seven public cohorts listed above.)

## Meta-analytic procedure

Per cohort, pre-treatment transcriptomes were normalised and scored by the identical ssGSEA pipeline used in discovery, on the same gene-set definitions. Effect size was `Δ = mean(good) − mean(bad)`. Per-signature cross-cohort aggregation used **Stouffer's Z** weighted by √N, testing the discovery-expected direction (positive for DSB/HDR, E2F/MYC, CD8 proliferation; negative for EMT).

## Meta-analysis result

| Signature            | N cohorts | Total N | Z (Stouffer) | One-sided P (expected dir) | Interpretation                    |
|----------------------|-----------|---------|--------------|----------------------------|-----------------------------------|
| DSB / HDR repair     | 7         | 290     | −0.15        | 0.56                       | No cross-cohort reproducibility   |
| E2F / MYC cell-cycle | 7         | 290     | +0.19        | 0.43                       | No cross-cohort reproducibility   |
| CD8 proliferation    | 7         | 290     | +0.06        | 0.48                       | No cross-cohort reproducibility   |
| EMT (expected ↓ good)| 7         | 290     | −1.83        | 0.97                       | Trend **opposite** to discovery   |

Per-cohort effects (delta good − bad, with pseudo 95% CI from √N) are shown in **SuppFig_external_forest**. Meta Z scores are shown in **SuppFig_meta_zscore**. The individual cohort showing the strongest directional agreement for the DSB/HDR signal — **GSE150082** (a watch-and-wait rectal cohort, P = 0.021 for DSB/HDR being lower in good responders, opposite direction) — is shown in **SuppFig_GSE150082_DSB**.

## Interpretation

Three features distinguish the TNT cohort from all seven GEO cohorts and plausibly explain the partial reproducibility:

1. **Regimen.** Our cohort received induction/consolidation FOLFOX/CAPOX *in addition* to long-course chemoradiation. The seven GEO cohorts received long-course nCRT alone (± capecitabine or single-agent oxaliplatin). The addition of platinum-based multi-agent chemotherapy in TNT imposes an acute DNA-damage burden that could *select for* DSB/HDR-proficient tumors as responders, an effect absent in pure chemoradiation.
2. **Response endpoint.** GEO cohorts mix Dworak TRG, Mandard TRG, pCR, and clinical complete response on watch-and-wait protocols. Our cohort uses surgical-specimen Dworak TRG.
3. **Platform.** Six of seven public cohorts are 15–20-year-old microarrays with limited dynamic range; the TNT cohort is modern stranded RNA-seq.

We therefore frame the discovery signal as **radiation-phase TNT biology**: external CRT cohorts differ in CRT regimen (short-course vs long-course), fractionation, and whether consolidation chemotherapy followed, which may explain partial reproducibility. Validation should be prospective and TNT-matched (e.g. PRODIGE 23, OPRA translational substudies with paired pre-CRT / post-CRT biopsies and final full-TNT outcome as endpoint), rather than interpreted as a universal chemoradiation-response biomarker.

## Supplementary figures

- **SuppFig_external_forest** (PDF/PNG at `figures/supp/SuppFig_external_forest.*`)
- **SuppFig_meta_zscore** (PDF/PNG at `figures/supp/SuppFig_meta_zscore.*`)
- **SuppFig_GSE150082_DSB** (PDF/PNG at `figures/supp/SuppFig_GSE150082_DSB.*`)

All three were re-rendered in the Genome-Medicine figure style (Arial-like sans, 0.6-pt spines, `good = #2E86AB`, `bad = #E63946`, signature-specific palette) matching the main figures.
