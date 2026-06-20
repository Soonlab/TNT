# Convergence test repath-1 (2026-05-02)

Tests two suspicions about the original 36-pair convergence test (manuscript §3.9 / Supp Fig S20), which reported 0/36 P<0.05 after response-group-adjusted partial Spearman + BH/36, and was interpreted in the manuscript as supporting a "two-read clinical algorithm" (static and dynamic layers independent).

The clinically more implementable thesis being probed: **pre-treatment levels mechanistically predict radiation-phase dynamics**, justifying a single pre-treatment biopsy ("one-read algorithm"). Two suspected reasons the original test missed dependence:
- **Suspicion 1**: response-group adjustment kills the between-group dependence the thesis predicts (good: high baseline + strong cascade; bad: low both).
- **Suspicion 2**: 36-pair multiple testing burden dilutes mechanistic signal that 5 a-priori pairs would survive.

All four configurations run on the **same 12 RNA-paired subjects** (intersection of RNA-paired and legacy paired sets): subjects 1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14 (good=6 [2,4,6,8,9,14] / bad=6 [1,5,7,10,12,13]). Per-pair n recorded in TSV (some WES-derived Δ have <12 due to missing values).

Inputs (unchanged from manuscript pipeline):
- `260418_add/integrated_subject_master_v2.tsv` — baseline features
- `09_integration/paired_delta/paired_feature_long.tsv` — legacy paired Δ (SBS5, missense, neo_binders, neo_sites, MHC_II, IGH_n, TRB_shannon, CD8_exhaustion, Treg)
- `260418_add/paired_immune_delta_per_subject.tsv` — RNA Δ (CD8_cytotoxic, Tcell_infiltration, Bcell_infiltration, TLS_Cabrita)

## Reproduction status (Configuration ORIGINAL)

Headline pair `DSB Repair × Δ CD8-cytotoxic`:
- **Plain Spearman** r=**−0.070**, P=**0.829** (n=12) — matches manuscript-quoted "r=−0.07, P=0.83" exactly
- **Partial Spearman, response-adjusted** r=−0.168, P=0.601 (n=12) — matches prior `260418_add/targeted_convergence_test.tsv` exactly

Note: the manuscript text quotes the **plain** r/P (−0.07/0.83), while the test machinery is the partial-Spearman group-adjusted. The two are consistent (both null at this pair).

Whole-grid summary (partial Spearman, group-adjusted, BH/36):
- 0/36 partial_P<0.05
- 0/36 BH q<0.05
- 0/36 BH q<0.10

→ Manuscript reproduction confirmed.

## Configuration A — pooled correlation (no group adjustment), 36 pairs, BH/36

Plain Spearman | Plain Pearson:
- 0/36 spearman_P<0.05, **0/36 BH q_spearman<0.05**
- 1/36 pearson_P<0.05 (TGFb-Mariathasan-like null spike), 0/36 BH q_pearson<0.05

Top 5 strongest |spearman_r| (none survive BH):

| baseline | cascade | r | P | BH q |
|---|---|---|---|---|
| Myc Targets V2 | IGH_n_delta | +0.559 | 0.0586 | 0.807 |
| DNA Double-Strand Break Repair R-HSA-5693532 | IGH_n_delta | +0.545 | 0.0666 | 0.807 |
| DNA Repair R-HSA-73894 | IGH_n_delta | +0.538 | 0.0709 | 0.807 |
| frac_amp | Treg_delta | −0.495 | 0.1017 | 0.807 |
| MHC_II | CD8_cytotoxic_delta | −0.483 | 0.1121 | 0.807 |

→ Removing response-group adjustment does **not** rescue significance. The three near-trend pairs (Myc V2 / DSB / DNA Repair × IGH_n) match Supp Fig S20A's pre-existing "near-0.05" rows, none of which clear BH/36.

## Configuration B — 5 pre-specified mechanistic pairs, BH/5

| # | baseline (column used) | cascade Δ (column used) | n | plain r | plain P | plain q | partial r | partial P | partial q |
|---|---|---|---|---|---|---|---|---|---|
| 1 | DNA Double-Strand Break Repair R-HSA-5693532 | SBS5_delta          | 12 | **−0.595** | **0.0411** | 0.205 | −0.365 | 0.244 | 0.739 |
| 2 | DNA Double-Strand Break Repair R-HSA-5693532 | neo_binders_delta   |  9 | −0.417 | 0.265 | 0.441 | −0.249 | 0.519 | 0.739 |
| 3 | E2F Targets                                  | missense_delta      | 10 | −0.438 | 0.206 | 0.441 | −0.312 | 0.380 | 0.739 |
| 4 | CD8_cytotoxic                                | CD8_exhaustion_delta| 12 | −0.126 | 0.697 | 0.697 | −0.149 | 0.645 | 0.739 |
| 5 | E2F Targets                                  | Treg_delta          | 12 | +0.154 | 0.633 | 0.697 | −0.108 | 0.739 | 0.739 |

**Column-name mappings** (no derived data; existing columns used as-is):
- "DSB/HDR repair" → `DNA Double-Strand Break Repair R-HSA-5693532` (Reactome composite, Thread-1 strongest single feature: top LASSO coef +0.89, external Z=+3.17)
- "Tumor cell cycle" → `E2F Targets` (Thread-1 cellcycle representative; G2-M Checkpoint and Myc Targets V2 give similar r at trend level — explored in Config A)
- "Δ TMB" → `missense_delta` (per-tumor missense count = TMB proxy)
- "Δ MHC-I neoantigen binders" → `neo_binders_delta`
- "Δ CD8 exhaustion" → `CD8_exhaustion_delta`
- "CD8-cytotoxic" baseline → `CD8_cytotoxic` (RNA score)

Verdict for B: 1/5 plain P<0.05 (DSB × ΔSBS5, P=0.041), but **0/5 BH q<0.05** under either method. Three pairs (rows 1–3) show negative r consistent with the "high baseline → more radiation-induced clearance" direction; magnitude ranges |r|=0.42–0.60 but only the DSB × ΔSBS5 pair reaches nominal significance. Group-adjusted partials are uniformly weaker (|r| down by ~0.1–0.2), confirming part of the signal *is* between-group and dies under partial correlation — but not enough to clear even the BH/5 threshold.

## Configuration C — between-group direction concordance, sign test (5 pairs)

| # | baseline | baseline good | baseline bad | baseline dir | cascade good | cascade bad | cascade dir | concordant |
|---|---|---|---|---|---|---|---|---|
| 1 | DSB Repair | +0.827 | +0.804 | good>bad | −29.3 | −81.3 | bad>good | **NO** |
| 2 | DSB Repair | +0.827 | +0.804 | good>bad | −346.2 | −164.2 | bad>good | **NO** |
| 3 | E2F Targets | +1.025 | +0.989 | good>bad | −69.0 | −12.7 | bad>good | **NO** |
| 4 | CD8_cytotoxic | −0.273 | −0.184 | bad>good | +0.97 | +0.17 | good>bad | **NO** |
| 5 | E2F Targets | +1.025 | +0.989 | good>bad | +1.16 | +0.18 | good>bad | YES |

**Direction concordant: 1/5; one-sided sign test (binomial p=0.5, H1: p>0.5) P=0.969** — strongly anti-concordant.

The pattern is interpretable: for the three "mutation clearance" pairs (1–3), bad responders show **more negative cascade Δ** (i.e., more "clearance" in absolute terms) — not because their clearance machinery is better, but because their baseline mutation burden is larger; absolute Δ ∝ baseline magnitude, not response quality. This is exactly the "regression-to-the-mean" / scaling artifact a between-group concordance test is designed to flush out.

## Verdict

Repath-1 is **null on every test that matters**. Removing the response-group adjustment (Config A) does not rescue any pair — the 3 near-trend correlations remain at P=0.06–0.07 and lose all their already-thin BH cushion. Pre-specifying 5 mechanistic pairs (Config B) yields 1/5 nominal P<0.05 (DSB × ΔSBS5, plain), but 0/5 after BH/5 correction; group-adjusted partials are weaker still. The direct between-group test (Config C) actively contradicts the thesis: 1/5 concordant, sign-test P=0.97, with the three mutation-clearance pairs flipping in the opposite direction because absolute Δ scales with baseline magnitude rather than response. The data do **not** support pre-treatment level → radiation-phase dynamics; the static and dynamic layers are observationally independent on this cohort, and the manuscript's two-read framing remains the supportable position. The "one-read algorithm" thesis would need a substantially larger or differently-stratified cohort, or a reframing in terms of relative (not absolute) cascade magnitude.

## Files

```
convergence_repath1_260502/
├── README.md
├── scripts/01_repath1_analysis.py
├── tables/
│   ├── config_original_36pair.tsv
│   ├── config_A_pooled_36pair.tsv
│   ├── config_B_mechanistic_5pair.tsv
│   ├── config_C_directionconc_5pair.tsv
│   └── _summary.json
├── figures/
│   ├── FigA_pooled_heatmap.{pdf,png}
│   └── FigB_mechanistic_scatter.{pdf,png}
└── logs/01_run.log
```
