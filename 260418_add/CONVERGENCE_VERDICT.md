# 2026-04-18 Convergence test — VERDICT

## The question

User: "If DSB repair / HDR / E2F / G2M / Myc V2 are high at baseline, do those subjects
show a stronger post-treatment CD8 cascade?"

If yes → Discovery LASSO predictor (Thread 1) → drives external-validated CD8 axis (Thread 2)
        → manifested as paired cascade (Thread 3). One causal chain.

If no  → Three threads are three independent observations of the same patients. No causal chain.

## The test (script `09_targeted_convergence_test.py`)

36 pre-specified pairs:
- 9 baseline tumor-intrinsic features (the LASSO winners + key GSEA hits): DSB Repair,
  HDR, DNA Repair, E2F Targets, G2-M Checkpoint, Myc Targets V2, frac_amp, MHC_II, MSI_pct
- 4 cascade Δ features: CD8_cytotoxic_delta, Treg_delta, MHC_II_delta, IGH_n_delta
- n = 12 paired subjects (6 good + 6 bad), Spearman + response-adjusted partial Spearman
- BH FDR

## The result — `targeted_convergence_test.tsv`

| Test | Hits | Expected by chance |
|---|---|---|
| partial p < 0.05 | **0** | 1.8 |
| partial p < 0.01 | 0 | 0.4 |
| BH-q < 0.10 | 0 | — |

**The headline pair** — DSB Repair (baseline) → CD8 cytotoxic Δ — gives Spearman r = **−0.07** (p = 0.83). Essentially zero correlation.

All other "predictor → cascade" pairs are equally null:
- HDR baseline → CD8_cyt Δ: r = +0.13, p = 0.70
- E2F Targets baseline → CD8_cyt Δ: r = +0.13, p = 0.70
- G2-M baseline → CD8_cyt Δ: r = +0.13, p = 0.70
- Myc V2 baseline → CD8_cyt Δ: r = +0.03, p = 0.93
- DNA Repair baseline → CD8_cyt Δ: r = −0.01, p = 0.98
- MSI% baseline → CD8_cyt Δ: r = −0.37, p = 0.26
- frac_amp baseline → CD8_cyt Δ: r = +0.18, p = 0.59

The strongest signal in the entire 36-pair table is MHC_II baseline ↔ CD8_cyt Δ at r = −0.48, p = 0.11 (negative direction, still NS).

## Verdict

**Convergence hypothesis FAILS.** The discovery LASSO predictor and the paired post-CRT
cascade are **observationally independent** in this cohort. The "tumor-intrinsic state →
effective CRT → cell killing → antigen release → CD8 axis activation" causal chain is **not
supported** by the n = 12 paired data.

This is not a power problem in the strict sense: at n = 12 we can detect Spearman r ≥ 0.55
at p < 0.05; the observed r values are mostly in [-0.2, +0.2]. There is no signal to detect.

## What this means for the manuscript narrative

### Current v0.7.3 framing (becomes untenable)

§3.10 "Assembled cascade":
> mutation clearance → neoantigen clearance → HLA-LOH clone elimination → Treg/MHC-II/CD8
> exhaustion reprogramming → B-cell infiltration

This is presented as a single mechanistic chain. The convergence test shows there is no
quantitative subject-level link from baseline tumor-intrinsic state to the magnitude of
this cascade. The cascade may still happen (within-good Δ are robust for several features),
but it is **not predicted by the LASSO winners**.

### Two honest re-framings

**Option A — CD8-cytotoxic-centered paper (recommended).**
- Lead with external N = 721 + Akiyoshi N = 298 → CD8-cytotoxic axis as the primary,
  reproducible, mechanism-grounded finding (>1,000 patients across 10 cohorts).
- Discovery cohort: confirms direction of CD8 axis (good > bad medians for CD8_cyt and
  Tcell_infiltration in the paired Δ table) but is underpowered to reach univariate
  significance at N = 33 pre-treatment subjects.
- Tumor-intrinsic predictor (DSB / HDR / E2F / G2M / Myc V2) → demoted to "discovery-stage
  observation in this single cohort whose external transferability is unresolved" — moves
  from main Fig 5 to Supplementary.
- Paired cascade → reframed as "post-CRT immune phenotype consistent with effective
  treatment in good responders", **without** the causal arrow from baseline tumor state.
  Treg Δ (the only between-group robust paired finding) becomes the main paired result.

**Option B — Two parallel observations explicitly (less ambitious but defensible).**
- Discovery predictor (Thread 1) and external CD8 axis (Thread 2) are presented as two
  independent windows on CRT response biology, validated on different cohorts at different
  scales.
- Paired cascade (Thread 3) is described as "phenomenology of post-CRT immune
  reprogramming in good responders" without claiming it bridges 1 and 2.
- Conclusion: MSS LARC TNT response has at least two independent molecular axes —
  cohort-local tumor-intrinsic predictability + universal CD8-cytotoxic engagement — that
  do not collapse into a single causal model on the available paired data.

### Why Option A is recommended

1. The CD8-cytotoxic axis has the **strongest evidence** in the entire study (N > 1,000
   external + Akiyoshi convergent at OR 3.81, all on independent platforms).
2. The convergence test we just ran is a **clean negative**: even the directionally
   "favourable" pairs land at r ≈ 0 (e.g. DSB → CD8_cyt: r = −0.07). This is not a
   marginal "trend" we can soften — it's null.
3. The current Discussion already concedes that the tumor-intrinsic axis is cohort-
   heterogeneous externally (§3.11, §4). Removing the implicit "but it drives the cascade
   anyway" subtext aligns the discussion with the data.
4. A single-narrative paper centred on a 1,000-patient externally validated finding is
   a better fit for Genome Medicine than a three-narrative paper whose unifying causal
   model has been falsified.

## Files

| File | Purpose |
|---|---|
| `07_paired_immune_delta.py` | Score CD8_cyt/Tcell/Bcell/TLS Δ for the 12 paired subjects |
| `08_baseline_vs_cascade_corr.py` | Exhaustive 312-pair correlation (also null after FDR) |
| `09_targeted_convergence_test.py` | Pre-specified 36-pair targeted test (decisive null) |
| `paired_immune_delta_per_subject.tsv` | Per-subject paired Δ for new immune sigs |
| `baseline_vs_delta_corr_long.tsv` | Full 312-pair correlation table |
| `targeted_convergence_test.tsv` | The 36-pair targeted test result |
| `Fig_baseline_vs_cascade_heatmap.{pdf,png}` | All 312 partial-r heatmap |
| `Fig_targeted_convergence_heatmap.{pdf,png}` | 36-pair targeted heatmap (the key figure) |
| `Fig_targeted_DSB_vs_CD8cyt_scatter.{pdf,png}` | The 4-panel scatter showing the null |

## Decision needed from the user

1. Adopt **Option A** (CD8-axis-centered paper) → I will draft v0.7.4 with §3.10 demoted
   to a "phenomenology of post-CRT immune state" subsection and the cascade narrative
   replaced by a CD8-axis-led narrative spanning Threads 2 and 3.
2. Adopt **Option B** (two parallel observations) → I will draft v0.7.4 with explicit
   parallel framing and a more cautious Discussion.
3. Stay with v0.7.3 narrative + add this convergence test to Supplementary as a "we
   tested this and didn't find it" disclaimer → least disruption but the assembled cascade
   in §3.10 becomes harder to defend if a reviewer asks for exactly this analysis.
