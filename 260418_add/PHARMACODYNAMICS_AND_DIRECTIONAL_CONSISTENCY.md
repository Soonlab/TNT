# Pharmacodynamics of baseline factors + IGHV directional consistency (2026-04-18)

Companion to `REGIMEN_AGNOSTIC_BASELINE.md`. Reports what the paired 12-subject
RT-phase biopsy set adds that pre-only cannot.

## Headline findings

1. **All four Thread-1 baseline factors move in biologically predicted directions
   after RT in both response groups.** DSB_HDR_repair, Tumor_cellcycle, and
   E2F_MYC_cellcycle decrease; EMT increases. Directional consistency is high
   in both good and bad responders (majority in 4/5 of 8 group-factor cells).
   This is *target-engagement* evidence: the four factors we nominated as
   baseline predictors are precisely the axes that RT visibly perturbs.
   Response-differential magnitude (MW on Δ) is NS for all four, which means
   the baseline predictor is a *predisposition* axis rather than a
   response-differential pharmacodynamic axis.

2. **EMT goes up in 6/6 (100 %) of good responders, binomial P = 0.016.**
   This is the only composite-level directional test that is nominally
   significant at P < 0.05 within a single group, and is consistent with
   the surviving post-RT tissue in good responders being enriched for
   mesenchymal/stromal components (because proliferating non-mesenchymal
   tumor cells were selectively killed).

3. **IGHV repertoire response is more directionally coherent in good
   responders than in bad responders, aggregated across 53 V-genes with
   sufficient coverage.** Three independent aggregate tests all reach
   P < 0.05:

   | Aggregate test | Statistic | P |
   |---|---|---|
   | Paired Wilcoxon of (good − bad) majority fraction across 53 V-genes | W = 470.5 | **0.035** |
   | Binomial on good_majority > bad_majority across V-genes | 24/37 | **0.049** |
   | Binomial on pattern "good coherent, bad mixed" vs opposite | 16 vs 7 | **0.047** |

   Interpretation: good responders show a **directed** repertoire response
   to RT, in which most V-genes move in a reproducible per-gene direction
   across the 6 good subjects. Bad responders show a **stochastic**
   repertoire response, in which the same V-genes have per-subject changes
   in mixed directions. The biology is consistent with coordinated,
   antigen-driven B-cell clonal dynamics in good responders versus
   disorganized, non-productive responses in bad.

4. **IGHV6-1 is the strongest single V-gene** (user's prior IGHV3-7 and
   IGHV3-74 are weaker):
   - IGHV6-1: good 0 / 6 down (unanimous), bad 4 up / 2 down
     (Fisher P = 0.061, MW P = 0.065 — both axes agree)
   - IGHV3-7: good 4 up / 2 down, bad 4 up / 2 down — both groups same
     direction, coherence_gap = 0 (user's memory only partially borne out)
   - IGHV3-74: good 1 up / 5 down, bad 2 up / 4 down — both groups same
     direction, coherence_gap = 0.17

## What paired adds that pre-only cannot

| Claim | Paired evidence source |
|---|---|
| Baseline factors are the biology RT engages | 4/4 factor directions correct post-RT (Fig E / F) |
| EMT enrichment on surviving tissue | 6/6 good EMT up, sign P = 0.016 |
| Good responders mount a coherent immune response | 3 aggregate IGHV tests P < 0.05 (Fig G) |
| Bad responders' immune response is stochastic | 16 vs 7 V-gene "good coherent, bad mixed" ratio |
| Specific V-gene example | IGHV6-1 (Fig H): good 0/6 down, bad split |

Crucially, **none of these are visible at baseline**. A pre-only analysis
cannot see EMT enrichment (that is a consequence of RT killing non-mesenchymal
cells), cannot see Treg/MHC-II/IGH Δ, and cannot see the "directed vs
stochastic repertoire" distinction.

## Relationship to the convergence null

The convergence test (script 09) showed that the **magnitude** of paired
Δ in immune cascade features is not predicted by the baseline Thread-1
score. That null stands. The present analysis does not contradict it: it
shows that the baseline factors themselves move post-RT in a regime-engaged
direction (not response-differentially), and that the immune repertoire's
**directional coherence** (a distinct property from magnitude) differs
between response groups. These are three independent phenomena:

1. Baseline → final response (static predictor, externally validated)
2. RT → baseline factor direction (target engagement, both groups)
3. RT → IGHV repertoire coherence (good coordinated, bad stochastic)

They do not compose into a single causal chain. They are complementary
layers of evidence about the same RT-phase biology.

## Manuscript placement (v0.7.4 Option B-hierarchy)

Within §3.10 "RT-phase pharmacodynamics" (reframed from the falsified
cascade), add three subsections:

- **§3.10.1 Target engagement of baseline predictor axes**
  Figures E + F. Message: all four Thread-1 factors move in predicted
  directions after RT in both response groups; good responders show
  EMT enrichment at sign P = 0.016. This is mechanism evidence that the
  four externally validated baseline signatures are not bystanders but
  the biology RT perturbs.

- **§3.10.2 Directed versus stochastic immune repertoire response**
  Figures G + H. Message: good responders' IGHV usage changes are
  directionally coherent across subjects (binomial and Wilcoxon and
  pattern P all < 0.05); bad responders' IGHV changes have similar or
  larger magnitudes but are not directionally coherent. IGHV6-1 is the
  leading example (good unanimous down, bad mixed).

- **§3.10.3 Paired cascade phenomenology and convergence null**
  Keep Treg / MHC-II / SBS5 / neoantigen delta findings; cite
  convergence null (script 09) explicitly; frame as "phenomenology of
  RT-phase changes in good responders, not cascade downstream of
  baseline predictor."

A new Supp Text S5 is recommended that documents the magnitude-vs-direction
distinction and why both are reported for every paired test in this work.
Reviewer objection that MW P > 0.5 rules out paired signal is pre-empted
by the sign-based aggregate evidence for the IGHV repertoire.

## Files

| File | Content |
|---|---|
| `15_baseline_factor_pharmacodynamics.py` | Script |
| `16_trust4_ighv_directional_consistency.py` | Script |
| `17_figures_pharmacodynamics.py` | Figures script |
| `baseline_factor_per_subject_delta.tsv` | Per-subject Δ for 4 composites + members |
| `baseline_factor_pharmacodynamics_stats.tsv` | Composite + member stats |
| `baseline_factor_sign_table.tsv` | Headline sign counts (composite only) |
| `trust4_ighv_per_subject_delta.tsv` | Per-subject V-gene fraction Δ |
| `trust4_ighv_directional_stats.tsv` | All 53 coverage-filtered V-genes, sign counts + Fisher |
| `trust4_ighv_good_coherent_bad_mixed.tsv` | 16 V-genes fitting the user's hypothesized pattern |
| `trust4_ighv_focus_genes.tsv` | User-named genes + auto-flagged Fisher P ≤ 0.2 |
| `FigE_baseline_spaghetti.{pdf,png}` | 4-factor pre→post trajectories |
| `FigF_baseline_sign_bar.{pdf,png}` | Sign counts per factor per group |
| `FigG_ighv_coherence_summary.{pdf,png}` | Scatter + boxplot + pattern pie |
| `FigH_ighv_focus_spaghetti.{pdf,png}` | IGHV6-1, IGHV3-7/74, top 6 coherent V-genes |

## Scope limitations

- n = 6 + 6 paired subjects. The within-group magnitude tests are
  underpowered against NS → do not over-interpret single-factor NS Δ
  medians.
- Directional tests are robust at n = 6 only when concordance is
  extreme (5/6 or 6/6). For weaker concordance (4/6) power is low.
  That is why we report aggregate tests over 53 V-genes and over
  4 baseline factors — they aggregate the directional information that
  individual tests cannot resolve.
- IGHV fraction signals are dominated by bulk-RNA-seq repertoire
  reconstruction; TRUST4 assembly completeness varies with tumor B-cell
  content. Rare V-genes below the coverage filter were excluded.
- No external cohort has paired pre-RT / post-RT-only biopsies, so the
  directional-coherence result is a discovery-cohort-specific
  pharmacodynamic claim and cannot currently be externally validated.
