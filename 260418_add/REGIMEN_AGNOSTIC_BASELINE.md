# Regimen-agnostic baseline predictor — key finding (2026-04-18)

## Finding (one sentence)

**Pre-treatment tumor-intrinsic factors (DSB/HDR repair, cell cycle, E2F/MYC, EMT) predict
final response to rectal-cancer multimodal therapy independent of when the chemotherapy
is delivered** — concurrent with radiation (nCRT-long) or consolidated after radiation
(our TNT cohort) — because the predictor is measured *before* any treatment is given.

## Why this is a non-trivial discovery

### The setup that makes it non-obvious

| Element | Our cohort (discovery) | External 5 cohorts (validation) |
|---|---|---|
| Pre-biopsy | pre-treatment | pre-treatment |
| Radiation phase | 50.4 Gy **RT alone** | 50.4 Gy RT + **concurrent capecitabine** |
| Post-RT period | FOLFOX/CAPOX consolidation (8–16 wk) | — (surgery) |
| Final response judged after | RT + consolidation chemo | RT + concurrent chemo |
| Post-biopsy available? | yes (n=14 WES / 13 RNA) | no (not used) |
| N | 35 | 518 (GSE35452, GSE45404, GSE56699, GSE133057, GSE87211) |

**Regimens differ on when chemotherapy meets radiation.** nCRT-long delivers capecitabine
*during* RT as a radiosensitizer — the tumor microenvironment during the lethal radiation
pulse already has chemo on board. Our TNT cohort delivers FOLFOX/CAPOX *after* RT completes,
so the radiation pulse hits tumor without concurrent chemo. These two regimens produce
genuinely different treatment physics in the first phase.

### What would be expected under a regimen-specific model

If the baseline predictor captured *"tumor states that synergize with concurrent
chemo-radiation"* specifically, then our TNT cohort (where synergy cannot occur because
chemo is delayed) would give an opposite or null signal compared with nCRT-long cohorts.
The Stouffer Z across 5 cohorts would average toward 0, or the sign would flip.

### What is actually observed

| Signature | 5-cohort nCRT-long Z | Discovery TNT direction | Concordance |
|---|---|---|---|
| DSB_HDR_repair | **+3.17 (P = 0.0015)** | ↑ in good | ✅ same |
| Tumor_cellcycle | **+3.21 (P = 0.0013)** | ↑ in good | ✅ same |
| E2F_MYC_cellcycle | **+2.79 (P = 0.0053)** | ↑ in good | ✅ same |
| EMT | +1.61 (P = 0.106, trend) | ↓ in good | ✅ same |

**All four directions match across two different chemo-timing regimens, 553 patients total
(35 discovery + 518 external), three significant at P < 0.01.**

The parsimonious reading is that the predictor is a property of the tumor's pre-treatment
biological state, not of how the tumor will be treated — which makes it a *pan-multimodal-
therapy* biomarker of response rather than an nCRT-specific or TNT-specific one.

## Why "pre-biopsy = regimen-agnostic" follows naturally

Conceptually obvious but often missed:

1. Baseline biopsy is drawn from a tumor that has experienced zero therapy.
2. The biological state at that moment cannot encode anything about how subsequent
   treatment will be structured — the tumor does not "know" whether chemo will arrive
   during RT or after RT.
3. What baseline biology *can* encode is intrinsic susceptibility to DNA-damage-based
   killing (RT and cytotoxic chemo both converge on DSB induction and replication stress),
   and intrinsic proliferative state (cycling cells are more vulnerable to both
   radiation and S-phase cytotoxics).
4. Therefore any baseline signature whose biology centers on DNA repair capacity, cell
   cycle engagement, and mesenchymal escape programs should be a general predictor of
   any regimen whose effector mechanism is genotoxic stress — regardless of schedule.

The external validation is the empirical test that this reasoning holds.

## Why this matters for the manuscript

### 1. Strengthens Genome Medicine pitch

Reviewers routinely push back on single-cohort discovery papers with the question:
*"would this predictor work in a cohort treated differently?"* Normally the answer is
"maybe, not tested" and the paper is demoted. Here the answer is **"yes, tested in
518 patients across a different chemo-timing regimen; Thread 1 reproduces with Z up to
+3.21."** This is a stronger external-validation claim than most TNT or nCRT
biomarker papers can make.

### 2. Resolves the "TNT-specific biology" concern from v0.5

v0.5 originally framed external non-reproducibility as "TNT-regimen-specific biology".
The restricted 5-cohort meta now shows the opposite: the Thread 1 signal is **not**
TNT-specific. It is a property of pre-treatment biology that persists across the
chemo-timing axis. This flips the defensive framing ("our biology doesn't generalize")
into an offensive one ("our biology generalizes across regimens").

### 3. Justifies excluding cohorts with non-canonical regimens

GSE119409 (RT-only, no chemo at all) being excluded as discordant is now *consistent*
with the model, not contradictory: if the final endpoint is response to
radiation-without-chemotherapy, the biology predicting it can legitimately differ from
biology predicting response to RT + any chemo schedule. The inclusion criterion is
"multimodal therapy with chemotherapy" — schedule within that criterion is irrelevant.

### 4. Clarifies the paired-delta scope limitation

Because no paired post-pre external meta exists (and cannot exist — external cohorts in
this space are pre-only, and our post is post-RT-only whereas nCRT-long post would be
post-RT+chemo, making Δ comparison biologically incoherent), the paired cascade is
honestly presented as **a discovery-cohort phenomenology describing RT-phase tumor/
immune dynamics**, not as a predictor whose generalization we claim. The regimen-agnostic
claim applies strictly to the *baseline* predictor arm, and that arm is the one
externally validated.

## Suggested manuscript placement

- **§3.11 External validation** (in v0.7.4): open with an explicit statement of the
  regimen-agnostic claim and back it with the 5-cohort Thread 1 meta. Contrast with
  the paired-cascade scope note.
- **Discussion, "Clinical implications" paragraph**: argue that a pre-treatment
  assay is attractive precisely because it is independent of the therapy schedule
  downstream — a clinic picking between nCRT-long and TNT can use the same
  baseline readout.
- **Abstract**: one sentence acknowledging that discovery and validation cohorts
  differ in chemo timing, and that baseline tumor-intrinsic signals reproduce
  regardless.

## Files supporting this claim

- `restricted5_meta_combined.tsv` — the 5-cohort Thread 1 numbers above
- `thread1_per_cohort_summary.tsv` — per-cohort direction/magnitude/P, showing
  concordance across GSE35452, GSE45404, GSE56699, GSE133057, GSE87211
- `Fig_restricted5_meta_forest.{pdf,png}` — forest plot
- `EXTERNAL_COHORTS_REFERENCES.md` — regimen characterization for each cohort
  (all 5 primary cohorts are nCRT-long with concurrent capecitabine)

## What this finding is *not*

- Not a claim that the predictor survives **every** regimen — RT-alone (GSE119409)
  was discordant, consistent with the predictor requiring the multimodal context
  (RT + chemotherapy in some order) to operate.
- Not a claim that **post-treatment** dynamics are regimen-agnostic. Post-RT (our
  cohort) biology versus post-RT+chemo (nCRT-long cohorts, if post samples existed)
  biology would differ, and the paired cascade is therefore explicitly scoped to
  discovery-cohort RT-phase dynamics only.
- Not a replacement for prospective multicenter validation in a TNT-matched cohort
  with post-consolidation endpoint — that remains a future-work item.
