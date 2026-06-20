# Sanity-check diagnosis (v1 → v2)

This document records exactly **what changed**, **why**, and **what disagreements
between the user-supplied directive and the actual v0.7.6 / v0.7.7 manuscript
text** were resolved, so that a future reviewer (or future me) can re-derive
the same decisions without re-reading the prompt history.

---

## 1. Problem 1 — headline-pair value mismatch (v1 showed −0.17, manuscript says −0.07)

### Diagnosis

The v1 heatmap displayed **partial Spearman ρ adjusted for response group**
(column `partial_r` of `config_original_36pair.tsv`) — value `−0.168` for the
DSB × Δ CD8-cytotoxic pair.

The manuscript (v0.7.6 §3.8 / v0.7.7 §3.8) quotes the same headline pair as
**r = −0.07, P = 0.83**. This is the **plain Spearman** value (column
`spearman_r` / `spearman_p` of the same TSV).

The manuscript text contains an internal inconsistency that propagates into
this confusion:

> "An a-priori convergence test (36-pair targeted Spearman of baseline LASSO
> winners versus cascade Δ, **partial-adjusted for response and BH-corrected
> across all 36 pairs**) found 0/36 hits at P < 0.05 (1.8 expected by chance),
> 0 at FDR < 0.10; the headline pair DSB-repair baseline → CD8-cytotoxic Δ
> gives **r = −0.07, P = 0.83** (n = 12; …)."   — v0.7.6 line 146 / v0.7.7 line 166

i.e. the manuscript describes the method as partial-adjusted, but the
headline number it quotes is the **unadjusted plain** ρ. Both interpretations
yield 0/36 at P < 0.05, so the high-level claim is robust to the choice, but
only the plain value reproduces the quoted −0.07.

### Fresh re-computation (v2 script)

`scripts/compute_convergence_v2.py` re-runs Spearman from raw inputs
(`260418_add/integrated_subject_master_v2.tsv` for baselines,
`09_integration/paired_delta/paired_feature_long.tsv` for paired Δ,
`260418_add/paired_immune_delta_per_subject.tsv` for RNA Δ). Method matches
the original `260418_add/09_targeted_convergence_test.py` machinery exactly
(rank residualisation against group dummy, then Pearson on residuals).

**Sanity row** (`tables/sanity_check_headline.tsv`):

| pair | n | plain r | plain P | partial r | partial P |
|---|---|---|---|---|---|
| DSB × Δ CD8-cytotoxic | 12 | **−0.070** | **0.829** | −0.169 | 0.601 |

→ **Plain Spearman reproduces the manuscript headline (−0.07 / 0.83) exactly.**

### Resolution

v2 cells display **plain Spearman r** (and BH q on plain P). The figure's
colorbar is now labelled "Plain Spearman r" to remove the previous ambiguity.
The partial-Spearman values are still computed for every pair and stored in
`tables/convergence_36pair_used_v2.tsv` for any reviewer who wants to inspect
the adjusted estimates — the 0/36 verdict holds in both frameworks.

---

## 2. Problem 2 — cascade Δ feature list

### Discrepancy with the user-supplied directive

The user-supplied directive states:

> "Manuscript Supp Fig S9 legend 의 4 cascade Δ features 는:
>  Δ SBS5 / Δ MHC-I neoantigen binders / Δ Treg infiltration / Δ IGH clonotype count"

Cross-checking against the actual manuscript:

- **Supp Fig S9 is *not* the convergence-test figure.** Supp Fig S9 is the GEO
  cohorts overview / 9 → 5 CONSORT diagram (v0.7.7 line 366: *"Supp Fig S9.
  GEO cohorts overview + CONSORT-style 9 → 5 primary + 4 excluded exclusion
  diagram"*).
- **The convergence-test figure is Supp Fig S20.** Its caption (v0.7.7 line 377)
  reads: *"Supp Fig S20. Convergence-null detail (Fig 8 companion; §3.11). A
  36-pair convergence-test scatter sorted by |Spearman r| descending… B
  Purity-adjusted paired Δ sensitivity — Δ before vs after CNVkit tumor-purity
  correction for the **4 key cascade features (Treg, MHC-II, CD8 exhaustion,
  IGH_n)**"*.
- The manuscript's actual convergence-test cascade list is therefore
  **CD8_cytotoxic, Treg, MHC-II, IGH_n** (all RNA-paired n=12 features).
- This matches the cascade list in `260418_add/09_targeted_convergence_test.py`
  line 67: `CASC = ['CD8_cytotoxic_delta', 'Treg_delta', 'MHC_II_delta', 'IGH_n_delta']`.
- The user-specified v2 set (SBS5 / neo_binders / Treg / IGH_n) *includes* two
  WES-paired features (SBS5 n=14, neo_binders n=11) that the original
  convergence test deliberately excluded, with this rationale in the
  manuscript Methods (v0.7.7 line 71): *"The 36-pair convergence test (§3.8)
  uses the RNA-paired n = 12 because all cascade Δ quantities it uses are
  RNA-derived."*

### Decision taken

The user's directive is the authoritative request. Per the explicit
instruction *"4 cascade Δ features 로 계산: Δ SBS5 / Δ MHC-I neoantigen
binders / Δ Treg infiltration / Δ IGH clonotype count"*, v2 is built on the
**user-specified** list, **not** the manuscript-actual list. This means:

- v2 is a **new** convergence test that does not appear in the manuscript
  text as written. It draws on the same 9 baselines as the manuscript test
  (Reactome DSB/DNA Repair/HRR + Hallmark E2F/G2-M/Myc V2 + MHC-II/MSI/amp)
  but a different 4-cascade column set.
- Per-pair n is **variable** (n = 10 – 14): SBS5 uses WES-paired (n = 13–14),
  neo_binders uses pVACseq-complete (n = 10–11), Treg / IGH_n use RNA-paired
  (n = 11–12). Each cell of the heatmap annotates its own n in the
  lower-right corner.
- The manuscript's headline pair (DSB × Δ CD8-cytotoxic, r = −0.07) is
  **not** in the v2 column set, so it cannot be highlighted as a cell. A
  separate gold-bordered call-out box below the heatmap records the
  headline-pair value verbatim so that a reader following the §3.8 narrative
  can still locate the manuscript-quoted number.

### Headline-pair visual treatment

- v1: gold border around the DSB × Δ CD8-cyt cell (in-grid emphasis).
- v2: gold-bordered call-out *below* the heatmap with the verbatim manuscript
  sentence ("DSB repair × Δ CD8-cytotoxic, r = −0.07, P = 0.83 [n = 12].
  Δ CD8-cytotoxic is not one of the four v2 cascade columns; see Supp Fig
  S20A for the full 36-pair forest including CD8-cytotoxic Δ.").

This is the option the user labelled "권장" (recommended) in the directive:
"4 columns + 별도 textbox 로 headline pair 정보 표시".

---

## 3. Verdict robustness across configurations

Whichever cascade list and Spearman variant is chosen, the convergence
null holds (per pre-existing analyses + v2 re-run, all on the same n = 12
RNA-paired baseline cohort):

| Configuration | Pairs P < 0.05 | BH q < 0.05 |
|---|---|---|
| v0 (manuscript actual: 4 RNA cascades, plain) | 0/36 | 0/36 |
| v0 (manuscript actual: 4 RNA cascades, partial) | 0/36 | 0/36 |
| **v2 (user-specified: 4 cascades, plain)** | **0/36** | **0/36** |
| v2 (user-specified: 4 cascades, partial) | 0/36 | 0/36 |

The strongest single |plain r| in v2 is `Myc V2 × Δ IGH = +0.559, P = 0.059`
(BH q = 0.457), comparable to the strongest in v0 (`MHC-II × Δ CD8-cyt =
−0.483, P = 0.112`). v2 has more near-trend (P 0.06 – 0.13) pairs visually
because the WES-paired SBS5 column behaves coherently with the tumor-intrinsic
baselines — a pattern that did not survive BH and that the manuscript
de-emphasises by sticking with the RNA-paired n = 12 framework.

---

## 4. Files emitted by this v2 build

- `tables/convergence_36pair_used_v2.tsv` — long-form 36-row table with both
  plain and partial Spearman + both BH-q columns; sorted by ascending plain P.
- `tables/sanity_check_headline.tsv` — single-row computation of DSB × Δ
  CD8-cytotoxic showing plain & partial values, to anchor the −0.07 / 0.83
  reproduction.
- `tables/sanity_check_diagnosis.md` — this file.
