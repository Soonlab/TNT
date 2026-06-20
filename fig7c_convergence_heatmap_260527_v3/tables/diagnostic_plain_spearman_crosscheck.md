# Diagnostic — Plain Spearman cross-check & partial implementation validation (v3)

Three numeric questions had to be settled before the v3 figure could be
trusted. Each was answered with fresh, in-script computation from raw inputs
and (where possible) compared against an independent prior result on disk.

Numbers below are produced by `scripts/diagnostic_plain_vs_partial.py`
(rerunnable; outputs to `tables/diagnostic_shared_pairs_r_match.tsv` and
`tables/diagnostic_partial_impl_match.tsv`).

---

## Q1. Does v2 / v3 plain Spearman reproduce the prior length-1 Configuration A?

**Yes — 18/18 shared pairs match exactly.**

Prior file: `convergence_repath1_260502/tables/config_A_pooled_36pair.tsv`
(Config A — pooled correlation, no group adjustment, 36 pairs).

Configurations share the *baselines × {Treg_delta, IGH_n_delta}* slice (the
2 cascades that appear in both the manuscript-actual 4-cascade list and the
user-specified v2/v3 4-cascade list). That's 9 × 2 = 18 pairs.

| Metric | Value |
|---|---|
| Shared pairs evaluated | **18 / 18** |
| Pairs matching within |Δr| < 0.01 **and** |ΔP| < 0.01 | **18 / 18** |
| max |Δr| (v3 fresh vs Config A prior) | **0.0005** |
| max |ΔP| | **< 1e-4** |

→ v3's plain Spearman pipeline is byte-identical (to rounding) with the
2026-05-02 length-1 analysis. The input TSVs and the Spearman call agree.

---

## Q2. Is "0/36 plain P < 0.05" numerically exact?

**Yes — confirmed by fresh `scipy.stats.spearmanr` on the top 8 |r| cells.**

The user flagged that v2's strongest |r| pairs (r ≈ 0.55–0.56 at n = 12, 13)
are close enough to nominal significance that the 0/36 P<0.05 claim might be
loose. Direct recomputation, no rounding:

| baseline | cascade | n | fresh r | fresh P |
|---|---|---|---|---|
| Myc Targets V2 | Δ IGH | 12 | +0.5594 | **0.05859** |
| DNA Repair (Reactome) | Δ SBS5 | 13 | −0.5365 | **0.05876** |
| DSB repair (Reactome) | Δ IGH | 12 | +0.5455 | **0.06661** |
| DNA Repair (Reactome) | Δ IGH | 12 | +0.5385 | **0.07089** |
| DSB repair (Reactome) | Δ SBS5 | 13 | −0.5117 | 0.07386 |
| frac_amp | Δ Treg | 12 | −0.4951 | 0.10172 |
| G2-M Checkpoint | Δ SBS5 | 13 | −0.4787 | 0.09796 |
| HRR (Reactome) | Δ SBS5 | 13 | −0.4649 | 0.10942 |

Distribution across the 36 user-specified cascade pairs:

| Threshold | Count |
|---|---|
| Plain P < 0.05 | **0 / 36**  ← the claim |
| Plain P < 0.06 | 2 / 36 |
| Plain P < 0.07 | 3 / 36 |

The two near-misses (P = 0.0586 and P = 0.0588) are within 1 standard
error of nominal significance, but both sit **above** 0.05. The "0 / 36 at
P < 0.05" statement is accurate; the two pairs would clear a *one-sided*
test (P_one_sided ≈ 0.029) but the manuscript uses two-sided throughout.

After BH-correction across 36 pairs, the smallest q is q = 0.457 (Myc V2 ×
Δ IGH), well above 0.05. **0 / 36 BH q < 0.05** is also correct.

---

## Q3. Is the partial Spearman implementation correct?

**Yes — manual rank-residualise vs sklearn-LR rank-residualise agree on 36/36 pairs (|Δ| < 1e-3).**

The headline pair (DSB × Δ CD8-cytotoxic, n = 12) computed by three
independent paths:

| Source | r | P |
|---|---|---|
| `convergence_repath1_260502/tables/config_original_36pair.tsv` (prior) | −0.1680 | 0.6007 |
| Manual rank-residualise-then-Pearson (script-09 style) | **−0.1685** | **0.6007** |
| sklearn LinearRegression rank-residualise-then-Pearson (independent) | **−0.1685** | **0.6007** |
| **Plain Spearman** (manuscript-quoted convention) | **−0.0700** | **0.8295** |

Three sources of truth agree on the partial value. The manuscript-quoted
**r = −0.07** is **the plain (unadjusted) Spearman**, not the partial; the
v0.7.6 / v0.7.7 text is internally inconsistent on this point (method
described as "partial-adjusted", number cited is the unadjusted one). The
0 / 36 BH q < 0.05 verdict is identical in both frameworks, so the
manuscript's high-level claim is robust to the choice.

For the full 9 × 4 = 36 user-cascade pairs, manual and sklearn agree on
**36/36** (max |Δr| < 1e-3, max |ΔP| < 1e-3; long form:
`tables/diagnostic_partial_impl_match.tsv`).

---

## Decision derived from this diagnostic

1. **Main figure (v3) uses partial Spearman.** That's the method the
   manuscript text describes, and the BH framework BH-corrects 36 pairs.
   This is `figures/Fig7C_convergence_heatmap_partial.pptx`.
2. **Diagnostic plain figure is kept** but explicitly labelled "DIAGNOSTIC
   ONLY" and bordered for archive use. This is
   `figures/diagnostic_plain_spearman_heatmap.pptx`.
3. **The headline pair has both values displayed**: the partial r = −0.17
   in the gold-bordered cell of the main figure; the plain r = −0.07
   verbatim in the gold-bordered call-out below the heatmap (so anyone
   following the §3.8 manuscript narrative can verify both numbers reproduce).
4. **No partial-Spearman implementation can yield the −0.07 the manuscript
   quotes** — confirmed by two independent reimplementations. If the
   manuscript wants the figure to display −0.07 as the headline-cell value,
   the manuscript itself must clarify that its convergence statistic is the
   *plain* Spearman ρ rather than the partial. The §3.8 / §3.10 prose would
   need a small wording edit accordingly.
