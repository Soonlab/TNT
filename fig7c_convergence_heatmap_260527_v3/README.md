# Fig 7C v3 — Convergence test heatmap (partial Spearman, manuscript-consistent)

Third iteration of the convergence-test heatmap. Two issues raised in the v2
review are settled here with a fresh diagnostic + a method switch:

1. **v1 displayed −0.17 (partial), manuscript quotes −0.07 (plain).**
   Identified in v2 as a manuscript internal inconsistency: method text says
   "partial-adjusted", headline value cited is the plain. v3 cross-validates
   two independent partial-Spearman implementations on the same input and
   confirms no partial computation can yield −0.07.
2. **v2 switched to plain Spearman to make the headline cell read −0.07 —
   but plain is NOT the manuscript method.**
   v3 reverts the cell statistic to **partial Spearman** (manuscript's stated
   method) while keeping the headline-pair plain value −0.07 visible verbatim
   in a gold-bordered call-out below the heatmap. A separate diagnostic
   plain-Spearman heatmap is retained for reviewer cross-check.

Build date: 2026-05-27 · Project: TNT (SC-RT rectal cancer, N = 35).

---

## 1. Diagnostic findings (Step 1)

Full text in `tables/diagnostic_plain_spearman_crosscheck.md`. Summary:

| Question | Answer |
|---|---|
| Does v2/v3 plain Spearman reproduce the 2026-05-02 length-1 Configuration A on the 18 shared baseline × {Treg, IGH_n} pairs? | **Yes — 18/18 match within |Δr| < 0.01 and |ΔP| < 1e-4.** |
| Is "0/36 nominal P < 0.05" numerically exact for plain Spearman? | **Yes — fresh `scipy.stats.spearmanr` confirms the top 2 cells are P = 0.05859 and P = 0.05876, both above 0.05.** |
| Does the manual rank-residualise partial-Spearman implementation agree with an independent sklearn-LR reimplementation? | **Yes — 36/36 pairs agree within |Δ| < 1e-3.** |
| Does any partial-Spearman computation yield the manuscript-quoted r = −0.07 for the headline pair? | **No.** Partial r = −0.169, P = 0.601 (three sources agree). −0.07 is the plain unadjusted ρ. |

→ The "0/36 plain P < 0.05" claim is **correct**, not a rounding artefact.
→ The manuscript's headline "r = −0.07, P = 0.83" is the **plain** Spearman, not the partial — an internal inconsistency in the v0.7.6 / v0.7.7 prose (method says "partial", number says "plain").
→ Both frameworks yield 0/36 BH q < 0.05, so the high-level convergence-null verdict is robust.

---

## 2. Method decision

| Variant | Used where | Reason |
|---|---|---|
| **Partial Spearman ρ (response-group adjusted)** | **Main figure: `figures/Fig7C_convergence_heatmap_partial.pptx`** | Manuscript §3.8 / §3.10 method text: *"36-pair targeted Spearman of baseline LASSO winners versus cascade Δ, partial-adjusted for response and BH-corrected across all 36 pairs"*. |
| Plain Spearman ρ (no adjustment) | Diagnostic only: `figures/diagnostic_plain_spearman_heatmap.pptx` | Matches the manuscript's *quoted headline value* (r = −0.07). Useful for reviewer cross-check; not the formal test. |

The main figure is the partial heatmap. The diagnostic heatmap is kept in
the same `figures/` directory with a bright "DIAGNOSTIC ONLY" banner at the
top so it cannot be accidentally mistaken for the main figure.

---

## 3. Figure design — main partial heatmap

### Grid
- 9 baseline tumor-intrinsic features (rows) × **5 cascade Δ columns** (option 2 of v2 prompt).
- First 4 cascade columns (`Δ SBS5`, `Δ MHC-I neoantigen binders`, `Δ Treg`,
  `Δ IGH clonotypes`) = the user-specified formal 36-pair convergence test.
  These column headers are bold.
- 5th column (`Δ CD8-cytotoxic`) = manuscript headline-pair anchor. Its
  column-header background is light grey to mark it as **outside the formal
  36-pair test** (shown for cross-reference only). The header is non-bold.

### Headline pair highlight
- DSB repair × Δ CD8-cytotoxic cell carries a **gold border** (2.5 pt).
- The cell displays its partial r value (−0.17). A gold-bordered call-out
  below the heatmap reproduces the manuscript-quoted plain-Spearman value
  verbatim:

> Headline pair (gold-bordered cell, 5th col.): DSB repair × Δ CD8-cytotoxic.
> Partial Spearman r = −0.17, P = 0.60 (n = 12). Manuscript text §3.8 quotes
> the plain-Spearman value r = −0.07, P = 0.83 (see Supp Fig S20A forest,
> and v3 diagnostic plain heatmap).

### Color scale
- Diverging ColorBrewer RdBu_r: deep blue (33, 102, 172) → white →
  deep red (178, 24, 43). Anchored at 0.
- Limits: **±0.4** (compressed from the plain figure's ±0.6, because partial
  values are smaller — strongest observed partial |r| is 0.48 for
  MHC-II × Δ Treg; ±0.4 makes the null-dominated cell distribution visually
  obvious without saturating).
- Cell text colour: white when |partial r| ≥ 0.30, black otherwise.

### Bottom annotation block
```
n = 10–14 paired subjects (varies by cascade Δ type) · partial Spearman ρ
(response-group-adjusted; manuscript-consistent) · BH-corrected across 36
pairs (4 user-specified cascades)

Cell value = partial Spearman r;  small n=… in lower-right corner of each
cell.   * nominal P < 0.05    ** BH q < 0.05.   5th column (grey header)
shown for headline-pair cross-reference; outside the 36-pair test.

Convergence result: 0/36 partial P < 0.05  ·  0/36 BH q < 0.05    (baseline
level and radiation-phase dynamics statistically independent).
```

---

## 4. Sanity check (must match manuscript exactly)

| Pair | n | Source | r | P |
|---|---|---|---|---|
| **DSB × Δ CD8-cytotoxic** (headline) | 12 | Plain Spearman, fresh | **−0.070** | **0.829** |
|  | 12 | Plain Spearman, manuscript quote | **−0.07** | **0.83** |
|  | 12 | Partial Spearman, fresh | −0.169 | 0.601 |
|  | 12 | Partial Spearman, manuscript stated method | (not quoted, but should be −0.17 / 0.60) |

→ **Plain Spearman reproduces the manuscript-quoted headline (−0.07 / 0.83) exactly.**
→ Partial Spearman gives −0.17 / 0.60. The manuscript's *method statement*
  says "partial", so −0.17 is what v3 displays in the gold-bordered cell;
  the *manuscript-quoted number* (−0.07) is then displayed in the call-out
  for transparency.
→ **Either way, 0/36 BH q < 0.05 holds; the convergence-null claim is unchanged.**

### 36-pair / 45-cell grid statistics

For the **formal 36-pair convergence test** (4 user-specified cascades):

| Quantity | Value |
|---|---|
| Pairs evaluated | 36 / 36 |
| Per-pair n range | 10 – 14 |
| Partial r range | [−0.478, +0.462] |
| Partial P range | [0.116, 0.982] |
| **Partial P < 0.05** | **0 / 36** |
| **BH q (partial) < 0.05** | **0 / 36** |
| Plain r range (diagnostic) | [−0.536, +0.559] |
| Plain P range (diagnostic) | [0.059, 0.965] |
| Plain P < 0.05 (diagnostic) | **0 / 36** |
| BH q (plain) < 0.05 (diagnostic) | **0 / 36** |

For the 5th column (Δ CD8-cytotoxic, 9 cells, displayed for reference only):

| Quantity | Value |
|---|---|
| Partial r range | [−0.478, +0.231] |
| Partial P range | [0.116, 0.914] |
| Partial P < 0.05 | 0 / 9 |
| BH q (partial, within 9-cell set) < 0.05 | 0 / 9 |
| Plain r range | [−0.483, +0.175] |
| Plain P < 0.05 | 0 / 9 |

---

## 5. PPT-native verification

| | Main (partial) | Diagnostic (plain) |
|---|---|---|
| Slides | 1 | 1 |
| Total shapes | **269** | **268** |
| Distinct fonts | `{'Arial'}` | `{'Arial'}` |
| Shapes with non-empty `<a:effectLst>` | **0** | **0** |
| Slide dimensions | 10.00 × 7.50 in | 10.00 × 7.50 in |

**How to verify in PowerPoint**

1. Open `figures/Fig7C_convergence_heatmap_partial.pptx`.
2. Double-click the headline cell (gold border in column 5) — value
   "−0.17" is selected as editable text, no ungroup needed.
3. Double-click "DSB repair (Reactome)" row label — same behaviour.
4. Double-click bottom gold-bordered call-out — editable in place.
5. Select any cell rectangle → Format Shape → confirm no shadow / glow /
   soft-edge enabled (Rule 6).
6. Font Pane on any text shape → font "Arial" (Rule 5).

Repeat for the diagnostic deck.

---

## 6. Visual cell-color distribution (what to expect when you open each PPT)

**Main partial heatmap.** Range ±0.4 + colours compressed: only 4 cells reach
white-text threshold (|r| ≥ 0.30) — MHC-II × Δ Treg (−0.48, deep blue), MSI %
× Δ Treg (+0.45, mid red), MHC-II × Δ CD8-cyt (−0.48, deep blue), Genomic
amp × Δ Treg (−0.36, mid blue). The rest of the grid is pale (|r| < 0.30),
making the null-dominated pattern visually unambiguous. No asterisks
anywhere (0/45 cells significant). DNA-repair-axis baselines × Δ IGH
clonotypes column shows a soft warm band (r ≈ +0.18–0.46) but no cell clears
nominal P. Δ MHC-I neoantigen binders column is essentially all near-white
(|r| ≤ 0.22 throughout, after partial adjustment).

**Diagnostic plain heatmap.** Range ±0.6: same near-cells appear darker
because plain r is unadjusted. DSB / DNA Repair / HRR / E2F / G2-M / Myc V2
rows show coherent **negative** band in the Δ SBS5 column (r ≈ −0.38 to
−0.54) and coherent **positive** band in the Δ IGH clonotypes column
(r ≈ +0.20 to +0.56). Headline cell (col 5) reads "−0.07" matching the
manuscript verbatim. Top 2 cells (Myc V2 × Δ IGH, DNA Repair × Δ SBS5) are
P = 0.0586 / 0.0588 — visually striking but not nominally significant.

---

## 7. Output files

```
fig7c_convergence_heatmap_260527_v3/
├── README.md                                  ← this file
├── scripts/
│   ├── diagnostic_plain_vs_partial.py         ← cross-check (Q1, Q2, Q3 above)
│   └── make_fig7c_partial.py                  ← main + diagnostic figure builder (single pass)
├── tables/
│   ├── convergence_36pair_partial.tsv         ← long-form 45 rows (36 formal + 9 reference); partial + plain values; BH-q both
│   ├── convergence_36pair_plain.tsv           ← same 45 rows sorted by plain P (diagnostic)
│   ├── diagnostic_plain_spearman_crosscheck.md← human-readable diagnostic (Q1/Q2/Q3)
│   ├── diagnostic_shared_pairs_r_match.tsv    ← 18-row plain Spearman cross-check vs prior Config A
│   └── diagnostic_partial_impl_match.tsv      ← 36-row manual vs sklearn partial Spearman implementation match
└── figures/
    ├── Fig7C_convergence_heatmap_partial.pptx ★ MAIN deliverable (partial Spearman, manuscript-method-consistent)
    ├── Fig7C_convergence_heatmap_partial.pdf
    ├── Fig7C_partial_preview-1.png            ← 200-dpi preview
    ├── diagnostic_plain_spearman_heatmap.pptx ← DIAGNOSTIC ONLY (plain Spearman, manuscript-quoted-value-consistent)
    ├── diagnostic_plain_spearman_heatmap.pdf
    └── diagnostic_plain_preview-1.png         ← 200-dpi preview
```

---

## 8. Manuscript-side TODOs (not modified here)

1. §3.10 last sentence: change forward reference *"(Fig 7E)"* → *"(Fig 7C)"*.
2. The §3.8 / §3.10 prose has an internal inconsistency on the convergence
   test — *"partial-adjusted for response"* (method) vs *"DSB-repair →
   CD8-cytotoxic Δ gives r = −0.07, P = 0.83"* (plain unadjusted value). The
   manuscript should resolve this either by:
   - **(a)** Changing the quoted headline to the partial value: *"DSB-repair
     → CD8-cytotoxic Δ partial r = −0.17, P = 0.60 (plain r = −0.07,
     P = 0.83 reported in Supp Fig S20A)"*. Recommended — matches Fig 7C v3.
   - **(b)** Changing the method statement to *"plain Spearman of baseline
     × cascade Δ, BH-corrected across 36 pairs"* and dropping the partial-
     adjustment language. Slightly less rigorous but matches the quoted
     headline directly.
   Either resolution is statistically fine; the 0/36 BH q<0.05 verdict
   holds in both frameworks.
3. The user's earlier prompts referenced *"Supp Fig S9 legend"* as the
   convergence figure — in v0.7.6 / v0.7.7 manuscript, Supp Fig S9 is the
   GEO cohorts CONSORT diagram and the convergence figure is Supp Fig S20.
   The cross-reference labels in the figure caption / forward references
   should be checked against the current manuscript draft (v0.7.7+).
4. The v3 figure uses 5 cascade columns (manuscript's 4 + reference Δ
   CD8-cyt). If kept this way the Supp Fig S20A forest plot caption should
   note that the main-figure version (Fig 7C) shows 5 cascades while the
   forest shows all 36 user-specified pairs.
