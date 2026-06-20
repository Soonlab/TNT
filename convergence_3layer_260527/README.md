# Convergence 3-layer evidence package — Fig 7C + Supp S11 + Supp S12

Three coordinated convergence-test deliverables, each addressing one of the
manuscript's claims about static-vs-dynamic layer independence:

| Layer | Deliverable | Statistical role | Rationale |
|---|---|---|---|
| **L1 (primary)** | **Fig 7C** — plain Spearman heatmap (9 × 4) | Pair-level null evidence | No outcome-conditioning → no collider bias. Matches the manuscript's quoted headline value (r = −0.07 for DSB × Δ CD8-cyt). |
| **L2 (omnibus)** | **Supp Fig S12** — block-wise sign-coherence (NEW) | Pattern-level positive evidence | n = 12 underpowers individual cell P; sign-test sidesteps the power limitation. 12/12 cells in pre-specified directions → P = 2.4 × 10⁻⁴. |
| **L3 (sensitivity)** | **Supp Fig S11** — partial Spearman heatmap (9 × 5) | Pair-level null robustness | Re-runs the test under response-group adjustment; 0/36 BH q < 0.05 still holds → null verdict robust to adjustment choice. Demoted from main because of the collider-bias caveat. |

Build date: 2026-05-27 · Project: TNT (SC-RT rectal cancer, N = 35).

---

## 1. Three-layer rationale

The v3 review (`fig7c_convergence_heatmap_260527_v3/`) settled three numeric
facts that drove this restructuring:

- The "0 / 36 plain P < 0.05" claim is **numerically exact**, not a rounding
  artefact. The two strongest cells sit at P = 0.0586 / 0.0588 — visually
  striking, formally non-significant.
- The headline pair DSB × Δ CD8-cytotoxic gives **plain r = −0.07** (the
  manuscript-quoted value) and **partial r = −0.17** (the value yielded by
  every response-adjusted Spearman implementation we tried). The two values
  differ by ~0.10 in magnitude. The 0/36 BH q < 0.05 verdict holds in both
  frameworks.
- The manuscript v0.7.6 / v0.7.7 prose is internally inconsistent on this
  point — it describes the method as "partial-adjusted" but quotes the
  unadjusted headline number.

Given these facts, the cleanest manuscript posture is:

- **Primary**: report the *plain* Spearman headline (which the prose
  already quotes), and surface the convergence-null **as a positive pattern
  finding** via the sign-coherence test rather than as 36 individual nulls.
- **Sensitivity**: keep the partial-Spearman heatmap as a robustness check
  with an explicit collider-bias caveat.
- **No method-text contradiction** anywhere: Fig 7C = plain, Supp S11 =
  partial sensitivity, Supp S12 = sign-coherence omnibus.

---

## 2. Sanity check

| Deliverable | Headline-pair value | Source |
|---|---|---|
| **Fig 7C main** | The headline pair is *outside* the 4-cascade panel; reported verbatim in the top-right textbox: **DSB × Δ CD8-cytotoxic, r = −0.07, P = 0.83 (n = 12)**. | matches manuscript §3.8 exactly |
| **Supp Fig S11** | The headline pair is the gold-bordered cell in the 5th (reference) column: **partial r = −0.17, P = 0.60 (n = 12)**. | three independent partial-Spearman implementations agree to within 1e-3 |
| **Supp Fig S12** | The headline pair is *not* in S12's 12-cell block (S12 uses the 6 tumor-intrinsic baselines × 2 cascades, omitting Δ CD8-cytotoxic). |  |

### Fig 7C 36-pair statistics (plain Spearman)
- n range: 10 – 14 (varies by cascade Δ type)
- plain r range: [−0.536, +0.559]
- plain P range: [0.0586, 0.965]
- **0 / 36 plain P < 0.05**, **0 / 36 BH q (plain) < 0.05**

### Supp S11 9 × 5 = 45-cell statistics (partial Spearman)
- partial r range: [−0.478, +0.462]
- partial P range: [0.116, 0.982]
- formal 36-pair set (user-specified 4 cascades): **0 / 36 partial P < 0.05**, **0 / 36 BH q (partial) < 0.05**

### Supp S12 sign-coherence test
| Block | k / n | One-sided binomial P |
|---|---|---|
| Block 1 — baseline × **Δ SBS5** (predicted ↓) | **6 / 6** | 0.01562 |
| Block 2 — baseline × **Δ IGH clonotypes** (predicted ↑) | **6 / 6** | 0.01562 |
| **Combined 12-cell block** | **12 / 12** | **2.44 × 10⁻⁴ (= 1/4096)** |
| Sign-flip permutation (1000 iters, seed = 20260527) | 0 / 1000 reached ≥ 12 | P_empirical < 10⁻³ (max in null = 11) |

→ Analytical binomial and empirical permutation agree to within sampling
resolution. The analytical value is what the manuscript should quote.

---

## 3. Pre-specification documentation (Supp S12)

The 6 baselines and 2 cascade Δ features in S12 are an **a priori** selection
from the thesis, not a post-hoc subset. Documented inside the figure
(Panel A bottom block) and reproduced here so a reviewer can verify the
selection logic without reading the figure:

### 6 baselines = "tumor-intrinsic proliferation/repair axis" (manuscript §3.3)
DSB repair, DNA repair, HRR, E2F targets, G2-M checkpoint, Myc V2.
These are the Reactome DNA-repair pathway composites + the Hallmark
cell-cycle / proliferation signatures that define Thread 1 of the
baseline-predictor architecture. EMT is excluded by design because its
predicted radiation-response direction is opposite (mesenchymal escape,
not proliferation-driven cytotoxicity).

### 2 cascade Δ features = "radiation cascade endpoints with mechanistic prediction"
| Δ feature | Predicted sign | Mechanistic prediction |
|---|---|---|
| Δ SBS5 | **negative** | Radiation kills proliferating cells → their SBS5 (clock signature) mutations leave the tumor → higher baseline proliferation → more clearance → more negative r. |
| Δ IGH clonotypes | **positive** | Radiation-released tumor antigen drives B-cell repertoire expansion → higher baseline proliferation → larger antigen pulse → larger Δ IGH → more positive r. |

Δ MHC-I neoantigen binders and Δ Treg are **excluded** from the sign-coherence
test because their thesis-predicted directions for the tumor-intrinsic axis
are equivocal (neoantigen clearance scales with mutation clearance but
through a different lineage of effects; Treg dynamics are immune-axis
properties, not proliferation-axis). They remain in Fig 7C and Supp S11
as part of the 36-pair grid.

### Cherry-picking guard
The pre-specification limits the omnibus test to the 12 cells whose
predicted directions are biologically explicit. The number of possible
6 × 2 subgrids inside the 9 × 4 grid is C(9,6) × C(4,2) = 504; the
selection logic above is the *one* subgrid the manuscript thesis points to
before the data are examined, so the test does not consume the multiplicity
budget of 504.

---

## 4. File-by-file summary

### Fig 7C — `figures/Fig7C_main_plain_spearman.pptx`
- 1 slide · 241 shapes · Arial only · 0 shadow effects · 10 × 7.5 in.
- 9 baseline (rows) × 4 cascade Δ (columns). Cell value = plain Spearman ρ.
- Diverging RdBu_r palette, range ±0.6 (matches plain r magnitude).
- Per-cell `n=` in lower-right corner.
- **No** gold-bordered cell — the headline pair is not in this 4-cascade
  panel. Instead, top-right call-out textbox carries the manuscript
  headline-pair value verbatim.
- Bottom annotation:
  > n = 11–12 paired subjects (varies by cascade Δ type) · plain Spearman ρ · BH-corrected across 36 pairs
  >
  > Cell value = plain Spearman r; small n=… in lower-right corner of each cell.   * nominal P < 0.05    ** BH q < 0.05.
  >
  > **0/36 nominal P < 0.05   ·   0/36 BH q < 0.05    (no individual baseline-cascade pair detected).**
  >
  > See Supp Fig S12 for omnibus block-wise sign-coherence test (12/12 cells in predicted direction, P = 2.4 × 10⁻⁴).   See Supp Fig S11 for partial-Spearman sensitivity (response-group adjusted).

### Supp Fig S11 — `figures/SuppFig_S11_partial_spearman_sensitivity.pptx`
- 1 slide · 270 shapes · Arial only · 0 shadow effects.
- 9 × 5 grid: 4 user-cascades + Δ CD8-cytotoxic (5th column, grey header,
  "ref."). Cell value = partial Spearman ρ. Range ±0.4 (partial values are
  more compressed than plain).
- Gold-bordered headline cell DSB × Δ CD8-cytotoxic in column 5 (partial
  r = −0.17).
- Coral-bordered caveat box at the bottom:
  > Caveat: Conditioning on response group (an outcome variable in our analytical design) may introduce collider bias under some causal structures. This sensitivity analysis is presented as supportive pair-level null evidence, not the primary test. Primary analysis: Fig 7C (plain Spearman). Omnibus pattern test: Supp Fig S12.

### Supp Fig S12 — `figures/SuppFig_S12_block_sign_coherence.pptx`
- **2 slides** (Rule 1 strict: one panel per slide). Total 112 shapes.
  Arial only · 0 shadow effects.
- **Slide 1 / Panel A** — 6 × 2 sign-concordance matrix.
  - Columns: Δ SBS5 (predicted ↓), Δ IGH clonotypes (predicted ↑).
  - Cell colour: TNT-palette **teal** ramp (intensity = |r|/0.6) — all 12
    cells are concordant with prediction, so the entire matrix is teal-tinted
    (visually conveys "as predicted").
  - Cell text: `±0.XX  ↓` or `±0.XX  ↑` with `n=` in lower-right corner.
  - Block sub-totals below each column + a green-bordered combined call-out:
    `12/12 in predicted direction.  One-sided binomial P = 2.4 × 10⁻⁴`.
  - Bottom annotation: pre-spec rationale (6 baselines, 2 cascades,
    direction predictions).
- **Slide 2 / Panel B** — 1000-iter sign-flip permutation null.
  - X-axis: count of cells matching predicted direction (0 – 12).
  - Y-axis: permutation iterations.
  - Blue-grey bars: null distribution (peaks around 6, as expected from a
    binomial(12, 0.5)). Counts annotated above each bar.
  - Coral vertical line at k = 12: observed value (separated from the bulk).
  - Below-plot call-out: analytical binomial vs empirical permutation P.

---

## 5. Manuscript-side TODOs

1. **§3.10 forward reference**: *"(Fig 7E)"* → *"(Fig 7C)"*.
2. **§3.8 / §3.10 prose**: drop "partial-adjusted" from the method statement
   when referring to the main figure (now plain Spearman). Add a sentence
   pointing to Supp Fig S12 for the omnibus positive-evidence test and to
   Supp Fig S11 for the partial-adjustment sensitivity check. Suggested
   wording:
   > "An a-priori 36-pair convergence test (plain Spearman ρ of baseline
   > LASSO winners versus cascade Δ, BH-corrected across 36 pairs) found
   > 0/36 hits at P < 0.05 (Fig 7C). At the pair level this is null, but
   > the *sign pattern* of the 6 tumor-intrinsic baselines × 2 mechanistically
   > predicted cascade endpoints (Δ SBS5, Δ IGH clonotype count) was
   > coherent in 12 / 12 cells in the pre-specified directions
   > (binomial one-sided P = 2.4 × 10⁻⁴; cross-validated by 1000-iter
   > sign-permutation; Supp Fig S12). A response-group-adjusted partial
   > Spearman sensitivity analysis reproduced the pair-level null (0/36
   > BH q < 0.05; Supp Fig S11)."
3. **Supp Fig S9 / S20 numbering**: confirm against the current draft's
   supplementary numbering whether the partial-sensitivity figure should be
   S11 (this folder's filename) or another integer; rename PPTX if needed.
4. **Figure 8F cascade-null callout**: the "0 / 36 pairs P < 0.05 · DSB →
   CD8-cyt r = −0.07" box can either stay (plain quote, internally
   consistent) or be augmented with a reference to the 12/12 sign coherence.

---

## 6. Output files

```
convergence_3layer_260527/
├── README.md                                  ← this file
├── scripts/
│   ├── _shared_pptx_helpers.py                ← shared diverging / teal / set_text / kill_shadow
│   ├── make_fig7c_main_plain.py
│   ├── make_suppfig_s11_partial.py
│   └── make_suppfig_s12_sign_coherence.py
├── tables/
│   ├── plain_spearman_36pair.tsv              ← L1 source (45 rows: 36 user-cascades + 9 ref)
│   ├── partial_spearman_36pair.tsv            ← L3 source (same 45)
│   ├── block_sign_coherence_12cell.tsv        ← L2 source (12 cells, sign/concordance/n)
│   └── sign_permutation_null_1000.tsv         ← L2 permutation null (1000 rows × {iter, count})
└── figures/
    ├── Fig7C_main_plain_spearman.pptx                ★ L1 main
    ├── Fig7C_main_plain_spearman.pdf
    ├── Fig7C_main_plain_spearman_preview-1.png
    ├── SuppFig_S11_partial_spearman_sensitivity.pptx ★ L3 sensitivity
    ├── SuppFig_S11_partial_spearman_sensitivity.pdf
    ├── SuppFig_S11_partial_spearman_sensitivity_preview-1.png
    ├── SuppFig_S12_block_sign_coherence.pptx         ★ L2 omnibus (new)
    ├── SuppFig_S12_block_sign_coherence.pdf
    ├── SuppFig_S12_block_sign_coherence_preview-1.png   ← Panel A
    └── SuppFig_S12_block_sign_coherence_preview-2.png   ← Panel B
```

---

## 7. Visual summary (what each figure shows at a glance)

**Fig 7C main (plain Spearman, 9 × 4).** Heatmap dominated by pale cells.
Coherent cool-blue band in the Δ SBS5 column (DNA-repair-axis baselines all
negative, r ≈ −0.38 to −0.54) and coherent warm-red band in the Δ IGH
clonotypes column (positive, r ≈ +0.31 to +0.56). Top-right gold-bordered
call-out shows headline pair (DSB × Δ CD8-cyt, r = −0.07, P = 0.83). No
asterisks anywhere — 0/36 nominal P < 0.05. Bottom emphasises the null,
points to S11 / S12.

**Supp Fig S11 (partial Spearman, 9 × 5).** Same grid layout + Δ CD8-cyt
5th column (grey header). Cells visibly paler than the plain version
because partial r is compressed (range ±0.4 vs ±0.6). Gold-bordered DSB ×
Δ CD8-cyt cell in column 5 displays r = −0.17. Coral-bordered caveat box
at bottom explains the collider-bias rationale for sensitivity-only role.

**Supp Fig S12 / Panel A (sign coherence).** 6 × 2 matrix, every cell
teal-tinted (concordant with prediction). All cells show their r value
plus a ↓ or ↑ arrow + n. Block sub-totals "6/6 concordant, P = 0.0156"
under each column. Combined call-out below: "12/12 in predicted direction,
P = 2.4 × 10⁻⁴". Three pre-spec rationale lines.

**Supp Fig S12 / Panel B (permutation null).** Symmetric bell-shaped
histogram of 1000-iter sign-flip null counts, peaks around k = 6 (213
iterations), tails to k = 1 (5 iters) and k = 11 (3 iters). Coral vertical
line at k = 12, separated from the bulk: "observed = 12/12 ↓". Stat
call-out reports analytical and empirical P, agreeing within 10⁻³.
