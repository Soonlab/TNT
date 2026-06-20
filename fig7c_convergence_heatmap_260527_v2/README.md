# Fig 7C v2 — Convergence test heatmap (PPT-native, rebuilt)

Rebuild of `fig7c_convergence_heatmap_260527/` (v1) after two issues raised
by the user:

1. **Headline-pair value mismatch.** v1 displayed −0.17 (partial Spearman),
   but the manuscript quotes −0.07 (plain Spearman) for DSB × Δ CD8-cytotoxic.
2. **Cascade Δ feature list mismatch.** v1 used the manuscript's actual
   convergence cascades (CD8_cytotoxic / IGH_n / MHC-II / Treg); v2 switches
   to the user-specified list (SBS5 / neo_binders / Treg / IGH_n).

Full diagnosis of both issues, including a transparent record of
disagreements between the user's directive and the v0.7.6 / v0.7.7
manuscript text as found on disk, is in `tables/sanity_check_diagnosis.md`.

Build date: 2026-05-27 · Project: TNT (SC-RT rectal cancer, N = 35).

---

## 1. Sanity-check resolution

| Item | v1 outcome | v2 outcome |
|---|---|---|
| Headline pair (DSB × Δ CD8-cytotoxic) cell value | r = −0.17 (partial Spearman) | r = −0.07, P = 0.83 (plain Spearman) — reproduced from raw inputs, matches manuscript §3.8 / §3.10 verbatim |
| Method displayed in colorbar / footnote | partial Spearman, response-adjusted | **plain Spearman ρ** (manuscript-quoted convention) |
| Headline cell location | inside the heatmap (gold-bordered) | **not** in the v2 grid — Δ CD8-cytotoxic is not a v2 cascade column. Conveyed via gold-bordered call-out box below the heatmap. |

Verification was a fresh re-computation from raw inputs
(`scripts/compute_convergence_v2.py`) which prints, at start:

```
Sanity check: headline pair (DSB × CD8_cytotoxic_delta)
  n             : 12
  plain Spearman: r = -0.070, P = 0.829
  partial Spear : r = -0.169, P = 0.601
  manuscript    : r = −0.07, P = 0.83  (plain) — should match plain row above
```

→ **Headline-pair sanity now matches manuscript exactly (r = −0.07, P = 0.83).**

---

## 2. Feature-list correction

### What changed

| Slot | v1 (manuscript-actual) | v2 (user-specified) |
|---|---|---|
| Cascade 1 | Δ CD8-cytotoxic (RNA, n = 12) | **Δ SBS5** (WES-paired, n = 13–14) |
| Cascade 2 | Δ IGH clonotypes (RNA, n = 12) | **Δ MHC-I neoantigen binders** (pVACseq, n = 10–11) |
| Cascade 3 | Δ MHC-II (RNA, n = 12) | **Δ Treg infiltration** (RNA, n = 11–12) |
| Cascade 4 | Δ Treg (RNA, n = 12) | **Δ IGH clonotype count** (RNA, n = 12) |

The 9 baseline features are unchanged from v1 (the same 9 the manuscript's
convergence test uses): DSB repair, DNA repair, HRR, E2F targets, G2-M
checkpoint, Myc V2, MHC-II, MSI %, genomic amp fraction.

### Per-pair n is now variable (10 – 14)

Because the v2 cascade list mixes WES-paired SBS5/neo_binders with RNA-paired
Treg/IGH_n, per-pair n varies. The heatmap **annotates `n=…` in the
lower-right corner of every cell** so this is locally visible. Distribution
across the 36 pairs:

| n | Pair count |
|---|---|
| 10 | 6 (E2F/G2M/Myc V2/HRR/MHC-II × neo_binders, and DNA Repair × neo_binders) |
| 11 | 4 (frac_amp & MSI × neo_binders, MSI × IGH_n/Treg) |
| 12 | 16 (all Treg / IGH_n pairs) |
| 13 | 8 (all SBS5 pairs except frac_amp) |
| 14 | 1 (frac_amp × SBS5) |

### Why the user-specified list differs from the manuscript-actual list

Documented in detail in `tables/sanity_check_diagnosis.md` §2. In short:

- The user's directive cites *"Supp Fig S9 legend"*, but Supp Fig S9 is the
  GEO cohorts CONSORT diagram, **not** the convergence figure (which is
  Supp Fig S20).
- The manuscript-actual convergence test, per Supp Fig S20B caption (v0.7.6
  line 288; v0.7.7 line 377), uses the four RNA-paired cascade Δ features
  *Treg, MHC-II, CD8 exhaustion, IGH_n*.
- The Methods §2.x explicitly states *"the 36-pair convergence test uses the
  RNA-paired n = 12 because all cascade Δ quantities it uses are RNA-derived"*
  (v0.7.7 line 71). Adding SBS5 (WES) and neo_binders (pVACseq) breaks that
  uniform-n constraint.
- v2 was nevertheless built on the user-specified set per direct directive.
  The README and the diagnosis file flag this so the user can decide whether
  to (a) merge v2 into the manuscript with a revised §3.x narrative, or (b)
  keep v1's CD8_cyt-inclusive list and ignore v2.

### Headline-pair handling

Per the user's recommended ("권장") option of *"4 columns + 별도 textbox 로
headline pair 정보 표시"*: the heatmap shows **no gold-bordered cell**, and
the manuscript headline pair is conveyed by a gold-bordered call-out below
the heatmap:

> **Manuscript headline pair (§3.10, plain Spearman, n = 12):** DSB repair ×
> Δ CD8-cytotoxic, r = −0.07, P = 0.83. Δ CD8-cytotoxic is not one of the
> four v2 cascade columns; see Supp Fig S20A for the full 36-pair forest
> including CD8-cytotoxic Δ.

---

## 3. Full 36-pair r / P / BH-q values

Long form: `tables/convergence_36pair_used_v2.tsv` (sorted by ascending plain P).

Summary:

| Quantity | Value |
|---|---|
| Pairs evaluated | **36 / 36** |
| Per-pair n range | 10 – 14 |
| Plain r range | [−0.536, +0.559] |
| Plain P range | [0.059, 0.965] |
| Partial r range | [−0.355, +0.462] |
| Partial P range | [0.131, 0.982] |
| **Plain P < 0.05** | **0 / 36** |
| **Partial P < 0.05** | **0 / 36** |
| **BH q (plain) < 0.05** | **0 / 36** |
| **BH q (partial) < 0.05** | **0 / 36** |

Top 6 pairs by ascending plain P (none survive BH):

| baseline | cascade | n | plain r | plain P | partial r | partial P | BH q (plain) |
|---|---|---|---|---|---|---|---|
| Myc Targets V2 | Δ IGH | 12 | +0.559 | 0.059 | +0.462 | 0.131 | 0.457 |
| DNA Repair (Reactome) | Δ SBS5 | 13 | −0.536 | 0.059 | −0.329 | 0.273 | 0.457 |
| DSB repair (Reactome) | Δ IGH | 12 | +0.545 | 0.067 | +0.422 | 0.172 | 0.457 |
| DNA Repair (Reactome) | Δ IGH | 12 | +0.538 | 0.071 | +0.427 | 0.167 | 0.457 |
| DSB repair (Reactome) | Δ SBS5 | 13 | −0.512 | 0.074 | −0.248 | 0.413 | 0.457 |
| G2-M Checkpoint | Δ SBS5 | 13 | −0.479 | 0.098 | −0.335 | 0.264 | 0.457 |

Note the two near-trend patterns in the v2 grid: (1) DNA-repair-axis baselines
(DSB / DNA repair / HRR / E2F / G2-M / Myc V2) all show **negative r vs Δ SBS5**
— interpretable as "higher pre-CRT DNA-repair score → more SBS5 clearance",
consistent with target engagement of the radiation-induced clock signature.
(2) The same DNA-repair-axis baselines show **positive r vs Δ IGH clonotypes**
— "higher pre-CRT DNA-repair score → larger post-RT IGH expansion". Neither
pattern clears BH q at n = 12 – 13, but the directional coherence is visible
on the heatmap (clusters of cool-blue cells in the SBS5 column, warm-red in
the IGH column among the DNA-repair rows). The user / manuscript can keep
these as "directionally consistent but underpowered trend" wording or treat
the convergence as flat; either reading is statistically defensible at this n.

---

## 4. Color scale

- **Palette**: ColorBrewer RdBu_r diverging — blue (33, 102, 172) → white
  → red (178, 24, 43); anchored at 0.
- **Limits**: **[−0.6, +0.6]** (v2 widened from v1's ±0.5 to accommodate
  observed |plain r| up to 0.559).
- **Cell text colour**: white when |r| ≥ 0.45, black otherwise — WCAG-legible
  against the cell fill.
- **No cell-level border highlight**: the headline pair is no longer
  expressible as a single cell, so headline emphasis moved to a call-out box
  below the heatmap (Rule 1: one panel per slide still satisfied — call-out
  is on the same slide but visually separated and unambiguously labelled as
  external annotation, not as a heatmap cell).

---

## 5. PPT-native verification

| Element | v1 | v2 |
|---|---|---|
| Slides | 1 | 1 |
| Total shapes | 202 | **240** |
| Distinct fonts | `{'Arial'}` | `{'Arial'}` |
| Shapes with non-empty `<a:effectLst>` (shadows / glows) | 0 | **0** |
| Slide dimensions | 10.00 × 7.50 in | 10.00 × 7.50 in |
| Text-bearing shapes | ~70 | **235** (the per-cell `n=` annotation adds 36 small textboxes) |

**How to verify in PowerPoint**

1. Open `figures/Fig7C_convergence_heatmap_v2.pptx`.
2. Double-click any cell value, e.g. "+0.55" — the caret appears immediately,
   no ungroup needed.
3. Double-click "DSB repair (Reactome)" row label — same behaviour.
4. Double-click the bottom gold-bordered call-out — text editable in place.
5. Select any cell rectangle → Format Shape → confirm no shadow / glow /
   soft-edge enabled (Rule 6).
6. Format Pane → Text Options on any selected textbox → font reads "Arial"
   (Rule 5).

LibreOffice headless render at `figures/Fig7C_convergence_heatmap_v2.pdf`;
150-dpi PNG preview at `figures/Fig7C_v2_preview-1.png`; 220-dpi at
`figures/Fig7C_v2_HIRES-1.png`.

---

## 6. Output files

```
fig7c_convergence_heatmap_260527_v2/
├── README.md                                  ← this file
├── scripts/
│   ├── compute_convergence_v2.py              ← Spearman re-computation (sanity-row + 36-pair grid)
│   └── make_fig7c_v2.py                       ← python-pptx native heatmap builder
├── tables/
│   ├── convergence_36pair_used_v2.tsv         ← 36-row long form: n / plain r,P / partial r,P / BH q ×2
│   ├── sanity_check_headline.tsv              ← DSB × Δ CD8-cyt single-row reproduction (−0.07 / 0.83)
│   └── sanity_check_diagnosis.md              ← transparent record of v1→v2 changes & user-prompt cross-check
└── figures/
    ├── Fig7C_convergence_heatmap_v2.pptx      ★ deliverable (PPT-native, single slide, single panel)
    ├── Fig7C_convergence_heatmap_v2.pdf       ← LibreOffice render
    ├── Fig7C_v2_preview-1.png                 ← 150-dpi PNG
    └── Fig7C_v2_HIRES-1.png                   ← 220-dpi PNG
```

---

## 7. Manuscript-side TODOs (not modified here)

Two items the user must reconcile depending on whether v2 is adopted:

1. If v2 is adopted as the manuscript's convergence figure (replacing the
   actual v0 design):
   - §3.8 / §3.10 must update *"36-pair targeted Spearman of baseline LASSO
     winners versus cascade Δ"* to clarify that the 4 cascades are
     SBS5 / neo_binders / Treg / IGH_n, and that per-pair n now varies
     (10 – 14).
   - The Methods sentence *"The 36-pair convergence test uses the RNA-paired
     n = 12 because all cascade Δ quantities it uses are RNA-derived"* would
     need to be edited (SBS5 / neo_binders are WES-paired).
   - Supp Fig S20B's *"4 key cascade features (Treg, MHC-II, CD8 exhaustion,
     IGH_n)"* would need updating too — or S20B could be kept on the v0
     cascade list while S20A / Fig 7C use v2 (slightly awkward, would need a
     footnote).

2. If v2 is **not** adopted and the manuscript keeps the v0 (CD8_cyt / Treg /
   MHC-II / IGH_n) convergence design:
   - Use **v1** of this folder (`fig7c_convergence_heatmap_260527/`) for the
     main figure, with the change that the headline cell shows the **plain
     Spearman value (−0.07)** instead of the partial value (−0.17). This is
     a 1-line edit in `fig7c_convergence_heatmap_260527/scripts/make_fig7c.py`
     (switch `row.partial_r` → `row.spearman_r`, `row.partial_P` →
     `row.spearman_p`). Happy to do this on request.

In either case: update §3.10 last-sentence forward reference from
*"(Fig 7E)"* → *"(Fig 7C)"* per the user's prompt note.
