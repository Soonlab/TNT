# v0.7.6 → v0.7.7 Change Summary

**Date**: 2026-04-24
**Scope**: Integrate NanoString PanCancer Immune orthogonal validation findings (both pre-spec NULL and post-hoc extreme-phenotype immune axis replication) into the submission manuscript.

## What changed

### Title / preamble
- Kept same title (the NanoString addition is a supplementary orthogonal layer, not a re-framing of the two load-bearing axes).
- Added v0.7.7 preamble bullet listing all additions.
- Kept v0.7.6 and v0.7.5 preamble bullets for traceability.

### Abstract
- Added Results sentence on NanoString orthogonal validation — 6-subject pre-spec NULL + 23/23 pre-treatment extreme-phenotype good > bad + 6 ceiling hits including Ayers TIS and IFN-γ 6/10, platform concordance median r = + 0.75.
- Added Conclusions sentence reiterating in-house NanoString extreme-phenotype reproduction of the Thread 2 immune axis on regulatory-grade signatures.

### Keywords
- Added: NanoString nCounter; tumor inflammation signature; extreme phenotype; regulatory-grade signature.

### Abbreviations
- Added: IAE (Inflamed-Active-Effective); IBI (Inflamed-but-Ineffective); IMPRES; TIS (Tumor Inflammation Signature).

### Methods
- Added new subsection **NanoString PanCancer Immune orthogonal validation (§3.14, new in v0.7.7; Supp Text S7)** describing:
  - 6-subject extreme-phenotype selection (s2/4/14 pCR vs s10/11/13 poor)
  - nSolver housekeeping-normalised counts, log₂ + z + mean-composite scoring
  - Pre-registered primary (P1–P4), secondary (S1–S5), tertiary (T1–T3) with one-sided good > bad direction lock and 3-vs-3 MW ceiling (P_min = 0.050)
  - Decision rule (rescue / partial / null) tied to `PRE_SPEC.md` commit 235021f
  - Post-hoc exploratory extensions: 11 additional composites (Ayers TIS, IFN-γ 6/10, IMPRES, effector/suppressor balance, M1/M2, NK_a/i) × 3 timepoints (pre / post / Δ), 4 ratios, full 730-probe pre-only and post-only MW
  - IBI vs IAE classification (ΔAyers-TIS > 0), descriptive-only given n = 2 vs n = 3
- Code-and-data footer updated to include `260424_nanostring/` scripts.

### Results
- Added new full section **§3.14 Orthogonal NanoString PanCancer Immune platform test** with 7 subsections (§3.14.1 pre-spec Arrow 5 Δ null; §3.14.2 decision — Fig 8F arrow 5 remains dashed; §3.14.3 post-hoc pre-treatment 23/23 good > bad with 6 ceiling hits; §3.14.4 post-RT maintenance; §3.14.5 IBI vs IAE molecular dichotomy; §3.14.6 subject-level resolution for s4 and s11; §3.14.7 integrated interpretation and caveats).
- Added inline Table with 6 ceiling-hit composites and two ratio ceiling hits.
- Added inline Table with IAE > IBI and IBI > IAE top effect sizes.
- §3.13 DFS/OS deferred section kept unchanged at its previous position; §3.14 added after it (so §3.14 is the final Results section).

### Discussion
- Changed "Three interlocking findings" → "Four interlocking findings".
- Added fourth finding paragraph describing NanoString extreme-phenotype replication + IBI/IAE dichotomy.
- Added NanoString-specific limitation paragraph: 3-vs-3 MW ceiling, effective-independence adjustment for gene-sharing, n = 2 vs n = 3 descriptive nature of IBI vs IAE.

### Conclusion
- Added phrase on NanoString orthogonal platform test reproducing the Thread 2 immune axis on regulatory-grade signatures (23/23 composites, 6 ceiling, Ayers TIS + IFN-γ).

### Figures
- No new main figures (Fig 1–9 unchanged). Fig 8F arrow 5 remains dashed / qualified as in v0.7.6 — the NanoString pre-spec does not rescue it.

### Supplementary figures (7 new: S22–S28)
- S22 NanoString composite 6-panel figure (pre/post/Δ + canonical + IBI/IAE + cascade T2)
- S23 pre/post/Δ heatmap (enlarged)
- S24 canonical signatures (enlarged)
- S25 IBI vs IAE fingerprint
- S26 subject radar (s2/s4/s11)
- S27 pre-spec Arrow 5 primary paired paired Δ
- S28 platform concordance histogram

### Supplementary tables (4 new: S12–S15)
- S12 pre-spec Arrow 5 rescue (P1–P4 + S1–S3 + T1–T3 with MW + Wilcoxon P)
- S13 exploratory pre/post/Δ composite MW (23 composites × 3 axes + 4 ratios × 3 axes = 81 rows)
- S14 IBI vs IAE descriptive fingerprint (composite + top-40 genes = 67 rows)
- S15 Subject deep-dive (s4 and s11 × 23 composites × pre/post/Δ = 46 rows)

### Supplementary text (1 new: S7)
- S7 NanoString extreme-phenotype orthogonal validation — frozen pre-spec + post-hoc extensions + platform concordance + IBI vs IAE classification + honest caveats.

### References
- Added: [63] Ayers TIS JCI 2017; [64] Ayers IFN-γ 6/10 JCI 2017; [65] Auslander IMPRES Nat Med 2018.

### Data availability
- Added NanoString RCC + count matrix GEO deposition statement.

## What did NOT change (for reviewer transparency)

- All Fig 1–9 main figures: unchanged.
- All v0.7.6 Results sections §3.1–§3.13: unchanged.
- All 35-patient discovery RNA-seq and WES analyses: unchanged.
- Nested-LOOCV LASSO AUC 0.745: unchanged.
- External LC-CRT 9-cohort + Akiyoshi + GSE254249 meta (§3.12): unchanged.
- Treg Δ MW P = 0.026 between-group BCa CI (§3.11 cascade terminal in Path B framing): unchanged.
- Arrow 5 in Fig 8F: dashed/qualified as in v0.7.6; the NanoString pre-spec does not rescue it.
- Supp Fig S1–S21 and Supp Table S1–S11: unchanged numbering and content.
- Abstract word count approximately matches v0.7.6 (remains within GM ~350 guidance).

## File layout (under `260424_nanostring/manuscript/`)

```
TNT_manuscript_v0.7.7_GenomeMedicine.md        # full body, 517 lines
TNT_manuscript_v0.7.7_GenomeMedicine.docx      # rendered (pandoc)
v07_to_v077_DIFF.md                             # this file
figures/
  FigS22_NanoString_composite.{pdf,png}        # 6-panel summary
  FigS23_NanoString_prepostdelta_heatmap.pdf
  FigS24_NanoString_canonical_signatures.pdf
  FigS25_IBI_vs_IAE_fingerprint.pdf
  FigS26_subject_radar.pdf
  FigS27_NanoString_prespec_primary_paired.pdf
  FigS28_NanoString_platform_concordance.pdf
  build_composite_fig.py
tables/
  TableS12_NanoString_prespec_Arrow5_rescue.tsv
  TableS13_NanoString_exploratory_pre_post_delta.tsv
  TableS14_IBI_vs_IAE_fingerprint.tsv
  TableS15_subject_deepdive_NanoString.tsv
  build_supp_tables.py
patches/                                        # (currently empty; diff is documented here)
```
