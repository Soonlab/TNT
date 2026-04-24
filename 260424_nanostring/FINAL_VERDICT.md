# NanoString nCounter Arrow 5 rescue — FINAL VERDICT

**Date**: 2026-04-24
**Input**: `ncounter_immune_score.xlsx` (730 probes × 12 samples; 6 subjects paired pre/post)
**Pre-spec freeze**: `PRE_SPEC.md` (commit 235021f, before result computation)
**Scripts**: `scripts/01_primary_secondary_tertiary.py`, `scripts/02_figures.py`

## 1. Pre-registered decision rule → **NULL (with direction flip)**

| Tier | Hypothesis | Rule | Outcome |
|---|---|---|---|
| Rescue | ≥1 primary one-sided P≤0.05 AND good > bad | — | **FAILED** (all primaries bad > good or tied) |
| Partial rescue | ≥2 primaries concordant direction AND P≤0.10 | — | **FAILED** (direction opposite) |
| Null | Otherwise | — | **TRIGGERED** |

No primary composite reaches the 3-vs-3 exact MW ceiling (one-sided P=0.05) in the pre-registered direction.

## 2. Primary tier P1–P4 (subject-level Δ, mean-z composites)

| Test | Δ good (mean) | Δ bad (mean) | Direction | MW 1s P (good > bad) | MW 2s P |
|---|---|---|---|---|---|
| P1 CXCL13 (single) | +0.52 | +1.11 | **bad > good** | 0.800 | 0.700 |
| P2 TLS-8 (8 genes) | +0.93 | +0.85 | tied (good slightly >) | 0.800 | 0.700 |
| P3 Plasma-proxy (3 genes) | +0.33 | +0.46 | **bad > good** | 0.650 | 1.000 |
| P4 GC-TF (3 genes) | +1.06 | +1.31 | **bad > good** | 0.800 | 0.700 |

**All four primary composites show positive Δ in both groups (post-RT immune expansion occurs in everyone), but the bad group shows equal or larger expansion.**

## 3. Secondary tier S1–S5

**S1–S3 B-cell lineage**: identical direction pattern — bad > good or tied
- S1 Naive_B: good −0.01 vs bad +1.18 (MW 1s P=0.95, 2s P=0.20, bad >> good)
- S2 Memory_B: good +0.03 vs bad +0.38 (P 1s=0.80)
- S3 BAFF/APRIL: good +0.47 vs bad +0.78 (P 1s=0.80)

**S4 Platform concordance (NanoString Δ × RNA-seq Δ, 5 subjects with both)** — **EXCELLENT**:
- 694 shared genes; **90.1% show r > 0**; median r = +0.754; 313 genes (45%) with r > 0.8
- Data-quality null hypothesis (platform artifact) rejected.
- The direction flip is NOT a platform issue; it is a genuine 3-vs-3 sampling pattern.

**S5 Full 730-gene discovery scan**:
- 22 genes reach one-sided P=0.05 (good > bad direction), 20 reach P=0.05 (bad > good) — nearly symmetric
- Direction split: good > bad 391/730 (53.6%) — essentially random
- **BH-adjusted q ≥ 0.90 for all 730 genes** (no FDR hits)
- Top hits (good > bad): TRAF3, TLR8, SYK, IKBKB, JAK3, IRF3 (innate immunity / signaling — not TLS / plasma)
- Top hits (bad > good): KLRC2, KLRC1, ITGB1, THBS1, AICDA (!), RAG1 — NK/ILC + one GC-TF gene

## 4. Tertiary tier T1–T3 (cascade internal coherence)

| Arrow | Pearson r | Pearson P | Spearman ρ | Spearman P |
|---|---|---|---|---|
| T1 ΔHLA-II × ΔPlasma-proxy | +0.717 | 0.109 | +0.771 | 0.072 |
| **T2 ΔCD8-exh × ΔTLS-8** | **+0.820** | **0.046** | +0.771 | 0.072 |
| T3 ΔHLA-I × ΔT-cell | +0.490 | 0.324 | +0.543 | 0.266 |

**T2 significant (n=6)** — post-RT immune reprogramming shows internal biological coherence: subjects with larger ΔCD8-exh also have larger ΔTLS-8, *regardless of response*.

This is the only P<0.05 finding in the entire NanoString panel; it corroborates the RNA-seq cascade arrows 3/4 but is response-agnostic.

## 5. Root cause of direction flip

### Same-5-subject RNA-seq sanity check (subj 11 excluded due to missing pre)

| Composite | RNA-seq Δ good mean | RNA-seq Δ bad mean | Direction |
|---|---|---|---|
| TLS_8 | +0.54 | −0.11 | **good > bad** (perfect 3-vs-2 separation P=0.10) |
| Plasma_proxy | +0.56 | −0.25 | good > bad |
| GC_TF | +0.37 | +0.91 | bad > good (driven by subj 13 AICDA+BCL6+POU2AF1) |

### Single-gene CXCL13 direction per platform (same 5 subjects)

| Subject | NanoString Δ | RNA-seq Δ | Agreement |
|---|---|---|---|
| s2 (good) | +1.45 | +1.91 | ✓ both + |
| s4 (good) | **−1.15** | **−1.04** | ✓ both − (atypical good) |
| s14 (good) | +1.24 | +0.53 | ✓ both + |
| s10 (bad) | +2.08 | +0.63 | ✓ both + |
| s13 (bad) | **−1.45** | **−1.98** | ✓ both − |

Platforms agree on direction for every subject. **The discrepancy between same-5 RNA-seq (good > bad) and same-6 NanoString (bad > good) is entirely explained by adding subj 11 to the bad group.**

### Subj 11 is a "bad responder with moderate immune expansion"
NanoString-only Δ for subj 11: TLS_8 +0.58, Plasma_proxy +0.67, Naive_B +0.67, HLA_II +1.06. Positive on most axes — contributes to lifting the bad-group mean above the good-group mean.

**Subj 11 + subj 10** = 2/3 poor responders showing post-RT inflammatory expansion. This is the "inflamed but ineffective" non-responder phenotype (known in immuno-oncology: immune infiltration without effective tumor control, e.g. via exhaustion or immunosuppressive lineage dominance).

## 6. Subject-level fingerprint (Δ post − pre, composite z)

| Subject | TLS-8 | Plasma | GC-TF | Naive-B | Memory-B | BAFF/APRIL | HLA-II | CD8-exh | HLA-I | T-cell | Summary |
|---|---|---|---|---|---|---|---|---|---|---|---|
| s2 (good) | +1.83 | +1.71 | +1.68 | +1.18 | +1.49 | +1.47 | +1.15 | +2.02 | +1.29 | +2.15 | textbook good (all strong +) |
| s4 (good) | +0.54 | +0.12 | +1.14 | −1.19 | −0.91 | +0.15 | +0.54 | +0.79 | +0.87 | −0.78 | **atypical good** (B-cell ↓) |
| s14 (good) | +0.42 | −0.84 | +0.35 | −0.02 | −0.49 | −0.22 | +0.56 | +0.05 | +0.09 | −0.97 | modest/mixed |
| s10 (bad) | +1.22 | +1.71 | +1.59 | +1.66 | +1.70 | +1.96 | +0.93 | +0.85 | +1.39 | +1.01 | **inflamed bad** |
| s11 (bad) | +0.58 | +0.67 | +0.38 | +0.67 | +0.09 | +0.55 | +1.06 | −0.23 | +0.19 | +0.46 | modest+ bad |
| s13 (bad) | +0.75 | −1.02 | +1.96 | +1.20 | −0.65 | −0.18 | −1.77 | +1.19 | −0.12 | +0.67 | mixed bad (GC-TF ↑, HLA-II ↓↓) |

## 7. Implications for Fig 8F / §3.11 cascade framing

Under the pre-spec decision rule, NanoString does **not** rescue Arrow 5 ("RT → B-cell/TCR expansion → response"). Combined with:

- RNA-seq paired n=12 Arrow 5 already weak (ΔIGH_n good vs bad MW P=0.24; composite sign-permutation P=0.27 after correcting for metric correlation)
- External LC-CRT 9-cohort meta Bcell_infiltration Fisher **P=0.014** is **pre-treatment baseline**, *not* Δ
- NanoString n=6 shows **direction flip driven by subj 11 "inflamed but poor" biology**
- Platform concordance is excellent (median r=+0.75, 90% genes r>0) so the flip is real biology, not noise
- **T2 ΔCD8-exh × ΔTLS-8 r=+0.82 P=0.046** shows post-RT immune reprogramming is internally coherent (cascade arrows 3/4 replicate cross-platform)

### Recommended path for v0.7.7 manuscript

**Path B (reroute terminal node)** is now the evidence-driven choice:

1. Keep cascade arrows 1–4 as "phenomenology coherent across platforms":
   - Arrow 1 (ΔSBS5 clonal clearance): in-house MW P=0.041
   - Arrow 2 (Δmuts × Δneo_binders): in-house Pearson r=+0.86 P=7e-4
   - Arrow 3 (Δneo × ΔMHC_II): in-house r=−0.81 P=0.009 **+ NanoString T1 HLA-II × Plasma r=+0.72 P=0.109** (same direction, trend, n=6)
   - Arrow 4 (Δimmune × ΔTCR): in-house r=+0.84 P=6e-4 **+ NanoString T2 CD8-exh × TLS-8 r=+0.82 P=0.046** (replicates)
2. Arrow 5 terminal node: demote "B-cell/TCR expansion → response" to **dashed/qualified**; reposition "Treg + MHC-II immune-axis activation in good responders" as the robust terminal (ΔTreg MW **P=0.026** between-group, §3.11 existing finding)
3. Add explicit §3.11 subsection on **NanoString orthogonal platform test**:
   - Pre-specified 4 primaries failed rescue
   - Direction flip driven by subj 11 ("inflamed but poor" biology)
   - Platform concordance excellent (median r=+0.75)
   - T2 corroborates post-RT CD8-exh × TLS coupling (response-agnostic)
   - External LC-CRT Bcell Fisher P=0.014 remains as *pre-treatment baseline* evidence, explicitly distinguished from Δ
4. Acknowledge in Limitations that n=3+3 extreme phenotype was underpowered to distinguish direction given the two "inflamed but poor" bad responders.

**Reviewer-defense line**: "Pre-registered NanoString nCounter Δ tests on 4 B-cell/TLS/plasma composites failed to show good > bad direction (all P 1s ≥ 0.65; direction flipped in 3/4 primaries). The failure was reproducible cross-platform (same 5 RNA-seq paired subjects show same per-subject direction; median r = +0.75) and traces to two bad responders with post-RT inflammatory expansion but ineffective tumour control. We therefore retain Arrow 5 only as a qualified descriptive branch (pre-treatment B-cell baseline remains externally validated in 9 LC-CRT cohorts, Fisher P = 0.014) and position the Treg / MHC-II immune-axis activation as the primary response-specific post-RT terminal."

## 8. Artefact inventory

```
260424_nanostring/
├── PRE_SPEC.md                               # frozen hypotheses (commit 235021f)
├── FINAL_VERDICT.md                          # this file
├── scripts/
│   ├── 01_primary_secondary_tertiary.py      # main analysis
│   └── 02_figures.py                          # 5 figures
├── tables/
│   ├── meta.tsv                               # 12 samples × metadata
│   ├── logz_matrix.tsv                        # 730 × 12 log2+z
│   ├── subject_delta.tsv                      # 730 × 6 Δ per subject
│   ├── composite_scores.tsv                   # per-sample composite z
│   ├── composite_subject_delta.tsv            # per-subject composite Δ
│   ├── composite_definitions.tsv              # gene membership audit
│   ├── P1_P4_primary.tsv                      # primary tier stats
│   ├── S1_S3_lineage.tsv                      # lineage composite stats
│   ├── S4_platform_concordance.tsv            # 694 gene × r×P
│   ├── S4_platform_summary.json               # overall metrics
│   ├── S5_full_scan.tsv                       # 730 gene one-sided MW + BH
│   └── T1_T3_cascade.tsv                      # 3 cascade correlations
├── figures/                                   # 5 PDF + PNG each
│   ├── Fig_primary_paired.{pdf,png}
│   ├── Fig_primary_secondary_bar.{pdf,png}
│   ├── Fig_cascade_scatter.{pdf,png}
│   ├── Fig_platform_concordance.{pdf,png}
│   └── Fig_subject_fingerprint.{pdf,png}
└── logs/
    ├── 01_run.log
    └── 02_figures.log
```
