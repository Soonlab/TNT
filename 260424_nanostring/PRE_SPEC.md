# NanoString nCounter PanCancer Immune — Pre-specified Analysis Plan

**Freeze date**: 2026-04-24 (frozen *before* running the full analysis script; only probe coverage and CXCL13 per-sample values were inspected during sanity-check on the same date)
**Data file**: `/mnt/sda1/data/TNT/analysis/ncounter_immune_score.xlsx` (730 probes × 12 samples, housekeeping-normalized counts from nSolver-style output)

## Cohort

Paired pre / post-SC-RT × 6 subjects = 12 samples. Extreme phenotype design (pCR vs poor; intermediate grades excluded).

| Subj | Pre | Post | Response | Bin | RNA-seq paired |
|---|---|---|---|---|---|
| 2 | TNT_RNA_5 | TNT_RNA_6 | CR | good | ✓ |
| 4 | TNT_RNA_11 | TNT_RNA_12 | CR | good | ✓ |
| 14 | TNT_RNA_41 | TNT_RNA_42 | CR | good | ✓ |
| 10 | TNT_RNA_29 | TNT_RNA_30 | poor | bad | ✓ |
| 11 | TNT_RNA_32 | TNT_RNA_33 | poor | bad | **✗ (pre missing)** |
| 13 | TNT_RNA_38 | TNT_RNA_39 | poor | bad | ✓ |

## Mathematical ceiling

- 3-vs-3 MW rank-sum: C(6,3)=20 → **one-sided P_min = 0.050**, **two-sided P_min = 0.100**.
- Within-group Wilcoxon signed-rank (n=3): two-sided P_min = 0.250, one-sided P_min = 0.125.
- **One-sided tests registered below, justified by prior directional hypothesis** (RNA-seq n=12 cascade direction + external LC-CRT meta Fisher Bcell_infiltration P=0.014, Tcell_infiltration P=0.048, CD8_cytotoxic P=0.013 — all positive direction in good; see memory `project_tnt.md` §D external meta).

## Pre-registered hypotheses (one-sided: good Δ > bad Δ)

### Primary tier (Arrow 5 direct rescue)

- **P1** — Δ CXCL13 (single-gene, TLS hub)
- **P2** — Δ TLS-8 composite: mean-z(CXCL13, CCL19, CCL21, CXCR5, CCR7, SELL, LAMP3, BCL6)
- **P3** — Δ Plasma-proxy composite: mean-z(TNFRSF17, CD38, IRF4)
- **P4** — Δ Germinal-center TF composite: mean-z(BCL6, AICDA, POU2AF1)

**Decision rule for Arrow 5 rescue**:
- **Rescue**: ≥1 primary reaches one-sided P ≤ 0.05 (i.e., 3-vs-3 perfect separation) AND direction = good > bad.
- **Partial rescue**: 0 primaries reach ceiling, but ≥2 show direction concordant with RNA-seq / external meta AND one-sided P ≤ 0.10.
- **Null**: 0 primaries concordant direction with P ≤ 0.10.

### Secondary tier (lineage decomposition + platform)

- **S1** — Naive-B composite: mean-z(MS4A1, CD19, CD22, PAX5)
- **S2** — Memory-B composite: mean-z(CD27, CD79B, TNFRSF13B)
- **S3** — BAFF/APRIL axis: mean-z(TNFSF13, TNFSF13B, TNFRSF13B, TNFRSF13C, TNFRSF17)
- **S4** — NanoString Δ × RNA-seq Δ Pearson per gene (platform concordance; 5 subjects with both: 2, 4, 14, 10, 13)
- **S5** — 730-gene one-sided MW on Δ + BH-FDR (discovery, not pre-specified directional)

### Tertiary tier (cascade arrows 3/4 cross-platform replication)

- **T1** — Δ HLA-II composite × Δ Plasma-proxy Pearson (n=6); HLA-II = mean-z(HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1) where probes present
- **T2** — Δ CD8-exhaustion composite × Δ TLS-8 Pearson; CD8-exh = mean-z(PDCD1, HAVCR2, LAG3, TIGIT, CTLA4)
- **T3** — Δ NLRC5/HLA-I axis × Δ T-cell composite Pearson; HLA-I = mean-z(NLRC5, HLA-A, HLA-B, HLA-C, TAP1, TAP2, PSMB8, PSMB9); T-cell = mean-z(CD3D, CD3E, CD3G, CD8A, CD8B)

## Normalization

- Input values are already nSolver housekeeping-normalized counts (float, decimals).
- Transform: **log2(count + 1)** before z-score across the 12-sample matrix (per-gene z).
- Composite score = arithmetic mean of z-scores across constituent genes (missing gene dropped, min ≥ 2 genes required for composite).
- Δ = post_z_mean − pre_z_mean at **subject level** (per-subject Δ).

## Statistical testing

- **Between-group (3-vs-3)**: Mann-Whitney U, one-sided (good Δ > bad Δ), exact distribution (`method='exact'`).
- **Within-good signed-rank**: Wilcoxon signed-rank, two-sided and one-sided (paired; reported for context, P_min=0.125).
- **Correlation (T1-T3)**: Pearson and Spearman (both reported; n=6 small, emphasize direction).
- **FDR (S5)**: Benjamini-Hochberg across 730 probes.
- **BCa bootstrap**: 5000× for Δ effect sizes where applicable (same scheme as Table S8).

## Gene composition guarantees

- If any constituent gene of a composite is missing from the panel, document explicitly and compute the composite with remaining genes (min 2).
- NanoString 730-panel known absences (per 2026-04-22 audit): IGHG/A/M, IGKC, IGLC, MZB1, JCHAIN, XBP1, PRDM1, IGHV/K/Lv. Composite definitions above deliberately use surrogates (TNFRSF17+CD38+IRF4 for plasma).

## Output schema

All result TSVs in `tables/`:
- `meta.tsv` — 12 samples × sample_id, subject_id, timepoint, response_bin, score
- `logz_matrix.tsv` — 730 genes × 12 samples (log2+z)
- `subject_delta.tsv` — 6 subjects × gene Δ (post_z − pre_z)
- `composite_scores.tsv` — 12 rows × composite z-scores
- `P1_P4_primary.tsv` — composite | direction | good Δ mean | bad Δ mean | MW_1s_P | MW_2s_P | Wilcoxon_good_2s_P | Wilcoxon_bad_2s_P
- `S1_S3_lineage.tsv` — same schema for S1-S3
- `S4_platform_concordance.tsv` — per-gene NanoString Δ × RNA-seq Δ Pearson r / P (5 subjects)
- `S5_full_scan.tsv` — 730 genes × direction + one-sided MW P + BH q
- `T1_T3_cascade.tsv` — 3 cascade arrow correlations (n=6)
- `FINAL_VERDICT.md` — Arrow 5 path decision

## Analysts' reviewer-defense statements

- **Direction lock**: Δ > 0 for good is pre-registered based on (1) RNA-seq paired n=12 ΔTreg/ΔMHC-II/ΔCD8-exh all good > bad, (2) external 9-cohort LC-CRT meta Bcell_infiltration Fisher P=0.014 / Tcell P=0.048 / CD8 P=0.013 all in good-up direction. This pre-specifies the one-sided alternative before unblinding the NanoString data.
- **No multiple testing burden on primary**: 4 primaries share a single biological hypothesis (Arrow 5 "RT → B-cell/plasma lineage expansion in good responders"). Report all 4 without Bonferroni; interpret as corroboration pattern, not independent tests.
- **Secondary/tertiary**: no family-wise correction except BH on S5 scan.
- **Data-snooping guard**: only CXCL13 row was inspected during sanity check; composite-level P-values and platform concordance were NOT computed prior to this freeze.
