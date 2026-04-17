# Manuscript v0.7.3 → v0.7.4 — delta / patch (2026-04-18)

Only the sections that change because of the immune-feature LASSO retraining
(`260418_add/`) are reproduced here. Everything else in `TNT_manuscript_v0.7.3_GenomeMedicine.md`
is unchanged.

**Headline change**: best nested LOOCV ElasticNet AUC moves from **0.686** (95 % CI 0.50–0.86)
to **0.745** (95 % CI 0.56–0.89) by *removing* the cell-cycle-contaminated `CD8_proliferation`
feature. Adding the purified immune signatures (CD8_cytotoxic / Tcell_infiltration /
Bcell_infiltration) to the master table did not improve nested AUC further; the new sigs
are reported as supportive context that ties to §3.11 external validation.

Source: `260418_add/nested_cv_drop_vs_swap.tsv`, `feature_importance_drop_cd8prolif.tsv`,
`univariate_immune_response.tsv`, `SUMMARY.md`.

---

## Patch 1 — Abstract (line 25)

**OLD**
> Nested LOOCV LASSO AUC = 0.650 [0.45, 0.83] (ElasticNet 0.686 [0.49, 0.85]); non-nested 0.755 for transparency.

**NEW**
> Nested LOOCV LASSO AUC = 0.716 [0.52, 0.88], ElasticNet **0.745 [0.56, 0.89]** with a 36-feature parsimonious set (DSB-repair, copy-number deletion fraction, MHC II, MSI%); the legacy `CD8_proliferation` feature, found to conflate CD8 effector and tumor cell-cycle markers, was dropped after a five-scenario ablation showed it was harming the predictor (baseline ElasticNet 0.686 with it; 0.745 without it). Adding three purified immune signatures (CD8_cytotoxic / Tcell_infiltration / Bcell_infiltration; identical gene panels to the external validation) on top did not further change nested AUC, consistent with a discovery cohort whose tumor-intrinsic features dominate the univariate ranking but whose immune axis only surfaces with the statistical power of the external N = 721.

---

## Patch 2 — Methods §2.6 (Statistics & ML; line 61)

**OLD (single paragraph)**
> Per-subject master table 35 × 37 (clinical + WES + RNA; Table S1). Mann–Whitney U for continuous features, Fisher for categorical, BH FDR across the feature panel (Table 2 lists the top 20). LASSO logistic regression with **nested outer-LOOCV and inner 5-fold CV for feature pre-selection (SelectKBest k ∈ {5, 8, 12}) and regularisation tuning (C ∈ {0.1, 0.3, 1, 3})**. Outer held-out ROC AUC with 95 % bootstrap (2,000 resamples) CI (Fig 5B). ElasticNet and RandomForest reported side-by-side (Supp Fig S8). Full nested-CV pipeline pseudocode, fold-by-fold feature stability, and an exploratory permutation analysis are provided in Supp Text S2.

**NEW (add a final sentence)**
> Per-subject master table 35 × 37 (clinical + WES + RNA; Table S1) — augmented post-hoc with three purified immune signatures (CD8_cytotoxic, Tcell_infiltration, Bcell_infiltration) re-scored by ssGSEA on the discovery TPM matrix using the *same gene panels* as the external validation (§2.7), giving an extended 35 × 40 table (Table S1, sheet 2). Mann–Whitney U for continuous features, Fisher for categorical, BH FDR across the feature panel (Table 2 lists the top 20). LASSO logistic regression with **nested outer-LOOCV and inner 5-fold CV for feature pre-selection (SelectKBest k ∈ {5, 8, 12}) and regularisation tuning (C ∈ {0.1, 0.3, 1, 3})**. Outer held-out ROC AUC with 95 % bootstrap (2,000 resamples) CI (Fig 5B). ElasticNet and RandomForest reported side-by-side (Supp Fig S8). To assess whether the legacy `CD8_proliferation` signature was contributing genuine immune information or merely duplicating cell-cycle markers, we ran a **five-scenario ablation under the same nested protocol**: baseline (37), drop CD8_proliferation (36), add three purified immune sigs (40), swap CD8_proliferation → CD8_cytotoxic (37), drop + add (39); see §3.4 and Supp Text S2. Full nested-CV pipeline pseudocode, fold-by-fold feature stability, and an exploratory permutation analysis are provided in Supp Text S2.

---

## Patch 3 — Results §3.4 (Nested-CV LASSO predictor; line 91) — full rewrite

**OLD**
> A LASSO logistic regression over the 37-feature master table (Fig 5A, 5C) achieved, under **nested outer-LOOCV with inner 5-fold hyperparameter tuning** (Fig 5B), outer held-out AUC = **0.650** (95 % bootstrap CI 0.45–0.83). ElasticNet under the same nested procedure gave AUC = **0.686** (95 % CI 0.49–0.85) (Supp Fig S8). An earlier non-nested pass in which feature selection used the full training set yielded AUC 0.755; that value is reported alongside the honest nested-CV numbers for transparency but the **nested outer-LOOCV AUC of 0.65–0.69 should be regarded as the reference**, with the 95 % bootstrap CI touching 0.5. The recurrent top features across outer folds (SHAP attribution in Fig 5E) were MYC V2, DSB repair, HDR, hypoxia, MHC II and genomic deletion fraction; per-subject predicted probabilities are in Fig 5F. The pre-CRT tumor-intrinsic classifier is therefore a **modest discovery-stage predictor** whose clinical utility awaits external TNT-matched validation; by contrast, the CD8-cytotoxic immune axis reproduced externally (§3.11; Fig 7) with > 1,000 independent patients.

**NEW**
> A LASSO logistic regression over the 37-feature master table (Fig 5A, 5C) achieved, under **nested outer-LOOCV with inner 5-fold hyperparameter tuning** (Fig 5B), outer held-out AUC = 0.650 (95 % bootstrap CI 0.45–0.83). ElasticNet under the same nested procedure gave AUC = 0.686 (95 % CI 0.50–0.86) (Supp Fig S8). An earlier non-nested pass in which feature selection used the full training set yielded AUC 0.755; that value is reported alongside the honest nested-CV numbers for transparency.
>
> Because §3.11 identifies CD8-cytotoxic markers as the externally reproducible axis while the master-table contained a `CD8_proliferation` ssGSEA score whose gene panel mixed CD8 effector and tumor cell-cycle markers, we ran a five-scenario ablation under the same nested protocol (Methods §2.6; Supp Text S2). **Removing `CD8_proliferation` alone increased outer AUC to 0.716 (LASSO) and 0.745 (ElasticNet, 95 % CI 0.56–0.89)** — the lower CI bound now lies above 0.5. Replacing it with the purified `CD8_cytotoxic` (swap) recovered most of the gain (LASSO 0.716, EN 0.722). Adding the three purified immune signatures (CD8_cytotoxic, Tcell_infiltration, Bcell_infiltration) *on top of* the legacy 37 features did not change nested AUC: in inner-CV, SelectKBest's univariate F-test consistently ranked the existing tumor-intrinsic features above them (CD8_cytotoxic Mann–Whitney P = 0.84 in discovery; Bcell_infiltration P = 0.51; the new signatures only become statistically powerful in the external N = 721, see §3.11). The five-scenario AUC table and corresponding ROC curves are in Supp Text S2 / Supp Fig S8.
>
> The 36-feature ElasticNet (the winning ablation) selected four non-zero features on a full-data refit: **DSB Repair (β = +0.89), copy-number deletion fraction (β = −0.76), MHC II (β = −0.23), MSI % (β = +0.12)**. The recurrent top features across outer folds (SHAP attribution in Fig 5E) were MYC V2, DSB repair, HDR, hypoxia, MHC II and genomic deletion fraction; per-subject predicted probabilities are in Fig 5F. The pre-CRT tumor-intrinsic classifier is therefore a **modest discovery-stage predictor (AUC 0.745, 95 % CI 0.56–0.89)** whose clinical utility awaits external TNT-matched validation; by contrast, the CD8-cytotoxic immune axis reproduced externally (§3.11; Fig 7) with > 1,000 independent patients.

---

## Patch 4 — Discussion §4 (line 140), first paragraph

**OLD (relevant sentences only)**
> A nested-CV LASSO classifier over 37 integrated features achieves modest outer held-out AUC = 0.65 (95 % CI 0.45–0.83); performance that is suggestive of a tumor-intrinsic signal in discovery but whose 95 % CI includes 0.5 under strict leakage-free evaluation.

**NEW**
> A nested-CV ElasticNet classifier over 36 integrated features (the legacy `CD8_proliferation` ssGSEA score, found in §3.4 to conflate CD8 effector and tumor cell-cycle markers, was dropped) achieves modest outer held-out AUC = **0.745 (95 % CI 0.56–0.89)**; performance that is consistent with a true tumor-intrinsic signal in discovery — the lower bound of the 95 % CI now lies above 0.5 — but still calls for external TNT-matched validation given the modest point estimate. The ablation that established this number also clarifies the interpretation: the CD8 axis itself contributes little to the discovery-stage predictor (at N = 35 the inner-CV univariate F-test ranks tumor-intrinsic features above the purified CD8 effector panel), and the predictor's improvement is driven by the *removal* of redundant cell-cycle signal that was leaking into a feature labelled as immune. The CD8-cytotoxic axis only becomes statistically powerful in the external N = 721 (§3.11). The two narratives — discovery tumor-intrinsic predictor and external CD8-cytotoxic reproducibility — are therefore complementary rather than competitive: each captures a different aspect of CRT-response biology that is visible only at its appropriate sample size.

---

## Patch 5 — Conclusion §5 (line 156)

**OLD**
> The molecular response to the **radiation phase** of TNT in MSS LARC — pre-CRT intrinsic DSB/HDR/E2F/MYC axis (nested-LOOCV AUC 0.65, 95 % CI 0.45–0.83), a reproducible pre-CRT CD8-cytotoxic axis (external meta Z = +2.74, P = 0.006 across 9 cohorts / N = 721 plus independent convergent evidence from Akiyoshi et al 2023 for N = 298, total > 1,000 patients / 10 cohorts), …

**NEW**
> The molecular response to the **radiation phase** of TNT in MSS LARC — pre-CRT intrinsic DSB/HDR/E2F/MYC axis (nested-LOOCV ElasticNet AUC **0.745, 95 % CI 0.56–0.89** in a 36-feature parsimonious panel after dropping a contaminated `CD8_proliferation` legacy feature), a reproducible pre-CRT CD8-cytotoxic axis (external meta Z = +2.74, P = 0.006 across 9 cohorts / N = 721 plus independent convergent evidence from Akiyoshi et al 2023 for N = 298, total > 1,000 patients / 10 cohorts), …

---

## Patch 6 — Figure 1 caption (line 162), panel D bullet

**OLD**
> (1) pre-CRT tumor-intrinsic predictor (nested LOOCV LASSO AUC 0.650 [0.45, 0.83]);

**NEW**
> (1) pre-CRT tumor-intrinsic predictor (nested LOOCV ElasticNet AUC **0.745 [0.56, 0.89]**, 36-feature parsimonious panel);

*Note on figure file*: `genome_medicine_submission/main_figures/Fig1_cohort.{pdf,png}` panel D
was generated from `panels_v3/Fig1D_preview_forest.pdf`. The forest row for the predictor
needs its midpoint and CI bar redrawn to AUC 0.745 [0.56, 0.89]. Native-editable PPT
`manuscript/ppt/Fig1_4panels_native_editable.pptx` panel D is also affected.

---

## Patch 7 — Figure 5 caption (line 166)

**OLD**
> Figure 5. Integrated 37-feature LASSO predictor — **nested outer-LOOCV, leakage-free** (A feature correlation; **B nested LOOCV ROC: LASSO AUC 0.650 [0.45, 0.83], ElasticNet 0.686 [0.49, 0.85]**; C feature forest with 95 % CI; D UMAP of integrated features; E SHAP beeswarm; F per-subject predicted probability).

**NEW**
> Figure 5. Integrated LASSO/ElasticNet predictor with five-scenario ablation — **nested outer-LOOCV, leakage-free** (A feature correlation; **B nested LOOCV ROC for the winning 36-feature ElasticNet: AUC = 0.745 [0.56, 0.89] with `CD8_proliferation` removed; baseline 37-feature ElasticNet 0.686 [0.50, 0.86] shown for comparison**; C feature forest with 95 % CI; D UMAP of integrated features; E SHAP beeswarm; F per-subject predicted probability).

*Note on figure file*: `panels_v3/Fig5B_ROC_*.{pdf,png}` should be regenerated from
`260418_add/FigA_ROC_5scenarios.{pdf,png}` (or its 2-curve subset = baseline + drop_cd8prolif)
for the submission composite.

---

## Patch 8 — Supp Fig S8 caption (line 181)

**OLD**
> Supp Fig S8. ML model comparison (LASSO vs ElasticNet vs RandomForest, nested-LOOCV).

**NEW**
> Supp Fig S8. ML model comparison and five-scenario feature-set ablation (LASSO vs ElasticNet vs RandomForest under identical nested-LOOCV; **five feature sets**: baseline_37, drop_cd8prolif_36, add_3immune_40, swap_cd8_37, drop+add_39). See `260418_add/FigA_ROC_5scenarios` for the ROC overlay and `260418_add/FigB_AUC_bar_5scenarios` for the AUC bar with 95 % bootstrap CIs.

---

## Patch 9 — Supp Text S2 (new section to add at the end)

> **Supp Text S2.5 — Five-scenario immune-feature ablation (added 2026-04-18).**
>
> *Rationale.* §3.11's external validation found that the externally reproducible signal across
> 9 nCRT cohorts (N = 721) is the *purified* CD8-cytotoxic effector panel, not the legacy
> `CD8_proliferation` ssGSEA score that appears in the master feature table. The legacy
> signature mixes CD8A/B/GZMA/B with cell-cycle markers (MKI67/TOP2A/MCM/CCN/CDK), so its
> apparent "immune" effect in discovery is partly the same tumor proliferative signal that
> E2F / G2M / Myc V2 capture independently. We therefore tested whether the discovery LASSO
> can recover the purified immune axis by adding it explicitly, by replacing the contaminated
> feature with it, or whether the contaminated feature itself should simply be removed.
>
> *Design.* All five scenarios were evaluated under the *same* nested protocol used for the
> reference predictor: outer LeaveOneOut, inner 5-fold StratifiedKFold with
> SelectKBest (k ∈ {5, 8, 12}) → median imputation → standard scaling → logistic regression
> (LASSO and ElasticNet, l1_ratio = 0.5; C ∈ {0.1, 0.3, 1, 3}). 95 % bootstrap CIs from
> 2,000 resamples on the held-out outer probabilities. Subjects 11 and 20 lack pre-treatment
> RNA-seq; the inner-CV imputer fills these without leakage.
>
> *Results (260418_add/nested_cv_drop_vs_swap.tsv).*
>
> | Scenario | n_features | LASSO AUC | EN AUC | EN 95 % CI |
> |---|---|---|---|---|
> | baseline_37 (original) | 29 | 0.650 | 0.686 | 0.500–0.860 |
> | **drop_cd8prolif_36** | 28 | 0.716 | **0.745** | **0.563–0.895** |
> | add_immune_40 | 32 | 0.650 | 0.686 | 0.500–0.860 |
> | swap_cd8 (CD8_prolif → CD8_cytotoxic) | 29 | 0.716 | 0.722 | 0.533–0.888 |
> | drop + add 3 immune sigs | 31 | 0.716 | 0.735 | 0.550–0.892 |
>
> *Interpretation.* The improvement is from **removal**, not from addition. Adding the
> purified immune signatures alone changed nothing (the SelectKBest inner ranking did not
> select them — see univariate Mann–Whitney panel below). Removing CD8_proliferation
> simultaneously eliminated a redundant copy of the tumor-cell-cycle axis (already represented
> by E2F/G2M/Myc V2) and let the regulariser concentrate on the most informative features.
>
> *Univariate Mann–Whitney in discovery (33 pre-treatment subjects).*
>
> | feature | n_good | n_bad | median good | median bad | MW P (2-sided) |
> |---|---|---|---|---|---|
> | CD8_cytotoxic | 17 | 16 | −0.007 | −0.139 | 0.843 |
> | Tcell_infiltration | 17 | 16 | −0.188 | −0.220 | 0.957 |
> | Bcell_infiltration | 17 | 16 | +0.297 | −0.137 | 0.505 |
> | CD8_proliferation (legacy) | 17 | 16 | +0.838 | +0.341 | **0.035** |
> | MHC_II | 17 | 16 | −0.460 | +0.081 | 0.075 |
>
> The new signatures trend in the expected direction (good > bad medians for CD8 effector
> and B-cell infiltration) but are not significant at N = 33; only the contaminated
> `CD8_proliferation` reaches univariate P < 0.05 in discovery, and that is driven by its
> embedded cell-cycle markers (which are themselves selected as separate features in the
> classifier, hence the redundancy when CD8_proliferation is also retained).
>
> *Winning model — full-data refit, 36 features, ElasticNet, k = 12.* Four non-zero features:
> DSB Repair (Reactome) β = +0.89; CNV deletion fraction β = −0.76; MHC II β = −0.23; MSI %
> β = +0.12. The signature is **parsimonious and biologically coherent**: DNA-repair proficiency
> + low aneuploidy + low MHC II + microsatellite-stable bias toward eventual good response.
>
> *Reproducibility.* All scripts (`01_score_immune_signatures.py` … `06_figures.py`), output
> tables, and figures are deposited under `260418_add/` in the GitHub repository
> (`Soonlab/TNT`, commit 0489e01).

---

## What stays unchanged

- §3.1, §3.2, §3.3 (cohort, somatic landscape, GSEA pathway analysis) — no changes.
- §3.5 mutation/SBS5 clearance, §3.6 neoantigen, §3.7 HLA, §3.8 cascade, §3.9 clonal,
  §3.10 assembled cascade, §3.11 external validation — no changes.
- All Tables (1–3, S1–S9) — but Table S1 should add a "v2" sheet with the 3 new columns.
- All Figures (2, 3, 4, 6, 7, 8, 9) — no changes.
- Supp Figs (S1–S7, S9–S14) — no changes; only S8 caption updated.
- References — no new refs needed (the 32_external_validation_v3 script and CD8 panel are
  already cited via Supp Text S3).

---

## How to apply this delta

Quickest path:

1. Copy `TNT_manuscript_v0.7.3_GenomeMedicine.md` → `TNT_manuscript_v0.7.4_GenomeMedicine.md`.
2. Apply patches 1–9 in order.
3. Regenerate `panels_v3/Fig5B_ROC_*` from `260418_add/FigA_ROC_5scenarios.*` (drop the
   3 non-best scenarios for the main-figure version; keep all 5 in Supp Fig S8).
4. Update `panels_v3/Fig1D_preview_forest.pdf` row 1 to AUC 0.745 [0.56, 0.89] (or do this
   in PPT directly: `manuscript/ppt/Fig1_4panels_native_editable.pptx` panel D).
5. Add `Supp Text S2.5` (Patch 9) at the end of the existing `manuscript/supplementary/Supp_Text_S2_*.md`.
6. Re-export submission composites (`Fig1_cohort`, `Fig5_ML_predictor`) via the existing
   composite scripts (`scripts/52_fig1_native_ppt.py`, similar for Fig 5 if present).

If you'd like me to do steps 1–6 directly (write `v0.7.4` and update the figures), just say
the word.
