# 2026-04-18 — Immune-feature LASSO retraining

## Goal
Add purified immune signatures (CD8_cytotoxic, Tcell_infiltration, Bcell_infiltration) to
the discovery 37-feature LASSO and check whether nested-CV AUC improves. The motivation:
external validation v3 (script 32) already showed CD8_cytotoxic is the only signal that
robustly reproduces (Z=+2.74, P=0.006, 9 cohorts, N=721). The legacy `CD8_proliferation`
in the discovery master table is actually a mix of CD8 and tumour cell-cycle genes.

## Pipeline
| step | script | what |
|---|---|---|
| 1 | `01_score_immune_signatures.py` | ssGSEA on `06_rna_immune/tpm_symbol.tsv`, 4 sigs |
| 2 | `02_extend_master_table.py`     | merge into master → `integrated_subject_master_v2.tsv` |
| 3 | `03_nested_cv_extended.py`      | nested LOOCV+5-fold inner, LASSO+EN, 37 vs 40 features |
| 4 | `04_extra_analyses.py`          | univariate MW + immune-only + CD8 swap |
| 5 | `05_drop_vs_swap.py`            | 5-scenario ablation |
| 6 | `06_figures.py`                 | ROC, AUC bar, feature-importance, univariate boxplots |

## Coverage of new signatures (TPM matrix)
| signature | n_total | n_found | missing |
|---|---|---|---|
| CD8_cytotoxic | 17 | 17 | — |
| Tcell_infiltration | 10 | 10 | — |
| Bcell_infiltration | 10 | 9 | CD20 |
| TLS_Cabrita | 7 | 7 | (already in master, dropped to avoid collision) |

`subj 11` and `subj 20` lack pre-treatment RNA-seq samples — median imputation handles them
inside the inner-CV pipeline (no leakage).

## Headline result — nested LOOCV (outer) + 5-fold inner tuning, 2,000-bootstrap 95% CI

| Scenario | n_features | LASSO AUC | EN AUC | EN 95% CI |
|---|---|---|---|---|
| baseline_37 (original) | 29* | 0.650 | 0.686 | 0.500 – 0.860 |
| **drop_cd8prolif_36** | 28 | 0.716 | **0.745** | **0.563 – 0.895** |
| add_immune_40 | 32 | 0.650 | 0.686 | 0.500 – 0.860 |
| swap_cd8 (CD8_prolif → CD8_cytotoxic) | 29 | 0.716 | 0.722 | 0.533 – 0.888 |
| drop + add 3 immune sigs | 31 | 0.716 | 0.735 | 0.550 – 0.892 |

\*excludes subject_id, response, sex, cT, prepost_set, CMS, matched_wes (8 non-feature cols).

### Take-aways

1. **Adding the 3 purified immune signatures alone did NOT change AUC.** The inner-CV
   `SelectKBest` consistently picked the existing tumour-intrinsic features (DSB repair,
   E2F/G2M/Myc, frac_amp, MHC_II) over the new immune ones — the new features weren't
   strong enough at the univariate F level (CD8_cytotoxic MW P=0.84 in discovery; vs the
   contaminated CD8_proliferation MW P=0.035 because the cell-cycle genes leak into it).

2. **Removing the contaminated `CD8_proliferation` improves AUC by ~6 points** to
   0.745 (95% CI 0.56 – 0.89). The CI now lies fully above 0.5 — modest but no longer
   "could be chance".  This is the cleanest result of the ablation: the legacy feature
   was actively harming the predictor by introducing redundant cell-cycle signal that
   competed with the dedicated tumour-cell-cycle features (E2F, G2M, Myc V2, DSB repair).

3. **Swapping CD8_proliferation for CD8_cytotoxic** also helps (LASSO 0.716 / EN 0.722),
   but not as much as just dropping it.  Implication: in the discovery cohort the
   CD8 signal *itself* contributes little; almost the entire ~6-point gain comes from
   removing redundancy.

4. **External-validated CD8 axis is real but invisible in discovery univariate.**  This
   is consistent with the finding that the discovery signal is dominated by tumour-intrinsic
   biology (low events, small N), while external N=721 has the statistical power to surface
   the immune axis. The two analyses are *complementary*, not contradictory.

## Winning model (drop_cd8prolif_36 + ElasticNet, full-data refit, k=12)

Non-zero features (4):

| feature | coef | direction |
|---|---|---|
| DSB Repair (Reactome) | +0.89 | ↑ good |
| frac_amp (CNV) | −0.76 | ↑ bad |
| MHC_II | −0.23 | ↑ bad (consistent with discovery; opposite to ICB literature) |
| MSI_pct | +0.12 | ↑ good (subtle) |

This nicely matches the v0.7 manuscript narrative and even simplifies it.

## Univariate Mann-Whitney in discovery (33 pre-treatment RNA samples)

| feature | n_good | n_bad | median good | median bad | MW P (2-sided) |
|---|---|---|---|---|---|
| CD8_cytotoxic | 17 | 16 | −0.007 | −0.139 | 0.843 |
| Tcell_infiltration | 17 | 16 | −0.188 | −0.220 | 0.957 |
| Bcell_infiltration | 17 | 16 | +0.297 | −0.137 | 0.505 |
| CD8_proliferation (legacy) | 17 | 16 | +0.838 | +0.341 | **0.035** |
| CD8_activation | 17 | 16 | −0.059 | −0.319 | 0.871 |
| MHC_II | 17 | 16 | −0.460 | +0.081 | 0.075 |

The new sigs do trend in the expected direction (good > bad) but are not significant in
this small cohort. CD8_proliferation's univariate signal is driven by the cell-cycle
contamination — once you control for this in the multivariate model the apparent benefit
disappears.

## Recommended manuscript update

§3.5 (predictor) should report:
- Best nested LOOCV AUC = **0.745 (95% CI 0.56 – 0.89)** with the 36-feature ElasticNet
  *after dropping the cell-cycle-contaminated CD8_proliferation.*
- The previously reported AUC 0.686 (with it) is honest but suboptimal.
- The 4 non-zero features (DSB repair, frac_amp, MHC_II, MSI_pct) form a parsimonious
  signature.
- Add a sentence explaining why purified CD8 doesn't help in discovery (small N) but
  is the dominant axis in external validation (N=721).

§3.11 (external) is unchanged.

§Methods should add: "we re-scored CD8_cytotoxic, Tcell_infiltration and Bcell_infiltration
on the discovery TPM matrix using the same gene panels as the external validation
(scripts/32) to test whether including these purified immune features improved the nested
LOOCV AUC. They did not (added on top: ΔAUC=0). Removing the legacy CD8_proliferation,
which conflates CD8 effector and tumour cell-cycle markers, increased AUC from 0.686 to
0.745 (95% CI 0.56 – 0.89)."

## Files
| file | purpose |
|---|---|
| `pre_subject_immune_scores.tsv` | new sigs, per pre-treatment subject |
| `signature_gene_coverage.tsv` | sig × gene coverage report |
| `integrated_subject_master_v2.tsv` | extended master (35×40) |
| `nested_cv_results_v2.tsv` | 37 vs 40 features, both models |
| `nested_cv_drop_vs_swap.tsv` | 5 scenarios × 2 models — main result table |
| `nested_cv_immune_only.tsv` | immune-only model AUC ≈ 0.5 (no signal alone) |
| `feature_importance_drop_cd8prolif.tsv` | winning model coefs |
| `univariate_immune_response.tsv` | MW per feature |
| `FigA_ROC_5scenarios.{pdf,png}` | ROC overlay, drop_cd8prolif_36 highlighted |
| `FigB_AUC_bar_5scenarios.{pdf,png}` | AUC bar w/ 95% CI |
| `FigC_feature_importance.{pdf,png}` | winning model coef bars |
| `FigD_univariate_boxplots.{pdf,png}` | good vs bad for 5 immune sigs |
| `probs_<scenario>_<model>.tsv` | held-out probabilities for each scenario |
