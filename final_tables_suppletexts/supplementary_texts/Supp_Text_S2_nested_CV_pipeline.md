# Supplementary Text S2 — Nested-CV LASSO Pipeline

This document details the leakage-free cross-validation procedure underlying the AUC values reported in Section 3.4 (LASSO outer LOOCV AUC = 0.650; ElasticNet 0.686).

## S2.1 Motivation

An earlier internal analysis selected features on the full training set and then evaluated a LASSO classifier by leave-one-out cross-validation, yielding AUC = 0.755. Because the feature-selection step had access to the held-out fold, this estimate is optimistically biased. We re-implemented the predictor with full nesting to obtain a leakage-free estimate; the honest nested-CV AUC of 0.65–0.69 is the reference value used throughout the manuscript.

## S2.2 Pseudocode

```
Input:
  X  : 35 × 37 feature matrix (clinical + WES + RNA features, Table S1)
  y  : 35-vector binary response (good = 1, bad = 0)
  K_outer = leave-one-out (n_outer = 35)
  K_inner = 5-fold stratified
  k_grid = {5, 8, 12}      # number of features to retain
  C_grid = {0.1, 0.3, 1, 3} # inverse L2 strength

For i in 1..n_outer:
  X_train, X_test = X[-i], X[i]
  y_train, y_test = y[-i], y[i]

  best_score, best_k, best_C = -inf, NA, NA
  For k in k_grid:
    For C in C_grid:
      inner_scores = []
      For each inner fold (X_inner_tr, y_inner_tr, X_inner_va, y_inner_va):
        # FEATURE SELECTION inside inner training set only
        feats = SelectKBest(f_classif, k=k).fit(X_inner_tr, y_inner_tr)
        X_inner_tr_sel = feats.transform(X_inner_tr)
        X_inner_va_sel = feats.transform(X_inner_va)
        clf = LogisticRegression(penalty='l1', C=C, solver='liblinear')
        clf.fit(X_inner_tr_sel, y_inner_tr)
        inner_scores.append(roc_auc_score(y_inner_va, clf.predict_proba(X_inner_va_sel)[:,1]))
      mean_inner = mean(inner_scores)
      If mean_inner > best_score:
        best_score, best_k, best_C = mean_inner, k, C

  # Refit on full outer training set with selected hyperparameters
  feats = SelectKBest(f_classif, k=best_k).fit(X_train, y_train)
  X_train_sel = feats.transform(X_train)
  X_test_sel  = feats.transform(X_test)
  clf = LogisticRegression(penalty='l1', C=best_C, solver='liblinear')
  clf.fit(X_train_sel, y_train)
  outer_pred[i] = clf.predict_proba(X_test_sel)[:,1]
  outer_features[i] = list of selected features

Output:
  outer_auc = roc_auc_score(y, outer_pred)            # = 0.650
  outer_auc_CI = bootstrap_2000(outer_auc) = [0.45, 0.83]
  feature_stability = frequency of each feature across n_outer folds
```

## S2.3 Feature stability across outer folds

Across the 35 outer LOOCV folds, the ten features most frequently selected (frequency in parentheses) were:
MYC_targets_V2 (35/35), DSB_repair (33/35), HDR (31/35), Hypoxia (29/35), MHC_II (28/35), genomic_deletion_fraction (26/35), DNA_repair (25/35), TLS_Cabrita (23/35), G2M_checkpoint (22/35), CD8_proliferation_legacy (19/35).

The recurrence of cell-cycle/DNA-repair features (MYC V2, DSB, HDR, G2M, DNA repair) across nearly every outer fold is the basis of the discovery-stage tumor-intrinsic predictor described in §3.4. The Hypoxia and MHC II features are stable secondary contributors. Per-fold selected feature lists are tabulated in `10_ml_predictor/nested_cv_features_per_fold.tsv`.

## S2.4 Comparator models

ElasticNet (`l1_ratio` ∈ {0.3, 0.5, 0.7, 0.9}; same nested protocol) gave outer AUC = 0.686 (95% CI 0.49–0.85). RandomForest (`n_estimators` = 500, default depth, same nested protocol) gave outer AUC = 0.581. SVM-RBF was not tested. LASSO and ElasticNet show overlapping CIs and are reported jointly; RandomForest's lower performance suggests the signal is approximately linear in the selected feature space.

## S2.5 Permutation null

A 1,000-permutation null was attempted in which `y` labels were shuffled inside the outer loop and the nested pipeline re-run. Wall-clock was prohibitive (estimated 14 hours on the available compute) and the 95% bootstrap CI of the unpermuted AUC ([0.45, 0.83]) already touches 0.5, transparently bounding the predictor's discriminative power. Per the project's pragmatism principle, the permutation test was deferred to an external multi-cohort validation effort.

## S2.6 Code

`scripts/28_nested_cv_permutation.py` (final nested pipeline).
`10_ml_predictor/nested_cv_results.tsv` (per-fold predictions).
`figures/panels_v3/Fig5B_nested_ROC.{pdf,png}` (rendered ROC).

