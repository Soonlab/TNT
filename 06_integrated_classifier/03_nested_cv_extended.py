"""
Re-run nested LOOCV (outer) + 5-fold inner CV for hyperparameter tuning, on the extended
40-feature master table that now includes purified immune signatures
(CD8_cytotoxic / Tcell_infiltration / Bcell_infiltration).

For comparability we run the SAME pipeline on:
  (A) the original 37-feature table  (integrated_subject_master.tsv)
  (B) the extended 40-feature table  (integrated_subject_master_v2.tsv)

Models: LASSO and ElasticNet logistic regression. RandomForest skipped (legacy run already
showed it was inferior to penalised LR with these features).

Outputs:
  nested_cv_results_v2.tsv             — A vs B side-by-side (AUC, 95% CI)
  nested_outer_probs_<set>_<model>.tsv — held-out probabilities per subject
  feature_importance_v2.tsv            — abs(coef) averaged across LOOCV folds for the
                                         winning ElasticNet/LASSO model on set B
"""
import os, warnings, json
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score

OUT = '/data/data/TNT/analysis/260418_add'
os.makedirs(OUT, exist_ok=True)

DROP_COLS = ['subject_id','response_bin','response_num','sex','cT','prepost_set',
             'CMS','matched_wes']

def load_xy(path):
    df = pd.read_csv(path, sep='\t')
    y  = (df['response_bin'] == 'good').astype(int).values
    feats = [c for c in df.columns if c not in DROP_COLS]
    X = df[feats].apply(pd.to_numeric, errors='coerce').values
    return df, X, y, feats

def make_pipe(penalty):
    common = dict(solver='saga', max_iter=20000)
    if penalty == 'elasticnet':
        clf = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, C=0.5, **common)
    else:
        clf = LogisticRegression(penalty='l1', C=0.5, **common)
    return Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc',  StandardScaler()),
        ('sel', SelectKBest(score_func=f_classif, k=8)),
        ('clf', clf),
    ])

GRID = {'sel__k': [5, 8, 12], 'clf__C': [0.1, 0.3, 1.0, 3.0]}

def nested_outer(make_pipe_fn, X, y, seed=0):
    loo = LeaveOneOut()
    probs = np.zeros(len(y))
    k_sel, C_sel = [], []
    for tr, te in loo.split(X):
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        gs = GridSearchCV(make_pipe_fn(), GRID, cv=inner, scoring='roc_auc',
                          n_jobs=4, refit=True)
        gs.fit(X[tr], y[tr])
        probs[te[0]] = gs.predict_proba(X[te])[:,1][0]
        k_sel.append(gs.best_params_['sel__k'])
        C_sel.append(gs.best_params_['clf__C'])
    return probs, np.median(k_sel), np.median(C_sel)

def boot_ci(y, p, n=2000, seed=0):
    rng = np.random.RandomState(seed); aucs = []
    for _ in range(n):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2: continue
        aucs.append(roc_auc_score(y[idx], p[idx]))
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)

results = []
for label, src in [('orig37', '/data/data/TNT/analysis/tables/integrated_subject_master.tsv'),
                   ('ext40',  f'{OUT}/integrated_subject_master_v2.tsv')]:
    df, X, y, feats = load_xy(src)
    print(f'\n=== {label}: N={len(y)} good={y.sum()} bad={(1-y).sum()} P={X.shape[1]}')
    for model in ('LASSO', 'ElasticNet'):
        pen = 'l1' if model == 'LASSO' else 'elasticnet'
        probs, kmed, cmed = nested_outer(lambda p=pen: make_pipe(p), X, y)
        auc = roc_auc_score(y, probs)
        lo, hi = boot_ci(y, probs)
        print(f'  {model:10s}  AUC={auc:.3f}  95%CI=[{lo:.3f},{hi:.3f}]  k_med={int(kmed)}  C_med={cmed}')
        pd.DataFrame({'subject_id': df['subject_id'], 'y': y, 'prob': probs}
                    ).to_csv(f'{OUT}/nested_outer_probs_{label}_{model}.tsv', sep='\t', index=False)
        results.append({'feature_set': label, 'model': model,
                        'n_features_total': X.shape[1],
                        'outer_AUC': round(auc,4),
                        'CI_low': round(lo,4), 'CI_high': round(hi,4),
                        'k_median': int(kmed), 'C_median': cmed})

res = pd.DataFrame(results)
res.to_csv(f'{OUT}/nested_cv_results_v2.tsv', sep='\t', index=False)
print('\n=== summary ===')
print(res.to_string(index=False))

# Feature importance: full-data LASSO + ElasticNet on extended set, abs-coef ranking
print('\nFeature importance (full-data refit on ext40 set, k=12 most permissive)...')
df, X, y, feats = load_xy(f'{OUT}/integrated_subject_master_v2.tsv')
imp_rows = []
for model, pen in [('LASSO','l1'), ('ElasticNet','elasticnet')]:
    pipe = make_pipe(pen)
    pipe.set_params(sel__k=12)
    pipe.fit(X, y)
    sel_mask = pipe.named_steps['sel'].get_support()
    sel_feats = [f for f, m in zip(feats, sel_mask) if m]
    coefs = pipe.named_steps['clf'].coef_.flatten()
    for f, c in zip(sel_feats, coefs):
        imp_rows.append({'model': model, 'feature': f, 'coef': float(c),
                         'abs_coef': float(abs(c))})
imp = pd.DataFrame(imp_rows).sort_values(['model','abs_coef'], ascending=[True,False])
imp.to_csv(f'{OUT}/feature_importance_v2.tsv', sep='\t', index=False)
print(imp.to_string(index=False))
print(f'\nWrote {OUT}/nested_cv_results_v2.tsv and feature_importance_v2.tsv')
