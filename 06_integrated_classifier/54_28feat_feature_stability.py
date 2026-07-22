"""Per-fold SelectKBest feature-stability sweep on the 28-feature reference set
(drop_cd8prolif_36 scenario). Reproduces the table in Supp Text S2.3.

Outputs (under 260418_add/):
  - nested_cv_28feat_features_per_fold.tsv      : selected-feature list per outer fold per learner
  - nested_cv_28feat_feature_selection_frequency.tsv : 35-fold selection counts per learner

Sanity check on AUC: should reproduce LASSO 0.7157 / ElasticNet 0.7451 from
260418_add/nested_cv_drop_vs_swap.tsv (drop_cd8prolif_36 row).
"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score

ORIG = '/data/data/TNT/analysis/tables/integrated_subject_master.tsv'
OUT  = '/data/data/TNT/analysis/260418_add'
DROP = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']

df = pd.read_csv(ORIG, sep='\t').drop(columns=['CD8_proliferation'])
y = (df['response_bin']=='good').astype(int).values
feats = [c for c in df.columns if c not in DROP]
X = df[feats].apply(pd.to_numeric, errors='coerce').values
print(f'28-feature reference: N={X.shape[0]}, P={X.shape[1]}')

def make_pipe(pen):
    common = dict(solver='saga', max_iter=20000)
    if pen == 'elasticnet':
        clf = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, C=0.5, **common)
    else:
        clf = LogisticRegression(penalty='l1', C=0.5, **common)
    return Pipeline([('imp', SimpleImputer(strategy='median')),
                     ('sc',  StandardScaler()),
                     ('sel', SelectKBest(score_func=f_classif, k=8)),
                     ('clf', clf)])

def nested_track(X, y, feats, pen):
    grid = {'sel__k':[5,8,12], 'clf__C':[0.1,0.3,1.0,3.0]}
    loo = LeaveOneOut(); probs = np.zeros(len(y))
    fold_features = []; best_ks = []
    for tr, te in loo.split(X):
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        gs = GridSearchCV(make_pipe(pen), grid, cv=inner, scoring='roc_auc',
                          n_jobs=4, refit=True)
        gs.fit(X[tr], y[tr])
        sel = gs.best_estimator_.named_steps['sel']
        chosen = [f for f, m in zip(feats, sel.get_support()) if m]
        fold_features.append(chosen)
        best_ks.append(gs.best_params_.get('sel__k'))
        probs[te[0]] = gs.predict_proba(X[te])[:,1][0]
    return probs, fold_features, best_ks

probs_e, ff_e, ks_e = nested_track(X, y, feats, 'elasticnet')
probs_l, ff_l, ks_l = nested_track(X, y, feats, 'l1')
print(f'ElasticNet outer AUC = {roc_auc_score(y, probs_e):.4f}  (expected 0.7451)')
print(f'LASSO outer AUC      = {roc_auc_score(y, probs_l):.4f}  (expected 0.7157)')

freq_e = pd.Series([f for fold in ff_e for f in fold]).value_counts()
freq_l = pd.Series([f for fold in ff_l for f in fold]).value_counts()
freq = pd.DataFrame({
    'selection_count_ElasticNet': freq_e.reindex(feats).fillna(0).astype(int),
    'selection_count_LASSO':      freq_l.reindex(feats).fillna(0).astype(int),
}).sort_values('selection_count_ElasticNet', ascending=False)
freq.to_csv(f'{OUT}/nested_cv_28feat_feature_selection_frequency.tsv', sep='\t')

rows = []
for i, (fe, fl, ke, kl) in enumerate(zip(ff_e, ff_l, ks_e, ks_l)):
    rows.append({'fold':i, 'model':'ElasticNet', 'best_k':ke, 'features':';'.join(fe)})
    rows.append({'fold':i, 'model':'LASSO',      'best_k':kl, 'features':';'.join(fl)})
pd.DataFrame(rows).to_csv(f'{OUT}/nested_cv_28feat_features_per_fold.tsv',
                           sep='\t', index=False)
print(f'Wrote {OUT}/nested_cv_28feat_features_per_fold.tsv')
print(f'Wrote {OUT}/nested_cv_28feat_feature_selection_frequency.tsv')
