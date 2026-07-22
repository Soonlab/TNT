"""
Disentangle the swap-model improvement: is the AUC gain from REMOVING the cell-cycle-contaminated
CD8_proliferation, or from ADDING the purified CD8_cytotoxic?

Comparisons (all nested LOOCV + 5-fold inner tuning):
  baseline_37            : original master, CD8_proliferation kept
  drop_cd8prolif_36      : drop CD8_proliferation entirely
  add_immune_40          : original + 3 new immune sigs (CD8_proliferation kept)
  swap_cd8_37            : CD8_proliferation -> CD8_cytotoxic (rest of new sigs dropped)
  drop_prolif_add_3_39   : drop CD8_proliferation AND add CD8_cyt + Tcell + Bcell
"""
import os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score

OUT = '/data/data/TNT/analysis/260418_add'
ORIG = '/data/data/TNT/analysis/tables/integrated_subject_master.tsv'
EXT  = f'{OUT}/integrated_subject_master_v2.tsv'
DROP = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']

def make_pipe(pen):
    common = dict(solver='saga', max_iter=20000)
    if pen=='elasticnet':
        clf = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, C=0.5, **common)
    else:
        clf = LogisticRegression(penalty='l1', C=0.5, **common)
    return Pipeline([('imp', SimpleImputer(strategy='median')),
                     ('sc',  StandardScaler()),
                     ('sel', SelectKBest(score_func=f_classif, k=8)),
                     ('clf', clf)])

def nested(X, y, pen):
    grid = {'sel__k':[5,8,12], 'clf__C':[0.1,0.3,1.0,3.0]}
    loo = LeaveOneOut(); probs = np.zeros(len(y))
    for tr, te in loo.split(X):
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        gs = GridSearchCV(make_pipe(pen), grid, cv=inner, scoring='roc_auc',
                          n_jobs=4, refit=True)
        gs.fit(X[tr], y[tr])
        probs[te[0]] = gs.predict_proba(X[te])[:,1][0]
    return probs

def boot_ci(y, p, n=2000, seed=0):
    rng = np.random.RandomState(seed); a=[]
    for _ in range(n):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx]))<2: continue
        a.append(roc_auc_score(y[idx], p[idx]))
    return np.percentile(a,2.5), np.percentile(a,97.5)

def get_xy(df_in):
    df = df_in.copy()
    y = (df['response_bin']=='good').astype(int).values
    feats = [c for c in df.columns if c not in DROP]
    X = df[feats].apply(pd.to_numeric, errors='coerce').values
    return X, y, feats

orig = pd.read_csv(ORIG, sep='\t')
ext  = pd.read_csv(EXT,  sep='\t')

scenarios = {
    'baseline_37'         : orig,
    'drop_cd8prolif_36'   : orig.drop(columns=['CD8_proliferation']),
    'add_immune_40'       : ext,
    'swap_cd8_37'         : ext.drop(columns=['CD8_proliferation','Tcell_infiltration','Bcell_infiltration']
                                ).rename(columns={'CD8_cytotoxic':'CD8_proliferation'}),
    'drop_prolif_add_3_39': ext.drop(columns=['CD8_proliferation']),
}

rows = []
for label, df in scenarios.items():
    X, y, feats = get_xy(df)
    print(f'\n=== {label} : P={X.shape[1]} ===')
    for model, pen in [('LASSO','l1'), ('ElasticNet','elasticnet')]:
        p = nested(X, y, pen)
        auc = roc_auc_score(y, p)
        lo, hi = boot_ci(y, p)
        print(f'  {model:10s}  AUC={auc:.3f}  95%CI=[{lo:.3f},{hi:.3f}]')
        rows.append({'scenario':label, 'model':model, 'n_features':X.shape[1],
                     'AUC':round(auc,4), 'CI_low':round(lo,4), 'CI_high':round(hi,4)})
        pd.DataFrame({'subject_id':df['subject_id'].values, 'y':y, 'prob':p}
                    ).to_csv(f'{OUT}/probs_{label}_{model}.tsv', sep='\t', index=False)

res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/nested_cv_drop_vs_swap.tsv', sep='\t', index=False)
print('\n=== summary ===')
print(res.to_string(index=False))
