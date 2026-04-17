"""
Negative-result follow-ups:

The base nested-CV (script 03) showed adding purified immune sigs did NOT improve outer-AUC
(LASSO 0.650, EN 0.686 — unchanged). This is because SelectKBest on the inner training folds
ranks the new sigs below the existing tumor-intrinsic features by univariate F.

To make this finding honest and complete, we add three follow-up analyses:

1. UNIVARIATE response association of the new features in the discovery cohort
   (Mann-Whitney U, good vs bad), with sample sizes documented.

2. IMMUNE-ONLY model: nested-CV using ONLY the new 3 immune features + age/sex/cT covariates.
   If immune axis predicts at all, it should at least be > 0.5 on its own.

3. SWAP model: replace the cell-cycle-contaminated CD8_proliferation with the purified
   CD8_cytotoxic, then re-run nested-CV. Tests whether the contamination is helping or
   hurting the legacy 37-feature predictor.

Outputs:
  univariate_immune_response.tsv
  nested_cv_immune_only.tsv
  nested_cv_cd8_swap.tsv
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score

OUT = '/data/data/TNT/analysis/260418_add'
SRC = f'{OUT}/integrated_subject_master_v2.tsv'
df = pd.read_csv(SRC, sep='\t')
y  = (df['response_bin'] == 'good').astype(int).values
print(f'N={len(y)} good={y.sum()} bad={(1-y).sum()}')

# ---------- 1. Univariate Mann-Whitney for the new sigs ----------
new_feats = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration']
rows = []
for f in new_feats:
    g = df[df.response_bin=='good'][f].dropna().values
    b = df[df.response_bin=='bad'][f].dropna().values
    u = stats.mannwhitneyu(g, b, alternative='two-sided')
    rows.append({'feature': f, 'n_good': len(g), 'n_bad': len(b),
                 'median_good': round(np.median(g),3), 'median_bad': round(np.median(b),3),
                 'delta': round(np.median(g)-np.median(b),3),
                 'mw_U': round(u.statistic,1), 'mw_p_2sided': round(u.pvalue,4)})
# Plus reference: existing CD8_proliferation
for f in ['CD8_proliferation','CD8_activation','MHC_II']:
    g = df[df.response_bin=='good'][f].dropna().values
    b = df[df.response_bin=='bad'][f].dropna().values
    u = stats.mannwhitneyu(g, b, alternative='two-sided')
    rows.append({'feature': f, 'n_good': len(g), 'n_bad': len(b),
                 'median_good': round(np.median(g),3), 'median_bad': round(np.median(b),3),
                 'delta': round(np.median(g)-np.median(b),3),
                 'mw_U': round(u.statistic,1), 'mw_p_2sided': round(u.pvalue,4)})
uni = pd.DataFrame(rows)
uni.to_csv(f'{OUT}/univariate_immune_response.tsv', sep='\t', index=False)
print('\n=== univariate response association ===')
print(uni.to_string(index=False))

# ---------- common nested-CV helpers ----------
def make_pipe(penalty, k=8, C=0.5):
    common = dict(solver='saga', max_iter=20000)
    if penalty == 'elasticnet':
        clf = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, C=C, **common)
    else:
        clf = LogisticRegression(penalty='l1', C=C, **common)
    return Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc',  StandardScaler()),
        ('sel', SelectKBest(score_func=f_classif, k=k)),
        ('clf', clf),
    ])

def nested(X, y, penalty, ks):
    grid = {'sel__k': ks, 'clf__C': [0.1, 0.3, 1.0, 3.0]}
    loo = LeaveOneOut()
    probs = np.zeros(len(y))
    for tr, te in loo.split(X):
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        gs = GridSearchCV(make_pipe(penalty), grid, cv=inner, scoring='roc_auc', n_jobs=4, refit=True)
        gs.fit(X[tr], y[tr])
        probs[te[0]] = gs.predict_proba(X[te])[:,1][0]
    return probs

def boot_ci(y, p, n=2000, seed=0):
    rng = np.random.RandomState(seed); a = []
    for _ in range(n):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2: continue
        a.append(roc_auc_score(y[idx], p[idx]))
    return np.percentile(a, 2.5), np.percentile(a, 97.5)

# ---------- 2. Immune-only model ----------
imm_feats = new_feats + ['age']  # add age as a non-omics covariate
sex_dummy = (df['sex']=='M').astype(int).values
cT_levels = pd.factorize(df['cT'])[0]  # encode T-stage
X_imm = df[imm_feats].apply(pd.to_numeric, errors='coerce').values
X_imm = np.column_stack([X_imm, sex_dummy, cT_levels])
print(f'\n=== immune-only model (P={X_imm.shape[1]}) ===')
res2 = []
for model, pen in [('LASSO','l1'), ('ElasticNet','elasticnet')]:
    probs = nested(X_imm, y, pen, ks=[3, 4, 5])
    auc = roc_auc_score(y, probs)
    lo, hi = boot_ci(y, probs)
    print(f'  {model:10s}  AUC={auc:.3f}  95%CI=[{lo:.3f},{hi:.3f}]')
    res2.append({'model': model, 'feature_set': 'immune_only_5feat',
                 'n_features': X_imm.shape[1],
                 'AUC': round(auc,4), 'CI_low': round(lo,4), 'CI_high': round(hi,4)})
pd.DataFrame(res2).to_csv(f'{OUT}/nested_cv_immune_only.tsv', sep='\t', index=False)

# ---------- 3. CD8 swap model (replace CD8_proliferation with CD8_cytotoxic) ----------
df_swap = df.drop(columns=['CD8_proliferation']).rename(columns={'CD8_cytotoxic':'CD8_proliferation'})
# (This is just a rename trick; we then drop the still-present new sigs to keep total feature count = 37)
df_swap = df_swap.drop(columns=['Tcell_infiltration','Bcell_infiltration'])
DROP_COLS = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']
feats_swap = [c for c in df_swap.columns if c not in DROP_COLS]
X_swap = df_swap[feats_swap].apply(pd.to_numeric, errors='coerce').values
print(f'\n=== CD8 swap model (P={X_swap.shape[1]}, CD8_proliferation := CD8_cytotoxic) ===')
res3 = []
for model, pen in [('LASSO','l1'), ('ElasticNet','elasticnet')]:
    probs = nested(X_swap, y, pen, ks=[5, 8, 12])
    auc = roc_auc_score(y, probs)
    lo, hi = boot_ci(y, probs)
    print(f'  {model:10s}  AUC={auc:.3f}  95%CI=[{lo:.3f},{hi:.3f}]')
    res3.append({'model': model, 'feature_set': 'cd8_swap_37feat',
                 'n_features': X_swap.shape[1],
                 'AUC': round(auc,4), 'CI_low': round(lo,4), 'CI_high': round(hi,4)})
pd.DataFrame(res3).to_csv(f'{OUT}/nested_cv_cd8_swap.tsv', sep='\t', index=False)

print(f'\nWrote univariate_immune_response.tsv, nested_cv_immune_only.tsv, nested_cv_cd8_swap.tsv')
