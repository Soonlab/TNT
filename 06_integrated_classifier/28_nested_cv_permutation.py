"""
Task 1 (v0.6 revision): Nested LOOCV + permutation test for the TNT response predictor.

Fixes the leakage in the original LOOCV (which selected features / lambda on the full
cohort). Outer loop = LOOCV (34 train / 1 test). Inner loop = 5-fold CV on the 34
training subjects for BOTH feature pre-selection (SelectKBest, k in {5,8,12}) and
regularisation parameter tuning.

Three model classes:
  - LASSO logistic regression (C grid)
  - ElasticNet logistic regression (C grid, l1_ratio=0.5)
  - RandomForest (top-10 by univariate ANOVA F, n_estimators grid)

For each model we record:
  - held-out outer-CV AUC
  - 95% bootstrap CI from outer-CV held-out probas
  - permutation p-value (n=1000 label shuffles, outer-CV AUC from a lighter
    non-nested but leakage-free pipeline -- nested perm*1000 is infeasible).
    We use a pipeline-wrapped LOOCV (no inner CV) with a fixed modest k=8 and
    default hyper-parameters as the permutation reference; this is standard
    practice for small-N perm tests and still tests the AUC statistic against
    random labels.
  - median number of features selected across the outer folds

Output:
  10_ml_predictor/nested_cv_results.tsv
  10_ml_predictor/nested_outer_probs_<model>.tsv
  figures/panels_v3/Fig5B_ROC_CI_nested.{pdf,png}
"""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RNG = np.random.RandomState(42)
OUT = '/mnt/sda1/data/TNT/analysis/10_ml_predictor'
FIG = '/mnt/sda1/data/TNT/analysis/figures/panels_v3'
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

df = pd.read_csv('/mnt/sda1/data/TNT/analysis/tables/integrated_subject_master.tsv', sep='\t')
y = (df['response_bin'] == 'good').astype(int).values
drop_cols = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']
feats = [c for c in df.columns if c not in drop_cols]
X = df[feats].apply(pd.to_numeric, errors='coerce').values
N, P = X.shape
print(f'N={N} good={y.sum()} bad={(1-y).sum()} features={P}')

# ---------- model definitions ----------
def pipe_lasso():
    return Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc', StandardScaler()),
        ('sel', SelectKBest(score_func=f_classif, k=8)),
        ('clf', LogisticRegression(penalty='l1', solver='saga', max_iter=20000, C=0.5)),
    ])
def pipe_enet():
    return Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc', StandardScaler()),
        ('sel', SelectKBest(score_func=f_classif, k=8)),
        ('clf', LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=20000, C=0.5)),
    ])
def pipe_rf():
    return Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc', StandardScaler()),
        ('sel', SelectKBest(score_func=f_classif, k=10)),
        ('clf', RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4)),
    ])

GRID = {
    'LASSO': (pipe_lasso, {'sel__k':[5,8,12], 'clf__C':[0.1,0.3,1.0,3.0]}),
    'ElasticNet': (pipe_enet, {'sel__k':[5,8,12], 'clf__C':[0.1,0.3,1.0,3.0]}),
    'RandomForest': (pipe_rf, {'sel__k':[8,10,15], 'clf__n_estimators':[300,500]}),
}

# ---------- Nested outer-LOOCV with inner 5-fold tuning ----------
def nested_outer_probs(pipe_fn, grid, Xm, ym):
    loo = LeaveOneOut()
    probs = np.zeros(len(ym)); k_sel = []
    for train_idx, test_idx in loo.split(Xm):
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        gs = GridSearchCV(pipe_fn(), grid, cv=inner, scoring='roc_auc', n_jobs=4, refit=True)
        gs.fit(Xm[train_idx], ym[train_idx])
        p = gs.predict_proba(Xm[test_idx])[:,1][0]
        probs[test_idx[0]] = p
        k_sel.append(gs.best_params_.get('sel__k', np.nan))
    return probs, np.nanmedian(k_sel)

def bootstrap_ci_auc(y_true, probs, n=1000, seed=0):
    rng = np.random.RandomState(seed)
    aucs = []
    n_obs = len(y_true)
    for _ in range(n):
        idx = rng.randint(0, n_obs, n_obs)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], probs[idx]))
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), np.array(aucs)

# ---------- Permutation test: LOOCV (non-nested but leakage-free) ----------
def loocv_auc(pipe_fn, Xm, ym):
    loo = LeaveOneOut()
    p = np.zeros(len(ym))
    for tr, te in loo.split(Xm):
        m = pipe_fn().fit(Xm[tr], ym[tr])
        p[te[0]] = m.predict_proba(Xm[te])[:,1][0]
    try:
        return roc_auc_score(ym, p), p
    except Exception:
        return 0.5, p

def perm_null(pipe_fn, Xm, ym, n_perm=1000, seed=1):
    rng = np.random.RandomState(seed)
    null = np.zeros(n_perm)
    for i in range(n_perm):
        yp = rng.permutation(ym)
        auc, _ = loocv_auc(pipe_fn, Xm, yp)
        null[i] = auc
        if (i+1) % 100 == 0: print(f'  perm {i+1}/{n_perm}')
    return null

# ---------- Run ----------
results = []
roc_data = {}
perm_cache = {}
for name, (fn, grid) in GRID.items():
    print(f'=== {name} nested outer-LOOCV ===')
    probs, kmed = nested_outer_probs(fn, grid, X, y)
    auc = roc_auc_score(y, probs)
    lo, hi, _ = bootstrap_ci_auc(y, probs, n=1000)
    print(f'  outer AUC={auc:.3f}  95%CI=[{lo:.3f},{hi:.3f}]  median_k={kmed}')
    # permutation
    print(f'  permutation null (1000x LOOCV)...')
    null = perm_null(fn, X, y, n_perm=1000, seed=1)
    p_perm = (np.sum(null >= auc) + 1) / (len(null) + 1)
    perm_cache[name] = null
    pd.DataFrame({'subject_id':df['subject_id'].values,'y':y,'prob':probs}).to_csv(
        f'{OUT}/nested_outer_probs_{name}.tsv', sep='\t', index=False)
    roc_data[name] = (probs, auc, lo, hi)
    results.append({'model':name, 'outer_AUC':round(auc,4),
                    '95%CI_low':round(lo,4), '95%CI_high':round(hi,4),
                    'perm_p_value':round(p_perm,4),
                    'n_features_selected_median':int(kmed) if np.isfinite(kmed) else -1})

res = pd.DataFrame(results)
res.to_csv(f'{OUT}/nested_cv_results.tsv', sep='\t', index=False)
print(res)

# ---------- Fig 5B: nested-CV ROC with CI band + permutation inset ----------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})
fig = plt.figure(figsize=(5.4, 4.4))
ax = fig.add_axes([0.13, 0.14, 0.82, 0.78])
palette = {'LASSO':'#1f77b4','ElasticNet':'#2ca02c','RandomForest':'#d62728'}
best = max(roc_data.items(), key=lambda kv: kv[1][1])[0]

for name,(probs,auc,lo,hi) in roc_data.items():
    fpr, tpr, _ = roc_curve(y, probs)
    # bootstrap CI band for best model only
    if name == best:
        rng = np.random.RandomState(0); boot_tprs = []
        grid_fpr = np.linspace(0,1,101)
        for _ in range(500):
            idx = rng.randint(0, len(y), len(y))
            if len(np.unique(y[idx])) < 2: continue
            f_i, t_i, _ = roc_curve(y[idx], probs[idx])
            boot_tprs.append(np.interp(grid_fpr, f_i, t_i))
        boot_tprs = np.array(boot_tprs)
        tpr_lo = np.percentile(boot_tprs, 2.5, axis=0)
        tpr_hi = np.percentile(boot_tprs, 97.5, axis=0)
        ax.fill_between(grid_fpr, tpr_lo, tpr_hi, color=palette[name], alpha=0.15,
                        label=f'{name} 95% CI')
    ax.plot(fpr, tpr, color=palette[name], lw=1.8,
            label=f'{name}  AUC={auc:.2f} [{lo:.2f}, {hi:.2f}]')
ax.plot([0,1],[0,1],'--',color='gray',lw=0.8)
ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
ax.set_title(f'Nested outer-LOOCV ROC (N={N})', fontsize=10)
ax.legend(loc='lower right', frameon=False, fontsize=8)
ax.set_xlim(-0.01,1.01); ax.set_ylim(-0.01,1.02)
for s in ['top','right']: ax.spines[s].set_visible(False)

# inset: permutation null for best model
ax_in = fig.add_axes([0.22, 0.58, 0.28, 0.30])
null_best = perm_cache[best]
ax_in.hist(null_best, bins=30, color='lightgray', edgecolor='gray', lw=0.3)
ax_in.axvline(roc_data[best][1], color=palette[best], lw=1.8)
p_b = (np.sum(null_best >= roc_data[best][1]) + 1) / (len(null_best) + 1)
ax_in.set_title(f'Perm null ({best})\nP={p_b:.3f}', fontsize=7)
ax_in.set_xlabel('LOOCV AUC', fontsize=7); ax_in.set_ylabel('Count', fontsize=7)
ax_in.tick_params(labelsize=6)
for s in ['top','right']: ax_in.spines[s].set_visible(False)

for ext in ('png','pdf'):
    fig.savefig(f'{FIG}/Fig5B_ROC_CI_nested.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'Wrote {FIG}/Fig5B_ROC_CI_nested.pdf/png')
