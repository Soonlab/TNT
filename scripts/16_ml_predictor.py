"""
ML predictor for TNT response using integrated features.
- LOOCV on 35 subjects
- ElasticNet, Random Forest, XGBoost
- Compare univariate to multivariate
"""
import pandas as pd, numpy as np, os, warnings
warnings.filterwarnings('ignore')
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, accuracy_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT='/mnt/sda1/data/TNT/analysis/10_ml_predictor'; os.makedirs(OUT, exist_ok=True)

I = pd.read_csv('/mnt/sda1/data/TNT/analysis/tables/integrated_subject_master.tsv', sep='\t')
y = (I['response_bin']=='good').astype(int).values
print(f'N={len(I)} (good={y.sum()}, bad={(1-y).sum()})')

# Feature set — drop non-numeric/id/label
drop_cols = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']
feats = [c for c in I.columns if c not in drop_cols]
X = I[feats].apply(pd.to_numeric, errors='coerce')
print(f'Features: {len(feats)}')

# Pipeline: impute median + scale + classifier
def make_pipe(clf):
    return Pipeline([('imp',SimpleImputer(strategy='median')),
                     ('sc',StandardScaler()),
                     ('clf',clf)])

models = {
    'ElasticNetLR': make_pipe(LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5, max_iter=10000)),
    'LassoLR': make_pipe(LogisticRegression(penalty='l1', solver='saga', C=0.5, max_iter=10000)),
    'RF': make_pipe(RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4)),
}

loo = LeaveOneOut()
results = []
for name, mdl in models.items():
    try:
        proba = cross_val_predict(mdl, X.values, y, cv=loo, method='predict_proba')[:,1]
        auc = roc_auc_score(y, proba)
        acc = accuracy_score(y, proba>=0.5)
        results.append((name, auc, acc))
        print(f'{name}: LOOCV AUC={auc:.3f}, Acc={acc:.3f}')
    except Exception as e:
        print(f'{name} err: {e}')

res = pd.DataFrame(results, columns=['model','LOOCV_AUC','LOOCV_Accuracy'])
res.to_csv(f'{OUT}/ml_loocv_results.tsv', sep='\t', index=False)

# Feature importance from RF (full-fit)
rf = make_pipe(RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))
rf.fit(X.values, y)
imp = pd.DataFrame({'feature':feats, 'importance':rf.named_steps['clf'].feature_importances_}).sort_values('importance', ascending=False)
imp.to_csv(f'{OUT}/rf_feature_importance.tsv', sep='\t', index=False)
print('\n=== RF top 15 features ===')
print(imp.head(15).to_string(index=False))

# Top-feature-only model (best for small n)
topN = 8
top_feats = imp['feature'].head(topN).tolist()
Xtop = X[top_feats]
for name, mdl in models.items():
    try:
        proba = cross_val_predict(mdl, Xtop.values, y, cv=loo, method='predict_proba')[:,1]
        auc = roc_auc_score(y, proba)
        print(f'{name} (top{topN}): AUC={auc:.3f}')
    except Exception as e:
        print(f'{name} top err: {e}')

# Univariate best feature
from scipy import stats as st
uni=[]
for f in feats:
    x = X[f].values
    mask = ~np.isnan(x)
    if mask.sum()<10: continue
    try:
        u = st.mannwhitneyu(x[mask & (y==1)], x[mask & (y==0)])
        uni.append((f, float(u.pvalue)))
    except: pass
uni = pd.DataFrame(uni, columns=['feature','p']).sort_values('p')
uni.to_csv(f'{OUT}/univariate_p.tsv', sep='\t', index=False)
print('\n=== Univariate top 10 ===')
print(uni.head(10).to_string(index=False))

# Simple plot: LOOCV AUC bar
fig, ax = plt.subplots(figsize=(6,3))
ax.bar(res['model'], res['LOOCV_AUC'], color='steelblue')
ax.axhline(0.5, color='gray', ls='--')
for i, v in enumerate(res['LOOCV_AUC']): ax.text(i, v+0.01, f'{v:.2f}', ha='center')
ax.set_ylabel('LOOCV AUC'); ax.set_ylim(0,1); ax.set_title(f'Response prediction (n={len(y)})')
plt.tight_layout(); plt.savefig(f'{OUT}/Fig_ML_LOOCV_AUC.png', dpi=150, bbox_inches='tight')
print('\nSaved to', OUT)
