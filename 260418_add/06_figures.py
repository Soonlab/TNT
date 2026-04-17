"""
Figures for the 260418 immune-feature LASSO analysis.

Fig A — ROC curves for all 5 scenarios (baseline / drop / add / swap / drop+add)
        with bootstrap 95% CI band on the winning model.
Fig B — Bar chart of nested-CV AUC + 95% CI for the 5 scenarios x 2 models.
Fig C — Feature-importance ranking for the winning model (drop_cd8prolif_36 + ElasticNet)
        based on a full-data refit at k=12 most permissive.
Fig D — Univariate Mann-Whitney response association of the new vs existing immune sigs
        (boxplot good vs bad).
"""
import os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/data/data/TNT/analysis/260418_add'
df_all = pd.read_csv(f'{OUT}/integrated_subject_master_v2.tsv', sep='\t')
res = pd.read_csv(f'{OUT}/nested_cv_drop_vs_swap.tsv', sep='\t')
print(res.to_string(index=False))

# ---------------- Fig A : ROC overlay ----------------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8})

scen_files = {
    'baseline_37':           'probs_baseline_37_ElasticNet.tsv',
    'drop_cd8prolif_36':     'probs_drop_cd8prolif_36_ElasticNet.tsv',
    'add_immune_40':         'probs_add_immune_40_ElasticNet.tsv',
    'swap_cd8_37':           'probs_swap_cd8_37_ElasticNet.tsv',
    'drop_prolif_add_3_39':  'probs_drop_prolif_add_3_39_ElasticNet.tsv',
}
scen_label = {
    'baseline_37':          'Baseline (37)',
    'drop_cd8prolif_36':    'Drop CD8_prolif (36) **',
    'add_immune_40':        'Add 3 immune sigs (40)',
    'swap_cd8_37':          'Swap CD8_prolif→cyt (37)',
    'drop_prolif_add_3_39': 'Drop + add (39)',
}
palette = {
    'baseline_37':          '#888888',
    'drop_cd8prolif_36':    '#E63946',
    'add_immune_40':        '#1f77b4',
    'swap_cd8_37':          '#2ca02c',
    'drop_prolif_add_3_39': '#9467bd',
}

fig, ax = plt.subplots(figsize=(5.5, 4.6))
best = 'drop_cd8prolif_36'
for sc, fn in scen_files.items():
    p = pd.read_csv(f'{OUT}/{fn}', sep='\t')
    fpr, tpr, _ = roc_curve(p['y'], p['prob'])
    auc = roc_auc_score(p['y'], p['prob'])
    row = res[(res.scenario==sc) & (res.model=='ElasticNet')].iloc[0]
    lw = 2.0 if sc==best else 1.2
    ls = '-' if sc==best else ('-' if sc in ('baseline_37',) else '--')
    ax.plot(fpr, tpr, color=palette[sc], lw=lw, ls=ls,
            label=f'{scen_label[sc]}  AUC={auc:.3f} [{row.CI_low:.2f}–{row.CI_high:.2f}]')

# bootstrap CI band on best model
p = pd.read_csv(f'{OUT}/{scen_files[best]}', sep='\t')
y, probs = p['y'].values, p['prob'].values
rng = np.random.RandomState(0)
grid_fpr = np.linspace(0,1,101); boot_tprs = []
for _ in range(800):
    idx = rng.randint(0, len(y), len(y))
    if len(np.unique(y[idx])) < 2: continue
    f_i, t_i, _ = roc_curve(y[idx], probs[idx])
    boot_tprs.append(np.interp(grid_fpr, f_i, t_i))
boot_tprs = np.array(boot_tprs)
ax.fill_between(grid_fpr, np.percentile(boot_tprs,2.5,axis=0),
                np.percentile(boot_tprs,97.5,axis=0),
                color=palette[best], alpha=0.15, label=None)

ax.plot([0,1],[0,1],'--',color='lightgray',lw=0.8)
ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
ax.set_title('Nested outer-LOOCV ROC — ElasticNet (N=35)\nimmune-feature ablation', fontsize=10)
ax.legend(loc='lower right', frameon=False, fontsize=7.5)
ax.set_xlim(-0.01,1.01); ax.set_ylim(-0.01,1.02)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/FigA_ROC_5scenarios.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote FigA_ROC_5scenarios')

# ---------------- Fig B : AUC bar chart with CI ----------------
fig, ax = plt.subplots(figsize=(6.4, 4.0))
order = ['baseline_37','drop_cd8prolif_36','add_immune_40','swap_cd8_37','drop_prolif_add_3_39']
xpos = np.arange(len(order))
w = 0.36
for i, model in enumerate(['LASSO','ElasticNet']):
    aucs, los, his = [], [], []
    for sc in order:
        r = res[(res.scenario==sc) & (res.model==model)].iloc[0]
        aucs.append(r.AUC); los.append(r.AUC - r.CI_low); his.append(r.CI_high - r.AUC)
    x = xpos + (i-0.5)*w
    color = '#1f77b4' if model=='LASSO' else '#E63946'
    ax.bar(x, aucs, width=w, color=color, alpha=0.85, label=model, edgecolor='black', lw=0.5)
    ax.errorbar(x, aucs, yerr=[los,his], fmt='none', color='black', lw=0.8, capsize=3)
    for xi, a in zip(x, aucs):
        ax.text(xi, a+0.015, f'{a:.3f}', ha='center', fontsize=7)

ax.axhline(0.5, color='gray', ls=':', lw=0.7)
ax.set_xticks(xpos)
ax.set_xticklabels([scen_label[s].split(' (')[0] for s in order], rotation=15, ha='right', fontsize=8.5)
ax.set_ylabel('Nested outer-LOOCV AUC')
ax.set_ylim(0.30, 0.95)
ax.set_title('AUC across feature-set ablations (95% bootstrap CI, n=2000)', fontsize=10)
ax.legend(frameon=False, loc='upper left')
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/FigB_AUC_bar_5scenarios.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote FigB_AUC_bar_5scenarios')

# ---------------- Fig C : feature importance for winning model ----------------
DROP = ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']
df_w = df_all.drop(columns=['CD8_proliferation','Tcell_infiltration','Bcell_infiltration','CD8_cytotoxic'])
y = (df_w['response_bin']=='good').astype(int).values
feats = [c for c in df_w.columns if c not in DROP]
X = df_w[feats].apply(pd.to_numeric, errors='coerce').values
pipe = Pipeline([('imp', SimpleImputer(strategy='median')),
                 ('sc',  StandardScaler()),
                 ('sel', SelectKBest(score_func=f_classif, k=12)),
                 ('clf', LogisticRegression(penalty='elasticnet', solver='saga',
                                            l1_ratio=0.5, C=0.5, max_iter=20000))])
pipe.fit(X, y)
mask = pipe.named_steps['sel'].get_support()
sel  = [f for f,m in zip(feats, mask) if m]
coefs = pipe.named_steps['clf'].coef_.flatten()
imp = pd.DataFrame({'feature':sel, 'coef':coefs})
imp['abs'] = imp['coef'].abs()
imp = imp.sort_values('abs', ascending=True)
imp_nz = imp[imp['abs'] > 0]
imp_nz.to_csv(f'{OUT}/feature_importance_drop_cd8prolif.tsv', sep='\t', index=False)

fig, ax = plt.subplots(figsize=(6.5, max(2.0, 0.32*len(imp))))
colors = ['#E63946' if c<0 else '#2E86AB' for c in imp['coef']]
ax.barh(imp['feature'], imp['coef'], color=colors, edgecolor='black', lw=0.4)
ax.axvline(0, color='black', lw=0.6)
ax.set_xlabel('Coefficient (good vs bad)')
ax.set_title('ElasticNet coefs — drop_cd8prolif_36 (full-data refit, k=12)', fontsize=10)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/FigC_feature_importance.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote FigC_feature_importance')

# ---------------- Fig D : univariate boxplots good vs bad ----------------
sigs = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration',
        'CD8_proliferation','MHC_II']
uni = pd.read_csv(f'{OUT}/univariate_immune_response.tsv', sep='\t')
uni_dict = uni.set_index('feature')['mw_p_2sided'].to_dict()

fig, axes = plt.subplots(1, len(sigs), figsize=(2.0*len(sigs), 3.4), sharey=False)
for ax, sig in zip(axes, sigs):
    g = df_all[df_all.response_bin=='good'][sig].dropna().values
    b = df_all[df_all.response_bin=='bad'][sig].dropna().values
    parts = ax.boxplot([g, b], positions=[0, 1], widths=0.6, patch_artist=True,
                        medianprops=dict(color='black', lw=1.2),
                        boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7),
                        capprops=dict(lw=0.7), flierprops=dict(marker='o', ms=3))
    for patch, c in zip(parts['boxes'], ['#2E86AB','#E63946']):
        patch.set_facecolor(c); patch.set_alpha(0.55)
    rng = np.random.RandomState(0)
    ax.scatter(rng.normal(0, 0.06, len(g)), g, color='#2E86AB', s=14, edgecolor='black', lw=0.3)
    ax.scatter(rng.normal(1, 0.06, len(b)), b, color='#E63946', s=14, edgecolor='black', lw=0.3)
    ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'])
    ax.set_title(f'{sig}\nMW P={uni_dict.get(sig, float("nan")):.3f}', fontsize=8.5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.suptitle('Univariate response association — discovery cohort (N=33 pre-treatment)', fontsize=10)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{OUT}/FigD_univariate_boxplots.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote FigD_univariate_boxplots')

print('\nAll 4 figures written to', OUT)
