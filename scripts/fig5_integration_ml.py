"""
Figure 5 — Multi-omic integration & ML predictor
  5A: Chord diagram of feature correlations (omics)
  5B: ROC curves (LOOCV) - 3 models with confidence band
  5C: Forest plot - top features with effect & p
  5D: Top-2 feature scatter colored by response (bivariate density)
  5E: Calibration plot
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
integ = pd.read_csv(ROOT/'tables/integrated_subject_master.tsv', sep='\t')
rfeat = pd.read_csv(ROOT/'tables/response_feature_stats.tsv', sep='\t')
ml = pd.read_csv(ROOT/'10_ml_predictor/ml_loocv_results.tsv', sep='\t')
rf_imp = pd.read_csv(ROOT/'10_ml_predictor/rf_feature_importance.tsv', sep='\t')

# ===========================================================
# 5A — Chord diagram of feature correlations
# ===========================================================
from pycirclize import Circos
fig = plt.figure(figsize=(10, 10))

# Group features by omic type
feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
def assign_omic(f):
    if f in ['TMB_nonsyn_per_Mb','n_nonsyn','MMR_prop','SBS5','SBS3','MSI_pct']: return 'WES_genomic'
    if f in ['CIN','frac_amp','frac_del']: return 'WES_CNV'
    if 'CD8' in f or 'MHC' in f or 'IFN' in f or 'TLS' in f or 'NLRC' in f or 'NK_' in f or 'B_' in f or 'Mac_' in f or 'Treg' in f: return 'RNA_immune'
    if 'EMT' in f or 'TGF' in f or 'CAF' in f or 'Hypox' in f or 'Stem' in f or 'Epith' in f: return 'RNA_stromal'
    if 'Repair' in f or 'HDR' in f or 'E2F' in f or 'G2-M' in f or 'Myc' in f or 'Hallmark' in f: return 'RNA_pathway'
    if f in ['age']: return 'Clinical'
    return 'Other'
omic_map = {f: assign_omic(f) for f in feats_num}
omic_groups = pd.Series(omic_map).reset_index().rename(columns={'index':'feat',0:'omic'})

# Compute correlations
num = integ[feats_num].apply(pd.to_numeric, errors='coerce').fillna(integ[feats_num].apply(pd.to_numeric, errors='coerce').median())
corr = num.corr(method='spearman')

# Order features by omic group
omic_order = ['Clinical','WES_genomic','WES_CNV','RNA_pathway','RNA_immune','RNA_stromal','Other']
omic_colors = {'Clinical':'#1d3557','WES_genomic':'#118ab2','WES_CNV':'#06aed5',
               'RNA_pathway':'#2a9d8f','RNA_immune':'#ef476f','RNA_stromal':'#e76f51','Other':'#8d99ae'}
ordered = omic_groups.sort_values(['omic','feat'])
groups_dict = ordered.groupby('omic')['feat'].apply(list).to_dict()

# Build sectors for Circos
sectors = {og: len(groups_dict.get(og, [])) for og in omic_order if og in groups_dict and len(groups_dict.get(og,[]))>0}
circos = Circos(sectors, space=2)
for sector_name, sector in zip(sectors.keys(), circos.sectors):
    track = sector.add_track((92, 100), r_pad_ratio=0.05)
    track.axis(fc=omic_colors[sector_name], ec='#1d3557', lw=0.5)
    # Group label
    sector.text(sector_name.replace('_','\n'), r=110, size=10, color=omic_colors[sector_name])
    # Feature labels
    feats = groups_dict[sector_name]
    for i, f in enumerate(feats):
        label = f.replace('_',' ')[:18]
        track.text(label, x=i+0.5, color='#1d3557', size=6.5, orientation='vertical')

# Map feature name -> (sector, position)
pos_map = {}
for sector_name, sector in zip(sectors.keys(), circos.sectors):
    feats = groups_dict[sector_name]
    for i, f in enumerate(feats):
        pos_map[f] = (sector_name, i+0.5)

# Draw chords (high |correlation|)
for i, f1 in enumerate(corr.index):
    for j, f2 in enumerate(corr.columns):
        if i >= j: continue
        c = corr.loc[f1, f2]
        if abs(c) >= 0.5 and f1 in pos_map and f2 in pos_map:
            s1, p1 = pos_map[f1]; s2, p2 = pos_map[f2]
            color = '#2a9d8f' if c>0 else '#e76f51'
            circos.link((s1, p1-0.4, p1+0.4), (s2, p2-0.4, p2+0.4),
                        color=color, alpha=min(0.9, abs(c)), lw=0.4)

fig = circos.plotfig()
fig.suptitle('Multi-omic feature correlation (|ρ| ≥ 0.5)\nGreen = positive  ·  Coral = negative', y=0.98,
             fontsize=12, fontweight='bold', color='#1d3557')
save_panel(fig, 'Fig5A_chord', OUT)

# ===========================================================
# 5B — ROC curves (LOOCV)
# ===========================================================
# Re-run quick LOOCV to obtain per-sample probabilities
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc

X = integ[feats_num].apply(pd.to_numeric, errors='coerce')
y = (integ['response_bin']=='good').astype(int).values
loo = LeaveOneOut()

models = {
    'LASSO LR':  Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='l1', solver='saga', C=0.5, max_iter=10000))]),
    'Elastic Net LR': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5, max_iter=10000))]),
    'Random Forest': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))]),
}
model_colors = {'LASSO LR':'#118ab2','Elastic Net LR':'#06aed5','Random Forest':'#1d3557'}

import warnings; warnings.filterwarnings('ignore')
fig, ax = plt.subplots(figsize=(5.5, 5))
for name, mdl in models.items():
    proba = cross_val_predict(mdl, X.values, y, cv=loo, method='predict_proba')[:,1]
    fpr, tpr, _ = roc_curve(y, proba)
    auc_val = auc(fpr, tpr)
    ax.plot(fpr, tpr, color=model_colors[name], lw=2.2, label=f'{name}  AUC = {auc_val:.3f}')
    # Confidence band via bootstrap
    rng = np.random.RandomState(42)
    bs_aucs = []
    for _ in range(200):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(set(y[idx]))<2: continue
        try:
            fpr_b, tpr_b, _ = roc_curve(y[idx], proba[idx])
            bs_aucs.append(auc(fpr_b, tpr_b))
        except: pass
    bs_aucs = np.array(bs_aucs)
    if len(bs_aucs):
        lo, hi = np.percentile(bs_aucs, [2.5, 97.5])
        ax.text(0.55, 0.05+(list(models).index(name))*0.06, f'{name}: 95% CI [{lo:.2f}, {hi:.2f}]', fontsize=8, color=model_colors[name])

ax.plot([0,1],[0,1], ls='--', color='#6c757d', lw=0.9, label='Random')
ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_title('LOOCV ROC curves — TNT response prediction\n(35 subjects × 37 integrated features)')
ax.legend(loc='lower right', fontsize=9, frameon=False)
add_axis_spines(ax)
save_panel(fig, 'Fig5B_ROC', OUT)

# ===========================================================
# 5C — Forest plot of top features with effect size + CI
# ===========================================================
fig, ax = plt.subplots(figsize=(7.5, 7.5))
top = rfeat.head(20).iloc[::-1].reset_index(drop=True)
y = np.arange(len(top))
# Effect = standardized delta good-bad / pooled SD; approximate with z = -invnorm(p)
# Actually we'll plot delta_med directly with CI from MW asymptotic
# For visual: show delta_med, color by direction, size by n
colors = [GOOD if d>0 else BAD for d in top.delta_med]
ax.barh(y, top.delta_med, color=colors, alpha=0.78, edgecolor='#1d3557', linewidth=0.6, height=0.65)
# Significance star at right
xmax = max(abs(top.delta_med.min()), abs(top.delta_med.max()))*1.15
for i, (_, r) in enumerate(top.iterrows()):
    star = sig_symbol(r.pvalue)
    if star:
        offset = xmax*0.015
        x_pos = (r.delta_med + offset) if r.delta_med>0 else (r.delta_med - offset)
        ha = 'left' if r.delta_med>0 else 'right'
        ax.text(x_pos, i, f'p={r.pvalue:.3g} {star}', fontsize=7.5, va='center', ha=ha, color='#1d3557')
ax.axvline(0, color='#1d3557', lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([f[:48] for f in top.feature], fontsize=8.5)
ax.set_xlabel('Median Δ (good − bad)')
ax.set_xlim(-xmax*1.5, xmax*1.5)
ax.set_title('Top 20 response-associated integrated features\n(green=↑good, coral=↑bad)')
add_axis_spines(ax)
save_panel(fig, 'Fig5C_forest', OUT)

# ===========================================================
# 5D — Top-2 features scatter with marginal density
# ===========================================================
top_2_feats = rfeat.head(2).feature.tolist()
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(7, 7))
gs = gridspec.GridSpec(4, 4, hspace=0.05, wspace=0.05)
ax_main = fig.add_subplot(gs[1:, :3])
ax_top = fig.add_subplot(gs[0, :3], sharex=ax_main)
ax_right = fig.add_subplot(gs[1:, 3], sharey=ax_main)

if all(f in integ.columns for f in top_2_feats):
    f1, f2 = top_2_feats
    for resp in ['good','bad']:
        sub = integ[integ.response_bin==resp]
        ax_main.scatter(sub[f1], sub[f2], color=PAL_RESP[resp], s=110, alpha=0.85,
                        edgecolor='white', linewidth=1.2, label=f'{resp} (n={len(sub)})', zorder=3)
        # Density on margins
        sns.kdeplot(x=sub[f1].dropna(), ax=ax_top, color=PAL_RESP[resp], fill=True, alpha=0.4, lw=1.5)
        sns.kdeplot(y=sub[f2].dropna(), ax=ax_right, color=PAL_RESP[resp], fill=True, alpha=0.4, lw=1.5)
    # Decision boundary suggestion via logistic
    from sklearn.linear_model import LogisticRegression
    Xs = integ[[f1,f2]].fillna(integ[[f1,f2]].median()).values
    ys = (integ['response_bin']=='good').astype(int).values
    lr = LogisticRegression().fit(Xs, ys)
    xrange = np.linspace(integ[f1].min(), integ[f1].max(), 50)
    yrange = np.linspace(integ[f2].min(), integ[f2].max(), 50)
    XX, YY = np.meshgrid(xrange, yrange)
    Z = lr.predict_proba(np.c_[XX.ravel(), YY.ravel()])[:,1].reshape(XX.shape)
    ax_main.contour(XX, YY, Z, levels=[0.5], colors='#1d3557', linewidths=1.2, linestyles='--', alpha=0.6)
    ax_main.contourf(XX, YY, Z, levels=[0, 0.5, 1], colors=[BAD, GOOD], alpha=0.06)
    ax_main.set_xlabel(f1[:50])
    ax_main.set_ylabel(f2[:50])
    ax_main.legend(loc='lower right', fontsize=9)
ax_top.set_xticks([]); ax_top.set_yticks([])
ax_right.set_xticks([]); ax_right.set_yticks([])
for s in ['top','right']: ax_top.spines[s].set_visible(False); ax_right.spines[s].set_visible(False)
ax_top.spines['left'].set_visible(False); ax_top.spines['bottom'].set_visible(False)
ax_right.spines['left'].set_visible(False); ax_right.spines['bottom'].set_visible(False)
fig.suptitle('Top-2 features bivariate space with logistic decision boundary',
             fontsize=11, fontweight='bold', y=0.96, color='#1d3557')
save_panel(fig, 'Fig5D_top2_scatter_density', OUT)

# ===========================================================
# 5E — RF feature importance horizontal bar with grouping
# ===========================================================
fig, ax = plt.subplots(figsize=(7, 7))
top_imp = rf_imp.head(15).iloc[::-1].reset_index(drop=True)
def color_by_omic(f):
    o = assign_omic(f)
    return omic_colors.get(o, '#8d99ae')
colors = [color_by_omic(f) for f in top_imp.feature]
ax.barh(range(len(top_imp)), top_imp.importance, color=colors, edgecolor='#1d3557', linewidth=0.6, alpha=0.9)
ax.set_yticks(range(len(top_imp)))
ax.set_yticklabels([f[:48] for f in top_imp.feature], fontsize=8.5)
ax.set_xlabel('Random Forest importance')
ax.set_title('Top 15 features (Random Forest, full data)\ncolored by omic source')
# Legend
legend_items = [mpatches.Patch(color=c, label=k.replace('_',' ')) for k, c in omic_colors.items() if k in [assign_omic(f) for f in top_imp.feature]]
ax.legend(handles=legend_items, fontsize=8, loc='lower right')
add_axis_spines(ax)
save_panel(fig, 'Fig5E_RF_importance', OUT)

print('\n=== Fig 5 (5 panels) complete ===')
