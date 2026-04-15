"""
Figure 5 v3 — journal-style multi-omic integration & ML predictor
Motifs:
  5A: Vasaikar Cell 2018 + Bagaev Cancer Cell 2021 — feature × feature correlation heatmap with omic-group annotation
  5B: Mariathasan Nature 2018 + Litchfield Cell 2021 — ROC curves with bootstrap CI bands
  5C: Mariathasan Nature 2018 — Forest plot with effect-size CI + p-value annotations
  5D: Capper Nature 2018 + Combes Cell 2022 — UMAP/PCA of subjects with response coloring + KDE density
  5E: Lundberg Nat Mach Intell 2020 + Krishna Cancer Cell 2024 — SHAP-style beeswarm of feature importance
  5F: Capper Nature 2018 + Bagaev — per-subject prediction confidence waterfall

Style: no titles, deep colors.
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch, Ellipse
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import warnings; warnings.filterwarnings('ignore')

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}

OMIC_COLORS = {
    'Clinical':'#0e2a47','WES_genomic':'#118ab2','WES_CNV':'#06aed5',
    'RNA_pathway':'#057a64','RNA_immune':'#c11456','RNA_stromal':'#b03219','Other':'#5a6772',
}

def assign_omic(f):
    if f in ['TMB_nonsyn_per_Mb','n_nonsyn','MMR_prop','SBS5','SBS3','MSI_pct']: return 'WES_genomic'
    if f in ['CIN','frac_amp','frac_del']: return 'WES_CNV'
    if 'CD8' in f or 'MHC' in f or 'IFN' in f or 'TLS' in f or 'NLRC' in f or 'NK_' in f or 'B_' in f or 'Mac_' in f or 'Treg' in f: return 'RNA_immune'
    if 'EMT' in f or 'TGF' in f or 'CAF' in f or 'Hypox' in f or 'Stem' in f or 'Epith' in f: return 'RNA_stromal'
    if 'Repair' in f or 'HDR' in f or 'E2F' in f or 'G2-M' in f or 'Myc' in f or 'Hallmark' in f: return 'RNA_pathway'
    if f in ['age']: return 'Clinical'
    return 'Other'

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'
clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
integ = pd.read_csv(ROOT/'tables/integrated_subject_master.tsv', sep='\t')
rfeat = pd.read_csv(ROOT/'tables/response_feature_stats.tsv', sep='\t')
ml = pd.read_csv(ROOT/'10_ml_predictor/ml_loocv_results.tsv', sep='\t')
rf_imp = pd.read_csv(ROOT/'10_ml_predictor/rf_feature_importance.tsv', sep='\t')

feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
omic_map = {f: assign_omic(f) for f in feats_num}

X_raw = integ[feats_num].apply(pd.to_numeric, errors='coerce')
X_imputed = X_raw.fillna(X_raw.median())
y = (integ['response_bin']=='good').astype(int).values

# ============================================================
# 5A — Multi-omic feature correlation heatmap (Vasaikar/Bagaev style)
# ============================================================
from scipy.cluster.hierarchy import linkage, leaves_list

corr = X_imputed.corr(method='spearman')
# Cluster within each omic group
omic_order = ['Clinical','WES_genomic','WES_CNV','RNA_pathway','RNA_immune','RNA_stromal','Other']
ordered_features = []
for og in omic_order:
    in_og = [f for f in feats_num if omic_map[f]==og]
    if len(in_og) > 1:
        sub_corr = corr.loc[in_og, in_og].fillna(0)
        Z = linkage(sub_corr.values, method='average')
        idx = leaves_list(Z)
        ordered_features.extend([in_og[i] for i in idx])
    else:
        ordered_features.extend(in_og)

corr_o = corr.loc[ordered_features, ordered_features]

# Build figure
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(3, 4,
    height_ratios=[0.16, 8, 1],
    width_ratios=[0.16, 8, 0.5, 1.6],
    hspace=0.04, wspace=0.04)

# Top color bar (omic group)
ax_top = fig.add_subplot(gs[0, 1])
top_colors = [OMIC_COLORS[omic_map[f]] for f in ordered_features]
arr_top = np.array([[to_rgb(c) for c in top_colors]])
ax_top.imshow(arr_top, aspect='auto', interpolation='nearest', extent=[0, len(ordered_features), 0, 1])
ax_top.set_xticks([]); ax_top.set_yticks([])
for s in ['top','right','left','bottom']: ax_top.spines[s].set_visible(False)

# Left color bar
ax_left = fig.add_subplot(gs[1, 0])
left_arr = np.array([[to_rgb(c) for c in top_colors]])
ax_left.imshow(left_arr.transpose(1,0,2), aspect='auto', interpolation='nearest',
               extent=[0, 1, 0, len(ordered_features)])
ax_left.set_xticks([]); ax_left.set_yticks([])
for s in ['top','right','left','bottom']: ax_left.spines[s].set_visible(False)
ax_left.invert_yaxis()

# Main heatmap
ax_h = fig.add_subplot(gs[1, 1])
im = ax_h.imshow(corr_o.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto', interpolation='nearest',
                 extent=[0, len(ordered_features), 0, len(ordered_features)])
ax_h.set_xticks(np.arange(len(ordered_features))+0.5)
ax_h.set_xticklabels([f[:25] for f in ordered_features], rotation=90, fontsize=6.5, color='#0e2a47')
ax_h.set_yticks(np.arange(len(ordered_features))+0.5)
ax_h.set_yticklabels([f[:25] for f in ordered_features[::-1]], fontsize=6.5, color='#0e2a47')
ax_h.invert_yaxis()
ax_h.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# Group dividers
group_starts = []
cur = 0
for og in omic_order:
    in_og = [f for f in ordered_features if omic_map[f]==og]
    if len(in_og)==0: continue
    group_starts.append(cur)
    cur += len(in_og)
for gs_pos in group_starts[1:]:
    ax_h.axvline(gs_pos, color='white', lw=2)
    ax_h.axhline(len(ordered_features)-gs_pos, color='white', lw=2)

# Right col: colorbar (top half) + omic legend (bottom half)
cax = fig.add_subplot(gs[1, 2])
# Inset colorbar: small at top
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cbar_ax = inset_axes(cax, width="80%", height="40%", loc='upper center',
                     bbox_to_anchor=(0, 0, 1, 1), bbox_transform=cax.transAxes)
cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
cb.set_label('Spearman ρ', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)
cax.axis('off')

# Omic legend (right column)
ax_leg = fig.add_subplot(gs[1, 3])
ax_leg.axis('off')
ax_leg.text(0.0, 0.99, 'Omic group', fontsize=10, color='#0e2a47', fontweight='bold',
            ha='left', va='top', transform=ax_leg.transAxes)
omic_present = [og for og in omic_order if any(omic_map[f]==og for f in ordered_features)]
y_step = 0.075
for i, og in enumerate(omic_present):
    y_pos = 0.92 - i*y_step
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.025), 0.18, 0.05, color=OMIC_COLORS[og],
                                transform=ax_leg.transAxes))
    ax_leg.text(0.27, y_pos, og.replace('_',' '), transform=ax_leg.transAxes,
                fontsize=9, va='center', color='#0e2a47')

save_panel(fig, 'Fig5A_correlation', OUT)

# ============================================================
# 5B — ROC curves with bootstrap CI bands (Mariathasan/Litchfield style)
# ============================================================
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc

models = {
    'LASSO LR':  Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='l1', solver='saga', C=0.5, max_iter=10000))]),
    'Elastic Net LR': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5, max_iter=10000))]),
    'Random Forest': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))]),
}
model_colors = {'LASSO LR':'#118ab2','Elastic Net LR':'#06aed5','Random Forest':'#7a3aad'}

X = X_raw.values
loo = LeaveOneOut()

# Compute predictions + bootstrap CI per model
n_boot = 200
fig, ax = plt.subplots(figsize=(6.5, 6))
np.random.seed(42)

mean_fpr = np.linspace(0, 1, 100)

for name, mdl in models.items():
    proba = cross_val_predict(mdl, X, y, cv=loo, method='predict_proba')[:, 1]
    fpr, tpr, _ = roc_curve(y, proba)
    auc_val = auc(fpr, tpr)
    # Bootstrap TPR at fixed FPR grid
    bs_tprs = []
    bs_aucs = []
    for _ in range(n_boot):
        idx = np.random.choice(len(y), len(y), replace=True)
        if len(set(y[idx])) < 2: continue
        try:
            f_b, t_b, _ = roc_curve(y[idx], proba[idx])
            interp_t = np.interp(mean_fpr, f_b, t_b)
            interp_t[0] = 0
            bs_tprs.append(interp_t)
            bs_aucs.append(auc(f_b, t_b))
        except: pass
    bs_tprs = np.array(bs_tprs)
    if len(bs_tprs) > 0:
        tpr_lo = np.percentile(bs_tprs, 2.5, axis=0)
        tpr_hi = np.percentile(bs_tprs, 97.5, axis=0)
        ax.fill_between(mean_fpr, tpr_lo, tpr_hi, color=model_colors[name], alpha=0.18)
    auc_lo, auc_hi = np.percentile(bs_aucs, [2.5, 97.5])
    ax.plot(fpr, tpr, color=model_colors[name], lw=2.6,
            label=f'{name}\n  AUC = {auc_val:.3f}  [95% CI: {auc_lo:.2f}–{auc_hi:.2f}]')

ax.plot([0,1],[0,1], ls='--', color='#5a6772', lw=1.0, label='Random (AUC = 0.5)')
ax.set_xlabel('False positive rate (1 − specificity)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('True positive rate (sensitivity)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.legend(loc='lower right', fontsize=9, frameon=True, framealpha=0.95, edgecolor='#0e2a47')
add_axis_spines(ax)
save_panel(fig, 'Fig5B_ROC_CI', OUT)

# ============================================================
# 5C — Forest plot with CI (Mariathasan Nature 2018 style)
# ============================================================
fig, ax = plt.subplots(figsize=(9, 8))
top = rfeat.head(20).iloc[::-1].reset_index(drop=True)

# Bootstrap CI per feature
np.random.seed(42)
def bootstrap_delta(g_vals, b_vals, n=1000):
    d = []
    for _ in range(n):
        gi = np.random.choice(g_vals, len(g_vals), replace=True)
        bi = np.random.choice(b_vals, len(b_vals), replace=True)
        d.append(np.median(gi) - np.median(bi))
    return np.percentile(d, [2.5, 97.5])

cis = []
for _, r in top.iterrows():
    if r.feature not in integ.columns:
        cis.append((0,0)); continue
    vals = pd.to_numeric(integ[r.feature], errors='coerce')
    g_vals = vals[integ.response_bin=='good'].dropna().values
    b_vals = vals[integ.response_bin=='bad'].dropna().values
    if len(g_vals)>=3 and len(b_vals)>=3:
        cis.append(bootstrap_delta(g_vals, b_vals))
    else:
        cis.append((0,0))
top['ci_low'] = [c[0] for c in cis]; top['ci_high'] = [c[1] for c in cis]

y_pos = np.arange(len(top))
xmax_data = max(top.ci_high.max(), top.delta_med.max())
xmin_data = min(top.ci_low.min(), top.delta_med.min())
right_margin = (xmax_data - xmin_data) * 0.55

# Effect size CI bars
for i, (_, r) in enumerate(top.iterrows()):
    color = GOOD_DEEP if r.delta_med>0 else BAD_DEEP
    ax.plot([r.ci_low, r.ci_high], [i, i], color=color, lw=2.2, alpha=0.85, solid_capstyle='round')
    ax.scatter(r.delta_med, i, s=180 if r.pvalue<0.05 else 100,
               color=color, edgecolor='white', linewidth=1.4, zorder=3)
    star = sig_symbol(r.pvalue)
    if star == 'ns': star = ''
    label = f'p = {r.pvalue:.3g} {star}'
    ax.text(xmax_data + right_margin*0.04, i, label, va='center', ha='left',
            fontsize=9, color='#0e2a47')

ax.axvline(0, color='#0e2a47', lw=1.0)
ax.axvspan(-0.03, 0.03, color='#dee2e6', alpha=0.4)
ax.set_yticks(y_pos)
ax.set_yticklabels([f[:50] for f in top.feature], fontsize=9.5, color='#0e2a47')
ax.set_xlabel('Median Δ (good − poor),  95% bootstrap CI',
              fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(xmin_data*1.1, xmax_data + right_margin)
add_axis_spines(ax)

# Direction guide top-right
ax.text(0.98, 0.98, '↑ Good      ↑ Poor →',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        color='#0e2a47', fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.4'))

save_panel(fig, 'Fig5C_forest_CI', OUT)

# ============================================================
# 5D — Dimensionality reduction (UMAP) of subjects (Capper / Combes)
# ============================================================
from sklearn.decomposition import PCA
import umap

X_scaled = StandardScaler().fit_transform(X_imputed.values)

# Try UMAP, fall back to PCA
try:
    reducer = umap.UMAP(n_neighbors=8, min_dist=0.3, random_state=42)
    embedding = reducer.fit_transform(X_scaled)
    method_label = 'UMAP'
except Exception as e:
    pca = PCA(n_components=2)
    embedding = pca.fit_transform(X_scaled)
    method_label = 'PCA (PC1 vs PC2)'

fig, ax = plt.subplots(figsize=(7.5, 6.5))

# KDE density background
for resp in ['good','bad']:
    mask = (integ.response_bin == resp).values
    try:
        sns.kdeplot(x=embedding[mask, 0], y=embedding[mask, 1], ax=ax,
                    color=PAL[resp], levels=4, alpha=0.4, linewidths=1.2)
    except: pass

# Scatter
for resp in ['good','bad']:
    mask = (integ.response_bin == resp).values
    ax.scatter(embedding[mask, 0], embedding[mask, 1],
               color=PAL[resp], s=140, alpha=0.92,
               edgecolor='white', linewidth=1.4, label=f'{resp} (n={mask.sum()})', zorder=4)

# Annotate each subject
for i in range(len(integ)):
    sid = integ.iloc[i].subject_id
    ax.annotate(f'S{int(sid)}', (embedding[i, 0], embedding[i, 1]),
                xytext=(5, 5), textcoords='offset points', fontsize=7,
                color='#0e2a47', alpha=0.6)

ax.set_xlabel(f'{method_label} dim 1', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel(f'{method_label} dim 2', fontsize=11, fontweight='bold', color='#0e2a47')
ax.legend(loc='upper right', fontsize=10, frameon=True, framealpha=0.92, edgecolor='#0e2a47')
add_axis_spines(ax)
save_panel(fig, 'Fig5D_UMAP', OUT)

# ============================================================
# 5E — SHAP-style beeswarm (Lundberg / Krishna Cancer Cell 2024 style)
# ============================================================
import shap
# Use Random Forest for interpretability
rf = Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
               ('clf',RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))])
rf.fit(X, y)
# Get processed X
X_proc = rf.named_steps['sc'].transform(rf.named_steps['imp'].transform(X))
explainer = shap.TreeExplainer(rf.named_steps['clf'])
shap_values = explainer.shap_values(X_proc)
# RandomForestClassifier returns shape (n_samples, n_features, n_classes)
if isinstance(shap_values, list):
    sv = shap_values[1]  # class 1 (good)
elif shap_values.ndim == 3:
    sv = shap_values[:, :, 1]
else:
    sv = shap_values

# Top features by mean(|SHAP|)
mean_abs = np.abs(sv).mean(axis=0)
top_idx = np.argsort(mean_abs)[-15:]  # top 15 features
top_feat_names = [feats_num[i] for i in top_idx]

fig, ax = plt.subplots(figsize=(9, 7))
np.random.seed(0)
for plot_i, fi in enumerate(top_idx):
    fname = feats_num[fi]
    sv_f = sv[:, fi]
    val_f = X_proc[:, fi]
    # Color by feature value (high=red, low=blue) — classic SHAP style
    cmap_shap = LinearSegmentedColormap.from_list('shap', ['#118ab2','#dde2e8','#c11456'])
    val_norm = (val_f - val_f.min()) / (val_f.max() - val_f.min() + 1e-9)
    colors = cmap_shap(val_norm)
    # Jitter on y axis
    jitter_y = plot_i + np.random.uniform(-0.3, 0.3, len(sv_f))
    ax.scatter(sv_f, jitter_y, s=45, c=colors, alpha=0.85, edgecolor='#0e2a47', linewidth=0.4, zorder=3)

ax.axvline(0, color='#0e2a47', lw=0.9)
ax.set_yticks(range(len(top_idx)))
ax.set_yticklabels([f[:48] for f in top_feat_names], fontsize=9.5, color='#0e2a47')
ax.set_xlabel('SHAP value (impact on predicted "good responder" probability)',
              fontsize=10.5, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)

# Color bar
from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mcolors
norm = mcolors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap_shap, norm=norm)
sm.set_array([])
cax_sb = fig.add_axes([0.97, 0.3, 0.012, 0.4])
cb = fig.colorbar(sm, cax=cax_sb, ticks=[0, 1])
cb.ax.set_yticklabels(['low', 'high'], fontsize=8.5, color='#0e2a47')
cb.set_label('Feature value\n(scaled)', fontsize=9.5, color='#0e2a47', fontweight='bold')
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)

save_panel(fig, 'Fig5E_SHAP_beeswarm', OUT)

# ============================================================
# 5F — Per-subject prediction confidence waterfall (Capper Nature 2018)
# ============================================================
# Use LASSO LR (best AUC) for per-subject prediction
mdl_best = models['LASSO LR']
proba_best = cross_val_predict(mdl_best, X, y, cv=loo, method='predict_proba')[:, 1]
pred_df = pd.DataFrame({'subject_id': integ.subject_id, 'true':y, 'true_label': integ.response_bin,
                         'prob_good': proba_best})
pred_df['correct'] = (pred_df.prob_good >= 0.5) == pred_df.true.astype(bool)
# Sort by predicted prob descending
pred_df = pred_df.sort_values('prob_good', ascending=False).reset_index(drop=True)
# Center bar around 0.5
pred_df['height'] = pred_df.prob_good - 0.5

fig, ax = plt.subplots(figsize=(13, 5))
x_pos = np.arange(len(pred_df))
colors = [PAL[lbl] for lbl in pred_df.true_label]
bars = ax.bar(x_pos, pred_df.height, color=colors, edgecolor='white', linewidth=0.8, width=0.85)

ax.axhline(0, color='#0e2a47', lw=1.0)
ax.set_ylim(-0.55, 0.55)

# Y axis: show true probability
ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
ax.set_yticklabels(['0\n(predicted poor)', '0.25', '0.5\n(decision)', '0.75', '1.0\n(predicted good)'],
                   fontsize=8.5, color='#0e2a47')

# Subject labels
ax.set_xticks(x_pos)
ax.set_xticklabels([f'S{int(s)}' for s in pred_df.subject_id], fontsize=7, rotation=90, color='#0e2a47')
ax.set_xlim(-0.6, len(pred_df)-0.4)

# Mark misclassified with × marker
mis = pred_df[~pred_df.correct]
for _, r in mis.iterrows():
    i = pred_df.index[pred_df.subject_id==r.subject_id][0]
    h = r.height
    y_marker = h + (0.04 if h>=0 else -0.04)
    ax.scatter(i, y_marker, marker='x', s=70, color='#d62828', linewidth=2.2, zorder=5)

# Legend
handles = [
    mpatches.Patch(color=GOOD_DEEP, label='True good (n={})'.format(int(pred_df.true.sum()))),
    mpatches.Patch(color=BAD_DEEP, label='True poor (n={})'.format(int((1-pred_df.true).sum()))),
    matplotlib.lines.Line2D([0],[0], marker='x', color='#d62828', markersize=10, lw=0,
                             label='Misclassified (n={})'.format(int((~pred_df.correct).sum()))),
]
ax.legend(handles=handles, loc='upper right', fontsize=9.5, frameon=True, framealpha=0.95,
          edgecolor='#0e2a47')

ax.set_xlabel('Subject (ranked by predicted good-responder probability)',
              fontsize=10.5, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)

# Accuracy callout
acc = pred_df.correct.mean()
ax.text(0.01, 0.96, f'LASSO LR LOOCV accuracy = {acc:.1%}',
        transform=ax.transAxes, ha='left', va='top', fontsize=10, color='#0e2a47',
        fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.4'))

save_panel(fig, 'Fig5F_per_subject_prediction', OUT)

print('\n=== Fig 5 v3 (6 journal-style panels) saved ===')
