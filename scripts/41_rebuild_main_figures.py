"""Rebuild submission main figures (Fig5, Fig6, Fig7) incorporating v0.7 additions.

Keeps previously-validated subpanels from panels_v3 and patches in:
  - Fig 5: replace ROC panel with nested-LOOCV leakage-free AUC (0.65) result
  - Fig 6: prepend BCa bootstrap-forest cascade summary panel
  - Fig 7: replace with CD8-cytotoxic external meta (9 cohorts, N = 721) + Akiyoshi 2023 convergent row

Output written to genome_medicine_submission/main_figures/ (PDF + PNG at 600 dpi).
Fig 1–Fig 4 left unchanged from the current submission.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PANELS = Path('/data/data/TNT/analysis/figures/panels_v3')
SUPP   = Path('/data/data/TNT/analysis/figures/supp')
MAIN   = Path('/data/data/TNT/analysis/genome_medicine_submission/main_figures')
FIGMAIN= Path('/mnt/sda1/data/TNT/analysis/figures/panels')
MAIN.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({'font.family':'Arial', 'font.size':9, 'axes.linewidth':0.7})

def place_image(ax, path, title=None, letter=None):
    if not Path(path).exists():
        ax.text(0.5, 0.5, f'missing: {Path(path).name}', ha='center', va='center',
                transform=ax.transAxes, color='red')
        ax.axis('off'); return
    img = mpimg.imread(path)
    ax.imshow(img); ax.axis('off')
    if letter:
        ax.text(-0.02, 1.03, letter, transform=ax.transAxes,
                fontsize=14, fontweight='bold', color='#1F3B5C', va='top')
    if title:
        ax.text(0.5, -0.04, title, transform=ax.transAxes,
                ha='center', va='top', fontsize=9, color='#333')

# ============================================================
# Fig 5 — ML predictor (nested-LOOCV honest version)
# ============================================================
# Layout: 2 rows × 3 cols
#   A correlation    B nested ROC     C forest
#   D UMAP           E SHAP           F per-subject prediction
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.08, wspace=0.06)
panels = [
    ('A', PANELS/'Fig5A_correlation.png',         'Feature correlation (37 features)'),
    ('B', Path('/mnt/sda1/data/TNT/analysis/figures/panels_v3/Fig4B_ROC_nested.png'),
                                                  'Nested outer-LOOCV ROC (LASSO AUC 0.65, ElasticNet 0.69)'),
    ('C', PANELS/'Fig5C_forest_CI.png',           'Top-feature forest (log OR with 95 % CI)'),
    ('D', PANELS/'Fig5D_UMAP.png',                'UMAP of integrated features'),
    ('E', PANELS/'Fig5E_SHAP_beeswarm.png',       'SHAP feature contribution'),
    ('F', PANELS/'Fig5F_per_subject_prediction.png', 'Per-subject predicted probability'),
]
for i,(L,p,t) in enumerate(panels):
    ax = fig.add_subplot(gs[i//3, i%3])
    place_image(ax, p, title=t, letter=L)
fig.suptitle('Figure 5. Nested LOOCV LASSO predictor of final TNT response (n = 35).',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig5_ML_predictor.{ext}', dpi=600, bbox_inches='tight')
print('saved Fig5')
plt.close(fig)

# ============================================================
# Fig 6 — Paired delta + BCa forest summary
# ============================================================
# Layout: top row = BCa forest (wide), bottom row = 3 panels (slope, waterfall, cascade)
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0], hspace=0.12, wspace=0.05)
# Panel A: BCa forest across top row (spans 3 columns)
ax = fig.add_subplot(gs[0, :])
place_image(ax, FIGMAIN/'Fig6_cascade_BCa_forest.png',
            title='Within-group (left) and between-group (right) Δ with BCa 95 % CI — only Treg has between-group CI excluding zero',
            letter='A')
# Panel B: Treg/MHC-II slope
ax = fig.add_subplot(gs[1, 0])
place_image(ax, PANELS/'Fig6B_slope_fancy.png',
            title='Paired pre → post slope (Treg, MHC-II, CD8 exhaustion, IGH)',
            letter='B')
# Panel C: Per-feature Δ waterfall
ax = fig.add_subplot(gs[1, 1])
place_image(ax, PANELS/'Fig6D_waterfall.png',
            title='Per-subject Δ waterfall across cascade features',
            letter='C')
# Panel D: Cascade schematic
ax = fig.add_subplot(gs[1, 2])
place_image(ax, PANELS/'Fig6F_cascade.png',
            title='Radiation-induced cascade model (exploratory)',
            letter='D')
fig.suptitle('Figure 6. Radiation-induced cascade in eventual good responders (n = 14 paired).',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig6_paired_delta.{ext}', dpi=600, bbox_inches='tight')
print('saved Fig6')
plt.close(fig)

# ============================================================
# Fig 7 — External validation (CD8 axis + Akiyoshi convergence)
# ============================================================
# Layout: top wide = forest v4 with Akiyoshi row, bottom row = 3 panels
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.12, wspace=0.05)
# Panel A: forest v4 (9-cohort + Akiyoshi) spans all columns
ax = fig.add_subplot(gs[0, :])
place_image(ax, FIGMAIN/'Fig7_external_CD8_validation_v4.png',
            title='9-cohort nCRT meta (N = 721) + Akiyoshi 2023 (N = 298) convergent row — > 1,000 patients / 10 cohorts',
            letter='A')
# Panel B: per-cohort heatmap
ax = fig.add_subplot(gs[1, 0])
place_image(ax, PANELS/'Fig7B_heatmap.png',
            title='Signature × cohort Δ heatmap',
            letter='B')
# Panel C: cohort concordance
ax = fig.add_subplot(gs[1, 1])
place_image(ax, PANELS/'Fig7D_cohort_concordance.png',
            title='Per-cohort concordance with discovery direction',
            letter='C')
# Panel D: discovery vs validation summary
ax = fig.add_subplot(gs[1, 2])
place_image(ax, PANELS/'Fig7F_discovery_validation.png',
            title='Discovery vs external-validation effect sizes by signature',
            letter='D')
fig.suptitle('Figure 7. External validation: CD8-cytotoxic axis reproducible pan-CRT (Stouffer Z = +2.74, P = 0.006).',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig7_external_validation.{ext}', dpi=600, bbox_inches='tight')
print('saved Fig7')
plt.close(fig)

print('\nMain figures rebuilt:')
for f in ('Fig5_ML_predictor','Fig6_paired_delta','Fig7_external_validation'):
    print(f'  {MAIN}/{f}.pdf + .png')
