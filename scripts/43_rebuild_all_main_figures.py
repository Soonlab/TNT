"""Full rebuild of submission main figures (Fig 1 – Fig 7) from panels_v3 subpanels.

Panel selection rationale:
  Fig 1 Cohort (A sankey · B waterfall · C clinical · D sample matrix · E design)
  Fig 2 WES landscape (A oncoprint · B TMB raincloud · C SBS signature attribution ·
                       D TMB/MSI waterfall · E CNV landscape · F HRD breakdown)
  Fig 3 RNA signatures (A TME radar · B signature heatmap · C volcano ·
                        D forest lollipop · E CD8 biaxial · F TLS Cabrita)
  Fig 4 GSEA (A running ES · B Hallmark bubble · C pathway dotplot ·
              D enrichment map · E ssGSEA heatmap · F category NES box)
  Fig 5 ML predictor (A correlation · B NESTED ROC [v0.7] · C forest · D UMAP · E SHAP · F per-subject)
  Fig 6 Paired cascade (A BCa forest [v0.7] · B slope · C waterfall · D cascade schema)
  Fig 7 External validation (A CD8 meta + Akiyoshi row [v0.7] · B heatmap · C concordance · D discovery-vs-validation)
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PANELS = Path('/data/data/TNT/analysis/figures/panels_v3')
MAIN   = Path('/data/data/TNT/analysis/genome_medicine_submission/main_figures')
V07    = Path('/mnt/sda1/data/TNT/analysis/figures/panels')
MAIN.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({'font.size':9, 'axes.linewidth':0.6})

def place(ax, path, letter=None, title=None):
    if not Path(path).exists():
        ax.text(0.5, 0.5, f'missing: {Path(path).name}', ha='center', va='center',
                transform=ax.transAxes, color='red'); ax.axis('off'); return
    img = mpimg.imread(path)
    ax.imshow(img); ax.axis('off')
    if letter:
        ax.text(-0.02, 1.05, letter, transform=ax.transAxes,
                fontsize=15, fontweight='bold', color='#1F3B5C', va='top')
    if title:
        ax.text(0.5, -0.05, title, transform=ax.transAxes,
                ha='center', va='top', fontsize=8.5, color='#444')

def build(fig_id, layout, panels_list, supertitle, outname):
    nrows, ncols = layout
    fig = plt.figure(figsize=(16, 5.0*nrows))
    gs = fig.add_gridspec(nrows, ncols, hspace=0.12, wspace=0.05)
    for i, p in enumerate(panels_list):
        if p is None: continue
        L, path, title = p
        ax = fig.add_subplot(gs[i//ncols, i%ncols])
        place(ax, path, letter=L, title=title)
    fig.suptitle(supertitle, fontsize=13, fontweight='bold', y=0.995)
    for ext in ('png','pdf'):
        fig.savefig(MAIN/f'{outname}.{ext}', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f'saved {outname}')

# ---------- Fig 1 ----------
build('1', (2,3),
    [('A', PANELS/'Fig1A_sankey.png', 'Patient-level sankey (sex → cT → response)'),
     ('B', PANELS/'Fig1B_waterfall.png', 'Clinical waterfall — per-patient variables'),
     ('C', PANELS/'Fig1C_clinical.png', 'Clinical characteristics by response'),
     ('D', PANELS/'Fig1D_sample_matrix.png', 'Sample matrix (WES × RNA-seq × timepoint)'),
     ('E', PANELS/'Fig1E_design.png', 'Study design — sampling around radiation phase of TNT'),
     None],
    'Figure 1. Cohort overview — 35 MSS LARC patients receiving TNT (good n=18 vs bad n=17).',
    'Fig1_cohort')

# ---------- Fig 2 ----------
build('2', (2,3),
    [('A', PANELS/'Fig2A_oncoprint_journal.png', 'Somatic driver oncoprint (APC, TP53, KRAS, FBXW7, KMT2D)'),
     ('B', PANELS/'Fig2B_TMB_raincloud.png', 'TMB distribution by response — good 1.85 vs bad 1.40 /Mb (P = 0.186)'),
     ('C', PANELS/'Fig2C_signature_attribution.png', 'SBS signature attribution (COSMIC v3.3)'),
     ('D', PANELS/'Fig2D_TMB_MSI_waterfall.png', 'Per-sample TMB and MSI waterfall — all tumors MSS (MSI < 0.19 %)'),
     ('E', PANELS/'Fig2E_CNV_HRD.png', 'CNV landscape and HRD proxy (LST/TAI/LOH)'),
     ('F', PANELS/'Fig2F_HRD_breakdown.png', 'HRD breakdown — LST marginally higher in bad responders (P = 0.037)')],
    'Figure 2. WES landscape — MSS, TMB-low, classic CRC driver distribution.',
    'Fig2_WES_landscape')

# ---------- Fig 3 ----------
build('3', (2,3),
    [('A', PANELS/'Fig3A_TME_radar.png', 'Tumor microenvironment signature radar (good vs bad)'),
     ('B', PANELS/'Fig3B_signature_heatmap.png', '22-signature z-score heatmap'),
     ('C', PANELS/'Fig3C_volcano_journal.png', 'Differential expression volcano (DESeq2)'),
     ('D', PANELS/'Fig3D_forest_lollipop.png', 'Signature-response forest lollipop'),
     ('E', PANELS/'Fig3E_CD8_biaxial.png', 'CD8 effector / exhaustion biaxial'),
     ('F', PANELS/'Fig3F_TLS_Cabrita.png', 'TLS signature (Cabrita) by response')],
    'Figure 3. RNA-seq immune and pathway signatures — pre-CRT.',
    'Fig3_RNA_signatures')

# ---------- Fig 4 ----------
build('4', (2,3),
    [('A', PANELS/'Fig4A_running_ES.png', 'GSEA running enrichment scores (top Hallmark)'),
     ('B', PANELS/'Fig4B_Hallmark_bubble.png', 'Hallmark NES × FDR bubble map'),
     ('C', PANELS/'Fig4C_pathway_dotplot.png', 'Reactome pathway dotplot'),
     ('D', PANELS/'Fig4D_enrichment_map.png', 'GSEA enrichment map network'),
     ('E', PANELS/'Fig4E_ssGSEA_heatmap.png', 'ssGSEA 95-set heatmap by response'),
     ('F', PANELS/'Fig4F_category_NES_box.png', 'Category-level NES box (DNA repair vs EMT vs immune)')],
    'Figure 4. GSEA and pathway integration — pre-CRT good responders enriched in E2F / G2M / DSB-HDR and depleted in EMT.',
    'Fig4_GSEA')

# ---------- Fig 5 — ML predictor (nested ROC [v0.7]) ----------
build('5', (2,3),
    [('A', PANELS/'Fig5A_correlation.png', 'Feature correlation (37 features)'),
     ('B', PANELS/'Fig4B_ROC_nested.png', 'Nested outer-LOOCV ROC [v0.7] — LASSO AUC 0.650 [0.45, 0.83]'),
     ('C', PANELS/'Fig5C_forest_CI.png', 'Feature forest — log OR with 95 % CI'),
     ('D', PANELS/'Fig5D_UMAP.png', 'UMAP of integrated features'),
     ('E', PANELS/'Fig5E_SHAP_beeswarm.png', 'SHAP feature contribution'),
     ('F', PANELS/'Fig5F_per_subject_prediction.png', 'Per-subject predicted probability')],
    'Figure 5. Integrated 37-feature LASSO predictor of final TNT response — nested outer-LOOCV [v0.7, leakage-free].',
    'Fig5_ML_predictor')

# ---------- Fig 6 — Paired cascade (BCa forest [v0.7]) ----------
# wide-top + 3 bottom
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0], hspace=0.12, wspace=0.05)
ax = fig.add_subplot(gs[0, :])
place(ax, V07/'Fig6_cascade_BCa_forest.png', letter='A',
      title='Within- and between-group paired Δ with BCa 95 % CI (v0.7) — only Treg has CI excluding 0')
ax = fig.add_subplot(gs[1, 0])
place(ax, PANELS/'Fig6B_slope_fancy.png', letter='B',
      title='Paired pre → post slopes (Treg, MHC-II, CD8 exhaustion, IGH)')
ax = fig.add_subplot(gs[1, 1])
place(ax, PANELS/'Fig6D_waterfall.png', letter='C',
      title='Per-subject Δ waterfall across cascade features')
ax = fig.add_subplot(gs[1, 2])
place(ax, PANELS/'Fig6F_cascade.png', letter='D',
      title='Radiation-induced cascade schematic (exploratory)')
fig.suptitle('Figure 6. Radiation-induced cascade in eventual good responders (n = 14 paired).',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig6_paired_delta.{ext}', dpi=600, bbox_inches='tight')
plt.close(fig); print('saved Fig6')

# ---------- Fig 7 — External validation (v0.7 CD8 + Akiyoshi) ----------
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.12, wspace=0.05)
ax = fig.add_subplot(gs[0, :])
place(ax, V07/'Fig7_external_CD8_validation_v4.png', letter='A',
      title='9-cohort nCRT meta (N = 721, Stouffer Z = +2.74, P = 0.006) + Akiyoshi 2023 convergent row (N = 298) — > 1,000 patients / 10 cohorts')
ax = fig.add_subplot(gs[1, 0])
place(ax, PANELS/'Fig7B_heatmap.png', letter='B', title='Signature × cohort Δ heatmap')
ax = fig.add_subplot(gs[1, 1])
place(ax, PANELS/'Fig7D_cohort_concordance.png', letter='C', title='Per-cohort concordance with discovery direction')
ax = fig.add_subplot(gs[1, 2])
place(ax, PANELS/'Fig7F_discovery_validation.png', letter='D', title='Discovery vs external effect sizes by signature')
fig.suptitle('Figure 7. External validation — CD8-cytotoxic axis reproducible pan-CRT.',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig7_external_validation.{ext}', dpi=600, bbox_inches='tight')
plt.close(fig); print('saved Fig7')

print('\nAll submission main figures rebuilt in:', MAIN)
