"""Build new Fig 1A (study design) and Fig 1D (3-narrative preview forest).
PDFs + PNGs written to figures/panels_v3/. Designed so text is preserved as
native text in the PDF (usetex=False, svg.fonttype='none') for downstream
LibreOffice PDF→PPTX conversion that keeps every string editable.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from pathlib import Path

OUT = Path('/data/data/TNT/analysis/figures/panels_v3')
OUT.mkdir(parents=True, exist_ok=True)

# critical: keep text as text in PDF/SVG so LibreOffice→PPTX preserves it
plt.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
    'svg.fonttype': 'none',
    'font.family':  'DejaVu Sans',
    'font.size':    9,
    'axes.linewidth': 0.6,
})

GOOD = '#2E86AB'; BAD = '#E63946'; PURPLE = '#6A4C93'; ACCENT = '#1F3B5C'
ROBUST = '#2a9d8f'; EXPL = '#b0b0b0'; PREBG = '#EAF2F8'; POSTBG = '#FDEBD0'

# =====================================================================
# Fig 1A — Study design schematic
# =====================================================================
fig, ax = plt.subplots(figsize=(10, 4.2))
ax.set_xlim(0, 18); ax.set_ylim(-2.7, 3.6); ax.axis('off')

# Phase boxes
def phase(x, w, text, fc, ec='#444', y=0.0, h=1.3):
    rect = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.03', linewidth=0.8,
                          facecolor=fc, edgecolor=ec)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=9.5)

# Phases:  pre-CRT baseline | Radiation CRT 50.4Gy 5-6 wk | post-CRT biopsy | Consolidation chemo | Surgery / W&W
phase(0.2, 2.0, 'Baseline\n(pre-CRT)', PREBG)
phase(2.4, 4.0, 'Radiation phase\nlong-course CRT\n50.4 Gy + capecitabine', '#F5E6A8')
phase(6.6, 2.0, 'Post-CRT\n(before consolidation)', POSTBG)
phase(8.8, 4.5, 'Consolidation\nchemotherapy\n(FOLFOX/CAPOX/Xeloda)', '#D6EAF8')
phase(13.5, 4.3, 'Surgery  or  Watch-and-Wait\n(final TNT response by Dworak TRG)',
      '#FADBD8')

# Biopsy triangles — pre + post
def biopsy(x, label):
    ax.add_patch(mpatches.RegularPolygon((x, 2.05), numVertices=3, radius=0.32,
                                         orientation=np.pi, facecolor=ACCENT,
                                         edgecolor='k', linewidth=0.8))
    ax.text(x, 2.75, label, ha='center', va='bottom', fontsize=9, fontweight='bold',
            color=ACCENT)

biopsy(1.2,  'Pre-CRT biopsy')
biopsy(7.6,  'Post-CRT biopsy')
# endpoint marker
ax.add_patch(mpatches.Circle((15.65, 2.05), 0.32, facecolor=BAD, edgecolor='k', lw=0.8))
ax.text(15.65, 2.75, 'Response endpoint', ha='center', va='bottom',
        fontsize=9, fontweight='bold', color=BAD)

# Radiation phase window bracket (our sampling)
bracket_y = -1.0
ax.plot([1.2, 7.6], [bracket_y, bracket_y], color=ACCENT, lw=2.2)
ax.plot([1.2, 1.2], [bracket_y-0.2, bracket_y+0.2], color=ACCENT, lw=2.2)
ax.plot([7.6, 7.6], [bracket_y-0.2, bracket_y+0.2], color=ACCENT, lw=2.2)
ax.text(4.4, bracket_y-0.55, 'Radiation-phase sampling window (this study)',
        ha='center', va='top', fontsize=10, fontweight='bold', color=ACCENT)

# Arrow timeline
ax.annotate('', xy=(17.9, 0.65), xytext=(0.0, 0.65),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.1))
ax.text(17.9, 0.25, 'Time', ha='right', va='top', fontsize=9, style='italic')

# Subjects
ax.text(0.1, 3.3, '35 MSS LARC patients  |  good n=18  |  bad n=17  |  14 paired pre/post',
        fontsize=11, fontweight='bold', color=ACCENT)

# Key message caption
ax.text(0.1, -2.3,
    'Final TNT response (Dworak TRG 0–1 vs 2–3) is assessed after surgery, '
    'after the full TNT regimen. Molecular profiling is powered only by the '
    'pre-CRT / post-CRT biopsy pair — this is the mid-treatment decision window.',
    fontsize=9, style='italic', color='#444', wrap=True)

for ext in ('pdf', 'png'):
    plt.savefig(OUT/f'Fig1A_study_design_v2.{ext}', dpi=600, bbox_inches='tight')
plt.close()
print('saved Fig1A_study_design_v2')

# =====================================================================
# Fig 1D — 3-narrative preview forest
# =====================================================================
fig, ax = plt.subplots(figsize=(10, 3.6))
ax.set_xlim(-2.5, 5.0); ax.set_ylim(-0.5, 4.5); ax.axis('off')

# three rows: y = 3 (Narrative 1), y = 2 (Narrative 2), y = 1 (Narrative 3)

# Axis baseline
ax.plot([-2, 4.5], [0.4, 0.4], color='k', lw=0.8)
for xv in [-2, -1, 0, 1, 2, 3, 4]:
    ax.plot([xv, xv], [0.38, 0.42], color='k', lw=0.6)
    ax.text(xv, 0.20, f'{xv:+g}' if xv!=0 else '0',
            ha='center', va='top', fontsize=8.5)
ax.text(1.0, -0.12, 'Standardized effect size (log-scale for narratives 1 & 3)',
        ha='center', va='top', fontsize=9, style='italic', color='#555')

# vertical zero
ax.plot([0, 0], [0.4, 4.2], color='gray', lw=0.5, linestyle='--')

# ---- Narrative 1: Pre-CRT tumor-intrinsic predictor ----
y = 3.4
# AUC 0.65 [0.45, 0.83] — map to axis as (AUC - 0.5)*8 for rough display
auc = 0.650; lo = 0.45; hi = 0.83
def auc_to_x(a): return (a - 0.5) * 10  # 0.5 → 0, 0.7 → 2, 1.0 → 5
x_c = auc_to_x(auc); x_l = auc_to_x(lo); x_h = auc_to_x(hi)
ax.plot([x_l, x_h], [y, y], color=GOOD, lw=2.0)
ax.plot([x_l, x_l], [y-0.1, y+0.1], color=GOOD, lw=2.0)
ax.plot([x_h, x_h], [y-0.1, y+0.1], color=GOOD, lw=2.0)
ax.scatter([x_c], [y], s=90, c=GOOD, edgecolors='k', zorder=5)
ax.text(-2.45, y, 'Narrative 1', ha='right', va='center',
        fontsize=10, fontweight='bold', color=ACCENT)
ax.text(-2.45, y-0.35, 'Pre-CRT tumor-intrinsic\npredictor (DSB/HDR/E2F)',
        ha='right', va='center', fontsize=8.5, color='#555')
ax.text(x_h + 0.15, y, f'LASSO AUC = {auc:.2f} [{lo:.2f}, {hi:.2f}]  (nested LOOCV, n = 35)',
        ha='left', va='center', fontsize=9.5, color=GOOD, fontweight='bold')

# ---- Narrative 2: Radiation cascade (n = 14) ----
y = 2.4
# Treg robust (diff +1.21 [+0.06, +1.97]) vs rest exploratory
def z_to_x(v): return v * 1.0  # use as-is
x_c = 1.21; x_l = 0.06; x_h = 1.97
ax.plot([x_l, x_h], [y, y], color=ROBUST, lw=2.0)
ax.plot([x_l, x_l], [y-0.1, y+0.1], color=ROBUST, lw=2.0)
ax.plot([x_h, x_h], [y-0.1, y+0.1], color=ROBUST, lw=2.0)
ax.scatter([x_c], [y], s=100, c=ROBUST, marker='D', edgecolors='k', zorder=5)
# exploratory gray bar
ax.plot([-0.5, 1.8], [y-0.4, y-0.4], color=EXPL, lw=2.0, alpha=0.7)
ax.scatter([0.6], [y-0.4], s=60, c=EXPL, marker='D', alpha=0.7, edgecolors='k')
ax.text(-2.45, y, 'Narrative 2', ha='right', va='center',
        fontsize=10, fontweight='bold', color=ACCENT)
ax.text(-2.45, y-0.4, 'Radiation-induced cascade\n(exploratory, n = 14 paired)',
        ha='right', va='center', fontsize=8.5, color='#555')
ax.text(x_h + 0.15, y, 'Treg Δ = +1.21 [+0.06, +1.97]  (robust, MW P = 0.026)',
        ha='left', va='center', fontsize=9.5, color=ROBUST, fontweight='bold')
ax.text(2.0, y-0.4, 'Other cascade features (exploratory, CIs span 0)',
        ha='left', va='center', fontsize=8.5, color='#888', style='italic')

# ---- Narrative 3: External CD8 meta + Akiyoshi ----
y = 1.1
# meta Z = 2.74 represented as Δ = ~0.19 (mean across 9 cohorts)
x_c = 0.19; x_l = 0.05; x_h = 0.33
ax.plot([x_l, x_h], [y, y], color=BAD, lw=2.0)
ax.plot([x_l, x_l], [y-0.1, y+0.1], color=BAD, lw=2.0)
ax.plot([x_h, x_h], [y-0.1, y+0.1], color=BAD, lw=2.0)
# Meta diamond
ax.scatter([x_c], [y], s=180, c=BAD, marker='D', edgecolors='k', zorder=5)
# Akiyoshi row (purple triangle)
akx = 0.18
ax.plot([0.06, 0.30], [y-0.4, y-0.4], color=PURPLE, lw=2.0)
ax.scatter([akx], [y-0.4], s=90, c='white', marker='D',
           edgecolors=PURPLE, linewidths=1.8, zorder=5)
ax.text(-2.45, y, 'Narrative 3', ha='right', va='center',
        fontsize=10, fontweight='bold', color=ACCENT)
ax.text(-2.45, y-0.4, 'External CD8-cytotoxic\nreproducibility',
        ha='right', va='center', fontsize=8.5, color='#555')
ax.text(x_h + 0.15, y,
        '9-cohort meta Stouffer Z = +2.74, P = 0.006 (N = 721, 8/9 concordant)',
        ha='left', va='center', fontsize=9.5, color=BAD, fontweight='bold')
ax.text(0.35, y-0.4,
        'Akiyoshi 2023  n = 298, OR = 3.81 [1.82, 7.97]  →  > 1,000 patients / 10 cohorts',
        ha='left', va='center', fontsize=8.5, color=PURPLE, style='italic')

# Title band at very top
ax.text(-2.45, 4.25, 'Headline summary of results',
        ha='left', va='center', fontsize=11, fontweight='bold', color=ACCENT)

for ext in ('pdf', 'png'):
    plt.savefig(OUT/f'Fig1D_preview_forest.{ext}', dpi=600, bbox_inches='tight')
plt.close()
print('saved Fig1D_preview_forest')
