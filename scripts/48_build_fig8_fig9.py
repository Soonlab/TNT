"""Build Fig 8 (HLA/neoantigen) and Fig 9 (clonal evolution) composite figures
from panels_v3 subpanels, to match the v0.7.1 submission layout.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PV3  = Path('/data/data/TNT/analysis/figures/panels_v3')
MAIN = Path('/data/data/TNT/analysis/genome_medicine_submission/main_figures')

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

# ---------- Fig 8 — HLA / neoantigen ----------
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.12, wspace=0.05)
panels = [
    ('A', PV3/'Fig8A_HLA_alleles.png', 'HLA class I allele frequency (35 subjects)'),
    ('B', PV3/'Fig8B_HLA_homozygosity.png', 'HLA class I homozygosity by response'),
    ('C', PV3/'Fig8C_HLA_LOH.png', 'HLA class I LOH prevalence by response (strict/lite)'),
    ('D', PV3/'Fig8D_neoantigen_pre.png', 'Pre-CRT neoantigen burden by response'),
    ('E', PV3/'Fig8E_neoantigen_paired.png', 'Paired Δ neoantigen binders across radiation phase'),
    ('F', PV3/'Fig8F_neoantigen_lollipop.png', 'Per-subject neoantigen lollipop by response'),
]
for i, (L, p, t) in enumerate(panels):
    ax = fig.add_subplot(gs[i//3, i%3])
    place(ax, p, letter=L, title=t)
fig.suptitle('Figure 8. HLA class I landscape and MHC-I neoantigen cascade across radiation phase.',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig8_HLA_neoantigen.{ext}', dpi=600, bbox_inches='tight')
plt.close(fig); print('saved Fig8')

# ---------- Fig 9 — Clonal evolution ----------
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.12, wspace=0.05)
panels = [
    ('A', PV3/'Fig9A_clone_trajectories.png', 'Per-subject clone trajectories (PyClone-VI)'),
    ('B', PV3/'Fig9B_CCF_pre_post.png', 'Cellular cancer fraction pre vs post-CRT'),
    ('C', PV3/'Fig9C_cluster_stacked.png', 'Cluster composition stacked by subject'),
    ('D', PV3/'Fig9D_dominant_shrink.png', 'Dominant-clone shrinkage Δ by response'),
    ('E', PV3/'Fig9E_shrink_expand_scatter.png', 'Shrink vs expand CCF scatter'),
    ('F', PV3/'Fig9F_fate_composition.png', 'Clone-fate composition by response'),
]
for i, (L, p, t) in enumerate(panels):
    ax = fig.add_subplot(gs[i//3, i%3])
    place(ax, p, letter=L, title=t)
fig.suptitle('Figure 9. Clonal evolution during the radiation phase of TNT (n = 12 paired, PyClone-VI).',
             fontsize=13, fontweight='bold', y=0.995)
for ext in ('png','pdf'):
    fig.savefig(MAIN/f'Fig9_clonal_evolution.{ext}', dpi=600, bbox_inches='tight')
plt.close(fig); print('saved Fig9')
