"""
Shared style module for TNT publication figures.
Professional aesthetic inspired by Nature/Cell/NEJM.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats

# ---------- Colors ----------
GOOD = '#2a9d8f'      # teal (good responder)
BAD = '#e76f51'       # coral (bad responder)
PRE = '#457b9d'       # steel blue (pre-treatment)
POST = '#f4a261'      # orange (post-treatment)
NORMAL = '#8d99ae'    # grey-blue
BLACK = '#1d3557'     # navy
HIGHLIGHT = '#e9c46a' # gold

PAL_RESP = {'good': GOOD, 'bad': BAD}
PAL_TIME = {'pre': PRE, 'post': POST, 'normal': NORMAL}
PAL_STAGE = {'T2': '#a8dadc', 'T2/T3': '#457b9d', 'T3': '#1d3557', 'T4': '#e63946'}
PAL_MUT = {
    'missense_variant': '#2a9d8f',
    'frameshift_variant': '#1d3557',
    'stop_gained': '#e63946',
    'inframe_insertion': '#8338ec',
    'inframe_deletion': '#b15dff',
    'splice_acceptor_variant': '#f77f00',
    'splice_donor_variant': '#f4a261',
    'start_lost': '#d62828',
    'stop_lost': '#dc2f02',
    'protein_altering_variant': '#6a0572',
    'multi': '#264653',
}

# Category-based palette for GSEA
GSEA_CAT = {
    'Proliferation': '#2a9d8f',
    'DNA repair': '#118ab2',
    'Immune': '#ef476f',
    'Stromal/EMT': '#e76f51',
    'Metabolism': '#f4a261',
    'Signaling': '#8338ec',
    'Apoptosis': '#6a4c93',
    'Other': '#8d99ae',
}

# Gradient colormaps
CMAP_DIVERGING = LinearSegmentedColormap.from_list('resp_div',
    ['#e76f51','#f4a261','#ffffff','#a8dadc','#2a9d8f'])
CMAP_HEAT = 'RdBu_r'

def setup_style():
    """Apply professional rcParams."""
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'axes.titleweight': 'bold',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.9,
        'axes.edgecolor': '#1d3557',
        'axes.labelcolor': '#1d3557',
        'axes.axisbelow': True,
        'xtick.color': '#1d3557',
        'ytick.color': '#1d3557',
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'legend.frameon': False,
        'legend.fontsize': 8.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.dpi': 110,
        'figure.facecolor': 'white',
        'grid.color': '#e0e5ec',
        'grid.linewidth': 0.5,
    })

def save_panel(fig, name, outdir):
    """Save both PDF and high-res PNG."""
    from pathlib import Path
    p = Path(outdir); p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(p/f'{name}.png', dpi=400, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  ✓ {name}')

def sig_symbol(p):
    if pd.isna(p): return ''
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    if p < 0.1: return '†'
    return 'ns'

def panel_letter(ax, letter, x=-0.15, y=1.05, size=14):
    """Place bold panel letter (A, B, ...)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight='bold', va='top', ha='left', color='#1d3557')

def raincloud(ax, data, x, y, order=None, palette=None, width=0.4, alpha=0.7):
    """Raincloud plot: half-violin + box + jitter."""
    if order is None: order = sorted(data[x].unique())
    if palette is None: palette = {k: BLACK for k in order}
    for i, g in enumerate(order):
        vals = pd.to_numeric(data[data[x]==g][y], errors='coerce').dropna().values
        if len(vals) < 2: continue
        color = palette.get(g, BLACK)
        # Half violin (left of center)
        parts = ax.violinplot(vals, positions=[i-0.15], widths=width, showmeans=False,
                              showextrema=False, showmedians=False)
        for pc in parts['bodies']:
            # Clip to left half
            pc.set_facecolor(color); pc.set_alpha(0.55); pc.set_edgecolor(color); pc.set_linewidth(0.8)
            m = np.mean(pc.get_paths()[0].vertices[:,0])
            pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], -np.inf, m)
        # Box in middle
        bp = ax.boxplot([vals], positions=[i+0.08], widths=0.15, patch_artist=True,
                        showfliers=False, medianprops=dict(color='white', linewidth=1.5),
                        boxprops=dict(facecolor=color, alpha=0.9, edgecolor=color, linewidth=0.8),
                        whiskerprops=dict(color=color, linewidth=0.8),
                        capprops=dict(color=color, linewidth=0.8))
        # Jitter points on right
        jitter_x = i + 0.28 + np.random.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(jitter_x, vals, s=14, color=color, alpha=alpha,
                   edgecolor='white', linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)

def stat_bracket(ax, x1, x2, y, p, h=None, lw=0.9, fs=8):
    """Draw a statistical bracket with p-value."""
    if h is None:
        y1, y2 = ax.get_ylim(); h = (y2-y1)*0.03
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], color='black', lw=lw, clip_on=False)
    label = f'p = {p:.3g}' if p >= 0.001 else f'p = {p:.1e}'
    if p < 0.05: label += ' ' + sig_symbol(p)
    ax.text((x1+x2)/2, y+h*1.1, label, ha='center', va='bottom', fontsize=fs)

def add_axis_spines(ax, sides=('left','bottom'), color='#1d3557', lw=0.9):
    for side in ['top','right','left','bottom']:
        ax.spines[side].set_visible(side in sides)
        if side in sides:
            ax.spines[side].set_color(color)
            ax.spines[side].set_linewidth(lw)

def annotate_count(ax, data, x, y, order=None, dy=0.08):
    """Add n= labels below boxplot."""
    if order is None: order = sorted(data[x].unique())
    ymin = ax.get_ylim()[0]
    for i, g in enumerate(order):
        n = data[data[x]==g][y].dropna().shape[0]
        ax.text(i, ymin - dy, f'n={n}', ha='center', va='top', fontsize=8, color='#6c757d')
