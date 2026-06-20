"""Fig 2E (alt) v2 — genome-wide CNV heatmap with explicit Poor/Good grouping.

Redesign vs. previous:
  - Keep poor (top) / good (bottom) ordering consistent with Fig 2B, 2D
  - Add a prominent horizontal black separator between Poor and Good groups
  - Add vertical colored sidebar (Poor red / Good teal) spanning the row range
    of each group, plus bold bracket text "Poor (TRG 2-3) n=17" and
    "Good (TRG 0-1) n=18" to the left of sample labels
  - Keep 5 Mb binning and divergent CN cmap (0 homo-del → 4+ amp)

Output: panels_v3/Fig2E_CNV_genome.{pdf,png}  (overwrites)
"""
import sys
sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)
OUT_SUB = ROOT/'genome_medicine_submission/main_figures'; OUT_SUB.mkdir(parents=True, exist_ok=True)

wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
pre_samples = wes_inv[wes_inv.timepoint=='pre'].sample_id.tolist()

GOOD='#2E86AB'; BAD='#E63946'; SEP='#1d3557'
PAL_RESP = {'good':GOOD, 'bad':BAD}

CHR_SIZE = {'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,
            '6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,
            '11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,
            '16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,
            '21':46709983,'22':50818468,'X':156040895}
chr_list = list(CHR_SIZE.keys())

BIN_SIZE = 5_000_000
chr_starts = {}
cum = 0
for c in chr_list:
    chr_starts[c] = cum
    cum += CHR_SIZE[c]
total_len = cum
n_bins = total_len // BIN_SIZE + 1

CNV_DIR = ROOT/'04_wes_cnv_clonal/cnvkit'
cn_mat = []
sample_used = []
for s in pre_samples:
    f = CNV_DIR/f'{s}_DNA.call.cns'
    if not f.exists(): continue
    df = pd.read_csv(f, sep='\t')
    df = df[df.chromosome.astype(str).isin(chr_list)].copy()
    bins = np.full(n_bins, 2.0)
    for _, r in df.iterrows():
        chrom = str(r.chromosome)
        gs = chr_starts[chrom] + int(r.start)
        ge = chr_starts[chrom] + int(r.end)
        b1 = gs // BIN_SIZE; b2 = ge // BIN_SIZE + 1
        bins[b1:b2] = r.cn
    cn_mat.append(bins)
    sample_used.append(s)
cn_mat = np.array(cn_mat)

# ---- Order: bad (top) → good (bottom); within group, subject_id asc
meta = pd.DataFrame({
    'sample': sample_used,
    'response': [wes_inv[wes_inv.sample_id==s].response_bin.iloc[0] for s in sample_used],
    'subject': [int(wes_inv[wes_inv.sample_id==s].subject_id.iloc[0]) for s in sample_used],
})
meta['resp_rank'] = meta.response.map({'bad':0, 'good':1})
meta = meta.sort_values(['resp_rank','subject']).reset_index(drop=True)
order_idx = [sample_used.index(s) for s in meta['sample']]
cn_mat = cn_mat[order_idx]
sample_used = meta['sample'].tolist()
n_bad = int((meta.response=='bad').sum())
n_good = int((meta.response=='good').sum())
N = len(sample_used)

fig = plt.figure(figsize=(15.5, 8.2))
# 3 cols: group sidebar (0.25) | main (rest) | colorbar-ish (none; use inset)
# rows: chromosome ideogram (0.35) | main (main height)
gs = fig.add_gridspec(2, 2,
    height_ratios=[0.30, N*0.20],
    width_ratios=[0.9, n_bins*0.022],
    hspace=0.04, wspace=0.03)

# Top: chromosome ideogram
ax_chr = fig.add_subplot(gs[0,1])
chr_colors_alt = ['#cccccc','#888888']*15
for i, c in enumerate(chr_list):
    start_b = chr_starts[c]//BIN_SIZE
    width_b = CHR_SIZE[c]//BIN_SIZE
    ax_chr.add_patch(Rectangle((start_b, 0), width_b, 1,
                               facecolor=chr_colors_alt[i],
                               edgecolor='white', linewidth=0.4))
    ax_chr.text(start_b + width_b/2, 0.5, c, ha='center', va='center',
                fontsize=8.8, color='black', fontweight='bold')
ax_chr.set_xlim(0, n_bins); ax_chr.set_ylim(0, 1)
ax_chr.set_xticks([]); ax_chr.set_yticks([])
for s in ['top','right','left','bottom']: ax_chr.spines[s].set_visible(False)

# Left: group sidebar with Poor on top, Good on bottom
ax_side = fig.add_subplot(gs[1,0])
ax_side.set_xlim(0, 1); ax_side.set_ylim(0, N)
# Poor block (top)
ax_side.add_patch(Rectangle((0.55, N - n_bad), 0.25, n_bad,
                            facecolor=BAD, edgecolor='white', linewidth=0.6))
# Good block (bottom)
ax_side.add_patch(Rectangle((0.55, 0), 0.25, n_good,
                            facecolor=GOOD, edgecolor='white', linewidth=0.6))
# Labels rotated on sidebar
ax_side.text(0.40, N - n_bad/2, f'Poor\n(TRG 2–3)\nn = {n_bad}',
             ha='right', va='center', fontsize=11, fontweight='bold', color=BAD,
             rotation=0)
ax_side.text(0.40, n_good/2, f'Good\n(TRG 0–1)\nn = {n_good}',
             ha='right', va='center', fontsize=11, fontweight='bold', color=GOOD,
             rotation=0)
ax_side.set_xticks([]); ax_side.set_yticks([])
for s in ['top','right','left','bottom']: ax_side.spines[s].set_visible(False)

# Main CNV heatmap
ax_main = fig.add_subplot(gs[1,1])
# Diverging cmap
cnv_cmap = LinearSegmentedColormap.from_list('cnv',
    ['#08306b','#2171b5','#deebf7','white','#fee0d2','#fb6a4a','#67000d'])
im = ax_main.imshow(cn_mat, cmap=cnv_cmap, vmin=0, vmax=4,
                    aspect='auto', interpolation='nearest',
                    extent=[0, n_bins, 0, N], origin='upper')
# Sample labels on right side (so they don't collide with side bar)
ax_main.set_yticks([N - i - 0.5 for i in range(N)])
ax_main.set_yticklabels(sample_used, fontsize=7.2)
ax_main.yaxis.tick_right(); ax_main.yaxis.set_label_position('right')
ax_main.tick_params(length=0)
ax_main.set_xticks([]); ax_main.set_xlim(0, n_bins)

# Chromosome dividers
for i, c in enumerate(chr_list):
    if i == 0: continue
    ax_main.axvline(chr_starts[c]//BIN_SIZE, color='white', lw=0.6)

# Horizontal solid separator between Poor and Good
ax_main.axhline(n_good, color='black', lw=2.0, clip_on=False)

# Colorbar inset
cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.45])
cbar = fig.colorbar(im, cax=cbar_ax, ticks=[0,1,2,3,4])
cbar.set_label('Copy number', fontsize=9, color=SEP)
cbar.set_ticklabels(['0\n(homo del)','1\n(loss)','2\n(neutral)','3\n(gain)','4+\n(amp)'])
cbar.ax.tick_params(labelsize=7, color=SEP, labelcolor=SEP)

fig.suptitle('Genome-wide copy number landscape — Poor (top, n=17) / Good (bottom, n=18); 5 Mb bins, chr 1–22 + X',
             fontsize=12, fontweight='bold', y=0.965, color=SEP)

for d in (OUT, OUT_SUB):
    plt.savefig(d/'Fig2E_CNV_genome.pdf', bbox_inches='tight')
    plt.savefig(d/'Fig2E_CNV_genome.png', dpi=600, bbox_inches='tight')
print(f'Saved: {OUT}/Fig2E_CNV_genome.pdf/png  (+submission mirror)')
print(f'  Order: Poor top (n={n_bad}) / Good bottom (n={n_good}); total {N} pre-CRT samples')
plt.close()
