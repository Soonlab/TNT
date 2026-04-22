"""
Regenerate Fig 3D / Fig 4D (v0.7.5) — signature-response forest lollipop.
Adds 3 externally-validated immune signatures as new rows:
  - CD8-cytotoxic      (GZMA/B, PRF1, IFNG, CD8A/B, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES)
  - T-cell infiltration (pure CD3 panel)
  - B-cell infiltration (CD19/MS4A1 panel, CD20 missing)

Rationale: manuscript §3.4 / Table 2 cite these three sigs (discovery P = 0.843 /
0.957 / 0.505; externally rescued in §3.12). Original fig3_v32 forest only
includes the legacy 22-signature panel, so they were invisible on the figure.

Output (overwrites):
  figures/panels_v3/Fig3D_forest_lollipop.{pdf,png}
Then the submission composite Fig 4 (physical filename Fig3_RNA_signatures.*) is
regenerated downstream via 43_rebuild_all_main_figures.py (user action or follow-up).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, '/data/data/TNT/analysis/scripts')
from _fig_style import setup_style, save_panel, sig_symbol, add_axis_spines  # noqa: E402
setup_style()

GOOD_DEEP = '#0a7d6e'
BAD_DEEP  = '#c53e1f'
BLACK_DEEP = '#0e2a47'
NEW_ACCENT = '#b58900'   # gold — marks the externally-validated rows

ROOT = Path('/data/data/TNT/analysis')
OUT  = ROOT/'figures/panels_v3'

# ---- 22-signature per-sample z-scores (already z-scored across all 56 samples)
sigs22 = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t')
sig_stats = pd.read_csv(ROOT/'06_rna_immune/sig_response_stats.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')

# ---- 3 additional purified immune signatures (raw ssGSEA), z-score across 56 samples
snew = pd.read_csv(ROOT/'260418_add/ssgsea_immune_all_samples.tsv', sep='\t')
NEW_SIGS = ['CD8_cytotoxic', 'Tcell_infiltration', 'Bcell_infiltration']
for c in NEW_SIGS:
    snew[c] = (snew[c] - snew[c].mean()) / snew[c].std()

sigs_ext = sigs22.merge(snew[['sample_id'] + NEW_SIGS], on='sample_id', how='left')
sigs_m = sigs_ext.merge(
    rna_inv[['sample_id','subject_id','timepoint','response_bin']],
    on='sample_id',
)

# ---- per-sig stats helper (mean good − mean bad, MW two-sided, 1000× bootstrap CI)
def sig_stats_row(name, values_good, values_bad, timepoint='pre'):
    g = np.asarray(values_good); b = np.asarray(values_bad)
    g = g[~np.isnan(g)]; b = b[~np.isnan(b)]
    mean_g, mean_b = g.mean(), b.mean()
    delta = mean_g - mean_b
    U, p = stats.mannwhitneyu(g, b, alternative='two-sided')
    rng = np.random.default_rng(42)
    diffs = [rng.choice(g, len(g), replace=True).mean() - rng.choice(b, len(b), replace=True).mean()
             for _ in range(1000)]
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return dict(timepoint=timepoint, signature=name, n_good=len(g), n_bad=len(b),
                mean_good=mean_g, mean_bad=mean_b, delta_good_minus_bad=delta,
                U=U, pvalue=p, ci_low=lo, ci_high=hi)

# ---- pre rows for the existing 22 sigs (re-compute CI the same way)
pre = sigs_m[sigs_m.timepoint=='pre']
rows = []
for _, r in sig_stats[sig_stats.timepoint=='pre'].iterrows():
    sig_n = r.signature
    g = pre[pre.response_bin=='good'][sig_n].dropna().values
    b = pre[pre.response_bin=='bad'][sig_n].dropna().values
    rows.append(sig_stats_row(sig_n, g, b))

# ---- pre rows for the 3 new externally-validated sigs
for sig_n in NEW_SIGS:
    g = pre[pre.response_bin=='good'][sig_n].dropna().values
    b = pre[pre.response_bin=='bad'][sig_n].dropna().values
    rows.append(sig_stats_row(sig_n, g, b))

df = pd.DataFrame(rows)
df = df.sort_values('pvalue').reset_index(drop=True)
df_plot = df.iloc[::-1].reset_index(drop=True)   # lowest P at top after flip

# ---- plot
LABEL_FIX = {
    'CD8_cytotoxic':      'CD8 cytotoxic ★',
    'Tcell_infiltration': 'T-cell infiltration ★',
    'Bcell_infiltration': 'B-cell infiltration ★',
}

fig, ax = plt.subplots(figsize=(9.3, 9.0))
y = np.arange(len(df_plot))
xmax_data = max(df_plot.ci_high.max(), df_plot.delta_good_minus_bad.max())
xmin_data = min(df_plot.ci_low.min(),  df_plot.delta_good_minus_bad.min())
right_margin = (xmax_data - xmin_data) * 0.55
ax.set_xlim(xmin_data*1.1, xmax_data + right_margin)

for i, (_, r) in enumerate(df_plot.iterrows()):
    is_new = r.signature in NEW_SIGS
    color = GOOD_DEEP if r.delta_good_minus_bad > 0 else BAD_DEEP
    # CI line
    ax.plot([r.ci_low, r.ci_high], [i, i], color=color, lw=2.0, alpha=0.85,
            solid_capstyle='round')
    # point estimate
    ax.scatter(r.delta_good_minus_bad, i,
               s=180 if r.pvalue < 0.05 else 100,
               color=color, edgecolor=NEW_ACCENT if is_new else 'white',
               linewidth=1.8 if is_new else 1.4, zorder=3)
    # p-value label
    star = sig_symbol(r.pvalue)
    if star == 'ns': star = ''
    label = f'p = {r.pvalue:.3g} {star}'
    ax.text(xmax_data + right_margin*0.04, i, label,
            va='center', ha='left', fontsize=9, color=BLACK_DEEP)

ax.axvline(0, color=BLACK_DEEP, lw=1.0)
ax.axvspan(-0.05, 0.05, color='#dee2e6', alpha=0.4)
ax.set_yticks(y)

# ---- styled y-tick labels: bold+gold for the 3 new externally-validated sigs
ytick_labels = []
for s in df_plot.signature:
    disp = LABEL_FIX.get(s, s.replace('_', ' '))
    ytick_labels.append(disp)
ax.set_yticklabels(ytick_labels, fontsize=10, color=BLACK_DEEP)
for tick, sig_n in zip(ax.get_yticklabels(), df_plot.signature):
    if sig_n in NEW_SIGS:
        tick.set_color(NEW_ACCENT)
        tick.set_fontweight('bold')

ax.set_xlabel('Δ z-score (good − poor), 95 % bootstrap CI',
              fontsize=11, fontweight='bold', color=BLACK_DEEP)

# legend note for ★ rows
ax.text(0.985, -0.06,
        '★ purified immune signatures (externally validated, §3.12)',
        transform=ax.transAxes, ha='right', va='top',
        fontsize=8.5, color=NEW_ACCENT, fontstyle='italic')

add_axis_spines(ax)
save_panel(fig, 'Fig3D_forest_lollipop', OUT)

# ---- also dump the augmented stats table for manuscript cross-check
out_tsv = ROOT/'260418_add/fig4D_forest_stats_with_immune3.tsv'
df.to_csv(out_tsv, sep='\t', index=False, float_format='%.5f')
print(f'wrote {out_tsv}')
print(f'wrote panel to {OUT}/Fig3D_forest_lollipop.{{pdf,png}} ({len(df_plot)} rows)')
