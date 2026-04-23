"""Build Fig_36 — MHC-I neoantigen landscape + radiation-phase clearance.
One panel per slide. No titles. Arial. Native shapes/textboxes only.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd
from _native_pptx import (NativePanel, make_prs_blank, panel_letter, footer_note,
                          GOOD, BAD, GOOD_FILL, BAD_FILL, GREY, GREY_LT,
                          NEUTRAL, BLACK, WHITE)

SRC = Path('/mnt/sda1/data/TNT/analysis/260423_added/source_data')
OUT = Path('/mnt/sda1/data/TNT/analysis/260423_added/Fig_36_neoantigen')
OUT.mkdir(parents=True, exist_ok=True)

per_sample = pd.read_csv(SRC/'neo_v2_per_sample.tsv', sep='\t')
pre_summary = pd.read_csv(SRC/'neo_v2_preCRT_summary.tsv', sep='\t')
paired = pd.read_csv(SRC/'neo_v2_paired_delta.tsv', sep='\t')
paired_summary = pd.read_csv(SRC/'neo_v2_paired_delta_summary.tsv', sep='\t')
bca = pd.read_csv(SRC/'neo_v2_bca_ci.tsv', sep='\t')
waterfall = pd.read_csv(SRC/'neo_v2_waterfall.tsv', sep='\t')

# Standard layout: plot area (2.0, 1.1) size (6.2, 4.9) → plot right edge 8.2
ORIG = (2.0, 1.1)
SIZE = (6.2, 4.9)


def fmt_p(p):
    if pd.isna(p): return 'P = NA'
    if p < 0.001: return 'P < 0.001'
    return f'P = {p:.3f}'


def nice_yticks_int(ymin, ymax, n=5):
    raw = np.linspace(ymin, ymax, n)
    return [(round(v), f'{int(round(v))}') for v in raw]


def legend_top_right(panel, items, x_start_in=7.7):
    """Put legend in upper-right of the plot area footer (above panel but below panel letter)."""
    y = 0.55
    for (color, label, is_outline) in items:
        if is_outline:
            panel.rect_in(x_start_in, y, 0.16, 0.11, fill=WHITE, line=color, line_width_pt=1.0)
        else:
            panel.rect_in(x_start_in, y, 0.16, 0.11, fill=color, line=None)
        panel.text(x_start_in + 0.22, y - 0.03, label, fontsize=8, color=color,
                   w_in=1.8, h_in=0.2)
        y += 0.22


# -------------------- Panel A/B/C: pre-CRT box + jitter --------------------
def build_pre_box_panel(prs, letter, metric, label_y, footer_stat):
    p = NativePanel(prs)
    panel_letter(p.slide, letter)
    pre = per_sample[per_sample.timepoint == 'pre']
    g = pre[pre.response == 'good'][metric].dropna().values
    b = pre[pre.response == 'bad'][metric].dropna().values
    vals_all = np.concatenate([g, b])
    ymin = min(0, float(np.min(vals_all)) * 1.05)
    ymax = float(np.max(vals_all)) * 1.12
    ticks = nice_yticks_int(ymin, ymax, 5)
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(0.2, 1.8), ylim=(ymin, ymax),
           xticks=[(0.6, f'good\n(n={len(g)})'), (1.2, f'bad\n(n={len(b)})')],
           yticks=ticks,
           xlabel='TNT response', ylabel=label_y,
           y_grid_values=[t[0] for t in ticks])
    p.boxplot(0.6, g, width=0.3, fill=GOOD_FILL, edge=GOOD)
    p.boxplot(1.2, b, width=0.3, fill=BAD_FILL, edge=BAD)
    p.jitter_scatter(0.6, g, width=0.08, color=GOOD, d_in=0.10, seed=3)
    p.jitter_scatter(1.2, b, width=0.08, color=BAD,  d_in=0.10, seed=5)
    p.diamond(0.6, float(np.median(g)), d_in=0.14, color=BLACK, edge=WHITE)
    p.diamond(1.2, float(np.median(b)), d_in=0.14, color=BLACK, edge=WHITE)
    # P-value top-right (outside plot area)
    p_val = pre_summary.set_index('metric').loc[metric, 'MW_p_twosided']
    p.text(6.2, 0.40, fmt_p(p_val), fontsize=11, bold=True, color=NEUTRAL,
           w_in=3.4, h_in=0.3, align='right')
    p.text(6.2, 0.66, 'Mann–Whitney U, two-sided', fontsize=8, italic=True, color=GREY,
           w_in=3.4, h_in=0.25, align='right')
    footer_note(p.slide, footer_stat)
    return p


# -------------------- Panel D/E/F: paired slope + meta --------------------
def build_paired_slope_panel(prs, letter, metric, ylabel, footer_stat):
    p = NativePanel(prs)
    panel_letter(p.slide, letter)
    g_vals = paired[paired.response == 'good'][metric].dropna().values
    b_vals = paired[paired.response == 'bad'][metric].dropna().values
    vals_all = np.concatenate([g_vals, b_vals])
    pad = 0.15 * (np.max(vals_all) - np.min(vals_all) + 1)
    ymin = float(np.min(vals_all)) - pad
    ymax = max(float(np.max(vals_all)) + pad, pad)
    ticks = nice_yticks_int(ymin, ymax, 6)
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(0.2, 1.8), ylim=(ymin, ymax),
           xticks=[(0.6, f'good\n(n={len(g_vals)})'), (1.2, f'bad\n(n={len(b_vals)})')],
           yticks=ticks,
           xlabel='TNT response', ylabel=ylabel,
           y_grid_values=[t[0] for t in ticks])
    p.hline(0, color=GREY, width_pt=0.6, dash=True)
    for v in g_vals:
        p.line(0.5, 0, 0.7, v, color=GOOD, width_pt=1.0)
        p.dot(0.5, 0, d_in=0.075, color=GREY, edge=WHITE)
        p.dot(0.7, v, d_in=0.13, color=GOOD, edge=WHITE)
    for v in b_vals:
        p.line(1.1, 0, 1.3, v, color=BAD,  width_pt=1.0)
        p.dot(1.1, 0, d_in=0.075, color=GREY, edge=WHITE)
        p.dot(1.3, v, d_in=0.13, color=BAD,  edge=WHITE)
    p.diamond(0.6, float(np.median(g_vals)), d_in=0.17, color=BLACK, edge=WHITE)
    p.diamond(1.2, float(np.median(b_vals)), d_in=0.17, color=BLACK, edge=WHITE)
    # Stats annotations (top-right header, above plot)
    ps = paired_summary.set_index('metric').loc[metric]
    p.text(6.2, 0.32, fmt_p(ps['MW_p_twosided']),
           fontsize=11, bold=True, color=NEUTRAL, w_in=3.4, h_in=0.3, align='right')
    p.text(6.2, 0.58, 'good vs bad (Mann–Whitney, two-sided)',
           fontsize=7.5, italic=True, color=GREY, w_in=3.4, h_in=0.25, align='right')
    p.text(6.2, 0.80, f'within-good Wilcoxon: {fmt_p(ps["wilcoxon_good_p"])}',
           fontsize=7.5, color=GOOD, w_in=3.4, h_in=0.22, align='right')
    p.text(6.2, 1.00, f'within-bad Wilcoxon: {fmt_p(ps["wilcoxon_bad_p"])}',
           fontsize=7.5, color=BAD, w_in=3.4, h_in=0.22, align='right')
    # Inline legend markers (top-left, in-inches so they stay OUTSIDE plot)
    p.ellipse_in(0.90, 0.58, 0.075, fill=GREY,  line=WHITE, line_width_pt=0.5)
    p.text(1.02, 0.51, 'pre-CRT baseline (0)', fontsize=7.5, italic=True,
           color=GREY, w_in=2.0, h_in=0.2)
    p.ellipse_in(0.90, 0.82, 0.12, fill=BLACK, line=WHITE, line_width_pt=0.5)
    p.text(1.02, 0.75, 'post-CRT Δ', fontsize=7.5, italic=True,
           color=BLACK, w_in=2.0, h_in=0.2)
    footer_note(p.slide, footer_stat)
    return p


# -------------------- Panel G: BCa CI forest --------------------
def build_forest_panel(prs, letter='G'):
    p = NativePanel(prs)
    panel_letter(p.slide, letter)
    metrics = [('delta_binders', 'Δ total binders'),
               ('delta_sites',   'Δ binder sites'),
               ('delta_strong',  'Δ strong binders'),
               ('delta_PCN',     'Δ PCN')]
    # Layout: plot area (2.8, 1.2) size (4.4, 4.8)
    p.axes(origin_in=(2.8, 1.2), size_in=(4.4, 4.8),
           xlim=(-1.1, 1.1), ylim=(0.0, 8.0),
           xticks=[(-1.0, '−max'), (0, '0'), (1.0, '+max')],
           yticks=None,
           xlabel='Δ normalized to per-metric |CI max|',
           ylabel='', y_grid_values=[0])
    p.vline(0, color=GREY, width_pt=0.6, dash=True)
    for yv in [2, 4, 6]:
        p.line(-1.1, yv, 1.1, yv, color=GREY_LT, width_pt=0.4)
    block_y = [7, 5, 3, 1]   # centers for 4 blocks (top to bottom)
    for (m, lbl), yc in zip(metrics, block_y):
        # metric label to the left of axis
        y_in = p.ax.y2in(yc) - 0.14
        p.text(1.4, y_in, lbl, fontsize=9.5, bold=True, color=NEUTRAL,
               w_in=1.35, h_in=0.3, align='right', anchor='middle')
        # rows
        good_row = bca[(bca.label == lbl) & (bca.group == 'good')].iloc[0]
        bad_row  = bca[(bca.label == lbl) & (bca.group == 'bad' )].iloc[0]
        scale = max(abs(good_row['ci_lo']), abs(good_row['ci_hi']),
                    abs(bad_row['ci_lo']),  abs(bad_row['ci_hi']), 1e-9)
        for row, color, yoff in [(good_row, GOOD, +0.45), (bad_row, BAD, -0.45)]:
            med = row['median'] / scale
            lo  = row['ci_lo']  / scale
            hi  = row['ci_hi']  / scale
            yv  = yc + yoff
            p.line(lo, yv, hi, yv, color=color, width_pt=2.0)
            p.line(lo, yv - 0.12, lo, yv + 0.12, color=color, width_pt=1.2)
            p.line(hi, yv - 0.12, hi, yv + 0.12, color=color, width_pt=1.2)
            p.dot(med, yv, d_in=0.14, color=color, edge=WHITE)
            # raw values printed to the right of axis
            p.text(7.4, p.ax.y2in(yv) - 0.10,
                   f'{row["group"]}: {row["median"]:.0f}  [{row["ci_lo"]:.0f}, {row["ci_hi"]:.0f}]',
                   fontsize=7.5, color=color, w_in=2.5, h_in=0.22, align='left', anchor='middle')
    # Legend (right of panel letter)
    p.text(1.4, 0.42, '●  median Δ with 95% BCa CI (bootstrap n=2000)',
           fontsize=8.5, italic=True, color=GREY, w_in=6.0, h_in=0.25)
    footer_note(p.slide,
                'Zero-crossing CI indicates non-significance. '
                'Normalization: each metric scaled to its own max |CI bound|; raw values shown on right.')
    return p


# -------------------- Panel H: per-subject waterfall --------------------
def build_waterfall_panel(prs, letter='H'):
    p = NativePanel(prs)
    panel_letter(p.slide, letter)
    n = len(waterfall)
    vals = waterfall.delta_binders.values
    ymin = min(0, float(vals.min())) * 1.1
    ymax = max(0, float(vals.max())) * 1.8
    if ymax == 0: ymax = 50
    yticks = nice_yticks_int(ymin, ymax, 6)
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(-0.5, n - 0.5), ylim=(ymin, ymax),
           xticks=[(i, f'#{int(s)}') for i, s in enumerate(waterfall.subject_id)],
           yticks=yticks,
           xlabel='Subject (sorted by Δ binders)',
           ylabel='Δ total MHC-I binder peptides (post − pre)',
           y_grid_values=[t[0] for t in yticks])
    p.hline(0, color=BLACK, width_pt=0.8, dash=False)
    for i, row in waterfall.iterrows():
        color = GOOD if row.response == 'good' else BAD
        bar_top = max(0, row.delta_binders)
        bar_bot = min(0, row.delta_binders)
        p.rect_data(i - 0.35, bar_bot, i + 0.35, bar_top,
                    fill=color, line=color, line_width_pt=0.6)
        ylbl = row.delta_binders + (0.045 * (ymax - ymin) if row.delta_binders >= 0
                                    else -0.045 * (ymax - ymin))
        p.text(2.0 + (i + 0.5) / n * 6.2 - 0.4,
               p.ax.y2in(ylbl) - 0.10,
               f'{int(row.delta_binders):+d}', fontsize=7.5, color=color,
               w_in=0.8, h_in=0.22, align='center')
    # Legend — top-left outside plot
    p.rect_in(2.0, 0.45, 0.18, 0.12, fill=GOOD, line=None)
    p.text(2.24, 0.42, 'good (pCR / near-CR)', fontsize=9, color=GOOD, w_in=2.4, h_in=0.2)
    p.rect_in(4.2, 0.45, 0.18, 0.12, fill=BAD,  line=None)
    p.text(4.44, 0.42, 'bad (PR / poor)',       fontsize=9, color=BAD,  w_in=2.0, h_in=0.2)
    footer_note(p.slide,
                f'N = {n} paired subjects (subject 13 excluded per manuscript '
                f'tumor-only rule). Negative bars = post-CRT neoantigen clearance.')
    return p


# -------------------- assemble --------------------
prs = make_prs_blank()
build_pre_box_panel(prs, 'A', 'n_binder_sites',
                    'Pre-CRT MHC-I binder sites',
                    'Full v2 cohort (good n=18 vs bad n=17). Each dot = one pre-CRT tumor. '
                    'Box = Q1/median/Q3, whiskers = 1.5×IQR. Black diamond = median.')
build_pre_box_panel(prs, 'B', 'n_binders',
                    'Pre-CRT total binder peptides',
                    'IC50 < 500 nM across all HLA–peptide predictions (v2 pvacseq re-run, MHCflurry).')
build_pre_box_panel(prs, 'C', 'PCN',
                    'Pre-CRT peptide-copy number (PCN)',
                    'PCN = Σ (2·VAF × binder count) across mutation sites.')
build_paired_slope_panel(prs, 'D', 'delta_binders',
                         'Δ total binder peptides (post − pre)',
                         'N = 11 paired subjects. Lines go from pre-baseline (0) to post-Δ. '
                         'Black diamond = group median.')
build_paired_slope_panel(prs, 'E', 'delta_sites',
                         'Δ MHC-I binder sites (post − pre)',
                         'Mutation sites with ≥1 binder after CRT minus those before CRT.')
build_paired_slope_panel(prs, 'F', 'delta_PCN',
                         'Δ PCN (post − pre)',
                         'PCN integrates allele-frequency-weighted binder count across sites.')
build_forest_panel(prs, 'G')
build_waterfall_panel(prs, 'H')

out_pptx = OUT / 'Fig_36_neoantigen_panels.pptx'
prs.save(str(out_pptx))
print('Saved', out_pptx)
print('  panels:', len(list(prs.slides)))
