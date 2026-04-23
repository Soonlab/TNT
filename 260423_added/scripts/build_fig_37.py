"""Build Fig_37 — HLA class I typing + HLA-LOH clone clearance.
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
OUT = Path('/mnt/sda1/data/TNT/analysis/260423_added/Fig_37_HLA_LOH')
OUT.mkdir(parents=True, exist_ok=True)

hom    = pd.read_csv(SRC/'hla_homozygosity.tsv', sep='\t')
hom_t  = pd.read_csv(SRC/'hla_homozygosity_tests.tsv', sep='\t')
subj   = pd.read_csv(SRC/'hla_loh_subject_summary.tsv', sep='\t')
tests  = pd.read_csv(SRC/'hla_loh_tests.tsv', sep='\t')
locus  = pd.read_csv(SRC/'hla_loh_per_locus.tsv', sep='\t')
clear  = pd.read_csv(SRC/'hla_loh_clearance.tsv', sep='\t')
conc   = pd.read_csv(SRC/'hla_loh_concordance.tsv', sep='\t')

ORIG = (2.0, 1.1)
SIZE = (6.2, 4.9)


def fmt_p(p):
    if pd.isna(p): return 'P = NA'
    if p < 0.001: return 'P < 0.001'
    return f'P = {p:.3f}'


# -------------------- Panel A: Homozygosity --------------------
def build_panel_A(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'A')
    order = ['homozygous_A','homozygous_B','homozygous_C','any_homozygous']
    labels = ['HLA-A', 'HLA-B', 'HLA-C', 'Any locus']
    yticks = [(0.0,'0%'),(0.1,'10'),(0.2,'20'),(0.3,'30')]
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(-0.5, len(order)-0.5), ylim=(0, 0.34),
           xticks=list(zip(range(len(order)), labels)),
           yticks=yticks,
           xlabel='HLA locus', ylabel='Homozygosity frequency (%)',
           y_grid_values=[t[0] for t in yticks])
    bar_w = 0.32
    for i, m in enumerate(order):
        r = hom_t[hom_t.metric == m].iloc[0]
        g_freq = r['good_freq']; b_freq = r['bad_freq']
        p.rect_data(i - bar_w, 0, i, g_freq, fill=GOOD, line=GOOD, line_width_pt=0.4)
        p.text(ORIG[0] + (i - bar_w/2 + 0.5) / len(order) * SIZE[0] - 0.5,
               p.ax.y2in(g_freq) - 0.26,
               f'{g_freq*100:.1f}%', fontsize=7.5, color=GOOD,
               w_in=1.0, h_in=0.22, align='center')
        p.rect_data(i, 0, i + bar_w, b_freq, fill=BAD, line=BAD, line_width_pt=0.4)
        p.text(ORIG[0] + (i + bar_w/2 + 0.5) / len(order) * SIZE[0] - 0.5,
               p.ax.y2in(b_freq) - 0.26,
               f'{b_freq*100:.1f}%', fontsize=7.5, color=BAD,
               w_in=1.0, h_in=0.22, align='center')
        # P-value above the tallest bar
        y_top = max(g_freq, b_freq) + 0.025
        p.text(ORIG[0] + (i + 0.5) / len(order) * SIZE[0] - 0.75,
               p.ax.y2in(y_top) - 0.11,
               fmt_p(r['fisher_p']), fontsize=7.5, italic=True, color=NEUTRAL,
               w_in=1.5, h_in=0.22, align='center')
    # Legend — above plot area, centered
    p.rect_in(3.5, 0.42, 0.18, 0.12, fill=GOOD, line=None)
    p.text(3.74, 0.40, f'good (n={(hom.response_bin=="good").sum()})',
           fontsize=9, color=GOOD, w_in=1.5, h_in=0.22)
    p.rect_in(5.3, 0.42, 0.18, 0.12, fill=BAD, line=None)
    p.text(5.54, 0.40, f'bad (n={(hom.response_bin=="bad").sum()})',
           fontsize=9, color=BAD, w_in=1.5, h_in=0.22)
    footer_note(p.slide,
                'Homozygosity per-locus Fisher exact (two-sided). '
                'Homozygous HLA reduces neopeptide repertoire diversity.')
    return p


# -------------------- Panel B/C: 2×2 Fisher mosaic --------------------
def build_mosaic_panel(prs, letter, call_label, footer_text):
    p = NativePanel(prs)
    panel_letter(p.slide, letter)
    row = tests[tests.call == call_label].iloc[0]
    g_n = row.good_pos + row.good_neg
    b_n = row.bad_pos + row.bad_neg
    total = g_n + b_n
    gw = g_n / total
    bw = b_n / total
    # Plot area with left margin for y-label
    x0, y0, w, h = 2.4, 1.2, 5.6, 4.7
    # x-axis baseline
    p.rect_in(x0, y0 + h, w, 0.03, fill=BLACK, line=None)
    good_w = gw * w
    bad_w  = bw * w
    good_top_frac = row.good_pos / g_n if g_n > 0 else 0
    bad_top_frac  = row.bad_pos  / b_n if b_n > 0 else 0
    good_top_h = good_top_frac * h
    good_bot_h = h - good_top_h
    p.rect_in(x0, y0,                 good_w, good_top_h, fill=GOOD,      line=WHITE, line_width_pt=1.2)
    p.rect_in(x0, y0 + good_top_h,    good_w, good_bot_h, fill=GOOD_FILL, line=WHITE, line_width_pt=1.2)
    bad_top_h = bad_top_frac * h
    bad_bot_h = h - bad_top_h
    p.rect_in(x0 + good_w, y0,                bad_w, bad_top_h, fill=BAD,      line=WHITE, line_width_pt=1.2)
    p.rect_in(x0 + good_w, y0 + bad_top_h,    bad_w, bad_bot_h, fill=BAD_FILL, line=WHITE, line_width_pt=1.2)
    def add_count_label(cx_in, cy_in, n, frac, bg_is_dark):
        color = WHITE if bg_is_dark else BLACK
        p.text(cx_in - 0.5, cy_in - 0.16,
               f'{n}\n({frac*100:.1f}%)', fontsize=10, bold=True, color=color,
               w_in=1.0, h_in=0.4, align='center', anchor='middle')
    if row.good_pos > 0:
        add_count_label(x0 + good_w/2, y0 + good_top_h/2, row.good_pos, good_top_frac, True)
    add_count_label(x0 + good_w/2, y0 + good_top_h + good_bot_h/2, row.good_neg, 1-good_top_frac, False)
    if row.bad_pos > 0:
        add_count_label(x0 + good_w + bad_w/2, y0 + bad_top_h/2, row.bad_pos, bad_top_frac, True)
    add_count_label(x0 + good_w + bad_w/2, y0 + bad_top_h + bad_bot_h/2, row.bad_neg, 1-bad_top_frac, False)
    # Group axis labels
    p.text(x0 + good_w/2 - 1.0, y0 + h + 0.10,
           f'good (n={g_n})', fontsize=10.5, bold=True, color=GOOD,
           w_in=2.0, h_in=0.26, align='center')
    p.text(x0 + good_w + bad_w/2 - 1.0, y0 + h + 0.10,
           f'bad (n={b_n})', fontsize=10.5, bold=True, color=BAD,
           w_in=2.0, h_in=0.26, align='center')
    p.text(x0, y0 + h + 0.42, 'TNT response  (block width ∝ group size)',
           fontsize=8.5, italic=True, color=NEUTRAL, w_in=w, h_in=0.25, align='center')
    # Y-axis label — clearer, placed further left with no overlap
    p.text_rotated(x0 - 1.8, y0 + h/2 - 0.9,
                   'HLA-LOH status   (top band = LOH+,  bottom = no LOH)',
                   fontsize=9, bold=True, w_in=1.8, h_in=0.3)
    # P-value — above plot area top-right
    p.text(5.5, 0.40, fmt_p(row.fisher_p),
           fontsize=12, bold=True, color=NEUTRAL, w_in=2.7, h_in=0.3, align='right')
    or_str = 'OR = ∞' if np.isinf(row.OR) else f'OR = {row.OR:.2f}'
    p.text(5.5, 0.68, f'Fisher exact, two-sided  —  {or_str}',
           fontsize=8.5, italic=True, color=GREY, w_in=2.7, h_in=0.25, align='right')
    footer_note(p.slide, footer_text)
    return p


# -------------------- Panel D: per-locus LOH freq --------------------
def build_panel_D(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'D')
    loci = ['HLA-A','HLA-B','HLA-C']
    yticks = [(0.0,'0%'),(0.05,'5'),(0.10,'10'),(0.15,'15'),(0.20,'20')]
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(-0.5, len(loci)-0.5), ylim=(0, 0.24),
           xticks=list(zip(range(len(loci)), loci)),
           yticks=yticks,
           xlabel='HLA locus', ylabel='HLA-LOH frequency (%) — strict criteria',
           y_grid_values=[t[0] for t in yticks])
    bar_w = 0.32
    for i, L in enumerate(loci):
        r_strict = locus[(locus.locus == L) & (locus.call == 'strict')].iloc[0]
        r_lite   = locus[(locus.locus == L) & (locus.call == 'lite')  ].iloc[0]
        for grp, color, xoff in [('good', GOOD, -bar_w/2 - 0.01),
                                  ('bad',  BAD,  +bar_w/2 + 0.01)]:
            freq = r_strict[f'{grp}_freq']
            npos = r_strict[f'{grp}_pos']; ntot = r_strict[f'{grp}_n_het']
            p.rect_data(i + xoff - bar_w/2, 0, i + xoff + bar_w/2, freq,
                        fill=color, line=color, line_width_pt=0.4)
            label = f'{int(npos)}/{int(ntot)}'
            p.text(ORIG[0] + (i + xoff + 0.5) / len(loci) * SIZE[0] - 0.55,
                   p.ax.y2in(freq) - 0.26,
                   label, fontsize=7.5, color=color,
                   w_in=1.1, h_in=0.22, align='center')
            # lite outline (where lite > strict)
            freq_lite = r_lite[f'{grp}_freq']
            if freq_lite > freq:
                p.rect_data(i + xoff - bar_w/2, freq,
                            i + xoff + bar_w/2, freq_lite,
                            fill=None, line=color, line_width_pt=1.0)
    # Legend above plot, top-right
    p.rect_in(6.0, 0.42, 0.18, 0.12, fill=GOOD, line=None)
    p.text(6.24, 0.40, 'good — strict', fontsize=8.5, color=GOOD, w_in=1.4, h_in=0.22)
    p.rect_in(7.6, 0.42, 0.18, 0.12, fill=BAD, line=None)
    p.text(7.84, 0.40, 'bad — strict',  fontsize=8.5, color=BAD, w_in=1.4, h_in=0.22)
    p.rect_in(6.0, 0.66, 0.18, 0.12, fill=WHITE, line=GREY, line_width_pt=1.0)
    p.text(6.24, 0.64, 'outline = lite only',
           fontsize=8, italic=True, color=GREY, w_in=2.0, h_in=0.22)
    footer_note(p.slide,
                'Strict: |Δratio| ≥ 0.20, Bonferroni-corrected Fisher P < 0.01, depth ≥ 30. '
                'Counts annotated as LOH+/het-in-normal.')
    return p


# -------------------- Panel E: paired pre→post clearance trajectory --------------------
def build_panel_E(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'E')
    candidates = clear[clear.pre_strict | clear.pre_lite].copy()
    candidates = candidates.sort_values('pre_imbalance', ascending=False).reset_index(drop=True)
    n = len(candidates)
    if n == 0:
        p.text(3, 3, 'No paired subject with pre-CRT LOH', fontsize=12)
        return p
    xticks = [(i, f'S{int(r.subject_id)} / {r.locus.replace("HLA-","")}')
              for i, r in candidates.iterrows()]
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(-0.5, n-0.5), ylim=(0, 0.50),
           xticks=xticks,
           yticks=[(0.0,'0'),(0.1,'0.1'),(0.2,'0.2'),(0.3,'0.3'),(0.4,'0.4'),(0.5,'0.5')],
           xlabel='Subject × HLA locus (sorted by pre-CRT allelic imbalance)',
           ylabel='Allelic imbalance  |tumor ratio − 0.5|',
           y_grid_values=[0.1,0.2,0.3,0.4,0.5])
    # Threshold line
    p.hline(0.20, color=GREY, width_pt=0.6, dash=True)
    p.text(ORIG[0] + SIZE[0] + 0.05, p.ax.y2in(0.20) - 0.10,
           'strict threshold  |Δ| ≥ 0.20',
           fontsize=7.5, italic=True, color=GREY, w_in=1.6, h_in=0.22)
    # Dumbbells
    for i, r in candidates.iterrows():
        color = GOOD if r.response_bin == 'good' else BAD
        p.line(i, r.pre_imbalance, i, r.post_imbalance, color=color, width_pt=1.8)
        p.dot(i, r.pre_imbalance,  d_in=0.17, color=color, edge=WHITE)
        p.dot(i, r.post_imbalance, d_in=0.15, color=WHITE, edge=color)
    # Legend (top-right, above plot)
    p.ellipse_in(6.4, 0.52, 0.14, fill=BLACK, line=WHITE, line_width_pt=0.5)
    p.text(6.55, 0.45, 'pre-CRT (filled)', fontsize=9, color=BLACK, w_in=1.5, h_in=0.22)
    p.ellipse_in(7.8, 0.52, 0.14, fill=WHITE, line=BLACK, line_width_pt=1.0)
    p.text(7.95, 0.45, 'post-CRT (open)', fontsize=9, color=BLACK, w_in=1.5, h_in=0.22)
    # Strict-clearance callout arrow overhead — ONE simple annotation above the plot,
    # not individual boxes (avoids x-axis label collision)
    s3s4_idx = candidates.index[candidates.subject_id.isin([3, 4])].tolist()
    if s3s4_idx:
        imin, imax = min(s3s4_idx), max(s3s4_idx)
        x_lo = ORIG[0] + (imin - 0.5 + 0.5) / n * SIZE[0]
        x_hi = ORIG[0] + (imax + 0.5 + 0.5) / n * SIZE[0]
        bracket_y = ORIG[1] + 0.05
        # bracket above strict-LOH events
        p.line_in(x_lo, bracket_y + 0.12, x_lo, bracket_y, color=NEUTRAL, width_pt=0.9)
        p.line_in(x_hi, bracket_y + 0.12, x_hi, bracket_y, color=NEUTRAL, width_pt=0.9)
        p.line_in(x_lo, bracket_y, x_hi, bracket_y, color=NEUTRAL, width_pt=0.9)
        p.text((x_lo + x_hi)/2 - 1.5, bracket_y - 0.28,
               'S3 & S4: 2-allele strict-LOH → full post-CRT clearance',
               fontsize=8, italic=True, color=NEUTRAL, w_in=3.0, h_in=0.22, align='center')
    footer_note(p.slide,
                'Imbalance = |tumor allele ratio − 0.5|. Filled = pre-CRT, open = post-CRT. '
                f'N = {n} (subject × locus) pre-CRT LOH-positive events among 13 paired subjects.')
    return p


# -------------------- Panel F: strict vs lite concordance --------------------
def build_panel_F(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'F')
    d = conc.copy()
    d['neglog10_p'] = -np.log10(d.fisher_p.clip(lower=1e-80))
    xmax = 0.55
    ymax = float(min(d.neglog10_p.max() * 1.1, 80))
    p.axes(origin_in=ORIG, size_in=SIZE,
           xlim=(-0.02, xmax), ylim=(0, ymax),
           xticks=[(0,'0'),(0.1,'0.1'),(0.2,'0.2'),(0.3,'0.3'),(0.4,'0.4'),(0.5,'0.5')],
           yticks=[(0,'0'),(2,'2'),(4,'4'),(10,'10'),(20,'20'),(40,'40'),(ymax,f'{int(ymax)}')],
           xlabel='|Δ ratio|  (|normal_ratio − tumor_ratio|)',
           ylabel='−log10 Fisher exact P (allele count)',
           y_grid_values=[2, 4, 10, 20, 40])
    p.line(0.20, 0, 0.20, ymax, color=GREY, width_pt=0.5, dash=True)
    p.line(-0.02, 2, xmax, 2, color=GREY, width_pt=0.5, dash=True)
    p.text(ORIG[0] + 0.20/xmax * SIZE[0] + 0.03,
           ORIG[1] + SIZE[1] + 0.05,
           '|Δ|=0.20', fontsize=7, italic=True, color=GREY, w_in=0.9, h_in=0.18)
    p.text(ORIG[0] + SIZE[0] + 0.03, p.ax.y2in(2) - 0.10,
           'P=0.01', fontsize=7, italic=True, color=GREY, w_in=0.7, h_in=0.18, align='left')
    # Strict region frame
    p.rect_data(0.20, 2, xmax, ymax, fill=None, line=NEUTRAL, line_width_pt=1.2)
    p.text(ORIG[0] + (0.20 + xmax)/2/xmax * SIZE[0] - 0.6, ORIG[1] + 0.05,
           'strict region', fontsize=8, italic=True, color=NEUTRAL,
           w_in=1.2, h_in=0.22, align='center')
    # Points
    for _, r in d.iterrows():
        color = GOOD if r.response_bin == 'good' else BAD
        if r.concordance == 'both':
            fill = color; edge = BLACK; dsize = 0.14
        elif r.concordance == 'lite_only':
            fill = WHITE; edge = color; dsize = 0.12
        else:
            fill = GREY_LT; edge = GREY; dsize = 0.09
        x_plot = max(-0.015, min(xmax - 0.005, r.delta_ratio))
        y_plot = max(0.05, min(ymax - 0.5, r.neglog10_p))
        p.dot(x_plot, y_plot, d_in=dsize, color=fill, edge=edge)
    # Legend above plot (top of slide)
    legend_items = [
        (GOOD, BLACK, 'good — strict', False),
        (BAD,  BLACK, 'bad — strict',  False),
        (WHITE, GOOD, 'good — lite only', True),
        (WHITE, BAD,  'bad — lite only',  True),
        (GREY_LT, GREY, 'neither (het)', False),
    ]
    lx = 0.5
    for fill, edge, text, is_open in legend_items:
        p.ellipse_in(lx, 0.58, 0.14, fill=fill, line=edge, line_width_pt=1.0)
        text_color = edge if is_open else (BLACK if fill == GREY_LT else fill)
        p.text(lx + 0.15, 0.51, text, fontsize=8.5,
               color=text_color, w_in=1.6, h_in=0.22)
        lx += 1.75
    footer_note(p.slide,
                'Each point = one (subject, locus) heterozygous event. '
                'Strict region defined by both x-threshold and P-threshold.')
    return p


# -------------------- assemble --------------------
prs = make_prs_blank()
build_panel_A(prs)
build_mosaic_panel(prs, 'B', 'strict any-LOH',
                   f'Per-subject strict HLA-LOH prevalence (pre-CRT). '
                   f'Good: {tests.iloc[0].good_pos}/{tests.iloc[0].good_pos+tests.iloc[0].good_neg}, '
                   f'Bad: {tests.iloc[0].bad_pos}/{tests.iloc[0].bad_pos+tests.iloc[0].bad_neg}.')
build_mosaic_panel(prs, 'C', 'lite any-LOH',
                   f'Per-subject LOHHLA-lite HLA-LOH (|Δratio| ≥ 0.15 + raw Fisher < 0.05). '
                   f'Good: {tests.iloc[1].good_pos}/{tests.iloc[1].good_pos+tests.iloc[1].good_neg}, '
                   f'Bad: {tests.iloc[1].bad_pos}/{tests.iloc[1].bad_pos+tests.iloc[1].bad_neg}.')
build_panel_D(prs)
build_panel_E(prs)
build_panel_F(prs)

out_pptx = OUT / 'Fig_37_HLA_LOH_panels.pptx'
prs.save(str(out_pptx))
print('Saved', out_pptx)
print('  panels:', len(list(prs.slides)))
