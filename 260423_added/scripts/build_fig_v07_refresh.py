"""Build Fig_8_v07_refresh — 6 panels matching Fig 8 A-F style of the v0.7
manuscript PPTX, rebuilt with v2 pvacseq cohort (46 tumors full) and all
new HLA/LOH stats. One panel per slide. No titles. Arial. Native shapes.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd
from _native_pptx import (NativePanel, make_prs_blank, panel_letter, footer_note,
                          half_violin, stacked_hbar,
                          GOOD, BAD, GOOD_FILL, BAD_FILL, GREY, GREY_LT,
                          NEUTRAL, BLACK, WHITE)
from pptx.dml.color import RGBColor

SRC = Path('/mnt/sda1/data/TNT/analysis/260423_added/source_data')
OUT = Path('/mnt/sda1/data/TNT/analysis/260423_added/Fig_8_v07_refresh')
OUT.mkdir(parents=True, exist_ok=True)

DARK_TEAL = RGBColor(0x1F, 0x4E, 0x5B)   # Fig 8A dark-teal bars
LIGHT_TEAL = RGBColor(0x3C, 0xA6, 0xB0)  # Fig 8A C-locus lighter
GOLD       = RGBColor(0xD6, 0xA4, 0x17)  # Fig 8A "from homozygous" overlay


def fmt_p(p):
    if pd.isna(p): return 'P = NA'
    if p < 0.001: return 'P < 0.001'
    return f'P = {p:.3f}'


# ========== Panel A: HLA allele top-10 (3 subplots, horizontal bar) ==========
def build_panel_A(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'A')
    allele = pd.read_csv(SRC/'hla_allele_topN.tsv', sep='\t')
    # three subplots side-by-side. Slide width 10", leave 0.5" margin both sides,
    # 3 plot columns of width 2.8", gap 0.3"
    colors = {'HLA-A': DARK_TEAL, 'HLA-B': DARK_TEAL, 'HLA-C': LIGHT_TEAL}
    for i, locus in enumerate(['HLA-A','HLA-B','HLA-C']):
        sub = allele[allele.locus == locus].sort_values('total', ascending=False).head(10)
        sub = sub.iloc[::-1].reset_index(drop=True)   # reverse so top-1 at top of chart
        ox = 0.9 + i * 3.05
        oy = 1.1
        pw = 2.6; ph = 4.8
        # Determine xmax per subplot
        xmax = float(sub.total.max()) * 1.2
        xticks_raw = np.linspace(0, xmax, 5)
        xticks = [(round(v), f'{int(round(v))}') for v in xticks_raw]
        yticks = list(zip(range(len(sub)),
                          [a.split(':')[0] + ':' + a.split(':')[1] if ':' in a else a
                           for a in sub.allele]))
        # Use panel.axes
        p.ax = None   # reset before next axes
        p.axes(origin_in=(ox, oy), size_in=(pw, ph),
               xlim=(0, xmax), ylim=(-0.5, len(sub)-0.5),
               xticks=xticks, yticks=yticks,
               xlabel='Allele count (2 × n_patients)', ylabel='',
               xlabel_fontsize=8, tick_fontsize=7, axis_lw_pt=0.6,
               y_grid_values=None)
        # y-axis allele labels are VERY long; redo them manually with narrow tickless left
        # Build bars
        for yi, r in sub.iterrows():
            stacked_hbar(p, 0, yi, 0,
                         [r.het_count, r.hom_count], [colors[locus], GOLD],
                         y_width=0.7)
            # total label
            p.text(ox + (r.total) / xmax * pw + 0.02,
                   p.ax.y2in(yi) - 0.10,
                   str(int(r.total)), fontsize=8, bold=True, color=NEUTRAL,
                   w_in=0.4, h_in=0.22, align='left')
        # locus header above subplot
        p.text(ox, 0.58, f'{locus}  (top 10)', fontsize=10.5, bold=True,
               color=colors[locus], w_in=pw, h_in=0.3, align='left')
    # Legend (slide top-right)
    p.rect_in(7.95, 0.40, 0.20, 0.12, fill=DARK_TEAL, line=None)
    p.text(8.20, 0.38, 'heterozygous allele count',
           fontsize=8, color=NEUTRAL, w_in=1.6, h_in=0.22)
    p.rect_in(7.95, 0.58, 0.20, 0.12, fill=GOLD, line=None)
    p.text(8.20, 0.56, 'from homozygous genotype',
           fontsize=8, color=NEUTRAL, w_in=1.6, h_in=0.22)
    footer_note(p.slide,
                'Allele counts per locus (2 alleles/patient); stacked yellow = contributed by '
                'homozygous genotype. n=35 patients (full v2 cohort).')
    return p


# ========== Panel B: Homozygosity strip-dot plot + MW P ==========
def build_panel_B(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'B')
    hom = pd.read_csv(SRC/'hla_homozygosity_per_subject.tsv', sep='\t')
    mw = pd.read_csv(SRC/'hla_homozygosity_mw.tsv', sep='\t').iloc[0]
    # x axis: n_homozygous_loci (0..3). Each dot = one patient; stack vertically within x bin
    # Separate columns for good (slight jitter left) and bad (jitter right)
    p.axes(origin_in=(2.0, 1.2), size_in=(6.2, 4.9),
           xlim=(-0.5, 3.5), ylim=(0, 14),
           xticks=[(0,'0'),(1,'1'),(2,'2'),(3,'3')],
           yticks=[(0,'0'),(2,'2'),(4,'4'),(6,'6'),(8,'8'),(10,'10'),(12,'12'),(14,'14')],
           xlabel='# homozygous HLA class-I loci',
           ylabel='Number of patients',
           y_grid_values=[2,4,6,8,10,12,14])
    # Stack dots vertically per (x, group)
    for grp, color, xoff in [('good', GOOD, -0.18), ('bad', BAD, +0.18)]:
        sub = hom[hom.response_bin == grp]
        counts_by_x = sub.groupby('n_homozygous_loci').size().to_dict()
        for x_bin, n in counts_by_x.items():
            for k in range(int(n)):
                yv = 0.5 + k
                p.dot(x_bin + xoff, yv, d_in=0.16, color=color, edge=WHITE)
            # count label above the column
            p.text(2.0 + (x_bin + xoff + 0.5) / 4 * 6.2 - 0.4,
                   p.ax.y2in(0.5 + n + 0.3) - 0.10,
                   str(int(n)), fontsize=9, bold=True, color=color,
                   w_in=0.8, h_in=0.22, align='center')
    # Group mean lines (vertical dashed)
    mean_good = float(mw['good_mean']); mean_bad = float(mw['bad_mean'])
    p.vline(mean_good, color=GOOD, width_pt=1.0, dash=True)
    p.vline(mean_bad,  color=BAD,  width_pt=1.0, dash=True)
    # Legend (left) + P-value (right), all outside plot at top
    p.rect_in(2.1, 0.40, 0.18, 0.12, fill=GOOD, line=None)
    p.text(2.32, 0.38, f'good (n={int(mw["n_good"])})',
           fontsize=9, color=GOOD, w_in=1.4, h_in=0.22)
    p.rect_in(3.6, 0.40, 0.18, 0.12, fill=BAD, line=None)
    p.text(3.82, 0.38, f'bad (n={int(mw["n_bad"])})',
           fontsize=9, color=BAD, w_in=1.4, h_in=0.22)
    p.text(6.2, 0.32, f'Mann-Whitney {fmt_p(mw["MW_p"])}',
           fontsize=11, bold=True, color=NEUTRAL, w_in=2.0, h_in=0.28, align='right')
    p.text(6.2, 0.58,
           f'good mean {mean_good:.2f} · bad mean {mean_bad:.2f}',
           fontsize=7.5, italic=True, color=GREY, w_in=2.0, h_in=0.22, align='right')
    footer_note(p.slide,
                f'Each dot = one patient. v2 full cohort n=35 (good {int(mw["n_good"])} / '
                f'bad {int(mw["n_bad"])}). Dashed lines = per-group mean.')
    return p


# ========== Panel C: LOH heatmap + side bar + per-allele Fisher ==========
def build_panel_C(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'C')
    heat = pd.read_csv(SRC/'hla_loh_heatmap_long.tsv', sep='\t')
    fisher = pd.read_csv(SRC/'hla_loh_per_allele_fisher.tsv', sep='\t')
    fisher_lite = fisher[fisher.call == 'LOHHLA-lite'].iloc[0]
    # Only include patients with ≥1 het locus (already in heatmap_long)
    patients = heat['subject_id'].unique().tolist()
    # Sort: bad first then good (like original Fig 8C), and within group by subject id asc
    patients_sorted = sorted(
        patients,
        key=lambda s: (0 if heat[heat.subject_id==s].response_bin.iloc[0]=='bad' else 1, s))
    # Plot area: heatmap on the left, bar on the right
    # Heatmap: 28 cols (patients) × 4 rows (response + A/B/C)
    # Use axes with no spines; draw cells as rect_in
    # heat area: (1.0, 1.5) to (7.3, 4.0), each cell ~0.22" × 0.35"
    hm_x0, hm_y0 = 1.0, 1.5
    n_pat = len(patients_sorted)
    cell_w = 6.4 / n_pat
    cell_h = 0.35
    row_y = {'response': hm_y0, 'HLA-A': hm_y0 + 0.55, 'HLA-B': hm_y0 + 0.92, 'HLA-C': hm_y0 + 1.29}
    # Response header row
    for i, sid in enumerate(patients_sorted):
        resp = heat[heat.subject_id==sid].response_bin.iloc[0]
        color = GOOD if resp == 'good' else BAD
        p.rect_in(hm_x0 + i*cell_w, row_y['response'], cell_w - 0.01, 0.28,
                  fill=color, line=WHITE, line_width_pt=0.6)
    # HLA-A/B/C rows
    for locus in ['HLA-A','HLA-B','HLA-C']:
        for i, sid in enumerate(patients_sorted):
            sub = heat[(heat.subject_id==sid) & (heat.locus==locus)]
            if len(sub) == 0:
                fill = GREY_LT   # no het data → grey
            else:
                is_loh = bool(sub.loh_lite.iloc[0])
                fill = BAD if is_loh else GREY_LT
            p.rect_in(hm_x0 + i*cell_w, row_y[locus], cell_w - 0.02, cell_h - 0.05,
                      fill=fill, line=WHITE, line_width_pt=0.4)
    # Row labels
    for label, yv in [('response', row_y['response']),
                       ('HLA-A',  row_y['HLA-A']),
                       ('HLA-B',  row_y['HLA-B']),
                       ('HLA-C',  row_y['HLA-C'])]:
        p.text(0.3, yv + 0.04, label, fontsize=9, bold=True, color=NEUTRAL,
               w_in=0.65, h_in=0.22, align='right')
    # Subject ID labels below (tight)
    for i, sid in enumerate(patients_sorted):
        p.text(hm_x0 + i*cell_w - 0.1, row_y['HLA-C'] + cell_h - 0.02,
               f'S{int(sid)}', fontsize=6.5, color=NEUTRAL,
               w_in=cell_w + 0.2, h_in=0.20, align='center')
    # Legend below heatmap
    p.rect_in(1.0, 3.45, 0.18, 0.11, fill=BAD, line=None)
    p.text(1.22, 3.43, 'LOH', fontsize=8, color=NEUTRAL, w_in=1.0, h_in=0.22)
    p.rect_in(1.9, 3.45, 0.18, 0.11, fill=GREY_LT, line=BLACK, line_width_pt=0.4)
    p.text(2.12, 3.43, 'retained / biallelic', fontsize=8, color=NEUTRAL, w_in=1.6, h_in=0.22)
    # Side bar: % LOH per group (horizontal bar)
    bar_x0 = 7.6; bar_w_total = 2.0
    p.text(bar_x0, 0.60, f'{fmt_p(fisher_lite.fisher_p)}',
           fontsize=11, bold=True, color=NEUTRAL, w_in=2.0, h_in=0.28)
    p.text(bar_x0, 0.87, 'Fisher exact, two-sided (per-allele)',
           fontsize=7, italic=True, color=GREY, w_in=2.0, h_in=0.22)
    # Bar axis: % LOH (0-30%)
    bar_y = {'bad': 1.60, 'good': 2.70}
    bar_max_pct = 30.0
    # axis line
    p.line_in(bar_x0, 3.10, bar_x0 + bar_w_total, 3.10, BLACK, width_pt=0.6)
    for pct in [0, 10, 20, 30]:
        xp = bar_x0 + pct / bar_max_pct * bar_w_total
        p.line_in(xp, 3.10, xp, 3.14, BLACK, width_pt=0.5)
        p.text(xp - 0.25, 3.17, f'{pct}', fontsize=7, color=NEUTRAL,
               w_in=0.5, h_in=0.2, align='center')
    p.text(bar_x0, 3.34, '% loci with LOH', fontsize=8, bold=True, color=NEUTRAL,
           w_in=bar_w_total, h_in=0.22, align='left')
    # bad bar
    bad_w_in = fisher_lite.bad_pct / bar_max_pct * bar_w_total
    p.rect_in(bar_x0, bar_y['bad'], bad_w_in, 0.40,
              fill=BAD, line=BAD, line_width_pt=0.4)
    p.text(bar_x0 - 0.55, bar_y['bad'] + 0.12, 'bad',
           fontsize=9, bold=True, color=BAD, w_in=0.5, h_in=0.22, align='right')
    p.text(bar_x0 + bad_w_in + 0.05, bar_y['bad'] + 0.12,
           f'{fisher_lite.bad_pos}/{fisher_lite.bad_total}  ({fisher_lite.bad_pct:.1f}%)',
           fontsize=8, bold=True, color=BAD, w_in=1.5, h_in=0.22)
    # good bar
    good_w_in = fisher_lite.good_pct / bar_max_pct * bar_w_total
    p.rect_in(bar_x0, bar_y['good'], good_w_in, 0.40,
              fill=GOOD, line=GOOD, line_width_pt=0.4)
    p.text(bar_x0 - 0.55, bar_y['good'] + 0.12, 'good',
           fontsize=9, bold=True, color=GOOD, w_in=0.5, h_in=0.22, align='right')
    p.text(bar_x0 + good_w_in + 0.05, bar_y['good'] + 0.12,
           f'{fisher_lite.good_pos}/{fisher_lite.good_total}  ({fisher_lite.good_pct:.1f}%)',
           fontsize=8, bold=True, color=GOOD, w_in=1.5, h_in=0.22)
    footer_note(p.slide,
                f'Per-allele Fisher (LOHHLA-lite): good {fisher_lite.good_pos}/{fisher_lite.good_total} '
                f'vs bad {fisher_lite.bad_pos}/{fisher_lite.bad_total}, P = {fisher_lite.fisher_p:.3f}. '
                f'Heatmap rows = HLA-A/B/C per patient (het normal only). Grey = retained or non-het.')
    return p


# ========== Panel D: Pre-CRT raincloud (3 sub-panels) ==========
def build_panel_D(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'D')
    per_sample = pd.read_csv(SRC/'neo_v2_per_sample.tsv', sep='\t')
    pre_summary = pd.read_csv(SRC/'neo_v2_preCRT_summary.tsv', sep='\t').set_index('metric')
    metrics = [('n_binder_sites',    'mutation sites with MHC-I binder (<500 nM)'),
               ('n_strong_binders',  'strong MHC-I binders (<50 nM)'),
               ('PCN',               'PCN score (HLA-LOH adjusted proxy)')]
    pre = per_sample[per_sample.timepoint == 'pre']
    for i, (metric, ylab) in enumerate(metrics):
        ox = 1.4 + i * 2.85
        oy = 1.35
        pw = 2.2; ph = 4.6
        g = pre[pre.response=='good'][metric].dropna().values
        b = pre[pre.response=='bad' ][metric].dropna().values
        vmin = min(float(min(g.min(), b.min())), 0)
        vmax = float(max(g.max(), b.max())) * 1.15
        ticks_raw = np.linspace(vmin, vmax, 5)
        ticks = [(round(t,1), f'{int(round(t))}' if abs(t)>=10 else f'{t:.1f}') for t in ticks_raw]
        p.ax = None
        p.axes(origin_in=(ox, oy), size_in=(pw, ph),
               xlim=(-0.5, 2.0), ylim=(vmin, vmax),
               xticks=[(0.3, 'good'), (1.3, 'bad')],
               yticks=ticks,
               xlabel='', ylabel=ylab,
               xlabel_fontsize=8, ylabel_fontsize=8, tick_fontsize=7,
               y_grid_values=[t[0] for t in ticks])
        # Raincloud: half-violin on LEFT, open box+scatter on RIGHT
        half_violin(p, 0.3, g, side='left', max_width=0.30,
                    color=GOOD, fill=GOOD_FILL, edge_width_pt=0.7)
        half_violin(p, 1.3, b, side='left', max_width=0.30,
                    color=BAD,  fill=BAD_FILL,  edge_width_pt=0.7)
        # Open box (no fill, thin)
        for xc, vals, color in [(0.3, g, GOOD), (1.3, b, BAD)]:
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            p.rect_data(xc + 0.18, q1, xc + 0.40, q3,
                        fill=None, line=color, line_width_pt=1.0)
            p.line(xc + 0.18, med, xc + 0.40, med, color=color, width_pt=1.2)
        # Jittered dots on far right
        p.jitter_scatter(0.68, g, width=0.08, color=GOOD, d_in=0.09, seed=3)
        p.jitter_scatter(1.68, b, width=0.08, color=BAD,  d_in=0.09, seed=7)
        # P-value above plot
        p_val = pre_summary.loc[metric, 'MW_p_twosided']
        p.text(ox, 0.78, fmt_p(p_val), fontsize=10, bold=True,
               color=NEUTRAL, w_in=pw, h_in=0.26, align='center')
        p.text(ox, 1.02, f'good n={len(g)} / bad n={len(b)}',
               fontsize=7, italic=True, color=GREY, w_in=pw, h_in=0.2, align='center')
    footer_note(p.slide,
                'Half-violin = kernel density, open box = IQR with median line, '
                'jittered points = individual pre-CRT tumors. v2 full cohort (n=35).')
    return p


# ========== Panel E: Paired pre→post slopes (2 metric) ==========
def build_panel_E(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'E')
    paired = pd.read_csv(SRC/'neo_v2_paired_per_subject.tsv', sep='\t')
    paired_summary = pd.read_csv(SRC/'neo_v2_paired_delta_summary.tsv', sep='\t').set_index('metric')
    metrics = [('sites',  'n_binder_sites',    'delta_sites',
                '# sites with MHC-I binder'),
               ('strong', 'n_strong_binders',  'delta_strong',
                '# strong binders (<50 nM)')]
    for i, (metric_key, _, delta_metric, ylab) in enumerate(metrics):
        ox = 1.0 + i * 4.3
        oy = 1.35
        pw = 3.8; ph = 4.6
        pre_col = f'pre_{metric_key}'; post_col = f'post_{metric_key}'
        all_vals = pd.concat([paired[pre_col], paired[post_col]]).dropna().values
        vmin = 0
        vmax = float(all_vals.max()) * 1.12
        ticks_raw = np.linspace(0, vmax, 6)
        ticks = [(round(t), f'{int(round(t))}') for t in ticks_raw]
        p.ax = None
        p.axes(origin_in=(ox, oy), size_in=(pw, ph),
               xlim=(-0.2, 1.2), ylim=(vmin, vmax),
               xticks=[(0,'pre'), (1,'post')],
               yticks=ticks,
               xlabel='', ylabel=ylab,
               xlabel_fontsize=9, ylabel_fontsize=9, tick_fontsize=8,
               y_grid_values=[t[0] for t in ticks])
        # Per-subject lines
        for _, r in paired.iterrows():
            color = GOOD if r.response == 'good' else BAD
            p.line(0, r[pre_col], 1, r[post_col], color=color, width_pt=0.8)
            p.dot(0, r[pre_col], d_in=0.09, color=color, edge=WHITE)
            p.dot(1, r[post_col], d_in=0.09, color=color, edge=WHITE)
            # subject label on left of pre
            p.text(ox - 0.18, p.ax.y2in(r[pre_col]) - 0.09,
                   f'S{int(r.subject_id)}', fontsize=5.5, color=color,
                   w_in=0.25, h_in=0.18, align='right')
        # Group mean thick line + big dots
        for grp, color in [('good', GOOD), ('bad', BAD)]:
            sub = paired[paired.response == grp]
            pre_mean = float(sub[pre_col].mean())
            post_mean = float(sub[post_col].mean())
            p.line(0, pre_mean, 1, post_mean, color=color, width_pt=2.8)
            p.dot(0, pre_mean, d_in=0.18, color=color, edge=WHITE)
            p.dot(1, post_mean, d_in=0.18, color=color, edge=WHITE)
            # Δ label near post dot
            delta_label = float(sub[post_col].mean() - sub[pre_col].mean())
            p.text(ox + (1.05 + 0.2) / 1.4 * pw - 0.5,
                   p.ax.y2in(post_mean) - 0.25,
                   f'{grp}  Δ = {delta_label:+.1f}',
                   fontsize=8, bold=True, color=color,
                   w_in=1.6, h_in=0.22)
        # Stats box top-left of subplot
        ps = paired_summary.loc[delta_metric]
        p.text(ox, 0.78, f'good (n=5): {fmt_p(ps["wilcoxon_good_p"])}',
               fontsize=8, color=GOOD, w_in=pw, h_in=0.2, align='left')
        p.text(ox, 1.00, f'bad  (n=6): {fmt_p(ps["wilcoxon_bad_p"])}',
               fontsize=8, color=BAD, w_in=pw, h_in=0.2, align='left')
        p.text(ox + pw - 1.8, 0.78,
               f'MW {fmt_p(ps["MW_p_twosided"])}',
               fontsize=8, italic=True, color=NEUTRAL, w_in=1.8, h_in=0.2, align='right')
    footer_note(p.slide,
                'Per-subject thin lines, group-mean thick line. '
                'Wilcoxon signed-rank within each group; Mann-Whitney U between groups. '
                'N=11 paired subjects (good 5, bad 6; subj 13 excluded).')
    return p


# ========== Panel F: Per-subject lollipop + HLA annotation columns ==========
def build_panel_F(prs):
    p = NativePanel(prs)
    panel_letter(p.slide, 'F')
    df = pd.read_csv(SRC/'neo_v2_lollipop.tsv', sep='\t')
    # Match original Fig 8F convention: bad at TOP, good at BOTTOM, each sorted
    # by pre binder sites ASCENDING (smallest at group top, largest at group bottom).
    # Use explicit response_order so pandas groups correctly regardless of dtype.
    df['_order'] = df.response.map({'bad': 0, 'good': 1})
    df = df.sort_values(['_order','n_binder_sites'], ascending=[True, True]).reset_index(drop=True)
    df = df.drop(columns='_order')
    n = len(df)
    xmax = float(df.n_binder_sites.max()) * 1.05
    # Layout: 4 annotation cols on left (response, HLA-A zygosity, B, C), lollipop axis on right
    # slide width 10" — leave 0.4 left margin.
    anno_w = 0.35
    x_resp   = 0.5
    x_A      = x_resp + anno_w + 0.05
    x_B      = x_A + anno_w + 0.05
    x_C      = x_B + anno_w + 0.05
    hm_x0    = x_C + anno_w + 0.25
    hm_w     = 10.0 - hm_x0 - 0.5
    row_h    = min(0.18, 5.3 / n)
    y_top    = 0.95
    # Annotation column headers
    for hx, hlabel in [(x_resp, 'resp'),
                        (x_A,    'HLA-A'),
                        (x_B,    'HLA-B'),
                        (x_C,    'HLA-C')]:
        p.text(hx - 0.02, y_top - 0.28, hlabel, fontsize=7, italic=True, bold=True,
               color=NEUTRAL, w_in=anno_w, h_in=0.22, align='center')
    # Lollipop axis region — create axes
    p.axes(origin_in=(hm_x0, y_top), size_in=(hm_w, row_h * n + 0.1),
           xlim=(0, xmax), ylim=(-0.5, n - 0.5),
           xticks=[(v, f'{int(v)}') for v in np.linspace(0, xmax, 6).round().astype(int)],
           yticks=None,
           xlabel='# mutation sites with MHC-I binder (pre)',
           ylabel='',
           xlabel_fontsize=9, tick_fontsize=7,
           y_grid_values=None)
    # Annotation rectangles + labels + lollipops: i=0 at TOP of slide.
    # df is sorted with bad first (i=0..16), then good (i=17..34).
    # Within each group, ascending by binder sites.
    # To render lollipop with consistent y-axis, map row i → data y = (n-1-i)
    # so annotation row i=0 (top) aligns with lollipop data y = n-1 (also axis-top).
    for i, r in df.iterrows():
        y_data = (len(df) - 1 - i)   # data-space y (axis-top when i=0)
        y_in_center = y_top + i * row_h + row_h / 2
        # response column
        resp_color = GOOD if r.response == 'good' else BAD
        p.rect_in(x_resp, y_in_center - row_h*0.35, anno_w, row_h*0.7,
                  fill=resp_color, line=None)
        # zygosity columns
        for hx, is_hom in [(x_A, r.homozygous_A),
                           (x_B, r.homozygous_B),
                           (x_C, r.homozygous_C)]:
            if bool(is_hom):
                p.rect_in(hx, y_in_center - row_h*0.35, anno_w, row_h*0.7,
                          fill=GOLD, line=None)
                p.text(hx, y_in_center - 0.06, 'hom', fontsize=5.5, bold=True,
                       color=WHITE, w_in=anno_w, h_in=0.14, align='center')
            else:
                p.rect_in(hx, y_in_center - row_h*0.35, anno_w, row_h*0.7,
                          fill=GREY_LT, line=None)
                p.text(hx, y_in_center - 0.06, 'het', fontsize=5.5,
                       color=NEUTRAL, w_in=anno_w, h_in=0.14, align='center')
        # Subject label just left of lollipop axis
        p.text(hm_x0 - 0.45, y_in_center - 0.06,
               f'S{int(r.subject_id)}', fontsize=6.5, color=NEUTRAL,
               w_in=0.40, h_in=0.16, align='right')
        # Lollipop stem (use y_data = (n-1-i) so i=0 lands at axis-top)
        p.line(0, y_data, r.n_binder_sites, y_data, color=resp_color, width_pt=1.0)
        # Strong binders small marker (inner)
        p.dot(r.n_strong_binders, y_data, d_in=0.08, color=GOLD, edge=NEUTRAL)
        # Main marker (outer) — X if any LOH, circle if retained
        marker_size = 0.16
        if bool(r.any_loh_strict):
            x_em = p.ax.x2in(r.n_binder_sites)
            y_em = p.ax.y2in(y_data)
            s = 0.11
            p.line_in(x_em - s, y_em - s, x_em + s, y_em + s, color=resp_color, width_pt=2.0)
            p.line_in(x_em - s, y_em + s, x_em + s, y_em - s, color=resp_color, width_pt=2.0)
        else:
            p.dot(r.n_binder_sites, y_data, d_in=marker_size, color=resp_color, edge=NEUTRAL)
        # Value label right of axis
        p.text(10.0 - 0.6, p.ax.y2in(y_data) - 0.08,
               f'{int(r.n_binder_sites)}', fontsize=6.5, bold=True, color=resp_color,
               w_in=0.5, h_in=0.16, align='right')
    # Legend bottom (well below plot area)
    leg_x = 0.6; leg_y = 7.30
    p.dot_radius = True
    p.ellipse_in(leg_x, leg_y, 0.15, fill=GOOD, line=NEUTRAL, line_width_pt=0.5)
    p.text(leg_x + 0.20, leg_y - 0.05, 'HLA retained (good)', fontsize=8,
           color=NEUTRAL, w_in=1.8, h_in=0.2)
    p.ellipse_in(leg_x + 2.2, leg_y, 0.15, fill=BAD, line=NEUTRAL, line_width_pt=0.5)
    p.text(leg_x + 2.4, leg_y - 0.05, 'HLA retained (bad)', fontsize=8,
           color=NEUTRAL, w_in=1.8, h_in=0.2)
    # X icon
    p.line_in(leg_x + 4.3 - 0.09, leg_y - 0.09, leg_x + 4.3 + 0.09, leg_y + 0.09,
              color=NEUTRAL, width_pt=1.8)
    p.line_in(leg_x + 4.3 - 0.09, leg_y + 0.09, leg_x + 4.3 + 0.09, leg_y - 0.09,
              color=NEUTRAL, width_pt=1.8)
    p.text(leg_x + 4.5, leg_y - 0.05, 'HLA strict LOH', fontsize=8,
           color=NEUTRAL, w_in=1.6, h_in=0.2)
    p.ellipse_in(leg_x + 6.2, leg_y, 0.10, fill=GOLD, line=NEUTRAL, line_width_pt=0.5)
    p.text(leg_x + 6.37, leg_y - 0.05, 'strong binders (<50 nM)',
           fontsize=8, color=NEUTRAL, w_in=2.0, h_in=0.2)
    return p


# ========== assemble ==========
prs = make_prs_blank()
build_panel_A(prs)
build_panel_B(prs)
build_panel_C(prs)
build_panel_D(prs)
build_panel_E(prs)
build_panel_F(prs)

out_pptx = OUT / 'Fig_8_v07_refresh_panels.pptx'
prs.save(str(out_pptx))
print('Saved', out_pptx)
print('  panels:', len(list(prs.slides)))
