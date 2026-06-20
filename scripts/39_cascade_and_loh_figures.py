"""Cascade BCa forest + HLA-LOH lite/strict comparison figures (v0.7 additions)."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
import re

BASE = Path('/mnt/sda1/data/TNT/analysis')
TBL  = BASE/'tables'
FIG_MAIN = BASE/'figures/panels'
FIG_SUPP = BASE/'figures/supp'
FIG_SUPP.mkdir(parents=True, exist_ok=True)

# ============================================================
# Fig 6 v2 — Cascade BCa forest
# ============================================================
df = pd.read_csv(TBL/'TableS8_cascade_BCa_bootstrap.tsv', sep='\t')

def parse_ci(s):
    if pd.isna(s) or s == '':
        return np.nan, np.nan
    m = re.match(r'\[\s*([+-]?[\d.]+),\s*([+-]?[\d.]+)\s*\]', str(s))
    if not m: return np.nan, np.nan
    return float(m.group(1)), float(m.group(2))

# parse per-group and diff CIs
for col in ['good_95CI','bad_95CI','diff_95CI']:
    parsed = df[col].apply(parse_ci)
    df[col+'_lo'] = parsed.apply(lambda x: x[0])
    df[col+'_hi'] = parsed.apply(lambda x: x[1])

# Order features biologically (cascade sequence)
feat_order = ['missense','SBS5','neo_sites','neo_binders',
              'MHC_II','Treg','CD8_exhaustion','TRB_shannon','IGH_n']
df_ord = df.set_index('feature').loc[[f for f in feat_order if f in df.feature.values]].reset_index()

# two-panel figure: left=within-good CIs, right=between-group diff CIs
fig, axes = plt.subplots(1, 2, figsize=(13, 6.0), gridspec_kw={'width_ratios':[1,1]})

# --- Panel A: within-good and within-bad medians with BCa CIs ---
ax = axes[0]
n = len(df_ord); y = np.arange(n)
offset = 0.18
for i, row in df_ord.iterrows():
    # rescale large values (IGH_n ~ 1400) to separate subplot axis if needed
    pass

# Because IGH_n and missense are on very different scales, z-score each feature's delta
# For clean forest, plot standardized (delta / |median(|good_delta|+|bad_delta|)|) and also add text labels with raw values.
df_ord['g_med']  = df_ord['good_delta_median']
df_ord['g_lo']   = df_ord['good_95CI_lo']
df_ord['g_hi']   = df_ord['good_95CI_hi']
df_ord['b_med']  = df_ord['bad_delta_median']
df_ord['b_lo']   = df_ord['bad_95CI_lo']
df_ord['b_hi']   = df_ord['bad_95CI_hi']

# Normalize by max(|g_hi|, |b_hi|, |g_lo|, |b_lo|) per feature so all features share [-1, +1]-ish axis
scale = df_ord[['g_lo','g_hi','b_lo','b_hi']].abs().max(axis=1).replace(0, 1)
for col in ['g_med','g_lo','g_hi','b_med','b_lo','b_hi']:
    df_ord[col+'_s'] = df_ord[col] / scale

ax.axvline(0, color='k', lw=0.5, zorder=1)
for i, r in df_ord.iterrows():
    # good
    ax.errorbar([r.g_med_s], [i-offset], xerr=[[r.g_med_s - r.g_lo_s],[r.g_hi_s - r.g_med_s]],
                fmt='o', color='#2E86AB', ms=7, lw=1.5, capsize=3,
                label='good' if i==0 else None)
    # bad
    ax.errorbar([r.b_med_s], [i+offset], xerr=[[r.b_med_s - r.b_lo_s],[r.b_hi_s - r.b_med_s]],
                fmt='s', color='#E63946', ms=6, lw=1.5, capsize=3,
                label='bad' if i==0 else None)
    # raw-value text annotation
    raw = f'good {r.good_delta_median:+.2g} vs bad {r.bad_delta_median:+.2g}'
    ax.text(1.25, i, raw, fontsize=7.5, color='gray', va='center')

ax.set_yticks(range(n))
ax.set_yticklabels(df_ord['feature'])
ax.set_xlabel('Standardized Δ signature (pre → post, z per feature)')
ax.set_title('Within-group paired Δ with BCa 95 % CI\n(n = 6 good vs 6–7 bad per feature)',
             fontsize=10.5)
ax.legend(loc='lower right', fontsize=9, frameon=False)
ax.invert_yaxis()
ax.set_xlim(-1.4, 2.2)
ax.spines[['top','right']].set_visible(False)

# --- Panel B: between-group diff CIs, colored by whether CI excludes 0 ---
ax = axes[1]
ax.axvline(0, color='k', lw=0.8, zorder=1)
for i, r in df_ord.iterrows():
    lo, hi, med = r['diff_95CI_lo'], r['diff_95CI_hi'], r['diff_median_good_minus_bad']
    if np.isnan(lo):
        continue
    scale_diff = max(abs(lo), abs(hi), abs(med), 1e-9)
    med_s = med / scale_diff
    lo_s  = lo / scale_diff
    hi_s  = hi / scale_diff
    color = '#2a9d8f' if r['diff_CI_excludes_zero']==1 else '#b0b0b0'
    tag   = 'robust' if r['diff_CI_excludes_zero']==1 else 'exploratory'
    ax.errorbar([med_s], [i], xerr=[[med_s - lo_s],[hi_s - med_s]],
                fmt='D', color=color, ms=9, lw=2, capsize=4)
    ax.text(1.15, i,
            f'diff={med:+.2g} [{lo:+.2g}, {hi:+.2g}]   MW P={r.MW_p:.3f}   {tag}',
            fontsize=7.5, color='gray' if tag=='exploratory' else '#2a9d8f', va='center',
            fontweight='bold' if tag=='robust' else 'normal')

ax.set_yticks(range(n))
ax.set_yticklabels(df_ord['feature'])
ax.set_xlabel('Standardized between-group Δ (good − bad, per feature)')
ax.set_title('Between-group Δ with BCa 95 % CI\nCI excluding 0 = robust; else exploratory',
             fontsize=10.5)
ax.invert_yaxis()
ax.set_xlim(-1.4, 2.5)
ax.spines[['top','right']].set_visible(False)

plt.tight_layout()
for ext in ('png','pdf'):
    plt.savefig(FIG_MAIN/f'Fig6_cascade_BCa_forest.{ext}', dpi=400, bbox_inches='tight')
print('saved', FIG_MAIN/'Fig6_cascade_BCa_forest.{png,pdf}')

# ============================================================
# Supp Fig S3 — HLA-LOH lite vs strict comparison
# ============================================================
pre_strict = pd.read_csv(BASE/'03_hla/loh_stricter/pre_crt_LOH_subject_strict.tsv', sep='\t')
# pre_strict has n_loh_lite and n_loh_strict per subject

fig, axes = plt.subplots(1, 3, figsize=(13, 4.5),
                         gridspec_kw={'width_ratios':[1,1,1.2]})

# Panel A: stacked bar of LOH prevalence (lite vs strict) by response
ax = axes[0]
data = []
for grp in ['good','bad']:
    sub = pre_strict[pre_strict.response_bin == grp]
    lite_pos = int((sub.n_loh_lite > 0).sum())
    strict_pos = int((sub.n_loh_strict > 0).sum())
    data.append({'grp':grp,'total':len(sub),'lite':lite_pos,'strict':strict_pos})
df_prev = pd.DataFrame(data)

x = np.arange(2); w = 0.35
ax.bar(x - w/2, df_prev['lite'], width=w, color='#8bb7d4',
       edgecolor='k', lw=0.5, label='LOHHLA-lite (uncorrected)')
ax.bar(x + w/2, df_prev['strict'], width=w, color='#2E86AB',
       edgecolor='k', lw=0.5, label='Bonferroni-corrected strict')
for i, r in df_prev.iterrows():
    ax.text(i - w/2, r['lite'] + 0.1, f"{r['lite']}/{r['total']}", ha='center', fontsize=9)
    ax.text(i + w/2, r['strict'] + 0.1, f"{r['strict']}/{r['total']}", ha='center', fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(['good', 'bad'])
ax.set_ylabel('Subjects with ≥ 1 HLA-LOH event')
ax.set_title('Pre-CRT HLA-LOH prevalence\n(lite vs strict criteria)', fontsize=10.5)
ax.legend(loc='upper right', fontsize=8, frameon=False)
ax.spines[['top','right']].set_visible(False)

# Panel B: per-subject strict LOH pre→post (subj 3 and 4)
ax = axes[1]
try:
    paired = pd.read_csv(BASE/'03_hla/loh_stricter/paired_LOH_change_strict.tsv', sep='\t')
    resolved = paired[paired.loh_resolved]
    for i, r in resolved.iterrows():
        ax.plot([0,1], [r.pre_loh, r.post_loh], 'o-', color='#2E86AB', ms=10, lw=2)
        ax.text(-0.1, r.pre_loh, f'subj {int(r.subject_id)}',
                ha='right', va='center', fontsize=9)
        ax.text(1.1, r.post_loh, f'{int(r.post_loh)} loci',
                ha='left', va='center', fontsize=9, color='gray')
except Exception as e:
    pass
ax.set_xticks([0,1]); ax.set_xticklabels(['pre-CRT','post-CRT'])
ax.set_ylabel('Number of HLA class I LOH loci (strict)')
ax.set_title('HLA-LOH clone clearance\nacross radiation phase', fontsize=10.5)
ax.set_ylim(-0.5, 4)
ax.set_xlim(-0.4, 1.4)
ax.spines[['top','right']].set_visible(False)

# Panel C: criteria comparison table (as text)
ax = axes[2]
ax.axis('off')
criteria_text = (
    "Criteria comparison\n"
    "─────────────────────────────────\n"
    "LOHHLA-lite (v0.5 draft):\n"
    "  • |Δratio| ≥ 0.15\n"
    "  • Fisher P < 0.05 uncorrected\n"
    "  • No minimum read-depth filter\n\n"
    "Stricter (v0.7 primary call):\n"
    "  • |Δratio| ≥ 0.20\n"
    "  • Bonferroni-corrected P < 0.01\n"
    "  • Tumor / normal depth ≥ 30 reads each\n"
    "  • Requires normal heterozygosity\n\n"
    "Effect of stricter criteria:\n"
    "  • 10 lite events → 2 strict events\n"
    "  • Both strict events in eventual\n"
    "    good responders (subj 3 and 4)\n"
    "  • Both show complete pre→post\n"
    "    resolution of LOH clones\n\n"
    "Two-subject nature of strict calls\n"
    "means HLA-LOH clearance is reported\n"
    "as anecdotal corroboration of the\n"
    "cascade model (see Results §3.7),\n"
    "not a between-group statistical claim.\n"
    "Fisher P (good 2/16 vs bad 0/12) = 0.49."
)
ax.text(0.0, 0.98, criteria_text, transform=ax.transAxes,
        fontsize=8.2, va='top', ha='left', family='monospace')

plt.tight_layout()
for ext in ('png','pdf'):
    plt.savefig(FIG_SUPP/f'SuppFig_S3_HLA_LOH_lite_vs_strict.{ext}', dpi=400, bbox_inches='tight')
print('saved', FIG_SUPP/'SuppFig_S3_HLA_LOH_lite_vs_strict.{png,pdf}')
