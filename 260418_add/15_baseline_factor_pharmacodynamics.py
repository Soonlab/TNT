#!/usr/bin/env python3
"""
15_baseline_factor_pharmacodynamics.py

Question
--------
The four baseline Thread-1 factors (DSB_HDR_repair, Tumor_cellcycle,
E2F_MYC_cellcycle, EMT) were externally validated as pre-treatment predictors
of final TNT response. If these factors genuinely index the biology that
radiation targets, then the paired pre/post biopsies of our 12 paired
subjects (RT-alone window) should show the factors themselves moving in
biologically predictable directions:

    DSB_HDR_repair        : DOWN post-RT (repair-proficient clones killed)
    Tumor_cellcycle       : DOWN post-RT (cycling clones killed)
    E2F_MYC_cellcycle     : DOWN post-RT (proliferation axis killed)
    EMT                   : UP   post-RT (surviving tissue enriches for mesenchymal)

We report two independent axes of evidence for each factor:

    1. Magnitude  (signed Wilcoxon, MW between groups) --- how big is the
       change, and is it different between good and bad responders?
    2. Direction  (sign / binomial, Fisher between groups) --- how
       consistently do subjects move in the predicted direction, regardless
       of the per-subject magnitude?

The second axis is necessary because small-N paired designs often show
reproducible direction but the magnitude varies subject-to-subject (so
Wilcoxon/MW lose power), and because a "good = coherent response, bad =
stochastic response" pattern is invisible to magnitude tests.

Inputs
------
    08_rna_pathway/ssgsea_scores.tsv          sample x pathway (95 cols)
    00_cohort/rna_inventory.tsv               sample->subject/timepoint/response

Outputs (to 260418_add/)
------------------------
    baseline_factor_per_subject_delta.tsv     per-subject Δ for each factor/composite
    baseline_factor_pharmacodynamics_stats.tsv  summary per factor+metric
    baseline_factor_sign_table.tsv            sign counts per group per factor
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/data/data/TNT/analysis"
OUT = f"{ROOT}/260418_add"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# 1) Load ssGSEA + inventory, restrict to paired subjects
# ---------------------------------------------------------------------------
ssgsea = pd.read_csv(f"{ROOT}/08_rna_pathway/ssgsea_scores.tsv", sep="\t")
inv = pd.read_csv(f"{ROOT}/00_cohort/rna_inventory.tsv", sep="\t")

# paired = subjects with both pre and post RNA
pairs = (inv[inv.timepoint.isin(["pre", "post"])]
         .groupby("subject_id")["timepoint"]
         .nunique())
paired_subjects = sorted(pairs[pairs == 2].index.tolist())
assert len(paired_subjects) == 12, paired_subjects

inv_p = inv[(inv.subject_id.isin(paired_subjects)) &
            (inv.timepoint.isin(["pre", "post"]))].copy()

# ---------------------------------------------------------------------------
# 2) Define 4 baseline factors as composites of representative pathways
#    (and keep individual members for transparency)
# ---------------------------------------------------------------------------
FACTOR_DEF = {
    "DSB_HDR_repair": {
        "members": [
            "DNA Double-Strand Break Repair R-HSA-5693532",
            "HDR Thru Homologous Recombination (HRR) R-HSA-5685942",
            "Homology Directed Repair R-HSA-5693538",
            "Processing Of DNA Double-Strand Break Ends R-HSA-5693607",
        ],
        "predicted_direction": "down",  # in both groups post-RT
    },
    "Tumor_cellcycle": {
        "members": [
            "G2-M Checkpoint",
            "Cell Cycle Checkpoints R-HSA-69620",
            "M Phase R-HSA-68886",
            "Mitotic G2-G2/M Phases R-HSA-453274",
            "S Phase R-HSA-69242",
        ],
        "predicted_direction": "down",
    },
    "E2F_MYC_cellcycle": {
        "members": [
            "E2F Targets",
            "Myc Targets V1",
            "Myc Targets V2",
        ],
        "predicted_direction": "down",
    },
    "EMT": {
        "members": [
            "Epithelial Mesenchymal Transition",
        ],
        "predicted_direction": "up",  # surviving tissue is mesenchymal-enriched
    },
}

# sanity check all members exist in ssGSEA columns
all_members = [m for d in FACTOR_DEF.values() for m in d["members"]]
missing = [m for m in all_members if m not in ssgsea.columns]
assert not missing, f"missing ssGSEA columns: {missing}"

# ---------------------------------------------------------------------------
# 3) Build per-subject factor values (composite = mean of member z-scores)
# ---------------------------------------------------------------------------
# First z-score each pathway across ALL samples in ssgsea (not just paired)
# to get a stable scale, then mean the members.
z = ssgsea.set_index("sample_id").copy()
z = (z - z.mean()) / z.std(ddof=0)
z = z.reset_index()

def composite(df, members):
    return df[members].mean(axis=1)

long_rows = []   # (subject, response, timepoint, factor, value)
for _, row in inv_p.iterrows():
    sid, subj, tp, resp = row.sample_id, row.subject_id, row.timepoint, row.response_bin
    zrow = z[z.sample_id == sid]
    if zrow.empty:
        continue
    for fname, fdef in FACTOR_DEF.items():
        val = zrow[fdef["members"]].mean(axis=1).values[0]
        long_rows.append((subj, resp, tp, fname, "composite", val))
        # also log each member individually
        for mem in fdef["members"]:
            long_rows.append((subj, resp, tp, fname, mem, zrow[mem].values[0]))

long_df = pd.DataFrame(long_rows, columns=["subject_id", "response_bin",
                                           "timepoint", "factor", "member",
                                           "value"])

# wide: pre vs post per subject
wide = (long_df
        .pivot_table(index=["subject_id", "response_bin", "factor", "member"],
                     columns="timepoint", values="value")
        .reset_index())
wide["delta"] = wide["post"] - wide["pre"]

wide.to_csv(f"{OUT}/baseline_factor_per_subject_delta.tsv", sep="\t", index=False)
print(f"wrote per-subject deltas: n={len(wide)} rows ({wide.factor.nunique()} factors)")

# ---------------------------------------------------------------------------
# 4) Statistical tests --- per factor, for composite + each member
# ---------------------------------------------------------------------------
def tests_for_group(deltas, predicted_dir):
    """magnitude = one-sample Wilcoxon of Δ; direction = binomial sign test."""
    arr = np.asarray(deltas, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) < 2:
        return dict(n=len(arr), median=np.nan, wilcox_p=np.nan,
                    n_predicted=np.nan, n_total=len(arr), sign_p=np.nan)
    # wilcoxon: two-sided against zero
    try:
        w_p = stats.wilcoxon(arr).pvalue
    except ValueError:
        w_p = np.nan
    # sign test: count subjects moving in predicted direction, binomial 0.5
    if predicted_dir == "down":
        n_pred = int((arr < 0).sum())
    else:
        n_pred = int((arr > 0).sum())
    n_total = int((arr != 0).sum())
    if n_total > 0:
        # one-sided binomial: probability of >= n_pred successes under p=0.5
        sign_p_one = stats.binom.sf(n_pred - 1, n_total, 0.5)
    else:
        sign_p_one = np.nan
    return dict(n=len(arr), median=float(np.median(arr)), wilcox_p=w_p,
                n_predicted=n_pred, n_total=n_total, sign_p_one=sign_p_one)


def tests_between(deltas_good, deltas_bad, predicted_dir):
    """magnitude = Mann-Whitney on Δ; direction = Fisher on predicted-sign 2x2."""
    g = np.asarray(deltas_good, dtype=float); g = g[~np.isnan(g)]
    b = np.asarray(deltas_bad, dtype=float);  b = b[~np.isnan(b)]
    # MW magnitude
    if len(g) >= 1 and len(b) >= 1:
        mw_p = stats.mannwhitneyu(g, b, alternative="two-sided").pvalue
    else:
        mw_p = np.nan
    # Fisher on directional 2x2: rows=group, cols=predicted/not
    def count_pred(arr):
        if predicted_dir == "down":
            return int((arr < 0).sum()), int((arr >= 0).sum())
        return int((arr > 0).sum()), int((arr <= 0).sum())
    gp, gn = count_pred(g); bp, bn = count_pred(b)
    try:
        _, fisher_p = stats.fisher_exact([[gp, gn], [bp, bn]])
    except Exception:
        fisher_p = np.nan
    return dict(mw_p=mw_p,
                good_predicted=gp, good_not=gn,
                bad_predicted=bp, bad_not=bn,
                fisher_p=fisher_p)


stats_rows = []
for fname, fdef in FACTOR_DEF.items():
    for lvl, label in [("composite", "composite")] + \
                      [(m, m) for m in fdef["members"]]:
        sub = wide[(wide.factor == fname) & (wide.member == lvl)]
        g = sub[sub.response_bin == "good"]["delta"].values
        b = sub[sub.response_bin == "bad"]["delta"].values
        pdir = fdef["predicted_direction"]
        row = dict(factor=fname, level=lvl, predicted_direction=pdir)
        # within-good
        row |= {f"good_{k}": v for k, v in tests_for_group(g, pdir).items()}
        row |= {f"bad_{k}": v for k, v in tests_for_group(b, pdir).items()}
        row |= tests_between(g, b, pdir)
        stats_rows.append(row)

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(f"{OUT}/baseline_factor_pharmacodynamics_stats.tsv",
                sep="\t", index=False)
print(f"wrote stats: n={len(stats_df)} rows")

# ---------------------------------------------------------------------------
# 5) Sign-table (composite only) --- the headline for the manuscript
# ---------------------------------------------------------------------------
sign_rows = []
for fname, fdef in FACTOR_DEF.items():
    sub = wide[(wide.factor == fname) & (wide.member == "composite")]
    pdir = fdef["predicted_direction"]
    for grp in ["good", "bad"]:
        deltas = sub[sub.response_bin == grp]["delta"].values
        if pdir == "down":
            n_pred = int((deltas < 0).sum())
        else:
            n_pred = int((deltas > 0).sum())
        n_tot = len(deltas)
        sign_p = stats.binom.sf(n_pred - 1, n_tot, 0.5) if n_tot else np.nan
        sign_rows.append({
            "factor": fname, "predicted_direction": pdir, "group": grp,
            "n_predicted": n_pred, "n_total": n_tot,
            "fraction_predicted": n_pred / n_tot if n_tot else np.nan,
            "sign_binomial_one_sided_P": sign_p,
        })

sign_df = pd.DataFrame(sign_rows)
sign_df.to_csv(f"{OUT}/baseline_factor_sign_table.tsv", sep="\t", index=False)
print("wrote sign-table (composite-level):")
print(sign_df.to_string(index=False))

# ---------------------------------------------------------------------------
# 6) Headline printout --- composite-only summary
# ---------------------------------------------------------------------------
print("\n=== Composite-level summary (paired n=6+6) ===")
comp = stats_df[stats_df.level == "composite"].copy()
cols = ["factor", "predicted_direction",
        "good_median", "good_wilcox_p", "good_n_predicted",
        "bad_median", "bad_wilcox_p", "bad_n_predicted",
        "mw_p", "fisher_p"]
print(comp[cols].to_string(index=False))
