#!/usr/bin/env python3
"""
16_trust4_ighv_directional_consistency.py

Question
--------
User prior observation: paired pre->post TRUST4 repertoire analysis showed
IGHV3-7 consistently up and IGHV3-74 consistently down in good responders,
while bad responders had large but non-directional shifts. Formalize this
observation across ALL IGHV genes with sufficient coverage on the 12 paired
subjects (6 good + 6 bad), using sign-consistency as the primary statistic
(magnitude tests conflate per-subject variance with direction).

Method
------
  1. Parse TNT_RNA_{i}_airr.tsv for each of the 24 paired samples. Filter
     locus==IGH, productive==T. Extract V-gene (strip allele suffix after *).
  2. Per sample: count productive IGH clonotypes per V-gene; compute fraction
     of total productive IGH clonotypes.
  3. Per subject: Δ_fraction = fraction_post - fraction_pre.
  4. Per V-gene with >=6 subjects having non-trivial coverage (fraction>=0.001
     in >=1 of {pre, post} in at least 6/12 subjects): for each group, count
     how many subjects moved up vs down. Sign-binomial within-group,
     Fisher 2x2 between-group on up/down counts.

Output
------
    trust4_ighv_per_subject_delta.tsv   V-gene x subject Δ table
    trust4_ighv_directional_stats.tsv   per V-gene sign stats (ranked)
    trust4_ighv_focus_genes.tsv         table for IGHV3-7 + IGHV3-74 and
                                        any V-gene with good_frac>=5/6 and
                                        good-vs-bad directional Fisher P<=0.1
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/data/data/TNT/analysis"
OUT = f"{ROOT}/260418_add"
AIRR_ROOT = f"{ROOT}/06_rna_immune/trust4"

inv = pd.read_csv(f"{ROOT}/00_cohort/rna_inventory.tsv", sep="\t")
pairs = (inv[inv.timepoint.isin(["pre", "post"])]
         .groupby("subject_id")["timepoint"].nunique())
paired_subjects = sorted(pairs[pairs == 2].index.tolist())
assert len(paired_subjects) == 12, paired_subjects

inv_p = inv[(inv.subject_id.isin(paired_subjects)) &
            (inv.timepoint.isin(["pre", "post"]))].copy()

# ---------------------------------------------------------------------------
# 1) parse AIRR files, build V-gene fraction per sample
# ---------------------------------------------------------------------------
def parse_airr(sample_id):
    path = f"{AIRR_ROOT}/{sample_id}/{sample_id}_airr.tsv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", usecols=["productive", "locus", "v_call"],
                     dtype=str, low_memory=False)
    df = df[(df.productive == "T") & (df.locus == "IGH")]
    if df.empty:
        return pd.Series(dtype=float)
    df["v_gene"] = df.v_call.str.split("*").str[0]
    # allow multi-gene v_call (comma-separated); take first
    df["v_gene"] = df.v_gene.str.split(",").str[0]
    counts = df.v_gene.value_counts()
    total = counts.sum()
    return counts / total

frac_by_sample = {}
for sid in inv_p.sample_id:
    s = parse_airr(sid)
    if s is None:
        print(f"  MISSING airr: {sid}")
        continue
    frac_by_sample[sid] = s

# long table: (subject, response, timepoint, v_gene, fraction)
rows = []
for _, row in inv_p.iterrows():
    sid, subj, tp, resp = row.sample_id, row.subject_id, row.timepoint, row.response_bin
    s = frac_by_sample.get(sid)
    if s is None:
        continue
    for vg, f in s.items():
        rows.append((subj, resp, tp, vg, float(f)))

long_df = pd.DataFrame(rows, columns=["subject_id", "response_bin",
                                      "timepoint", "v_gene", "fraction"])
print(f"samples parsed: {long_df.sample_id if 'sample_id' in long_df.columns else long_df.subject_id.nunique()} subjects x ~V-genes")

# ---------------------------------------------------------------------------
# 2) per-subject Δ_fraction per V-gene (fill missing with 0)
# ---------------------------------------------------------------------------
vgenes = sorted(long_df.v_gene.unique())
# wide: subject x timepoint per v_gene
wide = (long_df
        .pivot_table(index=["subject_id", "response_bin", "v_gene"],
                     columns="timepoint", values="fraction", fill_value=0.0)
        .reset_index())
wide["delta"] = wide["post"] - wide["pre"]

wide.to_csv(f"{OUT}/trust4_ighv_per_subject_delta.tsv", sep="\t", index=False)
print(f"wrote per-subject V-gene deltas: {len(wide)} rows, {wide.v_gene.nunique()} V-genes")

# ---------------------------------------------------------------------------
# 3) coverage filter: keep V-genes where >=6 subjects have non-trivial
#    expression (fraction >= 0.001 in either pre or post) AND data is present
#    for both time points for >=5 good and >=5 bad subjects
# ---------------------------------------------------------------------------
def subject_has_nontrivial(row):
    return (row["pre"] >= 0.001) or (row["post"] >= 0.001)

wide["nontrivial"] = wide.apply(subject_has_nontrivial, axis=1)

keep = []
for vg, g in wide.groupby("v_gene"):
    n_nontrivial = g.nontrivial.sum()
    n_good = g[g.response_bin == "good"].shape[0]
    n_bad = g[g.response_bin == "bad"].shape[0]
    if n_nontrivial >= 6 and n_good >= 5 and n_bad >= 5:
        keep.append(vg)
print(f"V-genes passing coverage filter: {len(keep)} / {len(vgenes)}")

# ---------------------------------------------------------------------------
# 4) directional stats per V-gene
# ---------------------------------------------------------------------------
def sign_table(deltas):
    arr = np.asarray(deltas, dtype=float)
    n_up = int((arr > 0).sum())
    n_down = int((arr < 0).sum())
    n_zero = int((arr == 0).sum())
    n_tot = len(arr) - n_zero  # ignore zero for sign test
    # two-sided binomial: P of seeing >=max(n_up,n_down) extremes under H0=0.5
    if n_tot > 0:
        k = max(n_up, n_down)
        # exact two-sided binomial test
        sign_p_two = min(1.0,
                         2 * stats.binom.sf(k - 1, n_tot, 0.5))
    else:
        sign_p_two = np.nan
    return dict(n_up=n_up, n_down=n_down, n_zero=n_zero,
                n_nontrivial=n_tot, sign_p_two=sign_p_two,
                median_delta=float(np.median(arr)) if len(arr) else np.nan)


stat_rows = []
for vg in keep:
    sub = wide[wide.v_gene == vg]
    g = sub[sub.response_bin == "good"]["delta"].values
    b = sub[sub.response_bin == "bad"]["delta"].values

    g_sign = sign_table(g)
    b_sign = sign_table(b)

    # between-group Fisher on up/down 2x2, ignoring zeros
    fisher_tbl = [[g_sign["n_up"], g_sign["n_down"]],
                  [b_sign["n_up"], b_sign["n_down"]]]
    if min(sum(fisher_tbl[0]), sum(fisher_tbl[1])) == 0:
        fisher_p = np.nan
    else:
        _, fisher_p = stats.fisher_exact(fisher_tbl)

    # MW on Δ magnitude
    try:
        mw_p = stats.mannwhitneyu(g, b, alternative="two-sided").pvalue
    except ValueError:
        mw_p = np.nan

    # directional concordance score: fraction of good in the majority direction
    # minus fraction of bad in that direction (positive => good more coherent)
    g_tot = g_sign["n_up"] + g_sign["n_down"]
    b_tot = b_sign["n_up"] + b_sign["n_down"]
    if g_tot > 0 and b_tot > 0:
        g_maj = max(g_sign["n_up"], g_sign["n_down"]) / g_tot
        b_maj = max(b_sign["n_up"], b_sign["n_down"]) / b_tot
        coherence_gap = g_maj - b_maj
    else:
        coherence_gap = np.nan

    stat_rows.append({
        "v_gene": vg,
        "good_n_up": g_sign["n_up"], "good_n_down": g_sign["n_down"],
        "good_median_delta": g_sign["median_delta"],
        "good_sign_P_two": g_sign["sign_p_two"],
        "bad_n_up": b_sign["n_up"], "bad_n_down": b_sign["n_down"],
        "bad_median_delta": b_sign["median_delta"],
        "bad_sign_P_two": b_sign["sign_p_two"],
        "fisher_P_updown": fisher_p,
        "mw_P_delta": mw_p,
        "coherence_gap": coherence_gap,
    })

stats_df = pd.DataFrame(stat_rows).sort_values("coherence_gap", ascending=False)
stats_df.to_csv(f"{OUT}/trust4_ighv_directional_stats.tsv",
                sep="\t", index=False)

# ---------------------------------------------------------------------------
# 5) focus V-genes: IGHV3-7, IGHV3-74, plus any with good majority >= 5/6 and
#    Fisher P <= 0.2
# ---------------------------------------------------------------------------
focus_names = {"IGHV3-7", "IGHV3-74"}
auto_focus = stats_df[
    ((stats_df.good_n_up >= 5) | (stats_df.good_n_down >= 5)) &
    (stats_df.fisher_P_updown <= 0.2)
]["v_gene"].tolist()
focus = sorted(focus_names.union(auto_focus))
focus_df = stats_df[stats_df.v_gene.isin(focus)].copy().sort_values(
    ["coherence_gap", "fisher_P_updown"], ascending=[False, True])
focus_df.to_csv(f"{OUT}/trust4_ighv_focus_genes.tsv", sep="\t", index=False)

print("\n=== top V-genes by good-vs-bad directional coherence_gap ===")
show_cols = ["v_gene", "good_n_up", "good_n_down", "good_median_delta",
             "bad_n_up", "bad_n_down", "bad_median_delta",
             "fisher_P_updown", "mw_P_delta", "coherence_gap"]
print(stats_df[show_cols].head(15).to_string(index=False))

print("\n=== focus genes (user-named + auto) ===")
print(focus_df[show_cols].to_string(index=False))

# user-prior gene check
for g in ["IGHV3-7", "IGHV3-74"]:
    row = stats_df[stats_df.v_gene == g]
    if row.empty:
        print(f"\n!! {g}: not in coverage-filtered set; raw wide row:")
        raw = wide[wide.v_gene == g]
        if not raw.empty:
            print(raw.sort_values(["response_bin", "subject_id"]).to_string(index=False))
        else:
            print("   not detected in any sample")
    else:
        print(f"\n{g}: good Δ up/down = "
              f"{int(row.good_n_up.values[0])}/{int(row.good_n_down.values[0])}, "
              f"bad = {int(row.bad_n_up.values[0])}/{int(row.bad_n_down.values[0])}, "
              f"Fisher P = {row.fisher_P_updown.values[0]:.3f}")

# ---------------------------------------------------------------------------
# 6) Repertoire-level aggregate tests
#    "Is good responder's V-gene repertoire response more directionally coherent
#     than bad responder's overall?"
# ---------------------------------------------------------------------------

# per-V-gene majority fraction (coherence) in each group
def majority_frac(n_up, n_down):
    tot = n_up + n_down
    if tot == 0:
        return np.nan
    return max(n_up, n_down) / tot

stats_df["good_majority_frac"] = stats_df.apply(
    lambda r: majority_frac(r.good_n_up, r.good_n_down), axis=1)
stats_df["bad_majority_frac"] = stats_df.apply(
    lambda r: majority_frac(r.bad_n_up, r.bad_n_down), axis=1)

# Classify each V-gene by direction pattern
def classify(r):
    def dom(n_up, n_down, thr=5):
        tot = n_up + n_down
        if tot == 0:
            return "none"
        frac = max(n_up, n_down) / tot
        if frac < 0.7:
            return "mixed"
        return "up" if n_up > n_down else "down"

    g_cls = dom(r.good_n_up, r.good_n_down)
    b_cls = dom(r.bad_n_up, r.bad_n_down)
    if g_cls in ("up", "down") and b_cls == "mixed":
        return "good_coherent_bad_mixed"
    if b_cls in ("up", "down") and g_cls == "mixed":
        return "bad_coherent_good_mixed"
    if g_cls == b_cls and g_cls in ("up", "down"):
        return "both_coherent_same"
    if g_cls in ("up", "down") and b_cls in ("up", "down") and g_cls != b_cls:
        return "both_coherent_opposite"
    return "both_mixed"

stats_df["pattern"] = stats_df.apply(classify, axis=1)
pattern_counts = stats_df.pattern.value_counts()
print("\n=== V-gene directional pattern breakdown (70% majority threshold) ===")
print(pattern_counts.to_string())

# Aggregate test 1: Wilcoxon signed-rank of (good_majority - bad_majority) across V-genes
paired_gap = (stats_df.good_majority_frac - stats_df.bad_majority_frac).dropna()
if len(paired_gap) >= 6:
    w_stat, w_p = stats.wilcoxon(paired_gap, alternative="greater")
    print(f"\nAggregate Wilcoxon (good majority_frac > bad majority_frac, "
          f"paired across {len(paired_gap)} V-genes): W={w_stat:.1f}, P={w_p:.4f}")

# Aggregate test 2: binomial --- how many V-genes have good_majority > bad_majority?
n_good_higher = int((stats_df.good_majority_frac > stats_df.bad_majority_frac).sum())
n_comparable = int(stats_df.good_majority_frac.notna().sum() -
                   (stats_df.good_majority_frac == stats_df.bad_majority_frac).sum())
binom_p = stats.binom.sf(n_good_higher - 1, n_comparable, 0.5) if n_comparable else np.nan
print(f"Aggregate binomial: good_majority_frac > bad_majority_frac in "
      f"{n_good_higher}/{n_comparable} V-genes, one-sided P = {binom_p:.4f}")

# Aggregate test 3: "good coherent, bad mixed" vs opposite pattern
ratio = pattern_counts.get("good_coherent_bad_mixed", 0)
opp_ratio = pattern_counts.get("bad_coherent_good_mixed", 0)
tot_directional = ratio + opp_ratio
if tot_directional > 0:
    bp = stats.binom.sf(ratio - 1, tot_directional, 0.5)
    print(f"\ngood_coherent_bad_mixed vs bad_coherent_good_mixed = "
          f"{ratio} vs {opp_ratio}, one-sided binomial P = {bp:.4f}")

# save augmented stats
stats_df.to_csv(f"{OUT}/trust4_ighv_directional_stats.tsv",
                sep="\t", index=False)

# top "good_coherent_bad_mixed" V-genes
gcbm = stats_df[stats_df.pattern == "good_coherent_bad_mixed"].copy()
gcbm = gcbm.sort_values("fisher_P_updown")
print("\n=== V-genes classified 'good_coherent_bad_mixed' (n={}) ===".format(len(gcbm)))
print(gcbm[["v_gene", "good_n_up", "good_n_down",
            "bad_n_up", "bad_n_down",
            "fisher_P_updown", "mw_P_delta",
            "good_median_delta", "bad_median_delta"]].to_string(index=False))
gcbm.to_csv(f"{OUT}/trust4_ighv_good_coherent_bad_mixed.tsv", sep="\t", index=False)
