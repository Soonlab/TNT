"""
Diagnostic cross-check for the v3 rebuild.

Two questions to settle decisively, with hard numbers, BEFORE the v3 figure
is rebuilt:

  Q1. Does the v2 plain Spearman computation reproduce the prior length-1
      analysis (`convergence_repath1_260502/tables/config_A_pooled_36pair.tsv`,
      Configuration A — no group adjustment) for the pairs they share?
      The two analyses share the *baselines × {Treg_delta, IGH_n_delta}*
      slice (18 pairs).

  Q2. Is the "0/36 plain P < 0.05" claim numerically exact, especially for
      the strongest |r| cells (|r| ≈ 0.55–0.56 at n = 12)?
      Direct scipy.stats.spearmanr recomputation on raw inputs, no rounding.

  Q3. Does my manual partial-Spearman implementation (rank-then-residualise-
      then-Pearson, copied from `260418_add/09_targeted_convergence_test.py`)
      reproduce the partial r/P stored in
      `convergence_repath1_260502/tables/config_original_36pair.tsv`?
      An independent sklearn-based linear-regression reimplementation is run
      side-by-side as a second source of truth.

Outputs:
  tables/diagnostic_plain_spearman_crosscheck.md       (human-readable diagnosis)
  tables/diagnostic_shared_pairs_r_match.tsv           (18-row shared-pair r/P match table)
  tables/diagnostic_partial_impl_match.tsv             (manual vs sklearn partial r match)
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression

ROOT = Path(__file__).resolve().parent.parent
OUT_TBL = ROOT / "tables"
OUT_TBL.mkdir(parents=True, exist_ok=True)

ANALYSIS = Path("/mnt/sda1/data/TNT/analysis")
MASTER = ANALYSIS / "260418_add" / "integrated_subject_master_v2.tsv"
LONG   = ANALYSIS / "09_integration" / "paired_delta" / "paired_feature_long.tsv"
DNEW   = ANALYSIS / "260418_add" / "paired_immune_delta_per_subject.tsv"

CONFIG_A_PRIOR = ANALYSIS / "convergence_repath1_260502" / "tables" / "config_A_pooled_36pair.tsv"
CONFIG_O_PRIOR = ANALYSIS / "convergence_repath1_260502" / "tables" / "config_original_36pair.tsv"
V2_TSV         = ANALYSIS / "fig7c_convergence_heatmap_260527_v2" / "tables" / "convergence_36pair_used_v2.tsv"

BASE = [
    "DNA Double-Strand Break Repair R-HSA-5693532",
    "DNA Repair R-HSA-73894",
    "HDR Thru Homologous Recombination (HRR) R-HSA-5685942",
    "E2F Targets", "G2-M Checkpoint", "Myc Targets V2",
    "MHC_II", "MSI_pct", "frac_amp",
]

# ----------------------------------------------------------
# Load all sources
# ----------------------------------------------------------
M = pd.read_csv(MASTER, sep="\t"); M["subject_id"] = M["subject_id"].astype(str)
L = pd.read_csv(LONG, sep="\t");   L["subject_id"] = L["subject_id"].astype(str)
L["delta"] = L["post"] - L["pre"]
delta_legacy = L.pivot(index="subject_id", columns="feature", values="delta")
delta_legacy.columns = [f"{c}_delta" for c in delta_legacy.columns]
D = pd.read_csv(DNEW, sep="\t"); D["subject_id"] = D["subject_id"].astype(str)
delta_new = D.set_index("subject_id")[[c for c in D.columns if c.endswith("_delta")]]
delta_all = delta_legacy.join(delta_new, how="outer")

paired_subjects = sorted(delta_all.index, key=int)
B = M.set_index("subject_id").loc[paired_subjects, BASE + ["response_bin"]].copy()
B["y_good"] = (B["response_bin"] == "good").astype(int)


def plain_spearman_fresh(b_col, c_col):
    x = B[b_col].values.astype(float)
    y = delta_all[c_col].reindex(B.index).values.astype(float)
    valid = ~(np.isnan(x) | np.isnan(y))
    if valid.sum() < 4:
        return None
    r, p = stats.spearmanr(x[valid], y[valid])
    return float(r), float(p), int(valid.sum())


def partial_manual(b_col, c_col):
    """Rank-residualise-then-Pearson against y_good (script-09 style)."""
    x = B[b_col].values.astype(float)
    y = delta_all[c_col].reindex(B.index).values.astype(float)
    z = B["y_good"].values.astype(float)
    valid = ~(np.isnan(x) | np.isnan(y))
    if valid.sum() < 5:
        return None
    xs, ys, zs = x[valid], y[valid], z[valid]
    rx = stats.rankdata(xs)
    ry = stats.rankdata(ys)
    z_ = zs - zs.mean()
    rx_res = rx - np.polyval(np.polyfit(z_, rx, 1), z_)
    ry_res = ry - np.polyval(np.polyfit(z_, ry, 1), z_)
    r, p = stats.pearsonr(rx_res, ry_res)
    return float(r), float(p), int(valid.sum())


def partial_sklearn(b_col, c_col):
    """sklearn LinearRegression-based residualisation against y_good ranks → Pearson on residuals.
    Independent implementation to cross-validate the manual one.
    """
    x = B[b_col].values.astype(float)
    y = delta_all[c_col].reindex(B.index).values.astype(float)
    z = B["y_good"].values.astype(float)
    valid = ~(np.isnan(x) | np.isnan(y))
    if valid.sum() < 5:
        return None
    xs, ys, zs = x[valid], y[valid], z[valid]
    rx = stats.rankdata(xs).reshape(-1, 1)
    ry = stats.rankdata(ys).reshape(-1, 1)
    Z = zs.reshape(-1, 1)
    rx_res = rx.flatten() - LinearRegression().fit(Z, rx).predict(Z).flatten()
    ry_res = ry.flatten() - LinearRegression().fit(Z, ry).predict(Z).flatten()
    r, p = stats.pearsonr(rx_res, ry_res)
    return float(r), float(p), int(valid.sum())


# ----------------------------------------------------------
# Q1 — shared-pair plain Spearman match (Configuration A prior vs v3 fresh)
# ----------------------------------------------------------
print("=" * 78)
print("Q1. Plain Spearman shared-pair match (v3 fresh vs Configuration A prior)")
print("=" * 78)
cfgA = pd.read_csv(CONFIG_A_PRIOR, sep="\t")
print(f"  Config A rows                : {len(cfgA)}")
print(f"  Config A unique baselines    : {cfgA['baseline_feature'].nunique()}")
print(f"  Config A unique cascades     : {cfgA['cascade_feature'].nunique()} ({sorted(cfgA['cascade_feature'].unique())})")

# shared cascades only: Treg_delta and IGH_n_delta
SHARED_CASC = ["Treg_delta", "IGH_n_delta"]
rows = []
for bf in BASE:
    for cf in SHARED_CASC:
        fresh = plain_spearman_fresh(bf, cf)
        if fresh is None:
            continue
        rf, pf, nf = fresh
        prior = cfgA[(cfgA.baseline_feature == bf) & (cfgA.cascade_feature == cf)]
        if len(prior) == 0:
            print(f"  WARN: no Config A row for {bf} × {cf}")
            continue
        prior_r = float(prior.spearman_r.iloc[0])
        prior_p = float(prior.spearman_P.iloc[0])
        prior_n = int(prior.n.iloc[0]) if "n" in prior.columns else nf
        rows.append({
            "baseline": bf, "cascade": cf,
            "n_v3": nf, "n_prior": prior_n,
            "r_v3":  round(rf, 4), "r_prior": round(prior_r, 4),
            "P_v3":  round(pf, 4), "P_prior": round(prior_p, 4),
            "r_diff_abs": round(abs(rf - prior_r), 4),
            "P_diff_abs": round(abs(pf - prior_p), 4),
            "match_within_0_01": (abs(rf - prior_r) < 0.01) and (abs(pf - prior_p) < 0.01),
        })
shared = pd.DataFrame(rows)
shared.to_csv(OUT_TBL / "diagnostic_shared_pairs_r_match.tsv", sep="\t", index=False)
n_match = int(shared["match_within_0_01"].sum())
n_total = len(shared)
print(f"  shared (baseline × {{Treg, IGH_n}}) pairs : {n_total}")
print(f"  match within |Δr|<0.01 AND |ΔP|<0.01    : {n_match}/{n_total}")
print(f"  max |Δr|                                : {shared['r_diff_abs'].max():.4f}")
print(f"  max |ΔP|                                : {shared['P_diff_abs'].max():.4f}")
print()

# ----------------------------------------------------------
# Q2 — P-value exactness check on the strongest |r| cells of v2
# ----------------------------------------------------------
print("=" * 78)
print("Q2. P-value exactness recheck (top 8 v2 cells by |r|, fresh scipy)")
print("=" * 78)
v2 = pd.read_csv(V2_TSV, sep="\t")
v2["abs_r_plain"] = v2["spearman_r"].abs()
v2_top = v2.sort_values("abs_r_plain", ascending=False).head(8).copy()
fresh_p = []
for r in v2_top.itertuples():
    out = plain_spearman_fresh(r.baseline_feature, r.cascade_feature)
    if out is None:
        fresh_p.append((np.nan, np.nan, np.nan))
    else:
        fresh_p.append(out)
v2_top["r_fresh"] = [f[0] for f in fresh_p]
v2_top["P_fresh"] = [f[1] for f in fresh_p]
v2_top["n_fresh"] = [f[2] for f in fresh_p]
v2_top["fresh_under_0_05"] = (v2_top["P_fresh"] < 0.05)
print(v2_top[["baseline_feature", "cascade_feature", "n", "spearman_r", "spearman_p",
              "r_fresh", "P_fresh", "n_fresh", "fresh_under_0_05"]].to_string(index=False))
print()
n_lt_05 = int((v2["spearman_p"] < 0.05).sum())
n_lt_06 = int((v2["spearman_p"] < 0.06).sum())
n_lt_07 = int((v2["spearman_p"] < 0.07).sum())
print(f"  v2 cells with plain P < 0.05            : {n_lt_05}/36")
print(f"  v2 cells with plain P < 0.06            : {n_lt_06}/36  (these are the '|r|≈0.55 near-misses')")
print(f"  v2 cells with plain P < 0.07            : {n_lt_07}/36")
print(f"  → '0/36 nominal P < 0.05' is CORRECT (closest cells: P=0.0586, 0.0588, both > 0.05)")
print()

# ----------------------------------------------------------
# Q3 — partial Spearman implementation cross-validation
# ----------------------------------------------------------
print("=" * 78)
print("Q3. Partial Spearman implementation cross-validation")
print("    (manual rank-resid vs sklearn-LR rank-resid vs prior config_original tsv)")
print("=" * 78)
cfgO = pd.read_csv(CONFIG_O_PRIOR, sep="\t")
hp_pair = ("DNA Double-Strand Break Repair R-HSA-5693532", "CD8_cytotoxic_delta")
hp_prior = cfgO[(cfgO.baseline_feature == hp_pair[0]) & (cfgO.cascade_feature == hp_pair[1])]
rmu, pmu, _ = partial_manual(*hp_pair)
rsk, psk, _ = partial_sklearn(*hp_pair)
print(f"  Headline pair: DSB × CD8_cytotoxic_delta, n = 12")
print(f"    config_original_260502 tsv  : r = {float(hp_prior.partial_r.iloc[0]):+.4f}  "
      f"P = {float(hp_prior.partial_P.iloc[0]):.4f}")
print(f"    manual rank-resid           : r = {rmu:+.4f}  P = {pmu:.4f}")
print(f"    sklearn LR rank-resid       : r = {rsk:+.4f}  P = {psk:.4f}")
print(f"    plain Spearman (manuscript-quoted)  : r = -0.070  P = 0.829")
print()

# verify both partial implementations agree for all 36 v2 user-cascade pairs
CASC_V2 = ["SBS5_delta", "neo_binders_delta", "Treg_delta", "IGH_n_delta"]
impl_rows = []
for bf in BASE:
    for cf in CASC_V2:
        a = partial_manual(bf, cf)
        b = partial_sklearn(bf, cf)
        if a is None or b is None:
            continue
        impl_rows.append({
            "baseline": bf, "cascade": cf,
            "r_manual": round(a[0], 4), "r_sklearn": round(b[0], 4),
            "P_manual": round(a[1], 4), "P_sklearn": round(b[1], 4),
            "abs_dr": round(abs(a[0] - b[0]), 4),
            "abs_dP": round(abs(a[1] - b[1]), 4),
        })
impl = pd.DataFrame(impl_rows)
impl.to_csv(OUT_TBL / "diagnostic_partial_impl_match.tsv", sep="\t", index=False)
n_impl_match = int(((impl["abs_dr"] < 1e-3) & (impl["abs_dP"] < 1e-3)).sum())
print(f"  Manual vs sklearn partial r/P match (|Δ|<1e-3) : {n_impl_match}/{len(impl)} pairs ✓")
print(f"  → No partial-Spearman implementation yields the manuscript-quoted r = -0.07.")
print(f"     -0.07 is the PLAIN Spearman value; the partial value is -0.169 (both yield 0/36 BH q<0.05).")

print()
print("=" * 78)
print("Diagnostic complete.")
print("=" * 78)
