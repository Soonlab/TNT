"""
Compute v2 9 × 4 convergence test using the user-specified cascade Δ list
(Δ SBS5, Δ MHC-I neoantigen binders, Δ Treg, Δ IGH clonotype count).

Reproduces the manuscript's analysis machinery (script 09_targeted_convergence_test.py):
- partial Spearman with response-group adjustment via rank residualisation
- plain Spearman for the manuscript-quoted headline value (r = −0.07, P = 0.83)
- BH-correction across the 36 pairs (separately for plain & partial)

Outputs:
- tables/convergence_36pair_used_v2.tsv   (long form: baseline × cascade × n × plain r/P × partial r/P × BH q)
- tables/sanity_check_headline.tsv        (DSB × CD8_cytotoxic_delta sanity row, computed fresh)
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
OUT_TBL = ROOT / "tables"
OUT_TBL.mkdir(parents=True, exist_ok=True)

ANALYSIS = Path("/mnt/sda1/data/TNT/analysis")
MASTER = ANALYSIS / "260418_add" / "integrated_subject_master_v2.tsv"
LONG   = ANALYSIS / "09_integration" / "paired_delta" / "paired_feature_long.tsv"
DNEW   = ANALYSIS / "260418_add" / "paired_immune_delta_per_subject.tsv"

# --- 9 baselines (same as the manuscript's 36-pair convergence test) ---
BASE = [
    "DNA Double-Strand Break Repair R-HSA-5693532",
    "DNA Repair R-HSA-73894",
    "HDR Thru Homologous Recombination (HRR) R-HSA-5685942",
    "E2F Targets",
    "G2-M Checkpoint",
    "Myc Targets V2",
    "MHC_II",
    "MSI_pct",
    "frac_amp",
]

# --- 4 cascade Δ features (user-specified v2 list) ---
CASC_V2 = ["SBS5_delta", "neo_binders_delta", "Treg_delta", "IGH_n_delta"]

# --- load master ---
M = pd.read_csv(MASTER, sep="\t")
M["subject_id"] = M["subject_id"].astype(str)

# --- build paired Δ table (pivot of paired_feature_long; identical to script 09) ---
L = pd.read_csv(LONG, sep="\t")
L["subject_id"] = L["subject_id"].astype(str)
L["delta"] = L["post"] - L["pre"]
delta_legacy = L.pivot(index="subject_id", columns="feature", values="delta")
delta_legacy.columns = [f"{c}_delta" for c in delta_legacy.columns]

# include the per-subject immune Δ table for CD8_cytotoxic_delta (sanity)
D = pd.read_csv(DNEW, sep="\t")
D["subject_id"] = D["subject_id"].astype(str)
delta_new = D.set_index("subject_id")[[c for c in D.columns if c.endswith("_delta")]]
delta_all = delta_legacy.join(delta_new, how="outer")

paired_subjects = sorted(delta_all.index, key=int)
B = M.set_index("subject_id").loc[paired_subjects, BASE + ["response_bin"]].copy()
B["y_good"] = (B["response_bin"] == "good").astype(int)


def partial_spearman(x, y, z):
    """Rank-residualise x and y against z, then Pearson on residuals."""
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    if len(d) < 5:
        return np.nan, np.nan, len(d)
    rx = stats.rankdata(d.x.values)
    ry = stats.rankdata(d.y.values)
    z_ = d.z.values - d.z.values.mean()
    rx_res = rx - np.polyval(np.polyfit(z_, rx, 1), z_)
    ry_res = ry - np.polyval(np.polyfit(z_, ry, 1), z_)
    if rx_res.std() == 0 or ry_res.std() == 0:
        return np.nan, np.nan, len(d)
    r, p = stats.pearsonr(rx_res, ry_res)
    return r, p, len(d)


def both_spearman(b_col: str, c_col: str):
    x = B[b_col].values
    y = delta_all[c_col].reindex(B.index).values
    valid = ~(np.isnan(x.astype(float)) | np.isnan(y.astype(float)))
    if valid.sum() < 5:
        return None
    xs, ys = x[valid].astype(float), y[valid].astype(float)
    plain_r, plain_p = stats.spearmanr(xs, ys)
    partial_r, partial_p, n_part = partial_spearman(xs, ys, B["y_good"].values[valid].astype(float))
    return {
        "baseline_feature": b_col,
        "cascade_feature": c_col,
        "n": int(valid.sum()),
        "spearman_r": round(float(plain_r), 4),
        "spearman_p": round(float(plain_p), 4),
        "partial_r": round(float(partial_r), 4) if not np.isnan(partial_r) else np.nan,
        "partial_P": round(float(partial_p), 4) if not np.isnan(partial_p) else np.nan,
    }


# ------------------------------------------------------------------
# (1) Sanity-check headline pair (DSB × CD8_cytotoxic_delta)
# ------------------------------------------------------------------
print("=== Sanity check: headline pair (DSB × CD8_cytotoxic_delta) ===")
sanity_row = both_spearman(
    "DNA Double-Strand Break Repair R-HSA-5693532", "CD8_cytotoxic_delta"
)
print(f"  n             : {sanity_row['n']}")
print(f"  plain Spearman: r = {sanity_row['spearman_r']:+.3f}, P = {sanity_row['spearman_p']:.3f}")
print(f"  partial Spear : r = {sanity_row['partial_r']:+.3f}, P = {sanity_row['partial_P']:.3f}")
print("  manuscript    : r = −0.07, P = 0.83  (plain) — should match plain row above")
print()
pd.DataFrame([sanity_row]).to_csv(OUT_TBL / "sanity_check_headline.tsv", sep="\t", index=False)

# ------------------------------------------------------------------
# (2) Full 9 × 4 v2 convergence grid (user-specified cascades)
# ------------------------------------------------------------------
rows = []
for bf in BASE:
    for cf in CASC_V2:
        r = both_spearman(bf, cf)
        if r is not None:
            rows.append(r)
R = pd.DataFrame(rows)
# BH on plain and partial separately
R["BH_q_plain"] = multipletests(R["spearman_p"].fillna(1.0), method="fdr_bh")[1].round(4)
R["BH_q_partial"] = multipletests(R["partial_P"].fillna(1.0), method="fdr_bh")[1].round(4)
R = R.sort_values("spearman_p")
R.to_csv(OUT_TBL / "convergence_36pair_used_v2.tsv", sep="\t", index=False)

print(f"=== v2 36-pair convergence test ({len(R)} pairs) ===")
print(f"  n range            : {R['n'].min()} – {R['n'].max()}")
print(f"  plain r range      : [{R['spearman_r'].min():+.3f}, {R['spearman_r'].max():+.3f}]")
print(f"  plain P range      : [{R['spearman_p'].min():.3f}, {R['spearman_p'].max():.3f}]")
print(f"  partial r range    : [{R['partial_r'].min():+.3f}, {R['partial_r'].max():+.3f}]")
print(f"  partial P range    : [{R['partial_P'].min():.3f}, {R['partial_P'].max():.3f}]")
print(f"  plain P < 0.05     : {(R['spearman_p'] < 0.05).sum()}/{len(R)}")
print(f"  partial P < 0.05   : {(R['partial_P'] < 0.05).sum()}/{len(R)}")
print(f"  BH q (plain) <0.05 : {(R['BH_q_plain'] < 0.05).sum()}/{len(R)}")
print(f"  BH q (partial)<0.05: {(R['BH_q_partial'] < 0.05).sum()}/{len(R)}")
print()
print(R.to_string(index=False))
