#!/usr/bin/env python
"""Build formalized supplementary tables S12-S14 for v0.7.7 from NanoString TSVs."""
from pathlib import Path
import pandas as pd
import numpy as np

SRC = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring/tables")
OUT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring/manuscript/tables")
OUT.mkdir(parents=True, exist_ok=True)


# ========== Table S12: Pre-spec Arrow 5 rescue result ==========
p1p4 = pd.read_csv(SRC / "P1_P4_primary.tsv", sep="\t")
s1s3 = pd.read_csv(SRC / "S1_S3_lineage.tsv", sep="\t")
s4_sum = pd.read_csv(SRC / "S4_platform_summary.json" if False else SRC / "S4_platform_concordance.tsv",
                     sep="\t")
t1t3 = pd.read_csv(SRC / "T1_T3_cascade.tsv", sep="\t")

# Build narrative Table S12
rows = []
tier_cols = ["target","n_good","n_bad","good_median","bad_median","good_mean","bad_mean",
             "MW_P_1s_good_gt_bad","MW_P_2s","Wgood_Wilcoxon_P_2s","Wbad_Wilcoxon_P_2s"]
for _, r in p1p4.iterrows():
    rows.append(dict(tier="Primary", test=r["target"],
                      good_median=round(r["good_median"], 3), bad_median=round(r["bad_median"], 3),
                      good_mean=round(r["good_mean"], 3), bad_mean=round(r["bad_mean"], 3),
                      direction=("good>bad" if r["good_mean"] > r["bad_mean"] else "bad>good"),
                      MW_1sP_good_gt_bad=round(r["MW_P_1s_good_gt_bad"], 3),
                      MW_2sP=round(r["MW_P_2s"], 3),
                      Wilcoxon_good_2sP=round(r["Wgood_Wilcoxon_P_2s"], 3),
                      Wilcoxon_bad_2sP=round(r["Wbad_Wilcoxon_P_2s"], 3)))
for _, r in s1s3.iterrows():
    rows.append(dict(tier="Secondary S1-S3", test=r["target"],
                      good_median=round(r["good_median"], 3), bad_median=round(r["bad_median"], 3),
                      good_mean=round(r["good_mean"], 3), bad_mean=round(r["bad_mean"], 3),
                      direction=("good>bad" if r["good_mean"] > r["bad_mean"] else "bad>good"),
                      MW_1sP_good_gt_bad=round(r["MW_P_1s_good_gt_bad"], 3),
                      MW_2sP=round(r["MW_P_2s"], 3),
                      Wilcoxon_good_2sP=round(r["Wgood_Wilcoxon_P_2s"], 3),
                      Wilcoxon_bad_2sP=round(r["Wbad_Wilcoxon_P_2s"], 3)))
# Append tertiary T1-T3 + S4 summary as footer rows
for _, r in t1t3.iterrows():
    rows.append(dict(tier="Tertiary cascade", test=r["test"],
                      good_median=np.nan, bad_median=np.nan,
                      good_mean=np.nan, bad_mean=np.nan, direction="-",
                      MW_1sP_good_gt_bad=np.nan, MW_2sP=np.nan,
                      Wilcoxon_good_2sP=round(r["pearson_r"], 3),  # repurpose field
                      Wilcoxon_bad_2sP=round(r["pearson_P"], 3)))   # repurpose field

df_s12 = pd.DataFrame(rows)
df_s12.to_csv(OUT / "TableS12_NanoString_prespec_Arrow5_rescue.tsv",
              sep="\t", index=False, float_format="%.4g")
print(f"TableS12 written ({len(df_s12)} rows)")

# ========== Table S13: Exploratory pre/post/Δ composite MW ==========
pre = pd.read_csv(SRC / "v2_pre_MW.tsv", sep="\t")
post = pd.read_csv(SRC / "v2_post_MW.tsv", sep="\t")
delta = pd.read_csv(SRC / "v2_delta_MW.tsv", sep="\t")
ratios = pd.read_csv(SRC / "v2_ratios_MW.tsv", sep="\t")

def tidy(df, axis_label):
    return pd.DataFrame({
        "composite_or_ratio": df["composite"],
        "axis": axis_label,
        "good_mean": df["good_mean"].round(3),
        "bad_mean": df["bad_mean"].round(3),
        "good_minus_bad": (df["good_mean"] - df["bad_mean"]).round(3),
        "direction": df["direction"],
        "MW_1sP_good_gt_bad": df["MW_P_1s_good_gt_bad"].round(3),
        "MW_1sP_good_lt_bad": df["MW_P_1s_good_lt_bad"].round(3),
        "MW_2sP": df["MW_P_2s"].round(3),
        "ceiling_hit_good_gt_bad_le_0p05": (df["MW_P_1s_good_gt_bad"] <= 0.05),
        "good_values": df["good_values"],
        "bad_values": df["bad_values"],
    })

s13_parts = [tidy(pre, "pre-treatment"),
             tidy(post, "post-treatment"),
             tidy(delta, "Δ (post − pre)")]
# Ratios: tidy to same schema
ratios_tidy = pd.DataFrame({
    "composite_or_ratio": ratios["composite"],
    "axis": ratios["tag"].str.replace("_ratio", " ratio"),
    "good_mean": ratios["good_mean"].round(3),
    "bad_mean": ratios["bad_mean"].round(3),
    "good_minus_bad": (ratios["good_mean"] - ratios["bad_mean"]).round(3),
    "direction": ratios["direction"],
    "MW_1sP_good_gt_bad": ratios["MW_P_1s_good_gt_bad"].round(3),
    "MW_1sP_good_lt_bad": ratios["MW_P_1s_good_lt_bad"].round(3),
    "MW_2sP": ratios["MW_P_2s"].round(3),
    "ceiling_hit_good_gt_bad_le_0p05": (ratios["MW_P_1s_good_gt_bad"] <= 0.05),
    "good_values": ratios["good_values"],
    "bad_values": ratios["bad_values"],
})
s13_parts.append(ratios_tidy)

df_s13 = pd.concat(s13_parts, ignore_index=True)
df_s13.to_csv(OUT / "TableS13_NanoString_exploratory_pre_post_delta.tsv",
              sep="\t", index=False)
print(f"TableS13 written ({len(df_s13)} rows)")

# ========== Table S14: IBI vs IAE descriptive fingerprint ==========
iae = pd.read_csv(SRC / "v2_IAE_vs_IBI_descriptive.tsv", sep="\t")
genes = pd.read_csv(SRC / "v2_IAE_vs_IBI_gene_descriptive.tsv", sep="\t").head(40)

# composite section
df_s14a = iae.rename(columns={
    "feature":"feature_or_gene", "source":"kind",
    "IAE_mean":"IAE_n2_mean", "IBI_mean":"IBI_n3_mean",
    "IAE_minus_IBI":"IAE_minus_IBI",
    "IAE_values":"IAE_values", "IBI_values":"IBI_values"
})[["feature_or_gene","kind","IAE_n2_mean","IBI_n3_mean","IAE_minus_IBI",
    "IAE_values","IBI_values"]].round(3)

# gene section
df_s14b = pd.DataFrame({
    "feature_or_gene": genes["gene"],
    "kind": "gene",
    "IAE_n2_mean": genes["IAE_mean"].round(3),
    "IBI_n3_mean": genes["IBI_mean"].round(3),
    "IAE_minus_IBI": genes["IAE_minus_IBI"].round(3),
    "IAE_values": "",
    "IBI_values": "",
})

df_s14 = pd.concat([df_s14a, df_s14b], ignore_index=True)
df_s14.to_csv(OUT / "TableS14_IBI_vs_IAE_fingerprint.tsv", sep="\t", index=False)
print(f"TableS14 written ({len(df_s14)} rows: {len(df_s14a)} composites/ratios + {len(df_s14b)} top genes)")

# ========== Subject deep-dive as TableS15 ==========
dd = pd.read_csv(SRC / "v2_subject_deepdive.tsv", sep="\t")
pheno = pd.read_csv(SRC / "v2_phenotype_classification.tsv", sep="\t")
dd_full = dd.merge(pheno, on="subject", how="left", suffixes=("","_pheno"))
dd_full.to_csv(OUT / "TableS15_subject_deepdive_NanoString.tsv", sep="\t", index=False,
               float_format="%.3f")
print(f"TableS15 written ({len(dd_full)} rows)")

print("\nAll 4 tables written to:", OUT)
