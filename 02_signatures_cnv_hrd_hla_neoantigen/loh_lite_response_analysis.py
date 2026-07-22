#!/usr/bin/env python3
"""Annotate LOH-lite HLA LOH calls with response and compute per-subject summary."""
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact, mannwhitneyu

R = Path("/mnt/sda1/data/TNT/analysis")
loh = pd.read_csv(R/"03_hla/loh_lite/hla_loh_lite_results.tsv", sep="\t")
meta = pd.read_excel("/mnt/sda1/data/TNT/TNT_WES/meta_WES.xlsx")

# map subject -> response_bin
resp = meta[["subject_id","response_bin"]].drop_duplicates().set_index("subject_id")["response_bin"]
loh["response_bin"] = loh["subject_id"].map(resp)

# keep tumor samples only (those with LOH_call column meaningful; skip where allele2 blank = homozygous)
het = loh[loh["allele2"].fillna("").astype(str).str.strip() != ""].copy()

# per-subject any-locus LOH (take pre-treatment if available)
pre = het[het["sample"].str.contains("-P|-PR", regex=True) & ~het["sample"].str.contains("PO")]
per_subj = pre.groupby(["subject_id","response_bin"]).agg(
    n_het_loci=("LOH_call","size"),
    n_loh=("LOH_call", lambda x: int(x.sum())),
    min_tumor_ratio=("tumor_ratio","min"),
    min_fisher_p=("fisher_p","min"),
).reset_index()
per_subj["any_LOH"] = per_subj["n_loh"]>0
per_subj.to_csv(R/"03_hla/loh_lite/per_subject_pre_LOH.tsv", sep="\t", index=False)

# Fisher: any_LOH vs response
tab = pd.crosstab(per_subj["response_bin"], per_subj["any_LOH"])
print("Any-locus HLA LOH (pre) vs response:")
print(tab)
if tab.shape == (2,2):
    odds, p = fisher_exact(tab.values)
    print(f"  Fisher OR={odds:.2f}, p={p:.3f}")

# Mann-Whitney on min_tumor_ratio (lower = more LOH)
g = per_subj[per_subj["response_bin"]=="good"]["min_tumor_ratio"].dropna()
b = per_subj[per_subj["response_bin"]=="bad"]["min_tumor_ratio"].dropna()
if len(g)>1 and len(b)>1:
    u,p = mannwhitneyu(g,b,alternative="two-sided")
    print(f"min_tumor_ratio good(n={len(g)}) med={g.median():.3f} vs bad(n={len(b)}) med={b.median():.3f}  p={p:.3f}")

# per-locus breakdown
print("\nPer-locus LOH by response (pre):")
loc = pre.groupby(["locus","response_bin"])["LOH_call"].agg(["sum","size"]).reset_index()
print(loc.to_string(index=False))

# save summary
with open(R/"03_hla/loh_lite/LOH_response_summary.txt","w") as f:
    f.write(f"Any-locus HLA LOH (pre-treatment) vs TNT response\n")
    f.write(tab.to_string()+"\n")
    if tab.shape==(2,2):
        f.write(f"Fisher exact OR={odds:.3f} p={p:.4f}\n")
    f.write(f"\nmin tumor_ratio: good median={g.median():.3f}, bad median={b.median():.3f}\n")
print("\nwrote per_subject_pre_LOH.tsv and LOH_response_summary.txt")
