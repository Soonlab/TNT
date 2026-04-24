#!/usr/bin/env python
"""
Exploratory round A-F (H):
A. Pre-treatment baseline (in-house complement of external LC-CRT Bcell Fisher P=0.014)
B. Canonical clinical signatures (Ayers TIS, IFN-γ 6/10, IMPRES)
C. Inflamed-but-ineffective dissection (subj 10,11 vs 2,14 subgroup)
D. Subject-level deep-dive (subj 4 atypical good, subj 11 NanoString-only)
E. Post-RT absolute state (post-only MW)
F. Checkpoint/exhaustion post-RT landscape
H. HLA-I machinery per-subject heatmap (visualization only)

G (HK QC) skipped: nSolver output already housekeeping-normalized, no HK probes in output.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring")
TABLES = ROOT / "tables"
INPUT = Path("/mnt/sda1/data/TNT/analysis/ncounter_immune_score.xlsx")

COHORT = [
    (2,  "TNT RNA 5",  "TNT RNA 6",  "good"),
    (4,  "TNT RNA 11", "TNT RNA 12", "good"),
    (14, "TNT RNA 41", "TNT RNA 42", "good"),
    (10, "TNT RNA 29", "TNT RNA 30", "bad"),
    (11, "TNT RNA 32", "TNT RNA 33", "bad"),
    (13, "TNT RNA 38", "TNT RNA 39", "bad"),
]

# ---------- Expanded composite panel (verified against 730-probe panel) ----------
COMPOSITES = {
    # Already in pre-spec (recomputed for completeness across pre/post/Δ)
    "TLS_8":        ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL","LAMP3","BCL6"],
    "Plasma_proxy": ["TNFRSF17","CD38","IRF4"],
    "GC_TF":        ["BCL6","AICDA","POU2AF1"],
    "Naive_B":      ["MS4A1","CD19","CD22","PAX5"],
    "Memory_B":     ["CD27","CD79B","TNFRSF13B"],
    "BAFF_APRIL":   ["TNFSF13","TNFSF13B","TNFRSF13B","TNFRSF13C","TNFRSF17"],
    "HLA_II":       ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"],
    "CD8_exh":      ["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"],
    "HLA_I_axis":   ["NLRC5","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"],
    "T_cell":       ["CD3D","CD3E","CD3G","CD8A","CD8B"],

    # B. Canonical signatures (regulatory/clinical standards)
    "Ayers_TIS":     ["CD274","CXCR6","TIGIT","CD27","CD8A","IDO1","LAG3","PDCD1LG2","PSMB10",
                      "STAT1","IFNG","HLA-E","CMKLR1","CXCL9","CCL5","HLA-DQA1"],  # 16/18 (NKG7, HLA-DRB1 drop)
    "IFNg_6":        ["IFNG","STAT1","CXCL9","CXCL10","GZMB","HLA-A"],
    "IFNg_10_Ayers": ["IFNG","STAT1","CCR5","CXCL9","CXCL10","CXCL11","GZMA","GZMB","HLA-DRA","PRF1"],
    "IMPRES_pos":    ["PDCD1","CD28","TNFRSF4","TNFRSF14","CTLA4","TNFRSF9","CD80","CD86","CD40",
                      "ICOS","IDO1","LAG3","TIGIT","CD27"],  # 14/15 (HLA-DRB1 drop)

    # C. Effector/suppressor balance composites
    "Teff_cytotoxic": ["GZMA","GZMB","GZMH","GZMK","PRF1","IFNG","CD8A","CD8B","GNLY"],  # 9/10
    "Treg":          ["FOXP3","IL10","TGFB1","TGFB2","CTLA4","TNFRSF4","TNFRSF18"],  # 7/8
    "CD8_cytotoxic": ["GZMA","GZMB","GZMH","GZMK","PRF1","IFNG","CD8A","CD8B"],  # for cyt/exh ratio
    "M1_macro":      ["CXCL9","CXCL10","CXCL11","IL12B","CD80","CD86","IFNG","TNF","IL1B"],  # 9/10 (no NOS2)
    "M2_macro":      ["ARG1","MRC1","CD163","IL10","TGFB1","VEGFA","CCL17","CCL22","MSR1"],
    "NK_activating": ["NCR1","KLRK1","KLRD1","FCGR3A","KLRF1","KIR_Activating_Subgroup_1",
                      "KIR_Activating_Subgroup_2"],  # panel-specific names
    "NK_inhibiting": ["KLRC1","KIR3DL1","KIR3DL2","KIR3DL3","KIR_Inhibiting_Subgroup_1",
                      "KIR_Inhibiting_Subgroup_2"],
    "DC_mature":     ["LAMP3","CCR7","CD83","IRF7","IL12B"],  # 5/6 (no FSCN1)
    "HLA_I_machinery_narrow": ["NLRC5","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"],
}

# ---------- Derived ratio composites (computed from composites above) ----------
RATIOS = [
    ("Teff_over_Treg",       "Teff_cytotoxic", "Treg"),
    ("CD8cyt_over_CD8exh",   "CD8_cytotoxic",  "CD8_exh"),
    ("M1_over_M2",           "M1_macro",       "M2_macro"),
    ("NK_activ_over_inhib",  "NK_activating",  "NK_inhibiting"),
]


# ================= helpers =================
def log2z(mat: pd.DataFrame) -> pd.DataFrame:
    x = np.log2(mat.values + 1.0)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    return pd.DataFrame((x - mu) / sd, index=mat.index, columns=mat.columns)


def composite(zmat: pd.DataFrame, genes: list[str]):
    present = [g for g in genes if g in zmat.index]
    missing = [g for g in genes if g not in zmat.index]
    return zmat.loc[present].mean(axis=0) if len(present) >= 2 else None, present, missing


def mw_one_two_sided(g: np.ndarray, b: np.ndarray) -> dict:
    u1s = stats.mannwhitneyu(g, b, alternative="greater", method="exact")
    u2s = stats.mannwhitneyu(g, b, alternative="two-sided", method="exact")
    u1s_ls = stats.mannwhitneyu(g, b, alternative="less", method="exact")  # for flip detection
    return dict(MW_P_1s_good_gt_bad=float(u1s.pvalue),
                MW_P_1s_good_lt_bad=float(u1s_ls.pvalue),
                MW_P_2s=float(u2s.pvalue))


# ================= main =================
def main():
    raw = pd.read_excel(INPUT).set_index("Probe Name")
    sample_cols = [c for pair in [(p, q) for (_, p, q, _) in COHORT] for c in pair]
    mat = raw[sample_cols].astype(float)
    mat = mat.loc[mat.sum(axis=1) > 0]
    zmat = log2z(mat)

    # ------- Compute composite scores per sample and per-subject pre/post/Δ -------
    sample_comp_rows = []
    comp_pre = {}   # {subj: {composite: pre_z}}
    comp_post = {}
    comp_delta = {}
    comp_def_rows = []
    for name, genes in COMPOSITES.items():
        ser, present, missing = composite(zmat, genes)
        if ser is None:
            print(f"SKIP {name}: only {len(present)} genes present")
            continue
        comp_def_rows.append(dict(composite=name, n_genes=len(genes),
                                   n_present=len(present),
                                   present=",".join(present),
                                   missing=",".join(missing)))
        for s in sample_cols:
            sample_comp_rows.append(dict(sample=s, composite=name, score=float(ser[s])))
        for subj, pre, post, _ in COHORT:
            comp_pre.setdefault(subj, {})[name] = float(ser[pre])
            comp_post.setdefault(subj, {})[name] = float(ser[post])
            comp_delta.setdefault(subj, {})[name] = float(ser[post] - ser[pre])

    pd.DataFrame(comp_def_rows).to_csv(TABLES / "v2_composite_definitions.tsv", sep="\t", index=False)

    # Build wide tables
    pre_df = pd.DataFrame(comp_pre).T.sort_index()
    post_df = pd.DataFrame(comp_post).T.sort_index()
    delta_df = pd.DataFrame(comp_delta).T.sort_index()
    pre_df.index.name = "subject"; pre_df.to_csv(TABLES / "v2_composite_pre.tsv", sep="\t")
    post_df.index.name = "subject"; post_df.to_csv(TABLES / "v2_composite_post.tsv", sep="\t")
    delta_df.index.name = "subject"; delta_df.to_csv(TABLES / "v2_composite_delta.tsv", sep="\t")

    # ------- Derived ratios -------
    ratio_rows = {}
    for rname, ca, cb in RATIOS:
        if ca not in pre_df.columns or cb not in pre_df.columns:
            continue
        for tag, df_ in [("pre", pre_df), ("post", post_df), ("delta", delta_df)]:
            ratio_rows.setdefault(tag, {})[rname] = (df_[ca] - df_[cb]).to_dict()
    ratio_pre = pd.DataFrame(ratio_rows["pre"])
    ratio_post = pd.DataFrame(ratio_rows["post"])
    ratio_delta = pd.DataFrame(ratio_rows["delta"])
    ratio_pre.to_csv(TABLES / "v2_ratio_pre.tsv", sep="\t")
    ratio_post.to_csv(TABLES / "v2_ratio_post.tsv", sep="\t")
    ratio_delta.to_csv(TABLES / "v2_ratio_delta.tsv", sep="\t")

    # ------- A/E/F/etc: MW pre-only, post-only, Δ for every composite and ratio -------
    def mw_panel(df_: pd.DataFrame, tag: str) -> pd.DataFrame:
        rows = []
        good_subj = [s for s, _p, _q, bn in COHORT if bn == "good"]
        bad_subj = [s for s, _p, _q, bn in COHORT if bn == "bad"]
        for col in df_.columns:
            g = df_.loc[good_subj, col].values
            b = df_.loc[bad_subj, col].values
            stats_ = mw_one_two_sided(g, b)
            rows.append(dict(composite=col, tag=tag,
                              good_mean=float(g.mean()), bad_mean=float(b.mean()),
                              good_median=float(np.median(g)), bad_median=float(np.median(b)),
                              good_values=",".join(f"{x:+.3f}" for x in g),
                              bad_values=",".join(f"{x:+.3f}" for x in b),
                              direction=("good>bad" if g.mean() > b.mean() else "bad>good"),
                              **stats_))
        return pd.DataFrame(rows)

    comp_pre_mw = mw_panel(pre_df, "pre")
    comp_post_mw = mw_panel(post_df, "post")
    comp_delta_mw = mw_panel(delta_df, "delta")
    comp_pre_mw.to_csv(TABLES / "v2_pre_MW.tsv", sep="\t", index=False)
    comp_post_mw.to_csv(TABLES / "v2_post_MW.tsv", sep="\t", index=False)
    comp_delta_mw.to_csv(TABLES / "v2_delta_MW.tsv", sep="\t", index=False)

    ratio_pre_mw = mw_panel(ratio_pre, "pre_ratio")
    ratio_post_mw = mw_panel(ratio_post, "post_ratio")
    ratio_delta_mw = mw_panel(ratio_delta, "delta_ratio")
    pd.concat([ratio_pre_mw, ratio_post_mw, ratio_delta_mw]).to_csv(
        TABLES / "v2_ratios_MW.tsv", sep="\t", index=False)

    # ------- Gene-level pre-only / post-only MW scan (expand S5 to pre + post) -------
    def gene_mw_scan(tag: str, sample_picker) -> pd.DataFrame:
        rows = []
        for g in zmat.index:
            good_vals = []
            bad_vals = []
            for subj, pre, post, bn in COHORT:
                v = zmat.loc[g, sample_picker(pre, post)]
                (good_vals if bn == "good" else bad_vals).append(float(v))
            good_vals = np.array(good_vals); bad_vals = np.array(bad_vals)
            if np.all(good_vals == bad_vals[0]) and np.all(bad_vals == bad_vals[0]) \
               and good_vals[0] == bad_vals[0]:
                continue
            try:
                u1s = stats.mannwhitneyu(good_vals, bad_vals, alternative="greater", method="exact")
                u2s = stats.mannwhitneyu(good_vals, bad_vals, alternative="two-sided", method="exact")
            except ValueError:
                continue
            rows.append(dict(gene=g, tag=tag,
                              good_mean=float(good_vals.mean()),
                              bad_mean=float(bad_vals.mean()),
                              MW_P_1s=float(u1s.pvalue), MW_P_2s=float(u2s.pvalue),
                              direction=("good>bad" if good_vals.mean() > bad_vals.mean() else "bad>good")))
        df_ = pd.DataFrame(rows).sort_values("MW_P_1s")
        # BH on one-sided
        m = len(df_)
        if m:
            sorted_p = df_["MW_P_1s"].sort_values().values
            bh = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
            bh = np.clip(bh, 0, 1)
            bh_order = {val: i for i, val in enumerate(df_["MW_P_1s"].sort_values().index)}
            df_["BH_q_1s"] = [bh[bh_order[i]] for i in df_.index]
        return df_

    pre_scan = gene_mw_scan("pre", lambda pre, post: pre)
    post_scan = gene_mw_scan("post", lambda pre, post: post)
    pre_scan.to_csv(TABLES / "v2_gene_pre_scan.tsv", sep="\t", index=False)
    post_scan.to_csv(TABLES / "v2_gene_post_scan.tsv", sep="\t", index=False)

    # ------- C: Inflamed-but-ineffective subgroup analysis -------
    # "Inflamed" = ΔAyers_TIS > 0 (regulatory-grade Tumor Inflammation Signature, broader than Teff)
    # classify subjects
    classes = []
    for subj in delta_df.index:
        d_tis = float(delta_df.loc[subj, "Ayers_TIS"])
        d_teff = float(delta_df.loc[subj, "Teff_cytotoxic"])
        inflamed = d_tis > 0
        bn = [b for s, _p, _q, b in COHORT if s == subj][0]
        classes.append(dict(subject=subj, response_bin=bn,
                             dAyers_TIS=d_tis, dTeff=d_teff,
                             inflamed=bool(inflamed),
                             phenotype=("inflamed_effective" if inflamed and bn == "good" else
                                        "inflamed_ineffective" if inflamed and bn == "bad" else
                                        "cold_effective" if not inflamed and bn == "good" else
                                        "cold_ineffective")))
    pheno_df = pd.DataFrame(classes).sort_values(["response_bin", "dAyers_TIS"])
    pheno_df.to_csv(TABLES / "v2_phenotype_classification.tsv", sep="\t", index=False)

    infl_bad = pheno_df[pheno_df["phenotype"] == "inflamed_ineffective"]["subject"].tolist()
    infl_good = pheno_df[pheno_df["phenotype"] == "inflamed_effective"]["subject"].tolist()
    print(f"Phenotypes (using ΔAyers_TIS>0 as 'inflamed'):")
    print(pheno_df.to_string(index=False))
    print(f"IAE (inflamed & good): {infl_good}, IBI (inflamed & bad): {infl_bad}")

    # IAE vs IBI descriptive comparison (requires >=2 in each group)
    ibi_iae_rows = []
    if len(infl_bad) >= 2 and len(infl_good) >= 2:
        for col in list(delta_df.columns) + list(ratio_delta.columns):
            if col in delta_df.columns:
                g = delta_df.loc[infl_good, col].values
                b = delta_df.loc[infl_bad, col].values
                src = "composite"
            else:
                g = ratio_delta.loc[infl_good, col].values
                b = ratio_delta.loc[infl_bad, col].values
                src = "ratio"
            ibi_iae_rows.append(dict(feature=col, source=src,
                                      IAE_mean=float(g.mean()), IBI_mean=float(b.mean()),
                                      IAE_minus_IBI=float(g.mean() - b.mean()),
                                      IAE_values=",".join(f"{x:+.3f}" for x in g),
                                      IBI_values=",".join(f"{x:+.3f}" for x in b),
                                      n_IAE=len(g), n_IBI=len(b)))
        ibi_iae_df = pd.DataFrame(ibi_iae_rows).sort_values(
            "IAE_minus_IBI", key=lambda s: s.abs(), ascending=False)
        ibi_iae_df.to_csv(TABLES / "v2_IAE_vs_IBI_descriptive.tsv", sep="\t", index=False)

        # gene-level IBI vs IAE
        delta_genes = {}
        for g in zmat.index:
            delta_genes[g] = {subj: float(zmat.loc[g, post] - zmat.loc[g, pre])
                              for subj, pre, post, _ in COHORT}
        dg_df = pd.DataFrame(delta_genes).T
        gene_ibi = []
        for g in dg_df.index:
            iae_vals = dg_df.loc[g, infl_good].values
            ibi_vals = dg_df.loc[g, infl_bad].values
            gene_ibi.append(dict(gene=g,
                                  IAE_mean=float(iae_vals.mean()),
                                  IBI_mean=float(ibi_vals.mean()),
                                  IAE_minus_IBI=float(iae_vals.mean() - ibi_vals.mean())))
        gibi_df = pd.DataFrame(gene_ibi).sort_values(
            "IAE_minus_IBI", key=lambda s: s.abs(), ascending=False)
        gibi_df.to_csv(TABLES / "v2_IAE_vs_IBI_gene_descriptive.tsv", sep="\t", index=False)
    else:
        print(f"WARN: IAE={len(infl_good)} IBI={len(infl_bad)} -- skipping IAE vs IBI comparison tables")

    # ------- D: Subject-level deep-dive (subj 4, 11) -------
    focus_subjects = [4, 11]
    dd_rows = []
    for subj in focus_subjects:
        bn = [b for s, _p, _q, b in COHORT if s == subj][0]
        # Subject's rank among cohort per composite
        for col in delta_df.columns:
            subj_val = delta_df.loc[subj, col]
            all_vals = delta_df[col].values
            rank = int(np.sum(all_vals <= subj_val))  # rank (1-6) low to high
            dd_rows.append(dict(subject=subj, response_bin=bn, composite=col,
                                  delta_z=float(subj_val), rank_in_6=rank,
                                  pre_z=float(pre_df.loc[subj, col]),
                                  post_z=float(post_df.loc[subj, col])))
    dd_df = pd.DataFrame(dd_rows)
    dd_df.to_csv(TABLES / "v2_subject_deepdive.tsv", sep="\t", index=False)

    # ------- Summary JSON -------
    summary = {}
    for tag, df_ in [("pre", comp_pre_mw), ("post", comp_post_mw), ("delta", comp_delta_mw)]:
        hits_1s_05 = df_[df_["MW_P_1s_good_gt_bad"] <= 0.05]
        hits_1s_ls_05 = df_[df_["MW_P_1s_good_lt_bad"] <= 0.05]
        summary[tag] = dict(
            total=int(len(df_)),
            good_gt_bad_direction=int((df_["direction"] == "good>bad").sum()),
            hits_good_gt_bad_1s_P_le_0p05=hits_1s_05["composite"].tolist(),
            hits_bad_gt_good_1s_P_le_0p05=hits_1s_ls_05["composite"].tolist(),
        )
    summary["phenotypes"] = pheno_df["phenotype"].value_counts().to_dict()
    summary["inflamed_ineffective_subjects"] = infl_bad
    summary["inflamed_effective_subjects"] = infl_good
    with open(TABLES / "v2_SUMMARY.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print("=== DONE ===")
    print("Summary:")
    print(json.dumps(summary, indent=2))
    print()
    print(f"Composites tested: {len(comp_def_rows)}")
    print(f"Pre MW hits (good>bad P<=0.05): {summary['pre']['hits_good_gt_bad_1s_P_le_0p05']}")
    print(f"Post MW hits (good>bad P<=0.05): {summary['post']['hits_good_gt_bad_1s_P_le_0p05']}")
    print(f"Δ MW hits (good>bad P<=0.05): {summary['delta']['hits_good_gt_bad_1s_P_le_0p05']}")
    print()
    print(f"Flipped hits (bad>good P<=0.05):")
    for tag in ["pre", "post", "delta"]:
        print(f"  {tag}: {summary[tag]['hits_bad_gt_good_1s_P_le_0p05']}")


if __name__ == "__main__":
    main()
