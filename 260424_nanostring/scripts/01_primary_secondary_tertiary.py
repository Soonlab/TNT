#!/usr/bin/env python
"""
NanoString PanCancer Immune — pre-specified analysis (P1-P4, S1-S5, T1-T3).
See ../PRE_SPEC.md for frozen hypotheses.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring")
INPUT = Path("/mnt/sda1/data/TNT/analysis/ncounter_immune_score.xlsx")
RNA_LOGTPM = Path("/mnt/sda1/data/TNT/analysis/06_rna_immune/logtpm_symbol.tsv")
TABLES = ROOT / "tables"

# ---------- Cohort table (frozen from memory + PRE_SPEC.md) ----------
COHORT = [
    # subj, pre_col, post_col, response_bin, score, rnaseq_paired
    (2,  "TNT RNA 5",  "TNT RNA 6",  "good", 0, True),
    (4,  "TNT RNA 11", "TNT RNA 12", "good", 0, True),
    (14, "TNT RNA 41", "TNT RNA 42", "good", 0, True),
    (10, "TNT RNA 29", "TNT RNA 30", "bad",  3, True),
    (11, "TNT RNA 32", "TNT RNA 33", "bad",  3, False),  # RNA-seq pre missing
    (13, "TNT RNA 38", "TNT RNA 39", "bad",  3, True),
]

# ---------- Composite definitions (frozen) ----------
COMPOSITES = {
    # Primary
    "TLS_8":          ["CXCL13", "CCL19", "CCL21", "CXCR5", "CCR7", "SELL", "LAMP3", "BCL6"],
    "Plasma_proxy":   ["TNFRSF17", "CD38", "IRF4"],
    "GC_TF":          ["BCL6", "AICDA", "POU2AF1"],
    # Secondary lineage
    "Naive_B":        ["MS4A1", "CD19", "CD22", "PAX5"],
    "Memory_B":       ["CD27", "CD79B", "TNFRSF13B"],
    "BAFF_APRIL":     ["TNFSF13", "TNFSF13B", "TNFRSF13B", "TNFRSF13C", "TNFRSF17"],
    # Tertiary
    "HLA_II":         ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"],
    "CD8_exh":        ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4"],
    "HLA_I_axis":     ["NLRC5", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "T_cell":         ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"],
}


# ================= helpers =================
def log2z(matrix: pd.DataFrame) -> pd.DataFrame:
    """log2(x+1) then z-score per gene across the 12 samples."""
    x = np.log2(matrix.values + 1.0)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    z = (x - mu) / sd
    return pd.DataFrame(z, index=matrix.index, columns=matrix.columns)


def composite(zmat: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    present = [g for g in genes if g in zmat.index]
    missing = [g for g in genes if g not in zmat.index]
    if len(present) < 2:
        raise ValueError(f"composite needs ≥2 genes, got {len(present)} present; missing {missing}")
    return zmat.loc[present].mean(axis=0), present, missing


def mw_exact_one_sided(good: np.ndarray, bad: np.ndarray) -> dict:
    """One-sided MW (good > bad) exact. Also report two-sided for context."""
    u1s = stats.mannwhitneyu(good, bad, alternative="greater", method="exact")
    u2s = stats.mannwhitneyu(good, bad, alternative="two-sided", method="exact")
    return {
        "MW_U_1s": float(u1s.statistic),
        "MW_P_1s_good_gt_bad": float(u1s.pvalue),
        "MW_P_2s": float(u2s.pvalue),
    }


def wilcoxon_signed(diffs: np.ndarray) -> dict:
    if np.all(diffs == 0):
        return {"Wilcoxon_stat": np.nan, "Wilcoxon_P_2s": 1.0, "Wilcoxon_P_1s_pos": 1.0}
    try:
        w2 = stats.wilcoxon(diffs, alternative="two-sided", mode="exact")
        wg = stats.wilcoxon(diffs, alternative="greater", mode="exact")
    except Exception:
        return {"Wilcoxon_stat": np.nan, "Wilcoxon_P_2s": np.nan, "Wilcoxon_P_1s_pos": np.nan}
    return {
        "Wilcoxon_stat": float(w2.statistic),
        "Wilcoxon_P_2s": float(w2.pvalue),
        "Wilcoxon_P_1s_pos": float(wg.pvalue),
    }


# ================= main =================
def main():
    TABLES.mkdir(parents=True, exist_ok=True)

    # ---- Load NanoString ----
    raw = pd.read_excel(INPUT)
    raw = raw.set_index("Probe Name")
    # Sort sample columns into our canonical order
    sample_cols = [c for pair in [(p, q) for (_, p, q, *_r) in COHORT] for c in pair]
    missing_cols = [c for c in sample_cols if c not in raw.columns]
    assert not missing_cols, f"missing sample columns: {missing_cols}"
    mat = raw[sample_cols].astype(float)
    mat.index.name = "gene"

    # drop all-zero genes (shouldn't happen with nSolver output, but safety)
    nz = mat.sum(axis=1) > 0
    mat = mat.loc[nz]

    # ---- log2 + per-gene z ----
    zmat = log2z(mat)
    zmat.to_csv(TABLES / "logz_matrix.tsv", sep="\t")

    # ---- Metadata ----
    meta_rows = []
    for subj, pre, post, bin_, score, rna in COHORT:
        meta_rows.append(dict(sample=pre,  subject=subj, timepoint="pre",  response_bin=bin_, score=score, rnaseq_paired=rna))
        meta_rows.append(dict(sample=post, subject=subj, timepoint="post", response_bin=bin_, score=score, rnaseq_paired=rna))
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(TABLES / "meta.tsv", sep="\t", index=False)

    # ---- Per-subject delta (post_z - pre_z) at gene level ----
    delta_rows = {}
    for subj, pre, post, *_ in COHORT:
        delta_rows[f"subj_{subj}"] = zmat[post] - zmat[pre]
    gene_delta = pd.DataFrame(delta_rows)  # genes × subjects
    gene_delta.to_csv(TABLES / "subject_delta.tsv", sep="\t")

    # ---- Composite scores per sample + per-subject delta ----
    comp_rows = []
    composite_meta_rows = []
    comp_delta = {}
    for name, genes in COMPOSITES.items():
        ser, present, missing = composite(zmat, genes)
        composite_meta_rows.append(dict(
            composite=name, n_genes=len(genes),
            n_present=len(present), present=",".join(present),
            missing=",".join(missing),
        ))
        for s in sample_cols:
            comp_rows.append(dict(sample=s, composite=name, score=float(ser[s])))
        # subject delta
        for subj, pre, post, *_ in COHORT:
            comp_delta.setdefault(f"subj_{subj}", {})[name] = float(ser[post] - ser[pre])

    pd.DataFrame(comp_rows).to_csv(TABLES / "composite_scores.tsv", sep="\t", index=False)
    pd.DataFrame(composite_meta_rows).to_csv(TABLES / "composite_definitions.tsv", sep="\t", index=False)
    comp_delta_df = pd.DataFrame(comp_delta).T  # subjects × composites
    comp_delta_df.index.name = "subject_key"
    comp_delta_df.to_csv(TABLES / "composite_subject_delta.tsv", sep="\t")

    # ---- Primary tier P1-P4 ----
    def tier_stats(target_name: str, series_subj: dict[int, float]) -> dict:
        good = np.array([series_subj[s] for (s, *_r, bn, _sc, _rn) in [(c[0], c[3], c[4], c[5]) for c in COHORT] if bn == "good"])  # noqa: E501
        # rebuild cleanly:
        good = np.array([series_subj[s] for (s, _p, _q, bn, _sc, _r) in COHORT if bn == "good"])
        bad  = np.array([series_subj[s] for (s, _p, _q, bn, _sc, _r) in COHORT if bn == "bad"])
        diffs_good = good  # Δ already post-pre; signed-rank vs 0
        diffs_bad  = bad
        out = dict(
            target=target_name,
            good_mean=float(np.mean(good)), good_median=float(np.median(good)),
            bad_mean=float(np.mean(bad)),   bad_median=float(np.median(bad)),
            good_min=float(good.min()), good_max=float(good.max()),
            bad_min=float(bad.min()),   bad_max=float(bad.max()),
            n_good=len(good), n_bad=len(bad),
        )
        out.update(mw_exact_one_sided(good, bad))
        out["Wilcoxon_good_"] = None  # placeholder replaced below
        # Signed-rank within-group
        wg = wilcoxon_signed(diffs_good); out.update({f"Wgood_{k}": v for k, v in wg.items()})
        wb = wilcoxon_signed(diffs_bad);  out.update({f"Wbad_{k}": v for k, v in wb.items()})
        # drop placeholder
        out.pop("Wilcoxon_good_")
        return out

    # Map composite→subject-Δ dict for testing
    def comp_series(name: str) -> dict[int, float]:
        return {subj: comp_delta_df.loc[f"subj_{subj}", name] for (subj, *_r) in COHORT}

    # Gene-level CXCL13
    def gene_series(gene: str) -> dict[int, float]:
        return {subj: float(gene_delta.loc[gene, f"subj_{subj}"]) for (subj, *_r) in COHORT}

    primary_rows = []
    primary_rows.append(tier_stats("P1_CXCL13",          gene_series("CXCL13")))
    primary_rows.append(tier_stats("P2_TLS_8",           comp_series("TLS_8")))
    primary_rows.append(tier_stats("P3_Plasma_proxy",    comp_series("Plasma_proxy")))
    primary_rows.append(tier_stats("P4_GC_TF",           comp_series("GC_TF")))
    pd.DataFrame(primary_rows).to_csv(TABLES / "P1_P4_primary.tsv", sep="\t", index=False)

    # ---- Secondary S1-S3 lineage ----
    sec_rows = []
    sec_rows.append(tier_stats("S1_Naive_B",      comp_series("Naive_B")))
    sec_rows.append(tier_stats("S2_Memory_B",     comp_series("Memory_B")))
    sec_rows.append(tier_stats("S3_BAFF_APRIL",   comp_series("BAFF_APRIL")))
    pd.DataFrame(sec_rows).to_csv(TABLES / "S1_S3_lineage.tsv", sep="\t", index=False)

    # ---- Secondary S4 platform concordance (5 subj with RNA-seq paired) ----
    rna_logtpm = pd.read_csv(RNA_LOGTPM, sep="\t", index_col=0)
    # shared genes
    shared = sorted(set(mat.index) & set(rna_logtpm.index))
    rna_sub = rna_logtpm.loc[shared]
    ns_genes = mat.loc[shared]  # raw nSolver counts (shared)

    # Compute per-subject delta on RNA-seq (log2-TPM already log-scale, so delta = post - pre)
    rna_subjects = [(s, p, q) for (s, p, q, _bn, _sc, rna) in COHORT if rna]
    # p, q are "TNT RNA N" -> rna_logtpm columns use "TNT_RNA_N"
    def col(s: str) -> str:
        return s.replace(" ", "_")
    rna_delta = {}
    for subj, pre, post in rna_subjects:
        rna_delta[f"subj_{subj}"] = rna_sub[col(post)] - rna_sub[col(pre)]
    rna_delta = pd.DataFrame(rna_delta)  # genes × 5 subjects

    # NanoString delta on shared genes (log2+1 of raw counts to match scale)
    ns_log = np.log2(ns_genes.values + 1.0)
    ns_log = pd.DataFrame(ns_log, index=ns_genes.index, columns=ns_genes.columns)
    ns_delta = {}
    for subj, pre, post, _bn, _sc, rna in COHORT:
        if not rna:
            continue
        ns_delta[f"subj_{subj}"] = ns_log[post] - ns_log[pre]
    ns_delta = pd.DataFrame(ns_delta)

    # Per-gene Pearson across 5 subjects
    concord_rows = []
    for g in shared:
        ns_d = ns_delta.loc[g].values
        rn_d = rna_delta.loc[g].values
        # both should be float vectors length 5
        if len(ns_d) != 5 or len(rn_d) != 5:
            continue
        if np.std(ns_d) == 0 or np.std(rn_d) == 0:
            continue
        r, p = stats.pearsonr(ns_d, rn_d)
        concord_rows.append(dict(gene=g, pearson_r=float(r), pearson_P=float(p),
                                  ns_delta_mean=float(ns_d.mean()), rna_delta_mean=float(rn_d.mean())))
    concord_df = pd.DataFrame(concord_rows).sort_values("pearson_r", ascending=False)
    concord_df.to_csv(TABLES / "S4_platform_concordance.tsv", sep="\t", index=False)

    # Summary: fraction positive r, median r
    pos_frac = float((concord_df["pearson_r"] > 0).mean())
    median_r = float(concord_df["pearson_r"].median())
    with open(TABLES / "S4_platform_summary.json", "w") as fh:
        json.dump(dict(n_genes=len(concord_df), fraction_positive_r=pos_frac,
                       median_r=median_r), fh, indent=2)

    # ---- S5 full 730 gene scan ----
    scan_rows = []
    for g in gene_delta.index:
        good = np.array([gene_delta.loc[g, f"subj_{s}"] for (s, _p, _q, bn, *_r) in COHORT if bn == "good"])
        bad  = np.array([gene_delta.loc[g, f"subj_{s}"] for (s, _p, _q, bn, *_r) in COHORT if bn == "bad"])
        if np.all(good == bad[0]) and np.all(bad == bad[0]) and good[0] == bad[0]:
            continue
        try:
            u1s = stats.mannwhitneyu(good, bad, alternative="greater", method="exact")
            u2s = stats.mannwhitneyu(good, bad, alternative="two-sided", method="exact")
        except ValueError:
            continue
        scan_rows.append(dict(
            gene=g,
            good_median=float(np.median(good)), bad_median=float(np.median(bad)),
            good_mean=float(good.mean()),      bad_mean=float(bad.mean()),
            MW_P_1s=float(u1s.pvalue), MW_P_2s=float(u2s.pvalue),
            direction=("good>bad" if good.mean() > bad.mean() else "bad>good"),
        ))
    scan_df = pd.DataFrame(scan_rows)
    # BH-FDR on one-sided P
    m = len(scan_df)
    if m:
        order = scan_df["MW_P_1s"].rank(method="first").astype(int).values
        sorted_p = scan_df["MW_P_1s"].sort_values().values
        bh = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
        bh = np.clip(bh, 0, 1)
        rank_of = {val: i for i, val in enumerate(scan_df["MW_P_1s"].sort_values().index)}
        q_by_idx = {idx: bh[rank_of[idx]] for idx in scan_df.index}
        scan_df["BH_q_1s"] = [q_by_idx[i] for i in scan_df.index]
    scan_df = scan_df.sort_values("MW_P_1s")
    scan_df.to_csv(TABLES / "S5_full_scan.tsv", sep="\t", index=False)

    # ---- Tertiary T1-T3 ----
    ter_rows = []
    def corr_pair(name: str, comp_a: str, comp_b: str):
        a = np.array([comp_delta_df.loc[f"subj_{s}", comp_a] for (s, *_r) in COHORT])
        b = np.array([comp_delta_df.loc[f"subj_{s}", comp_b] for (s, *_r) in COHORT])
        pr, pp = stats.pearsonr(a, b)
        sr, sp = stats.spearmanr(a, b)
        return dict(test=name, comp_A=comp_a, comp_B=comp_b, n=len(a),
                    pearson_r=float(pr), pearson_P=float(pp),
                    spearman_r=float(sr), spearman_P=float(sp),
                    A_values=",".join(f"{x:.3f}" for x in a),
                    B_values=",".join(f"{x:.3f}" for x in b))

    ter_rows.append(corr_pair("T1_HLAII_vs_Plasma",    "HLA_II",     "Plasma_proxy"))
    ter_rows.append(corr_pair("T2_CD8exh_vs_TLS8",     "CD8_exh",    "TLS_8"))
    ter_rows.append(corr_pair("T3_HLAI_vs_Tcell",      "HLA_I_axis", "T_cell"))
    pd.DataFrame(ter_rows).to_csv(TABLES / "T1_T3_cascade.tsv", sep="\t", index=False)

    print("=== DONE ===")
    print(f"Tables written to {TABLES}/")
    for p in sorted(TABLES.glob("*")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
