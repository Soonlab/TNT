# TNT Project — Analysis Progress Log

**Started:** 2026-04-13
**Workspace:** `/mnt/sda1/data/TNT/analysis/`
**Goal:** Validate 기존 4-step analysis + execute A/B/C/D/E/F + produce figures/tables/manuscript.

---

## Cohort summary (final)
- **N=35 subjects**, WES 77 samples (1 missing: 13-N), RNA-seq 56 samples
- **Good=18** (TRG0:15, TRG1:3), **Bad=17** (TRG2:12, TRG3:5)
- **Paired(Y)=14** pre+post; **Unpaired(N)=21** mostly pre±normal
- cT imbalance: T4 enriched in bad (7/17 vs 2/18, p=0.086)
- Sex imbalance: M-biased in good (15/18 vs 10/17, p=0.146)
- Files: `00_cohort/{clinical_master.tsv, wes_inventory.tsv, rna_inventory.tsv, table1.tsv}`

## Key finding: WES VCFs are GERMLINE, not somatic
- Macrogen `final.vcf` = GATK HaplotypeCaller/GenotypeGVCFs (germline)
- **Action needed:** For somatic analyses → either
  1. Re-run Mutect2 T-N (matched) / tumor-only + PoN (unmatched) on `.recal.bam`, OR
  2. T-N variant subtraction from annotated xlsx + gnomAD<0.001 + FILTER=PASS + AD thresholds (fast, acceptable for TMB/signature/driver)
- Recommendation: do BOTH (pseudo-somatic fast pass + Mutect2 gold-standard rerun)

---

## Status snapshot (by task ID)

| Task | Status | Notes |
|---|---|---|
| #1 Cohort + Table 1 | ✅ DONE | good=18/bad=17, T4 imbalance noted |
| #16 RNA immune signatures (22 sigs) | ✅ DONE | CD8_proliferation↑ good (pre, p=0.035) |
| #3 TMB | ⏳ Script ready, need pseudo-somatic pipeline |
| #2 SBS refit | ⏳ Need SigProfilerAssignment install |
| #4 MSI | ⏳ Need MSIsensor2 install + reference |
| #5 HLA typing | ⏳ Need OptiType/POLYSOLVER |
| #6 HLA LOH | ⏳ Blocked by HLA typing |
| #7 Neoantigen | ⏳ Blocked by HLA + somatic calls |
| #8 CNV | ⏳ Need FACETS/CNVkit + SNP pileups |
| #9 HRD | ⏳ Blocked by CNV |
| #10 Driver panel | ⏳ Blocked by somatic calls |
| #11 Clonal evolution | ⏳ Blocked by somatic calls |
| #12 DEG good vs bad | pending, use DESeq2 on count matrix |
| #13 GSEA | pending, chained after DEG |
| #14 ssGSEA matrix | partial — 22 sigs scored; Hallmark 50 still pending |
| #15 Full deconvolution | partial — coarse immune cell sigs done; CIBERSORTx/xCell/TIMER pending |
| #17 TCR (TRUST4) | pending |
| #18 CMS/CRIS | pending (CMScaller, R) |
| #19 Resistance sigs | partial — TGFb/EMT/Hypoxia/Stemness/CAF done |
| #20 Multi-omics integration | pending |
| #21 ML predictor | pending |
| #22 External validation | pending (GSE150082, TCGA-READ) |
| #23 Figures/Tables/Manuscript | initial figs only |

---

## Preliminary biological findings (2026-04-13 session)

### Pre-treatment (n=33 RNA)
| Signature | good mean z | bad mean z | Δ | p | q |
|---|---|---|---|---|---|
| CD8_proliferation | 0.73 | 0.15 | +0.58 | **0.035** | 0.72 |
| MHC_II | −0.39 | +0.15 | −0.54 | 0.075 | 0.72 |
| Stemness_mRNAsi | +0.17 | −0.01 | +0.18 | 0.14 | 0.72 |
| B_cell | −0.04 | −0.25 | +0.21 | 0.17 | 0.72 |

- **CD8_proliferation signal confirmed** (matches user's prior finding)
- MHC_II opposite direction — worth investigating (possibly high in bad because of immune evasion via checkpoint upregulation without effective cytotoxic axis)
- After BH correction across 22 sigs, none pass q<0.05 — n limited

### Post-treatment (n=13 RNA)
| Signature | good mean z | bad mean z | Δ | p |
|---|---|---|---|---|
| MHC_II | +0.88 | −0.11 | +0.98 | 0.051 |
| Treg | +0.92 | −0.04 | +0.96 | 0.14 |
| CD8_exhaustion | +0.67 | −0.00 | +0.67 | 0.23 |
| IFNg_Ayers_18 | +0.50 | −0.04 | +0.54 | 0.29 |
| Checkpoint_inhibitory | +0.71 | +0.07 | +0.64 | 0.37 |

- **Post-treatment inflamed microenvironment in good responders** (direction consistent across all immune axes)
- n=13 limits power; effect sizes large (Δ≈0.5-1.0 z)
- Pre (CD8 prolif) → Post (activated/exhausted inflamed TIL) — dynamic immune activation hypothesis

### Immediate interpretation
User's hypothesis (antigen presentation → CD8 proliferation → response) is **partially supported** by CD8 proliferation pre signal, but **class I antigen presentation** signature itself was not significant at pre. MHC_II shows **opposite direction at pre** but **concordant direction at post**. Suggests model should be refined: *effective response = pre-existing proliferative CD8 + post-treatment immune activation*.

---

## Next session priorities (ordered)

### P0 — complete WES somatic backbone
1. Pseudo-somatic variant extraction from annotated xlsx (fast, 78 samples)
2. TMB per sample + response association
3. SBS signature refit (SigProfilerAssignment) → verify MMR(1/9/12) and SBS3(5/14)
4. MSI calling (MSIsensor2 tumor-only as backup to tumor-normal)
5. Driver mutation oncoprint (CRC panel)

### P1 — deepen RNA
6. DESeq2 DEG good vs bad (pre), covariate-adjusted (sex, cT)
7. fgsea Hallmark + Reactome
8. CMS classification (CMScaller)
9. Full CIBERSORTx / xCell deconvolution
10. TCR TRUST4 run

### P2 — multi-omics integration
11. HLA typing + LOH + neoantigen (when WES somatic ready)
12. CNV/HRD
13. Mediation analysis (neoantigen → antigen_pres → CD8_prolif → response)
14. ML predictor (elastic net, LOOCV)

### P3 — manuscript
15. Fig1 cohort + Fig2 WES + Fig3 immune + Fig4 integration + Fig5 ML
16. External validation (GSE150082, TCGA-READ)
17. Manuscript draft (Intro/Methods/Results/Discussion)

---

## Environment / tools installed so far
- Python: pandas, scipy, statsmodels, matplotlib, seaborn, openpyxl (in `rnaseq_arabidopsis` conda env)

## Tools to install
- `SigProfilerAssignment` (SBS refit)
- `MSIsensor2` (+ GRCh38 microsatellite sites)
- `OptiType` or `POLYSOLVER` (HLA typing)
- `LOHHLA` (HLA LOH)
- `pVACtools` + NetMHCpan 4.1 (neoantigen; NetMHCpan requires registration)
- `FACETS` R package (CNV)
- `scarHRD` R package (HRD score)
- `dNdScv`, `maftools` (R)
- `PyClone-VI` (clonal evolution)
- `CMScaller` (R, CRC CMS)
- `TRUST4` (TCR)
- `CIBERSORTx` (web/docker, LM22)
- `xCell`, `immunedeconv` (R)
- `GATK4` Mutect2 (T-N somatic rerun)

## Figures produced so far
- `figures/Fig_pre_signatures_boxplot.png`
- `figures/Fig_pre_signature_heatmap.png`

## Scripts
- `scripts/01_build_tpm_and_signatures.py` — rebuild anytime

---

## 2026-04-15 — v0.4 Genome Medicine revision

**Tasks completed**
- **Task 1 (DFS/OS)**: audited `meta.xlsx`, `meta_WES.xlsx`, `meta_RNA.xlsx`, `clinical_master.tsv`; recursive search of `/mnt/sda1/data/TNT/` for `*surviv*|*clinic*|*followup*|*DFS*|*OS.*`. **No survival data present.** Documented in `analysis/clinical_survival_status.md`; KM/Cox deferred. Not fabricated.
- **Task 2 (Fig 5 neoantigen cascade — new main figure)**: `scripts/26_fig5_neoantigen_cascade.py` → `figures/panels/Fig5_neoantigen_cascade.{pdf,png}`. 4-panel, Arial-like sans, thin spines, `good=#2E86AB bad=#E63946`, 600 dpi PNG + vector PDF, 7.2×5.5 in canvas. Panel A pre missense (P=0.987), B pre strong binders (P=0.546), C paired Δ binders slopegraph+box (P=0.429), D HLA-LOH subj 3 & 4 paired allelic-imbalance pre→post.
- **Task 3 (external meta as supplementary)**: `scripts/27_supp_external_meta.py` → `figures/supp/SuppFig_external_forest.{pdf,png}`, `SuppFig_meta_zscore.{pdf,png}`, `SuppFig_GSE150082_DSB.{pdf,png}` (same style helper). Narrative in `manuscript/supplementary/Supp_external_meta.md` with cohort table (7 cohorts, N=290), Stouffer's Z recap, TNT-regimen-specific-biology framing.
- **Task 4 (manuscript v0.4)**: `manuscript/TNT_manuscript_v0.4_GenomeMedicine.md`. Restructured into 2 core narratives — (1) pre-treatment DNA-repair/cell-cycle predictor (AUC 0.755), (2) post-treatment mutation→neoantigen→HLA-LOH→immune→B-cell cascade. External validation in §3.10 as "TNT-regimen-specific biology". DFS/OS in §3.11 documented-as-deferred. Fig list updated (Fig 5 neoantigen cascade new main; former Fig 7A–C → supplementary).
- **Task 5 (GitHub sync)**: see commit block below.

**Scripts added**
- `scripts/26_fig5_neoantigen_cascade.py`
- `scripts/27_supp_external_meta.py`

**Files added / updated**
- `manuscript/TNT_manuscript_v0.4_GenomeMedicine.md` (new)
- `manuscript/supplementary/Supp_external_meta.md` (new)
- `clinical_survival_status.md` (new)
- `figures/panels/Fig5_neoantigen_cascade.{pdf,png}` (new main)
- `figures/supp/SuppFig_external_forest.{pdf,png}` (new)
- `figures/supp/SuppFig_meta_zscore.{pdf,png}` (new)
- `figures/supp/SuppFig_GSE150082_DSB.{pdf,png}` (new)
