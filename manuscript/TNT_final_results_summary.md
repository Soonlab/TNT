# TNT Project — Final Results Summary
**Date**: 2026-04-14  
**Workspace**: `/mnt/sda1/data/TNT/analysis/`

## Executive Summary

In 35 MSS locally advanced rectal cancer (LARC) patients receiving total neoadjuvant therapy (TNT), tumor-intrinsic DNA repair and proliferative programs are associated with good response in our discovery cohort. However, **external validation in 7 public microarray cohorts (n=290) did not broadly reproduce this association**, suggesting cohort-specific biology or treatment-regimen-specific effects (TNT vs conventional nCRT) that require larger prospective validation.

---

## 1. Cohort
- N=35, SNU-Hospital, Korean
- WES 77 samples (28 matched normals, 8 unmatched tumors), RNA-seq 56 samples
- TRG-based: good=18 (TRG0-1), bad=17 (TRG2-3)
- cT4 enriched in bad (p=0.086)

## 2. WES Somatic Landscape (Mutect2 T-N)
- **MSI: 0/41 MSI-H** (all MSS, max 0.19%)
- **TMB**: median 1.6/Mb; matched pre good 1.85 vs bad 1.40 (p=0.186)
- **CNV/CIN**: good 0.20 vs bad 0.23 (p=0.659, no difference)
- **HRD proxy LST**: good 4 vs bad 5 (p=0.037) — bad has more large transitions
- **Mutational signatures (COSMIC v3.3 refit)**:
  - Dominant: SBS5 (ageing, mean 248/sample), SBS1 (49)
  - MMR (SBS6/15/20/26): sporadic in both good/bad; NOT a response discriminator
  - SBS3 (HRD): **absent** (0 in all samples)
  - **Prior claim of MMR signature in subj 1/9/12 refuted**: subj 1 has SBS15=22%, subj 9 SBS15=14% (also present in several good responders), subj 12 has 0% MMR
  - **Prior claim of SBS3 in subj 5/14 refuted**: SBS3=0 in all samples
- **Drivers** (classic MSS-CRC pattern): APC 30, TP53 20, KRAS 14, FBXW7 7, KMT2D 4
  - FBXW7 trend ↑ in good (4/16 vs 1/12, OR=3.7, p=0.36)
- **HLA class I typing**: Korean-typical (A*24:02, A*33:03, A*02:01, B*51:01, C*01:02); homozygosity not associated with response

## 3. Transcriptome (Discovery: RNA-seq n=56)

### DEG (pre, good vs bad, adjusted for sex + cT)
- No single gene FDR<0.05
- Top: CPLX2, LRRC37A2, ACCS, LCN15, SHC2

### GSEA Hallmark (pre good vs bad)
| Pathway | NES | p | Direction |
|---|---|---|---|
| E2F_TARGETS | **+2.78** | **8×10⁻²⁶** | ↑good |
| G2M_CHECKPOINT | +2.46 | 2×10⁻¹⁵ | ↑good |
| MYC_TARGETS_V1 | +2.36 | 1×10⁻¹³ | ↑good |
| MYC_TARGETS_V2 | +2.23 | 8×10⁻⁷ | ↑good |
| MITOTIC_SPINDLE | +1.72 | 8×10⁻⁵ | ↑good |
| **EMT** | **−2.16** | **6×10⁻¹⁰** | ↓good |
| MYOGENESIS | −1.79 | 2×10⁻⁵ | ↓good |

### GSEA Reactome top
- Cell cycle checkpoints, M_PHASE, **HOMOLOGY_DIRECTED_REPAIR** (p=10⁻⁹), DSB_REPAIR, DNA_REPLICATION up in good
- ECM_ORGANIZATION, (and TGFβ-related) down in good

### ssGSEA (95 pathways, pre)
| Pathway | p | Δ | Direction |
|---|---|---|---|
| DSB Repair | **0.007** | +0.030 | ↑good |
| Myc Targets V2 | 0.018 | +0.029 | ↑good |
| HDR | 0.020 | +0.038 | ↑good |
| DNA Repair | 0.020 | +0.023 | ↑good |
| G2-M Checkpoint | 0.032 | +0.041 | ↑good |
| E2F Targets | 0.035 | +0.046 | ↑good |
| **CD8 proliferation** | **0.035** | **+0.497** | ↑good |
| MHC_II | 0.074 | −0.054 | ↓good (paradox) |

### TCR/BCR (TRUST4, 56 RNA)
- TRB Shannon, top1, Gini: response-invariant
- Post IGH clonotype count trend ↑ in good (3403 vs 1121, p=0.37) — borderline

### CMS classification
- Pre distribution not associated with response (Fisher p=1)
- CMS4 (mesenchymal) minor trend in bad (4 vs 3, nonsig)

## 4. Clonal Evolution (PyClone-VI, 12 paired pre+post)
- Good responders show larger dominant clone shrinkage (dominant_shrink median -0.67 vs bad -0.15, p=0.34 trend)
- 3/6 good responders had ≥2 shrinking clones vs 2/6 bad
- Consistent with effective chemoradiation eliminating clonal tumor populations

## 5. Multi-omics Integration + ML Predictor
- 37-feature per-subject master table (clinical + WES + RNA signatures + ssGSEA)
- Top response-associated features (Mann-Whitney): DSB Repair (p=0.007), Myc V2 (0.018), HDR (0.020), DNA Repair (0.020), G2M (0.032), CD8 prolif (0.035), E2F (0.035), MHC_II (0.074)
- **LOOCV ML AUC**:
  - **LassoLR full features: AUC=0.755** (best)
  - ElasticNetLR top-8: AUC=0.70
  - RF top-8: AUC=0.70
- RF feature importance top: Myc V2, DSB Repair, Hypoxia, HDR, MHC_II, frac_del, DNA Repair, TLS

## 6. External Validation ⚠️ (Key finding)

### Cohorts (manual response-column mapping)
| Cohort | N | n_good | n_bad | Platform |
|---|---|---|---|---|
| GSE150082 | 39 | 13 | 20 | LARC/TNT, PTRG |
| GSE35452 | 46 | 24 | 22 | LARC nCRT |
| GSE45404 | 80 | 35 | 45 | LARC Responder class |
| GSE68204 | 96 | 20 | 51 | LARC TRG |
| GSE69657 | 30 | 13 | 17 | Stage-IV CRT |
| GSE94104 | 40 | 12 | 28 | LARC nCRT biopsy |
| GSE119409 | 56 | 15 | 41 | LARC, sensitivity |
| **Total** | **387** | **132** | **224** | |

### Per-cohort signature direction (good vs bad delta z-score)
| Cohort | DSB | E2F/Myc | CD8 prolif | EMT |
|---|---|---|---|---|
| GSE35452 | **+0.21** | **+0.23** | **+0.26** | +0.15 |
| GSE45404 | **+0.34*** | +0.28 | +0.35 | +0.05 |
| GSE68204 | +0.07 | +0.06 | +0.02 | -0.14 |
| GSE94104 | -0.03 | +0.04 | -0.06 | -0.07 |
| GSE150082 | -0.41* | **-0.53*** | -0.48 | +0.02 |
| GSE119409 | -0.18 | -0.10 | -0.12 | **+0.34*** |
| GSE69657 | -0.32 | -0.37 | -0.37 | **+0.50*** |

*significant at p<0.05; expected sign for our hypothesis: DSB/E2F/CD8 positive, EMT negative*

### Meta-analysis (Stouffer, 7 cohorts, 290 samples)
| Signature | Z | p (one-sided, expected dir) |
|---|---|---|
| DSB_HDR_repair | -0.15 | 0.56 |
| E2F_MYC_cellcycle | +0.19 | 0.43 |
| CD8_proliferation | +0.06 | 0.48 |
| **EMT** | **-1.83** | **0.97** (reverse direction!) |

### Interpretation
- **2 cohorts agree** with discovery direction (GSE35452, GSE45404 — classical nCRT)
- **3 cohorts disagree** (GSE150082, GSE119409, GSE69657)
- **2 cohorts null** (GSE68204, GSE94104)
- EMT signal notably **reverses** in meta (↑good in validation, ↓good in discovery)
- Plausible explanations:
  1. TNT vs conventional nCRT regimen differences (our cohort used TNT; many validation cohorts used only nCRT)
  2. Platform differences (microarray probe coverage of key DSB/HDR genes variable)
  3. Response definition differences (TRG grading vs Responder/Non-Responder vs sensitivity)
  4. Ethnic/genetic background: Korean TNT vs Western nCRT
  5. True biological variability — our discovery finding may be cohort-specific
  6. Small n in many cohorts (n=30-80) → high variance

## 7. Limitations
- Discovery n=35 single-center
- External validation results mixed → limits generalizability claim
- No HLA LOH analysis (LOHHLA not yet run)
- No neoantigen prediction (pVACtools installed but not yet run — requires VEP cache)
- 8 unmatched WES samples may have residual germline (conservatively excluded in matched analyses)

## 8. Conclusion
In the TNT-treated MSS LARC discovery cohort, pre-treatment DNA repair and cell-cycle programs associate with response (LOOCV AUC=0.755). External validation in conventional nCRT cohorts shows **inconsistent reproduction**, suggesting TNT-specific biology or cohort heterogeneity. **A prospective TNT-specific multicenter cohort is needed** before clinical use. The observed reversal of the EMT axis in meta-analysis warrants careful mechanistic investigation.

---

## 9. Output Inventory

**Tables** (`tables/`):
- `integrated_subject_master.tsv` — 35 subjects × 37 features
- `response_feature_stats.tsv` — ranked feature-response associations

**WES** (`01_wes_signatures/`, `02_wes_tmb_msi/`, `04_wes_cnv_clonal/`, `03_hla/`):
- SBS refit activities, variant master, TMB, MSI, CNV segments, HRD proxy, HLA class I, driver oncoprint, PyClone clonal evolution

**RNA** (`05_rna_deg_gsea/`, `06_rna_immune/`, `07_rna_cms/`, `08_rna_pathway/`):
- DEG, GSEA (Hallmark, Reactome), 22 immune sigs, ssGSEA 95 pathways, CMS, TRUST4

**External** (`11_external_validation/`):
- 9 GEO signature score files + meta-analysis

**ML** (`10_ml_predictor/`):
- LOOCV results, RF importance

**Figures** (`figures/`):
- 8 key figures (signature boxplots/heatmaps, TMB, oncoprint, GSEA, correlation, ML AUC)

**Manuscript** (`manuscript/`):
- `TNT_manuscript_outline.md` — full draft v0.1
- `TNT_final_results_summary.md` — this document

**Scripts** (`scripts/`):
- 19 analysis scripts (Python + R + bash)

## 10. Outstanding (for future work)
- LOHHLA (HLA LOH) analysis
- pVACseq neoantigen (needs VEP GRCh38 cache ~30GB)
- Additional external validation cohorts once FTP issues resolved (GSE87211, GSE46862, GSE56699, GSE133057, GSE190826, GSE119174, GSE93375)
- Prospective TNT-specific validation cohort
- GitHub repo creation + sync
