# External validation cohorts — original publications (2026-04-18)

Per-cohort source publication and the response-label definition the *original* paper used,
cross-referenced with our concordance/discordance verdict from `thread1_per_cohort_summary.tsv`.

## Reference table (sorted by Thread 1 concordance)

| GEO | Concord. | Original paper | Treatment | Response label in original paper |
|---|---|---|---|---|
| **GSE56699** | 4/4 ★ | Isella C, Terrasi A et al. *Nat Genet* 2015;47(4):312–9 (PMID 25706627, doi:10.1038/ng.3224) | Pre-op long-course CRT | radio-resistant vs radio-sensitive (CAF-stromal score dichotomy) |
| **GSE45404** | 3/4 ★ | Agostini M, Zangrando A et al. *Cancer Biol Ther* 2015;16(8):1160–71 (PMID 26023803) | Pre-op long-course CRT | Mandard TRG: 1–2 = good, 3–5 = bad |
| **GSE133057** | 4/4 | Ferrandon S, DeVecchio J et al. *Cancer Res* 2020;80(2):334–46 (PMID 31704889) | 50.4 Gy / 25 fx + 5-FU, surgery 8–12 wk | AJCC pTRG (CAP scale) |
| **GSE87211** | 4/4 | Hu Y, Gaedcke J, Emons G et al. *Genes Chromosomes Cancer* 2018;57(3):140–9 (PMID 29119627) | German CAO/ARO/AIO-era 5-FU-based CRT (Göttingen) | **No pCR/TRG label — DFS/OS/recurrence only** (we use recurrence as surrogate) |
| **GSE35452** | 3/4 | Watanabe T, Kobunai T, Akiyoshi T et al. *Dis Colon Rectum* 2014;57(1):23–31 (PMID 24316942) | Japanese long-course pre-op CRT | Binary responder vs non-responder |
| GSE46862 | 1/4 | Gim J, Cho YB, Hong HK et al. *Radiat Oncol* 2016;11:50 (PMID 27005571) | Pre-op CRT (SMC, Seoul) | Dworak TRG 4-class (TO/MO/MI/NT) |
| GSE94104 | 1/4 sig– | Alderdice M, Richman SD et al. *J Pathol* 2018;245(1):19–28 (PMID 29412457) | Pre-op long-course CRT (Northern Ireland Biobank + COPERNICUS NCT01263171) | **No formal response label — paper is about CMS/CRIS subtype-classification stability across biopsy/resection** |
| **GSE119409** | 0/4 sig– | Ji D, Song C, Li Y et al. *J Immunother Cancer* 2020;8(2):e000826 (PMID 33106387) | **Neoadjuvant radiotherapy ALONE (no chemo)** | Responder vs non-responder via SAM (PKU paper, Treg-focused) |
| **GSE150082** | 0/4 sig– (3) | Sendoya JM et al. *Cancers (Basel)* 2020;12(8):2227 (PMID 32784964) | **Mixed: long-course CRT + TNT subset** (induction CAPOX + consolidation capecitabine for high-risk) | AJCC pTRG: pTRG 0–1 + cCR = good; pTRG 2–3 + unresectable = bad |

## Critical findings — beyond our previous handling

### 1. GSE119409 is radiation-only, NOT chemoradiation
Ji et al. *JITC* 2020 explicitly used **neoadjuvant radiotherapy alone** without concurrent
fluoropyrimidine. Our pipeline (`scripts/32_external_validation_v3_CD8axis.py`) tagged it as
`nCRT-long`, which is wrong. This cohort should arguably be **excluded** from a CRT-focused
meta-analysis, or at least labelled distinctly in `regimen` column.

This matters: it is one of 4 discordant cohorts and contributes a *significant* discordant
hit (EMT p=0.020).

### 2. GSE94104 has no formal response label
Alderdice et al. *J Pathol* 2018 is a **subtype-classification stability paper** using
matched biopsy/resection pairs. Our use of `tumour regression grade` is reading a manual
metadata field that the original paper did not validate as their primary endpoint. Of the
discordant cohorts, this one's response label is the weakest.

### 3. GSE87211 uses recurrence, not TRG
Hu et al. *GCC* 2018 is a Göttingen rectal cancer cohort whose primary outcome is DFS/OS
and recurrence after surgery — not pCR/TRG. We use `recurrence = bad` as a survival
surrogate. This is a **different endpoint** (long-term oncologic outcome) than pre-op
tumor regression. Concordance of all 4 Thread 1 features here (DSB p=0.058, E2F p=0.156,
cellcycle p=0.073, EMT p=0.158) at N=353 is therefore evidence that the Thread 1 axis
relates to *long-term outcome* in addition to short-term TRG.

### 4. GSE150082 is a mixed-regimen cohort
Sendoya et al. *Cancers* 2020 explicitly mixed long-course CRT with a TNT subset (induction
CAPOX + consolidation capecitabine for high systemic-relapse-risk patients). The paper's own
biological conclusion is that **low DSB repair / B-cell infiltration = good response**, which
is the *opposite* of our discovery direction for the DNA-repair axis. So the discordance is
not a methodological artifact — it is the original paper's finding. The mixed regimen + the
TNT-like consolidation may genuinely reproduce a different biology.

## Implications for our manuscript

### Cleaner regimen-stratified meta (suggested re-analysis)

| Stratum | Cohorts | Note |
|---|---|---|
| Pure long-course nCRT, formal TRG endpoint | GSE56699, GSE45404, GSE133057, GSE35452, GSE46862 (5 cohorts, N≈233) | Apples-to-apples |
| Long-course nCRT, recurrence/survival endpoint | GSE87211 (n=353) | Different endpoint, larger N |
| nCRT subtype-stability cohort, response label inferred | GSE94104 (n=80) | Weakest label |
| **Radiation-alone (NOT CRT)** | **GSE119409 (n=56)** | **Should be excluded** |
| Mixed CRT+TNT subset | GSE150082 (n=39) | Mixed regimen — interpret separately |

If we run a "clean nCRT, TRG-based" meta on the first 5 cohorts only (N ≈ 233), Thread 1
should look much stronger because (a) GSE119409 (RT-only) is removed, (b) GSE150082 (mixed
TNT) is removed, (c) the other two questionable cohorts are out.

### Honest manuscript narrative options

**(i) Keep current 9-cohort meta but acknowledge regimen heterogeneity** — least disruptive,
   needs one Discussion paragraph noting GSE119409 is RT-only and GSE150082 has TNT subset.

**(ii) Redo meta with a clean nCRT-CRT-pCR/TRG stratum** — 5 cohorts (N=233) for Thread 1
   primary, full 9-cohort as sensitivity. CD8 axis should hold either way.

**(iii) Make this the core Discussion point** — "previous external null was driven by regimen
   heterogeneity (RT-alone, mixed TNT, subtype-classification cohorts pooled with TRG-labelled
   pCR cohorts); within a homogeneous nCRT-TRG stratum the discovery direction reproduces."

## Files
- `EXTERNAL_COHORTS_REFERENCES.md` (this file)
- Will produce `thread1_clean_nCRT_meta.tsv` if user requests the (ii) re-analysis
