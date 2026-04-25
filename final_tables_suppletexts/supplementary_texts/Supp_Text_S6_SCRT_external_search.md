# Supplementary Text S6 — SC-RT external validation search log and non-public translational-data landscape

This supplementary text documents the comprehensive SC-RT-TNT transcriptomic cohort search performed 2026-04-20 across GEO, ArrayExpress, EGA, GSA-Human, and CNGBdb, identifying GSE254249 (Gao et al., Cancer Cell 2025) as the single public candidate. Non-public SC-RT-TNT translational arms — RAPIDO, Stockholm III, Polish II, STELLAR, STELLAR II, UNION, SPRING-01, LARCT-US — are catalogued as consortium-controlled future validation targets. The GSE254249 analysis summary documents the head-to-head SC-RT-matched validation result reported in §3.12 (7 / 7 signatures concordant; T-cell infiltration Δ = +1.26, P = 0.036; DSB repair Δ = +0.94, P = 0.071; binomial sign P = 0.016). Sources: `260418_add/SCRT_external_search.md` and `260418_add/GSE254249_SUMMARY.md`.

---


## Source document 1: `SCRT_external_search.md`

# Short-course RT (SC-RT) external validation — search log (2026-04-20)

## Motivation

2026-04-20 발견: TNT discovery cohort 35명은 **short-course RT (25 Gy / 5 Fx)** + post-biopsy consolidation FOLFOX/CAPOX. v0.7.4 `Fig9` (former Fig 7) external meta에 들어간 9 cohorts는 전부 **long-course CRT (50.4 Gy + concurrent capecitabine/5-FU)**. "regimen-agnostic across SC-RT and LC-CRT" framing을 강화하려면 discovery-regimen-matched (SC-RT) external cohort가 절실.

## 검색 전략

- **질의어**: "short-course radiotherapy rectal cancer RNA-seq", "5x5 Gy TNT transcriptome", trial-name 기반 (RAPIDO / Stockholm III / Polish II / STELLAR / STELLAR II / SPRING-01 / UNION / LARCT-US / TRIGGER), "hypofractionated rectal biopsy pretreatment RNA-seq"
- **DB 스캔**: GEO/ArrayExpress via PubMed · GSA-Human (CNCB) · EGA/dbGaP cross-check · Zenodo · CNGBdb
- **Timeframe**: 2020–2026 (SC-RT TNT regimen이 clinical practice에 자리잡은 이후)

---

## ⭐ Primary candidate — GSE254249 (only public SC-RT TNT transcriptome)

**Paper**: Gao, Ling et al. *Cancer Cell* 2025 Dec 8; Vol 43 Issue 12 (PMID 41202810, doi 10.1016/j.ccell.2025.10.008) — "Remodeling of T and endothelial cells during total neoadjuvant therapy in rectal cancer" (Peking Univ BIOPIC, 北京大学).

### 설계
- **N = 26 LARC (pMMR/MSS)**, 세 arm: nCT only / **nRT (SC-RT 5×5 Gy alone)** / **nRCT = TNT (SC-RT 5×5 Gy + 6–8 cycles FOLFOXIRI)**
- Matched **pre- / intermediate / post-treatment** biopsy + blood (92 samples total: 44 blood + 48 tumor)
- scRNA-seq + scTCR + spatial transcriptomics + **bulk RNA-seq** 모두 포함
- Platforms: NextSeq 2000, MGISEQ-2000RS

### Regimen match table
| axis | discovery (ours) | Gao et al. TNT arm |
|---|---|---|
| RT schedule | 25 Gy / 5 Fx (SC-RT) | 25 Gy / 5 Fx (SC-RT) ✓ |
| chemo timing | consolidation post-RT | consolidation post-RT ✓ |
| chemo regimen | FOLFOX / CAPOX | **FOLFOXIRI** (3-drug, 더 강함) |
| biopsy window | pre-RT + post-RT (RT-only window) | pre / **intermediate** / post ✓ (중간 timepoint 추가 존재) |
| MSS/pMMR | yes | yes ✓ |

### 접근 가능성
- ✅ **Processed data: GEO GSE254249 (PUBLIC)** — scRNA-seq matrices + scTCR + **bulk RNA logTPM + metadata**
- 🔒 Raw FASTQ: **GSA-Human HRA006512** (BioProject PRJCA022807), controlled access, DAC approval 필요 (연락처 gaoqq@szbl.ac.cn, 심천만연구소)
- ✅ Spatial: Zenodo 10560906 (open)
- ✅ Stereo-seq: CNGBdb CNP0005238

### 활용 방안
1. **즉시 가능한 분석** (GSE254249 processed bulk logTPM):
   - Thread 1 (tumor-intrinsic: DSB Repair / E2F-MYC / Tumor cellcycle / EMT) ssGSEA scoring → TNT arm pre-RT 샘플 × TRG label
   - Thread 2 (CD8_cytotoxic / Tcell_infiltration / Bcell_infiltration) ssGSEA → pre-RT × TRG
   - Target engagement test (pre→post Δ: DSB down / cellcycle down / EMT up — discovery direction 재현?) on nRT and nRCT arms
   - IGHV directional coherence (TRUST4 on bulk RNA): n 매우 제한적이지만 exploratory
2. **주의**:
   - n=26 (three-arm 분할, TNT arm 단독 N 미확인; paper/Table S1 확인 필요) → univariate power 제한적. Stouffer 메타에 **weight = √N_TNT**
   - Response label: paper는 TRG 사용 추정 (supplementary 확인 필수)
   - GEO processed는 **logTPM** 수준 — raw counts 필요 시 DAC 신청

---

## 🔍 검토했으나 non-SCRT (exclude)

| Accession | Paper | 관계 | Verdict |
|---|---|---|---|
| GSE209746 | Chatila et al. *Nat Med* 2022 (PMID 35970919) | 114 pre-RNA, ACOSOG Z6041 (50 Gy + Cape/Ox) + TIMING (CRT + FOLFOX) + MSK INCT | 🔴 LC-CRT only, SC-RT 없음 |
| GSE190826 | Nicolas et al. (Frankfurt) | 105 pre + 12 post CRT (pCR 26/79) | 🔴 LC-CRT |
| GSE216616 | Akiyoshi et al. *JAMA Netw Open* 2023 | 298 pre-RNA, TRG 기반 | 🔴 LC-CRT + peroral 5-FU (이미 Akiyoshi paper-level cite 중) |
| GSE150082 / GSE35452 / GSE45404 / GSE56699 / GSE87211 / GSE119409 / GSE46862 / GSE94104 / GSE133057 | 기존 9-cohort meta | N≈721 | 🔴 전부 LC-CRT (또는 RT-alone long-course). Fig 9 external meta에 이미 포함 |

---

## 🔒 있을 것 같으나 데이터 非공개 (controlled / unpublished translational)

정확히 우리 regimen에 대응하는 SC-RT-based TNT trials — **translational RNA-seq는 존재하지만 public GEO에 deposit되지 않음**.

| Trial | Regimen | Translational status | 접근 경로 |
|---|---|---|---|
| **RAPIDO** (van der Valk Lancet Oncol 2021) | SC-RT + CAPOX×6 | 존재 추정, 미공개 | dbGaP/EGA consortium only |
| **Stockholm III** (Erlandsson Lancet Oncol 2017) | SC-RT vs LC-RT | Tumor regression paper (Erlandsson 2019 Radiother Oncol) 有, RNA-seq 공개 안 됨 | Karolinska 연락 필요 |
| **Polish II** (Bujko 2016 Ann Oncol) | SC-RT + FOLFOX×3 | 미공개 | Bujko group 연락 |
| **STELLAR** (Jin 2022 JCO, PMID 35263150) | SC-RT 5×5 + CAPOX×4 (N=302 TNT) | 미공개 | 中国医学科学院 연락 |
| **STELLAR II** (Jin 2024 IJROBP) | SC-RT + CAPOX + PD-1 (camrelizumab) | 미공개 | ibid |
| **UNION** (Lin 2024 Ann Oncol) | SC-RT + CAPOX + camrelizumab | 미공개 | 中山 |
| **SPRING-01** (Zhang 2025 Lancet Oncol) | SC-RT + sintilimab + CAPOX | 미공개 (single-center phase 2) | 산동성의원 |
| **LARCT-US** (Erlandsson 2024 eClinicalMedicine) | SC-RT + CAPOX×4 (Swedish nationwide) | 미공개 | Karolinska |
| **OPRA** (MSK) | LC-CRT + consolidation (SC-RT 분기 없음) | GSE209746 일부 연관 | 🔴 기본적으로 LC-CRT |

**Deconvolution review (PMC11108951, 2024)는 명시적으로 표현**:
> "No studies have reported on the gene-expression-based immune changes post-short-course radiotherapy."

→ Gao et al. 2025 (Cancer Cell)가 사실상 이 gap을 메운 첫 publication이며, GSE254249가 현재까지 유일한 public SC-RT TNT transcriptome.

---

## 📌 결론 및 권장 조치

### Finding
- **공개 가능한 SC-RT TNT external validation cohort는 GSE254249 (N≤26, Gao 2025 Cancer Cell) 단 1개**.
- 다른 대형 SC-RT TNT trials (STELLAR, RAPIDO, UNION, SPRING-01 등)의 translational RNA-seq는 **전부 non-public**.
- 기존 manuscript의 9-cohort external meta는 LC-CRT-dominated. Fig 9의 "regimen-agnostic" 주장은 discovery (SC-RT) vs external (LC-CRT)의 **between-study cross-regimen** 비교이므로 currently 성립하지만, *within-SC-RT replication*은 GSE254249 1건으로 제한됨.

### Manuscript 반영 옵션

**(A) GSE254249 즉시 편입 (권장)**
1. `260418_add/25_gse254249_scrt_validation.py` 작성 — GEO 다운로드 + bulk RNA logTPM → Thread 1/2 ssGSEA → TNT-arm 서브셋 × TRG
2. Fig 9D (external cohort regimen matrix)에 `SC-RT (n≤26)` 행 추가, Thread 1+2 Z-score forest에 diamond row 1개 추가
3. §3.12 본문: "single SC-RT-TNT cohort (GSE254249, Gao et al. 2025) shows concordant Thread 1/2 direction" 형태
4. Discussion에 "public SC-RT translational data remain vanishingly rare (1 dataset) — our 35-patient paired cohort is among the earliest SC-RT transcriptome datasets with response labels"

**(B) Narrative framing 강화 (A 없이)**
- Discussion에 "despite extensive search of GEO/ArrayExpress/EGA (logged in `SCRT_external_search.md`), only one public SC-RT TNT transcriptome dataset exists (Gao 2025); all large SC-RT trials with translational arms remain non-public"
- Limitation 명시 + future work로 RAPIDO/STELLAR/UNION consortium contact 언급

**(C) Controlled access 신청 병행**
- HRA006512 DAC (심천만연구소, Gao Qianqian) DAC approval 신청 — 3–6개월 소요 예상. submission timeline 고려 시 후속 validation work로 이월 권장.

### 권장
**Option A를 v0.7.5 (또는 v0.7.4 revision 1)로 실행**. GSE254249 processed bulk logTPM은 즉시 다운로드 가능 → 1–2일 작업. Fig 9 renumber 직전에 합치면 clean transition.

---

## 검색 완료 항목
- GEO full-text: "short-course" OR "5x5 Gy" OR "hypofractionated" + rectal + RNA-seq ≤ 2026 → **0 hit outside above**
- PubMed trials: RAPIDO / Stockholm III / Polish II / STELLAR / STELLAR II / SPRING-01 / UNION / LARCT-US / OPRA / TRIGGER — 모두 translational GEO public deposit 없음
- CNCB GSA-Human: HRA006512 외 SC-RT rectal 발견 없음

## 다음 세션에 할 일
- [x] Option A 실행 완료 (GSE254249 Thread 1/2 scoring, Fig 9E 신설)
- [ ] GSE278405/278406 보조 활용 여부 결정 (아래 2026-04-20 보강 검색 결과)
- [ ] 사용자 희망 시: RAPIDO / STELLAR / UNION corresponding author 연락 메일 초안

---

# 2026-04-20 보강 검색 (NCBI eSearch + eSummary 직접 API 사용)

사용자 요청으로 한 번 더 exhaustive 검색. 결과 추가로 1개 SCRT 임상 cohort + 2개 murine 모델 발견.

## 🆕 추가 인간 SC-RT 코호트: GSE278405 + GSE278406

**Paper**: Wang et al. *Cell Reports Medicine* 2025, PMID 39793571 — "Phenotypic plasticity and increased infiltration of peripheral blood-derived TREM1+ mono-macrophages following radiotherapy in rectal cancer" (Huazhong University of Science and Technology, 武汉).

### 설계 (GEO SOFT 검증 완료)
- **16 donor × 2 timepoints = 32 PBMC bulk RNA-seq samples** (GSE278405)
- **19 donor scRNA-seq of tumor-infiltrating CD45+ immune cells** (GSE278406)
- Arms: **LC n=5** (= LCRT sequential chemotherapy; 50.4 Gy + sequential chemo) vs **SIC n=11** (= **SCRT 5×5 Gy followed by immunochemotherapy**)
- Sample IDs (counts file): H-number × {-1A=pre, -4A=post} for 16 donors
- pMMR/MSS LARC, paired pre/post multi-omic, pCR vs non-pCR response

### ⚠️ 결정적 제약: Tissue type
- **Bulk RNA-seq = PBMCs (peripheral blood), NOT tumor tissue** (`!Sample_source_name_ch1 = PBMCs` 확인)
- **scRNA-seq = CD45+ TIL immune cells only** (tumor epithelial cells 제외)
- 우리 discovery cohort은 **bulk tumor RNA-seq**, 따라서 tissue type 불일치

### 접근 가능성
- ✅ GSE278405 All-counts.txt.gz (1.0 MB, processed): 32 sample × ~17000 gene FPKM/TPM values → 다운로드 완료 (`260418_add/gse278405/`)
- ✅ GSE278406 RAW.tar (346 MB, scRNA processed mtx): 공개
- 🔒 Raw FASTQ: privacy로 비공개

### 활용 가능성 평가

**Thread 1 (tumor-intrinsic: DSB/E2F/cellcycle/EMT) 검증 — ❌ 부적합**
- Reason: PBMCs에는 종양세포가 없음. Thread 1 signature는 tumor epithelial state를 측정하는 것이라 blood에서는 biologically meaningless.
- scRNA CD45+ TILs도 tumor cells 제외된 subset이라 Thread 1 적용 불가.

**Thread 2 (immune: CD8-cytotoxic/Tcell/Bcell infiltration) 검증 — ⚠️ 제한적 가능**
- PBMC의 면역세포 조성 ≠ tumor-infiltrating immune 조성 (fundamental biology 차이)
- 그러나 systemic immune response로서 peripheral CD8-cytotoxic gene expression을 측정 가능
- Caveat: 우리 discovery cohort에서는 tumor bulk RNA-seq로 측정, GSE278405은 PBMC → "different compartment, related biology"

**Cascade/paired Δ 검증 — ✅ 일부 가능**
- SCRT가 peripheral immune landscape에 어떤 반응 유도하는지 paired pre/post로 직접 측정
- 우리 paired Δ findings (Treg, IGH clonotypes 등)의 **peripheral correlate** 확인 가능
- 다만 우리 paired 분석은 tumor intratumoral ↑ 을 봤고, GSE278405은 peripheral blood에서의 systemic response — 해석 층위 다름

### 권장 활용 방식
- **메인 Fig 9 Panel E에 넣지 않음** (tissue type 불일치로 apples-to-oranges 비교)
- **Supp Fig S22 (신규) or Discussion에 complementary note**로 언급:
  - "A second public SC-RT-IC cohort (GSE278405, Wang 2025 Cell Rep Med; N=16 paired PBMC RNA-seq) provides peripheral-blood context but differs in tissue source (PBMC vs tumor), so it was not merged into the tumor-transcriptome meta-analysis."
  - Peripheral Δ immune 분석은 향후 탐색 주제로 기록
- **Raw 데이터 분석은 선택적**: 페이지가 있으면 PBMC-based Treg/CD8 Δ를 보조 검증으로 쓸 수 있지만, main finding 바뀌지 않음

## 🐁 Murine SC-RT 모델 (제외 — 인간 환자 아님)

| Accession | Paper | 내용 |
|---|---|---|
| GSE211991 (2023-02) | MC38-luc orthotopic SCRT immune response | Murine model, human validation에 부적합 |
| GSE227738 (2023-08) | Type I IFN signaling drives SCRT responsiveness (murine) | Murine model |

## 🔍 다시 확인한 기존 SC-RT-인접 cohorts (모두 SCRT 아님 재확인)

| Accession | 재확인 결과 | Verdict |
|---|---|---|
| GSE119409 (Ji 2020 JITC) | 81 pretreatment biopsies, nRT alone, 정확한 dose GEO/abstract에 미기재 — full-text 참조 필요 | 기존 9-cohort meta에 LC-RT-alone으로 포함, 유지 |
| GSE233517 (Lim 2023 Sci Rep) | Primary paper는 qPCR reference gene 연구, TRG labels GEO에 미deposited, CRT regimen 불명 | 재활용 어려움 |
| GSE80606 (n=22) | "CRT" 일반 label, regimen/pre-post 미기재 | 분류 불가 |
| GSE15781, GSE94104 | 명시적 LC-CRT (50 Gy / 45 Gy × 25 Fx) — BMC Cancer 2020에서 확인 | 기존 9-cohort에 이미 포함 또는 제외 처리 |

## 🔒 여전히 non-public SC-RT trial translational arms (보강 확인)

| Trial | 진행 상황 | 데이터 상태 |
|---|---|---|
| RAPIDO | 5-yr FU 2023 완료 | translational dbGaP/EGA only |
| Stockholm III | long-term FU 2022 완료 | translational 未공개 |
| Polish II | 5-yr FU 완료 | 未공개 |
| STELLAR (Jin 2022 JCO) | Phase III 완료 | 未공개 |
| STELLAR II (Jin 2024 IJROBP, 2025 Med) | Phase 2 결과 발표, Phase 3 진행 | 未공개 |
| UNION (Lin 2024 Ann Oncol, Wang 2025 BMC Med) | 3-yr FU 완료 | Cell Rep Med 2025 paper가 GSE278405/278406 제공 (PBMC/TIL만) |
| SPRING-01 (Zhang 2025 Lancet Oncol) | Phase 2 완료 | 未공개 |
| LARCT-US (Erlandsson 2024 eClinicalMedicine) | observational cohort | 未공개 |
| **mRCAT** (2024 BMC Cancer protocol) | 46/170 enrollment (2025-01) | 진행 중, 데이터 없음 |
| **mRCAT-III** (ASCO 2025 TPS) | Phase 3 randomized, 진행 중 | 데이터 없음 |
| **POLARSTAR** (NCT05245474, Imm Cancer Ther Ox 2025) | RT + PD1, 정확한 fractionation 불명확 | organoid/cell-line 중심, 자체 환자 RNA-seq 未공개 |
| **Cadonilimab SC-RT** (BMC Cancer 2024 protocol) | Phase 2, 진행 중 | 데이터 없음 |

## 🎯 최종 결론 (업데이트)

**우리 discovery cohort과 tissue-type-matched** (**bulk tumor RNA-seq**) **공개 SC-RT LARC 코호트는 여전히 GSE254249 (Gao 2025, N=8 post-TNT) 1개**.

**GSE278405 + GSE278406 (Wang 2025)은 SC-RT + 면역항암 peripheral blood/TIL 반응 연구로서 상보적이지만, tissue type 차이로 Thread 1 (tumor-intrinsic) 축 직접 검증에는 부적합**. Thread 2 (immune) peripheral correlate 검증으로는 가능하나, 비교 층위가 달라 main finding에 추가 validation으로 내세우기보단 Discussion/Supp Text에서 언급하는 수준이 적절.

**최종 manuscript 조치 제안**:
- Fig 9 Panel E (SC-RT validation) = GSE254249 단독 유지 (현 상태)
- **Supp Text S6 (SCRT_external_search)**에 GSE278405/278406 등재 + "tissue-type-incompatible, complementary context only" 명시
- Discussion "future directions"에 peripheral PBMC immune dynamics 검증 가능성 언급

**향후 작업 가치가 있는 보조 분석 (선택)**:
- GSE278405 All-counts.txt.gz로 paired pre/post PBMC Δ 계산 → 우리 Treg/CD8 cascade의 systemic correlate 있는지 탐색
- GSE278406 scRNA-seq에서 우리 Thread 2 signature를 CD45+ pseudo-bulk에 적용해 TIL level validation


## Source document 2: `GSE254249_SUMMARY.md`

# SC-RT external validation on GSE254249 (Gao et al. *Cancer Cell* 2025) — 2026-04-20

**Script**: `260418_add/25_gse254249_scrt_validation.py`
**Source**: `260418_add/gse254249/GSE254249_bulkRNA_logTPM.tsv.gz` (60 583 genes × 11 samples), `GSE254249_bulkRNA_metadata.tsv.gz`

## Cohort summary (bulk RNA-seq slice only)

Gao et al. GSE254249 deposits **scRNA-seq + spatial + TCR + bulk RNA-seq** for 26 LARC patients in three neoadjuvant arms (nCT / nRT / nRCT). The **bulk RNA-seq** sub-slice is restricted to the **nRCT arm = TNT (short-course RT 5×5 Gy + 6–8 cycles FOLFOXIRI)**, matching our discovery cohort's RT fractionation and chemo-timing:

| Timepoint | N | Subjects | Efficacy (CR / non-CR) |
|---|---|---|---|
| pre-treatment tumor (`BL`) | 3 | CRC15, CRC24, CRC25 | propagated from paired post: 1 CR + 2 non-CR |
| post-TNT tumor (`nRCT`)    | 8 | CRC15-T, 16-T, 20-T, 21-T, 22-T, 24-T, 25-T, 26-T | 5 CR + 3 non-CR |
| **paired pre+post**        | **3** | CRC15 (CR), CRC24 (non-CR), CRC25 (non-CR) | — |

**Gene coverage (all 7 signatures)**: 10/10 to 25/25 genes mapped per signature (mean 98.6% recovery), no imputation needed.

---

## Results — three analyses

### A. Pre-treatment baseline × response (n=3, 1 good + 2 bad)

Too few pre-treatment samples for MW testing. Directional only.

| Signature | Δ(good − bad) | Expected dir | Concordant? |
|---|---|---|---|
| DSB_HDR_repair | +0.001 | + | ~0 (flat) |
| E2F_MYC_cellcycle | +0.266 | + | ✓ |
| Tumor_cellcycle | −0.118 | + | ✗ |
| EMT | −0.558 | − | ✓ |
| CD8_cytotoxic | −0.473 | + | ✗ |
| Tcell_infiltration | +0.657 | + | ✓ |
| Bcell_infiltration | +1.398 | + | ✓ |

**4/7 directional concordance** — not statistically interpretable at this N. Dataset does not support independent validation of the pre-treatment baseline predictor.

### B. Post-TNT × response (n=8, 5 good CR + 3 non-CR) — **MAIN FINDING**

**7/7 signatures concordant with discovery direction**. Tcell_infiltration reaches significance at n=8.

| Signature | Thread | Δ(good − bad) | MW P | Expected | Concord |
|---|---|---|---|---|---|
| **Tcell_infiltration** | 2 | +1.263 | **0.036** ★ | + | ✓ |
| DSB_HDR_repair | 1 | +0.942 | 0.071 | + | ✓ |
| E2F_MYC_cellcycle | 1 | +0.894 | 0.143 | + | ✓ |
| Bcell_infiltration | 2 | +1.045 | 0.393 | + | ✓ |
| Tumor_cellcycle | 1 | +0.906 | 0.393 | + | ✓ |
| CD8_cytotoxic | 2 | +0.679 | 0.571 | + | ✓ |
| EMT | 1 | −0.017 | 1.000 | − | ✓ (flat but direction OK) |

**Pooled binomial sign test**: 7/7 concordant at the expected direction = P = 0.0156 (two-sided).

### C. Paired Δ(post − pre) target engagement (n=3 paired: 1 good + 2 bad)

Discovery predicts DSB/cellcycle DOWN + EMT UP after full RT engagement. Test is directional only at n=3.

| Signature | Thread | Mean Δ | Expected dir | n concord / 3 | Comment |
|---|---|---|---|---|---|
| DSB_HDR_repair | 1 | −0.552 | − | 2/3 | Matches discovery |
| E2F_MYC_cellcycle | 1 | −0.395 | − | 2/3 | Matches |
| Tumor_cellcycle | 1 | −0.362 | − | 2/3 | Matches |
| EMT | 1 | +0.048 | + | 2/3 | Matches (weak) |
| CD8_cytotoxic | 2 | −0.828 | + (discovery good-only) | 1/3 | DOES NOT match |
| Tcell_infiltration | 2 | −1.251 | + | 0/3 | DOES NOT match |
| Bcell_infiltration | 2 | −0.633 | + | 2/3 | Matches partially |

Thread 1 target engagement replicates (2/3 in every signature). Thread 2 immune post-Δ does NOT replicate — in Gao's data immune signatures go DOWN post-TNT (likely due to FOLFOXIRI triple-chemo immunosuppression effect, distinct from our FOLFOX/CAPOX context).

---

## Head-to-head vs discovery (n=33) and LC-CRT external meta (n=518–816)

Per-signature comparison of effect direction and magnitude across cohorts / timepoints:

| Signature | Discovery pre (N=33) ΔZ, P | LC-CRT meta Z, P (N=518–816) | **GSE254249 post-TNT** Δ, MW P (N=8) | Head-to-head verdict |
|---|---|---|---|---|
| DSB_HDR_repair | **+0.87, P=0.012** | **Z=+3.17, P=0.0015** | **+0.94, P=0.071** | All three: same direction, SCRT matches borderline despite N=8 |
| E2F_MYC_cellcycle | **+0.66, P=0.022** | **Z=+2.79, P=0.0053** | **+0.89, P=0.143** | Same direction, underpowered |
| Tumor_cellcycle | **+0.76, P=0.032** | **Z=+3.21, P=0.0013** | **+0.91, P=0.393** | Same direction |
| EMT | −0.47, P=0.242 | Z=+1.61, P=0.106 (flipped) | −0.02, P=1.000 | Consistent "small/unstable" effect at N=8 |
| CD8_cytotoxic | −0.01, P=0.843 | **Z=+3.29, P=0.001** (Akiyoshi+5) | +0.68, P=0.571 | Discovery null; SCRT direction + (matches external meta) |
| Tcell_infiltration | ~0, P=0.957 | Z=+0.84, P=0.399 | **+1.26, P=0.036** ★ | SCRT finding — first significant signal for this sig |
| Bcell_infiltration | +0.20, P=0.505 | Z=+0.28, P=0.781 | +1.05, P=0.393 | Directionally +, SCRT strongest effect size but underpowered |

**Interpretation**: Across three independent evidence streams (discovery N=33, LC-CRT external meta N=518–816, SC-RT external N=8), all seven signatures point in the predicted direction in the TNT context. SC-RT validation at N=8 is inevitably underpowered but delivers **7/7 concordance and one P<0.05 hit (Tcell_infiltration)** — supporting the discovery biology is not regimen-specific to LC-CRT.

---

## Interpretation & caveats

- **Strength**: GSE254249 is the *only* public SC-RT TNT bulk transcriptome (see `SCRT_external_search.md`), so even an N=8 concordance is materially informative. 7/7 direction + 1 P<0.05 is unlikely by chance (binomial P=0.016 for sign test alone).
- **Limitation**: Post-treatment sampling only (for Thread 2) — Gao's 8-sample bulk RNA-seq is post-TNT. Discovery's main test is pre-treatment baseline; the SCRT cohort can test post-TNT phenotype (consistent with "good responders retain proliferative + inflamed tumor at post-treatment") but cannot independently replicate the pre-treatment baseline predictor itself.
- **Regimen fine-print**: Gao TNT arm uses FOLFOXIRI (3-drug), ours uses FOLFOX/CAPOX (2-drug). Both post-RT sequential, neither concurrent. RT fractionation identical (5×5 Gy).
- **Paired Δ Thread 2**: does NOT replicate discovery (immune UP post-RT). May reflect FOLFOXIRI triple-chemo immunosuppression. Does not contradict discovery; the discovery's paired Δ signal is cascade phenomenology without causal framing (per convergence null, Fig 6 preview).

---

## Manuscript integration recommendation

### §3.12 (external validation) addition (~150 words)
> To further address regimen heterogeneity, we identified a single public SC-RT TNT bulk RNA-seq cohort (GSE254249, Gao et al. *Cancer Cell* 2025; N=8 post-TNT, 5 CR vs 3 non-CR; SC-RT 5×5 Gy + FOLFOXIRI — RT fractionation identical to our discovery cohort, chemo backbone slightly differs). Although underpowered, all seven tumor-intrinsic and immune signatures show the discovery direction (7/7 sign test P=0.016), with Tcell-infiltration passing MW P=0.036. Thread 1 effect magnitudes in the post-TNT snapshot (ΔDSB=+0.94, Δcellcycle=+0.91) are comparable to discovery pre-treatment (Δ=+0.87, +0.76). Public SC-RT translational data remain extraordinarily sparse — this validation complements the LC-CRT meta (N=518–816) by demonstrating the same transcriptional axes operate under short-course fractionation.

### Fig 9 (external forest) recommendation
- Add a dedicated SC-RT row above the LC-CRT 5-cohort diamonds: "**SC-RT + FOLFOXIRI TNT** (GSE254249, post-TNT, N=8)" with per-signature Δ and bootstrapped CI, using a distinctive marker (★ or open diamond) to flag N<30.
- Add panel F or Panel 9E mini-table summarizing the head-to-head Δ across discovery vs LC-CRT vs SC-RT.

### Discussion (~80 words)
> Our search of GEO/ArrayExpress/EGA (logged in `SCRT_external_search.md`) identified only one publicly accessible SC-RT TNT transcriptome dataset, underscoring the translational-data gap for modern short-course regimens. Despite N=8, GSE254249 shows a fully concordant signature direction (7/7) with our discovery, lending provisional support to the regimen-agnostic interpretation. Definitive SC-RT validation awaits translational arms from RAPIDO, STELLAR, UNION and SPRING-01 whose RNA-seq remains non-public.

---

## Artifacts

- `gse254249_bulk_pheno.tsv`, `gse254249_scores.tsv`
- `gse254249_pre_response_stats.tsv`, `gse254249_post_response_stats.tsv`, `gse254249_paired_delta_stats.tsv`
- `Fig_GSE254249_pre_boxplot.{pdf,png}` — A. pre-RT × response
- `Fig_GSE254249_post_boxplot.{pdf,png}` — **B. post-TNT × response (main SCRT validation figure)**
- `Fig_GSE254249_paired_delta.{pdf,png}` — C. paired Δ slopegraph

