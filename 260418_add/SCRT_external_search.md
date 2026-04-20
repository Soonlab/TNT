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
- [ ] 사용자 Option A/B/C 선택
- [ ] Option A 선택 시: `260418_add/25_gse254249_scrt_validation.py` 작성 + Fig 9 업데이트
- [ ] Option B 선택 시: Discussion/Limitations 패러그래프 작성 + MEMORY 반영
- [ ] 사용자 희망 시: RAPIDO / STELLAR / UNION corresponding author 연락 메일 초안
