# Clinical Survival Data Status — TNT Cohort

**Date**: 2026-04-15
**Author**: Automated audit (Task 1 of GM-target revision)

## Summary

**No disease-free survival (DFS) or overall survival (OS) data are available in the currently curated TNT metadata.** Kaplan-Meier and Cox analyses are therefore not performed in manuscript v0.5 and are flagged as a limitation.

**Survival data are not yet mature — this is a recently-accrued cohort; DFS/OS analyses are planned for a follow-up report.**

## Files audited

| File | Sheet(s) | Columns |
|------|----------|---------|
| `/mnt/sda1/data/TNT/meta.xlsx` | Sheet1 | Macrogen Number, Chart #, Initial, 성별, 나이, Subject No., Clinical T Stage, TNT Response, TNT Response.1, Binomial, prepost_set |
| `/mnt/sda1/data/TNT/TNT_WES/meta_WES.xlsx` | Sheet1 | sample_id, subject_id, timepoint, clinical_T, response_num, response_bin, sex, age, prepost_set |
| `/mnt/sda1/data/TNT/TNT_RNAseq/meta_RNA.xlsx` | Sheet1 | sample_id, subject_id, timepoint, clinical_T, response_num, response_bin, sex, age, prepost_set |
| `/mnt/sda1/data/TNT/analysis/00_cohort/clinical_master.tsv` | – | subject_id, macrogen_chart, initial, sex, age, subject_no, cT, response_num, response_bin, prepost_set |

A recursive search under `/mnt/sda1/data/TNT/` for filenames matching `*surviv*`, `*clinic*`, `*followup*`, `*DFS*`, `*OS.*` returned no additional survival-bearing datasets.

## Implication

- Endpoint used throughout the manuscript remains **TNT response (good vs bad)** operationalised from `TNT Response` / `response_bin` per subject.
- Long-term outcome analyses (DFS, OS, locoregional recurrence) are deferred to a follow-up study when outcome data become available.
- We recommend the clinical team export an updated cohort table including `date_of_diagnosis`, `date_last_followup`, `recurrence_event`, `recurrence_date`, `death_event`, `death_date` so KM/Cox can be added at revision.
