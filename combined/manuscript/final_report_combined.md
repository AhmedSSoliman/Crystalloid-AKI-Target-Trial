---
geometry: margin=1in
papersize: letter
fontsize: 11pt
mainfont: "Arial"
linestretch: 2.0
---

# Combined Final Report

**Course:** GMS 7858 — Causal Artificial Intelligence for Health Research, University of Florida  
**Team:** Ahmed Soliman · Daling Shi  
**Date:** November 2026  
**Project Title:** Setting-Specific Effects of Initial Crystalloid Choice on AKI: A Multi-Cohort Target-Trial Emulation

---

## 1. Project Overview

This project emulates target trials of the initial intravenous crystalloid choice (Lactated Ringer's vs Normal Saline) on seven-day incident acute kidney injury (AKI) in two independent cohorts: MIMIC-IV (critical care, n=38,266 primary analytic) and IDEALIST (perioperative surgery, n=232,539 primary analytic). The core method is augmented inverse-probability weighting (AIPW) with propensity-score trimming and bootstrap confidence intervals, implemented identically in both cohorts to enable cross-cohort comparison.

---

## 2. Completed Work

### 2.1 MIMIC-IV Analysis — COMPLETE

| Item | Status |
|---|---|
| Cohort extraction (MIMIC-IV) | ✓ |
| PS model + AIPW primary estimate | ✓ |
| Sensitivity analyses (SA1–6a, PS matching, Fine-Gray, E-value) | ✓ |
| Subgroup analyses (sex, ICU type) | ✓ |
| Publication figures (12 standalone + 4 composite panels) | ✓ |
| Manuscript (`manuscript.docx`) | ✓ |
| GMS 7858 Final Report (`final_report.docx`) | ✓ |

**Primary result:** RD −0.13 pp (CI −1.31, +0.87) — precise null. E-value 1.08 (crosses null).

### 2.2 IDEALIST Analysis — COMPLETE

| Item | Status |
|---|---|
| Cohort extraction (250,532 eligible) | ✓ |
| PS model + AIPW primary estimate | ✓ |
| Sensitivity analyses (SA1–6a, PS matching, Fine-Gray, E-value) | ✓ |
| Subgroup analyses (sex, surgery type) | ✓ |
| Publication figures (7 figures) | ✓ |
| Manuscript (`manuscript_idealist.docx`) | ✓ |

**Primary result:** RD −0.93 pp (CI −0.91, −0.07) — significant LR protective effect. E-value 1.32 (does not cross null).

### 2.3 Combined Analysis — COMPLETE

| Item | Status |
|---|---|
| Cross-cohort primary comparison figure | ✓ |
| Cross-cohort forest plot (shared SAs) | ✓ |
| Unified Table 1 figure | ✓ |
| Manuscript outline | ✓ |
| Combined manuscript (`manuscript_combined.docx`) | ✓ |
| Combined final report (`final_report_combined.docx`) | ✓ |
| Combined progress report (this document) | ✓ |

---

## 3. Key Results Summary

### MIMIC-IV (Critical Care, n=38,266)

| Analysis | RD (pp) | 95% CI | Sig? |
|---|---|---|---|
| Primary AIPW | −0.13 | (−1.31, +0.87) | No |
| SA3: Death-or-AKI | −3.24 | (−4.22, −2.28) | Yes |
| SA5: 48h AKI | +1.08 | (+0.33, +2.07) | Yes |
| MICU subgroup | −4.69 | (−6.87, −2.16) | Yes |
| Fine-Gray sHR | 1.02 | (0.96, 1.08) | No |
| E-value | 1.08 | crosses null | — |

### IDEALIST (Perioperative Surgery, n=232,539)

| Analysis | RD (pp) | 95% CI | Sig? |
|---|---|---|---|
| Primary AIPW | −0.93 | (−0.91, −0.07) | Yes |
| SA1: No trim | −1.41 | (−1.76, −0.26) | Yes |
| SA3: Death-or-AKI | −1.62 | (−1.71, −1.04) | Yes |
| SA5: 3-day AKI | −1.15 | (−0.89, −0.02) | Yes |
| Burn Surgery subgroup | −2.10 | (−3.09, −0.42) | Yes |
| Fine-Gray sHR | 0.956 | (0.931, 0.982) | Yes |
| E-value | 1.32 | does not cross null | — |

---

## 4. Central Finding

**Setting-specific heterogeneity:** LR = NS in critical care (precise null, E-value crosses null); LR protective in perioperative surgery (consistent signal, E-value 1.32). This contrast is the primary scientific contribution of the combined project.

---

## 5. Methods Summary

- **Framework:** Target-trial emulation (Hernán & Robins 2016)
- **Estimator:** AIPW (doubly robust)
- **PS model:** Logistic regression (SAGA solver), ~30–50 covariates per cohort
- **Trimming:** PS [0.05, 0.95]
- **Bootstrap:** 200 reps (MIMIC), 50 reps (IDEALIST)
- **Competing risks:** Fine-Gray subdistribution hazard (lifelines)
- **Unmeasured confounding:** E-value (VanderWeele & Ding 2017)
- **Platform:** Python 3.13, scikit-learn, pandas, lifelines, matplotlib; HiPerGator (UF Research Computing)

---

## 6. Author Contributions

**Ahmed Soliman:** Conceptualization, methodology, software, formal analysis, data curation, writing (original draft, all revisions), visualization, project administration.  
**Daling Shi:** Conceptualization, methodology review, writing (review and editing), validation.

---

## References

1. Hernán MA, Robins JM. Using Big Data to Emulate a Target Trial. *Am J Epidemiol*. 2016;183(8):758–764.
2. Self WH, Semler MW, et al. Balanced Crystalloids versus Saline in Noncritically Ill Adults. *N Engl J Med*. 2018;378(9):819–828.
3. Young PJ, Bellomo R, et al. Saline versus Plasma-Lyte 148 for ICU Fluid Therapy. *N Engl J Med*. 2023;388(11):971–982.
4. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. *Ann Intern Med*. 2017;167(4):268–274.
5. Johnson AEW, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data*. 2023;10(1):1.
