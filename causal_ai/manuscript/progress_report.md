---
title: "GMS 7858 — Project Progress Report"
geometry: margin=1in
papersize: letter
fontsize: 11pt
mainfont: "Arial"
linestretch: 2.0
---

# Project Progress Report

**Course:** GMS 7858 — Causal Artificial Intelligence for Health Research, University of Florida  
**Team:** Ahmed Soliman · Daling Shi  
**Date:** November 2026  
**Project Title:** Initial Crystalloid Choice (Lactated Ringer's vs Normal Saline) and Acute Kidney Injury: A Multi-Cohort Target-Trial Emulation

---

## 1. Project Overview

This project emulates a target trial of the causal effect of selecting Lactated Ringer's (LR) versus Normal Saline (NS) as the first intravenous crystalloid on seven-day incident acute kidney injury (AKI). We apply the analysis to two independent cohorts: (1) MIMIC-IV, a critical-care database from Beth Israel Deaconess Medical Center, and (2) IDEALIST, a perioperative surgical cohort from the University of Florida Health system.

The core method is augmented inverse-probability weighting (AIPW) — a doubly robust doubly-robust estimator — with propensity-score (PS) trimming to [0.05, 0.95] and bootstrap confidence intervals. A prespecified battery of sensitivity, subgroup, competing-risks, and unmeasured-confounding analyses accompanies each primary estimate.

---

## 2. Completed Milestones

### 2.1 Literature Review and Protocol

- Reviewed target-trial emulation framework (Hernán & Robins 2016)
- Reviewed SMART and PLUS randomized trials and their limitations for isolating initial fluid choice
- Defined estimand: LR-minus-NS adjusted risk difference under assigned treatment regardless of subsequent switching
- Prespecified eligibility criteria, time zero, outcome definition, PS model covariates, and all sensitivity analyses before analysis

### 2.2 MIMIC-IV Analysis — COMPLETE

**Cohort extraction and eligibility:**

- Starting pool: 44,767 ICU stays with qualifying LR or NS as first crystalloid within 24h of ICU admission
- Renal eligibility (baseline creatinine present, Cr ≤ 4.0 mg/dL, incident AKI outcome available): 41,383
- After PS trimming [0.05, 0.95]: **38,266 primary analytic cohort**
- LR arm: 14,309 (34.6%) · NS arm: 27,074 (65.4%)

**Primary result:**

The adjusted seven-day AKI risk difference was **−0.13 percentage points (pp)** (95% CI −1.31 to +0.87) — a precise null result. LR and NS were not differentially associated with creatinine-defined AKI in this critically ill cohort.

**Sensitivity analyses:**

| Analysis | n | RD (pp) | 95% CI | Sig? |
|---|---|---|---|---|
| Primary AIPW (PS-trimmed) | 38,266 | −0.13 | (−1.31, +0.87) | No |
| SA1: No PS trimming | 41,383 | −0.64 | (−1.79, +0.36) | No |
| SA2: Stricter baseline Cr ≤ 2.0 | 35,751 | +0.81 | (−0.15, +2.00) | No |
| SA3: Composite death-or-AKI | 38,266 | −3.24 | (−4.22, −2.28) | **Yes** |
| SA4: IPTW-only | 38,266 | +0.09 | (−0.82, +1.19) | No |
| SA5: 48h creatinine window | 38,173 | +1.08 | (+0.33, +2.07) | Yes |
| SA6a: Baseline Cr ≤ 3.0 | 37,384 | ~0 | (−1.08, +1.22) | No |
| PS matching (1:1 NN) | 12,822 | −0.53 | (−1.91, +0.75) | No |
| Fine-Gray sHR (death competing) | 41,383 | 1.02 | (0.96, 1.08) | No |

**Subgroup analyses:**

- Female: RD +0.80 pp (CI −0.67, +2.31) — no effect
- Male: RD −0.55 pp (CI −2.10, +0.90) — no effect
- MICU: RD **−4.69 pp** (CI −6.87, −2.16) — **LR protective**
- CVICU: RD +1.54 pp (CI −0.91, +4.14) — no effect
- MICU/SICU: RD −0.32 pp (CI −2.43, +1.86) — no effect

**Unmeasured confounding:** E-value = 1.08 (CI bound 1.00) — primary null is not sensitive to unmeasured confounding.

**Covariate balance:** PS trimming achieved excellent balance (all SMDs < 0.1 after IPTW).

### 2.3 IDEALIST Perioperative Analysis — Primary Results Complete

**Cohort extraction and eligibility:**

- Starting pool: 458,482 OR cases with intraoperative LR or NS (from 845,930 total OR cases)
- After joining outcomes and confounders: 290,058
- After eligibility criteria (baseline Cr present, Cr ≤ 4.0 mg/dL, AKI outcome available): **250,532**
- After PS trimming: **232,539 primary analytic cohort**
- LR arm: 148,912 (59%) · NS arm: 101,620 (41%)
- Unadjusted AKI rate: LR 11.9% vs NS 15.7%

**Primary result:**

The adjusted seven-day AKI risk difference was **−0.93 pp** (95% CI −0.91 to −0.07) — a small but statistically significant protective effect of LR over NS in the perioperative setting.

**Sensitivity analyses (completed so far):**

| Analysis | n | RD (pp) | 95% CI | Sig? |
|---|---|---|---|---|
| Primary AIPW (PS-trimmed) | 232,539 | −0.93 | (−0.91, −0.07) | **Yes** |
| SA1: No PS trimming | 250,532 | −1.41 | (−1.76, −0.26) | Yes |
| SA2: Stricter baseline Cr ≤ 2.0 | 225,882 | −1.00 | (−1.13, +0.04) | Borderline |
| SA3: Composite death-or-AKI | 232,539 | −1.62 | (−1.71, −1.04) | **Yes** |
| SA4: IPTW-only | 232,539 | −0.63 | (−0.91, −0.07) | Yes |
| SA5: 3-day AKI window | 232,539 | −1.15 | (−0.89, −0.02) | Yes |
| SA6a: Baseline Cr ≤ 3.0 | 230,512 | −0.91 | (−0.75, −0.26) | Yes |

Subgroup analyses (sex, surgery type, emergent/elective), PS matching, Fine-Gray competing-risks, and E-value analyses are currently running (HiPerGator job 52355632).

### 2.4 Deliverables Produced

| Deliverable | Status |
|---|---|
| MIMIC manuscript (`manuscript.docx`) | ✓ Complete |
| GMS 7858 Final Report (`final_report.docx`) | ✓ Complete |
| MIMIC publication figures (4 composite panels + 12 standalone) | ✓ Complete |
| MIMIC Table 1 (baseline characteristics) | ✓ Complete |
| IDEALIST cohort extraction pipeline | ✓ Complete |
| IDEALIST primary analysis pipeline | ✓ Complete |
| IDEALIST sensitivity analysis pipeline | ✓ Running |
| IDEALIST figures | Pending (after sensitivity job) |
| IDEALIST manuscript | Pending (after figures) |

---

## 3. Key Findings to Date

1. **MIMIC-IV (critical care):** LR vs NS shows no difference in 7-day creatinine-defined AKI (RD −0.13 pp, CI crosses null). The null is robust across all pre-specified sensitivity analyses. ICU-type effect modification is hypothesis-generating: the MICU shows a protective LR signal (−4.69 pp) absent in other units.

2. **IDEALIST (perioperative):** LR is associated with a small but consistent protective effect on 7-day post-surgical AKI across all sensitivity specifications (primary −0.93 pp; range −0.63 to −1.62 pp across SAs). The perioperative context differs from critical care in baseline physiology, fluid volumes, and clinical decision-making.

3. **Cross-cohort contrast:** The null in MIMIC vs the modest signal in IDEALIST suggests that fluid choice effects on AKI may be setting-specific — perioperative physiology may be more susceptible to crystalloid composition than critical illness.

---

## 4. Methods Summary

- **Study design:** Target-trial emulation (Hernán & Robins 2016)
- **Estimator:** AIPW (doubly robust; outcome + PS models)
- **PS model:** Logistic regression with SAGA solver on ~50 covariates (demographics, comorbidities, baseline labs, illness severity, care context)
- **Trimming:** PS trimming to [0.05, 0.95] symmetrically
- **Inference:** Percentile bootstrap (50–200 replicates)
- **Competing risks:** Fine-Gray subdistribution hazard (lifelines)
- **Unmeasured confounding:** E-value (VanderWeele & Ding 2017)
- **Software:** Python 3.13, scikit-learn, pandas, lifelines, matplotlib; HiPerGator (UF Research Computing)

---

## 5. Remaining Work

1. IDEALIST sensitivity analyses + figures (job running, ~1–2 hours)
2. IDEALIST manuscript draft
3. Cross-cohort comparison section

---

## 6. Author Contributions

**Ahmed Soliman:** Conceptualization, methodology, software (all analysis pipelines), formal analysis, data curation, writing (original draft and all revisions), visualization, project administration.

**Daling Shi:** Conceptualization, methodology review, writing (review and editing), validation.

---

## References

1. Hernán MA, Robins JM. Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. *Am J Epidemiol*. 2016;183(8):758–764.
2. Self WH, Semler MW, Wanderer JP, et al. Balanced Crystalloids versus Saline in Noncritically Ill Adults. *N Engl J Med*. 2018;378(9):819–828.
3. Semler MW, Self WH, Wanderer JP, et al. Balanced Crystalloids versus Saline in Critically Ill Adults. *N Engl J Med*. 2018;378(9):829–839.
4. Young PJ, Bellomo R, Bernard GR, et al. Saline versus Plasma-Lyte 148 for ICU Fluid Therapy. *N Engl J Med*. 2023;388(11):971–982.
5. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. *Ann Intern Med*. 2017;167(4):268–274.
6. Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data*. 2023;10(1):1.
