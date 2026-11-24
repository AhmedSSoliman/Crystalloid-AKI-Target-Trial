---
title: "Manuscript Outline — Combined MIMIC + IDEALIST"
geometry: margin=1in
papersize: letter
fontsize: 11pt
mainfont: "Arial"
linestretch: 2.0
---

# Manuscript Outline

**Working Title:** Setting-Specific Effects of Initial Crystalloid Choice on Acute Kidney Injury: A Multi-Cohort Target-Trial Emulation in Critical Care and Perioperative Surgery

**Authors:** Ahmed Soliman · Daling Shi  
**Target Journal:** Critical Care Medicine / NEJM Evidence / JASN

---

## I. Abstract (~250 words)

- Background: crystalloid controversy, gap in comparative evidence across settings
- Methods: target-trial emulation, AIPW, two independent cohorts
- Results: MIMIC null (−0.13 pp); IDEALIST signal (−0.93 pp)
- Conclusion: setting-specific heterogeneity; perioperative surgery more sensitive

---

## II. Introduction

1. Clinical context: LR vs NS debate, SMART and PLUS trials, residual uncertainty
2. Gap: no direct comparison of initial fluid choice effect across care settings
3. Target-trial emulation as the methodological solution
4. **Contributions** (6 numbered points)
5. Hypothesis: effect differs between critical care and perioperative surgery

---

## III. Methods

### A. Data Sources
- MIMIC-IV: critical-care EHR, BIDMC 2008–2019
- IDEALIST: perioperative EHR, University of Florida Health 2012–2024

### B. Target Trial Specification (parallel for both cohorts)
- Eligibility criteria
- Time zero definition
- Exposure: first standalone LR or NS
- Outcome: 7-day incident AKI (KDIGO Stage ≥1)
- Causal contrast: LR vs NS as assigned

### C. Propensity Score Model
- Covariates (ICU-focused for MIMIC; perioperative for IDEALIST)
- Logistic regression, SAGA solver
- PS trimming [0.05, 0.95]

### D. Estimation
- AIPW doubly-robust estimator
- Bootstrap CIs (200 reps MIMIC; 50 reps IDEALIST)

### E. Sensitivity Analyses (pre-specified, identical structure)
- SA1: No PS trimming
- SA2: Stricter baseline Cr ≤ 2.0
- SA3: Composite death-or-AKI
- SA4: IPTW-only
- SA5: Shorter AKI window (48h / 3-day)
- SA6a: Baseline Cr ≤ 3.0
- PS matching (1:1 NN)
- Fine-Gray competing risks
- E-value unmeasured confounding

### F. Subgroup Analyses
- Sex (Male, Female)
- Care-unit / surgery-type subgroups

### G. Cross-Cohort Comparison
- Unified Table 1 (baseline characteristics side-by-side)
- Cross-cohort forest plot
- Formal heterogeneity discussion

---

## IV. Results

### A. Cohort Flow and Baseline Characteristics
- Table 1 (unified, both cohorts)
- Figure 1A: MIMIC cohort flow
- Figure 1B: IDEALIST cohort flow

### B. Covariate Balance
- Figure 2: Love plots (MIMIC | IDEALIST)
- Figure 3: PS overlap (MIMIC | IDEALIST)

### C. Primary Results
- **Figure 4: Primary comparison panel** (new — MIMIC vs IDEALIST head-to-head)
- MIMIC: −0.13 pp (CI −1.31 to +0.87), null, E-value 1.08
- IDEALIST: −0.93 pp (CI −0.91 to −0.07), significant, E-value 1.32

### D. Sensitivity Analyses
- **Figure 5: Cross-cohort forest plot** (new — shared SAs, both cohorts)
- Table 2: Full sensitivity results (both cohorts)

### E. Subgroup Analyses
- Figure 6: Subgroup forest (MIMIC ICU-type | IDEALIST surgery-type)
- Notable: MICU protective (−4.69 pp); ORTHO NS favored (+0.69 pp)

### F. Competing Risks and Unmeasured Confounding
- Figure 7: CIF plots (MIMIC | IDEALIST)
- Fine-Gray sHR: MIMIC 1.02 (null); IDEALIST 0.956 (significant)
- E-values: MIMIC 1.08 (crosses null); IDEALIST 1.32 (does not)

---

## V. Discussion

1. **Main finding:** setting-specific heterogeneity — null in critical care, signal in surgery
2. **Mechanistic explanation:** surgical AKI window, volume × composition interaction
3. **MICU subgroup:** hypothesis-generating protective LR signal in medical ICU
4. **ORTHO subgroup:** NS favored — tourniquet, blood loss, different physiology
5. **Consistency across estimators:** AIPW, IPTW, matching convergent in IDEALIST
6. **Comparison with RCTs:** SMART/PLUS — policy-level vs first-choice emulation
7. **Limitations:** residual confounding, AKI phenotype (creatinine only), generalizability
8. **Implications:** perioperative LR as default may reduce postoperative AKI

---

## VI. Conclusions

One paragraph: null in critical care, protective signal in perioperative surgery, robust across SAs, E-value supports causality in IDEALIST, clinical implication.

---

## VII. Supporting Sections

- Author Contributions (CRediT)
- Data and Code Availability
- Funding (GMS 7858)
- Declaration of Competing Interests
- References (JAMIA/Vancouver style, ~20–25 citations)

---

## Figure List

| Figure | Description |
|---|---|
| 1A/1B | Cohort flow diagrams (MIMIC / IDEALIST) |
| 2 | Love plots — covariate balance before/after IPTW |
| 3 | PS overlap histograms |
| **4** | **Primary estimates comparison (new cross-cohort panel)** |
| **5** | **Cross-cohort forest plot (shared SAs)** |
| 6 | Subgroup forest plots |
| 7 | CIF competing risks |
| **T1** | **Unified Table 1 (new)** |
| T2 | Full sensitivity results table |
