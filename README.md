# Crystalloid-AKI-Target-Trial

**Setting-Specific Effects of Initial Crystalloid Choice (Lactated Ringer's vs Normal Saline) on Acute Kidney Injury: A Multi-Cohort Target-Trial Emulation**

Ahmed Soliman · Daling Shi  
University of Florida — GMS 7858: Causal Artificial Intelligence for Health Research  
November 2026

---

## Overview

This repository contains the full analysis pipeline, results, figures, and manuscripts for a multi-cohort target-trial emulation of the initial intravenous crystalloid choice (LR vs NS) on seven-day incident acute kidney injury (AKI). We apply an identical causal inference framework — augmented inverse-probability weighting (AIPW) with propensity-score trimming — to two independent cohorts:

| Cohort | Setting | n (analytic) | Primary RD | Significant? |
|---|---|---|---|---|
| **MIMIC-IV** | Critical care ICU | 38,266 | −0.13 pp | No (null) |
| **IDEALIST** | Perioperative surgery | 232,539 | −0.93 pp | **Yes** |

**Central finding:** Setting-specific heterogeneity — a precise null in critical care (E-value 1.08, crosses null) and a consistent LR-protective signal in perioperative surgery (E-value 1.32, does not cross null).

---

## Repository Structure

```
Crystalloid-AKI-Target-Trial/
│
├── 01_mimic/                         # MIMIC-IV critical care analysis
│   ├── 00_exploratory/               # Initial EDA notebooks and figures
│   ├── analysis/
│   │   ├── full_analysis.py          # AIPW primary + SA1–6a + bootstrap
│   │   ├── sensitivity_analyses.py   # SA5–SA10, PS matching, Fine-Gray, E-value
│   │   ├── make_figures.py           # Individual publication figures
│   │   ├── make_paper_figures.py     # 4 composite panels (300 DPI)
│   │   ├── results_summary.csv       # 14 analysis results
│   │   ├── evalue_summary.csv        # E-value
│   │   └── figures/                  # All figures (panels/ + standalone/)
│   ├── docs/                         # Proposal and supplementary docs
│   └── manuscript/
│       ├── manuscript.docx           # Full MIMIC paper
│       ├── final_report.docx         # GMS 7858 final report
│       └── progress_report.docx      # Project progress report
│
├── 02_idealist/                      # IDEALIST perioperative analysis
│   ├── analysis/
│   │   ├── extract_cohort.py         # Cohort extraction
│   │   ├── full_analysis.py          # AIPW primary + SA1–6a
│   │   ├── sensitivity_analyses.py   # SA5, SA7–9, PS matching, Fine-Gray, E-value
│   │   ├── make_figures_idealist.py  # 7 publication figures
│   │   ├── results_summary_idealist.csv
│   │   ├── balance_idealist.csv
│   │   ├── evalue_summary_idealist.csv
│   │   └── figures/                  # 7 figures
│   ├── manuscript/
│   │   └── manuscript_idealist.docx  # Full IDEALIST paper
│   └── run_idealist.sh               # HiPerGator sbatch master script
│
└── 03_combined/                      # Combined multi-cohort paper
    ├── analysis/
    │   ├── make_combined_figures.py  # Cross-cohort figures
    │   └── figures/
    │       ├── fig_primary_comparison.png   # Head-to-head primary estimates
    │       ├── fig_cross_cohort_forest.png  # Shared SAs both cohorts
    │       └── fig_unified_table1.png       # Unified Table 1
    └── manuscript/
        ├── outline.docx                     # Manuscript outline
        ├── manuscript_combined.docx         # Full combined paper
        ├── final_report_combined.docx       # Combined final report
        └── progress_report_combined.docx    # Combined progress report
```

---

## Methods

### Framework
Target-trial emulation (Hernán & Robins 2016). The emulated trial assigns patients at time zero to LR or NS as first crystalloid, with subsequent care proceeding as usual.

### Estimator
Augmented inverse-probability weighting (AIPW) — doubly robust. PS trimming to [0.05, 0.95]. Bootstrap confidence intervals (200 reps MIMIC, 50 reps IDEALIST).

### Cohort Definitions

**MIMIC-IV:** Adults in first ICU stay, first qualifying standalone LR or NS within 24h of ICU admission, baseline creatinine ≤4.0 mg/dL, incident AKI outcome available.

**IDEALIST:** Adult surgical patients, first standalone intraoperative LR or NS in operative window, baseline creatinine ≤4.0 mg/dL, 7-day AKI outcome available. Fluid classification by `med_order_display_name` (NOT `iv_class`, which mislabels LR as "saline" in IDEALIST).

### Outcome
Incident AKI within 7 days: KDIGO Stage ≥1 (creatinine rise ≥0.3 mg/dL within 48h or ≥1.5× reference within 7 days).

### Sensitivity Analyses (pre-specified, identical structure)
SA1: No PS trimming · SA2: Stricter Cr ≤ 2.0 · SA3: Composite death-or-AKI · SA4: IPTW-only · SA5: Shorter AKI window · SA6a: Cr ≤ 3.0 · PS matching (1:1 NN) · Fine-Gray competing risks · E-value

---

## Key Results

### MIMIC-IV (Critical Care)

| Analysis | n | RD (pp) | 95% CI | Sig? |
|---|---|---|---|---|
| Primary AIPW | 38,266 | −0.13 | (−1.31, +0.87) | No |
| SA3: Death-or-AKI | 38,266 | −3.24 | (−4.22, −2.28) | Yes |
| MICU subgroup | 7,577 | −4.69 | (−6.87, −2.16) | Yes |
| Fine-Gray sHR | 41,383 | 1.02 | (0.96, 1.08) | No |
| E-value | — | 1.08 | crosses null | — |

### IDEALIST (Perioperative Surgery)

| Analysis | n | RD (pp) | 95% CI | Sig? |
|---|---|---|---|---|
| Primary AIPW | 232,539 | −0.93 | (−0.91, −0.07) | **Yes** |
| SA1: No trim | 250,532 | −1.41 | (−1.76, −0.26) | Yes |
| SA3: Death-or-AKI | 232,539 | −1.62 | (−1.71, −1.04) | Yes |
| Burn Surgery | 16,645 | −2.10 | (−3.09, −0.42) | Yes |
| Fine-Gray sHR | 232,539 | 0.956 | (0.931, 0.982) | Yes |
| E-value | — | 1.32 | does not cross null | — |

---

## Figures

### Cross-Cohort Comparison (03_combined)

**Primary Estimates: Head-to-Head**  
*AIPW risk difference (pp) with 95% bootstrap CI — MIMIC-IV vs IDEALIST*

![Primary Comparison](03_combined/analysis/figures/fig_primary_comparison.png)

---

**Cross-Cohort Sensitivity Analysis Forest Plot**  
*Shared sensitivity analyses run identically in both cohorts*

![Cross-Cohort Forest](03_combined/analysis/figures/fig_cross_cohort_forest.png)

---

**Unified Baseline Characteristics (Table 1)**  
*Demographics and comorbidities side-by-side across both cohorts*

![Unified Table 1](03_combined/analysis/figures/fig_unified_table1.png)

---

### MIMIC-IV Figures (01_mimic)

**Panel 1 — Study Design & Cohort Flow**

![MIMIC Panel 1](01_mimic/analysis/figures/panels/panel_1_design.png)

---

**Panel 2 — Propensity Score Balance**

![MIMIC Panel 2](01_mimic/analysis/figures/panels/panel_2_balance.png)

---

**Panel 3 — Treatment Effect Estimates**

![MIMIC Panel 3](01_mimic/analysis/figures/panels/panel_3_effects.png)

---

**Panel 4 — Competing Events (Fine-Gray)**

![MIMIC Panel 4](01_mimic/analysis/figures/panels/panel_4_competing.png)

---

**Causal DAG**

![Causal DAG](01_mimic/analysis/figures/fig_causal_dag.png)

---

**Propensity Score Overlap**

![PS Overlap MIMIC](01_mimic/analysis/figures/fig_ps_overlap.png)

---

**Love Plot — Covariate Balance**

![Love Plot MIMIC](01_mimic/analysis/figures/fig_love_plot.png)

---

**AIPW Primary Estimates + Sensitivity Analyses**

![AIPW Results MIMIC](01_mimic/analysis/figures/fig_aipw_results.png)

---

**Subgroup Forest Plot**

![Subgroup Forest MIMIC](01_mimic/analysis/figures/fig_subgroup_forest.png)

---

**Cumulative Incidence Functions (Competing Risks)**

![CIF MIMIC](01_mimic/analysis/figures/fig_cif_competing.png)

---

**E-Value Sensitivity**

![E-Value MIMIC](01_mimic/analysis/figures/fig_evalue.png)

---

**PS Matching Balance**

![Matching Balance MIMIC](01_mimic/analysis/figures/fig_matching_balance.png)

---

### IDEALIST Figures (02_idealist)

**Cohort Flow Diagram**

![Cohort Flow IDEALIST](02_idealist/analysis/figures/fig_cohort_flow_idealist.png)

---

**Propensity Score Overlap**

![PS Overlap IDEALIST](02_idealist/analysis/figures/fig_ps_overlap_idealist.png)

---

**Love Plot — Covariate Balance**

![Love Plot IDEALIST](02_idealist/analysis/figures/fig_love_plot_idealist.png)

---

**AIPW Primary Estimates + Sensitivity Analyses**

![AIPW Results IDEALIST](02_idealist/analysis/figures/fig_aipw_results_idealist.png)

---

**Subgroup Forest Plot**

![Subgroup Forest IDEALIST](02_idealist/analysis/figures/fig_subgroup_forest_idealist.png)

---

**Cumulative Incidence Functions (Competing Risks)**

![CIF IDEALIST](02_idealist/analysis/figures/fig_cif_idealist.png)

---

**E-Value Sensitivity**

![E-Value IDEALIST](02_idealist/analysis/figures/fig_evalue_idealist.png)

---

## Reproducing the Analysis

### Dependencies
```bash
conda create -n causal_aki python=3.11
pip install pandas numpy scikit-learn matplotlib lifelines pyarrow joblib
```

### MIMIC-IV
```bash
cd 01_mimic/analysis
python3 full_analysis.py
python3 sensitivity_analyses.py
python3 make_figures.py
python3 make_paper_figures.py
```

### IDEALIST (HiPerGator)
```bash
cd 02_idealist
sbatch run_idealist.sh
```

### Combined figures
```bash
cd 03_combined/analysis
python3 make_combined_figures.py
```

---

## Data Availability

- **MIMIC-IV:** Freely available at [physionet.org](https://physionet.org/content/mimiciv/) under a data use agreement.
- **IDEALIST:** Available to approved investigators through the University of Florida PRISMAP data governance framework. Raw data are not included in this repository.

---

## Citation

If you use this code or methodology, please cite:

> Soliman A, Shi D. Setting-Specific Effects of Initial Crystalloid Choice on Acute Kidney Injury: A Multi-Cohort Target-Trial Emulation. University of Florida GMS 7858, November 2026.

---

## License

MIT License — code only. Data access governed by respective data use agreements.

---

## Contact

Ahmed Soliman · University of Florida  
GMS 7858 — Causal Artificial Intelligence for Health Research
