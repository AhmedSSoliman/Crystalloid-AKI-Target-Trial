**Initial Lactated Ringer's versus Normal Saline and Subsequent Acute Kidney Injury in ICU Patients: A Causal Analysis Using MIMIC-IV**

**Background and motivation**

Intravenous fluid choice is a clinically relevant question in intensive care. Lactated Ringer's (LR) and normal saline (NS) are commonly used crystalloid solutions. Randomized trials, including SMART, have investigated balanced crystalloids versus saline and kidney-related outcomes. However, those comparisons do not exactly match the effect of choosing LR as the initial recorded ICU fluid, followed by subsequent usual care. [SMART trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC5846085/)

This project would use MIMIC-IV to examine that initial treatment decision. MIMIC-IV contains hospital and ICU records, including laboratory measurements and treatment information, which may support a temporally aligned observational analysis. Dataset access has been obtained; the suitability of the available records is being assessed. [MIMIC-IV documentation](https://physionet.org/content/mimiciv/3.1/)

The main methodological challenge is that fluid choice is not random. Patients receiving LR and NS may differ in kidney function, illness severity, and clinical circumstances before treatment. These differences may also influence subsequent acute kidney injury (AKI). The project would therefore combine an explicit study design with methods for adjusting measured confounding.

**Research question**

Among eligible adult ICU patients, what is the effect of choosing LR rather than NS at the first qualifying recorded ICU fluid initiation on subsequent creatinine-defined AKI during hospitalization within seven days?

The intended estimand is the average difference in outcome risk under the two initial choices in the eligible population, with subsequent treatment allowed to follow usual care. This interpretation will depend on adequate treatment overlap and the assumptions required for causal inference.

**Specific aims**

1. **Establish a reproducible cohort** with clearly defined fluid exposure, pretreatment renal status, and aligned eligibility, treatment assignment, and follow-up.
2. **Estimate the risk difference** for subsequent creatinine-defined AKI using inverse probability weighting and doubly robust estimation.
3. **Assess robustness** to reasonable changes in fluid classification, baseline creatinine definitions, and treatment overlap, while examining limitations in outcome monitoring.

**Proposed study design**

The study will use target trial principles: first specify the treatment decision and follow-up that an ideal trial would evaluate, then assess how closely the available observational records can represent that design. This approach helps avoid inconsistencies between treatment assignment and the start of follow-up. [Target trial design principles](https://pmc.ncbi.nlm.nih.gov/articles/PMC5124536/)

| **Component**               | **Preliminary specification**                                                                                              |
| --------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| Population                  | Adults aged ≥18 years at their first recorded ICU stay                                                                     |
| ---                         | ---                                                                                                                        |
| Treatment initiation window | First qualifying LR or NS initiation within 24 hours of ICU admission; narrower windows considered in sensitivity analyses |
| ---                         | ---                                                                                                                        |
| Treatment comparison        | Initial recorded LR versus initial recorded NS                                                                             |
| ---                         | ---                                                                                                                        |
| Time zero                   | Start of the first qualifying administration                                                                               |
| ---                         | ---                                                                                                                        |
| Renal eligibility           | No chronic dialysis and sufficient pretreatment information to assess prevalent AKI; operational criteria to be finalized  |
| ---                         | ---                                                                                                                        |
| Baseline information        | Measurements and clinical information preceding time zero                                                                  |
| ---                         | ---                                                                                                                        |
| Primary outcome             | New creatinine-defined AKI recorded during hospitalization within seven days after time zero                               |
| ---                         | ---                                                                                                                        |
| Primary effect measure      | Adjusted absolute risk difference with a 95% confidence interval                                                           |
| ---                         | ---                                                                                                                        |
| Subsequent treatment        | Later fluid switching will not change initial group assignment                                                             |
| ---                         | ---                                                                                                                        |

"Initial" will refer to the **first qualifying recorded ICU administration**, rather than the patient's first-ever fluid exposure. Earlier emergency department or ward fluids will be assessed where available.

Fluid classification will distinguish standalone fluid administration from medication carriers, diluents, and flushes as far as the records permit. Standalone administration will not automatically be interpreted as confirmed resuscitation intent. Classification rules will also be checked for dependence on information observed after infusion initiation, such as completed infusion volume.

**Outcome and baseline kidney function**

The proposed outcome is creatinine-defined AKI based on the KDIGO 2012 creatinine criteria, with explicit measurement comparisons and time windows. Because urine-output criteria may not be included, the outcome will be described as _creatinine-defined AKI_, rather than the full KDIGO endpoint. [KDIGO 2012 guideline](https://kdigo.org/wp-content/uploads/2016/10/KDIGO-2012-AKI-Guideline-English.pdf)

A recent pretreatment creatinine measurement does not, by itself, establish that a patient is free of AKI. The feasibility assessment will distinguish patients with evidence of existing AKI, patients with sufficient observed history without detected AKI, and patients whose baseline status cannot be determined. Reference creatinine definitions will use only information preceding treatment.

Missing follow-up measurements will not be treated as absence of AKI. Death and discharge before day seven will also be evaluated explicitly. The final analysis plan will specify how these events define or limit the inpatient endpoint, and whether competing-event or observation-process methods are needed.

**Analysis plan**

A clinically informed causal diagram will guide selection of pretreatment confounders. Candidate variables include demographics, historical kidney disease and other comorbidities, recent creatinine and electrolytes, blood pressure, pretreatment organ support, ICU type, and recorded earlier fluid exposure.

The analysis will proceed in three stages:

1. Describe cohort characteristics, missingness, and unadjusted outcome frequencies.
2. Estimate propensity scores—the probability of receiving LR given measured baseline characteristics—and assess treatment overlap and covariate balance.
3. Estimate adjusted risks and their difference using inverse probability weighting and augmented inverse probability weighting, a doubly robust method combining treatment and outcome models.

Machine learning may be used for the treatment and outcome models, with cross-fitting to reduce overfitting. Model complexity will be kept appropriate for the available sample size and course timeline.

Causal interpretation will require adequate overlap, sufficiently measured confounding, and a well-defined treatment comparison. Doubly robust estimation does not remove bias from unmeasured clinical decision factors.

**Feasibility assessment underway**

The local data audit is designed to answer three questions before the protocol is finalized:

| **Question**                                               | **Evidence needed**                                                                                             |
| ---------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| Can LR and NS fluid administration be identified reliably? | Item mappings, order information, ambiguity rates, and counts under defensible exposure definitions             |
| ---                                                        | ---                                                                                                             |
| Can pretreatment kidney function be established?           | Creatinine availability and timing, repeated historical measurements, and ascertainment of existing AKI         |
| ---                                                        | ---                                                                                                             |
| Are the groups sufficiently comparable?                    | Eligible counts by group, baseline differences, propensity-score overlap, and diagnostic effective sample sizes |
| ---                                                        | ---                                                                                                             |

The audit will also assess follow-up creatinine availability, death, and discharge. No minimum sample size or claim of adequate power is assumed at this stage. Precision will be evaluated once eligible counts and outcome information are available.

If comparable patients exist only within certain ICU types or clinical settings, a narrower target population may be considered. Any restriction will be justified clinically and defined before treatment-effect comparisons are examined.

**Expected contribution and scope**

The project would provide a transparent application of causal inference to a clinically meaningful treatment decision, with particular attention to timing, exposure definition, baseline renal status, and outcome ascertainment.

The core deliverables would be a reproducible cohort pipeline, a causal diagram, balance and overlap diagnostics, adjusted effect estimates, and sensitivity analyses. These align with the course topics of confounding, treatment-effect estimation, machine learning, and target trial emulation.

The primary scope will be average treatment effects. Subgroup analyses would be optional and pursued only if the main analysis is feasible and sufficiently precise. A null or imprecise result would still be informative; the project does not depend on demonstrating that one fluid is superior.

**Tentative timeline**

- **Before the October 5 proposal deadline:** Complete the feasibility audit, review related studies, and finalize the population, exposure, outcome, and causal diagram.
- **October:** Build and validate the cohort; complete baseline and overlap diagnostics.
- **November:** Conduct effect estimation and focused sensitivity analyses.
- **Before the final course deadline:** Complete the technical report, reproducible code, and presentation.

**Questions**

1. Is initial LR versus NS a clinically meaningful and sufficiently focused comparison?
2. Should the population be restricted to a particular ICU type or clinical indication?
3. What pretreatment history would be sufficient to assess baseline AKI reliably?
4. Is seven-day inpatient creatinine-defined AKI a suitable primary endpoint?
5. Which clinical factors are most important to capture when explaining why clinicians choose LR or NS?

The immediate decision is whether the clinical question and proposed design are worth pursuing, while the data audit determines whether MIMIC-IV can support them.