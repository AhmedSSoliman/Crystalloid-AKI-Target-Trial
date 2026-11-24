"""
IDEALIST LR vs NS -> 7-day post-surgical AKI  —  Cohort extraction
===================================================================
Target-trial emulation of the peri-operative initial crystalloid choice.

Time zero  : surgery_start_datetime (from or_cases_cleaned)
Exposure   : first standalone LR or NS given WITHIN the OR case
             (taken_datetime in [surgery_start - 2h, surgery_stop + 1h])
Outcome    : KDIGO-defined AKI within 7 days of surgery start
             (aki_surgery_7d_without_race_correction from outcome_final)
Baseline Cr: reference creatinine from Variable_Generation

Run:
    python3 extract_cohort.py --out cohort_idealist.parquet
"""

import argparse, glob, os, sys, warnings, shutil, tempfile
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA = ("/orange/prismap-data-core/AIDAIR/Cleaned_files/EHR/"
        "IDEALIST_RETRO_CLEANED_2025/miscellaneous/new_primary_logic_recleaning/cleaned")

# ── LR / NS identification ─────────────────────────────────────────────────────
# Use display name, NOT iv_class (iv_class mislabels LR as 'saline' in IDEALIST)
LR_PAT  = r'(?i)^LR[\s,]|^LACTATED|RINGER'
NS_PAT  = r'(?i)^0\.9\s*%\s*NACL|^NORMAL SALINE|^NS[\s,]'
# Exclude drug carriers (antibiotic-in-NS, etc.) and flushes/nebulizers
EXCL_PAT = r'(?i)\bFLUSH\b|NEBULIZER|3\s*%|7\s*%| IN 0\.9| IN LR| IN LACTATED|IRRIGATION|HEPARIN|PRIMING'

SEED = 42
np.random.seed(SEED)


def read_batched(pattern, columns=None):
    """Read all parquet batches matching pattern, return concatenated DataFrame."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files for pattern: {pattern}")
    frames = []
    for f in files:
        pf = pq.ParquetFile(f)
        for i in range(pf.num_row_groups):
            rg = pf.read_row_group(i, columns=columns)
            frames.append(rg.to_pandas())
    return pd.concat(frames, ignore_index=True)


def load_or_cases():
    """Load OR cases: one row per surgical case with timing."""
    print("Loading OR cases…")
    cols = ["patient_deiden_id", "encounter_deiden_id_primary",
            "surgery_start_datetime", "surgery_stop_datetime",
            "or_case_num", "sched_service", "sched_post_op_location"]
    df = read_batched(f"{DATA}/or_cases_cleaned_*.parquet", columns=cols)
    df["surgery_start_datetime"] = pd.to_datetime(df["surgery_start_datetime"], utc=False, errors="coerce")
    df["surgery_stop_datetime"]  = pd.to_datetime(df["surgery_stop_datetime"],  utc=False, errors="coerce")
    # Drop cases with no surgery time
    df = df.dropna(subset=["surgery_start_datetime"])
    # One row per case (or_case_num), keeping earliest start if duplicated
    df = df.sort_values("surgery_start_datetime").drop_duplicates(
        subset=["encounter_deiden_id_primary", "or_case_num"], keep="first"
    )
    print(f"  OR cases: {len(df):,}")
    return df


def load_fluids():
    """Load meds: standalone LR or NS administrations with timestamps."""
    print("Loading fluid medications…")
    cols = ["patient_deiden_id", "encounter_deiden_id_primary",
            "med_order_display_name", "taken_datetime",
            "fluid_volume", "imputed_fluid_bolus", "iv_class"]
    df = read_batched(f"{DATA}/meds_cleaned_*.parquet", columns=cols)
    df["taken_datetime"] = pd.to_datetime(df["taken_datetime"], utc=False, errors="coerce")
    df = df.dropna(subset=["taken_datetime"])

    # Classify fluid type by display name
    is_lr   = df.med_order_display_name.str.contains(LR_PAT,   na=False, regex=True)
    is_ns   = df.med_order_display_name.str.contains(NS_PAT,   na=False, regex=True)
    is_excl = df.med_order_display_name.str.contains(EXCL_PAT, na=False, regex=True)

    fluids = df[(is_lr | is_ns) & ~is_excl].copy()
    fluids["fluid_type"] = np.where(is_lr[fluids.index], "LR", "NS")
    print(f"  Fluid records (LR+NS standalone): {len(fluids):,}  "
          f"LR={fluids.fluid_type.eq('LR').sum():,}  NS={fluids.fluid_type.eq('NS').sum():,}")
    return fluids[["patient_deiden_id", "encounter_deiden_id_primary",
                   "taken_datetime", "fluid_type", "fluid_volume"]]


def load_outcomes():
    """Load outcome_final: AKI flag + hospital mortality."""
    print("Loading outcome_final…")
    cols = ["patient_deiden_id", "encounter_deiden_id_primary",
            "aki_surgery_7d_without_race_correction",
            "aki_surgery_3d_without_race_correction",
            "hospital_mortality",
            "postop_icu_yn",
            "surgery_death_7d",
            "admit_datetime", "discharge_datetime"]
    df = read_batched(f"{DATA}/Outcome_Generation/output/outcome_final_*.parquet", columns=cols)
    print(f"  Outcomes: {len(df):,}")
    return df


def load_variables():
    """Load Variable_Generation: 408 pre-built periop confounders."""
    print("Loading Variable_Generation (confounders)…")
    # Key confounder columns we need
    cov_cols = [
        "patient_deiden_id", "encounter_deiden_id_primary", "or_case_num",
        # demographics
        "sex", "race", "ethnicity", "smoking_status",
        "age", "bmi",
        # creatinine (pre-op window 0-7d and 8-365d before surgery)
        "serum_creatinine_min_0_7", "serum_creatinine_max_0_7",
        "serum_creatinine_mean_0_7", "count_serum_creatinine_0_7",
        "serum_creatinine_min_8_365", "serum_creatinine_max_8_365",
        "serum_creatinine_mean_8_365",
        "reference creatinine", "egfr",
        # comorbidities (Elixhauser-style flags in Variable_Generation)
        "diabetes", "ckd", "liverd", "ichf", "icancer", "icpd",
        "icvd", "imi", "ipvd", "htn_c", "hypothy", "depress",
        "anemia", "obese", "wghtloss", "coag", "alc_drug",
        # nephrotoxic meds
        "nephrotoxic", "nsaids", "aminoglycosides", "vancomycin",
        # illness severity / surgical context
        "emergent", "pressors_inotropes", "surgery_type", "site",
        "specialty", "payer", "admit_year", "admit_month",
        "postop_loc", "sched_post_op_location",
        "time_to_surgery_hours", "days_since_last_surgery",
        "prev_surgery_7d",
        # labs (pre-op)
        "serum_albumin_min_0_7", "serum_albumin_mean_0_7",
        "serum_lactate_min_0_7", "serum_lactate_max_0_7",
        "wbc_mean_0_7", "plt_mean_0_7", "hgb_mean_0_7",
        "serum_creatinine_min_8_365",
        # surgery timing
        "surgery_start_datetime", "surgery_stop_datetime",
        "night_admission",
    ]
    # Read available columns (Variable_Generation has 408; not all may match)
    files = sorted(glob.glob(f"{DATA}/Variable_Generation/output/all_generated_variables_*.parquet"))
    if not files:
        raise FileNotFoundError("Variable_Generation output not found")
    # Get actual columns from first file
    actual = [c for c,_ in [(f.name, str(f.type))
                             for f in pq.ParquetFile(files[0]).schema_arrow]]
    use_cols = [c for c in cov_cols if c in actual]
    missing  = [c for c in cov_cols if c not in actual]
    if missing:
        print(f"  [NOTE] Variable_Generation missing cols (will impute): {missing}")

    frames = []
    for f in files:
        pf = pq.ParquetFile(f)
        for i in range(pf.num_row_groups):
            frames.append(pf.read_row_group(i, columns=use_cols).to_pandas())
    df = pd.concat(frames, ignore_index=True)
    print(f"  Variable_Generation rows: {len(df):,}, cols used: {len(use_cols)}")
    return df


def merge_cohort(or_cases, fluids, outcomes, variables):
    """
    Build analysis cohort:
    1. For each OR case find the FIRST standalone LR or NS given during surgery
       (taken_datetime in [surgery_start - 2h, surgery_stop + 1h])
    2. Join outcomes (AKI, mortality) on encounter_deiden_id_primary
    3. Join confounders from Variable_Generation on encounter+or_case
    """
    print("Merging cohort…")

    # ── Step 1: match fluids to OR cases ──────────────────────────────────────
    merged = or_cases.merge(
        fluids, on=["patient_deiden_id", "encounter_deiden_id_primary"], how="inner"
    )
    # Keep only fluids given during (or just before/after) surgery
    window_start = merged["surgery_start_datetime"] - pd.Timedelta(hours=2)
    window_end   = merged["surgery_stop_datetime"]  + pd.Timedelta(hours=1)
    in_window = (merged["taken_datetime"] >= window_start) & \
                (merged["taken_datetime"] <= window_end)
    merged = merged[in_window].copy()

    # For each OR case, keep only the FIRST crystalloid (time zero)
    merged = merged.sort_values("taken_datetime")
    first_fluid = merged.drop_duplicates(
        subset=["encounter_deiden_id_primary", "or_case_num"], keep="first"
    )
    first_fluid = first_fluid.rename(columns={"taken_datetime": "time_zero"})
    print(f"  Cases with intraop LR or NS: {len(first_fluid):,}  "
          f"LR={first_fluid.fluid_type.eq('LR').sum():,}  "
          f"NS={first_fluid.fluid_type.eq('NS').sum():,}")

    # ── Step 2: join outcomes ──────────────────────────────────────────────────
    cohort = first_fluid.merge(outcomes, on=["patient_deiden_id", "encounter_deiden_id_primary"],
                               how="inner")
    print(f"  After joining outcomes: {len(cohort):,}")

    # ── Step 3: join confounders ───────────────────────────────────────────────
    # or_case_num may be int in Variable_Generation vs string in or_cases; coerce
    variables = variables.copy()
    if "or_case_num" in variables.columns and "or_case_num" in cohort.columns:
        variables["or_case_num"] = variables["or_case_num"].astype(str)
        cohort["or_case_num"]    = cohort["or_case_num"].astype(str)
        cohort = cohort.merge(variables,
                              on=["encounter_deiden_id_primary", "or_case_num"],
                              how="left", suffixes=("", "_vg"))
    else:
        cohort = cohort.merge(variables.drop(columns=["or_case_num"], errors="ignore"),
                              on="encounter_deiden_id_primary", how="left",
                              suffixes=("", "_vg"))
    print(f"  After joining confounders: {len(cohort):,}")
    return cohort


def apply_eligibility(cohort):
    """
    Eligibility criteria:
    - reference creatinine present and ≤ 4.0 mg/dL (consistent with MIMIC)
    - AKI outcome available (aki_surgery_7d not null)
    - Adult (assume all IDEALIST surgical patients are adult)
    """
    print("Applying eligibility criteria…")
    n0 = len(cohort)

    # Rename reference creatinine column (has a space)
    if "reference creatinine" in cohort.columns:
        cohort = cohort.rename(columns={"reference creatinine": "baseline_creatinine"})

    # Coerce to numeric
    cohort["baseline_creatinine"] = pd.to_numeric(
        cohort["baseline_creatinine"], errors="coerce"
    )
    cohort["aki_surgery_7d"] = pd.to_numeric(
        cohort["aki_surgery_7d_without_race_correction"], errors="coerce"
    )

    # Require baseline creatinine
    n_no_cr = cohort["baseline_creatinine"].isna().sum()
    cohort = cohort.dropna(subset=["baseline_creatinine"])
    print(f"  Missing baseline creatinine excluded: {n_no_cr:,}")

    # Require AKI outcome
    n_no_aki = cohort["aki_surgery_7d"].isna().sum()
    cohort = cohort.dropna(subset=["aki_surgery_7d"])
    print(f"  Missing AKI outcome excluded: {n_no_aki:,}")

    # Prevalent AKI / advanced CKD proxy: Cr > 4.0
    n_prev = (cohort["baseline_creatinine"] > 4.0).sum()
    cohort = cohort[cohort["baseline_creatinine"] <= 4.0].copy()
    print(f"  Prevalent AKI proxy (Cr > 4.0) excluded: {n_prev:,}")

    print(f"  Eligible cohort: {len(cohort):,} (from {n0:,})")
    return cohort.reset_index(drop=True)


def main(out_path):
    or_cases  = load_or_cases()
    fluids    = load_fluids()
    outcomes  = load_outcomes()
    variables = load_variables()

    cohort = merge_cohort(or_cases, fluids, outcomes, variables)
    cohort = apply_eligibility(cohort)

    print(f"\nFinal cohort: {len(cohort):,} rows, {cohort.shape[1]} cols")
    print(f"LR: {cohort.fluid_type.eq('LR').sum():,}  "
          f"NS: {cohort.fluid_type.eq('NS').sum():,}")
    print(f"AKI rate LR: {cohort[cohort.fluid_type=='LR'].aki_surgery_7d.mean():.3f}  "
          f"NS: {cohort[cohort.fluid_type=='NS'].aki_surgery_7d.mean():.3f}")

    cohort.to_parquet(out_path, index=False)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="cohort_idealist.parquet")
    args = p.parse_args()
    main(args.out)
