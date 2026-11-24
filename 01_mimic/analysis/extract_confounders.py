"""
Confounder enrichment: reads raw MIMIC-IV and joins new columns onto the analytic parquet.

New columns added:
  cr_48h           - max serum Cr in [t0, t0+48h] (KDIGO 48h absolute criterion)
  ckd              - Elixhauser CKD flag (whole-hosp ICD codes; pre-t0 audit needed)
  liver_disease    - Elixhauser liver disease flag
  heart_failure    - Elixhauser heart failure flag
  sepsis           - sepsis flag
  cancer           - cancer flag
  diabetes         - diabetes flag
  electrolyte_disorder - electrolyte disorder flag
  gcs_min_pre      - min GCS in 24h before t0 (severity proxy)
  hr_max_pre       - max heart rate in 24h before t0
  map_min_pre      - min MAP in 24h before t0
  spo2_min_pre     - min SpO2 in 24h before t0
  vasopressor_pre  - 1 if any vasopressor started before t0
  cr_slope_pre     - linear slope of Cr over 7d before baseline (mg/dL/day)
  n_cr_pre         - number of pre-baseline Cr measurements
  year             - calendar year of t0

Usage:
  python3 extract_confounders.py \\
    --parquet parquet/cohort_analytic.parquet \\
    --mimic /path/to/mimiciv/3.1 \\
    --out parquet/cohort_enriched.parquet
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

SEED = 42
np.random.seed(SEED)

# ── ICD code maps ─────────────────────────────────────────────────────────────
# NOTE: These are whole-hospitalization codes; pre-time-zero audit is needed
# before treating them as confirmed pre-exposure comorbidities.
COMORBIDITY_ICD = {
    "ckd": {
        9:  ["585", "V420", "V451", "V56"],
        10: ["N18", "Z940", "Z992"],
    },
    "liver_disease": {
        9:  ["571", "572", "070"],
        10: ["K70", "K71", "K72", "K73", "K74", "B18", "B19"],
    },
    "heart_failure": {
        9:  ["428", "39891", "40201", "40211", "40291"],
        10: ["I50", "I110", "I130", "I132"],
    },
    "sepsis": {
        9:  ["038", "9959", "78552"],
        10: ["A40", "A41", "R652"],
    },
    "cancer": {
        9:  [str(c) for c in range(140, 196)] + [str(c) for c in range(200, 209)],
        10: ["C"],
    },
    "diabetes": {
        9:  ["250"],
        10: ["E10", "E11", "E12", "E13", "E14"],
    },
    "electrolyte_disorder": {
        9:  ["276"],
        10: ["E86", "E87"],
    },
}

# ── vital/lab itemids ─────────────────────────────────────────────────────────
CREAT_ITEMID   = 50912
GCS_ITEMID     = 220739
HR_ITEMID      = 220045
MAP_ITEMID     = 220052
SPO2_ITEMID    = 220277
VASO_ITEMIDS   = [221906, 221289, 221662, 222315, 221749]
CHART_ITEMIDS  = [GCS_ITEMID, HR_ITEMID, MAP_ITEMID, SPO2_ITEMID]


def load_cohort(parquet_path):
    df = pd.read_parquet(parquet_path)
    for col in ["fluid_starttime", "icu_intime", "baseline_cr_time"]:
        df[col] = pd.to_datetime(df[col])
    print(f"  Cohort loaded: {len(df):,} rows")
    return df


# ── A. 48h creatinine ─────────────────────────────────────────────────────────
def extract_cr_48h(cohort, mimic_hosp):
    print("[A] 48h creatinine from labevents ...")
    lab_path = os.path.join(mimic_hosp, "labevents.csv")
    # Read in chunks — file is very large (~158M lines)
    stay_map = cohort[["subject_id", "stay_id", "fluid_starttime"]].copy()
    stay_map["t48"] = stay_map["fluid_starttime"] + pd.Timedelta(hours=48)

    chunks = []
    chunksize = 2_000_000
    for chunk in pd.read_csv(
        lab_path,
        usecols=["subject_id", "itemid", "charttime", "valuenum"],
        dtype={"subject_id": "int32", "itemid": "int32"},
        parse_dates=["charttime"],
        chunksize=chunksize,
    ):
        c = chunk[chunk["itemid"] == CREAT_ITEMID].copy()
        if c.empty:
            continue
        c = c.merge(stay_map, on="subject_id", how="inner")
        c = c[
            (c["charttime"] >= c["fluid_starttime"]) &
            (c["charttime"] <= c["t48"]) &
            c["valuenum"].notna()
        ]
        chunks.append(c[["stay_id", "valuenum"]])

    if not chunks:
        print("  No 48h creatinine found.")
        cohort["cr_48h"] = np.nan
        return cohort

    cr48 = pd.concat(chunks).groupby("stay_id")["valuenum"].max().reset_index()
    cr48.columns = ["stay_id", "cr_48h"]
    cohort = cohort.merge(cr48, on="stay_id", how="left")
    n = cohort["cr_48h"].notna().sum()
    print(f"  cr_48h: {n:,} / {len(cohort):,} non-null")
    return cohort


# ── B. Comorbidities from diagnoses_icd ───────────────────────────────────────
def extract_comorbidities(cohort, mimic_hosp):
    print("[B] Comorbidities from diagnoses_icd ...")
    dx_path = os.path.join(mimic_hosp, "diagnoses_icd.csv")
    dx = pd.read_csv(dx_path, usecols=["hadm_id", "icd_code", "icd_version"],
                     dtype={"hadm_id": "int32", "icd_version": "int8", "icd_code": str})
    dx["icd_code"] = dx["icd_code"].str.strip().str.upper()

    hadm_ids = set(cohort["hadm_id"].unique())
    dx = dx[dx["hadm_id"].isin(hadm_ids)]

    result = cohort[["hadm_id"]].drop_duplicates()
    for name, maps in COMORBIDITY_ICD.items():
        flags = []
        for version, prefixes in maps.items():
            sub = dx[dx["icd_version"] == version]
            mask = sub["icd_code"].str.startswith(tuple(p.upper() for p in prefixes))
            flagged = sub.loc[mask, "hadm_id"].unique()
            flags.append(flagged)
        all_flagged = np.concatenate(flags) if flags else np.array([], dtype=int)
        result[name] = result["hadm_id"].isin(all_flagged).astype(int)

    cohort = cohort.merge(result, on="hadm_id", how="left")
    for name in COMORBIDITY_ICD:
        cohort[name] = cohort[name].fillna(0).astype(int)
        print(f"  {name}: {cohort[name].sum():,} flagged ({cohort[name].mean()*100:.1f}%)")
    return cohort


# ── C. Vitals/severity from ICU chartevents ───────────────────────────────────
def extract_vitals(cohort, mimic_icu):
    print("[C] Vitals/severity from ICU chartevents ...")
    chart_path = os.path.join(mimic_icu, "chartevents.csv")

    stay_map = cohort[["stay_id", "fluid_starttime"]].copy()
    stay_map["pre_start"] = stay_map["fluid_starttime"] - pd.Timedelta(hours=24)
    stay_ids = set(stay_map["stay_id"].unique())

    agg = {sid: {"gcs": [], "hr": [], "map": [], "spo2": []} for sid in stay_ids}

    chunksize = 3_000_000
    for chunk in pd.read_csv(
        chart_path,
        usecols=["stay_id", "itemid", "charttime", "valuenum"],
        dtype={"stay_id": "int32", "itemid": "int32"},
        parse_dates=["charttime"],
        chunksize=chunksize,
    ):
        c = chunk[
            chunk["itemid"].isin(CHART_ITEMIDS) &
            chunk["stay_id"].isin(stay_ids) &
            chunk["valuenum"].notna()
        ].copy()
        if c.empty:
            continue
        c = c.merge(stay_map, on="stay_id", how="inner")
        c = c[(c["charttime"] >= c["pre_start"]) & (c["charttime"] < c["fluid_starttime"])]
        for _, row in c.iterrows():
            sid = row["stay_id"]
            if row["itemid"] == GCS_ITEMID:
                agg[sid]["gcs"].append(row["valuenum"])
            elif row["itemid"] == HR_ITEMID:
                agg[sid]["hr"].append(row["valuenum"])
            elif row["itemid"] == MAP_ITEMID:
                agg[sid]["map"].append(row["valuenum"])
            elif row["itemid"] == SPO2_ITEMID:
                agg[sid]["spo2"].append(row["valuenum"])

    rows = []
    for sid, vals in agg.items():
        rows.append({
            "stay_id": sid,
            "gcs_min_pre":  min(vals["gcs"])  if vals["gcs"]  else np.nan,
            "hr_max_pre":   max(vals["hr"])   if vals["hr"]   else np.nan,
            "map_min_pre":  min(vals["map"])  if vals["map"]  else np.nan,
            "spo2_min_pre": min(vals["spo2"]) if vals["spo2"] else np.nan,
        })
    vitals = pd.DataFrame(rows)
    cohort = cohort.merge(vitals, on="stay_id", how="left")
    for col in ["gcs_min_pre", "hr_max_pre", "map_min_pre", "spo2_min_pre"]:
        n = cohort[col].notna().sum()
        print(f"  {col}: {n:,} / {len(cohort):,} non-null")
    return cohort


# ── D. Vasopressor use ────────────────────────────────────────────────────────
def extract_vasopressors(cohort, mimic_icu):
    print("[D] Vasopressor use from inputevents ...")
    ie_path = os.path.join(mimic_icu, "inputevents.csv")
    stay_map = cohort[["stay_id", "fluid_starttime"]].copy()
    stay_ids = set(stay_map["stay_id"].unique())

    chunks = []
    for chunk in pd.read_csv(
        ie_path,
        usecols=["stay_id", "itemid", "starttime"],
        dtype={"stay_id": "int32", "itemid": "int32"},
        parse_dates=["starttime"],
        chunksize=1_000_000,
    ):
        c = chunk[
            chunk["itemid"].isin(VASO_ITEMIDS) &
            chunk["stay_id"].isin(stay_ids)
        ]
        chunks.append(c[["stay_id", "starttime"]])

    if not chunks:
        cohort["vasopressor_pre"] = 0
        return cohort

    ie = pd.concat(chunks)
    ie = ie.merge(stay_map, on="stay_id", how="inner")
    ie = ie[ie["starttime"] < ie["fluid_starttime"]]
    vaso = ie.groupby("stay_id").size().reset_index(name="n_vaso")
    vaso["vasopressor_pre"] = 1
    cohort = cohort.merge(vaso[["stay_id", "vasopressor_pre"]], on="stay_id", how="left")
    cohort["vasopressor_pre"] = cohort["vasopressor_pre"].fillna(0).astype(int)
    n = cohort["vasopressor_pre"].sum()
    print(f"  vasopressor_pre: {n:,} / {len(cohort):,} ({n/len(cohort)*100:.1f}%)")
    return cohort


# ── E. Creatinine trajectory ──────────────────────────────────────────────────
def extract_cr_trajectory(cohort, mimic_hosp):
    print("[E] Pre-baseline creatinine trajectory from labevents ...")
    lab_path = os.path.join(mimic_hosp, "labevents.csv")
    stay_map = cohort[["subject_id", "stay_id", "icu_intime", "baseline_cr_time"]].copy()
    stay_map["traj_start"] = stay_map["icu_intime"] - pd.Timedelta(days=7)

    all_cr = []
    for chunk in pd.read_csv(
        lab_path,
        usecols=["subject_id", "itemid", "charttime", "valuenum"],
        dtype={"subject_id": "int32", "itemid": "int32"},
        parse_dates=["charttime"],
        chunksize=2_000_000,
    ):
        c = chunk[chunk["itemid"] == CREAT_ITEMID].copy()
        if c.empty:
            continue
        c = c.merge(stay_map, on="subject_id", how="inner")
        c = c[
            (c["charttime"] >= c["traj_start"]) &
            (c["charttime"] <= c["baseline_cr_time"]) &
            c["valuenum"].notna()
        ]
        all_cr.append(c[["stay_id", "charttime", "valuenum"]])

    if not all_cr:
        cohort["cr_slope_pre"] = np.nan
        cohort["n_cr_pre"] = 0
        return cohort

    cr = pd.concat(all_cr)

    def compute_slope(grp):
        grp = grp.sort_values("charttime")
        n = len(grp)
        if n < 2:
            return pd.Series({"cr_slope_pre": np.nan, "n_cr_pre": n})
        t = (grp["charttime"] - grp["charttime"].iloc[0]).dt.total_seconds() / 86400.0
        y = grp["valuenum"].values
        if t.std() == 0:  # all measurements at same time — slope undefined
            return pd.Series({"cr_slope_pre": np.nan, "n_cr_pre": n})
        try:
            slope = np.polyfit(t, y, 1)[0]
        except np.linalg.LinAlgError:
            slope = np.nan
        return pd.Series({"cr_slope_pre": slope, "n_cr_pre": n})

    traj = cr.groupby("stay_id").apply(compute_slope).reset_index()
    cohort = cohort.merge(traj, on="stay_id", how="left")
    cohort["n_cr_pre"] = cohort["n_cr_pre"].fillna(0).astype(int)
    n = cohort["cr_slope_pre"].notna().sum()
    print(f"  cr_slope_pre: {n:,} / {len(cohort):,} non-null")
    print(f"  n_cr_pre: mean={cohort['n_cr_pre'].mean():.1f}")
    return cohort


# ── F. Calendar year ──────────────────────────────────────────────────────────
def add_calendar_year(cohort):
    cohort["year"] = pd.to_datetime(cohort["fluid_starttime"]).dt.year
    print(f"[F] year: {cohort['year'].value_counts().sort_index().to_dict()}")
    return cohort


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--mimic",   required=True)
    parser.add_argument("--out",     required=True)
    args = parser.parse_args()

    mimic_hosp = os.path.join(args.mimic, "hosp")
    mimic_icu  = os.path.join(args.mimic, "icu")

    print("=" * 60)
    print("Confounder enrichment pipeline")
    print("=" * 60)

    cohort = load_cohort(args.parquet)

    cohort = extract_cr_48h(cohort, mimic_hosp)
    cohort = extract_comorbidities(cohort, mimic_hosp)
    cohort = extract_vasopressors(cohort, mimic_icu)
    cohort = add_calendar_year(cohort)

    # Vitals and Cr trajectory are slow (huge files); run them too
    cohort = extract_vitals(cohort, mimic_icu)
    cohort = extract_cr_trajectory(cohort, mimic_hosp)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    cohort.to_parquet(args.out, index=False)
    print(f"\nEnriched parquet saved → {args.out}")
    print(f"Shape: {cohort.shape}")

    new_cols = ["cr_48h", "ckd", "liver_disease", "heart_failure", "sepsis",
                "cancer", "diabetes", "electrolyte_disorder",
                "gcs_min_pre", "hr_max_pre", "map_min_pre", "spo2_min_pre",
                "vasopressor_pre", "cr_slope_pre", "n_cr_pre", "year"]
    print("\nSummary of new columns:")
    for col in new_cols:
        if col in cohort.columns:
            nn = cohort[col].notna().sum()
            print(f"  {col:<25} {nn:>7,} non-null  ({nn/len(cohort)*100:.1f}%)")
        else:
            print(f"  {col:<25} MISSING")


if __name__ == "__main__":
    main()
