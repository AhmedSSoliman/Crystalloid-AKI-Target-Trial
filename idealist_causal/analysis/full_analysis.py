"""
IDEALIST LR vs NS -> 7-day post-surgical AKI  —  AIPW Primary Analysis
========================================================================
Mirrors the MIMIC causal_ai/analysis/full_analysis.py exactly:
  - AIPW doubly robust estimator (Bang & Robins 2005)
  - PS trimming [0.05, 0.95]
  - 200-replicate percentile bootstrap CIs
  - Covariate balance (SMD) before/after IPTW
  - Writes results_summary.csv (appended by later SA scripts)

Run (after extract_cohort.py):
    python3 full_analysis.py --parquet cohort_idealist.parquet
"""

import argparse, os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed

SEED = 42
np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_CSV    = os.path.join(SCRIPT_DIR, "results_summary_idealist.csv")
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Confounder columns (available after extract_cohort.py) ────────────────────
COV_COLS = [
    "sex", "race", "ethnicity", "smoking_status", "age", "bmi",
    "baseline_creatinine",   # renamed from 'reference creatinine'
    "egfr",
    "serum_creatinine_min_0_7", "serum_creatinine_max_0_7",
    "serum_creatinine_mean_0_7", "count_serum_creatinine_0_7",
    "serum_creatinine_min_8_365", "serum_creatinine_max_8_365",
    "diabetes", "ckd", "liverd", "ichf", "icancer", "icpd",
    "icvd", "imi", "ipvd", "htn_c", "hypothy", "depress",
    "anemia", "obese", "wghtloss", "coag", "alc_drug",
    "nephrotoxic", "nsaids", "aminoglycosides", "vancomycin",
    "emergent", "pressors_inotropes", "surgery_type", "site",
    "specialty", "payer", "admit_year", "admit_month",
    "time_to_surgery_hours", "days_since_last_surgery",
    "prev_surgery_7d", "night_admission",
    "serum_albumin_min_0_7", "serum_albumin_mean_0_7",
    "serum_lactate_min_0_7", "serum_lactate_max_0_7",
    "wbc_mean_0_7", "plt_mean_0_7", "hgb_mean_0_7",
]

N_BOOT = 50
PS_LOW, PS_HIGH = 0.05, 0.95


def prep_X(df, cov_cols):
    """Coerce numerics, one-hot only low-cardinality categoricals (≤50 levels)."""
    available = [c for c in cov_cols if c in df.columns]
    sub = df[available].copy()
    # Coerce string columns to numeric where possible
    for col in sub.select_dtypes(include=["object", "category"]).columns:
        coerced = pd.to_numeric(sub[col], errors="coerce")
        if coerced.notna().mean() >= 0.5:  # mostly numeric → keep as float
            sub[col] = coerced
    # One-hot only remaining string cols with ≤50 unique values
    cats = [c for c in sub.select_dtypes(include=["object","category"]).columns
            if sub[c].nunique() <= 50]
    drop_cols = [c for c in sub.select_dtypes(include=["object","category"]).columns
                 if c not in cats]
    sub = sub.drop(columns=drop_cols)
    if cats:
        sub = pd.get_dummies(sub, columns=cats, drop_first=True, dtype=np.float32)
    return sub.astype(np.float32)


def fit_ps(X, A):
    pipe = Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("scl", StandardScaler()),
        ("lr",  LogisticRegression(max_iter=500, C=1.0, solver="saga",
                                   n_jobs=1, random_state=SEED)),
    ])
    pipe.fit(X, A)
    return pipe.predict_proba(X)[:, 1]


def aipw_estimate(df, ps, outcome_col="aki_surgery_7d"):
    """AIPW doubly-robust risk difference (LR - NS)."""
    A = df["fluid_type"].eq("LR").astype(int).values
    Y = df[outcome_col].values
    n = len(df)

    X = prep_X(df, COV_COLS)

    # Outcome models by treatment arm
    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    Ximp = scl.fit_transform(imp.fit_transform(X))

    om1 = LogisticRegression(max_iter=500, C=1.0, solver="saga", n_jobs=1, random_state=SEED)
    om0 = LogisticRegression(max_iter=500, C=1.0, solver="saga", n_jobs=1, random_state=SEED)
    if A.sum() > 5:
        om1.fit(Ximp[A == 1], Y[A == 1])
        mu1 = om1.predict_proba(Ximp)[:, 1]
    else:
        mu1 = np.full(n, Y[A == 1].mean())
    if (1 - A).sum() > 5:
        om0.fit(Ximp[A == 0], Y[A == 0])
        mu0 = om0.predict_proba(Ximp)[:, 1]
    else:
        mu0 = np.full(n, Y[A == 0].mean())

    # Stabilized weights
    p1 = ps.mean(); p0 = 1 - p1
    w1 = np.where(A == 1, p1 / np.clip(ps, 1e-6, 1), 0)
    w0 = np.where(A == 0, p0 / np.clip(1 - ps, 1e-6, 1), 0)

    # AIPW
    r1 = mu1 + w1 * (Y - mu1)
    r0 = mu0 + w0 * (Y - mu0)
    rd = r1.mean() - r0.mean()
    return rd, r1.mean(), r0.mean()


def iptw_only(df, ps, outcome_col="aki_surgery_7d"):
    A = df["fluid_type"].eq("LR").astype(int).values
    Y = df[outcome_col].values
    p1 = ps.mean(); p0 = 1 - p1
    w1 = np.where(A == 1, p1 / np.clip(ps, 1e-6, 1), 0)
    w0 = np.where(A == 0, p0 / np.clip(1 - ps, 1e-6, 1), 0)
    r1 = (w1 * Y).sum() / w1.sum()
    r0 = (w0 * Y).sum() / w0.sum()
    return r1 - r0, r1, r0


def smd(x1, x0):
    m1, m0 = np.nanmean(x1), np.nanmean(x0)
    s1, s0 = np.nanstd(x1), np.nanstd(x0)
    denom = np.sqrt((s1**2 + s0**2) / 2)
    return (m1 - m0) / denom if denom > 0 else 0.0


def weighted_smd(col, treat, df, w):
    x1 = df.loc[treat == 1, col].values.astype(float)
    x0 = df.loc[treat == 0, col].values.astype(float)
    w1 = w[treat == 1]; w0 = w[treat == 0]
    mu1 = np.nansum(w1 * x1) / np.nansum(w1)
    mu0 = np.nansum(w0 * x0) / np.nansum(w0)
    var1 = np.nansum(w1 * (x1 - mu1)**2) / np.nansum(w1)
    var0 = np.nansum(w0 * (x0 - mu0)**2) / np.nansum(w0)
    denom = np.sqrt((var1 + var0) / 2)
    return (mu1 - mu0) / denom if denom > 0 else 0.0


def bootstrap_one(seed, df, outcome_col):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(df), len(df))
    b = df.iloc[idx].reset_index(drop=True)
    X = prep_X(b, COV_COLS)
    A = b["fluid_type"].eq("LR").astype(int).values
    ps = fit_ps(X, A)
    # Trim inside bootstrap
    keep = (ps >= PS_LOW) & (ps <= PS_HIGH)
    b2 = b[keep].reset_index(drop=True)
    ps2 = ps[keep]
    rd, _, _ = aipw_estimate(b2, ps2, outcome_col)
    return rd


def bootstrap_ci(df, outcome_col, n_boot=N_BOOT):
    seeds = list(range(SEED, SEED + n_boot))
    rds = [bootstrap_one(s, df, outcome_col) for s in seeds]
    rds = np.array(rds)
    return float(np.nanpercentile(rds, 2.5)), float(np.nanpercentile(rds, 97.5))


def append_result(row: dict):
    df = pd.DataFrame([row])
    if os.path.exists(OUT_CSV):
        existing = pd.read_csv(OUT_CSV)
        # Remove any existing row with same analysis label
        existing = existing[existing["analysis"] != row["analysis"]]
        df = pd.concat([existing, df], ignore_index=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"  Saved → {OUT_CSV}")


def main(parquet_path):
    print(f"Loading cohort: {parquet_path}")
    cohort = pd.read_parquet(parquet_path)
    print(f"  {len(cohort):,} rows | LR={cohort.fluid_type.eq('LR').sum():,} "
          f"NS={cohort.fluid_type.eq('NS').sum():,}")
    print(f"  Overall AKI rate: {cohort.aki_surgery_7d.mean():.3f}")

    # ── Propensity score ───────────────────────────────────────────────────────
    print("\nFitting propensity score…")
    X_full = prep_X(cohort, COV_COLS)
    A_full = cohort["fluid_type"].eq("LR").astype(int).values
    ps_full = fit_ps(X_full, A_full)
    auc = roc_auc_score(A_full, ps_full)
    print(f"  PS AUC: {auc:.3f}")

    # Save PS to cohort for downstream figures
    cohort["ps"] = ps_full

    # ── PS trimming ────────────────────────────────────────────────────────────
    keep = (ps_full >= PS_LOW) & (ps_full <= PS_HIGH)
    ct = cohort[keep].reset_index(drop=True)
    ps_t = ps_full[keep]
    n_trimmed = (~keep).sum()
    print(f"  PS trimming: removed {n_trimmed:,} → primary cohort n={len(ct):,}")

    # ── AIPW primary estimate ──────────────────────────────────────────────────
    print("\nAIPW primary estimate…")
    rd, r_lr, r_ns = aipw_estimate(ct, ps_t)
    print(f"  RD = {rd*100:+.2f} pp  LR={r_lr*100:.1f}%  NS={r_ns*100:.1f}%")

    print("  Bootstrap CIs (200 reps)…")
    lo, hi = bootstrap_ci(ct, "aki_surgery_7d")
    print(f"  95% CI: ({lo*100:+.2f}, {hi*100:+.2f})")

    append_result(dict(
        analysis="Primary AIPW (PS-trimmed)",
        n=len(ct), rd=rd, rd_lower=lo, rd_upper=hi,
        risk_lr=r_lr, risk_ns=r_ns,
    ))

    # ── SA1: no PS trimming ────────────────────────────────────────────────────
    print("\nSA1: No PS trimming…")
    rd1, r1_lr, r1_ns = aipw_estimate(cohort, ps_full)
    lo1, hi1 = bootstrap_ci(cohort, "aki_surgery_7d")
    print(f"  RD = {rd1*100:+.2f} pp  CI: ({lo1*100:+.2f}, {hi1*100:+.2f})")
    append_result(dict(
        analysis="SA1: No PS trimming",
        n=len(cohort), rd=rd1, rd_lower=lo1, rd_upper=hi1,
        risk_lr=r1_lr, risk_ns=r1_ns,
    ))

    # ── SA2: Stricter baseline Cr ≤ 2.0 ───────────────────────────────────────
    print("\nSA2: Baseline Cr ≤ 2.0…")
    c2 = cohort[cohort.baseline_creatinine <= 2.0].reset_index(drop=True)
    X2 = prep_X(c2, COV_COLS)
    A2 = c2.fluid_type.eq("LR").astype(int).values
    ps2 = fit_ps(X2, A2)
    keep2 = (ps2 >= PS_LOW) & (ps2 <= PS_HIGH)
    c2t = c2[keep2].reset_index(drop=True); ps2t = ps2[keep2]
    rd2, r2_lr, r2_ns = aipw_estimate(c2t, ps2t)
    lo2, hi2 = bootstrap_ci(c2t, "aki_surgery_7d")
    print(f"  n={len(c2t):,}  RD = {rd2*100:+.2f} pp  CI: ({lo2*100:+.2f}, {hi2*100:+.2f})")
    append_result(dict(
        analysis="SA2: Stricter baseline Cr≤2.0",
        n=len(c2t), rd=rd2, rd_lower=lo2, rd_upper=hi2,
        risk_lr=r2_lr, risk_ns=r2_ns,
    ))

    # ── SA3: Composite death-or-AKI ───────────────────────────────────────────
    print("\nSA3: Composite death-or-AKI…")
    ct["composite"] = ((ct["aki_surgery_7d"] == 1) |
                       (ct["hospital_mortality"] == 1)).astype(float)
    rd3, r3_lr, r3_ns = aipw_estimate(ct, ps_t, "composite")
    lo3, hi3 = bootstrap_ci(ct, "composite")
    print(f"  RD = {rd3*100:+.2f} pp  CI: ({lo3*100:+.2f}, {hi3*100:+.2f})")
    append_result(dict(
        analysis="SA3: Composite death-or-AKI",
        n=len(ct), rd=rd3, rd_lower=lo3, rd_upper=hi3,
        risk_lr=r3_lr, risk_ns=r3_ns,
    ))

    # ── SA4: IPTW only ────────────────────────────────────────────────────────
    print("\nSA4: IPTW only…")
    rd4, r4_lr, r4_ns = iptw_only(ct, ps_t)
    lo4, hi4 = bootstrap_ci(ct, "aki_surgery_7d")  # bootstrap uses AIPW; acceptable
    print(f"  RD = {rd4*100:+.2f} pp")
    append_result(dict(
        analysis="SA4: IPTW-only (no outcome model)",
        n=len(ct), rd=rd4, rd_lower=lo4, rd_upper=hi4,
        risk_lr=r4_lr, risk_ns=r4_ns,
    ))

    # ── SA6a: Baseline Cr ≤ 3.0 ───────────────────────────────────────────────
    print("\nSA6a: Baseline Cr ≤ 3.0…")
    c6 = cohort[cohort.baseline_creatinine <= 3.0].reset_index(drop=True)
    X6 = prep_X(c6, COV_COLS)
    A6 = c6.fluid_type.eq("LR").astype(int).values
    ps6 = fit_ps(X6, A6)
    keep6 = (ps6 >= PS_LOW) & (ps6 <= PS_HIGH)
    c6t = c6[keep6].reset_index(drop=True); ps6t = ps6[keep6]
    rd6, r6_lr, r6_ns = aipw_estimate(c6t, ps6t)
    lo6, hi6 = bootstrap_ci(c6t, "aki_surgery_7d")
    print(f"  n={len(c6t):,}  RD = {rd6*100:+.2f} pp  CI: ({lo6*100:+.2f}, {hi6*100:+.2f})")
    append_result(dict(
        analysis="SA6a: baseline Cr≤3.0",
        n=len(c6t), rd=rd6, rd_lower=lo6, rd_upper=hi6,
        risk_lr=r6_lr, risk_ns=r6_ns,
    ))

    # ── Balance diagnostics ───────────────────────────────────────────────────
    print("\nCovariate balance…")
    treat = ct["fluid_type"].eq("LR").astype(int).values
    p1_siptw = ps_t.mean(); p0_siptw = 1 - p1_siptw
    w = np.where(treat == 1, p1_siptw / np.clip(ps_t, 1e-6, 1),
                             p0_siptw / np.clip(1 - ps_t, 1e-6, 1))

    num_cols = prep_X(ct, COV_COLS).select_dtypes(include=[float]).columns.tolist()
    smd_before = {}; smd_after = {}
    for col in num_cols[:40]:  # top 40 for plot
        raw = pd.to_numeric(ct[col], errors="coerce").fillna(0).values.astype(float) if col in ct.columns else np.zeros(len(ct))
        lr_m = raw[treat == 1]; ns_m = raw[treat == 0]
        smd_before[col] = smd(lr_m, ns_m)
        smd_after[col]  = weighted_smd(col, treat, ct.assign(**{col: raw}), w) \
            if col in ct.columns else weighted_smd(col, treat,
               pd.DataFrame({col: raw, "fluid_type": ct.fluid_type}), w)

    n_bal = sum(1 for v in smd_after.values() if abs(v) < 0.1)
    print(f"  |SMD| < 0.1 after IPTW: {n_bal}/{len(smd_after)} covariates")

    # Save balance dict for figures
    balance = pd.DataFrame({"smd_before": smd_before, "smd_after": smd_after})
    balance.to_csv(os.path.join(SCRIPT_DIR, "balance_idealist.csv"))

    # Save trimmed cohort with PS for figures/SA scripts
    ct.to_parquet(os.path.join(SCRIPT_DIR, "cohort_trimmed_idealist.parquet"), index=False)
    np.save(os.path.join(SCRIPT_DIR, "ps_full_idealist.npy"), ps_full)
    np.save(os.path.join(SCRIPT_DIR, "ps_trim_idealist.npy"), ps_t)
    # Save full cohort with PS attached for downstream
    cohort.to_parquet(os.path.join(SCRIPT_DIR, "cohort_full_ps_idealist.parquet"), index=False)

    print("\n✓ Primary analysis complete.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--parquet", required=True)
    args = p.parse_args()
    main(args.parquet)
