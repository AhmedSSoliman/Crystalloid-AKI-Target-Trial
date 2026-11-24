"""
Full causal analysis: LR vs NS and 7-day incident AKI in MIMIC-IV ICU patients.

Methods: AIPW (augmented inverse probability weighting) with bootstrap 95% CI.
Data: analytic parquet built from MIMIC-IV (no raw CSV reads required).

Usage:
    python full_analysis.py --parquet /path/to/cohort_analytic.parquet
"""

import argparse
import os
import shutil
import tempfile
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score

# ── reproducibility ──────────────────────────────────────────────────────────
SEED = 42
np.random.seed(SEED)

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def savefig(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  saved → {path}")


def smd(col, treat, df):
    """Absolute standardized mean difference."""
    g1 = df.loc[treat == 1, col].dropna()
    g0 = df.loc[treat == 0, col].dropna()
    pooled_sd = np.sqrt((g1.var() + g0.var()) / 2)
    if pooled_sd == 0:
        return 0.0
    return abs((g1.mean() - g0.mean()) / pooled_sd)


def weighted_smd(col, treat, df, weights):
    """Weighted SMD using IPTW weights."""
    w1 = weights[treat == 1]
    w0 = weights[treat == 0]
    x1 = df.loc[treat == 1, col].fillna(df[col].median())
    x0 = df.loc[treat == 0, col].fillna(df[col].median())
    mu1 = np.average(x1, weights=w1)
    mu0 = np.average(x0, weights=w0)
    var1 = np.average((x1 - mu1) ** 2, weights=w1)
    var0 = np.average((x0 - mu0) ** 2, weights=w0)
    pooled = np.sqrt((var1 + var0) / 2)
    return abs((mu1 - mu0) / pooled) if pooled > 0 else 0.0


def fit_ps(X, A):
    """Fit logistic PS model, return propensity scores and AUC."""
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
    ])
    pipe.fit(X, A)
    ps = pipe.predict_proba(X)[:, 1]
    auc = roc_auc_score(A, ps)
    return ps, auc, pipe


def aipw_estimate(Y, A, ps, X_out, pipe_out=None):
    """
    AIPW risk difference estimator.
    If pipe_out is None, fits an outcome model on (A, X_out).
    Returns (rd, risk_lr, risk_ns).
    """
    n = len(Y)
    X_all = np.column_stack([A, X_out])

    if pipe_out is None:
        pipe_out = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
        ])
        pipe_out.fit(X_all, Y)

    X1 = np.column_stack([np.ones(n), X_out])
    X0 = np.column_stack([np.zeros(n), X_out])
    mu1 = pipe_out.predict_proba(X1)[:, 1]
    mu0 = pipe_out.predict_proba(X0)[:, 1]

    psi1 = mu1 + A * (Y - mu1) / np.clip(ps, 1e-6, 1 - 1e-6)
    psi0 = mu0 + (1 - A) * (Y - mu0) / np.clip(1 - ps, 1e-6, 1 - 1e-6)

    return psi1.mean() - psi0.mean(), psi1.mean(), psi0.mean()


def iptw_only_rd(Y, A, ps):
    """Stabilized IPTW risk difference (no outcome model)."""
    p1 = A.mean()
    p0 = 1 - p1
    w = np.where(A == 1, p1 / np.clip(ps, 1e-6, 1), p0 / np.clip(1 - ps, 1e-6, 1))
    r1 = np.sum(A * w * Y) / np.sum(A * w)
    r0 = np.sum((1 - A) * w * Y) / np.sum((1 - A) * w)
    return r1 - r0, r1, r0


def bootstrap_ci(fn, n_boot=200, alpha=0.05, n_jobs=-1):
    """Parallel bootstrap 95% CI for a scalar-returning function fn()."""
    estimates = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(fn)() for _ in range(n_boot)
    )
    rds = [e[0] for e in estimates]
    lo = np.percentile(rds, 100 * alpha / 2)
    hi = np.percentile(rds, 100 * (1 - alpha / 2))
    return lo, hi


# =============================================================================
# MAIN
# =============================================================================

def main(parquet_path, n_boot=200):
    print("=" * 65)
    print("LR vs NS → Incident AKI: Full Causal Analysis")
    print("=" * 65)

    # ── load data ─────────────────────────────────────────────────────────────
    print("\n[1/9] Loading data ...")
    # Copy to local /tmp first — /blue filesystem has ~50s latency for parquet reads
    _tmp_dir = tempfile.mkdtemp(prefix="causalai_")
    _local_parquet = os.path.join(_tmp_dir, "cohort.parquet")
    print(f"  Copying parquet to local tmp ({_tmp_dir}) ...")
    shutil.copy2(parquet_path, _local_parquet)
    df = pd.read_parquet(_local_parquet)
    print(f"  Raw cohort: {len(df):,} stays  (LR={df.fluid_type.eq('LR').sum():,}, NS={df.fluid_type.eq('NS').sum():,})")

    # ── Step 1: AKI outcome ───────────────────────────────────────────────────
    print("\n[2/9] Defining AKI outcome ...")
    eligible = df.query("has_baseline_cr == 1 and has_followup_cr == 1").copy()
    eligible["aki_abs"]   = (eligible["max_fu_creatinine"] - eligible["baseline_creatinine"]) >= 0.3
    eligible["aki_ratio"] = (eligible["max_fu_creatinine"] / eligible["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    eligible["aki_flag"]  = (eligible["aki_abs"] | eligible["aki_ratio"]).astype(int)

    for g, sub in eligible.groupby("fluid_type"):
        n, k = len(sub), sub["aki_flag"].sum()
        print(f"  {g}: AKI {k}/{n} ({100*k/n:.1f}%)")

    # ── Step 2: Cohort restriction ─────────────────────────────────────────────
    print("\n[3/9] Feasibility cohort restriction ...")
    n_all = len(df)
    n_cr  = len(eligible)
    # Exclude prevalent-AKI proxy: very high baseline Cr (>4 mg/dL)
    # This is a conservative proxy; true prevalent-AKI exclusion requires
    # time-stamped prior creatinine trajectories from raw labevents.
    cohort = eligible.query("baseline_creatinine <= 4.0").copy()
    n_primary = len(cohort)
    print(f"  All stays:              {n_all:>8,}")
    print(f"  Both creatinines avail: {n_cr:>8,}  (dropped {n_all-n_cr:,})")
    print(f"  Excl. baseline Cr>4:    {n_primary:>8,}  (dropped {n_cr-n_primary:,})")
    print(f"  Primary cohort:         {n_primary:>8,}")

    # ── Step 3: Covariate preparation ─────────────────────────────────────────
    print("\n[4/9] Preparing covariates ...")
    cohort["female"]          = (cohort["gender"] == "F").astype(int)
    cohort["race_white"]      = (cohort["race_cat"] == "WHITE").astype(int)
    cohort["race_black"]      = (cohort["race_cat"] == "BLACK").astype(int)
    cohort["log_baseline_cr"] = np.log1p(cohort["baseline_creatinine"])
    # Simple CKD proxy — severity scores (SOFA/APACHE) require re-extraction
    # from raw MIMIC-IV icustays/chartevents (not in this parquet).
    cohort["ckd_proxy"]       = (cohort["baseline_creatinine"] > 1.5).astype(int)
    cohort["treat"]           = (cohort["fluid_type"] == "LR").astype(int)

    # ICU type dummies (drop first careunit as reference)
    cu_dummies = pd.get_dummies(cohort["first_careunit"], prefix="cu", drop_first=True).astype(int)
    cohort = pd.concat([cohort.reset_index(drop=True), cu_dummies.reset_index(drop=True)], axis=1)

    # Enriched covariates — present when cohort_enriched.parquet is used
    ENRICHED_COLS = [
        "ckd", "liver_disease", "heart_failure", "sepsis", "cancer",
        "diabetes", "electrolyte_disorder", "vasopressor_pre",
        "gcs_min_pre", "hr_max_pre", "spo2_min_pre", "cr_slope_pre",
        "n_cr_pre", "year",
    ]
    extra_cov = [c for c in ENRICHED_COLS if c in cohort.columns]
    if extra_cov:
        print(f"  Using {len(extra_cov)} enriched covariates: {extra_cov}")
    cov_cols = (
        ["female", "anchor_age", "race_white", "race_black",
         "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"]
        + list(cu_dummies.columns)
        + extra_cov
    )

    X = cohort[cov_cols].values
    A = cohort["treat"].values
    Y = cohort["aki_flag"].values

    # ── Step 4: Propensity score ───────────────────────────────────────────────
    print("\n[5/9] Fitting propensity score model ...")
    ps_full, auc_full, pipe_ps = fit_ps(X, A)
    print(f"  PS model AUC (full cohort): {auc_full:.3f}")

    # Stabilized IPTW
    p_treat = A.mean()
    iptw = np.where(A == 1, p_treat / np.clip(ps_full, 1e-6, 1),
                    (1 - p_treat) / np.clip(1 - ps_full, 1e-6, 1))
    ess_lr = (iptw[A == 1].sum() ** 2) / (iptw[A == 1] ** 2).sum()
    ess_ns = (iptw[A == 0].sum() ** 2) / (iptw[A == 0] ** 2).sum()
    print(f"  ESS: LR={ess_lr:.0f}, NS={ess_ns:.0f}")

    # Trimming
    keep = (ps_full >= 0.05) & (ps_full <= 0.95)
    n_trimmed = (~keep).sum()
    print(f"  Trimmed (PS<0.05 or >0.95): {n_trimmed:,}  → {keep.sum():,} remain")
    ct = cohort[keep].copy()
    Xt = X[keep]; At = A[keep]; Yt = Y[keep]; ps_t = ps_full[keep]
    p_treat_t = At.mean()
    iptw_t = np.where(At == 1, p_treat_t / np.clip(ps_t, 1e-6, 1),
                      (1 - p_treat_t) / np.clip(1 - ps_t, 1e-6, 1))
    ps_t_new, auc_t, _ = fit_ps(Xt, At)
    print(f"  PS model AUC (trimmed): {auc_t:.3f}")

    # PS overlap figure
    print("\n  Generating PS overlap figure ...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for fluid, color, ax_idx in [("LR", "#2196F3", 0), ("NS", "#FF9800", 0)]:
        mask = cohort["fluid_type"] == fluid
        axes[0].hist(ps_full[mask.values], bins=60, alpha=0.5, label=fluid, color=color, density=True)
    axes[0].axvline(0.05, color="red", ls="--", lw=1.5, label="Trim 0.05/0.95")
    axes[0].axvline(0.95, color="red", ls="--", lw=1.5)
    axes[0].set_xlabel("Propensity Score"); axes[0].set_ylabel("Density")
    axes[0].set_title("PS Distribution (full cohort)"); axes[0].legend()

    for fluid, color in [("LR", "#2196F3"), ("NS", "#FF9800")]:
        mask = ct["fluid_type"] == fluid
        axes[1].hist(ps_t_new[mask.values], bins=60, alpha=0.5, label=fluid, color=color, density=True)
    axes[1].set_xlabel("Propensity Score"); axes[1].set_ylabel("Density")
    axes[1].set_title("PS Distribution (trimmed)"); axes[1].legend()
    plt.suptitle("Propensity Score Overlap by Fluid Group", fontweight="bold")
    savefig("fig_ps_overlap.png")

    # ── Step 5: Covariate balance / SMD ───────────────────────────────────────
    print("\n[6/9] Covariate balance (SMD) ...")
    smd_before, smd_after = {}, {}
    cov_df = pd.DataFrame(Xt, columns=cov_cols)
    treat_s = pd.Series(At)

    for col in cov_cols:
        smd_before[col] = smd(col, treat_s, cov_df)
        smd_after[col]  = weighted_smd(col, treat_s, cov_df, iptw_t)

    pct_below_01_before = np.mean([v < 0.1 for v in smd_before.values()])
    pct_below_01_after  = np.mean([v < 0.1 for v in smd_after.values()])
    print(f"  |SMD|<0.1 before IPTW: {100*pct_below_01_before:.0f}%")
    print(f"  |SMD|<0.1 after  IPTW: {100*pct_below_01_after:.0f}%")

    # Love plot
    fig, ax = plt.subplots(figsize=(8, max(5, len(cov_cols) * 0.4)))
    y = np.arange(len(cov_cols))
    ax.scatter([smd_before[c] for c in cov_cols], y, marker="o", color="#E53935", label="Before IPTW", zorder=3)
    ax.scatter([smd_after[c]  for c in cov_cols], y, marker="D", color="#1E88E5", label="After IPTW",  zorder=3)
    ax.axvline(0.1, color="gray", ls="--", lw=1.2, label="|SMD|=0.1")
    ax.set_yticks(y); ax.set_yticklabels(cov_cols, fontsize=8)
    ax.set_xlabel("Absolute SMD"); ax.set_title("Covariate Balance: Before and After IPTW")
    ax.legend(loc="lower right"); ax.set_xlim(left=0)
    savefig("fig_love_plot.png")

    # ── Step 6: AIPW primary estimate ─────────────────────────────────────────
    print(f"\n[7/9] AIPW primary estimation (bootstrap n={n_boot}) ...")

    rd_pt, r1_pt, r0_pt = aipw_estimate(Yt, At, ps_t_new, Xt)

    rng = np.random.default_rng(SEED)

    def boot_primary():
        idx = rng.integers(0, len(Yt), size=len(Yt))
        ps_b, _, _ = fit_ps(Xt[idx], At[idx])
        return aipw_estimate(Yt[idx], At[idx], ps_b, Xt[idx])

    if n_boot > 0:
        ci_lo, ci_hi = bootstrap_ci(boot_primary, n_boot=n_boot)
    else:
        ci_lo, ci_hi = float("nan"), float("nan")
    print(f"  Primary AIPW:  RD={rd_pt:.4f} (95% CI {ci_lo:.4f}, {ci_hi:.4f})")
    print(f"  LR risk: {r1_pt:.4f}   NS risk: {r0_pt:.4f}")

    results = [{
        "analysis": "Primary AIPW (PS-trimmed)",
        "n": int(len(Yt)),
        "rd": rd_pt, "rd_lower": ci_lo, "rd_upper": ci_hi,
        "risk_lr": r1_pt, "risk_ns": r0_pt,
    }]

    # ── Step 7: Sensitivity analyses ──────────────────────────────────────────
    print("\n[8/9] Sensitivity analyses ...")

    # SA1 — no trimming
    print("  SA1: No PS trimming ...")
    ps_all, _, _ = fit_ps(X, A)
    rd_sa1, r1_sa1, r0_sa1 = aipw_estimate(Y, A, ps_all, X)

    def boot_sa1():
        idx = rng.integers(0, len(Y), size=len(Y))
        ps_b, _, _ = fit_ps(X[idx], A[idx])
        return aipw_estimate(Y[idx], A[idx], ps_b, X[idx])

    lo1, hi1 = bootstrap_ci(boot_sa1, n_boot) if n_boot > 0 else (float("nan"), float("nan"))
    print(f"  SA1 RD={rd_sa1:.4f} (95% CI {lo1:.4f}, {hi1:.4f})")
    results.append({"analysis": "SA1: No PS trimming", "n": int(len(Y)),
                    "rd": rd_sa1, "rd_lower": lo1, "rd_upper": hi1,
                    "risk_lr": r1_sa1, "risk_ns": r0_sa1})

    # SA2 — stricter baseline Cr ≤ 2.0
    print("  SA2: Stricter baseline Cr≤2.0 ...")
    df_sa2  = eligible.query("baseline_creatinine <= 2.0").copy()
    df_sa2["female"]          = (df_sa2["gender"] == "F").astype(int)
    df_sa2["race_white"]      = (df_sa2["race_cat"] == "WHITE").astype(int)
    df_sa2["race_black"]      = (df_sa2["race_cat"] == "BLACK").astype(int)
    df_sa2["log_baseline_cr"] = np.log1p(df_sa2["baseline_creatinine"])
    df_sa2["ckd_proxy"]       = (df_sa2["baseline_creatinine"] > 1.5).astype(int)
    df_sa2["treat"]           = (df_sa2["fluid_type"] == "LR").astype(int)
    cu2 = pd.get_dummies(df_sa2["first_careunit"], prefix="cu", drop_first=True).astype(int)
    df_sa2 = pd.concat([df_sa2.reset_index(drop=True), cu2.reset_index(drop=True)], axis=1)
    cov2 = ["female", "anchor_age", "race_white", "race_black",
            "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"] + list(cu2.columns)
    X2 = df_sa2[cov2].values; A2 = df_sa2["treat"].values; Y2 = df_sa2["aki_flag"].values
    ps2, _, _ = fit_ps(X2, A2)
    keep2 = (ps2 >= 0.05) & (ps2 <= 0.95)
    X2t, A2t, Y2t, ps2t = X2[keep2], A2[keep2], Y2[keep2], ps2[keep2]
    ps2tr, _, _ = fit_ps(X2t, A2t)
    rd_sa2, r1_sa2, r0_sa2 = aipw_estimate(Y2t, A2t, ps2tr, X2t)

    def boot_sa2():
        idx = rng.integers(0, len(Y2t), size=len(Y2t))
        ps_b, _, _ = fit_ps(X2t[idx], A2t[idx])
        return aipw_estimate(Y2t[idx], A2t[idx], ps_b, X2t[idx])

    lo2, hi2 = bootstrap_ci(boot_sa2, n_boot) if n_boot > 0 else (float("nan"), float("nan"))
    print(f"  SA2 RD={rd_sa2:.4f} (95% CI {lo2:.4f}, {hi2:.4f})  n={len(Y2t):,}")
    results.append({"analysis": "SA2: Stricter baseline Cr≤2.0", "n": int(len(Y2t)),
                    "rd": rd_sa2, "rd_lower": lo2, "rd_upper": hi2,
                    "risk_lr": r1_sa2, "risk_ns": r0_sa2})

    # SA3 — composite outcome (death-or-AKI)
    print("  SA3: Composite death-or-AKI ...")
    ct2 = cohort[keep].copy()
    ct2["composite"] = ((ct2["hospital_expire_flag"] == 1) | (ct2["aki_flag"] == 1)).astype(int)
    Ycomp = ct2["composite"].values
    rd_sa3, r1_sa3, r0_sa3 = aipw_estimate(Ycomp, At, ps_t_new, Xt)

    def boot_sa3():
        idx = rng.integers(0, len(Ycomp), size=len(Ycomp))
        ps_b, _, _ = fit_ps(Xt[idx], At[idx])
        return aipw_estimate(Ycomp[idx], At[idx], ps_b, Xt[idx])

    lo3, hi3 = bootstrap_ci(boot_sa3, n_boot) if n_boot > 0 else (float("nan"), float("nan"))
    print(f"  SA3 RD={rd_sa3:.4f} (95% CI {lo3:.4f}, {hi3:.4f})")
    results.append({"analysis": "SA3: Composite death-or-AKI", "n": int(len(Ycomp)),
                    "rd": rd_sa3, "rd_lower": lo3, "rd_upper": hi3,
                    "risk_lr": r1_sa3, "risk_ns": r0_sa3})

    # SA4 — IPTW only (no outcome model)
    print("  SA4: IPTW-only (no outcome model) ...")
    rd_sa4, r1_sa4, r0_sa4 = iptw_only_rd(Yt, At, ps_t_new)

    def boot_sa4():
        idx = rng.integers(0, len(Yt), size=len(Yt))
        ps_b, _, _ = fit_ps(Xt[idx], At[idx])
        return iptw_only_rd(Yt[idx], At[idx], ps_b)

    lo4, hi4 = bootstrap_ci(boot_sa4, n_boot) if n_boot > 0 else (float("nan"), float("nan"))
    print(f"  SA4 RD={rd_sa4:.4f} (95% CI {lo4:.4f}, {hi4:.4f})")
    results.append({"analysis": "SA4: IPTW-only (no outcome model)", "n": int(len(Yt)),
                    "rd": rd_sa4, "rd_lower": lo4, "rd_upper": hi4,
                    "risk_lr": r1_sa4, "risk_ns": r0_sa4})

    # Save results CSV
    res_df = pd.DataFrame(results)
    csv_path = os.path.join(SCRIPT_DIR, "results_summary.csv")
    res_df.to_csv(csv_path, index=False, float_format="%.5f")
    print(f"\n  Results saved → {csv_path}")
    print(res_df.to_string(index=False))

    # Forest plot
    labels = [r["analysis"] for r in results]
    rds    = [r["rd"] for r in results]
    los    = [r["rd_lower"] for r in results]
    his    = [r["rd_upper"] for r in results]
    colors = ["#1E88E5", "#43A047", "#8E24AA", "#FB8C00", "#E53935"]

    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(labels))
    has_ci = not any(np.isnan(v) for v in los + his)
    for i in range(len(labels)):
        if has_ci:
            ax.plot([los[i], his[i]], [y[i], y[i]], color=colors[i], lw=2.5)
        ax.scatter(rds[i], y[i], color=colors[i], s=80, zorder=5)
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Risk Difference (LR − NS)")
    ci_note = "(95% bootstrap CI)" if has_ci else "(point estimates; run without --fast for CI)"
    ax.set_title(f"AIPW Risk Differences: Primary and Sensitivity Analyses\n{ci_note}", fontweight="bold")
    for i, r in enumerate(results):
        label = f'{r["rd"]:+.3f}'
        if has_ci:
            label += f' ({r["rd_lower"]:+.3f}, {r["rd_upper"]:+.3f})'
        ref_x = r["rd_upper"] if has_ci and not np.isnan(r["rd_upper"]) else r["rd"]
        ax.text(ref_x + 0.002, i, label, va="center", fontsize=7.5)
    xvals = [v for v in rds + (los if has_ci else []) + (his if has_ci else []) if not np.isnan(v)]
    ax.set_xlim(min(xvals) - 0.04, max(xvals) + 0.08)
    savefig("fig_aipw_results.png")

    # ── Step 8: Competing events ───────────────────────────────────────────────
    print("\n[9a/9] Competing events ...")
    ct3 = cohort[keep].copy()
    ct3["died_7d"] = (
        ct3["deathtime"].notna() &
        ((ct3["deathtime"] - ct3["fluid_starttime"]).dt.total_seconds() / 86400 <= 7)
    ).astype(int)
    ct3["discharged_alive_7d"] = (
        (ct3["hospital_expire_flag"] == 0) &
        ((ct3["dischtime"] - ct3["fluid_starttime"]).dt.total_seconds() / 86400 <= 7)
    ).astype(int)
    ct3["missing_fu_cr"] = (ct3["has_followup_cr"] == 0).astype(int)

    comp_rows = []
    for grp in ["LR", "NS"]:
        sub = ct3[ct3["fluid_type"] == grp]
        n_grp = len(sub)
        comp_rows.append({
            "Fluid": grp, "N": n_grp,
            "AKI (%)":           f'{100*sub["aki_flag"].mean():.1f}%',
            "Died ≤7d (%)":      f'{100*sub["died_7d"].mean():.1f}%',
            "Disch. alive ≤7d": f'{100*sub["discharged_alive_7d"].mean():.1f}%',
            "Missing fu Cr (%)": f'{100*sub["missing_fu_cr"].mean():.1f}%',
        })
    comp_df = pd.DataFrame(comp_rows)
    print(comp_df.to_string(index=False))

    # Competing events figure — stacked bar
    cats = ["AKI", "Died ≤7d", "Disch. alive ≤7d", "No event"]
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(2)
    bar_colors = ["#E53935", "#424242", "#43A047", "#90A4AE"]
    bottoms = np.zeros(2)
    for cat, color in zip(cats, bar_colors):
        if cat == "AKI":
            vals = [ct3[ct3.fluid_type == g]["aki_flag"].mean() for g in ["LR", "NS"]]
        elif cat == "Died ≤7d":
            vals = [ct3[ct3.fluid_type == g]["died_7d"].mean() for g in ["LR", "NS"]]
        elif cat == "Disch. alive ≤7d":
            vals = [ct3[ct3.fluid_type == g]["discharged_alive_7d"].mean() for g in ["LR", "NS"]]
        else:
            vals = [1 - ct3[ct3.fluid_type == g]["aki_flag"].mean()
                    - ct3[ct3.fluid_type == g]["died_7d"].mean()
                    - ct3[ct3.fluid_type == g]["discharged_alive_7d"].mean()
                    for g in ["LR", "NS"]]
        ax.bar(x, vals, bottom=bottoms, color=color, label=cat, edgecolor="white")
        bottoms += np.array(vals)
    ax.set_xticks(x); ax.set_xticklabels(["LR", "NS"])
    ax.set_ylabel("Proportion"); ax.set_ylim(0, 1.05)
    ax.set_title("7-Day Outcomes by Fluid Group (PS-trimmed cohort)", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    savefig("fig_competing_events.png")

    # ── Cohort flow figure ────────────────────────────────────────────────────
    print("\n[9b/9] Cohort flow diagram ...")
    fig, ax = plt.subplots(figsize=(7, 9))
    ax.axis("off")
    boxes = [
        (0.5, 0.92, f"All MIMIC-IV ICU stays with LR/NS ≤24h\n(from parquet)  n={len(df):,}", "#E3F2FD"),
        (0.5, 0.74, f"Has baseline AND follow-up creatinine\nn={n_cr:,}", "#E8F5E9"),
        (0.5, 0.56, f"Baseline Cr ≤4.0 mg/dL\n(excl. prevalent-AKI proxy)\nn={n_primary:,}", "#FFF8E1"),
        (0.5, 0.38, f"PS-trimmed (0.05–0.95)\nn={int(keep.sum()):,}", "#F3E5F5"),
        (0.3, 0.17, f"LR\nn={int(At.sum()):,}", "#BBDEFB"),
        (0.7, 0.17, f"NS\nn={int((1-At).sum()):,}", "#FFE0B2"),
    ]
    for (bx, by, text, color) in boxes:
        ax.text(bx, by, text, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor=color, edgecolor="gray"),
                transform=ax.transAxes)
    for i in range(len(boxes) - 3):
        ax.annotate("", xy=(0.5, boxes[i+1][1] + 0.06), xytext=(0.5, boxes[i][1] - 0.06),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="black"))
    ax.annotate("", xy=(0.3, 0.26), xytext=(0.5, boxes[-3][1] - 0.06),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="black"))
    ax.annotate("", xy=(0.7, 0.26), xytext=(0.5, boxes[-3][1] - 0.06),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="black"))
    ax.set_title("Cohort Flow Diagram", fontweight="bold", fontsize=12, y=0.98)
    savefig("fig_cohort_flow.png")

    # ── Final summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("RESULTS SUMMARY")
    print("=" * 65)
    print(res_df[["analysis", "n", "rd", "rd_lower", "rd_upper", "risk_lr", "risk_ns"]]
          .round(4).to_string(index=False))
    print(f"\nFigures saved to: {FIG_DIR}")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LR vs NS AKI causal analysis (AIPW)")
    parser.add_argument("--parquet", required=True, help="Path to cohort_analytic.parquet")
    parser.add_argument("--n-boot", type=int, default=200, help="Bootstrap iterations (default 200)")
    parser.add_argument("--fast", action="store_true", help="Skip bootstrap (point estimates only, for testing)")
    args = parser.parse_args()
    main(args.parquet, n_boot=0 if args.fast else args.n_boot)
