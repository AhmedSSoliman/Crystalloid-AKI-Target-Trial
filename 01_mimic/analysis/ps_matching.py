"""
PS matching secondary analysis: 1:1 nearest-neighbor matching without replacement.

Caliper = 0.2 × SD(logit(PS)).
Outputs:
  figures/fig_matching_ps.png      - PS before/after matching
  figures/fig_matching_balance.png - SMD love plot pre/post match
  Appends result row to results_summary.csv

Usage:
  python3 ps_matching.py --parquet parquet/cohort_analytic.parquet
  python3 ps_matching.py --parquet parquet/cohort_enriched.parquet
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score

SEED = 42
np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

ENRICHED_COMORBIDITY_COLS = [
    "ckd", "liver_disease", "heart_failure", "sepsis",
    "cancer", "diabetes", "electrolyte_disorder",
    "vasopressor_pre", "gcs_min_pre", "hr_max_pre", "map_min_pre",
]


def savefig(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  saved → {path}")


def smd(col, treat, df):
    treat_arr = np.asarray(treat)
    vals = df[col].values if isinstance(df, pd.DataFrame) else np.asarray(df[col])
    g1 = vals[treat_arr == 1]
    g0 = vals[treat_arr == 0]
    g1 = g1[~np.isnan(g1.astype(float))]
    g0 = g0[~np.isnan(g0.astype(float))]
    pooled = np.sqrt((g1.var() + g0.var()) / 2)
    return abs((g1.mean() - g0.mean()) / pooled) if pooled > 0 else 0.0


def fit_ps(X, A):
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale",  StandardScaler()),
        ("clf",    LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
    ])
    pipe.fit(X, A)
    ps  = pipe.predict_proba(X)[:, 1]
    auc = roc_auc_score(A, ps)
    return ps, auc, pipe


def match_1to1(ps_lr, ps_ns, lr_idx, ns_idx, caliper):
    """
    Greedy 1:1 NN matching without replacement.
    ps_lr, ps_ns: logit-transformed PS for LR and NS units.
    Returns arrays of matched LR and NS original indices.
    """
    logit_lr = ps_lr.reshape(-1, 1)
    logit_ns = ps_ns.reshape(-1, 1)

    nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
    nn.fit(logit_ns)
    dists, indices = nn.kneighbors(logit_lr)

    matched_lr, matched_ns = [], []
    used_ns = set()

    # Sort LR units by ascending distance to their nearest NS (greedy)
    order = np.argsort(dists[:, 0])
    for i in order:
        d = dists[i, 0]
        j = indices[i, 0]
        if d <= caliper and j not in used_ns:
            matched_lr.append(lr_idx[i])
            matched_ns.append(ns_idx[j])
            used_ns.add(j)

    return np.array(matched_lr), np.array(matched_ns)


def prepare_cohort(parquet_path):
    df = pd.read_parquet(parquet_path)
    eligible = df.query("has_baseline_cr == 1 and has_followup_cr == 1").copy()
    eligible["aki_abs"]   = (eligible["max_fu_creatinine"] - eligible["baseline_creatinine"]) >= 0.3
    eligible["aki_ratio"] = (eligible["max_fu_creatinine"] / eligible["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    eligible["aki_flag"]  = (eligible["aki_abs"] | eligible["aki_ratio"]).astype(int)
    cohort = eligible.query("baseline_creatinine <= 4.0").copy()

    cohort["female"]          = (cohort["gender"] == "F").astype(int)
    cohort["race_white"]      = (cohort["race_cat"] == "WHITE").astype(int)
    cohort["race_black"]      = (cohort["race_cat"] == "BLACK").astype(int)
    cohort["log_baseline_cr"] = np.log1p(cohort["baseline_creatinine"])
    cohort["ckd_proxy"]       = (cohort["baseline_creatinine"] > 1.5).astype(int)
    cohort["treat"]           = (cohort["fluid_type"] == "LR").astype(int)
    cu_dummies = pd.get_dummies(cohort["first_careunit"], prefix="cu", drop_first=True).astype(int)
    cohort = pd.concat([cohort.reset_index(drop=True), cu_dummies.reset_index(drop=True)], axis=1)

    base_cov = (
        ["female", "anchor_age", "race_white", "race_black",
         "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"]
        + list(cu_dummies.columns)
    )
    extra_cov = [c for c in ENRICHED_COMORBIDITY_COLS if c in cohort.columns]
    if extra_cov:
        print(f"  Using enriched covariates: {extra_cov}")
    cov_cols = base_cov + extra_cov
    return cohort.reset_index(drop=True), cov_cols


def main(parquet_path, n_boot=200):
    print("=" * 60)
    print("PS Matching Secondary Analysis")
    print("=" * 60)

    cohort, cov_cols = prepare_cohort(parquet_path)
    X = cohort[cov_cols].values
    A = cohort["treat"].values
    Y = cohort["aki_flag"].values

    print(f"  Cohort: {len(cohort):,}  LR={A.sum():,}  NS={(1-A).sum():,}")

    # PS trimming (same as primary)
    ps_full, auc_full, _ = fit_ps(X, A)
    print(f"  PS AUC (full): {auc_full:.3f}")
    keep = (ps_full >= 0.05) & (ps_full <= 0.95)
    cohort_t = cohort[keep].reset_index(drop=True)
    X_t = X[keep]; A_t = A[keep]; Y_t = Y[keep]; ps_t = ps_full[keep]
    print(f"  After PS trimming: {keep.sum():,}")

    # Logit-transform PS for caliper matching
    logit_ps = np.log(ps_t / (1 - np.clip(ps_t, 1e-6, 1 - 1e-6)))
    caliper = 0.2 * logit_ps.std()
    print(f"  Caliper = 0.2 × SD(logit PS) = {caliper:.4f}")

    lr_mask = A_t == 1
    ns_mask = A_t == 0
    lr_idx  = np.where(lr_mask)[0]
    ns_idx  = np.where(ns_mask)[0]
    logit_lr = logit_ps[lr_idx]
    logit_ns = logit_ps[ns_idx]

    print("  Running 1:1 NN matching ...")
    matched_lr_idx, matched_ns_idx = match_1to1(logit_lr, logit_ns, lr_idx, ns_idx, caliper)
    n_pairs = len(matched_lr_idx)
    print(f"  Matched pairs: {n_pairs:,}  (from LR={lr_mask.sum():,}, NS={ns_mask.sum():,})")

    all_matched = np.concatenate([matched_lr_idx, matched_ns_idx])
    matched_df = cohort_t.iloc[all_matched].copy()
    A_m = A_t[all_matched]
    Y_m = Y_t[all_matched]
    ps_m = ps_t[all_matched]

    # Crude matched RD
    r_lr = Y_m[A_m == 1].mean()
    r_ns = Y_m[A_m == 0].mean()
    rd_match = r_lr - r_ns
    print(f"  Matched crude RD: {rd_match:.4f}  (LR={r_lr:.3f}, NS={r_ns:.3f})")

    # PS distribution after matching
    ps_dist = np.abs(logit_ps[matched_lr_idx] - logit_ps[matched_ns_idx])
    print(f"  PS distance (logit): mean={ps_dist.mean():.4f}, max={ps_dist.max():.4f}")

    # Bootstrap CI
    print(f"  Bootstrap CI ({n_boot} iterations) ...")
    def boot_match():
        idx = np.random.choice(n_pairs, n_pairs, replace=True)
        Y_lr_b = Y_t[matched_lr_idx[idx]]
        Y_ns_b = Y_t[matched_ns_idx[idx]]
        return (Y_lr_b.mean() - Y_ns_b.mean(),)

    estimates = Parallel(n_jobs=-1, prefer="threads")(
        delayed(boot_match)() for _ in range(n_boot)
    )
    rds = [e[0] for e in estimates]
    ci_lo = np.percentile(rds, 2.5)
    ci_hi = np.percentile(rds, 97.5)
    print(f"  RD = {rd_match:.4f}  95% CI ({ci_lo:.4f}, {ci_hi:.4f})")

    # ── Figures ───────────────────────────────────────────────────────────────
    # Fig 1: PS overlap before/after matching
    print("  Generating PS overlap figure ...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for fluid, color in [("LR", "#2196F3"), ("NS", "#FF9800")]:
        mask = cohort_t["fluid_type"] == fluid
        axes[0].hist(ps_t[mask.values], bins=60, alpha=0.5, label=fluid, color=color, density=True)
    axes[0].set_title("Before matching"); axes[0].set_xlabel("PS"); axes[0].legend()

    for fluid, color in [("LR", "#2196F3"), ("NS", "#FF9800")]:
        mask = matched_df["fluid_type"] == fluid
        axes[1].hist(ps_m[matched_df["fluid_type"].eq(fluid).values], bins=40,
                     alpha=0.5, label=fluid, color=color, density=True)
    axes[1].set_title("After 1:1 matching"); axes[1].set_xlabel("PS"); axes[1].legend()
    plt.suptitle("PS Distribution Before/After Matching", fontweight="bold")
    savefig("fig_matching_ps.png")

    # Fig 2: Love plot before/after matching
    print("  Generating balance love plot ...")
    num_cov = ["anchor_age", "female", "race_white", "race_black",
               "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"]
    num_cov = [c for c in num_cov if c in cohort_t.columns]
    smd_before = [smd(c, A_t, cohort_t) for c in num_cov]
    smd_after  = [smd(c, A_m, matched_df) for c in num_cov]

    fig, ax = plt.subplots(figsize=(8, max(5, len(num_cov) * 0.45 + 1)))
    y = np.arange(len(num_cov))
    ax.scatter(smd_before, y, color="#FF9800", s=60, zorder=3, label="Before matching")
    ax.scatter(smd_after,  y, color="#2196F3", s=60, zorder=3, label="After matching")
    ax.axvline(0.1, color="red", ls="--", lw=1.2, label="|SMD|=0.1")
    ax.axvline(0.0, color="grey", ls="-",  lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(num_cov, fontsize=9)
    ax.set_xlabel("|SMD|"); ax.set_title("Covariate Balance: Before/After PS Matching", fontweight="bold")
    ax.legend(fontsize=9); plt.tight_layout()
    savefig("fig_matching_balance.png")

    # ── Append to results_summary.csv ─────────────────────────────────────────
    csv_path = os.path.join(SCRIPT_DIR, "results_summary.csv")
    new_row = pd.DataFrame([{
        "analysis": "PS Matching (1:1 NN caliper=0.2SD)",
        "n": n_pairs * 2,
        "rd": rd_match,
        "rd_lower": ci_lo,
        "rd_upper": ci_hi,
        "risk_lr": r_lr,
        "risk_ns": r_ns,
    }])
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        # Remove old matching row if present
        existing = existing[~existing["analysis"].str.startswith("PS Matching")]
        combined = pd.concat([existing, new_row], ignore_index=True)
    else:
        combined = new_row
    combined.to_csv(csv_path, index=False)
    print(f"  Appended to {csv_path}")

    print("\nMatching summary:")
    print(f"  Pairs: {n_pairs:,}")
    print(f"  RD = {rd_match:.4f}  95% CI ({ci_lo:.4f}, {ci_hi:.4f})")
    print(f"  Risk LR={r_lr:.3f}  NS={r_ns:.3f}")
    pct_bal = sum(s < 0.1 for s in smd_after) / len(smd_after) * 100
    print(f"  Balance: {pct_bal:.0f}% of covariates |SMD|<0.1 after matching")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args()
    main(args.parquet, args.n_boot)
