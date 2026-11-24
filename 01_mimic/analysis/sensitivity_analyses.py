"""
Additional sensitivity analyses (SA5–SA8) and causal DAG figure.

SA5: Shorter observation window (3d) — skipped if cr_48h not available
SA6a: Prevalent AKI exclusion with baseline Cr threshold > 3.0
SA7: Sex subgroup analysis (female vs male)
SA8: ICU-type subgroup analysis (top 3 ICU types)

Also generates: figures/fig_causal_dag.png

Usage:
  python3 sensitivity_analyses.py --parquet parquet/cohort_analytic.parquet
  python3 sensitivity_analyses.py --parquet parquet/cohort_enriched.parquet
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
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
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


def fit_ps(X, A):
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale",  StandardScaler()),
        ("clf",    LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
    ])
    pipe.fit(X, A)
    ps = pipe.predict_proba(X)[:, 1]
    return ps, roc_auc_score(A, ps), pipe


def aipw_estimate(Y, A, ps, X_out):
    n = len(Y)
    X_all = np.column_stack([A, X_out])
    pipe_out = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale",  StandardScaler()),
        ("clf",    LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
    ])
    pipe_out.fit(X_all, Y)
    X1 = np.column_stack([np.ones(n), X_out])
    X0 = np.column_stack([np.zeros(n), X_out])
    mu1 = pipe_out.predict_proba(X1)[:, 1]
    mu0 = pipe_out.predict_proba(X0)[:, 1]
    psi1 = mu1 + A * (Y - mu1) / np.clip(ps, 1e-6, 1 - 1e-6)
    psi0 = mu0 + (1 - A) * (Y - mu0) / np.clip(1 - ps, 1e-6, 1 - 1e-6)
    return psi1.mean() - psi0.mean(), psi1.mean(), psi0.mean()


def bootstrap_ci(fn, n_boot=200):
    estimates = Parallel(n_jobs=-1, prefer="threads")(
        delayed(fn)() for _ in range(n_boot)
    )
    rds = [e[0] for e in estimates]
    return np.percentile(rds, 2.5), np.percentile(rds, 97.5)


def prepare_base(df, cr_threshold=4.0, outcome="aki_flag"):
    eligible = df.query("has_baseline_cr == 1 and has_followup_cr == 1").copy()
    eligible["aki_abs"]   = (eligible["max_fu_creatinine"] - eligible["baseline_creatinine"]) >= 0.3
    eligible["aki_ratio"] = (eligible["max_fu_creatinine"] / eligible["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    eligible["aki_flag"]  = (eligible["aki_abs"] | eligible["aki_ratio"]).astype(int)
    eligible["death_7d"]  = (
        (eligible["deathtime"].notna()) &
        ((pd.to_datetime(eligible["deathtime"]) - pd.to_datetime(eligible["fluid_starttime"])).dt.total_seconds() <= 7 * 86400)
    ).astype(int)
    eligible["death_or_aki"] = ((eligible["aki_flag"] == 1) | (eligible["death_7d"] == 1)).astype(int)

    cohort = eligible.query(f"baseline_creatinine <= {cr_threshold}").copy()

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
    cov_cols = base_cov + extra_cov
    return cohort.reset_index(drop=True), cov_cols


def run_aipw_trimmed(cohort, cov_cols, outcome_col="aki_flag", n_boot=200, label=""):
    """Run AIPW with PS trimming and bootstrap CI. Returns result dict."""
    X = cohort[cov_cols].values
    A = cohort["treat"].values
    Y = cohort[outcome_col].values

    ps, auc, _ = fit_ps(X, A)
    keep = (ps >= 0.05) & (ps <= 0.95)
    Xt, At, Yt, pst = X[keep], A[keep], Y[keep], ps[keep]
    n = keep.sum()

    rd, r1, r0 = aipw_estimate(Yt, At, pst, Xt)

    def boot():
        idx = np.random.choice(n, n, replace=True)
        ps_b, _, _ = fit_ps(Xt[idx], At[idx])
        keep_b = (ps_b >= 0.05) & (ps_b <= 0.95)
        if keep_b.sum() < 20:
            return (np.nan, np.nan, np.nan)
        return aipw_estimate(Yt[idx][keep_b], At[idx][keep_b], ps_b[keep_b], Xt[idx][keep_b])

    lo, hi = bootstrap_ci(boot, n_boot)
    if label:
        print(f"  {label}: n={n:,}  RD={rd:.4f}  95%CI ({lo:.4f},{hi:.4f})  LR={r1:.3f} NS={r0:.3f}")
    return {"n": n, "rd": rd, "rd_lower": lo, "rd_upper": hi, "risk_lr": r1, "risk_ns": r0}


def sa5_short_window(cohort, cov_cols):
    """SA5: 3-day observation window — requires cr_48h from enriched parquet."""
    if "cr_48h" not in cohort.columns:
        print("  SA5 (3d window): SKIPPED — cr_48h not available; re-run after extract_confounders.py")
        return None
    # Use cr_48h as the follow-up creatinine (approximates 3d better than 7d max)
    c = cohort.copy()
    c["aki_abs_48h"]   = (c["cr_48h"] - c["baseline_creatinine"]) >= 0.3
    c["aki_ratio_48h"] = (c["cr_48h"] / c["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    c["aki_48h"]       = (c["aki_abs_48h"] | c["aki_ratio_48h"]).fillna(False).astype(int)
    c = c[c["cr_48h"].notna()]
    print(f"  SA5 (48h window): n={len(c):,} with cr_48h")
    return run_aipw_trimmed(c, cov_cols, outcome_col="aki_48h", label="SA5 48h window")


def sa6a_stricter_baseline(df, cov_cols_base):
    """SA6a: Exclude baseline Cr > 3.0 (instead of > 4.0)."""
    cohort, cov_cols = prepare_base(df, cr_threshold=3.0)
    print(f"  SA6a (baseline Cr≤3.0): n={len(cohort):,}")
    return run_aipw_trimmed(cohort, cov_cols, label="SA6a baseline Cr≤3.0")


def sa7_sex_subgroups(cohort, cov_cols, n_boot=200):
    """SA7: Female vs male subgroup."""
    results = []
    for sex_val, label in [(1, "Female"), (0, "Male")]:
        sub = cohort[cohort["female"] == sex_val].copy()
        if len(sub) < 200:
            print(f"  SA7 {label}: too few ({len(sub)}), skipping")
            continue
        res = run_aipw_trimmed(sub, cov_cols, n_boot=n_boot, label=f"SA7 {label}")
        res["analysis"] = f"SA7: {label} subgroup"
        results.append(res)
    return results


def sa8_icu_subgroups(cohort, cov_cols, n_boot=200):
    """SA8: Top 3 ICU-type subgroups."""
    top3 = cohort["first_careunit"].value_counts().head(3).index.tolist()
    results = []
    for icu in top3:
        sub = cohort[cohort["first_careunit"] == icu].copy()
        # Rebuild ICU dummies within subgroup
        sub2, cov2 = prepare_base_from_cohort(sub)
        if len(sub2) < 100:
            continue
        res = run_aipw_trimmed(sub2, cov2, n_boot=n_boot, label=f"SA8 {icu}")
        res["analysis"] = f"SA8: {icu}"
        results.append(res)
    return results


def prepare_base_from_cohort(sub):
    """Re-derive covariates from an already-restricted cohort subset."""
    cohort = sub.copy()
    cu_dummies = pd.get_dummies(cohort["first_careunit"], prefix="cu", drop_first=True).astype(int)
    cohort = pd.concat([cohort.reset_index(drop=True), cu_dummies.reset_index(drop=True)], axis=1)
    base_cov = (
        ["female", "anchor_age", "race_white", "race_black",
         "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"]
        + list(cu_dummies.columns)
    )
    extra = [c for c in ENRICHED_COMORBIDITY_COLS if c in cohort.columns]
    return cohort.reset_index(drop=True), base_cov + extra


# ── Causal DAG ────────────────────────────────────────────────────────────────
def draw_causal_dag():
    """Draw causal DAG using matplotlib patches and arrows."""
    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14); ax.set_ylim(0, 9)
    ax.axis("off")
    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("#FAFAFA")

    # ── node definitions: (x, y, label, shape, color) ────────────────────────
    confounders = [
        (1.5, 8.0, "Age / Sex / Race"),
        (1.5, 6.8, "ICU Type"),
        (1.5, 5.6, "Baseline Creatinine"),
        (1.5, 4.4, "CKD / Comorbidities"),
        (1.5, 3.2, "Illness Severity\n(SOFA / hemodynamics)"),
        (1.5, 2.0, "Pre-ICU Fluids"),
        (1.5, 0.8, "Calendar Period"),
    ]
    exposure  = (5.5, 4.5, "Initial ICU Fluid\n(LR vs NS)")
    outcome   = (10.5, 4.5, "7-day Incident AKI")
    mediator  = (7.5, 2.0, "Post-treatment\nFluid Management\n(not adjusted)")

    def draw_ellipse(ax, cx, cy, label, fc, ec, fs=9):
        ellipse = mpatches.Ellipse((cx, cy), 2.6, 0.95, fc=fc, ec=ec, lw=1.5, zorder=2)
        ax.add_patch(ellipse)
        ax.text(cx, cy, label, ha="center", va="center", fontsize=fs,
                fontweight="normal", zorder=3)

    def draw_rect(ax, cx, cy, label, fc, ec, fs=9.5, lw=2.0, ls="-", width=2.8, height=1.0):
        rect = mpatches.FancyBboxPatch(
            (cx - width / 2, cy - height / 2), width, height,
            boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=lw, ls=ls, zorder=2
        )
        ax.add_patch(rect)
        ax.text(cx, cy, label, ha="center", va="center", fontsize=fs,
                fontweight="bold", zorder=3)

    def arrow(ax, x0, y0, x1, y1, color="#444", lw=1.3, ls="-", style="->"):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle=style, color=color,
                                   lw=lw, linestyle=ls,
                                   connectionstyle="arc3,rad=0.0"),
                    zorder=1)

    # Draw confounders
    for cx, cy, label in confounders:
        draw_ellipse(ax, cx, cy, label, fc="#E8E8E8", ec="#888888", fs=8.5)

    # Draw exposure
    draw_rect(ax, *exposure[:2], exposure[2], fc="#BBDEFB", ec="#1565C0", fs=9.5)

    # Draw outcome
    draw_rect(ax, *outcome[:2], outcome[2], fc="#FFCDD2", ec="#B71C1C", fs=9.5)

    # Draw mediator (dashed)
    draw_rect(ax, *mediator[:2], mediator[2], fc="#FFE0B2", ec="#E65100",
              fs=8, lw=1.5, ls="--", width=2.8, height=1.1)
    ax.text(mediator[0], mediator[1] - 0.75, "(not adjusted)", ha="center",
            va="center", fontsize=7.5, color="#E65100", style="italic")

    # Arrows: confounder → exposure
    for cx, cy, _ in confounders:
        arrow(ax, cx + 1.3, cy, exposure[0] - 1.4, exposure[1],
              color="#888888", lw=1.0)

    # Arrows: confounder → outcome
    for cx, cy, _ in confounders:
        arrow(ax, cx + 1.3, cy, outcome[0] - 1.4, outcome[1],
              color="#888888", lw=1.0)

    # Exposure → outcome (main causal arrow)
    arrow(ax, exposure[0] + 1.4, exposure[1],
          outcome[0] - 1.4, outcome[1],
          color="#1565C0", lw=2.5, style="-|>")

    # Exposure → mediator (dashed)
    arrow(ax, exposure[0] + 0.3, exposure[1] - 0.5,
          mediator[0] - 1.0, mediator[1] + 0.3,
          color="#E65100", lw=1.5, ls="dashed", style="->")

    # Mediator → outcome (dashed)
    arrow(ax, mediator[0] + 1.4, mediator[1] + 0.3,
          outcome[0] - 0.1, outcome[1] - 0.5,
          color="#E65100", lw=1.5, ls="dashed", style="->")

    # Legend
    legend_items = [
        mpatches.Patch(fc="#E8E8E8", ec="#888888", label="Confounder (adjusted)"),
        mpatches.Patch(fc="#BBDEFB", ec="#1565C0", label="Exposure"),
        mpatches.Patch(fc="#FFCDD2", ec="#B71C1C", label="Outcome"),
        mpatches.Patch(fc="#FFE0B2", ec="#E65100", label="Mediator (not adjusted)"),
    ]
    ax.legend(handles=legend_items, loc="lower right", fontsize=8.5,
              framealpha=0.9, edgecolor="#aaa")

    ax.set_title("Proposed Causal DAG: Initial ICU Fluid and 7-day AKI",
                 fontsize=12, fontweight="bold", pad=10)
    plt.tight_layout()
    savefig("fig_causal_dag.png")


# ── Subgroup forest plot ──────────────────────────────────────────────────────
def plot_subgroup_forest(results, csv_path):
    """Forest plot for SA7 (sex) and SA8 (ICU type) subgroups."""
    if not results:
        print("  No subgroup results to plot.")
        return
    labels  = [r["analysis"] for r in results]
    rds     = [r["rd"]       for r in results]
    lows    = [r["rd_lower"] for r in results]
    highs   = [r["rd_upper"] for r in results]
    ns      = [r["n"]        for r in results]

    fig, ax = plt.subplots(figsize=(9, max(4, len(results) * 0.7 + 1.5)))
    y = np.arange(len(results))
    colors = ["#2196F3" if "Female" in l or "Male" in l else "#FF7043" for l in labels]

    for i in range(len(results)):
        lo, hi = lows[i], highs[i]
        ax.plot([lo, hi], [y[i], y[i]], color=colors[i], lw=2, zorder=2)
        ax.scatter([rds[i]], [y[i]], color=colors[i], s=70, zorder=3)

    ax.axvline(0, color="black", lw=1.0, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}  (n={n:,})" for l, n in zip(labels, ns)], fontsize=9)
    ax.set_xlabel("Risk Difference (LR − NS)", fontsize=10)
    ax.set_title("Subgroup Analyses (SA7: Sex, SA8: ICU Type)\nAIPW Risk Difference with 95% Bootstrap CI",
                 fontweight="bold")

    # Annotate RD values
    for i, (rd, lo, hi) in enumerate(zip(rds, lows, highs)):
        ax.text(max(highs) + 0.002, y[i],
                f"{rd:.3f} ({lo:.3f}, {hi:.3f})",
                va="center", ha="left", fontsize=8)

    ax.set_xlim(min(lows) - 0.01, max(highs) + 0.04)
    plt.tight_layout()
    savefig("fig_subgroup_forest.png")

    # Append rows to CSV
    new_rows = pd.DataFrame(results)
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        # Remove old SA7/SA8 rows
        existing = existing[~existing["analysis"].str.startswith("SA7") &
                            ~existing["analysis"].str.startswith("SA8")]
        combined = pd.concat([existing, new_rows], ignore_index=True)
    else:
        combined = new_rows
    combined.to_csv(csv_path, index=False)
    print(f"  Appended {len(results)} subgroup rows to {csv_path}")


# ── main ──────────────────────────────────────────────────────────────────────
def main(parquet_path, n_boot=200):
    print("=" * 60)
    print("Additional Sensitivity Analyses (SA5–SA8) + Causal DAG")
    print("=" * 60)

    df = pd.read_parquet(parquet_path)
    cohort, cov_cols = prepare_base(df)
    print(f"  Primary cohort: {len(cohort):,} rows")

    csv_path = os.path.join(SCRIPT_DIR, "results_summary.csv")
    all_new_rows = []

    # SA5: short window
    print("\n[SA5] 3-day / 48h observation window")
    res5 = sa5_short_window(cohort, cov_cols)
    if res5:
        res5["analysis"] = "SA5: 48h creatinine window"
        all_new_rows.append(res5)

    # SA6a: stricter baseline Cr threshold
    print("\n[SA6a] Stricter prevalent-AKI exclusion (baseline Cr ≤ 3.0)")
    res6a = sa6a_stricter_baseline(df, cov_cols)
    res6a["analysis"] = "SA6a: baseline Cr≤3.0"
    all_new_rows.append(res6a)

    # SA7: sex subgroups
    print("\n[SA7] Sex subgroups")
    res7 = sa7_sex_subgroups(cohort, cov_cols, n_boot=n_boot)
    all_new_rows.extend(res7)

    # SA8: ICU-type subgroups
    print("\n[SA8] ICU-type subgroups")
    res8 = sa8_icu_subgroups(cohort, cov_cols, n_boot=n_boot)
    all_new_rows.extend(res8)

    # Subgroup forest plot (SA7 + SA8)
    print("\n  Subgroup forest plot ...")
    subgroup_rows = [r for r in all_new_rows if r.get("analysis", "").startswith("SA7")
                     or r.get("analysis", "").startswith("SA8")]
    plot_subgroup_forest(subgroup_rows, csv_path)

    # Append SA5 and SA6a to CSV separately
    non_subgroup = [r for r in all_new_rows if not (
        r.get("analysis", "").startswith("SA7") or r.get("analysis", "").startswith("SA8")
    )]
    if non_subgroup:
        new_df = pd.DataFrame(non_subgroup)
        if os.path.exists(csv_path):
            existing = pd.read_csv(csv_path)
            for lab in [r["analysis"] for r in non_subgroup]:
                existing = existing[existing["analysis"] != lab]
            combined = pd.concat([existing, new_df], ignore_index=True)
        else:
            combined = new_df
        combined.to_csv(csv_path, index=False)
        print(f"  Wrote {len(non_subgroup)} non-subgroup SA rows to {csv_path}")

    # Causal DAG
    print("\n  Drawing causal DAG ...")
    draw_causal_dag()

    print("\n" + "=" * 60)
    print("Done. New figures:")
    for f in ["fig_subgroup_forest.png", "fig_causal_dag.png"]:
        fp = os.path.join(FIG_DIR, f)
        if os.path.exists(fp):
            print(f"  ✓ {fp}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args()
    main(args.parquet, args.n_boot)
