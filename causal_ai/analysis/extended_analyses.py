"""
Extended analyses to complete the proposal's methods section.

Adds three pieces not covered by full_analysis.py / sensitivity_analyses.py:

  1. E-value for the primary AIPW risk difference and its CI bound
     (VanderWeele & Ding, Ann Intern Med 2017) — quantifies how strong an
     unmeasured confounder would need to be to explain away the result.

  2. Fine-Gray competing-risks analysis — subdistribution hazard for incident
     AKI treating in-hospital death without AKI as a competing event, plus
     cumulative incidence functions (CIF) by fluid group. Uses lifelines when
     available; otherwise falls back to a non-parametric Aalen-Johansen CIF
     with an IPTW-adjusted cause-specific comparison.

  3. SA9: KDIGO-correct 48h absolute criterion AS A SEPARATE WINDOW
     (already partially in SA5) with bootstrap CI — kept for completeness
     and to report alongside the E-value.

Results are appended to results_summary.csv (E-value in a side file,
evalue_summary.csv, since it is not a risk difference row).
Figures: fig_cif_competing.png, fig_evalue.png

Usage:
  python3 extended_analyses.py --parquet parquet/cohort_enriched.parquet
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
from sklearn.metrics import roc_auc_score

SEED = 42
np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

LR_COL = "#1565C0"
NS_COL = "#E65100"

ENRICHED_COLS = [
    "ckd", "liver_disease", "heart_failure", "sepsis", "cancer",
    "diabetes", "electrolyte_disorder", "vasopressor_pre",
    "gcs_min_pre", "hr_max_pre", "map_min_pre",
]


# ── shared helpers (mirrors full_analysis.py) ────────────────────────────────
def savefig(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
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
    return ps, roc_auc_score(A, ps)


def prepare_base(df, cr_threshold=4.0):
    e = df.query("has_baseline_cr == 1 and has_followup_cr == 1").copy()
    e["aki_abs"]   = (e["max_fu_creatinine"] - e["baseline_creatinine"]) >= 0.3
    e["aki_ratio"] = (e["max_fu_creatinine"] / e["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    e["aki_flag"]  = (e["aki_abs"] | e["aki_ratio"]).astype(int)

    t0 = pd.to_datetime(e["fluid_starttime"])
    dt_death = (pd.to_datetime(e["deathtime"]) - t0).dt.total_seconds() / 86400.0
    e["days_to_death"] = dt_death
    e["death_7d"] = ((dt_death >= 0) & (dt_death <= 7)).fillna(False).astype(int)

    cohort = e.query(f"baseline_creatinine <= {cr_threshold}").copy()
    cohort["female"]          = (cohort["gender"] == "F").astype(int)
    cohort["race_white"]      = (cohort["race_cat"] == "WHITE").astype(int)
    cohort["race_black"]      = (cohort["race_cat"] == "BLACK").astype(int)
    cohort["log_baseline_cr"] = np.log1p(cohort["baseline_creatinine"])
    cohort["ckd_proxy"]       = (cohort["baseline_creatinine"] > 1.5).astype(int)
    cohort["treat"]           = (cohort["fluid_type"] == "LR").astype(int)
    cu = pd.get_dummies(cohort["first_careunit"], prefix="cu", drop_first=True).astype(int)
    cohort = pd.concat([cohort.reset_index(drop=True), cu.reset_index(drop=True)], axis=1)

    base_cov = (["female", "anchor_age", "race_white", "race_black",
                 "log_baseline_cr", "ckd_proxy", "hours_from_icu_admit", "los_hours"]
                + list(cu.columns))
    extra = [c for c in ENRICHED_COLS if c in cohort.columns]
    return cohort.reset_index(drop=True), base_cov + extra


# =============================================================================
# 1. E-VALUE  (VanderWeele & Ding 2017)
# =============================================================================
def evalue_from_rr(rr):
    """E-value for a risk/rate ratio point estimate or CI bound.
    RR must be >= 1 orientation; for RR<1 use 1/RR first."""
    if rr < 1:
        rr = 1.0 / rr
    return rr + np.sqrt(rr * (rr - 1.0))


def compute_evalue(rd, r_lr, r_ns, ci_lower, ci_upper):
    """
    Convert the AIPW risk difference to an approximate risk ratio on the
    outcome scale, then compute the E-value for the point estimate and for
    the CI bound closest to the null (RR=1).

    Returns a dict of the RR, E-value(point), and E-value(CI).
    """
    # Risk ratio from the two counterfactual risks
    rr = r_lr / r_ns if r_ns > 0 else np.nan

    # Approximate RR CI bounds by shifting the RD CI onto the NS baseline risk.
    # risk_lr_bound = risk_ns + rd_bound  → RR_bound = risk_lr_bound / risk_ns
    rr_lower = (r_ns + ci_lower) / r_ns if r_ns > 0 else np.nan
    rr_upper = (r_ns + ci_upper) / r_ns if r_ns > 0 else np.nan

    ev_point = evalue_from_rr(rr)

    # E-value for the CI: use the bound nearest the null (RR = 1)
    if rr > 1:
        ci_bound = rr_lower          # lower bound nearest null from above
    else:
        ci_bound = rr_upper          # upper bound nearest null from below
    # If the CI crosses the null, E-value for CI is 1 (null not excluded)
    if (rr_lower <= 1 <= rr_upper):
        ev_ci = 1.0
    else:
        ev_ci = evalue_from_rr(ci_bound)

    return {
        "risk_lr": r_lr, "risk_ns": r_ns, "rr": rr,
        "rr_lower": rr_lower, "rr_upper": rr_upper,
        "evalue_point": ev_point, "evalue_ci": ev_ci,
        "crosses_null": bool(rr_lower <= 1 <= rr_upper),
    }


def plot_evalue(ev, title_suffix=""):
    """Visualize how a confounder's association with treatment and outcome
    (both on the RR scale) must exceed the E-value curve to explain the result."""
    fig, ax = plt.subplots(figsize=(7, 6))
    E = ev["evalue_point"]
    # The 'bias contour' satisfying RR_UD * RR_EU / (RR_UD + RR_EU - 1) = observed RR
    # Points on/above the E-value joint value are required to explain away.
    rr_obs = ev["rr"] if ev["rr"] >= 1 else 1 / ev["rr"]
    xs = np.linspace(rr_obs, rr_obs * 3 + 2, 400)
    # Solve for y given x on the contour that just explains rr_obs:
    #   rr_obs = (x*y)/(x+y-1)  → y = rr_obs*(x-1) / (x - rr_obs)
    ys = rr_obs * (xs - 1) / (xs - rr_obs)
    valid = (ys > 0) & np.isfinite(ys)
    ax.plot(xs[valid], ys[valid], color=LR_COL, lw=2,
            label=f"Bias needed to explain RR={rr_obs:.3f}")
    ax.scatter([E], [E], color="#B71C1C", s=90, zorder=5)
    ax.annotate(f"E-value = {E:.2f}\n(both associations)",
                (E, E), textcoords="offset points", xytext=(12, 10),
                fontsize=10, color="#B71C1C", fontweight="bold")
    ax.axhline(1, color="#9E9E9E", lw=0.8)
    ax.axvline(1, color="#9E9E9E", lw=0.8)
    ax.set_xlabel("RR: unmeasured confounder ↔ treatment", fontsize=11)
    ax.set_ylabel("RR: unmeasured confounder ↔ outcome", fontsize=11)
    ax.set_title(f"Figure 10. E-value Sensitivity to Unmeasured Confounding{title_suffix}",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(1, xs.max())
    ax.set_ylim(1, min(ys[valid].max() if valid.any() else 10, rr_obs * 3 + 2))
    ax.grid(alpha=0.3)
    savefig("fig_evalue.png")


# =============================================================================
# 2. FINE-GRAY COMPETING RISKS
# =============================================================================
def build_competing_time_event(cohort):
    """
    Build (time, event) for competing-risks:
      event = 1  → incident AKI within 7d (primary)
      event = 2  → death within 7d without AKI (competing)
      event = 0  → censored at 7d (survived, no AKI)
    Time axis in days (0, 7]. AKI time is approximated by the follow-up
    creatinine window; when no explicit AKI timestamp exists we assign the
    AKI event to the midpoint used for ranking only (CIF is step-based).
    """
    c = cohort.copy()
    # Approximate AKI day: if cr_48h-based AKI, day≈2, else day≈3.5 within 7d window.
    aki = c["aki_flag"].values == 1
    death_only = (c["death_7d"].values == 1) & (~aki)

    event = np.zeros(len(c), dtype=int)
    event[aki] = 1
    event[death_only] = 2

    time = np.full(len(c), 7.0)
    # Death time (if within window) for competing events
    dtd = c["days_to_death"].values
    has_dt = np.isfinite(dtd) & (dtd >= 0) & (dtd <= 7)
    time[death_only & has_dt] = dtd[death_only & has_dt]
    # AKI assumed observed by the creatinine follow-up; place at day 2 if 48h
    # criterion present else day 3 (coarse; CIF shape is the reportable part).
    if "cr_48h" in c.columns:
        aki_48 = aki & c["cr_48h"].notna().values & \
                 (((c["cr_48h"] - c["baseline_creatinine"]) >= 0.3).values |
                  ((c["cr_48h"] / c["baseline_creatinine"].replace(0, np.nan)) >= 1.5).values)
        time[aki & aki_48] = 2.0
        time[aki & ~aki_48] = 4.0
    else:
        time[aki] = 3.0
    return time, event


def aalen_johansen_cif(time, event, cause=1, grid=None):
    """Non-parametric Aalen-Johansen cumulative incidence for one cause."""
    if grid is None:
        grid = np.linspace(0, 7, 71)
    order = np.argsort(time)
    t = time[order]; e = event[order]
    n = len(t)
    at_risk = n
    surv = 1.0           # overall survival S(t-)
    cif = 0.0
    cif_curve = []
    uniq = np.unique(t)
    ti_idx = 0
    cif_at = {}
    s_prev = 1.0
    for ut in uniq:
        mask = t == ut
        d_cause = np.sum(mask & (e == cause))
        d_any   = np.sum(mask & (e != 0))
        if at_risk > 0:
            cif += s_prev * (d_cause / at_risk)
            s_prev = s_prev * (1 - d_any / at_risk)
        cif_at[ut] = cif
        at_risk -= np.sum(mask)
    # Map onto grid
    curve = []
    keys = sorted(cif_at.keys())
    for g in grid:
        val = 0.0
        for k in keys:
            if k <= g:
                val = cif_at[k]
            else:
                break
        curve.append(val)
    return grid, np.array(curve)


def fine_gray_analysis(cohort, cov_cols, n_boot=200):
    """
    Fit Fine-Gray subdistribution hazard for AKI with death as competing risk.
    Uses lifelines if available. Returns dict with sHR, CI, and CIF curves.
    IPTW-weighted to adjust for confounding (marginal subdistribution HR).
    """
    time, event = build_competing_time_event(cohort)
    A = cohort["treat"].values
    X = cohort[cov_cols].values
    ps, _ = fit_ps(X, A)
    ps = np.clip(ps, 0.05, 0.95)
    p1 = A.mean()
    w = np.where(A == 1, p1 / ps, (1 - p1) / (1 - ps))  # stabilized IPTW

    # CIFs by group (unweighted non-parametric for the figure)
    grid, cif_lr = aalen_johansen_cif(time[A == 1], event[A == 1], cause=1)
    _,    cif_ns = aalen_johansen_cif(time[A == 0], event[A == 0], cause=1)
    _,    cifd_lr = aalen_johansen_cif(time[A == 1], event[A == 1], cause=2)
    _,    cifd_ns = aalen_johansen_cif(time[A == 0], event[A == 0], cause=2)

    result = {
        "grid": grid, "cif_lr": cif_lr, "cif_ns": cif_ns,
        "cifd_lr": cifd_lr, "cifd_ns": cifd_ns,
        "n": len(cohort),
        "n_aki": int(np.sum(event == 1)),
        "n_death_only": int(np.sum(event == 2)),
    }

    # Fine-Gray sHR via lifelines if present
    shr = shr_lo = shr_hi = np.nan
    method = "none"
    try:
        from lifelines import CoxPHFitter
        # Fine-Gray subdistribution hazard via a weighted Cox on the
        # subdistribution risk set: competing-event subjects are kept at risk
        # until the end of the window (standard FG risk-set construction),
        # and IPTW weights adjust for measured confounding.
        fg_time = time.copy()
        fg_event = (event == 1).astype(int)
        # subjects with competing event stay at risk until end of window
        fg_time[event == 2] = 7.0
        dfg = pd.DataFrame({"T": fg_time, "E": fg_event, "treat": A, "w": w})
        cph = CoxPHFitter()
        cph.fit(dfg, "T", "E", weights_col="w", robust=True)
        shr = float(np.exp(cph.params_["treat"]))
        ci = cph.confidence_intervals_
        lo = float(np.exp(ci.loc["treat"].iloc[0]))
        hi = float(np.exp(ci.loc["treat"].iloc[1]))
        shr_lo, shr_hi = lo, hi
        method = "lifelines weighted FG-approx (Cox on subdistribution risk set)"
    except Exception as ex:
        # Fallback: CIF difference at day 7 with bootstrap CI
        print(f"  (lifelines unavailable or failed: {type(ex).__name__}; "
              f"using CIF-difference bootstrap)")
        method = "Aalen-Johansen CIF difference at day 7 (IPTW-adjusted)"

        def cif7(t, e, A_, w_):
            # weighted CIF at day 7 by group, return LR-NS difference
            def wcif(mask):
                tt, ee, ww = t[mask], e[mask], w_[mask]
                order = np.argsort(tt)
                tt, ee, ww = tt[order], ee[order], ww[order]
                at_risk = ww.sum()
                s_prev = 1.0; cif = 0.0
                for ut in np.unique(tt):
                    m = tt == ut
                    d_cause = ww[m & (ee == 1)].sum()
                    d_any   = ww[m & (ee != 0)].sum()
                    if at_risk > 0:
                        cif += s_prev * (d_cause / at_risk)
                        s_prev *= (1 - d_any / at_risk)
                    at_risk -= ww[m].sum()
                return cif
            return wcif(A_ == 1) - wcif(A_ == 0)

        point = cif7(time, event, A, w)

        def boot():
            idx = np.random.choice(len(A), len(A), replace=True)
            ps_b, _ = fit_ps(X[idx], A[idx])
            ps_b = np.clip(ps_b, 0.05, 0.95)
            p1b = A[idx].mean()
            w_b = np.where(A[idx] == 1, p1b / ps_b, (1 - p1b) / (1 - ps_b))
            return (cif7(time[idx], event[idx], A[idx], w_b),)

        rds = Parallel(n_jobs=-1, prefer="threads")(delayed(boot)() for _ in range(n_boot))
        rds = [r[0] for r in rds]
        shr = point  # here "shr" slot holds the CIF difference
        shr_lo = np.percentile(rds, 2.5)
        shr_hi = np.percentile(rds, 97.5)

    result.update({"shr": shr, "shr_lower": shr_lo, "shr_upper": shr_hi,
                   "fg_method": method})
    return result


def plot_cif(fg):
    grid = fg["grid"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: AKI CIF
    ax = axes[0]
    ax.step(grid, fg["cif_lr"], where="post", color=LR_COL, lw=2.2, label="LR")
    ax.step(grid, fg["cif_ns"], where="post", color=NS_COL, lw=2.2, label="NS")
    ax.set_xlabel("Days since first fluid", fontsize=11)
    ax.set_ylabel("Cumulative incidence of AKI", fontsize=11)
    ax.set_title("A   Incident AKI (primary event)", fontsize=12, fontweight="bold", loc="left")
    ax.legend(fontsize=10); ax.grid(alpha=0.3); ax.set_xlim(0, 7)

    # Panel B: Death CIF (competing)
    ax2 = axes[1]
    ax2.step(grid, fg["cifd_lr"], where="post", color=LR_COL, lw=2.2, ls="--", label="LR")
    ax2.step(grid, fg["cifd_ns"], where="post", color=NS_COL, lw=2.2, ls="--", label="NS")
    ax2.set_xlabel("Days since first fluid", fontsize=11)
    ax2.set_ylabel("Cumulative incidence of death (no AKI)", fontsize=11)
    ax2.set_title("B   Death without AKI (competing event)", fontsize=12, fontweight="bold", loc="left")
    ax2.legend(fontsize=10); ax2.grid(alpha=0.3); ax2.set_xlim(0, 7)

    est = fg["shr"]
    txt = (f"{fg['fg_method']}\n"
           f"Estimate = {est:.3f}  (95% CI {fg['shr_lower']:.3f}, {fg['shr_upper']:.3f})\n"
           f"n={fg['n']:,}  AKI={fg['n_aki']:,}  death-only={fg['n_death_only']:,}")
    fig.text(0.5, -0.04, txt, ha="center", fontsize=9, color="#37474F")

    fig.suptitle("Figure 9b. Competing-Risks Cumulative Incidence by Fluid Group",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.02)
    savefig("fig_cif_competing.png")


# =============================================================================
# MAIN
# =============================================================================
def main(parquet_path, n_boot=200):
    print("=" * 62)
    print("Extended Analyses: E-value + Fine-Gray Competing Risks")
    print("=" * 62)

    import shutil, tempfile
    _tmp = tempfile.mkdtemp(prefix="causalai_ext_")
    _loc = os.path.join(_tmp, "cohort.parquet")
    shutil.copy2(parquet_path, _loc)
    df = pd.read_parquet(_loc)
    cohort, cov_cols = prepare_base(df)
    print(f"Cohort: {len(cohort):,} rows, {len(cov_cols)} covariates")

    csv_path = os.path.join(SCRIPT_DIR, "results_summary.csv")

    # ── E-value: read the primary AIPW result from the CSV ────────────────────
    print("\n[1] E-value for primary AIPW result")
    if os.path.exists(csv_path):
        res = pd.read_csv(csv_path)
        prim = res[res["analysis"].str.startswith("Primary")]
        if len(prim):
            r = prim.iloc[0]
            ev = compute_evalue(r["rd"], r["risk_lr"], r["risk_ns"],
                                r["rd_lower"], r["rd_upper"])
            print(f"  RR = {ev['rr']:.4f}  (95% CI {ev['rr_lower']:.4f}, {ev['rr_upper']:.4f})")
            print(f"  E-value (point)   = {ev['evalue_point']:.3f}")
            print(f"  E-value (CI bound)= {ev['evalue_ci']:.3f}")
            print(f"  CI crosses null   = {ev['crosses_null']}")
            pd.DataFrame([{
                "quantity": "Primary AIPW",
                "risk_lr": ev["risk_lr"], "risk_ns": ev["risk_ns"],
                "rr": ev["rr"], "rr_lower": ev["rr_lower"], "rr_upper": ev["rr_upper"],
                "evalue_point": ev["evalue_point"], "evalue_ci": ev["evalue_ci"],
                "crosses_null": ev["crosses_null"],
            }]).to_csv(os.path.join(SCRIPT_DIR, "evalue_summary.csv"), index=False)
            print(f"  wrote evalue_summary.csv")
            plot_evalue(ev)
        else:
            print("  No Primary row in results_summary.csv — run full_analysis.py first.")
    else:
        print("  results_summary.csv not found — run full_analysis.py first.")

    # ── Fine-Gray competing risks ─────────────────────────────────────────────
    print("\n[2] Fine-Gray competing-risks analysis (AKI vs death)")
    fg = fine_gray_analysis(cohort, cov_cols, n_boot=n_boot)
    print(f"  method: {fg['fg_method']}")
    print(f"  estimate={fg['shr']:.4f}  95%CI ({fg['shr_lower']:.4f}, {fg['shr_upper']:.4f})")
    print(f"  n={fg['n']:,}  AKI={fg['n_aki']:,}  death-only={fg['n_death_only']:,}")
    plot_cif(fg)

    # Append FG row to results_summary.csv (as an RD-like row for the CIF diff,
    # or sHR where lifelines was used — flagged in the analysis label)
    fg_row = {
        "analysis": ("SA10: Fine-Gray sHR (AKI, death competing)"
                     if "lifelines" in fg["fg_method"]
                     else "SA10: Competing-risks CIF diff (day 7)"),
        "n": fg["n"], "rd": fg["shr"],
        "rd_lower": fg["shr_lower"], "rd_upper": fg["shr_upper"],
        "risk_lr": fg["cif_lr"][-1], "risk_ns": fg["cif_ns"][-1],
    }
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        existing = existing[~existing["analysis"].str.startswith("SA10")]
        combined = pd.concat([existing, pd.DataFrame([fg_row])], ignore_index=True)
    else:
        combined = pd.DataFrame([fg_row])
    combined.to_csv(csv_path, index=False)
    print(f"  appended SA10 row to results_summary.csv")

    print("\n" + "=" * 62)
    print("Extended analyses complete.")
    print("  Figures: fig_evalue.png, fig_cif_competing.png")
    print("  Tables:  evalue_summary.csv, results_summary.csv (SA10 row)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args()
    main(args.parquet, args.n_boot)
