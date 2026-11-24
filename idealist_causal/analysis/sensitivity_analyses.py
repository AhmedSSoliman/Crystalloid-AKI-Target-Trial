"""
IDEALIST LR vs NS  —  Sensitivity / subgroup / PS-matching analyses
====================================================================
SA5  : 3-day AKI window (aki_surgery_3d)
SA7  : Sex subgroups
SA8  : Surgery-type subgroups
SA9  : Emergent vs elective
PS matching: 1:1 NN, caliper 0.2*SD(logit PS)
E-value for unmeasured confounding
Fine-Gray competing-risks (death without AKI)

Run after full_analysis.py:
    python3 sensitivity_analyses.py
"""

import os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed
from scipy.stats import norm

SEED = 42
np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_CSV    = os.path.join(SCRIPT_DIR, "results_summary_idealist.csv")
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")

PS_LOW, PS_HIGH = 0.05, 0.95
N_BOOT = 200

# Load COV_COLS from full_analysis
from full_analysis import COV_COLS, prep_X, fit_ps, aipw_estimate, bootstrap_ci, append_result


# ── PS matching ────────────────────────────────────────────────────────────────
def ps_match_1to1(cohort, ps_full):
    """1:1 nearest-neighbour matching, caliper = 0.2*SD(logit PS)."""
    logit_ps = np.log(np.clip(ps_full, 1e-6, 1 - 1e-6) /
                      np.clip(1 - ps_full, 1e-6, 1))
    caliper = 0.2 * np.std(logit_ps)

    A = cohort["fluid_type"].eq("LR").astype(int).values
    idx_lr = np.where(A == 1)[0]
    idx_ns = np.where(A == 0)[0]

    # Match each LR to nearest NS within caliper
    nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
    nn.fit(logit_ps[idx_ns].reshape(-1, 1))
    dist, ind = nn.kneighbors(logit_ps[idx_lr].reshape(-1, 1))

    matched_lr = []; matched_ns = []
    for i, (d, j) in enumerate(zip(dist[:, 0], ind[:, 0])):
        if d <= caliper:
            matched_lr.append(idx_lr[i])
            matched_ns.append(idx_ns[j])

    # Drop duplicate NS matches
    ns_seen = set(); final_lr = []; final_ns = []
    for lr_i, ns_i in zip(matched_lr, matched_ns):
        if ns_i not in ns_seen:
            final_lr.append(lr_i); final_ns.append(ns_i)
            ns_seen.add(ns_i)

    m_lr = cohort.iloc[final_lr].reset_index(drop=True)
    m_ns = cohort.iloc[final_ns].reset_index(drop=True)
    return m_lr, m_ns


def matched_rd_ci(m_lr, m_ns, outcome="aki_surgery_7d", n_boot=N_BOOT):
    """Risk difference in matched pairs with bootstrap CI."""
    matched = pd.concat([m_lr.assign(_arm="LR"),
                         m_ns.assign(_arm="NS")], ignore_index=True)
    rd = m_lr[outcome].mean() - m_ns[outcome].mean()

    rds_boot = []
    rng = np.random.default_rng(SEED)
    for _ in range(n_boot):
        pairs = rng.integers(0, len(m_lr), len(m_lr))
        rds_boot.append(m_lr[outcome].iloc[pairs].mean() -
                        m_ns[outcome].iloc[pairs].mean())
    lo = float(np.nanpercentile(rds_boot, 2.5))
    hi = float(np.nanpercentile(rds_boot, 97.5))
    return rd, lo, hi


# ── E-value ────────────────────────────────────────────────────────────────────
def evalue_from_rr(rr):
    if rr <= 0: return np.nan
    if rr >= 1:
        return rr + np.sqrt(rr * (rr - 1))
    else:
        rr_inv = 1.0 / rr
        return rr_inv + np.sqrt(rr_inv * (rr_inv - 1))


def compute_evalue(rd, r_lr, r_ns):
    """E-value on RR scale from marginal risks."""
    if r_ns <= 0 or r_ns >= 1:
        return np.nan, np.nan
    rr = r_lr / r_ns
    return evalue_from_rr(rr), rr


# ── Fine-Gray (IPTW-weighted subdistribution) ──────────────────────────────────
def fine_gray_analysis(cohort_full, ps_full):
    """
    Competing-risks analysis.
    Event 1 = AKI within 7 days
    Event 2 = death without AKI (competing event)
    Uses IPTW-weighted Cox on subdistribution risk set via lifelines.
    Falls back to CIF difference if lifelines unavailable.
    """
    try:
        from lifelines import CoxPHFitter
        HAS_LIFELINES = True
    except ImportError:
        HAS_LIFELINES = False

    # Keep trimmed cohort
    keep = (ps_full >= PS_LOW) & (ps_full <= PS_HIGH)
    ct   = cohort_full[keep].reset_index(drop=True)
    ps_t = ps_full[keep]

    A = ct["fluid_type"].eq("LR").astype(int).values
    p1 = ps_t.mean(); p0 = 1 - p1
    w = np.where(A == 1, p1 / np.clip(ps_t, 1e-6, 1),
                         p0 / np.clip(1 - ps_t, 1e-6, 1))

    # Build competing-event indicator
    # event=1 AKI, event=2 death-without-AKI, event=0 censored
    aki  = ct["aki_surgery_7d"].fillna(0).astype(int).values
    dead = ct["hospital_mortality"].fillna(0).astype(int).values
    event = np.where(aki == 1, 1, np.where(dead == 1, 2, 0))

    # Use coarse time: 3 for event, 7 for censored (calendar days not available exactly)
    time = np.where(event > 0, 3.5, 7.0)

    df_fg = pd.DataFrame({
        "T": time, "E": event, "A": A, "W": w,
    })

    if HAS_LIFELINES:
        # Fine-Gray via weighted Cox on sub-distribution risk set
        # event=2 → set time to 7 (keeps competing-event subjects at risk)
        df_fg["T_fg"] = np.where(df_fg["E"] == 2, 7.0, df_fg["T"])
        df_fg["E_fg"] = (df_fg["E"] == 1).astype(int)
        cph = CoxPHFitter()
        cph.fit(df_fg[["T_fg", "E_fg", "A", "W"]], duration_col="T_fg",
                event_col="E_fg", weights_col="W", robust=True)
        shr   = float(np.exp(cph.params_["A"]))
        ci_lo = float(np.exp(cph.confidence_intervals_.loc["A", "95% lower-bound"]))
        ci_hi = float(np.exp(cph.confidence_intervals_.loc["A", "95% upper-bound"]))
        return shr, ci_lo, ci_hi, df_fg
    else:
        # Fallback: IPTW-adjusted CIF difference
        from full_analysis import iptw_only
        rd_cif, r_lr, r_ns = iptw_only(ct, ps_t, "aki_surgery_7d")
        return 1.0 + rd_cif, None, None, df_fg


def main():
    print("Loading trimmed cohort and PS arrays…")
    ct        = pd.read_parquet(os.path.join(SCRIPT_DIR, "cohort_trimmed_idealist.parquet"))
    cohort    = pd.read_parquet(os.path.join(SCRIPT_DIR, "cohort_full_ps_idealist.parquet"))
    ps_full   = np.load(os.path.join(SCRIPT_DIR, "ps_full_idealist.npy"))
    ps_t      = np.load(os.path.join(SCRIPT_DIR, "ps_trim_idealist.npy"))

    # ── SA5: 3-day AKI window ─────────────────────────────────────────────────
    print("\nSA5: 3-day AKI window…")
    if "aki_surgery_3d_without_race_correction" in ct.columns:
        ct["aki_3d"] = ct["aki_surgery_3d_without_race_correction"].fillna(0).astype(float)
        rd5, r5_lr, r5_ns = aipw_estimate(ct, ps_t, "aki_3d")
        lo5, hi5 = bootstrap_ci(ct, "aki_3d")
        print(f"  RD = {rd5*100:+.2f} pp  CI: ({lo5*100:+.2f}, {hi5*100:+.2f})")
        append_result(dict(analysis="SA5: 3-day AKI window",
                           n=len(ct), rd=rd5, rd_lower=lo5, rd_upper=hi5,
                           risk_lr=r5_lr, risk_ns=r5_ns))
    else:
        print("  aki_surgery_3d not found — skipping SA5")

    # ── SA7: Sex subgroups ────────────────────────────────────────────────────
    if "sex" in ct.columns:
        for sex_val in ct["sex"].dropna().unique():
            sub = ct[ct["sex"] == sex_val].reset_index(drop=True)
            if len(sub) < 200: continue
            ps_sub = fit_ps(prep_X(sub, COV_COLS),
                            sub.fluid_type.eq("LR").astype(int).values)
            keep_s = (ps_sub >= PS_LOW) & (ps_sub <= PS_HIGH)
            sub2 = sub[keep_s].reset_index(drop=True); ps2 = ps_sub[keep_s]
            rd_s, r_lr_s, r_ns_s = aipw_estimate(sub2, ps2)
            lo_s, hi_s = bootstrap_ci(sub2, "aki_surgery_7d", n_boot=100)
            print(f"\nSA7 sex={sex_val}: n={len(sub2):,}  RD={rd_s*100:+.2f}  "
                  f"CI:({lo_s*100:+.2f},{hi_s*100:+.2f})")
            append_result(dict(analysis=f"SA7: {sex_val} subgroup",
                               n=len(sub2), rd=rd_s, rd_lower=lo_s, rd_upper=hi_s,
                               risk_lr=r_lr_s, risk_ns=r_ns_s))

    # ── SA8: Surgery type subgroups ───────────────────────────────────────────
    stype_col = next((c for c in ["surgery_type", "specialty", "sched_service"]
                      if c in ct.columns), None)
    if stype_col:
        top_types = ct[stype_col].value_counts().head(4).index.tolist()
        for stype in top_types:
            sub = ct[ct[stype_col] == stype].reset_index(drop=True)
            if len(sub) < 100: continue
            ps_sub = fit_ps(prep_X(sub, COV_COLS),
                            sub.fluid_type.eq("LR").astype(int).values)
            keep_s = (ps_sub >= PS_LOW) & (ps_sub <= PS_HIGH)
            sub2 = sub[keep_s].reset_index(drop=True); ps2 = ps_sub[keep_s]
            rd_s, r_lr_s, r_ns_s = aipw_estimate(sub2, ps2)
            lo_s, hi_s = bootstrap_ci(sub2, "aki_surgery_7d", n_boot=100)
            print(f"\nSA8 {stype_col}={stype}: n={len(sub2):,}  RD={rd_s*100:+.2f}  "
                  f"CI:({lo_s*100:+.2f},{hi_s*100:+.2f})")
            append_result(dict(analysis=f"SA8: {stype}",
                               n=len(sub2), rd=rd_s, rd_lower=lo_s, rd_upper=hi_s,
                               risk_lr=r_lr_s, risk_ns=r_ns_s))

    # ── SA9: Emergent vs elective ─────────────────────────────────────────────
    if "emergent" in ct.columns:
        for emerg_val, label in [(1, "Emergent"), (0, "Elective")]:
            sub = ct[ct["emergent"] == emerg_val].reset_index(drop=True)
            if len(sub) < 100: continue
            ps_sub = fit_ps(prep_X(sub, COV_COLS),
                            sub.fluid_type.eq("LR").astype(int).values)
            keep_s = (ps_sub >= PS_LOW) & (ps_sub <= PS_HIGH)
            sub2 = sub[keep_s].reset_index(drop=True); ps2 = ps_sub[keep_s]
            rd_s, r_lr_s, r_ns_s = aipw_estimate(sub2, ps2)
            lo_s, hi_s = bootstrap_ci(sub2, "aki_surgery_7d", n_boot=100)
            print(f"\nSA9 {label}: n={len(sub2):,}  RD={rd_s*100:+.2f}  "
                  f"CI:({lo_s*100:+.2f},{hi_s*100:+.2f})")
            append_result(dict(analysis=f"SA9: {label}",
                               n=len(sub2), rd=rd_s, rd_lower=lo_s, rd_upper=hi_s,
                               risk_lr=r_lr_s, risk_ns=r_ns_s))

    # ── PS matching ───────────────────────────────────────────────────────────
    print("\nPS matching (1:1 NN)…")
    keep = (ps_full >= PS_LOW) & (ps_full <= PS_HIGH)
    cohort_for_match = cohort[keep].reset_index(drop=True)
    ps_for_match = ps_full[keep]
    m_lr, m_ns = ps_match_1to1(cohort_for_match, ps_for_match)
    n_pairs = len(m_lr)
    rd_m, lo_m, hi_m = matched_rd_ci(m_lr, m_ns)
    r_lr_m = m_lr["aki_surgery_7d"].mean()
    r_ns_m = m_ns["aki_surgery_7d"].mean()
    print(f"  {n_pairs:,} matched pairs  RD={rd_m*100:+.2f}  CI:({lo_m*100:+.2f},{hi_m*100:+.2f})")
    append_result(dict(analysis=f"PS Matching (1:1 NN caliper=0.2SD)",
                       n=n_pairs, rd=rd_m, rd_lower=lo_m, rd_upper=hi_m,
                       risk_lr=r_lr_m, risk_ns=r_ns_m))

    # ── Fine-Gray competing risks ─────────────────────────────────────────────
    print("\nFine-Gray competing-risks analysis…")
    shr, ci_lo, ci_hi, df_fg = fine_gray_analysis(cohort, ps_full)
    r_lr_fg = cohort[keep]["aki_surgery_7d"].mean()
    r_ns_fg = cohort[keep][cohort[keep].fluid_type == "NS"]["aki_surgery_7d"].mean()
    print(f"  Fine-Gray sHR = {shr:.3f}  CI: ({ci_lo:.3f}, {ci_hi:.3f})" if ci_lo else
          f"  IPTW-adjusted CIF difference (lifelines fallback): {shr:.4f}")
    append_result(dict(
        analysis="SA10: Fine-Gray sHR (AKI, death competing)",
        n=len(cohort[keep]),
        rd=shr if ci_lo else shr - 1.0,
        rd_lower=ci_lo if ci_lo else np.nan,
        rd_upper=ci_hi if ci_hi else np.nan,
        risk_lr=r_lr_fg, risk_ns=r_ns_fg,
    ))
    df_fg.to_csv(os.path.join(SCRIPT_DIR, "fg_data_idealist.csv"), index=False)

    # ── E-value ───────────────────────────────────────────────────────────────
    print("\nE-value…")
    # Read primary estimate from results
    res = pd.read_csv(OUT_CSV)
    prim = res[res.analysis == "Primary AIPW (PS-trimmed)"].iloc[0]
    ev, rr = compute_evalue(prim.rd, prim.risk_lr, prim.risk_ns)
    print(f"  Primary RR = {rr:.4f}  E-value = {ev:.3f}")
    ev_df = pd.DataFrame([dict(
        quantity="Primary AIPW",
        risk_lr=prim.risk_lr, risk_ns=prim.risk_ns,
        rr=rr, evalue_point=ev, crosses_null=abs(prim.rd) < 0.005
    )])
    ev_df.to_csv(os.path.join(SCRIPT_DIR, "evalue_summary_idealist.csv"), index=False)

    print("\n✓ Sensitivity analyses complete.")
    print(f"Results: {OUT_CSV}")


if __name__ == "__main__":
    main()
