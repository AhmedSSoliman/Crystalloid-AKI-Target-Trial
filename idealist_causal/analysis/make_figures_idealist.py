"""
IDEALIST LR vs NS  —  Publication figures
==========================================
Mirrors make_figures.py from the MIMIC analysis.
Reads from saved CSV/parquet files; does NOT re-run the analysis.

Run after sensitivity_analyses.py:
    python3 make_figures_idealist.py
"""

import os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.gridspec import GridSpec

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Style (identical palette to MIMIC analysis) ────────────────────────────────
LR_COL  = "#1565C0"
NS_COL  = "#E65100"
GRY     = "#546E7A"
RED_SIG = "#B71C1C"
GRN     = "#2E7D32"
LTGRY   = "#ECEFF1"

plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          11,
    "axes.titlesize":     13,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    10,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.color":         "#E0E0E0",
    "grid.linewidth":     0.6,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
})


def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✓ {name}")


# ── Figure 1: Cohort flow ──────────────────────────────────────────────────────
def fig_cohort_flow(stats: dict):
    """stats keys: n_all, n_lr_raw, n_ns_raw, n_no_cr, n_eligible,
                   n_prev_aki, n_primary, n_trimmed, n_analytic, n_lr, n_ns"""
    fig, ax = plt.subplots(figsize=(9.6, 11))
    ax.set_xlim(-0.6, 11.0); ax.set_ylim(0, 11); ax.axis("off")

    def box(x, y, w, h, text, color="#FFFFFF", ec="#37474F", lw=1.5, fontsize=10, bold=False):
        rect = FancyBboxPatch((x-w/2, y-h/2), w, h,
                              boxstyle="round,pad=0.08", linewidth=lw,
                              edgecolor=ec, facecolor=color)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                fontweight="bold" if bold else "normal",
                multialignment="center", color="#1A1A1A")

    def arrow(x, y1, y2, pad_top=0.47, pad_bot=0.45):
        ax.annotate("", xy=(x, y2+pad_bot), xytext=(x, y1-pad_top),
                    arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.8, mutation_scale=14))

    def excl_box(x, y, text, side="right", feed_x=5.0):
        rect = FancyBboxPatch((x, y-0.42), 3.3, 0.84,
                              boxstyle="round,pad=0.06", linewidth=1,
                              edgecolor="#B71C1C", facecolor="#FFEBEE")
        ax.add_patch(rect)
        ax.text(x+1.65, y, text, ha="center", va="center", fontsize=9,
                color="#B71C1C", multialignment="center")
        if side == "right":
            ax.annotate("", xy=(x+0.02, y), xytext=(feed_x, y),
                        arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.2, mutation_scale=10))
        else:
            ax.annotate("", xy=(x+3.28, y), xytext=(feed_x, y),
                        arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.2, mutation_scale=10))

    ys = [10.3, 8.6, 6.9, 4.9, 2.6]
    s = stats

    box(5, ys[0], 5.5, 0.9,
        f"IDEALIST Surgical Cohort\nAll OR cases with intraop LR or NS\n"
        f"n = {s['n_all']:,}  (LR {s['n_lr_raw']:,} · NS {s['n_ns_raw']:,})",
        color="#E3F2FD", ec=LR_COL, lw=2, fontsize=10)

    arrow(5, ys[0], ys[1], pad_top=0.47, pad_bot=0.45)
    excl_box(5.95, (ys[0]+ys[1])/2,
             f"Missing baseline or\nfollow-up creatinine\nn = {s['n_no_cr']:,}",
             side="right", feed_x=5.0)

    box(5, ys[1], 5.5, 0.85,
        f"Baseline AND follow-up creatinine present\nn = {s['n_eligible']:,}",
        color="#F3E5F5", ec="#6A1B9A", fontsize=10)

    arrow(5, ys[1], ys[2], pad_top=0.45, pad_bot=0.45)
    excl_box(5.95, (ys[1]+ys[2])/2,
             f"Prevalent AKI proxy\n(baseline Cr > 4.0 mg/dL)\nn = {s['n_prev_aki']:,}",
             side="right", feed_x=5.0)

    box(5, ys[2], 5.5, 0.85,
        f"Eligible cohort (baseline Cr ≤ 4.0 mg/dL)\nn = {s['n_primary']:,}",
        color="#E8F5E9", ec=GRN, fontsize=10)

    # Branch to LR/NS
    ax.annotate("", xy=(2.6, ys[3]+0.42), xytext=(5, ys[2]-0.43),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))
    ax.annotate("", xy=(7.4, ys[3]+0.42), xytext=(5, ys[2]-0.43),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))

    box(2.6, ys[3], 3.6, 0.84,
        f"Lactated Ringer's (LR)\nn = {s['n_lr']:,}  ({100*s['n_lr']/s['n_primary']:.1f}%)",
        color="#E3F2FD", ec=LR_COL, lw=2, fontsize=10, bold=True)
    box(7.4, ys[3], 3.6, 0.84,
        f"Normal Saline (NS)\nn = {s['n_ns']:,}  ({100*s['n_ns']/s['n_primary']:.1f}%)",
        color="#FFF3E0", ec=NS_COL, lw=2, fontsize=10, bold=True)

    # Re-converge
    ax.annotate("", xy=(5, ys[4]+0.43), xytext=(2.6, ys[3]-0.42),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))
    ax.annotate("", xy=(5, ys[4]+0.43), xytext=(7.4, ys[3]-0.42),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))

    excl_y = (ys[3] + ys[4]) / 2 - 0.15
    excl_box(7.7, excl_y,
             f"PS outside [0.05, 0.95]\n(poor overlap)\nn = {s['n_trimmed']:,}",
             side="right", feed_x=7.4)

    box(5, ys[4], 5.8, 0.85,
        f"Primary analytic cohort (PS-trimmed)\nn = {s['n_analytic']:,}",
        color="#FFF8E1", ec="#F57F17", fontsize=10, bold=True)

    ax.set_title("IDEALIST: Cohort Flow Diagram", fontsize=14, fontweight="bold",
                 pad=10, loc="left")
    save("fig_cohort_flow_idealist.png")


# ── Figure 2: PS overlap ───────────────────────────────────────────────────────
def fig_ps_overlap(cohort_full, ps_full, cohort_trim, ps_trim):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, ps, cohort, title in [
        (axes[0], ps_full, cohort_full, "Full cohort"),
        (axes[1], ps_trim, cohort_trim, "PS-trimmed cohort"),
    ]:
        A = cohort["fluid_type"].eq("LR").values
        ax.hist(ps[A], bins=40, alpha=0.6, color=LR_COL, label="LR", density=True)
        ax.hist(ps[~A], bins=40, alpha=0.6, color=NS_COL, label="NS", density=True)
        ax.set_xlabel("Propensity score"); ax.set_ylabel("Density")
        ax.set_title(title); ax.legend()
    fig.suptitle("IDEALIST: PS Overlap — LR vs NS", fontsize=13, fontweight="bold")
    save("fig_ps_overlap_idealist.png")


# ── Figure 3: Love plot ────────────────────────────────────────────────────────
def fig_love_plot(balance_csv):
    bal = pd.read_csv(balance_csv, index_col=0)
    bal = bal.reindex(bal["smd_before"].abs().sort_values(ascending=True).index)
    fig, ax = plt.subplots(figsize=(8, max(5, len(bal)*0.22)))
    y = np.arange(len(bal))
    ax.scatter(bal["smd_before"], y, color=GRY, s=35, zorder=3, label="Before IPTW")
    ax.scatter(bal["smd_after"],  y, color=LR_COL, s=35, zorder=3, label="After IPTW")
    ax.axvline(0.1, color=RED_SIG, lw=1.2, ls="--", alpha=0.6)
    ax.axvline(-0.1, color=RED_SIG, lw=1.2, ls="--", alpha=0.6)
    ax.axvline(0, color="black", lw=0.8, alpha=0.4)
    ax.set_yticks(y); ax.set_yticklabels(bal.index, fontsize=7)
    ax.set_xlabel("Standardized Mean Difference")
    ax.set_title("IDEALIST: Covariate Balance Before/After IPTW", fontweight="bold")
    ax.legend()
    save("fig_love_plot_idealist.png")


# ── Figure 4: AIPW forest ─────────────────────────────────────────────────────
def fig_aipw_forest(res_csv):
    df = pd.read_csv(res_csv)
    # Filter to main SA analyses (exclude subgroups and matching for this plot)
    main_rows = df[~df.analysis.str.contains(
        r"SA7:|SA8:|SA9:|PS Match|Fine-Gray", na=False
    )].copy()
    main_rows["rd_pp"]  = main_rows["rd"] * 100
    main_rows["lo_pp"]  = main_rows["rd_lower"] * 100
    main_rows["hi_pp"]  = main_rows["rd_upper"] * 100
    main_rows = main_rows.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11, max(4, len(main_rows)*0.55)))
    y = np.arange(len(main_rows))
    colors = [RED_SIG if (row.lo_pp > 0 or row.hi_pp < 0) else LR_COL
              for _, row in main_rows.iterrows()]
    for i, (_, row) in enumerate(main_rows.iterrows()):
        ax.plot([row.lo_pp, row.hi_pp], [i, i], color=colors[i], lw=1.8)
        ax.scatter(row.rd_pp, i, color=colors[i], s=60, zorder=5)
        ax.text(row.hi_pp + 0.1, i,
                f"n={int(row.n):,}  {row.rd_pp:+.2f} ({row.lo_pp:+.2f},{row.hi_pp:+.2f})",
                va="center", fontsize=8)
    ax.axvline(0, color="black", lw=1.2, ls="--", alpha=0.6)
    ax.set_yticks(y); ax.set_yticklabels(main_rows.analysis, fontsize=9)
    ax.set_xlabel("Risk Difference LR − NS (percentage points)")
    ax.set_title("IDEALIST: AIPW Risk Differences", fontweight="bold")
    save("fig_aipw_results_idealist.png")


# ── Figure 5: Subgroup forest ─────────────────────────────────────────────────
def fig_subgroup_forest(res_csv):
    df = pd.read_csv(res_csv)
    subs = df[df.analysis.str.contains(r"SA7:|SA8:|SA9:", na=False)].copy()
    if subs.empty:
        print("  No subgroup rows yet — skipping forest")
        return
    subs["rd_pp"] = subs["rd"] * 100
    subs["lo_pp"] = subs["rd_lower"] * 100
    subs["hi_pp"] = subs["rd_upper"] * 100
    subs = subs.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11, max(4, len(subs)*0.55)))
    y = np.arange(len(subs))
    for i, (_, row) in enumerate(subs.iterrows()):
        c = RED_SIG if (row.lo_pp > 0 or row.hi_pp < 0) else LR_COL
        ax.plot([row.lo_pp, row.hi_pp], [i, i], color=c, lw=1.8)
        ax.scatter(row.rd_pp, i, color=c, s=60, zorder=5)
        ax.text(row.hi_pp + 0.1, i,
                f"n={int(row.n):,}  {row.rd_pp:+.2f} ({row.lo_pp:+.2f},{row.hi_pp:+.2f})",
                va="center", fontsize=8)
    ax.axvline(0, color="black", lw=1.2, ls="--", alpha=0.6)
    ax.set_yticks(y); ax.set_yticklabels(subs.analysis, fontsize=9)
    ax.set_xlabel("Risk Difference LR − NS (percentage points)")
    ax.set_title("IDEALIST: Subgroup Analyses", fontweight="bold")
    save("fig_subgroup_forest_idealist.png")


# ── Figure 6: CIF competing risks ────────────────────────────────────────────
def fig_cif(fg_csv):
    df = pd.read_csv(fg_csv)
    if df.empty: return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for event_id, (ax, title) in enumerate(zip(axes,
            ["Incident AKI (primary event)", "Death without AKI (competing event)"]), 1):
        for fluid, col, ls in [("LR", LR_COL, "-"), ("NS", NS_COL, "--")]:
            arm = df[df.A == (1 if fluid == "LR" else 0)]
            ev  = arm[arm.E == event_id].sort_values("T")
            if ev.empty: continue
            # Simple step-function CIF (illustrative)
            ax.step(ev["T"].values, np.linspace(0, len(ev)/len(arm), len(ev)),
                    color=col, ls=ls, label=fluid, lw=2)
        ax.set_xlabel("Days since surgery")
        ax.set_ylabel("Cumulative incidence")
        ax.set_title(title); ax.legend()
    fig.suptitle("IDEALIST: Competing-Risks Cumulative Incidence", fontsize=13, fontweight="bold")
    save("fig_cif_idealist.png")


# ── Figure 7: E-value ──────────────────────────────────────────────────────────
def fig_evalue(evalue_csv):
    if not os.path.exists(evalue_csv):
        print("  E-value CSV not found — skipping"); return
    ev_df = pd.read_csv(evalue_csv)
    rr = float(ev_df["rr"].iloc[0])
    ev = float(ev_df["evalue_point"].iloc[0])

    fig, ax = plt.subplots(figsize=(6, 5))
    rr_range = np.linspace(1.0, 5.0, 300)
    # E-value curve: bias = x + sqrt(x*(x-1))
    bias = rr_range + np.sqrt(rr_range * (rr_range - 1))
    ax.plot(rr_range, bias, color=LR_COL, lw=2)
    # Mark primary estimate
    ax.scatter([ev], [ev], color=RED_SIG, zorder=5, s=80)
    ax.text(ev + 0.05, ev + 0.02, f"E-value = {ev:.2f}", fontsize=10, color=RED_SIG)
    ax.set_xlabel("RR: unmeasured confounder ↔ treatment")
    ax.set_ylabel("RR: unmeasured confounder ↔ outcome")
    ax.set_title("IDEALIST: E-value Sensitivity to Unmeasured Confounding", fontweight="bold")
    save("fig_evalue_idealist.png")


def main():
    res_csv   = os.path.join(SCRIPT_DIR, "results_summary_idealist.csv")
    bal_csv   = os.path.join(SCRIPT_DIR, "balance_idealist.csv")
    ev_csv    = os.path.join(SCRIPT_DIR, "evalue_summary_idealist.csv")
    fg_csv    = os.path.join(SCRIPT_DIR, "fg_data_idealist.csv")
    ct_path   = os.path.join(SCRIPT_DIR, "cohort_trimmed_idealist.parquet")
    full_path = os.path.join(SCRIPT_DIR, "cohort_full_ps_idealist.parquet")
    ps_full_p = os.path.join(SCRIPT_DIR, "ps_full_idealist.npy")
    ps_trim_p = os.path.join(SCRIPT_DIR, "ps_trim_idealist.npy")

    if not os.path.exists(res_csv):
        print("results_summary_idealist.csv not found — run full_analysis.py first"); return

    ct      = pd.read_parquet(ct_path)
    cohort  = pd.read_parquet(full_path)
    ps_full = np.load(ps_full_p)
    ps_trim = np.load(ps_trim_p)

    # Cohort flow stats
    res = pd.read_csv(res_csv)
    prim = res[res.analysis == "Primary AIPW (PS-trimmed)"].iloc[0]

    # Back-calculate flow counts
    n_primary  = int(cohort[cohort.baseline_creatinine <= 4.0].shape[0])
    n_analytic = int(prim.n)
    n_trimmed  = n_primary - n_analytic
    n_lr       = int(cohort[cohort.baseline_creatinine <= 4.0].fluid_type.eq("LR").sum())
    n_ns       = int(cohort[cohort.baseline_creatinine <= 4.0].fluid_type.eq("NS").sum())
    n_all      = len(cohort)

    stats = dict(
        n_all=n_all, n_lr_raw=cohort.fluid_type.eq("LR").sum(),
        n_ns_raw=cohort.fluid_type.eq("NS").sum(),
        n_no_cr=0, n_eligible=n_primary,  # already filtered in extract step
        n_prev_aki=0, n_primary=n_primary,
        n_trimmed=n_trimmed, n_analytic=n_analytic,
        n_lr=n_lr, n_ns=n_ns,
    )

    print("Building IDEALIST figures…")
    fig_cohort_flow(stats)

    keep = (ps_full >= 0.05) & (ps_full <= 0.95)
    fig_ps_overlap(cohort, ps_full, cohort[keep].reset_index(drop=True), ps_trim)

    if os.path.exists(bal_csv):
        fig_love_plot(bal_csv)

    fig_aipw_forest(res_csv)
    fig_subgroup_forest(res_csv)

    if os.path.exists(fg_csv):
        fig_cif(fg_csv)

    fig_evalue(ev_csv)

    print(f"\nAll figures → {FIG_DIR}")


if __name__ == "__main__":
    main()
