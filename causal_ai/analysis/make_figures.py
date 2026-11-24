"""
Publication-quality figures for LR vs NS → 7-day AKI causal analysis.
Replaces all auto-generated figures with consistent, journal-ready versions.

Run:
    python3 make_figures.py --parquet ../../parquet/cohort_enriched.parquet
"""

import argparse, os, shutil, tempfile, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Ellipse
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score

# ── Global style ──────────────────────────────────────────────────────────────
FONT = "DejaVu Sans"
plt.rcParams.update({
    "font.family":        FONT,
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

# Palette
LR_COL  = "#1565C0"   # deep blue
NS_COL  = "#E65100"   # deep orange
GRY     = "#546E7A"   # blue-grey
RED_SIG = "#B71C1C"   # dark red  (significant result)
GRN     = "#2E7D32"   # dark green
LTGRY   = "#ECEFF1"

SEED = 42
np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✓ {name}")

# ── Human-readable variable labels ───────────────────────────────────────────
VAR_LABELS = {
    # Demographics
    "female":               "Female sex",
    "anchor_age":           "Age (years)",
    "race_white":           "Race: White",
    "race_black":           "Race: Black",
    # Kidney
    "log_baseline_cr":      "Log baseline creatinine",
    "ckd_proxy":            "CKD proxy (Cr >1.5)",
    "ckd":                  "CKD (ICD-coded)",
    # Timing
    "hours_from_icu_admit": "Hours to first fluid",
    "los_hours":            "ICU LOS at admission (h)",
    # Comorbidities
    "liver_disease":        "Liver disease",
    "heart_failure":        "Heart failure",
    "sepsis":               "Sepsis",
    "cancer":               "Cancer",
    "diabetes":             "Diabetes",
    "electrolyte_disorder": "Electrolyte disorder",
    # Severity / vitals
    "vasopressor_pre":      "Vasopressor use pre-t₀",
    "gcs_min_pre":          "GCS min (24h pre-t₀)",
    "hr_max_pre":           "HR max (24h pre-t₀)",
    "spo2_min_pre":         "SpO₂ min (24h pre-t₀)",
    # Creatinine trajectory
    "cr_slope_pre":         "Creatinine slope (7d pre-t₀)",
    "n_cr_pre":             "No. prior creatinine values",
    # Calendar
    "year":                 "Calendar year",
}

# ICU type dummy labels (strip "cu_" prefix, shorten)
ICU_SHORT = {
    "Cardiac Vascular Intensive Care Unit (CVICU)": "CVICU",
    "Coronary Care Unit (CCU)":                     "CCU",
    "Medical Intensive Care Unit (MICU)":           "MICU",
    "Medical/Surgical Intensive Care Unit (MICU/SICU)": "MICU/SICU",
    "Med/Surg":                                     "Med/Surg",
    "Medicine/Cardiology Intermediate":             "Med/Cardiol. Intermediate",
    "Medicine":                                     "Medicine",
    "Neuro Intermediate":                           "Neuro Intermediate",
    "Neuro Stepdown":                               "Neuro Stepdown",
    "Neuro Surgical Intensive Care Unit (Neuro SICU)": "Neuro SICU",
    "Neurology":                                    "Neurology",
    "PACU":                                         "PACU",
    "Surgical Intensive Care Unit (SICU)":          "SICU",
    "Surgery/Trauma":                               "Surgery/Trauma",
    "Surgery/Vascular/Intermediate":                "Surgery/Vasc. Intermediate",
    "Trauma SICU (TSICU)":                          "TSICU",
    "Intensive Care Unit (ICU)":                    "ICU (general)",
}

def clean_label(raw):
    if raw.startswith("cu_"):
        suffix = raw[3:]
        return "ICU type: " + ICU_SHORT.get(suffix, suffix)
    return VAR_LABELS.get(raw, raw.replace("_", " ").title())

# ── Domain grouping for love plot ─────────────────────────────────────────────
DOMAINS = {
    "Demographics":    ["female","anchor_age","race_white","race_black"],
    "Kidney":          ["log_baseline_cr","ckd_proxy","ckd","cr_slope_pre","n_cr_pre"],
    "Comorbidities":   ["liver_disease","heart_failure","sepsis","cancer","diabetes","electrolyte_disorder"],
    "Severity":        ["vasopressor_pre","gcs_min_pre","hr_max_pre","spo2_min_pre"],
    "Timing":          ["hours_from_icu_admit","los_hours","year"],
    "ICU Type":        [],  # filled dynamically
}

# ── Helper: fit PS model ──────────────────────────────────────────────────────
def fit_ps(X_, A_):
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale",  StandardScaler()),
        ("clf",    LogisticRegression(max_iter=200, C=1.0, solver="liblinear")),
    ])
    pipe.fit(X_, A_)
    ps = pipe.predict_proba(X_)[:, 1]
    return ps, roc_auc_score(A_, ps)

def smd_fn(x1, x0):
    s = np.sqrt((np.nanvar(x1) + np.nanvar(x0)) / 2)
    return abs(np.nanmean(x1) - np.nanmean(x0)) / s if s > 0 else 0.0

def weighted_smd_fn(col, treat, df_, w):
    w1, w0 = w[treat==1], w[treat==0]
    med = df_[col].median()
    x1 = df_.loc[treat==1, col].fillna(med).values
    x0 = df_.loc[treat==0, col].fillna(med).values
    mu1 = np.average(x1, weights=w1); mu0 = np.average(x0, weights=w0)
    v1  = np.average((x1-mu1)**2, weights=w1)
    v0  = np.average((x0-mu0)**2, weights=w0)
    s   = np.sqrt((v1+v0)/2)
    return abs(mu1-mu0)/s if s > 0 else 0.0


# =============================================================================
# FIGURE 1 — Cohort Flow Diagram (CONSORT-style)
# =============================================================================
def fig_cohort_flow(df):
    print("Figure 1: Cohort flow diagram")
    n_raw   = 44767
    n_all   = len(df)
    n_lr_raw = df.fluid_type.eq("LR").sum()
    n_ns_raw = df.fluid_type.eq("NS").sum()

    eligible = df.query("has_baseline_cr==1 and has_followup_cr==1")
    n_no_cr  = n_all - len(eligible)
    cohort   = eligible.query("baseline_creatinine<=4.0")
    n_prev_aki = len(eligible) - len(cohort)
    keep_ps  = (cohort.index.isin(cohort.index))  # placeholder; actual trim done below
    n_primary = len(cohort)

    # rough trim count from results CSV
    n_trimmed = n_primary - 38266
    n_analytic = 38266
    n_lr = cohort.fluid_type.eq("LR").sum()
    n_ns = cohort.fluid_type.eq("NS").sum()

    fig, ax = plt.subplots(figsize=(9.6, 11))
    ax.set_xlim(-0.6, 11.0); ax.set_ylim(0, 11)
    ax.axis("off")

    def box(x, y, w, h, text, color="#FFFFFF", ec="#37474F", lw=1.5, fontsize=11, bold=False):
        rect = FancyBboxPatch((x-w/2, y-h/2), w, h,
                              boxstyle="round,pad=0.08", linewidth=lw,
                              edgecolor=ec, facecolor=color)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                fontweight="bold" if bold else "normal",
                multialignment="center", color="#1A1A1A")

    def arrow(x, y1, y2, pad_top=0.45, pad_bot=0.43):
        # Draw between box EDGES: start just below the upper box, end just above
        # the lower box, so the arrow never pierces either box.
        ax.annotate("", xy=(x, y2 + pad_bot), xytext=(x, y1 - pad_top),
                    arrowprops=dict(arrowstyle="-|>", color="#37474F",
                                   lw=1.8, mutation_scale=14))

    def excl_box(x, y, text, side="right", feed_x=5.0):
        """Red exclusion box. side='right' (box to the right of the spine, arrow
        points rightward from feed_x) or 'left' (mirrored)."""
        rect = FancyBboxPatch((x, y-0.42), 3.3, 0.84,
                              boxstyle="round,pad=0.06", linewidth=1,
                              edgecolor="#B71C1C", facecolor="#FFEBEE")
        ax.add_patch(rect)
        ax.text(x+1.65, y, text, ha="center", va="center", fontsize=9,
                color="#B71C1C", multialignment="center")
        if side == "right":
            ax.annotate("", xy=(x+0.02, y), xytext=(feed_x, y),
                        arrowprops=dict(arrowstyle="-|>", color="#B71C1C",
                                       lw=1.2, mutation_scale=10))
        else:  # left: box sits left of the spine, arrow points leftward
            ax.annotate("", xy=(x+3.28, y), xytext=(feed_x, y),
                        arrowprops=dict(arrowstyle="-|>", color="#B71C1C",
                                       lw=1.2, mutation_scale=10))

    # Sanity check: the per-arm split below is the ELIGIBLE cohort (n_primary),
    # whose arms sum exactly to n_primary. The PS-trimmed count (n_analytic) has
    # no per-arm breakdown in the source data, so it is shown only as a spine
    # node, never split by arm — this keeps every number traceable to the data.
    assert n_lr + n_ns == n_primary, (
        f"Cohort-flow arm split {n_lr}+{n_ns} != eligible cohort {n_primary}"
    )

    # Row y-positions: raw, Cr-present, eligible, LR/NS branch, PS-trimmed
    ys = [10.3, 8.6, 6.9, 4.9, 2.6]

    box(5, ys[0], 5.5, 0.9,
        f"MIMIC-IV v3.1 — All ICU stays with\nfirst LR or NS within 24 h of ICU admission\n"
        f"n = {n_all:,}  (LR {n_lr_raw:,} · NS {n_ns_raw:,})",
        color="#E3F2FD", ec=LR_COL, lw=2, fontsize=10)

    arrow(5, ys[0], ys[1], pad_top=0.47, pad_bot=0.45)
    excl_box(5.95, (ys[0]+ys[1])/2,
             f"Missing baseline or\nfollow-up creatinine\nn = {n_no_cr:,}",
             side="right", feed_x=5.0)

    box(5, ys[1], 5.5, 0.85,
        f"Baseline AND follow-up creatinine present\nn = {len(eligible):,}",
        color="#F3E5F5", ec="#6A1B9A", fontsize=10)

    arrow(5, ys[1], ys[2], pad_top=0.45, pad_bot=0.45)
    excl_box(5.95, (ys[1]+ys[2])/2,
             f"Prevalent AKI proxy\n(baseline Cr > 4.0 mg/dL)\nn = {n_prev_aki:,}",
             side="right", feed_x=5.0)

    box(5, ys[2], 5.5, 0.85,
        f"Eligible cohort (baseline Cr ≤ 4.0 mg/dL)\nn = {n_primary:,}",
        color="#E8F5E9", ec=GRN, fontsize=10)

    # Branch the eligible cohort into its LR / NS arms (these sum to n_primary).
    ax.annotate("", xy=(2.6, ys[3]+0.42), xytext=(5, ys[2]-0.43),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))
    ax.annotate("", xy=(7.4, ys[3]+0.42), xytext=(5, ys[2]-0.43),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))

    box(2.6, ys[3], 3.6, 0.84,
        f"Lactated Ringer's (LR)\nn = {n_lr:,}  ({100*n_lr/n_primary:.1f}%)",
        color="#E3F2FD", ec=LR_COL, lw=2, fontsize=10, bold=True)
    box(7.4, ys[3], 3.6, 0.84,
        f"Normal Saline (NS)\nn = {n_ns:,}  ({100*n_ns/n_primary:.1f}%)",
        color="#FFF3E0", ec=NS_COL, lw=2, fontsize=10, bold=True)

    # Re-converge both arms into the PS-trimmed primary analytic cohort.
    ax.annotate("", xy=(5, ys[4]+0.43), xytext=(2.6, ys[3]-0.42),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))
    ax.annotate("", xy=(5, ys[4]+0.43), xytext=(7.4, ys[3]-0.42),
                arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.6, mutation_scale=12))
    # PS-trimming exclusion: placed lower-right, below the NS box, feeding from
    # the NS box bottom so it clears the converging arrows (which live at x<7.4).
    excl_y = (ys[3] + ys[4]) / 2 - 0.15
    excl_box(7.7, excl_y,
             f"PS outside [0.05, 0.95]\n(poor overlap)\nn = {n_trimmed:,}",
             side="right", feed_x=7.4)

    box(5, ys[4], 5.8, 0.85,
        f"Primary analytic cohort (PS-trimmed)\nn = {n_analytic:,}",
        color="#FFF8E1", ec="#F57F17", fontsize=10, bold=True)

    ax.set_title("Figure 1. Cohort Flow Diagram", fontsize=14, fontweight="bold",
                 pad=10, loc="left")
    save("fig_cohort_flow.png")


# =============================================================================
# FIGURE 2 — Causal DAG
# =============================================================================
def fig_causal_dag():
    print("Figure 2: Causal DAG")
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set_xlim(0, 13); ax.set_ylim(0, 7)
    ax.axis("off")
    ax.set_facecolor("white")

    def ellipse(cx, cy, w, h, text, fc="#ECEFF1", ec="#546E7A", lw=1.4, fs=9.5):
        e = Ellipse((cx, cy), w, h, facecolor=fc, edgecolor=ec, linewidth=lw)
        ax.add_patch(e)
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                multialignment="center", color="#1A1A1A")

    def rect(cx, cy, w, h, text, fc, ec, lw=2, fs=10.5, bold=True):
        r = FancyBboxPatch((cx-w/2, cy-h/2), w, h,
                           boxstyle="round,pad=0.12", facecolor=fc,
                           edgecolor=ec, linewidth=lw)
        ax.add_patch(r)
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal",
                multialignment="center", color="#1A1A1A")

    def arrow(x1, y1, x2, y2, color="#546E7A", ls="-", lw=1.3, hw=0.18, hl=0.25):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle=f"-|>,head_width={hw},head_length={hl}",
                                   color=color, lw=lw, linestyle=ls,
                                   connectionstyle="arc3,rad=0.0"))

    # Confounder nodes (left column)
    conf = [
        (1.6, 6.3, "Age / Sex / Race"),
        (1.6, 5.2, "ICU Type"),
        (1.6, 4.1, "Baseline Creatinine"),
        (1.6, 3.0, "CKD / Comorbidities"),
        (1.6, 1.9, "Illness Severity\n(GCS / HR / SpO₂)"),
        (1.6, 0.9, "Pre-ICU Fluids /\nCalendar Period"),
    ]
    for cx, cy, txt in conf:
        ellipse(cx, cy, 2.7, 0.75, txt, fc="#ECEFF1", ec="#546E7A")

    # Exposure
    rect(6.0, 3.5, 2.8, 1.0, "Initial ICU Fluid\n(LR vs NS)",
         fc="#BBDEFB", ec=LR_COL, lw=2.5)

    # Outcome
    rect(10.8, 3.5, 2.8, 1.0, "7-day\nIncident AKI",
         fc="#FFCDD2", ec="#C62828", lw=2.5)

    # Mediator
    rect(6.0, 1.3, 2.6, 0.8, "Post-treatment\nFluid Management\n(not adjusted)",
         fc="#FFF9C4", ec="#F9A825", lw=1.5, fs=9, bold=False)

    # Confounder → exposure arrows
    for cx, cy, _ in conf:
        arrow(cx+1.35, cy, 5.6, 3.5+np.sign(cy-3.5)*0.4, color="#78909C", lw=1.1)

    # Confounder → outcome arrows
    for cx, cy, _ in conf:
        arrow(cx+1.35, cy, 10.4, 3.5+np.sign(cy-3.5)*0.4, color="#78909C", lw=1.1)

    # Exposure → outcome (thick, main effect)
    arrow(6.9, 3.5, 10.3, 3.5, color=LR_COL, lw=2.5, hw=0.22, hl=0.3)

    # Exposure → mediator (dashed)
    arrow(6.0, 2.98, 6.0, 1.72, color="#F9A825", ls="dashed", lw=1.5, hw=0.15, hl=0.2)

    # Mediator → outcome (dashed)
    arrow(7.3, 1.3, 10.3, 3.1, color="#F9A825", ls="dashed", lw=1.5, hw=0.15, hl=0.2)

    # Legend
    handles = [
        mpatches.Patch(fc="#ECEFF1", ec="#546E7A", label="Confounder (adjusted)"),
        mpatches.Patch(fc="#BBDEFB", ec=LR_COL,   label="Exposure"),
        mpatches.Patch(fc="#FFCDD2", ec="#C62828", label="Outcome"),
        mpatches.Patch(fc="#FFF9C4", ec="#F9A825", label="Mediator (not adjusted)"),
        mlines.Line2D([], [], color=LR_COL, lw=2.5, label="Causal effect of interest"),
        mlines.Line2D([], [], color="#F9A825", lw=1.5, ls="--", label="Mediated path (excluded)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=9,
              framealpha=0.95, edgecolor="#BDBDBD", bbox_to_anchor=(1.0, 0.0))

    ax.set_title("Figure 2. Proposed Causal DAG: Initial ICU Fluid and 7-day AKI",
                 fontsize=13, fontweight="bold", loc="left", pad=8)
    save("fig_causal_dag.png")


# =============================================================================
# FIGURE 3 — PS Overlap
# =============================================================================
def fig_ps_overlap(cohort_full, ps_full, cohort_trim, ps_trim):
    """cohort_full/ps_full: all eligible; cohort_trim/ps_trim: PS-trimmed subset."""
    print("Figure 3: PS overlap")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=False)

    for fluid, color, label in [("LR", LR_COL, "LR"), ("NS", NS_COL, "NS")]:
        mask_full = cohort_full["fluid_type"].eq(fluid).values
        mask_trim = cohort_trim["fluid_type"].eq(fluid).values
        axes[0].hist(ps_full[mask_full], bins=70, alpha=0.55,
                     label=label, color=color, density=True, linewidth=0)
        axes[1].hist(ps_trim[mask_trim], bins=70, alpha=0.55,
                     label=label, color=color, density=True, linewidth=0)

    n_full = len(cohort_full)
    n_trim = len(cohort_trim)
    for ax, title in zip(axes, [f"A   Full cohort (n={n_full:,})",
                                  f"B   PS-trimmed (n={n_trim:,})"]):
        ax.axvline(0.05, color="#B71C1C", ls="--", lw=1.5, label="Trim 0.05 / 0.95")
        ax.axvline(0.95, color="#B71C1C", ls="--", lw=1.5)
        ax.set_xlabel("Propensity Score (P[LR])", fontsize=11)
        ax.set_ylabel("Density", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold", loc="left", pad=4)
        ax.legend(framealpha=0.9)
        ax.set_xlim(0, 1)

    axes[0].text(0.5, 0.97,
                 f"AUC = {roc_auc_score(cohort_full['treat'].values, ps_full):.3f}",
                 transform=axes[0].transAxes, ha="center", va="top",
                 fontsize=10, color="#37474F",
                 bbox=dict(fc="white", ec="#BDBDBD", pad=3, alpha=0.9))

    fig.suptitle("Figure 3. Propensity Score Distribution by Fluid Group",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.01)
    save("fig_ps_overlap.png")


# =============================================================================
# FIGURE 4 — Love Plot (SMD)
# =============================================================================
def fig_love_plot(cov_cols, smd_before, smd_after):
    print("Figure 4: Love plot")

    # Build ordered list by domain
    domain_order = ["Demographics", "Kidney", "Comorbidities", "Severity", "Timing", "ICU Type"]
    icu_cols = [c for c in cov_cols if c.startswith("cu_")]
    DOMAINS["ICU Type"] = icu_cols

    ordered = []
    for domain in domain_order:
        members = [c for c in DOMAINS.get(domain, []) if c in cov_cols]
        if domain == "ICU Type":
            members = [c for c in cov_cols if c.startswith("cu_")]
        for c in members:
            if c in smd_before:
                ordered.append((domain, c))
    # catch any leftover
    already = {c for _, c in ordered}
    for c in cov_cols:
        if c in smd_before and c not in already:
            ordered.append(("Other", c))

    n = len(ordered)
    fig_h = max(7, n * 0.32 + 2)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    # Domain color bands
    domain_colors = {
        "Demographics": "#E3F2FD", "Kidney": "#E8F5E9",
        "Comorbidities": "#FFF3E0", "Severity": "#FCE4EC",
        "Timing": "#F3E5F5", "ICU Type": "#ECEFF1", "Other": "#FAFAFA",
    }
    prev_domain = None
    band_start = 0
    band_map = {}
    for i, (domain, c) in enumerate(ordered):
        if domain != prev_domain:
            if prev_domain is not None:
                band_map[prev_domain] = (band_start, i)
            band_start = i
            prev_domain = domain
    if prev_domain:
        band_map[prev_domain] = (band_start, len(ordered))

    for domain, (i0, i1) in band_map.items():
        ax.axhspan(i0 - 0.5, i1 - 0.5, facecolor=domain_colors.get(domain, "#FAFAFA"),
                   alpha=0.4, zorder=0)
        ax.text(-0.005, (i0 + i1) / 2, domain, va="center", ha="right",
                fontsize=8.5, color="#546E7A", fontstyle="italic",
                transform=ax.get_yaxis_transform())

    ys = np.arange(n)
    before_vals = [smd_before[c] for _, c in ordered]
    after_vals  = [smd_after[c]  for _, c in ordered]

    ax.scatter(before_vals, ys, color="#E53935", marker="o", s=38, zorder=4,
               label="Before IPTW", linewidths=0)
    ax.scatter(after_vals,  ys, color=LR_COL,   marker="D", s=30, zorder=5,
               label="After IPTW", linewidths=0)

    # Connect before→after with thin line
    for i, (b, a) in enumerate(zip(before_vals, after_vals)):
        ax.plot([b, a], [i, i], color="#BDBDBD", lw=0.6, zorder=2)

    ax.axvline(0.1, color="#37474F", ls="--", lw=1.3, label="|SMD| = 0.1", zorder=3)
    ax.axvline(0.0, color="#9E9E9E", lw=0.8, zorder=1)

    labels = [clean_label(c) for _, c in ordered]
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Absolute Standardized Mean Difference (SMD)", fontsize=11)
    ax.set_xlim(-0.02, max(before_vals) * 1.08)
    ax.set_ylim(-0.8, n - 0.2)
    ax.invert_yaxis()

    pct_b = 100 * np.mean([v < 0.1 for v in smd_before.values()])
    pct_a = 100 * np.mean([v < 0.1 for v in smd_after.values()])
    ax.text(0.98, 0.01,
            f"|SMD| < 0.10: {pct_b:.0f}% → {pct_a:.0f}% after IPTW",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            color="#37474F", bbox=dict(fc="white", ec="#BDBDBD", pad=3, alpha=0.9))

    ax.legend(loc="lower right", framealpha=0.9, fontsize=9)
    ax.grid(axis="x", color="#E0E0E0", linewidth=0.6)
    ax.set_axisbelow(True)

    fig.suptitle("Figure 4. Covariate Balance Before and After IPTW",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.005)
    fig.tight_layout()
    save("fig_love_plot.png")


# =============================================================================
# FIGURE 5 — AIPW Forest Plot (Primary + SA1–SA4)
# =============================================================================
def fig_aipw_forest(res_csv):
    print("Figure 5: AIPW forest plot")
    res = pd.read_csv(res_csv)
    # Only primary + SA1–SA4
    keep = ["Primary AIPW (PS-trimmed)", "SA1: No PS trimming",
            "SA2: Stricter baseline Cr≤2.0", "SA3: Composite death-or-AKI",
            "SA4: IPTW-only (no outcome model)"]
    res = res[res["analysis"].isin(keep)].copy()
    # Canonical order: primary first, then SA1–4 bottom-up for readability
    order = keep[::-1]
    res["_ord"] = res["analysis"].map({v: i for i, v in enumerate(order)})
    res = res.sort_values("_ord")

    # Short labels
    short = {
        "Primary AIPW (PS-trimmed)":        "Primary AIPW\n(PS-trimmed, enriched covariates)",
        "SA1: No PS trimming":               "SA1  No PS trimming",
        "SA2: Stricter baseline Cr≤2.0":     "SA2  Stricter baseline Cr ≤ 2.0 mg/dL",
        "SA3: Composite death-or-AKI":       "SA3  Composite: death OR AKI",
        "SA4: IPTW-only (no outcome model)": "SA4  IPTW-only (no outcome model)",
    }

    fig = plt.figure(figsize=(13, 5))
    # 3-panel layout: label | forest | table
    gs = GridSpec(1, 3, figure=fig, width_ratios=[3.5, 4, 3], wspace=0.04)
    ax_lbl  = fig.add_subplot(gs[0])
    ax_frst = fig.add_subplot(gs[1])
    ax_tbl  = fig.add_subplot(gs[2])

    n = len(res)
    ys = np.arange(n)
    is_primary = res["analysis"].eq("Primary AIPW (PS-trimmed)").values
    is_sig = ((res["rd_lower"] > 0) | (res["rd_upper"] < 0)).values

    # --- Forest panel ---
    for i, row in enumerate(res.itertuples()):
        sig = (row.rd_lower > 0 or row.rd_upper < 0)
        col = RED_SIG if sig else ("#1565C0" if row.analysis.startswith("Primary") else GRY)
        lw  = 2.5 if row.analysis.startswith("Primary") else 1.5
        ms  = 10  if row.analysis.startswith("Primary") else 7

        ax_frst.plot([row.rd_lower, row.rd_upper], [i, i], color=col, lw=lw, solid_capstyle="round")
        ax_frst.scatter([row.rd], [i], color=col, s=ms**2, zorder=5,
                        marker="D" if row.analysis.startswith("Primary") else "o")

    ax_frst.axvline(0, color="#37474F", lw=1.5, ls="--", zorder=2)
    ax_frst.axhline(n - 1.5, color="#BDBDBD", lw=0.8, ls="-")  # separator after primary

    ax_frst.set_xlim(-0.065, 0.065)
    ax_frst.set_yticks(ys); ax_frst.set_yticklabels([])
    ax_frst.set_xlabel("Risk Difference  LR − NS  (percentage points × 1)", fontsize=10)
    ax_frst.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x*100:+.1f}"))
    ax_frst.set_ylim(-0.8, n - 0.2)
    ax_frst.invert_yaxis()
    ax_frst.grid(axis="x", color="#E0E0E0", lw=0.6)
    ax_frst.set_axisbelow(True)
    ax_frst.spines["left"].set_visible(False)
    ax_frst.set_title("Risk Difference (95% CI)", fontsize=10, loc="center", pad=4)

    # Annotation: "LR higher risk" / "NS higher risk"
    ax_frst.text(0.001, -0.6, "LR higher →", fontsize=8, color="#B71C1C", va="center")
    ax_frst.text(-0.001, -0.6, "← NS higher", fontsize=8, color=GRN, va="center", ha="right")

    # --- Label panel ---
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(-0.8, n - 0.2)
    ax_lbl.axis("off"); ax_lbl.invert_yaxis()
    ax_lbl.set_title("Analysis", fontsize=10, loc="right", pad=4)
    for i, row in enumerate(res.itertuples()):
        bold = row.analysis.startswith("Primary")
        ax_lbl.text(0.98, i, short[row.analysis], va="center", ha="right",
                    fontsize=9.5, fontweight="bold" if bold else "normal")
    ax_lbl.axhline(n - 1.5, color="#BDBDBD", lw=0.8)

    # --- Table panel ---
    ax_tbl.set_xlim(0, 1); ax_tbl.set_ylim(-0.8, n - 0.2)
    ax_tbl.axis("off"); ax_tbl.invert_yaxis()
    # Header
    ax_tbl.text(0.05, -0.6, "n",        fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.35, -0.6, "LR risk",  fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.62, -0.6, "NS risk",  fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.90, -0.6, "RD (95% CI)", fontsize=9, ha="center", fontweight="bold")
    ax_tbl.axhline(n - 1.5, color="#BDBDBD", lw=0.8)
    for i, row in enumerate(res.itertuples()):
        sig = (pd.notna(row.rd_lower) and (row.rd_lower > 0 or row.rd_upper < 0))
        col = RED_SIG if sig else "#1A1A1A"
        ax_tbl.text(0.05, i, f"{row.n:,}",           fontsize=9, ha="center", va="center", color=col)
        ax_tbl.text(0.35, i, f"{row.risk_lr*100:.1f}%", fontsize=9, ha="center", va="center", color=LR_COL)
        ax_tbl.text(0.62, i, f"{row.risk_ns*100:.1f}%", fontsize=9, ha="center", va="center", color=NS_COL)
        ci_txt = (f"{row.rd*100:+.2f} ({row.rd_lower*100:+.2f}, {row.rd_upper*100:+.2f})"
                  if pd.notna(row.rd_lower) else f"{row.rd*100:+.2f}")
        ax_tbl.text(0.90, i, ci_txt, fontsize=8.5, ha="center", va="center",
                    color=col, fontweight="bold" if sig else "normal")
    ax_tbl.set_title("n / Risk (%) / RD pp (95% CI)", fontsize=10, loc="center", pad=4)

    fig.suptitle("Figure 5. AIPW Risk Differences: Primary and Sensitivity Analyses (enriched covariates)",
                 fontsize=12, fontweight="bold", x=0.02, ha="left", y=1.02)
    save("fig_aipw_results.png")


# =============================================================================
# FIGURE 6 — Subgroup Forest (SA7 sex + SA8 ICU type)
# =============================================================================
def fig_subgroup_forest(res_csv):
    print("Figure 6: Subgroup forest plot")
    res = pd.read_csv(res_csv)
    keep_pat = ["SA7:", "SA8:"]
    sub = res[res["analysis"].str.startswith(tuple(keep_pat))].copy()

    short_map = {
        "SA7: Female subgroup":                                              "Female",
        "SA7: Male subgroup":                                                "Male",
        "SA8: Cardiac Vascular Intensive Care Unit (CVICU)":                "CVICU",
        "SA8: Medical Intensive Care Unit (MICU)":                          "MICU",
        "SA8: Medical/Surgical Intensive Care Unit (MICU/SICU)":            "MICU/SICU",
    }
    group_map = {k: ("Sex" if "SA7" in k else "ICU Type") for k in short_map}

    sub["label"] = sub["analysis"].map(short_map)
    sub["group"] = sub["analysis"].map(group_map)
    sub = sub.dropna(subset=["label"])

    # Order: ICU type first (most interesting), then sex
    order = ["MICU","CVICU","MICU/SICU","Female","Male"]
    sub["_ord"] = sub["label"].map({v:i for i,v in enumerate(order)})
    sub = sub.sort_values("_ord")

    fig = plt.figure(figsize=(13, 4.5))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[2.5, 4.5, 3], wspace=0.04)
    ax_lbl  = fig.add_subplot(gs[0])
    ax_frst = fig.add_subplot(gs[1])
    ax_tbl  = fig.add_subplot(gs[2])

    n = len(sub)
    ys = np.arange(n)

    # Color: ICU type = orange family, sex = blue family
    group_cols = {"MICU": RED_SIG, "CVICU": "#E65100", "MICU/SICU": GRY,
                  "Female": "#6A1B9A", "Male": "#1565C0"}

    for i, row in enumerate(sub.itertuples()):
        col = group_cols.get(row.label, GRY)
        sig = pd.notna(row.rd_lower) and (row.rd_lower > 0 or row.rd_upper < 0)
        lw  = 2.2 if sig else 1.4
        ax_frst.plot([row.rd_lower, row.rd_upper], [i, i], color=col, lw=lw, solid_capstyle="round")
        ax_frst.scatter([row.rd], [i], color=col, s=64, zorder=5)
        if sig:
            ax_frst.scatter([row.rd], [i], color=col, s=200, marker="*",
                            zorder=6, alpha=0.3)

    # Reference line = primary result
    ax_frst.axvline(0, color="#37474F", lw=1.5, ls="--")
    # Primary AIPW estimate
    primary_rd = -0.00134
    ax_frst.axvline(primary_rd, color="#1565C0", lw=1.0, ls=":", alpha=0.6,
                    label=f"Primary AIPW ({primary_rd*100:+.2f} pp)")

    # Group separator
    n_icu = sub["group"].eq("ICU Type").sum()
    ax_frst.axhline(n_icu - 0.5, color="#BDBDBD", lw=1.0)
    ax_frst.text(0.072, n_icu - 0.5, "Sex", fontsize=9, color="#546E7A",
                 va="center", style="italic")
    ax_frst.text(0.072, n_icu/2 - 0.5, "ICU Type", fontsize=9, color="#546E7A",
                 va="center", style="italic")

    ax_frst.set_xlim(-0.095, 0.095)
    ax_frst.set_yticks(ys); ax_frst.set_yticklabels([])
    ax_frst.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x*100:+.0f}"))
    ax_frst.set_xlabel("Risk Difference  LR − NS  (percentage points)", fontsize=10)
    ax_frst.set_ylim(-0.8, n - 0.2); ax_frst.invert_yaxis()
    ax_frst.grid(axis="x", color="#E0E0E0", lw=0.6); ax_frst.set_axisbelow(True)
    ax_frst.spines["left"].set_visible(False)
    ax_frst.set_title("Risk Difference (95% CI)", fontsize=10, loc="center", pad=4)
    ax_frst.text(0.001, -0.6, "LR higher →", fontsize=8, color="#B71C1C", va="center")
    ax_frst.text(-0.001, -0.6, "← LR protective", fontsize=8, color=GRN, va="center", ha="right")
    ax_frst.legend(fontsize=8, loc="lower right", framealpha=0.9)

    # --- Labels ---
    ax_lbl.set_xlim(0,1); ax_lbl.set_ylim(-0.8, n-0.2); ax_lbl.axis("off"); ax_lbl.invert_yaxis()
    ax_lbl.set_title("Subgroup", fontsize=10, loc="right", pad=4)
    ax_lbl.axhline(n_icu-0.5, color="#BDBDBD", lw=1.0)
    for i, row in enumerate(sub.itertuples()):
        ax_lbl.text(0.98, i, row.label, va="center", ha="right", fontsize=10,
                    color=group_cols.get(row.label, GRY), fontweight="bold")

    # --- Table ---
    ax_tbl.set_xlim(0,1); ax_tbl.set_ylim(-0.8, n-0.2); ax_tbl.axis("off"); ax_tbl.invert_yaxis()
    ax_tbl.axhline(n_icu-0.5, color="#BDBDBD", lw=1.0)
    ax_tbl.text(0.05, -0.6, "n",         fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.33, -0.6, "LR risk",   fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.58, -0.6, "NS risk",   fontsize=9, ha="center", fontweight="bold")
    ax_tbl.text(0.90, -0.6, "RD pp (95% CI)", fontsize=9, ha="center", fontweight="bold")
    ax_tbl.set_title("n / Risk / RD pp (95% CI)", fontsize=10, loc="center", pad=4)
    for i, row in enumerate(sub.itertuples()):
        sig = pd.notna(row.rd_lower) and (row.rd_lower > 0 or row.rd_upper < 0)
        col = group_cols.get(row.label, GRY)
        ax_tbl.text(0.05, i, f"{row.n:,}",             fontsize=9, ha="center", va="center")
        ax_tbl.text(0.33, i, f"{row.risk_lr*100:.1f}%", fontsize=9, ha="center", va="center", color=LR_COL)
        ax_tbl.text(0.58, i, f"{row.risk_ns*100:.1f}%", fontsize=9, ha="center", va="center", color=NS_COL)
        ci = (f"{row.rd*100:+.1f} ({row.rd_lower*100:+.1f}, {row.rd_upper*100:+.1f})"
              if pd.notna(row.rd_lower) else f"{row.rd*100:+.1f}")
        ax_tbl.text(0.90, i, ci, fontsize=8.5, ha="center", va="center",
                    color=col, fontweight="bold" if sig else "normal")

    fig.suptitle("Figure 6. Subgroup Analyses — Sex (SA7) and ICU Type (SA8)",
                 fontsize=12, fontweight="bold", x=0.02, ha="left", y=1.02)
    save("fig_subgroup_forest.png")


# =============================================================================
# FIGURE 7+8 — PS Matching
# =============================================================================
def fig_matching(cohort, ps_full):
    print("Figure 7: PS matching overlap")

    from sklearn.neighbors import NearestNeighbors

    treat = cohort["treat"].values
    logit_ps = np.log(np.clip(ps_full, 1e-6, 1-1e-6) /
                      np.clip(1-ps_full, 1e-6, 1-1e-6))
    caliper = 0.2 * logit_ps.std()

    lr_idx = np.where(treat == 1)[0]
    ns_idx = np.where(treat == 0)[0]
    lr_logit = logit_ps[lr_idx].reshape(-1, 1)
    ns_logit = logit_ps[ns_idx].reshape(-1, 1)

    nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
    nn.fit(ns_logit)
    dists, matches = nn.kneighbors(lr_logit)
    within = dists.flatten() <= caliper

    matched_lr = lr_idx[within]
    matched_ns = ns_idx[matches.flatten()[within]]
    n_pairs = len(matched_lr)

    ps_lr_m = ps_full[matched_lr]
    ps_ns_m = ps_full[matched_ns]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Before matching
    for fluid, col, ps_arr in [("LR", LR_COL, ps_full[treat==1]),
                                ("NS", NS_COL, ps_full[treat==0])]:
        axes[0].hist(ps_arr, bins=70, alpha=0.5, color=col, label=fluid, density=True, linewidth=0)
    axes[0].set_title("A   Before matching", fontsize=12, fontweight="bold", loc="left")
    axes[0].set_xlabel("Propensity Score"); axes[0].set_ylabel("Density")
    axes[0].set_xlim(0, 1); axes[0].legend()
    axes[0].text(0.98, 0.97, f"All  LR n={treat.sum():,}  NS n={(1-treat).sum():,}",
                 transform=axes[0].transAxes, ha="right", va="top", fontsize=9,
                 bbox=dict(fc="white", ec="#BDBDBD", pad=2, alpha=0.9))

    # After matching
    axes[1].hist(ps_lr_m, bins=50, alpha=0.5, color=LR_COL, label="LR (matched)", density=True, linewidth=0)
    axes[1].hist(ps_ns_m, bins=50, alpha=0.5, color=NS_COL, label="NS (matched)", density=True, linewidth=0)
    axes[1].set_title("B   After 1:1 matching", fontsize=12, fontweight="bold", loc="left")
    axes[1].set_xlabel("Propensity Score"); axes[1].set_ylabel("Density")
    axes[1].set_xlim(0, 1); axes[1].legend()
    axes[1].text(0.98, 0.97, f"Matched pairs n={n_pairs:,}\ncaliper = {caliper:.3f} (0.2×SD logit)",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=9,
                 bbox=dict(fc="white", ec="#BDBDBD", pad=2, alpha=0.9))

    fig.suptitle("Figure 7. Propensity Score Distribution Before and After 1:1 Nearest-Neighbour Matching",
                 fontsize=12, fontweight="bold", x=0.02, ha="left", y=1.01)
    save("fig_matching_ps.png")
    return matched_lr, matched_ns


def fig_matching_balance(cohort, cov_cols, ps_full, matched_lr, matched_ns):
    print("Figure 8: Matching balance")
    treat = cohort["treat"].values
    cov_df_full = pd.DataFrame(
        cohort[cov_cols].values, columns=cov_cols)

    matched_mask = np.zeros(len(cohort), dtype=bool)
    matched_mask[matched_lr] = True
    matched_mask[matched_ns] = True
    cov_df_match = cov_df_full[matched_mask].reset_index(drop=True)
    treat_match  = treat[matched_mask]

    # SMD unmatched and matched
    treat_s_full  = pd.Series(treat)
    treat_s_match = pd.Series(treat_match)

    smd_unmatched = {c: smd_fn(cov_df_full.loc[treat_s_full==1, c].dropna().values,
                                cov_df_full.loc[treat_s_full==0, c].dropna().values)
                     for c in cov_cols}
    smd_matched   = {c: smd_fn(cov_df_match.loc[treat_s_match==1, c].dropna().values,
                                cov_df_match.loc[treat_s_match==0, c].dropna().values)
                     for c in cov_cols}

    # Reuse love-plot layout but for matching
    domain_order = ["Demographics","Kidney","Comorbidities","Severity","Timing","ICU Type"]
    DOMAINS["ICU Type"] = [c for c in cov_cols if c.startswith("cu_")]
    ordered = []
    for domain in domain_order:
        members = [c for c in DOMAINS.get(domain, []) if c in smd_unmatched]
        if domain == "ICU Type":
            members = [c for c in cov_cols if c.startswith("cu_")]
        for c in members:
            if c in smd_unmatched:
                ordered.append(c)
    for c in cov_cols:
        if c in smd_unmatched and c not in ordered:
            ordered.append(c)

    n = len(ordered)
    fig_h = max(7, n * 0.32 + 2)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    bv = [smd_unmatched[c] for c in ordered]
    av = [smd_matched[c]   for c in ordered]
    ys = np.arange(n)

    for i, (b, a) in enumerate(zip(bv, av)):
        ax.plot([b, a], [i, i], color="#BDBDBD", lw=0.6, zorder=2)
    ax.scatter(bv, ys, color="#E53935", marker="o", s=38, zorder=4, label="Before matching")
    ax.scatter(av, ys, color=GRN,       marker="D", s=30, zorder=5, label="After 1:1 matching")
    ax.axvline(0.1, color="#37474F", ls="--", lw=1.3, label="|SMD| = 0.1", zorder=3)
    ax.axvline(0.0, color="#9E9E9E", lw=0.8)

    labels = [clean_label(c) for c in ordered]
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Absolute Standardized Mean Difference (SMD)", fontsize=11)
    ax.set_xlim(-0.02, max(bv) * 1.1)
    ax.set_ylim(-0.8, n - 0.2); ax.invert_yaxis()

    pct_b = 100 * np.mean([v < 0.1 for v in smd_unmatched.values()])
    pct_a = 100 * np.mean([v < 0.1 for v in smd_matched.values()])
    ax.text(0.98, 0.01, f"|SMD| < 0.10: {pct_b:.0f}% → {pct_a:.0f}% after matching",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(fc="white", ec="#BDBDBD", pad=3, alpha=0.9))
    ax.legend(loc="lower right", framealpha=0.9, fontsize=9)
    ax.grid(axis="x", color="#E0E0E0", lw=0.6); ax.set_axisbelow(True)

    fig.suptitle("Figure 8. Covariate Balance Before and After 1:1 PS Matching",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.005)
    fig.tight_layout()
    save("fig_matching_balance.png")


# =============================================================================
# FIGURE 9 — Competing Events
# =============================================================================
def fig_competing_events(cohort):
    print("Figure 9: Competing events")

    df = cohort.copy()
    df["days_to_death"] = (
        (pd.to_datetime(df["deathtime"]) - pd.to_datetime(df["fluid_starttime"]))
        .dt.total_seconds() / 86400
    )
    df["days_to_disch"] = (
        (pd.to_datetime(df["dischtime"]) - pd.to_datetime(df["fluid_starttime"]))
        .dt.total_seconds() / 86400
    )
    df["died_7d"]     = (df["days_to_death"] >= 0) & (df["days_to_death"] <= 7)
    df["disch_7d"]    = (~df["died_7d"]) & (df["days_to_disch"] >= 0) & (df["days_to_disch"] <= 7)
    df["missing_cr"]  = df["has_followup_cr"] == 0

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=False)

    groups = ["LR", "NS"]
    colors_ev = {
        "7-day AKI": "#E53935",
        "Died ≤7d":  "#37474F",
        "Disch. alive ≤7d": GRN,
        "Missing follow-up Cr": "#9E9E9E",
        "No event (obs. ≥7d)": "#B0BEC5",
    }

    # Panel A: stacked 100% bar with % labels
    ax = axes[0]
    bottoms = {g: 0.0 for g in groups}
    xs = np.arange(len(groups))
    bar_w = 0.45

    event_data = {}
    for g in groups:
        sub = df[df["fluid_type"] == g]
        n   = len(sub)
        aki  = sub["aki_flag"].mean() if "aki_flag" in sub else 0
        died = sub["died_7d"].mean()
        disc = sub["disch_7d"].mean()
        miss = sub["missing_cr"].mean()
        rest = max(0, 1 - aki - died - disc - miss)
        event_data[g] = {
            "7-day AKI": aki,
            "Died ≤7d":  died,
            "Disch. alive ≤7d": disc,
            "Missing follow-up Cr": miss,
            "No event (obs. ≥7d)": rest,
        }

    for ev, col in colors_ev.items():
        vals = [event_data[g][ev] for g in groups]
        bots = [bottoms[g] for g in groups]
        bars = ax.bar(xs, vals, bottom=bots, color=col, label=ev, width=bar_w,
                      edgecolor="white", linewidth=0.5)
        for j, (bar, v, bot) in enumerate(zip(bars, vals, bots)):
            if v > 0.025:
                ax.text(bar.get_x() + bar.get_width()/2, bot + v/2,
                        f"{v*100:.1f}%", ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")
        for g, v in zip(groups, vals):
            bottoms[g] += v

    ax.set_xticks(xs)
    ax.set_xticklabels([f"LR\n(n={df[df.fluid_type=='LR'].shape[0]:,})",
                         f"NS\n(n={df[df.fluid_type=='NS'].shape[0]:,})"], fontsize=11)
    ax.set_ylabel("Proportion of patients", fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
    ax.set_title("A   7-day outcomes (all eligible)", fontsize=12, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax.grid(axis="y", color="#E0E0E0", lw=0.6); ax.set_axisbelow(True)

    # Panel B: grouped bars for key events only (AKI + death)
    ax2 = axes[1]
    ev_keys = ["7-day AKI", "Died ≤7d", "Missing follow-up Cr"]
    ev_cols = ["#E53935", "#37474F", "#9E9E9E"]
    x2 = np.arange(len(ev_keys))
    bw = 0.3
    for i, g in enumerate(["LR", "NS"]):
        col = LR_COL if g == "LR" else NS_COL
        vals = [event_data[g][e] for e in ev_keys]
        offset = -bw/2 if g == "LR" else bw/2
        bars = ax2.bar(x2 + offset, [v*100 for v in vals], width=bw,
                       color=col, label=g, edgecolor="white", linewidth=0.5, alpha=0.85)
        for bar, v in zip(bars, vals):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.4,
                     f"{v*100:.1f}%", ha="center", va="bottom", fontsize=9, color=col,
                     fontweight="bold")

    ax2.set_xticks(x2); ax2.set_xticklabels(ev_keys, fontsize=10)
    ax2.set_ylabel("Percentage (%)", fontsize=11)
    ax2.set_title("B   Key events by group (grouped comparison)", fontsize=12,
                  fontweight="bold", loc="left")
    ax2.legend(fontsize=10, framealpha=0.9)
    ax2.grid(axis="y", color="#E0E0E0", lw=0.6); ax2.set_axisbelow(True)
    ax2.set_ylim(0, max(event_data["LR"]["7-day AKI"],
                        event_data["NS"]["7-day AKI"]) * 100 * 1.35)

    fig.suptitle("Figure 9. 7-Day Competing Events by Fluid Group (Primary Analytic Cohort)",
                 fontsize=12, fontweight="bold", x=0.02, ha="left", y=1.01)
    plt.tight_layout()
    save("fig_competing_events.png")


# =============================================================================
# TABLE 1 — Baseline Characteristics (saved as PNG)
# =============================================================================
def fig_table1(cohort):
    print("Table 1: Baseline characteristics")

    df = cohort.copy()
    lr = df[df["fluid_type"] == "LR"]
    ns = df[df["fluid_type"] == "NS"]
    N_lr, N_ns = len(lr), len(ns)

    def pct(sub, col):
        return f"{sub[col].mean()*100:.1f}%  ({sub[col].sum():.0f})"

    def med_iqr(sub, col):
        q1, med, q3 = sub[col].quantile([0.25, 0.5, 0.75])
        return f"{med:.1f}  [{q1:.1f}–{q3:.1f}]"

    def smd_val(col, binary=False):
        x1 = lr[col].dropna().values
        x0 = ns[col].dropna().values
        return smd_fn(x1, x0)

    rows = []

    def add(label, lr_val, ns_val, smd_v, indent=False):
        rows.append(("  " + label if indent else label, lr_val, ns_val, f"{smd_v:.2f}"))

    rows.append(("", f"LR  (n={N_lr:,})", f"NS  (n={N_ns:,})", "|SMD|"))

    add("Age, median [IQR]",   med_iqr(lr,"anchor_age"),   med_iqr(ns,"anchor_age"),   smd_val("anchor_age"))
    add("Female sex, % (n)",   pct(lr,"female"),            pct(ns,"female"),            smd_val("female"))
    add("Race: White, % (n)",  pct(lr,"race_white"),        pct(ns,"race_white"),         smd_val("race_white"))
    add("Race: Black, % (n)",  pct(lr,"race_black"),        pct(ns,"race_black"),         smd_val("race_black"))
    rows.append(("", "", "", ""))
    add("Baseline creatinine (mg/dL), median [IQR]",
        med_iqr(lr,"baseline_creatinine"), med_iqr(ns,"baseline_creatinine"),
        smd_val("baseline_creatinine"))
    add("CKD proxy (Cr >1.5), % (n)", pct(lr,"ckd_proxy"), pct(ns,"ckd_proxy"), smd_val("ckd_proxy"))

    if "ckd" in df.columns:
        rows.append(("", "", "", ""))
        for col, lbl in [("ckd","CKD (ICD-coded)"),("liver_disease","Liver disease"),
                          ("heart_failure","Heart failure"),("sepsis","Sepsis"),
                          ("cancer","Cancer"),("diabetes","Diabetes"),
                          ("electrolyte_disorder","Electrolyte disorder")]:
            if col in df.columns:
                add(f"{lbl}, % (n)", pct(lr,col), pct(ns,col), smd_val(col), indent=True)

    if "vasopressor_pre" in df.columns:
        rows.append(("", "", "", ""))
        add("Vasopressor pre-t₀, % (n)", pct(lr,"vasopressor_pre"), pct(ns,"vasopressor_pre"), smd_val("vasopressor_pre"))
    if "gcs_min_pre" in df.columns:
        add("GCS min (24h pre-t₀), median [IQR]", med_iqr(lr,"gcs_min_pre"), med_iqr(ns,"gcs_min_pre"), smd_val("gcs_min_pre"))
    if "hr_max_pre" in df.columns:
        add("HR max (24h pre-t₀), median [IQR]", med_iqr(lr,"hr_max_pre"), med_iqr(ns,"hr_max_pre"), smd_val("hr_max_pre"))

    rows.append(("", "", "", ""))
    add("Hours to first fluid, median [IQR]",
        med_iqr(lr,"hours_from_icu_admit"), med_iqr(ns,"hours_from_icu_admit"),
        smd_val("hours_from_icu_admit"))
    add("ICU LOS at admission (h), median [IQR]",
        med_iqr(lr,"los_hours"), med_iqr(ns,"los_hours"), smd_val("los_hours"))

    rows.append(("", "", "", ""))
    rows.append(("7-day AKI (crude)", "", "", ""))
    if "aki_flag" in df.columns:
        add("AKI event, % (n)", pct(lr,"aki_flag"), pct(ns,"aki_flag"), smd_val("aki_flag"), indent=True)

    # Draw table
    n_rows = len(rows)
    fig_h = max(8, n_rows * 0.38 + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax.axis("off")

    col_x = [0.00, 0.44, 0.65, 0.86, 0.97]
    col_w = [0.44, 0.21, 0.21, 0.11]
    row_h = 1.0 / (n_rows + 2)

    # Header background
    ax.add_patch(plt.Rectangle((0, 1 - row_h), 1, row_h,
                                transform=ax.transAxes, fc="#1565C0", zorder=0))
    for j, (txt, cx) in enumerate(zip(["Characteristic",
                                        f"LR  (n={N_lr:,})",
                                        f"NS  (n={N_ns:,})",
                                        "|SMD|"], col_x[:4])):
        ax.text(cx + col_w[j]/2, 1 - row_h/2, txt,
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="white")

    for i, (label, lr_v, ns_v, smd_v) in enumerate(rows):
        y = 1 - (i + 2) * row_h
        bg = "#F5F5F5" if i % 2 == 0 else "white"
        ax.add_patch(plt.Rectangle((0, y), 1, row_h,
                                    transform=ax.transAxes, fc=bg, zorder=0))

        is_header = (lr_v == "" and ns_v == "" and smd_v == "" and label != "")
        is_col_hdr = (i == 0)
        fw = "bold" if is_header else "normal"
        fs = 10 if is_header else 9.5
        col = "#1A237E" if is_header else "#1A1A1A"

        ax.text(col_x[0] + 0.005, y + row_h/2, label,
                transform=ax.transAxes, ha="left", va="center",
                fontsize=fs, fontweight=fw, color=col)

        for j, val in enumerate([lr_v, ns_v]):
            ax.text(col_x[j+1] + col_w[j+1]/2, y + row_h/2, val,
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9.5, color="#1A1A1A")

        # SMD: color-code
        if smd_v and smd_v != "|SMD|":
            try:
                sv = float(smd_v)
                sc = "#B71C1C" if sv >= 0.1 else GRN
            except:
                sc = "#1A1A1A"
            ax.text(col_x[3] + col_w[3]/2, y + row_h/2, smd_v,
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9.5, color=sc, fontweight="bold" if sv >= 0.1 else "normal")

    # Outer border
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                                fc="none", ec="#37474F", lw=1.5, zorder=5))
    # Column dividers
    for cx in col_x[1:4]:
        ax.plot([cx, cx], [0, 1], transform=ax.transAxes,
                color="#BDBDBD", lw=0.6, zorder=4)

    ax.set_title("Table 1. Baseline Characteristics by Initial ICU Fluid Group\n"
                 "(Primary analytic cohort, before IPTW)",
                 fontsize=13, fontweight="bold", loc="left", pad=8)

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    save("table_1_baseline.png")


# =============================================================================
# MAIN
# =============================================================================
def main(parquet_path):
    print(f"\n{'='*60}")
    print("Publication Figures — LR vs NS / AKI Causal Analysis")
    print(f"{'='*60}\n")

    # Load
    _tmp = tempfile.mkdtemp(prefix="causalai_")
    _loc = os.path.join(_tmp, "cohort.parquet")
    shutil.copy2(parquet_path, _loc)
    df = pd.read_parquet(_loc)
    print(f"Loaded: {len(df):,} rows, {df.shape[1]} columns")

    ENRICHED_COLS = ["ckd","liver_disease","heart_failure","sepsis","cancer",
                     "diabetes","electrolyte_disorder","vasopressor_pre",
                     "gcs_min_pre","hr_max_pre","spo2_min_pre","cr_slope_pre","n_cr_pre","year"]

    # Build cohort
    eligible = df.query("has_baseline_cr==1 and has_followup_cr==1").copy()
    eligible["aki_abs"]   = (eligible["max_fu_creatinine"] - eligible["baseline_creatinine"]) >= 0.3
    eligible["aki_ratio"] = (eligible["max_fu_creatinine"] / eligible["baseline_creatinine"].replace(0, np.nan)) >= 1.5
    eligible["aki_flag"]  = (eligible["aki_abs"] | eligible["aki_ratio"]).astype(int)
    cohort = eligible.query("baseline_creatinine<=4.0").copy()

    cohort["female"]          = (cohort["gender"] == "F").astype(int)
    cohort["race_white"]      = (cohort["race_cat"] == "WHITE").astype(int)
    cohort["race_black"]      = (cohort["race_cat"] == "BLACK").astype(int)
    cohort["log_baseline_cr"] = np.log1p(cohort["baseline_creatinine"])
    cohort["ckd_proxy"]       = (cohort["baseline_creatinine"] > 1.5).astype(int)
    cohort["treat"]           = (cohort["fluid_type"] == "LR").astype(int)

    cu_dummies = pd.get_dummies(cohort["first_careunit"], prefix="cu", drop_first=True).astype(int)
    cohort = pd.concat([cohort.reset_index(drop=True), cu_dummies.reset_index(drop=True)], axis=1)

    BASE_COVS = (["female","anchor_age","race_white","race_black",
                  "log_baseline_cr","ckd_proxy","hours_from_icu_admit","los_hours"]
                 + list(cu_dummies.columns))
    extra_cov = [c for c in ENRICHED_COLS if c in cohort.columns]
    cov_cols  = BASE_COVS + extra_cov
    print(f"Covariates: {len(cov_cols)} ({len(BASE_COVS)} base + {len(extra_cov)} enriched)")

    X = cohort[cov_cols].values
    A = cohort["treat"].values
    ps_full, auc = fit_ps(X, A)
    print(f"PS AUC: {auc:.3f}")

    keep = (ps_full >= 0.05) & (ps_full <= 0.95)
    ct = cohort[keep].copy()
    Xt, At = X[keep], A[keep]
    ps_t, _ = fit_ps(Xt, At)

    p_t = At.mean()
    iptw_t = np.where(At==1, p_t/np.clip(ps_t,1e-6,1), (1-p_t)/np.clip(1-ps_t,1e-6,1))
    cov_df_t = pd.DataFrame(Xt, columns=cov_cols)
    treat_s  = pd.Series(At)
    smd_before = {c: smd_fn(cov_df_t.loc[treat_s==1,c].dropna().values,
                              cov_df_t.loc[treat_s==0,c].dropna().values) for c in cov_cols}
    smd_after  = {c: weighted_smd_fn(c, treat_s, cov_df_t, iptw_t) for c in cov_cols}

    res_csv = os.path.join(SCRIPT_DIR, "results_summary.csv")

    print()
    fig_cohort_flow(df)
    fig_causal_dag()
    fig_ps_overlap(cohort, ps_full, ct.reset_index(drop=True), ps_t)
    fig_love_plot(cov_cols, smd_before, smd_after)
    fig_aipw_forest(res_csv)
    fig_subgroup_forest(res_csv)
    ml, mn = fig_matching(cohort, ps_full)
    fig_matching_balance(cohort, cov_cols, ps_full, ml, mn)
    fig_competing_events(ct)
    fig_table1(cohort)

    print(f"\nAll figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    import matplotlib.ticker
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    args = parser.parse_args()
    main(args.parquet)
