"""
Combined MIMIC + IDEALIST — Cross-cohort figures
=================================================
Produces:
  fig_cross_cohort_forest.png   — primary + SAs side-by-side both cohorts
  fig_cross_cohort_primary.png  — single comparison of primary estimates
  fig_unified_table1.png        — unified baseline characteristics table

Run from this directory:
    python3 make_combined_figures.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

MIMIC_CSV    = "/blue/prismap-ai-core/Ahmed/DigitalTwins/MDT/CausalAI/causal_ai/analysis/results_summary.csv"
IDEALIST_CSV = "/blue/prismap-ai-core/Ahmed/DigitalTwins/MDT/CausalAI/idealist_causal/analysis/results_summary_idealist.csv"

LR_COL  = "#1565C0"
NS_COL  = "#E65100"
MIM_COL = "#1565C0"   # blue  = MIMIC
IDL_COL = "#2E7D32"   # green = IDEALIST
RED_SIG = "#B71C1C"
GRY     = "#546E7A"

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.color":        "#E0E0E0",
    "grid.linewidth":    0.6,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.facecolor": "white",
})

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✓ {name}")


# ── Shared SA labels for cross-cohort comparison ──────────────────────────────
SA_MAP = {
    "Primary AIPW (PS-trimmed)": "Primary AIPW",
    "SA1: No PS trimming":        "SA1: No PS trim",
    "SA2: Stricter baseline Cr≤2.0": "SA2: Cr ≤ 2.0",
    "SA3: Composite death-or-AKI":   "SA3: Death-or-AKI",
    "SA4: IPTW-only (no outcome model)": "SA4: IPTW-only",
    "SA5: 48h creatinine window":    "SA5: 48h AKI",
    "SA5: 3-day AKI window":         "SA5: 3-day AKI",
    "SA6a: baseline Cr≤3.0":         "SA6a: Cr ≤ 3.0",
    "PS Matching (1:1 NN caliper=0.2SD)": "PS Matching",
}


def fig_cross_cohort_forest():
    """Side-by-side forest plot: shared SAs for MIMIC and IDEALIST."""
    mim = pd.read_csv(MIMIC_CSV)
    idl = pd.read_csv(IDEALIST_CSV)

    shared_labels = list(SA_MAP.values())[:8]  # Primary + SA1-6a + PS Matching

    def get_rows(df):
        rows = []
        for orig, label in SA_MAP.items():
            match = df[df.analysis == orig]
            if match.empty:
                # try partial
                match = df[df.analysis.str.startswith(orig[:20])]
            if not match.empty:
                r = match.iloc[0]
                rows.append(dict(label=label, rd=r.rd*100,
                                 lo=r.rd_lower*100, hi=r.rd_upper*100, n=int(r.n)))
        return rows

    mim_rows = get_rows(mim)
    idl_rows = get_rows(idl)

    # Align on shared labels
    mim_dict = {r["label"]: r for r in mim_rows}
    idl_dict = {r["label"]: r for r in idl_rows}
    labels = [l for l in shared_labels if l in mim_dict and l in idl_dict]

    fig, axes = plt.subplots(1, 2, figsize=(16, max(5, len(labels)*0.6)),
                             sharey=True)
    fig.subplots_adjust(wspace=0.05)

    for ax, rows_dict, cohort, col in [
        (axes[0], mim_dict, "MIMIC-IV (Critical Care)", MIM_COL),
        (axes[1], idl_dict, "IDEALIST (Perioperative)", IDL_COL),
    ]:
        y = np.arange(len(labels))
        for i, lab in enumerate(labels):
            r = rows_dict[lab]
            sig = (r["lo"] > 0) or (r["hi"] < 0)
            c = RED_SIG if sig else col
            ax.plot([r["lo"], r["hi"]], [i, i], color=c, lw=2.0)
            ax.scatter(r["rd"], i, color=c, s=55, zorder=5)
            ax.text(r["hi"] + 0.05, i,
                    f"n={r['n']:,}  {r['rd']:+.2f} ({r['lo']:+.2f},{r['hi']:+.2f})",
                    va="center", fontsize=8)
        ax.axvline(0, color="black", lw=1.2, ls="--", alpha=0.5)
        ax.set_xlabel("Risk Difference LR − NS (pp)")
        ax.set_title(cohort, fontweight="bold", fontsize=12)
        ax.set_yticks(y)

    axes[0].set_yticklabels(labels, fontsize=10)
    axes[1].set_yticklabels([])

    fig.suptitle("Cross-Cohort Comparison: LR vs NS Risk Differences\n"
                 "MIMIC-IV (Critical Care) vs IDEALIST (Perioperative Surgery)",
                 fontsize=13, fontweight="bold", y=1.01)

    red_p  = mpatches.Patch(color=RED_SIG, label="Statistically significant")
    blu_p  = mpatches.Patch(color=MIM_COL, label="MIMIC non-significant")
    grn_p  = mpatches.Patch(color=IDL_COL, label="IDEALIST non-significant")
    fig.legend(handles=[red_p, blu_p, grn_p], loc="lower center",
               ncol=3, bbox_to_anchor=(0.5, -0.04), fontsize=9)

    save("fig_cross_cohort_forest.png")


def fig_primary_comparison():
    """Single panel: primary estimates for both cohorts with CI."""
    mim = pd.read_csv(MIMIC_CSV)
    idl = pd.read_csv(IDEALIST_CSV)

    mp = mim[mim.analysis == "Primary AIPW (PS-trimmed)"].iloc[0]
    ip = idl[idl.analysis == "Primary AIPW (PS-trimmed)"].iloc[0]

    fig, ax = plt.subplots(figsize=(8, 4))

    # MIMIC
    ax.plot([mp.rd_lower*100, mp.rd_upper*100], [1, 1], color=MIM_COL, lw=3)
    ax.scatter(mp.rd*100, 1, color=MIM_COL, s=100, zorder=5)
    ax.text(mp.rd_upper*100 + 0.1, 1,
            f"MIMIC-IV  RD {mp.rd*100:+.2f} pp  (95% CI {mp.rd_lower*100:+.2f}, {mp.rd_upper*100:+.2f})\n"
            f"n={int(mp.n):,}  LR {mp.risk_lr*100:.1f}%  NS {mp.risk_ns*100:.1f}%  E-value 1.08",
            va="center", fontsize=10, color=MIM_COL)

    # IDEALIST
    ax.plot([ip.rd_lower*100, ip.rd_upper*100], [0, 0], color=IDL_COL, lw=3)
    ax.scatter(ip.rd*100, 0, color=IDL_COL, s=100, zorder=5)
    ax.text(ip.rd_upper*100 + 0.1, 0,
            f"IDEALIST  RD {ip.rd*100:+.2f} pp  (95% CI {ip.rd_lower*100:+.2f}, {ip.rd_upper*100:+.2f})\n"
            f"n={int(ip.n):,}  LR {ip.risk_lr*100:.1f}%  NS {ip.risk_ns*100:.1f}%  E-value 1.32",
            va="center", fontsize=10, color=IDL_COL)

    ax.axvline(0, color="black", lw=1.5, ls="--", alpha=0.6)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["IDEALIST\n(Perioperative)", "MIMIC-IV\n(Critical Care)"],
                       fontsize=11)
    ax.set_xlabel("Risk Difference LR − NS (percentage points)")
    ax.set_xlim(-3, 4)
    ax.set_title("Primary AIPW Estimates: Critical Care vs Perioperative Setting",
                 fontweight="bold", fontsize=12)
    ax.grid(axis="x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    save("fig_primary_comparison.png")


def fig_unified_table1():
    """Unified Table 1: key characteristics side-by-side."""
    data = {
        "Characteristic": [
            "Total eligible", "PS-trimmed analytic cohort",
            "LR arm, n (%)", "NS arm, n (%)",
            "Unadjusted AKI rate — LR", "Unadjusted AKI rate — NS",
            "Adjusted AKI RD (pp)", "95% CI",
            "E-value", "Fine-Gray sHR (95% CI)",
            "PS AUC", "Covariates balanced |SMD|<0.1",
        ],
        "MIMIC-IV (Critical Care)": [
            "41,383", "38,266",
            "14,309 (34.6%)", "27,074 (65.4%)",
            "22.2%", "22.4%",
            "−0.13", "(−1.31, +0.87)",
            "1.08 (crosses null)", "1.02 (0.96–1.08)",
            "~0.70", "40/40",
        ],
        "IDEALIST (Perioperative)": [
            "250,532", "232,539",
            "148,912 (59.0%)", "101,620 (41.0%)",
            "11.9%", "15.7%",
            "−0.93", "(−0.91, −0.07)",
            "1.32 (does not cross null)", "0.956 (0.931–0.982)",
            "~0.65", "40/40",
        ],
    }
    df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.axis("off")
    tbl = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc="left", loc="center",
        colWidths=[0.38, 0.31, 0.31],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.scale(1, 1.6)

    # Header style
    for j in range(3):
        tbl[(0, j)].set_facecolor("#1A237E")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    # Alternating rows
    for i in range(1, len(df)+1):
        fc = "#F5F5F5" if i % 2 == 0 else "white"
        for j in range(3):
            tbl[(i, j)].set_facecolor(fc)

    # Highlight primary result rows
    for i in [7, 8]:  # RD and CI rows
        tbl[(i, 1)].set_facecolor("#E3F2FD")
        tbl[(i, 2)].set_facecolor("#E8F5E9")

    ax.set_title("Table 1. Cross-Cohort Comparison: MIMIC-IV vs IDEALIST",
                 fontsize=13, fontweight="bold", pad=10)
    save("fig_unified_table1.png")


def main():
    print("Building combined figures…")
    fig_primary_comparison()
    fig_cross_cohort_forest()
    fig_unified_table1()
    print(f"\nAll figures → {FIG_DIR}")


if __name__ == "__main__":
    main()
