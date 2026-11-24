"""
Organize publication figures for the LR-vs-NS -> 7-day AKI target-trial paper.

This does NOT re-run the analysis. It renders from saved results:
  - the 12 standalone 300-DPI PNGs already in figures/ (produced by make_figures.py)
  - results_summary.csv / evalue_summary.csv (for the caption/annotation values)

Outputs, both at 300 DPI:
  figures/standalone/   the 12 individual figures (clean copies, stable names)
  figures/panels/       4 composite multi-panel figures for the paper + report:
      panel_1_design.png      A: causal DAG            B: cohort flow
      panel_2_balance.png     A: PS overlap            B: covariate love plot
      panel_3_effects.png     A: AIPW forest (SA)      B: subgroup forest
      panel_4_competing.png   A: competing-risks CIF   B: E-value curve
  (standalone Table 1 and Table 2 are kept separate; tables are not paneled.)

Run:
    python3 make_paper_figures.py
"""

import os, shutil, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR    = os.path.join(SCRIPT_DIR, "figures")
STANDALONE = os.path.join(FIG_DIR, "standalone")
PANELS     = os.path.join(FIG_DIR, "panels")
os.makedirs(STANDALONE, exist_ok=True)
os.makedirs(PANELS, exist_ok=True)

DPI = 300
PANEL_FONT = "DejaVu Sans"
plt.rcParams.update({
    "font.family":       PANEL_FONT,
    "savefig.dpi":       DPI,
    "savefig.facecolor": "white",
    "figure.facecolor":  "white",
})

# The 12 standalone figures already produced by make_figures.py
STANDALONE_FIGS = [
    "fig_cohort_flow.png",
    "fig_causal_dag.png",
    "fig_ps_overlap.png",
    "fig_love_plot.png",
    "fig_aipw_results.png",
    "fig_subgroup_forest.png",
    "fig_matching_ps.png",
    "fig_matching_balance.png",
    "fig_competing_events.png",
    "fig_cif_competing.png",
    "fig_evalue.png",
    "table_1_baseline.png",
]


def _require(name):
    p = os.path.join(FIG_DIR, name)
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"Expected source figure missing: {p}\n"
            f"Run make_figures.py first to regenerate the standalone figures."
        )
    return p


def copy_standalone():
    """Keep clean copies of the 12 individual figures under figures/standalone/."""
    print("Organizing standalone figures ->", STANDALONE)
    for name in STANDALONE_FIGS:
        src = _require(name)
        shutil.copy2(src, os.path.join(STANDALONE, name))
        print(f"  ✓ {name}")


def _place_image(ax, img_path, letter, title):
    """Draw one sub-panel: the source PNG and a bold corner letter.

    The source PNGs already carry their own descriptive titles, so we add only
    the panel letter (A/B) to avoid doubled headings.
    """
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.axis("off")
    # Panel letter, top-left, outside the image content
    ax.text(-0.02, 1.01, letter, transform=ax.transAxes,
            fontsize=24, fontweight="bold", va="bottom", ha="left")


def _two_panel(out_name, left, right, suptitle=None,
               width_ratios=(1, 1), figwidth=16):
    """
    Compose a 2-panel figure from two (path, letter, title) specs.
    Figure height is chosen from the taller source image's aspect ratio so
    nothing is distorted.
    """
    specs = [left, right]
    aspects = []
    for path, _, _ in specs:
        h, w = mpimg.imread(path).shape[:2]
        aspects.append(w / h)
    # allocate width to panels by width_ratios; height = max over panels
    wsum = sum(width_ratios)
    panel_w = [figwidth * wr / wsum for wr in width_ratios]
    panel_h = [pw / ar for pw, ar in zip(panel_w, aspects)]
    figheight = max(panel_h) + 0.25

    fig, axes = plt.subplots(
        1, 2, figsize=(figwidth, figheight),
        gridspec_kw={"width_ratios": list(width_ratios)},
    )
    for ax, (path, letter, title) in zip(axes, specs):
        _place_image(ax, path, letter, title)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.97,
                        bottom=0.02, wspace=0.04)
    out = os.path.join(PANELS, out_name)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✓ {out_name}")


def build_panels():
    print("Building composite panels ->", PANELS)

    _two_panel(
        "panel_1_design.png",
        (_require("fig_causal_dag.png"),  "A", "Prespecified causal DAG"),
        (_require("fig_cohort_flow.png"), "B", "Cohort flow"),
        suptitle="Figure 1. Study design",
        width_ratios=(1.55, 1.0),
    )

    _two_panel(
        "panel_2_balance.png",
        (_require("fig_ps_overlap.png"), "A", "Propensity-score overlap (full vs trimmed)"),
        (_require("fig_love_plot.png"),  "B", "Covariate balance before/after IPTW"),
        suptitle="Figure 2. Positivity and covariate balance",
        width_ratios=(1.5, 1.0),
    )

    _two_panel(
        "panel_3_effects.png",
        (_require("fig_aipw_results.png"),    "A", "AIPW risk differences (primary + SA1–SA4)"),
        (_require("fig_subgroup_forest.png"), "B", "Subgroups (sex, ICU type)"),
        suptitle="Figure 3. Treatment-effect estimates",
        width_ratios=(1.0, 1.0),
    )

    _two_panel(
        "panel_4_competing.png",
        (_require("fig_cif_competing.png"), "A", "Competing-risks cumulative incidence"),
        (_require("fig_evalue.png"),        "B", "E-value sensitivity"),
        suptitle="Figure 4. Competing risks and unmeasured confounding",
        width_ratios=(1.7, 1.0),
    )


def summarize():
    """Print the key numbers so figure captions stay consistent with results."""
    rs = os.path.join(SCRIPT_DIR, "results_summary.csv")
    ev = os.path.join(SCRIPT_DIR, "evalue_summary.csv")
    print("\nValue check (for captions):")
    if os.path.exists(rs):
        df = pd.read_csv(rs)
        prim = df.iloc[0]
        print(f"  Primary AIPW RD = {prim['rd']*100:+.2f} pp "
              f"({prim['rd_lower']*100:+.2f}, {prim['rd_upper']*100:+.2f}); n={int(prim['n'])}")
    if os.path.exists(ev):
        e = pd.read_csv(ev).iloc[0]
        print(f"  E-value point = {e['evalue_point']:.2f}, CI bound = {e['evalue_ci']:.2f}")


def main():
    copy_standalone()
    build_panels()
    summarize()
    print("\nDone. Panels in:", PANELS)
    print("      Standalone in:", STANDALONE)


if __name__ == "__main__":
    main()
