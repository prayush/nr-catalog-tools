#!/usr/bin/env python
"""
Generate in-depth figures for mismatch and phase error as a function of eccentricity.
Focuses on both eccentric ('a', 'c', 'e') and quasi-circular ('b', 'd', 'f') categories.
Produces:
- fig_eccentricity_mismatch.png
- fig_eccentricity_phase.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(_SCRIPTS_DIR)
_RESULTS_DIR = os.path.join(_PROJECT_DIR, "results")
_FIGS_DIR = os.path.join(_PROJECT_DIR, "figs")

# Aesthetics
SERIES_CONFIG = [
    {
        "name": "SXS (train)",
        "cat": "SXS",
        "cal": True,
        "color": "#2166ac",
        "marker": "o",
        "zorder": 4,
        "alpha": 0.4,
    },
    {
        "name": "SXS (test)",
        "cat": "SXS",
        "cal": False,
        "color": "#4393c3",
        "marker": "x",
        "zorder": 3,
        "alpha": 0.8,
    },
    {
        "name": "RIT",
        "cat": "RIT",
        "cal": None,
        "color": "#d6604d",
        "marker": "^",
        "zorder": 2,
        "alpha": 0.6,
    },
    {
        "name": "MAYA",
        "cat": "MAYA",
        "cal": None,
        "color": "#4dac26",
        "marker": "v",
        "zorder": 1,
        "alpha": 0.6,
    },
]
MODE_LABEL = {
    "22": r"$(2,2)$",
    "21": r"$(2,1)$",
    "33": r"$(3,3)$",
    "44": r"$(4,4)$",
    "32": r"$(3,2)$",
    "43": r"$(4,3)$",
}
MODES = ["22", "21", "33", "44", "32", "43"]

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
    }
)


def load_data():
    merged = os.path.join(_RESULTS_DIR, "batch_aligned_all.csv")
    if not os.path.exists(merged):
        raise FileNotFoundError(f"Missing {merged}")
    df = pd.read_csv(merged)
    if "error" in df.columns:
        df = df[df["error"].isna() | (df["error"] == "")]

    for m in MODES:
        match_col = f"match_{m}"
        if match_col in df.columns:
            df[f"mismatch_{m}"] = 1.0 - df[match_col]

    return df


def plot_eccentricity_dependence(
    df, metric="mismatch", outname="fig_eccentricity.png", log_y=True, ylabel=""
):
    n_modes = len(MODES)
    ncols = 3
    nrows = (n_modes + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows), sharex=True
    )

    # We want to distinguish catalogs. We can just use scatter.
    # We will only plot where eccentricity > 0
    df_plot = df[df["eccentricity"] > 0].copy()

    for idx, mode in enumerate(MODES):
        ax = axes[idx // ncols][idx % ncols]
        col = f"{metric}_{mode}"

        if col not in df_plot.columns:
            ax.set_visible(False)
            continue

        for config in SERIES_CONFIG:
            if config["cal"] is not None:
                sub = df_plot[
                    (df_plot["catalog"] == config["cat"])
                    & (df_plot["nrsur7dq4_calibration"] == config["cal"])
                ]
            else:
                sub = df_plot[df_plot["catalog"] == config["cat"]]

            xv = sub["eccentricity"]
            yv = sub[col]

            valid = np.isfinite(xv) & np.isfinite(yv) & (xv > 1e-5)
            if log_y:
                valid &= yv > 0

            if valid.sum() == 0:
                continue

            ax.scatter(
                xv[valid],
                yv[valid],
                c=config["color"],
                marker=config["marker"],
                s=15,
                alpha=config["alpha"],
                edgecolors="none" if config["marker"] != "x" else config["color"],
                zorder=config["zorder"],
                label=config["name"] if idx == 0 else "",
            )

        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_title(MODE_LABEL[mode], fontsize=11)
        ax.grid(True, lw=0.4, alpha=0.4)
        if idx >= n_modes - ncols:
            ax.set_xlabel("Eccentricity $e$", fontsize=10)
        if idx % ncols == 0:
            ax.set_ylabel(ylabel, fontsize=10)

    # Add legend to the first axis
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        axes[0][0].legend(handles, labels, loc="best", fontsize=9, framealpha=0.7)

    plt.tight_layout()
    os.makedirs(_FIGS_DIR, exist_ok=True)
    outpath = os.path.join(_FIGS_DIR, outname)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


def main():
    print("Loading data...")
    df = load_data()

    print("Generating mismatch vs eccentricity figure...")
    plot_eccentricity_dependence(
        df,
        metric="mismatch",
        outname="fig_eccentricity_mismatch.png",
        log_y=True,
        ylabel=r"Mismatch $1 - \mathcal{F}$",
    )

    print("Generating phase-difference vs eccentricity figure...")
    plot_eccentricity_dependence(
        df,
        metric="phase_diff",
        outname="fig_eccentricity_phase.png",
        log_y=True,
        ylabel=r"$\Delta\Phi$/cycle [rad]",
    )
    print("Done.")


if __name__ == "__main__":
    main()
