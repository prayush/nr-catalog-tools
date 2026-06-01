#!/usr/bin/env python
"""Compute per-mode NR vs NRSur7dq4 match as a function of total binary mass.

Runs all four pilot simulations over a grid of total masses and produces a
figure of log10(1 - match) vs M for each mode.

Usage
-----
    python mass_scan.py [--outdir results]
"""

from __future__ import annotations

import os
import sys
import csv
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

from nrcatalogtools.surrogate import generate_surrogate_modes, SURROGATE_MODES
from nrcatalogtools.waveform.matching import (
    compute_mode_match,
    compute_phase_diff_per_cycle,
    mode_f_lower,
)
from nrcatalogtools import load_catalog

# ── Configuration ─────────────────────────────────────────────────────────────

SIMS = [
    ("SXS", "SXS:BBH:0001", r"$q=1,\;\chi=0$"),
    ("SXS", "SXS:BBH:0005", r"$q=1,\;\chi_{1z}=+0.5$"),
    ("SXS", "SXS:BBH:0169", r"$q=2,\;\chi=0$"),
    ("SXS", "SXS:BBH:0162", r"$q=2,\;\chi_{1z}=+0.6$"),
]

MASS_GRID = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100]  # M_sun

PSD_NAME = "aLIGOZeroDetHighPower"
DELTA_T = 1.0 / 4096  # seconds
DISTANCE = 1.0  # Mpc

# Colours and line styles per simulation
SIM_STYLES = {
    "SXS:BBH:0001": dict(color="C0", ls="-", marker="o", ms=5),
    "SXS:BBH:0005": dict(color="C1", ls="--", marker="s", ms=5),
    "SXS:BBH:0169": dict(color="C2", ls="-", marker="^", ms=5),
    "SXS:BBH:0162": dict(color="C3", ls="--", marker="D", ms=5),
}


# ── Core per-mass comparison ──────────────────────────────────────────────────


def compare_one_mass(wfm, cat, sim_name, total_mass):
    """Return per-mode {match, phase_diff_per_cycle, f_lower_mode} for one mass.

    Parameters
    ----------
    wfm : WaveformModes
        Pre-loaded NR waveform object.
    cat : CatalogBase
        Pre-loaded catalog (needed for get_parameters).
    sim_name : str
    total_mass : float

    Returns
    -------
    dict  {(ell, em): {'match', 'phase_diff_per_cycle', 'f_lower_mode'}}
        or None if surrogate generation failed entirely.
    """
    params = cat.get_parameters(sim_name, total_mass=total_mass)
    f_lower = params["f_lower"]

    try:
        h_sur, f_lower_sur = generate_surrogate_modes(
            params,
            total_mass=total_mass,
            distance=DISTANCE,
            delta_t_seconds=DELTA_T,
        )
    except Exception as exc:
        print(f"      M={total_mass}: surrogate failed — {exc}")
        return None

    f_lower_match = max(f_lower, f_lower_sur)

    result = {}
    for (ell, em) in SURROGATE_MODES:
        try:
            h_nr_complex = wfm.get_mode(
                ell,
                em,
                total_mass=total_mass,
                distance=DISTANCE,
                delta_t_seconds=DELTA_T,
            )
        except Exception:
            result[(ell, em)] = {
                "match": float("nan"),
                "phase_diff_per_cycle": float("nan"),
                "f_lower_mode": float("nan"),
            }
            continue

        h_nr = h_nr_complex.real()

        if (ell, em) not in h_sur:
            result[(ell, em)] = {
                "match": float("nan"),
                "phase_diff_per_cycle": float("nan"),
                "f_lower_mode": float("nan"),
            }
            continue

        h_sur_mode = h_sur[(ell, em)].real()
        f_low_mode = mode_f_lower(f_lower_match, em)

        mm = compute_mode_match(h_nr, h_sur_mode, f_low_mode, psd_name=PSD_NAME)
        dp, _ = compute_phase_diff_per_cycle(h_nr_complex, h_sur[(ell, em)])

        result[(ell, em)] = {
            "match": mm,
            "phase_diff_per_cycle": dp,
            "f_lower_mode": f_low_mode,
        }

    return result


# ── Main ──────────────────────────────────────────────────────────────────────


def run_mass_scan(outdir=None, figsdir=None):
    if outdir is None:
        outdir = os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "results"))
    if figsdir is None:
        figsdir = os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "figs"))

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(figsdir, exist_ok=True)

    # results[sim_name][total_mass] = {(ell,em): {...}}
    all_results = {}

    for catalog_name, sim_name, label in SIMS:
        print(f"\n{'='*60}")
        print(f"Loading {sim_name} from {catalog_name} catalog...")
        cat = load_catalog(catalog_name)
        wfm = cat.get(sim_name)
        print("  Waveform loaded. Starting mass scan...")

        sim_results = {}
        for M in MASS_GRID:
            print(f"  M = {M:5.0f} M_sun ...", end="", flush=True)
            res = compare_one_mass(wfm, cat, sim_name, M)
            sim_results[M] = res
            if res is not None:
                mm22 = res.get((2, 2), {}).get("match", float("nan"))
                print(
                    f"  match(2,2) = {mm22:.4f}"
                    if not np.isnan(mm22)
                    else "  match(2,2) = N/A"
                )
            else:
                print("  FAILED")

        all_results[sim_name] = sim_results

    # ── Write CSV ─────────────────────────────────────────────────────────────
    csv_path = os.path.join(outdir, "mass_scan_results.csv")
    rows = []
    for catalog_name, sim_name, label in SIMS:
        for M, res in all_results[sim_name].items():
            if res is None:
                continue
            for (ell, em), vals in res.items():
                rows.append(
                    {
                        "sim_id": sim_name,
                        "total_mass": M,
                        "ell": ell,
                        "em": em,
                        "match": vals["match"],
                        "phase_diff_per_cycle": vals["phase_diff_per_cycle"],
                        "f_lower_mode": vals["f_lower_mode"],
                    }
                )
    if rows:
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nCSV written: {csv_path}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    _plot_mass_scan(all_results, figsdir)

    return all_results


def _plot_mass_scan(all_results, outdir):
    n_modes = len(SURROGATE_MODES)
    ncols = 3
    nrows = (n_modes + ncols - 1) // ncols  # = 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 8), sharex=True)
    axes = axes.flatten()

    for ax_idx, (ell, em) in enumerate(SURROGATE_MODES):
        ax = axes[ax_idx]
        mode_label = f"$({ell},{em:+d})$"

        for catalog_name, sim_name, sim_label in SIMS:
            masses, log_mismatches = [], []
            for M in MASS_GRID:
                res = all_results[sim_name].get(M)
                if res is None:
                    continue
                mm = res.get((ell, em), {}).get("match", float("nan"))
                if np.isnan(mm) or mm >= 1.0 or mm <= 0.0:
                    continue
                masses.append(M)
                log_mismatches.append(np.log10(1.0 - mm))

            if masses:
                style = SIM_STYLES[sim_name]
                ax.plot(
                    masses,
                    log_mismatches,
                    label=sim_label,
                    color=style["color"],
                    ls=style["ls"],
                    marker=style["marker"],
                    ms=style["ms"],
                    lw=1.5,
                    alpha=0.9,
                )

        ax.set_title(mode_label, fontsize=11)
        ax.set_ylabel(r"$\log_{10}(1 - \mathcal{F})$", fontsize=9)
        ybot = -4.0
        ytop = 0.5
        # Shaded bands (drawn before data so they sit behind the lines)
        ax.axhspan(
            np.log10(0.1), ytop, facecolor="#606060", alpha=0.18, zorder=0
        )  # match < 0.90
        ax.axhspan(
            np.log10(0.03), np.log10(0.1), facecolor="#909090", alpha=0.15, zorder=0
        )  # 0.90–0.97
        ax.axhspan(
            np.log10(0.01), np.log10(0.03), facecolor="#c0c0c0", alpha=0.15, zorder=0
        )  # 0.97–0.99
        ax.set_ylim(ybot, ytop)
        ax.grid(True, which="both", alpha=0.2)
        ax.tick_params(axis="both", labelsize=8)

    # x-axis labels on bottom row
    for ax in axes[ncols:]:
        ax.set_xlabel(r"$M_{\rm tot}\ [M_\odot]$", fontsize=10)

    # Remove unused axes if n_modes < nrows*ncols
    for ax_idx in range(n_modes, len(axes)):
        fig.delaxes(axes[ax_idx])

    # Band labels on the right edge of every panel
    for ax in axes[:n_modes]:
        ax.text(
            MASS_GRID[-1] + 1,
            (np.log10(0.1) + 0.5) / 2,
            "< 0.90",
            va="center",
            fontsize=6,
            color="#404040",
        )
        ax.text(
            MASS_GRID[-1] + 1,
            (np.log10(0.1) + np.log10(0.03)) / 2,
            "0.97",
            va="center",
            fontsize=6,
            color="#555555",
        )
        ax.text(
            MASS_GRID[-1] + 1,
            (np.log10(0.03) + np.log10(0.01)) / 2,
            "0.99",
            va="center",
            fontsize=6,
            color="#777777",
        )

    # Legend in top-left panel (2,2), which has the most free space
    axes[0].legend(fontsize=8, loc="upper left")
    # Remove legend from last panel where it was before
    axes[n_modes - 1].legend().remove()

    fig.suptitle(
        r"NR vs NRSur7dq4: per-mode mismatch as a function of $M_{\rm tot}$"
        "\n(aLIGO ZDHP noise curve)",
        fontsize=11,
    )
    plt.tight_layout(rect=[0, 0, 0.97, 1])

    fig_path = os.path.join(outdir, "mass_scan_mismatch.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure written: {fig_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--outdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "results")),
        help="Output directory for CSV results (default: project/results)",
    )
    p.add_argument(
        "--figsdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "figs")),
        help="Output directory for figures (default: project/figs)",
    )
    args = p.parse_args()
    run_mass_scan(outdir=args.outdir, figsdir=args.figsdir)
