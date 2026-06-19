#!/usr/bin/env python
"""Compare one NR simulation against NRSur7dq4 mode by mode.

For each of the 7 NR modes (NRSur7dq4 supports ell ≤ 4; the (5,5)
mode is unavailable from the surrogate and appears with match=N/A), this script:

1. Loads the NR waveform from the requested catalog.
2. Calls ``get_parameters()`` to extract intrinsic parameters.
3. Generates NRSur7dq4 modes at those parameters, rescaled to *total_mass*.
4. Computes ``pycbc.filter.match()`` per mode (maximised over time and phase).
5. Writes a CSV to ``results/<sim_id>_mode_matches.csv``.
6. Saves a 3-panel figure (amplitude comparison + phase merger + match bars).

Usage
-----
    python compare_one_sim_vs_surrogate.py \\
        --catalog SXS \\
        --sim SXS:BBH:0001 \\
        --total-mass 40 \\
        --psd aLIGOZeroDetHighPower \\
        --outdir results

    python compare_one_sim_vs_surrogate.py \\
        --catalog RIT \\
        --sim "RIT:BBH:0001-n100-id3" \\
        --total-mass 40

    # Use --rotate to also run SO(3)-optimized match (slower)
    python compare_one_sim_vs_surrogate.py --catalog SXS --sim SXS:BBH:0001 --rotate
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

from nrcatalogtools.comparisons import compare_sim_vs_surrogate, DELTA_T


def _build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--catalog",
        default="SXS",
        help="Catalog name: SXS, RIT, or MAYA (default: SXS)",
    )
    p.add_argument(
        "--sim", default="SXS:BBH:0001", help="Simulation ID (default: SXS:BBH:0001)"
    )
    p.add_argument(
        "--total-mass",
        type=float,
        default=40.0,
        help="Total binary mass in solar masses (default: 40)",
    )
    p.add_argument(
        "--psd",
        default="aLIGOZeroDetHighPower",
        help="PyCBC analytic PSD name (default: aLIGOZeroDetHighPower)",
    )
    p.add_argument(
        "--delta-t",
        type=float,
        default=DELTA_T,
        help=f"Sample spacing in seconds (default: {DELTA_T})",
    )
    p.add_argument(
        "--outdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "results")),
        help="Output directory for CSV matches (default: project/results)",
    )
    p.add_argument(
        "--figsdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "figs")),
        help="Output directory for figures (default: project/figs)",
    )
    p.add_argument(
        "--rotate",
        action="store_true",
        help="Compute SO(3)-rotation-optimized match (slow; auto-enabled for precessing systems)",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    compare_sim_vs_surrogate(
        catalog_name=args.catalog,
        sim_name=args.sim,
        total_mass=args.total_mass,
        psd_name=args.psd,
        outdir=args.outdir,
        figsdir=args.figsdir,
        delta_t=args.delta_t,
        rotate=args.rotate,
    )
