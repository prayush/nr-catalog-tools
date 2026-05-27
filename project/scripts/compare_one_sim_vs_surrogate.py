#!/usr/bin/env python
"""Compare one NR simulation against NRSur7dq4 mode by mode.

For each of the 7 NR modes (NRSur7dq4 supports ell ≤ 4; the (5,5)
mode is unavailable from the surrogate and appears with match=N/A), this script:

1. Loads the NR waveform from the requested catalog.
2. Calls ``get_parameters()`` to extract intrinsic parameters.
3. Generates NRSur7dq4 modes at those parameters, rescaled to *total_mass*.
4. Computes ``pycbc.filter.match()`` per mode (maximised over time and phase).
5. Writes a CSV to ``results/<sim_id>_mode_matches.csv``.
6. Saves a 2-panel figure (amplitude comparison + match summary).

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

    # Use --rotate to also run SO(3)-optimized match (Step 2, slower)
    python compare_one_sim_vs_surrogate.py --catalog SXS --sim SXS:BBH:0001 --rotate

Convention notes
----------------
* NRSur7dq4 is precessing and returns all modes with |m| ≤ ell ≤ 4 as
  complex h_lm arrays (same spin-weight -2 SWSH convention as NR catalogs).
* Matches are computed on the real part (h₊ component) of each mode.
* f_lower for mode (ell, m) is set to ``f_lower * |m| / 2`` (GW frequency
  ≈ |m| × orbital frequency).
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make local helpers importable regardless of cwd.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPTS_DIR)
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

from surrogate_utils import (
    generate_surrogate_modes,
    check_surrogate_prior,
    NR_MODES,
)
from match_utils import compute_mode_match, compute_phase_diff_per_cycle, mode_f_lower
from catalog_utils import load_catalog

# ── Constants ────────────────────────────────────────────────────────────────

DELTA_T = 1.0 / 4096  # seconds
DISTANCE = 1.0  # Mpc (amplitude-irrelevant for match, kept consistent)


# ── Helpers ──────────────────────────────────────────────────────────────────


def _safe_sim_id(sim_name: str) -> str:
    """Convert sim name to a filesystem-safe string."""
    return sim_name.replace(":", "_").replace("/", "_")


def _waveform_duration(wfm, total_mass: float) -> float:
    """Estimate waveform duration in seconds from the (2,2) mode."""
    from nrcatalogtools import utils

    m_secs = utils.time_to_physical(total_mass)
    return float(wfm.time[-1] - wfm.time[0]) * m_secs


# ── Core comparison ──────────────────────────────────────────────────────────


def compare_sim_vs_surrogate(
    catalog_name: str,
    sim_name: str,
    total_mass: float = 40.0,
    psd_name: str = "aLIGOZeroDetHighPower",
    outdir: str = "results",
    delta_t: float = DELTA_T,
    rotate: bool = False,
) -> dict:
    """Run the full NR vs NRSur7dq4 comparison for one simulation.

    Parameters
    ----------
    catalog_name : str
        One of ``'SXS'``, ``'RIT'``, ``'MAYA'``.
    sim_name : str
        Simulation identifier (catalog-specific, e.g. ``'SXS:BBH:0001'``).
    total_mass : float, optional
        Total mass in solar masses (default 40).
    psd_name : str, optional
        PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``).
    outdir : str, optional
        Output directory for CSV and figure (default ``'results'``).
    delta_t : float, optional
        Sample spacing in physical seconds (default 1/4096).
    rotate : bool, optional
        If True, also compute the SO(3)-rotation-optimized match via
        ``WaveformModes.match_sphere_averaged()`` (slow).

    Returns
    -------
    dict
        Per-mode match results::

            {
              (ell, em): {
                  'match': float,
                  'f_lower_mode': float,
                  'match_rotated': float or None,
              },
              ...
            }
    """
    os.makedirs(outdir, exist_ok=True)

    # ── 1. Load catalog and waveform ─────────────────────────────────────────
    print(f"\n[1/6] Loading {catalog_name} catalog...")
    cat = load_catalog(catalog_name)
    print(f"      Fetching {sim_name}...")
    wfm = cat.get(sim_name)

    # ── 2. Extract parameters ─────────────────────────────────────────────────
    print(f"[2/6] Extracting source parameters (M={total_mass} M☉)...")
    params = cat.get_parameters(sim_name, total_mass=total_mass)
    q = params["mass1"] / params["mass2"]
    chi1 = np.sqrt(
        params["spin1x"] ** 2 + params["spin1y"] ** 2 + params["spin1z"] ** 2
    )
    chi2 = np.sqrt(
        params["spin2x"] ** 2 + params["spin2y"] ** 2 + params["spin2z"] ** 2
    )
    f_lower = params["f_lower"]

    print(
        f"      q={q:.4f}  |χ₁|={chi1:.4f}  |χ₂|={chi2:.4f}  f_lower={f_lower:.2f} Hz"
    )

    if not check_surrogate_prior(params):
        print(
            "WARNING: parameters lie outside NRSur7dq4 prior (q > 4 or |χ| > 0.8). "
            "Proceeding, but surrogate extrapolation may be unreliable."
        )

    # ── 3. Generate surrogate modes ──────────────────────────────────────────
    print(f"[3/6] Generating NRSur7dq4 modes (M={total_mass} M☉, D={DISTANCE} Mpc)...")
    h_sur, f_lower_sur = generate_surrogate_modes(
        params, total_mass=total_mass, distance=DISTANCE, delta_t_seconds=delta_t
    )
    print(f"      Generated {len(h_sur)} surrogate modes.")
    if f_lower_sur > f_lower * 1.05:
        print(
            f"      NOTE: surrogate starts at f_GW={f_lower_sur:.1f} Hz "
            f"(NRSur7dq4 minimum for these params), above NR f_lower={f_lower:.1f} Hz. "
            f"Match f_lower raised to {f_lower_sur:.1f} Hz."
        )
    # Use the larger of the two as the effective lower cutoff for the match.
    f_lower_match = max(f_lower, f_lower_sur)

    # ── 4. (PSD is built per-mode inside compute_mode_match to guarantee
    #         delta_f consistency after length alignment.)
    print(f"[4/6] PSD: {psd_name} (built per-mode at matched frequency resolution)")

    # ── 5. Per-mode match ─────────────────────────────────────────────────────
    print("[5/6] Computing per-mode matches...")
    results = {}

    for (ell, em) in NR_MODES:
        # NR mode in physical units
        try:
            h_nr_complex = wfm.get_mode(
                ell,
                em,
                total_mass=total_mass,
                distance=DISTANCE,
                delta_t_seconds=delta_t,
            )
        except Exception as exc:
            print(f"      ({ell},{em:+d}): NR mode unavailable — {exc}")
            results[(ell, em)] = {
                "match": float("nan"),
                "f_lower_mode": float("nan"),
                "phase_diff_per_cycle": float("nan"),
                "n_cycles": float("nan"),
                "match_rotated": None,
            }
            continue

        h_nr = h_nr_complex.real()

        # Surrogate mode
        if (ell, em) not in h_sur:
            print(f"      ({ell},{em:+d}): surrogate mode unavailable (ell > 4?)")
            results[(ell, em)] = {
                "match": float("nan"),
                "f_lower_mode": float("nan"),
                "phase_diff_per_cycle": float("nan"),
                "n_cycles": float("nan"),
                "match_rotated": None,
            }
            continue

        h_sur_mode = h_sur[(ell, em)].real()

        f_low_mode = mode_f_lower(f_lower_match, em)
        mm = compute_mode_match(h_nr, h_sur_mode, f_low_mode, psd_name=psd_name)
        dphase, n_cycles = compute_phase_diff_per_cycle(h_nr_complex, h_sur[(ell, em)])
        results[(ell, em)] = {
            "match": mm,
            "f_lower_mode": f_low_mode,
            "phase_diff_per_cycle": dphase,
            "n_cycles": n_cycles,
            "match_rotated": None,
        }
        flag = "" if np.isnan(mm) else f"{mm:.6f}"
        dp_str = "N/A" if np.isnan(dphase) else f"{dphase:.4f} rad"
        print(
            f"      ({ell},{em:+d}): match = {flag}  phase_diff/cycle = {dp_str}"
            f"  [f_lower_mode={f_low_mode:.1f} Hz]"
        )

    # ── 5b. Optional SO(3)-rotation-optimized match ───────────────────────────
    if rotate:
        print("[5b] Computing SO(3)-rotation-optimized match (Nelder-Mead)...")
        # Build a surrogate WaveformModes object to pass to match_sphere_averaged.
        # For now we fall back to the sphere-averaged method on the WaveformModes
        # side, passing the surrogate modes as a plain dict.
        try:
            from match_utils import load_psd

            dur_nr = _waveform_duration(wfm, total_mass)
            dur_sur = len(next(iter(h_sur.values()))) * delta_t
            psd_rot = load_psd(
                f_lower_match, delta_t, max(dur_nr, dur_sur) * 1.1, psd_name=psd_name
            )
            mm_rot = wfm.match_sphere_averaged(
                _make_waveform_modes_from_dict(h_sur, wfm),
                psd=psd_rot,
                f_lower=f_lower_match,
                delta_t=delta_t,
            )
            print(f"      SO(3)-optimized match = {mm_rot:.6f}")
            for key in results:
                results[key]["match_rotated"] = mm_rot
        except Exception as exc:
            print(f"      SO(3)-optimized match failed: {exc}")

    # ── 6. Output ─────────────────────────────────────────────────────────────
    print(f"[6/6] Writing outputs to {outdir}/...")
    _write_csv(results, sim_name, catalog_name, total_mass, params, outdir)
    _plot(
        results,
        wfm,
        h_sur,
        sim_name,
        catalog_name,
        total_mass,
        params,
        delta_t,
        psd_name,
        outdir,
    )
    _print_table(results, sim_name)

    return results


def _make_waveform_modes_from_dict(h_sur_dict, wfm_ref):
    """Wrap surrogate mode dict into a WaveformModes object for rotation matching.

    Uses the NR waveform's time/ell grid as the reference frame.
    This is a best-effort wrapper; mismatched time grids are handled by
    ``match_sphere_averaged``'s mode-by-mode resampling.
    """
    # This is the surrogate_wrapper.py functionality (Step 2).
    # For now, return the raw dict — match_sphere_averaged accepts dicts directly
    # when the NR side calls get_mode() per mode in its objective function.
    # A proper WaveformModes wrapper will be added in surrogate_wrapper.py.
    return h_sur_dict


# ── Output helpers ────────────────────────────────────────────────────────────


def _write_csv(results, sim_name, catalog_name, total_mass, params, outdir):
    safe_id = _safe_sim_id(sim_name)
    path = os.path.join(outdir, f"{safe_id}_mode_matches.csv")
    rows = []
    for (ell, em), res in results.items():
        rows.append(
            {
                "sim_id": sim_name,
                "catalog": catalog_name,
                "total_mass": total_mass,
                "q": params["mass1"] / params["mass2"],
                "chi1z": params["spin1z"],
                "chi2z": params["spin2z"],
                "ell": ell,
                "em": em,
                "match": res["match"],
                "f_lower_mode": res["f_lower_mode"],
                "phase_diff_per_cycle": res["phase_diff_per_cycle"],
                "n_cycles": res["n_cycles"],
                "match_rotated": res["match_rotated"]
                if res["match_rotated"] is not None
                else "",
            }
        )
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"      CSV written: {path}")


def _plot(
    results,
    wfm,
    h_sur,
    sim_name,
    catalog_name,
    total_mass,
    params,
    delta_t,
    psd_name,
    outdir,
):
    """3-panel figure: (2,2) amplitude, (2,2) phase, match bar chart."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 11))
    fig.suptitle(
        f"{catalog_name}: {sim_name}\n"
        f"q={params['mass1']/params['mass2']:.3f}  "
        f"χ₁z={params['spin1z']:.3f}  χ₂z={params['spin2z']:.3f}  "
        f"M={total_mass} M☉",
        fontsize=11,
    )

    # ── Panel 1: (2,2) amplitude ─────────────────────────────────────────────
    ax = axes[0]
    try:
        h22_nr = wfm.get_mode(
            2, 2, total_mass=total_mass, distance=DISTANCE, delta_t_seconds=delta_t
        )
        t_nr = np.array(h22_nr.sample_times)
        amp_nr = np.abs(np.array(h22_nr))
        ax.plot(t_nr, amp_nr, color="C0", lw=1.5, label="NR", alpha=0.9)
    except Exception:
        pass

    if (2, 2) in h_sur:
        h22_sur = h_sur[(2, 2)]
        t_sur = np.array(h22_sur.sample_times)
        amp_sur = np.abs(np.array(h22_sur))
        ax.plot(
            t_sur, amp_sur, color="C1", lw=1.5, ls="--", label="NRSur7dq4", alpha=0.9
        )

    ax.set_xlabel("t - t_peak [s]")
    ax.set_ylabel("|h₂₂| [strain]")
    ax.set_title("(2,2) mode amplitude")
    ax.legend(fontsize=9)
    ax.set_xlim(left=max(-2.0, ax.get_xlim()[0]))  # zoom in to last 2 s before merger

    # ── Panel 2: (2,2) real part near merger ─────────────────────────────────
    ax = axes[1]
    try:
        t_nr_full = np.array(h22_nr.sample_times)
        re_nr = np.array(h22_nr.real())
        mask = t_nr_full >= -0.25  # last 250 ms
        ax.plot(t_nr_full[mask], re_nr[mask], color="C0", lw=1.2, label="NR")
    except Exception:
        pass

    if (2, 2) in h_sur:
        t_sur_full = np.array(h22_sur.sample_times)
        re_sur = np.array(h22_sur.real())
        mask = t_sur_full >= -0.25
        ax.plot(
            t_sur_full[mask],
            re_sur[mask],
            color="C1",
            lw=1.2,
            ls="--",
            label="NRSur7dq4",
        )

    ax.set_xlabel("t - t_peak [s]")
    ax.set_ylabel("Re(h₂₂) [strain]")
    ax.set_title("(2,2) real part — merger region")
    ax.legend(fontsize=9)

    # ── Panel 3: match bar chart ──────────────────────────────────────────────
    ax = axes[2]
    mode_labels = [f"({ell},{em:+d})" for (ell, em) in NR_MODES]
    match_vals = [
        results.get((ell, em), {}).get("match", float("nan")) for (ell, em) in NR_MODES
    ]

    colors = [
        "C2"
        if (not np.isnan(m) and m >= 0.99)
        else "C3"
        if (not np.isnan(m) and m >= 0.95)
        else "C4"
        if not np.isnan(m)
        else "lightgray"
        for m in match_vals
    ]

    bars = ax.bar(
        mode_labels,
        [m if not np.isnan(m) else 0 for m in match_vals],
        color=colors,
        edgecolor="k",
        linewidth=0.8,
    )

    for bar, mv in zip(bars, match_vals):
        if not np.isnan(mv):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.002,
                f"{mv:.4f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
        else:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                0.5,
                "N/A",
                ha="center",
                va="center",
                fontsize=8,
                color="gray",
            )

    ax.axhline(0.99, color="C2", ls=":", lw=1, alpha=0.6, label="0.99")
    ax.axhline(0.95, color="C3", ls=":", lw=1, alpha=0.6, label="0.95")
    ax.set_ylim(0, 1.04)
    ax.set_ylabel("Match")
    ax.set_title(f"Per-mode match (NR vs NRSur7dq4, PSD: {psd_name})")
    ax.legend(fontsize=8, loc="lower right")

    plt.tight_layout()
    safe_id = _safe_sim_id(sim_name)
    fig_path = os.path.join(outdir, f"{safe_id}_mode_matches.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"      Figure written: {fig_path}")


def _print_table(results, sim_name):
    w = 74
    print(f"\n{'─'*w}")
    print(f"  Mode-by-mode comparison: {sim_name} vs NRSur7dq4")
    print(f"{'─'*w}")
    print(
        f"  {'Mode':<10} {'Match':>10}  {'ΔΦ/cycle [rad]':>16}  {'N_cycles':>10}  {'f_lower_mode':>14}"
    )
    print(f"{'─'*w}")
    for (ell, em) in NR_MODES:
        res = results.get((ell, em), {})
        mm = res.get("match", float("nan"))
        fl = res.get("f_lower_mode", float("nan"))
        dp = res.get("phase_diff_per_cycle", float("nan"))
        nc = res.get("n_cycles", float("nan"))
        mm_str = f"{mm:.6f}" if not np.isnan(mm) else "    N/A  "
        dp_str = f"{dp:.4f}" if not np.isnan(dp) else "           N/A"
        nc_str = f"{nc:.1f}" if not np.isnan(nc) else "       N/A"
        fl_str = f"{fl:.1f} Hz" if not np.isnan(fl) else "    N/A"
        print(
            f"  ({ell},{em:+d})      {mm_str:>10}  {dp_str:>16}  {nc_str:>10}  {fl_str:>14}"
        )
    print(f"{'─'*w}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────


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
        default=os.path.join(_SCRIPTS_DIR, "results"),
        help="Output directory for CSV and figure",
    )
    p.add_argument(
        "--rotate",
        action="store_true",
        help="Also compute SO(3)-rotation-optimized match (Step 2, slow)",
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
        delta_t=args.delta_t,
        rotate=args.rotate,
    )
