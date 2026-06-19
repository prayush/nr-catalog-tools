#!/usr/bin/env python
"""Plot NR vs NRSur7dq4 mode-by-mode waveform comparison.

For precessing SXS systems the NR modes are frame-rotated (via Wigner
D-matrices) into the surrogate's inertial reference frame before plotting,
so both curves live in the same physical frame:

    Surrogate inertial frame:  ẑ = L̂(t_ref),  x̂ = n̂(t_ref)

The rotation matrix R is the one computed by Phase 2 (_epoch_align_spins).
Applying D^l_{mm'}(R) to the NR modes maps them from the NR simulation
coordinate frame to the surrogate's inertial frame.

Usage
-----
    python plot_mode_comparison.py --catalog SXS --sim SXS:BBH:0053 \\
        --total-mass 40 --outdir project/figs
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from nrcatalogtools import load_catalog
from nrcatalogtools.surrogate import (
    generate_surrogate_modes,
    check_surrogate_prior,
    NR_MODES,
    _epoch_align_spins,
    load_nrsur7dq4,
)
from nrcatalogtools import utils as nrutils

DISTANCE = 1.0  # Mpc
_ELLS = [2, 3, 4]


def _safe_sim_id(s: str) -> str:
    return s.replace(":", "_").replace("/", "_")


def _to_numpy_complex(pycbc_ts) -> tuple[np.ndarray, np.ndarray]:
    t = np.array(pycbc_ts.sample_times)
    h = np.array(pycbc_ts)
    return t, h


def _get_phase2_rotation(
    sim_name: str, params: dict, total_mass: float
) -> np.ndarray | None:
    """Return the final Phase 2 rotation matrix R (accounts for f_ref clipping).

    Replicates the same two-pass logic as generate_surrogate_modes so that R
    corresponds exactly to the epoch used when the surrogate was evaluated.
    Returns None if Phase 2 is not applicable (non-SXS or aligned-spin).
    """
    chi1_perp = np.sqrt(params["spin1x"] ** 2 + params["spin1y"] ** 2)
    chi2_perp = np.sqrt(params["spin2x"] ** 2 + params["spin2y"] ** 2)
    is_precessing = (chi1_perp > 1e-4) or (chi2_perp > 1e-4)
    if not is_precessing:
        return None

    import sxs as _sxs

    sim_obj = _sxs.load(sim_name, auto_supersede=True, download=False)
    m_secs = nrutils.time_to_physical(total_mass)
    dt_dimless = (1.0 / 4096) / m_secs

    sur = load_nrsur7dq4()
    q = params["mass1"] / params["mass2"]

    # --- first attempt at t_peak - 4500M ---
    try:
        chiA, chiB, f_ref, R = _epoch_align_spins(
            sim_obj, target_t_before_merger=4500.0
        )
    except Exception as e:
        print(f"  [Phase 2] initial epoch failed: {e}; R unavailable")
        return None

    # --- test the surrogate call to detect f_ref clip ---
    try:
        sur(q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=0, f_ref=f_ref)
    except Exception as exc:
        err_str = str(exc)
        if (
            "too small" in err_str.lower()
            or "omega_ref" in err_str
            or "omega_0" in err_str
        ):
            import re

            m_match = re.search(r"([0-9.]+)\s*=\s*omega_0", err_str)
            omega_min = float(m_match.group(1)) * 1.01 if m_match else 0.0170
            f_ref_clipped = omega_min / np.pi
            try:
                chiA, chiB, f_ref, R = _epoch_align_spins(
                    sim_obj, f_ref_target=f_ref_clipped
                )
            except Exception as e2:
                print(f"  [Phase 2] clipped epoch failed: {e2}; using first R")

    err = np.max(np.abs(R @ R.T - np.eye(3)))
    if err > 1e-6:
        print(f"  [Phase 2] R not orthogonal (err={err:.2e}); skipping rotation")
        return None

    print("\n  Phase 2 rotation matrix R (NR frame → surrogate frame):")
    print(f"    R = {np.round(R, 4)}")
    print(f"    det(R) = {np.linalg.det(R):.6f}  (should be +1)")
    return R


def rotate_nr_modes(h_nr_dict: dict, R: np.ndarray, ell_max: int = 4) -> dict:
    """Rotate NR modes from the NR frame into the surrogate's inertial frame.

    Applies the Wigner D-matrix corresponding to rotation R:
        h_lm^(sur) = Σ_{m'} D^l_{mm'}(R) · h_lm'^(NR)

    where R is the 3×3 passive rotation matrix (rows = new basis vectors in
    old basis), i.e. the Phase 2 matrix with rows [n̂, λ̂, L̂].

    Only positive-m modes present in h_nr_dict are returned (negative-m modes
    are filled from the h_{l,-m} = (-1)^m conj(h_{l,m}) symmetry, which holds
    in the NR frame for waveforms extracted via the Regge-Wheeler-Zerilli
    formalism that SXS uses).
    """
    import quaternionic
    import spherical

    # R is a passive rotation: to rotate modes we use the active equivalent R^T
    q_active = quaternionic.array.from_rotation_matrix(R.T)

    wig = spherical.Wigner(ell_max=ell_max)
    D = wig.D(q_active)  # D^l_{mm'}(R^T): shape (n_D_elements,)

    # Determine time grid from any available mode
    n_t = None
    t_ref = None
    for ts in h_nr_dict.values():
        t_ref, _ = _to_numpy_complex(ts)
        n_t = len(t_ref)
        break

    result = {}
    for ell in range(2, ell_max + 1):
        # Assemble all 2l+1 m-modes into h_in[t, m_idx]
        h_in = np.zeros((n_t, 2 * ell + 1), dtype=complex)
        for m_idx, m in enumerate(range(-ell, ell + 1)):
            if (ell, m) in h_nr_dict:
                _, h_arr = _to_numpy_complex(h_nr_dict[(ell, m)])
                h_in[:, m_idx] = h_arr
            elif (ell, -m) in h_nr_dict and m < 0:
                # Conjugation symmetry: h_{l,-m} = (-1)^m conj(h_{l,m})
                _, h_pos = _to_numpy_complex(h_nr_dict[(ell, -m)])
                h_in[:, m_idx] = ((-1) ** m) * np.conj(h_pos)

        # Rotate: h_out[t, m] = Σ_{m'} D^l_{mm'} h_in[t, m']
        h_out = np.zeros_like(h_in)
        for m_idx, m in enumerate(range(-ell, ell + 1)):
            for mp_idx, mp in enumerate(range(-ell, ell + 1)):
                d_val = D[wig.Dindex(ell, m, mp)]
                h_out[:, m_idx] += d_val * h_in[:, mp_idx]

        # Store only positive-m modes that were present in the original dict
        for m_idx, m in enumerate(range(-ell, ell + 1)):
            if (ell, m) in h_nr_dict:
                result[(ell, m)] = (t_ref, h_out[:, m_idx])

    return result


def _find_optimal_shift(
    h_nr: np.ndarray,
    h_sur: np.ndarray,
    delta_t: float,
    f_lower_hz: float,
    max_shift_s: float = 0.2,
) -> tuple[float, float]:
    """Return (t_shift_s, phi_shift_rad) maximising <Re(h_nr)|Re(e^{iφ} h_sur(t-τ))>.

    Uses a flat (white-noise) PSD with a hard low-frequency cutoff at f_lower_hz.
    Both arrays are assumed to be complex; only Re(h_nr) enters the objective.
    The search is limited to |τ| ≤ max_shift_s to suppress spurious wrap-around peaks.
    """
    n = max(len(h_nr), len(h_sur))
    n_fft = 1
    while n_fft < 2 * n:
        n_fft <<= 1

    freqs = np.fft.fftfreq(n_fft, d=delta_t)
    weight = (np.abs(freqs) >= f_lower_hz).astype(float)

    H1 = np.fft.fft(h_nr.real, n=n_fft)
    H2 = np.fft.fft(h_sur, n=n_fft)

    G = np.conj(H1) * H2 * weight
    C = np.fft.ifft(G)

    # Normalise so that max|C| = 1 for identical signals
    norm = np.sqrt(np.sum(np.abs(H1) ** 2 * weight) * np.sum(np.abs(H2) ** 2 * weight))
    if norm > 0:
        C /= norm

    # Restrict search to ±max_shift_s to avoid circular wrap-around peaks
    max_lag = min(int(np.ceil(max_shift_s / delta_t)), n_fft // 2)
    # Build an index array: [0, 1, ..., max_lag, n_fft-max_lag, ..., n_fft-1]
    pos_lags = np.arange(0, max_lag + 1)
    neg_lags = np.arange(n_fft - max_lag, n_fft)
    search_idx = np.concatenate([pos_lags, neg_lags])

    abs_C_search = np.abs(C[search_idx])
    best_in_search = int(np.argmax(abs_C_search))
    peak_idx = int(search_idx[best_in_search])

    t_shift = (peak_idx if peak_idx <= n_fft // 2 else peak_idx - n_fft) * delta_t
    phi_shift = float(np.angle(C[peak_idx]))

    return t_shift, phi_shift


def plot_mode_comparison(
    catalog_name: str,
    sim_name: str,
    total_mass: float = 40.0,
    delta_t: float = 1.0 / 4096,
    outdir: str = ".",
    t_window_s: tuple[float, float] = (-0.4, 0.06),
) -> None:
    os.makedirs(outdir, exist_ok=True)

    # --- load catalog & waveform ---
    print(f"Loading {catalog_name} catalog / {sim_name}...")
    cat = load_catalog(catalog_name)
    wfm = cat.get(sim_name)
    params = cat.get_parameters(sim_name, total_mass=total_mass)

    q = params["mass1"] / params["mass2"]
    chi1 = np.sqrt(
        params["spin1x"] ** 2 + params["spin1y"] ** 2 + params["spin1z"] ** 2
    )
    chi2 = np.sqrt(
        params["spin2x"] ** 2 + params["spin2y"] ** 2 + params["spin2z"] ** 2
    )

    if not check_surrogate_prior(params):
        print(f"  WARNING: outside NRSur7dq4 prior (q={q:.2f}, |χ₁|={chi1:.2f})")

    # --- generate surrogate (includes Phase 2 for precessing SXS) ---
    print("Generating NRSur7dq4 modes...")
    h_sur, f_lower_sur = generate_surrogate_modes(
        params,
        total_mass,
        DISTANCE,
        delta_t_seconds=delta_t,
        sim_name=sim_name,
        catalog=cat,
        nr_wfm=wfm,
    )

    # --- get Phase 2 rotation matrix (same R that was used above) ---
    R = None
    is_sxs = getattr(cat, "CATALOG_TYPE", None) == "SXS"
    if is_sxs:
        print("Computing Phase 2 rotation matrix R...")
        R = _get_phase2_rotation(sim_name, params, total_mass)

    # --- load NR modes and (if Phase 2) rotate into surrogate frame ---
    available = [(l, m) for (l, m) in NR_MODES if (l, m) in h_sur]
    h_nr_raw = {}
    for (ell, em) in available:
        try:
            ts = wfm.get_mode(
                ell,
                em,
                total_mass=total_mass,
                distance=DISTANCE,
                delta_t_seconds=delta_t,
            )
            h_nr_raw[(ell, em)] = ts
        except Exception:
            pass

    if R is not None and h_nr_raw:
        print(
            "  Rotating NR modes into surrogate inertial frame via Wigner D-matrices..."
        )
        h_nr_rotated = rotate_nr_modes(h_nr_raw, R, ell_max=4)
        frame_label = "NR (frame-rotated)"
    else:
        h_nr_rotated = None
        frame_label = "NR"

    # --- build figure ---
    n_modes = len(available)
    n_cols = 2
    n_rows = (n_modes + 1) // n_cols

    fig = plt.figure(figsize=(14, 3.5 * n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, hspace=0.45, wspace=0.32)

    fig.suptitle(
        f"{catalog_name}: {sim_name}  —  NR vs NRSur7dq4\n"
        f"M = {total_mass} M☉   q = {q:.3f}   "
        f"|χ₁| = {chi1:.3f}   |χ₂| = {chi2:.3f}\n"
        + (
            "NR modes rotated by Phase 2 R into surrogate frame (ẑ=L̂, x̂=n̂ at t_ref)"
            if R is not None
            else "NR modes in native simulation frame (unrotated)"
        ),
        fontsize=11,
        y=1.02,
    )

    for idx, (ell, em) in enumerate(available):
        row, col = divmod(idx, n_cols)
        ax_amp = fig.add_subplot(gs[row, col])

        # --- NR mode (rotated or raw) ---
        if h_nr_rotated is not None and (ell, em) in h_nr_rotated:
            t_nr, h_nr = h_nr_rotated[(ell, em)]
            t_nr = t_nr.copy()  # avoid mutating the shared array in rotate_nr_modes
        elif (ell, em) in h_nr_raw:
            t_nr, h_nr = _to_numpy_complex(h_nr_raw[(ell, em)])
        else:
            ax_amp.set_title(f"({ell},{em:+d})  NR unavailable")
            continue

        # --- surrogate mode ---
        t_sur, h_sur_arr = _to_numpy_complex(h_sur[(ell, em)])

        # Align amplitude peaks to t = 0 for both
        t_nr -= t_nr[np.argmax(np.abs(h_nr))]
        t_sur -= t_sur[np.argmax(np.abs(h_sur_arr))]

        # --- Optimal time and phase shift (white-noise, f_lower cutoff) ---
        # Use the surrogate effective start frequency as the low-frequency cutoff.
        # A single cutoff is used for all modes (no per-mode scaling) to keep
        # the shift search well above the noisy low-frequency floor.
        f_lower_shift = max(f_lower_sur, params.get("f_lower", 20.0))
        t_shift, phi_shift = _find_optimal_shift(
            h_nr, h_sur_arr, delta_t, f_lower_shift
        )

        # Apply shifts to surrogate
        h_sur_opt = h_sur_arr * np.exp(1j * phi_shift)
        t_sur_opt = t_sur - t_shift

        # Window and masks
        t0, t1 = t_window_s
        mask_nr = (t_nr >= t0) & (t_nr <= t1)
        mask_sur = (t_sur_opt >= t0) & (t_sur_opt <= t1)

        peak_amp = np.max(np.abs(h_nr))
        if peak_amp == 0:
            peak_amp = 1.0

        ax_re = ax_amp.twinx()

        lns = []
        lns += ax_amp.plot(
            t_nr[mask_nr],
            np.abs(h_nr[mask_nr]) / peak_amp,
            color="C0",
            lw=1.8,
            label=f"|h| {frame_label}",
        )
        lns += ax_amp.plot(
            t_sur_opt[mask_sur],
            np.abs(h_sur_opt[mask_sur]) / peak_amp,
            color="C1",
            lw=1.8,
            ls="--",
            label="|h| Surrogate",
        )
        lns += ax_re.plot(
            t_nr[mask_nr],
            h_nr[mask_nr].real / peak_amp,
            color="C0",
            lw=0.9,
            alpha=0.4,
            label=f"Re {frame_label}",
        )
        lns += ax_re.plot(
            t_sur_opt[mask_sur],
            h_sur_opt[mask_sur].real / peak_amp,
            color="C1",
            lw=0.9,
            ls="--",
            alpha=0.4,
            label="Re Surrogate (opt-aligned)",
        )

        ax_amp.set_title(f"({ell},{em:+d})", fontsize=11, fontweight="bold")
        ax_amp.set_xlabel("t − t_peak  [s]", fontsize=8)
        ax_amp.set_ylabel("|h| / peak", fontsize=8)
        ax_re.set_ylabel("Re(h) / peak", fontsize=8, color="gray")
        ax_re.tick_params(axis="y", labelcolor="gray", labelsize=7)
        ax_amp.tick_params(axis="both", labelsize=7)
        ax_amp.axvline(0, color="k", lw=0.6, alpha=0.3)
        ax_amp.set_xlim(t0, t1)
        ax_amp.grid(True, alpha=0.3, ls=":", lw=0.6)

        # Shift annotation
        ax_amp.text(
            0.03,
            0.97,
            f"Δt = {t_shift * 1e3:.2f} ms\nΔφ = {np.degrees(phi_shift):.1f}°",
            transform=ax_amp.transAxes,
            fontsize=7,
            va="top",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85),
        )

        labels = [ln.get_label() for ln in lns]
        ax_amp.legend(lns, labels, fontsize=6.5, loc="upper right", ncol=2)

    if n_modes % n_cols:
        fig.add_subplot(gs[n_rows - 1, n_cols - 1]).set_visible(False)

    safe_id = _safe_sim_id(sim_name)
    out_path = os.path.join(outdir, f"{safe_id}_waveform_modes.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out_path}")


def _build_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--catalog", default="SXS")
    p.add_argument("--sim", default="SXS:BBH:0053")
    p.add_argument("--total-mass", type=float, default=40.0)
    p.add_argument("--delta-t", type=float, default=1.0 / 4096)
    p.add_argument(
        "--outdir", default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "figs"))
    )
    p.add_argument(
        "--t-start",
        type=float,
        default=-0.4,
        help="Window start relative to peak (seconds)",
    )
    p.add_argument(
        "--t-end",
        type=float,
        default=0.06,
        help="Window end relative to peak (seconds)",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    plot_mode_comparison(
        catalog_name=args.catalog,
        sim_name=args.sim,
        total_mass=args.total_mass,
        delta_t=args.delta_t,
        outdir=args.outdir,
        t_window_s=(args.t_start, args.t_end),
    )
