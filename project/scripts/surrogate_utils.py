"""Surrogate waveform utilities for NR vs NRSur7dq4 comparisons.

Handles loading the surrogate, calling it with NR-compatible parameters,
and wrapping the output as physical-unit pycbc TimeSeries objects.

NRSur7dq4 notes
---------------
* Precessing model — takes full 3-vector spins, returns all modes up to ell=4.
* The ``ellMax`` parameter caps the output (max is 4); ``mode_list`` is not
  supported for precessing models.
* ``f_low`` argument is in dimensionless M·Ω units (orbital angular frequency),
  **not** Hz.  Convert via ``Momega = π · f_lower_hz · M · MTSUN_SI``.
* ``dt`` argument is in dimensionless M units.
* Output ``h[(ell,em)]`` is a complex numpy array representing the spin-weight
  -2 spherical harmonic mode ``h_lm`` (same convention as SXS/RIT/MAYA NR
  data: ``h_+ - i h_×`` decomposed as ``Σ h_lm Y_{-2,lm}``).
* Time array ``t_sur`` is in dimensionless M units; ``t_sur[-1] ≈ +100 M``
  is near merger (peak amplitude).
"""

from __future__ import annotations

import numpy as np
from pycbc.types import TimeSeries

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from nrcatalogtools import utils

# Module-level surrogate singleton to avoid reloading the large model file.
_nrsur7dq4 = None


def load_nrsur7dq4():
    """Load and cache the NRSur7dq4 surrogate model.

    Returns
    -------
    gwsurrogate surrogate object
    """
    global _nrsur7dq4
    if _nrsur7dq4 is None:
        import gwsurrogate as gws

        _nrsur7dq4 = gws.LoadSurrogate("NRSur7dq4")
    return _nrsur7dq4


# Modes available from NRSur7dq4 (ellMax ≤ 4; (5,5) is not available).
SURROGATE_MODES = [(2, 2), (2, 1), (3, 3), (4, 4), (3, 2), (4, 3)]


def generate_surrogate_modes(
    params: dict,
    total_mass: float,
    distance: float = 1.0,
    delta_t_seconds: float = 1.0 / 4096,
) -> dict:
    """Call NRSur7dq4 and return physical-unit modes as a pycbc TimeSeries dict.

    Parameters
    ----------
    params : dict
        PyCBC-compatible binary parameter dict as returned by
        ``CatalogBase.get_parameters()``.  Must contain:
        ``mass1``, ``mass2``, ``spin1x/y/z``, ``spin2x/y/z``, ``f_lower``.
    total_mass : float
        Total binary mass in solar masses (sets the physical time/frequency
        scale for the surrogate call).
    distance : float, optional
        Luminosity distance in Mpc for amplitude scaling (default 1).
    delta_t_seconds : float, optional
        Desired sample spacing in physical seconds (default 1/4096).

    Returns
    -------
    dict
        ``{(ell, em): pycbc.types.TimeSeries}`` — complex mode time series in
        physical units (same strain units as ``WaveformModes.get_mode()``),
        epoch set so that ``t=0`` coincides with the peak of the (2,2) mode.

    Raises
    ------
    ValueError
        If ``f_lower`` is not positive or the surrogate call fails.
    """
    q = params["mass1"] / params["mass2"]  # ≥ 1 by PyCBC convention
    chiA = [params["spin1x"], params["spin1y"], params["spin1z"]]
    chiB = [params["spin2x"], params["spin2y"], params["spin2z"]]
    f_lower_hz = params["f_lower"]

    if f_lower_hz <= 0:
        raise ValueError(
            f"f_lower = {f_lower_hz} Hz is not positive. "
            "Cannot determine surrogate start frequency."
        )

    m_secs = utils.time_to_physical(total_mass)  # M * MTSUN_SI  [seconds]

    # Dimensionless surrogate arguments.
    # gwsurrogate uses the "Mf" convention: f_low = M * f_GW where M is in
    # seconds (M_total * MTSUN_SI) and f_GW is the gravitational-wave frequency
    # in Hz.  f_lower_hz is the (2,2)-mode GW frequency, so f_low = m_secs * f_GW.
    # (NOT pi * m_secs * f_GW — that would be M*Omega_orbital, which is wrong.)
    f_low_sur = f_lower_hz * m_secs  # M * f_GW  (dimensionless Mf)
    dt_dimless = delta_t_seconds / m_secs  # dimensionless time step

    sur = load_nrsur7dq4()
    try:
        t_sur, h_sur, _ = sur(q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=f_low_sur)
    except Exception as exc:
        # f_low_sur may be below the surrogate's minimum starting orbital frequency
        # (NRSur7dq4 is limited to the extent of its training waveforms, ~4500 M).
        # Parse omega_min from the error message, convert to Mf, and retry.
        # Error format: "Got omega_ref = X < Y = omega_0, too small!"
        if "too small" in str(exc).lower() or "omega_ref" in str(exc):
            import re

            # Error format: "Got omega_ref = X < Y = omega_0, too small!"
            # We need Y (the minimum), which appears as "Y = omega_0".
            m = re.search(r"([0-9.]+)\s*=\s*omega_0", str(exc))
            if m:
                omega_min = float(m.group(1)) * 1.01  # 1 % above minimum
            else:
                omega_min = 0.0170  # safe fallback above known NRSur7dq4 bound
            f_low_clipped = omega_min / np.pi  # Mf = omega_orbital / pi
            t_sur, h_sur, _ = sur(
                q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=f_low_clipped
            )
        else:
            raise

    # Physical scaling factors
    amp_scale = utils.amp_to_physical(total_mass, distance)
    t_physical = t_sur * m_secs  # physical seconds (still relative to some epoch)

    # Set epoch so t=0 is at peak (2,2) amplitude — mirrors WaveformModes convention.
    peak_idx = int(np.argmax(np.abs(h_sur[(2, 2)])))
    peak_time_phys = t_physical[peak_idx]
    epoch = t_physical[0] - peak_time_phys

    # Compute the effective starting GW frequency from the (2,2) phase derivative.
    # This is the f_lower that callers should use in the match integral — it may be
    # higher than the requested f_lower_hz if the surrogate was clipped to its
    # minimum training extent.
    phase22 = np.unwrap(np.angle(h_sur[(2, 2)]))
    omega_gw_start = abs(phase22[1] - phase22[0]) / dt_dimless  # M * Omega_GW (dimless)
    f_lower_effective = omega_gw_start / (2.0 * np.pi) / m_secs  # Hz

    result = {}
    for (ell, em) in SURROGATE_MODES:
        if (ell, em) not in h_sur:
            continue
        h_phys = h_sur[(ell, em)] * amp_scale
        result[(ell, em)] = TimeSeries(
            h_phys.astype(np.complex128),
            delta_t=delta_t_seconds,
            epoch=epoch,
        )

    # Return both the mode dict and the effective starting frequency so callers can
    # set the match lower cutoff to max(f_lower_nr, f_lower_sur).
    return result, f_lower_effective


def check_surrogate_prior(
    params: dict, q_max: float = 4.0, chi_max: float = 0.8
) -> bool:
    """Return True if the binary parameters fall within the NRSur7dq4 prior volume.

    NRSur7dq4 is valid for:
    - ``q = m1/m2 ∈ [1, 4]``
    - ``|χ₁|, |χ₂| ≤ 0.8``
    - Waveform length ≥ ~4350 M (not checked here — handled by the surrogate call).

    Parameters
    ----------
    params : dict
        PyCBC-compatible parameter dict.
    q_max : float, optional
        Maximum mass ratio (default 4).
    chi_max : float, optional
        Maximum spin magnitude (default 0.8).

    Returns
    -------
    bool
    """
    q = params["mass1"] / params["mass2"]
    if not (1.0 <= q <= q_max):
        return False
    chi1 = np.sqrt(
        params["spin1x"] ** 2 + params["spin1y"] ** 2 + params["spin1z"] ** 2
    )
    chi2 = np.sqrt(
        params["spin2x"] ** 2 + params["spin2y"] ** 2 + params["spin2z"] ** 2
    )
    return chi1 <= chi_max and chi2 <= chi_max
