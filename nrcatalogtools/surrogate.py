"""Surrogate waveform utilities for NR vs NRSur7dq4 comparisons.

Handles loading the surrogate, calling it with NR-compatible parameters,
and wrapping the output as physical-unit pycbc TimeSeries objects.

NRSur7dq4 notes
---------------
* Precessing model — takes full 3-vector spins, returns all modes up to ell=4.
* The ``ellMax`` parameter caps the output (max is 4); ``mode_list`` is not
  supported for precessing models.
* ``f_low`` controls waveform truncation only (start frequency); for NRSur7dq4
  the recommended value is 0 (return the entire waveform).
* ``f_ref`` sets the reference epoch at which the spins are defined; it is in
  **cycles/M** (= M * f_GW where M is in seconds).  For the NR–surrogate
  comparison the correct choice is ``f_ref = f_lower_NR * M_seconds``, which
  aligns the surrogate spin epoch with the NR relaxation-time spin values.
* ``dt`` argument is in dimensionless M units.
* Output ``h[(ell,em)]`` is a complex numpy array representing the spin-weight
  -2 spherical harmonic mode ``h_lm`` (same convention as SXS/RIT/MAYA NR
  data: ``h_+ - i h_×`` decomposed as ``Σ h_lm Y_{-2,lm}``).
* Time array ``t_sur`` is in dimensionless M units; ``t_sur[-1] ≈ +100 M``
  is near merger (peak amplitude).

Spin epoch conventions
----------------------
Two cases arise depending on whether the NR waveform starts before or after
the surrogate's minimum training frequency (~0.0165 M·Ω_orbital):

1. **NR shorter than surrogate** (``f_lower_NR > f_min_sur``):
   Pass ``f_low=0`` (full waveform) and ``f_ref=f_lower_NR_dimless``.  The
   surrogate backward-evolves the spins from the NR epoch to its start.

2. **NR longer than surrogate** (``f_lower_NR < f_min_sur``):
   The surrogate cannot extrapolate before its minimum frequency, so
   ``f_ref`` is clipped to ``f_min_sur``.  For **aligned-spin or non-spinning
   systems** the spin components are constant (no precession), so the metadata
   spins remain valid regardless of epoch.  For **precessing systems** one
   would need to extract the instantaneous spins from the NR dynamics at
   ``f_min_sur`` — not yet implemented; a warning is emitted in that case.
"""

from __future__ import annotations

import numpy as np
from pycbc.types import TimeSeries

from . import utils

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


# Full set of NR modes to include in the comparison table.  Modes absent from
# the surrogate (e.g. (5,5)) will appear with match=NaN in the output.
NR_MODES = [(2, 2), (2, 1), (3, 3), (4, 4), (5, 5), (3, 2), (4, 3)]

# Positive-m modes produced by NRSur7dq4 (ellMax=4); subset of NR_MODES.
SURROGATE_MODES = [(ell, em) for ell, em in NR_MODES if ell <= 4]

# Approximate surrogate minimum frequency in cycles/M (= M * f_GW in seconds).
# Derived from the known NRSur7dq4 minimum M*Omega_orbital ≈ 0.0165:
#   f_min_cycles_per_M = omega_min / pi ≈ 0.00525
_SURROGATE_F_MIN_CYCLES_PER_M = 0.01717 / np.pi  # slightly above omega_0 ≈ 0.0165


def generate_surrogate_modes(
    params: dict,
    total_mass: float,
    distance: float = 1.0,
    delta_t_seconds: float = 1.0 / 4096,
) -> tuple[dict, float]:
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
    tuple[dict, float]
        ``({(ell, em): pycbc.types.TimeSeries}, f_lower_effective)`` —
        complex mode time series in physical units with epoch at peak (2,2)
        amplitude, plus the effective starting GW frequency in Hz.

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

    # f_ref in cycles/M (= M * f_GW_hz, where M is in seconds).
    f_ref_dimless = f_lower_hz * m_secs  # cycles/M
    dt_dimless = delta_t_seconds / m_secs  # dimensionless time step

    sur = load_nrsur7dq4()

    chi1_perp = np.sqrt(chiA[0] ** 2 + chiA[1] ** 2)
    chi2_perp = np.sqrt(chiB[0] ** 2 + chiB[1] ** 2)
    is_precessing = (chi1_perp > 1e-4) or (chi2_perp > 1e-4)

    try:
        t_sur, h_sur, _ = sur(
            q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=0, f_ref=f_ref_dimless
        )
    except Exception as exc:
        err_str = str(exc)
        if (
            "too small" in err_str.lower()
            or "omega_ref" in err_str
            or "omega_0" in err_str
        ):
            import re

            m_match = re.search(r"([0-9.]+)\s*=\s*omega_0", err_str)
            if m_match:
                omega_min = float(m_match.group(1)) * 1.01
            else:
                omega_min = 0.0170
            f_ref_clipped = omega_min / np.pi

            if is_precessing:
                print(
                    f"      WARNING: f_ref={f_ref_dimless:.5f} cycles/M is below "
                    f"surrogate minimum; clipping to {f_ref_clipped:.5f}.  "
                    "For precessing systems the NR spins should be extracted from "
                    "NR dynamics at f_min_sur — using metadata spins instead "
                    "(may introduce a spin-epoch error)."
                )
            t_sur, h_sur, _ = sur(
                q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=0, f_ref=f_ref_clipped
            )
        else:
            raise

    amp_scale = utils.amp_to_physical(total_mass, distance)
    t_physical = t_sur * m_secs

    peak_idx = int(np.argmax(np.abs(h_sur[(2, 2)])))
    peak_time_phys = t_physical[peak_idx]
    epoch = t_physical[0] - peak_time_phys

    phase22 = np.unwrap(np.angle(h_sur[(2, 2)]))
    omega_gw_start = abs(phase22[1] - phase22[0]) / dt_dimless
    f_lower_effective = omega_gw_start / (2.0 * np.pi) / m_secs

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


def surrogate_dict_to_waveform_modes(h_sur_dict: dict, ell_max: int = 4):
    """Wrap a dictionary of pycbc mode TimeSeries into a WaveformModes object.

    Parameters
    ----------
    h_sur_dict : dict
        Dict mapping (ell, em) -> pycbc.TimeSeries.
    ell_max : int, optional
        Maximum ell value to support (default 4).

    Returns
    -------
    WaveformModes
    """
    from nrcatalogtools.waveform.modes import WaveformModes

    if not h_sur_dict:
        raise ValueError("h_sur_dict is empty")

    class SurrogateWaveformModes(WaveformModes):
        def __init__(self, h_dict, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._h_dict = h_dict

        def get_mode(self, ell, em, **kwargs):
            return self._h_dict[(ell, em)]

    # Extract time array from any mode
    any_mode = next(iter(h_sur_dict.values()))
    time_array = np.array(any_mode.sample_times)

    # Calculate number of modes for ell_min=2
    num_modes = (ell_max + 1) ** 2 - 4
    data = np.zeros((len(time_array), num_modes), dtype=complex)

    wfm = SurrogateWaveformModes(
        h_sur_dict, data, time=time_array, ell_min=2, ell_max=ell_max
    )

    for (ell, m), ts in h_sur_dict.items():
        if ell > ell_max or ell < 2:
            continue
        try:
            idx = wfm.index(ell, m)
            wfm.data[:, idx] = ts.data
        except ValueError:
            pass  # Ignore modes that are not within the index range

    return wfm
