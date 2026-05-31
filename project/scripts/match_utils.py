"""Per-mode match computation helpers.

Wraps ``pycbc.filter.match()`` with the mode-specific frequency cutoff
``f_lower_mode = f_lower * |m| / 2`` (GW frequency ≈ |m| × orbital frequency)
and handles PSD loading, length alignment, and graceful NaN returns for
near-zero (degenerate) modes.
"""

from __future__ import annotations

import numpy as np
from pycbc.types import TimeSeries, FrequencySeries


def load_psd(
    f_lower: float,
    delta_t: float,
    waveform_length_seconds: float,
    psd_name: str = "aLIGOZeroDetHighPower",
) -> FrequencySeries:
    """Load a named analytic PSD sampled to match a waveform's frequency grid.

    Parameters
    ----------
    f_lower : float
        Low-frequency cutoff in Hz.
    delta_t : float
        Time step of the waveforms in seconds (sets the Nyquist limit).
    waveform_length_seconds : float
        Duration of the longest waveform in seconds (sets frequency resolution).
    psd_name : str, optional
        PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``).

    Returns
    -------
    pycbc.types.FrequencySeries
    """
    from pycbc.psd import from_string

    # Next power-of-two sample count for the FFT
    n_samples = int(waveform_length_seconds / delta_t)
    n_fft = 1
    while n_fft < n_samples:
        n_fft <<= 1

    delta_f = 1.0 / (n_fft * delta_t)
    length_f = n_fft // 2 + 1

    return from_string(psd_name, length_f, delta_f, low_freq_cutoff=f_lower)


def compute_mode_match(
    h_nr: TimeSeries,
    h_sur: TimeSeries,
    f_lower_mode: float,
    psd_name: str = "aLIGOZeroDetHighPower",
    f_upper: float | None = None,
) -> float:
    """Compute the noise-weighted match between one NR and one surrogate mode.

    Both inputs should be the *real part* of the complex strain mode
    (h₊ component), sampled at the same ``delta_t``.  The function pads to
    the next power-of-two, builds a PSD at the matching frequency resolution,
    and calls ``pycbc.filter.match()``.

    Parameters
    ----------
    h_nr : pycbc.types.TimeSeries
        Real-valued NR mode time series.
    h_sur : pycbc.types.TimeSeries
        Real-valued surrogate mode time series.
    f_lower_mode : float
        Low-frequency cutoff for this mode in Hz.
        Use ``f_lower * |m| / 2`` (GW frequency scales as |m| × f_orbital).
    psd_name : str, optional
        PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``).
    f_upper : float or None, optional
        Upper frequency cutoff in Hz (default: Nyquist).

    Returns
    -------
    float
        Match in [0, 1], or ``float('nan')`` if either waveform has zero norm.
    """
    from pycbc.filter import match as pycbc_match
    from pycbc.psd import from_string

    # Guard: skip degenerate (near-zero) modes.
    if (
        float(np.max(np.abs(np.array(h_nr)))) < 1e-50
        or float(np.max(np.abs(np.array(h_sur)))) < 1e-50
    ):
        return float("nan")

    t_start = max(float(h_nr.start_time), float(h_sur.start_time))
    t_end = min(float(h_nr.end_time), float(h_sur.end_time))

    if t_end <= t_start:
        return float("nan")

    # Time-slice both waveforms to the common time window
    h_nr_sliced = h_nr.time_slice(t_start, t_end)
    h_sur_sliced = h_sur.time_slice(t_start, t_end)

    # Apply a gentle taper to suppress spectral leakage from the hard cut
    from pycbc.waveform.utils import taper_timeseries

    h1_tapered = taper_timeseries(h_nr_sliced, tapermethod="startend", return_lal=False)
    h2_tapered = taper_timeseries(
        h_sur_sliced, tapermethod="startend", return_lal=False
    )

    # Pad to the next power-of-two >= max length.  PyCBC match() requires
    # the PSD delta_f to exactly equal the waveform frequency series delta_f,
    # so we must build the PSD *after* we know the final padded length.
    raw_len = max(len(h1_tapered), len(h2_tapered))
    n_fft = 1
    while n_fft < raw_len:
        n_fft <<= 1

    h1 = h1_tapered.copy()
    h1.resize(n_fft)
    h2 = h2_tapered.copy()
    h2.resize(n_fft)

    delta_f = 1.0 / (n_fft * h1.delta_t)
    length_f = n_fft // 2 + 1
    psd = from_string(psd_name, length_f, delta_f, low_freq_cutoff=f_lower_mode)

    mm, _ = pycbc_match(
        h1,
        h2,
        psd=psd,
        low_frequency_cutoff=f_lower_mode,
        high_frequency_cutoff=f_upper,
    )
    return float(mm)


def compute_phase_diff_per_cycle(
    h_nr: TimeSeries,
    h_sur: TimeSeries,
) -> tuple[float, float]:
    """Compute accumulated phase difference per GW cycle over the common window.

    Both inputs are the *complex* mode time series (h_lm = h+ - i h×).
    The two waveforms are trimmed to their shared time window (both should have
    epoch set so t=0 is at peak amplitude), then the total accumulated phase of
    each is computed from the unwrapped angle.

    The metric returned is::

        phase_diff_per_cycle = |ΔΦ_NR - ΔΦ_sur| / N_cycles_NR   [rad / cycle]

    where ``ΔΦ = |φ(t_end) - φ(t_start)|`` is the total phase evolved and
    ``N_cycles_NR = ΔΦ_NR / (2π)``.

    Parameters
    ----------
    h_nr : pycbc.types.TimeSeries
        Complex NR mode time series.
    h_sur : pycbc.types.TimeSeries
        Complex surrogate mode time series.

    Returns
    -------
    tuple[float, float]
        ``(phase_diff_per_cycle, n_cycles_nr)``.
        Returns ``(nan, nan)`` if either waveform has zero norm or the common
        window contains fewer than 2 samples.
    """
    arr_nr = np.array(h_nr)
    arr_sur = np.array(h_sur)

    if float(np.max(np.abs(arr_nr))) < 1e-50 or float(np.max(np.abs(arr_sur))) < 1e-50:
        return float("nan"), float("nan")

    dt = float(h_nr.delta_t)
    t_start = max(float(h_nr.start_time), float(h_sur.start_time))
    t_end = min(float(h_nr.end_time), float(h_sur.end_time))

    if t_end <= t_start:
        return float("nan"), float("nan")

    # Slice both to the common window.
    i_nr_s = max(0, int(round((t_start - float(h_nr.start_time)) / dt)))
    i_nr_e = min(len(arr_nr), int(round((t_end - float(h_nr.start_time)) / dt)) + 1)
    i_sur_s = max(0, int(round((t_start - float(h_sur.start_time)) / dt)))
    i_sur_e = min(len(arr_sur), int(round((t_end - float(h_sur.start_time)) / dt)) + 1)

    n = min(i_nr_e - i_nr_s, i_sur_e - i_sur_s)
    if n < 2:
        return float("nan"), float("nan")

    phi_nr = np.unwrap(np.angle(arr_nr[i_nr_s : i_nr_s + n]))
    phi_sur = np.unwrap(np.angle(arr_sur[i_sur_s : i_sur_s + n]))

    delta_phi_nr = abs(phi_nr[-1] - phi_nr[0])
    delta_phi_sur = abs(phi_sur[-1] - phi_sur[0])

    n_cycles_nr = delta_phi_nr / (2.0 * np.pi)
    if n_cycles_nr < 0.5:
        return float("nan"), float("nan")

    phase_diff_per_cycle = abs(delta_phi_nr - delta_phi_sur) / n_cycles_nr
    return float(phase_diff_per_cycle), float(n_cycles_nr)


def mode_f_lower(f_lower: float, em: int) -> float:
    """Return the GW frequency cutoff for mode (ell, m).

    GW frequency for the (ell, |m|) mode is approximately |m| times the
    orbital frequency: ``f_gw ≈ |m| * f_orbital = |m| * f_lower / 2``
    (since the (2,2) mode has ``f_gw = 2 * f_orbital``).

    Parameters
    ----------
    f_lower : float
        Orbital reference frequency in Hz (= half the (2,2) GW frequency).
    em : int
        Azimuthal mode number m.

    Returns
    -------
    float
        Mode-specific GW frequency cutoff in Hz.
    """
    return f_lower * abs(em) / 2.0 if em != 0 else f_lower
