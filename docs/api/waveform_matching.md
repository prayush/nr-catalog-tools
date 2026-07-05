---
title: nrcatalogtools.waveform.matching
parent: API Reference
nav_order: 7
---

<!-- GENERATED FILE — DO NOT EDIT. Regenerate with `python bin/generate_api_docs.py`. -->
{% raw %}

# `nrcatalogtools.waveform.matching`

Standalone waveform matching and rotation helpers.

These are module-level functions (not bound to WaveformModes) so they can
be unit-tested and used independently of the class.


## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### `apply_wigner_rotation_to_mode_dict`

```python
apply_wigner_rotation_to_mode_dict(mode_dict, R, ell_max=4)
```

Apply a Wigner rotation to a dictionary of spherical harmonic modes.

This is useful for rotating the output of ``gwsurrogate`` or
``pycbc.waveform.get_td_waveform_modes`` (which return dicts) into the
NR source frame before computing mode-by-mode matches.

The rotation is applied mode-by-mode via Wigner D-matrices:

    h'_{ℓm}(t) = Σ_{m'} D^{(ℓ)}_{m'm}(R) h_{ℓm'}(t)

where R ∈ SO(3) is a unit quaternion and D^{(ℓ)} is the (2ℓ+1)×(2ℓ+1)
Wigner D-matrix for angular momentum ℓ.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `mode_dict` | `dict` | Keys are ``(l, m)`` integer tuples; values are complex ``pycbc.types.TimeSeries`` objects (or 1-D numpy arrays of matching length). |
| `R` | `quaternionic.array` | Unit quaternion representing the rotation. |
| `ell_max` | `int` | Maximum ℓ to include (default 4). |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `dict` | Rotated mode dictionary with the same ``(l, m)`` keys. |

---

### `load_psd`

```python
load_psd(f_lower: float, delta_t: float, waveform_length_seconds: float, psd_name: str = 'aLIGOZeroDetHighPower')
```

Load a named analytic PSD sampled to match a waveform's frequency grid.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `f_lower` | `float` | Low-frequency cutoff in Hz. |
| `delta_t` | `float` | Time step of the waveforms in seconds (sets the Nyquist limit). |
| `waveform_length_seconds` | `float` | Duration of the longest waveform in seconds (sets frequency resolution). |
| `psd_name` | `str` | PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``). |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `pycbc.types.FrequencySeries` |  |

---

### `compute_mode_match`

```python
compute_mode_match(h_nr, h_sur, f_lower_mode: float, psd_name: str = 'aLIGOZeroDetHighPower', f_upper=None) -> float
```

Compute the noise-weighted match between one NR and one surrogate mode.

Both inputs should be the *real part* of the complex strain mode
(h₊ component), sampled at the same ``delta_t``.  The function pads to
the next power-of-two, builds a PSD at the matching frequency resolution,
and calls ``pycbc.filter.match()``.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `h_nr` | `pycbc.types.TimeSeries` | Real-valued NR mode time series. |
| `h_sur` | `pycbc.types.TimeSeries` | Real-valued surrogate mode time series. |
| `f_lower_mode` | `float` | Low-frequency cutoff for this mode in Hz. Use ``f_lower * \|m\| / 2`` (GW frequency scales as \|m\| × f_orbital). |
| `psd_name` | `str` | PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``). |
| `f_upper` | `float or None` | Upper frequency cutoff in Hz (default: Nyquist). |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `float` | Match in [0, 1], or ``float('nan')`` if either waveform has zero norm. |

---

### `compute_phase_diff_per_cycle`

```python
compute_phase_diff_per_cycle(h_nr, h_sur) -> tuple
```

Compute accumulated phase difference per GW cycle over the common window.

Both inputs are the *complex* mode time series (h_lm = h+ - i h×).
The two waveforms are trimmed to their shared time window (both should have
epoch set so t=0 is at peak amplitude), then the total accumulated phase of
each is computed from the unwrapped angle.

The metric returned is::

    phase_diff_per_cycle = |ΔΦ_NR - ΔΦ_sur| / N_cycles_NR   [rad / cycle]

where ``ΔΦ = |φ(t_end) - φ(t_start)|`` is the total phase evolved and
``N_cycles_NR = ΔΦ_NR / (2π)``.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `h_nr` | `pycbc.types.TimeSeries` | Complex NR mode time series. |
| `h_sur` | `pycbc.types.TimeSeries` | Complex surrogate mode time series. |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `tuple[float, float]` | ``(phase_diff_per_cycle, n_cycles_nr)``. Returns ``(nan, nan)`` if either waveform has zero norm or the common window contains fewer than 2 samples. |

---

### `mode_f_lower`

```python
mode_f_lower(f_lower: float, em: int) -> float
```

Return the GW frequency cutoff for mode (ell, m).

GW frequency for the (ell, |m|) mode is approximately |m| times the
orbital frequency: ``f_gw ≈ |m| * f_orbital = |m| * f_lower / 2``
(since the (2,2) mode has ``f_gw = 2 * f_orbital``).

#### Parameters

| Name | Type | Description |
|---|---|---|
| `f_lower` | `float` | GW frequency of the (2,2) mode in Hz (= 2 × orbital frequency). This is what ``CatalogBase.get_parameters()`` returns as ``f_lower``. |
| `em` | `int` | Azimuthal mode number m. For m=0 the mode carries no oscillatory GW power at a well-defined frequency; ``f_lower`` is returned as a conservative lower bound but the result should not be interpreted as a physically meaningful frequency cutoff for that mode. |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `float` | Mode-specific GW frequency cutoff in Hz. |

---

### `interpolate_in_amp_phase`

```python
interpolate_in_amp_phase(obj, new_time, k=3, kind=None)
```

Interpolate in amplitude and phase using a variety of methods.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `obj` | `sxs.TimeSeries` | Complex waveform time series. |
| `new_time` | `array_like` | New time axis to interpolate onto. |
| `k` | `int` | Spline order for ``InterpolatedUnivariateSpline`` (default 3). |
| `kind` | `str` | Alternative interpolation: ``'linear'``, ``'quadratic'``, ``'cubic'``, or ``'CubicSpline'``. When specified, ``k`` is ignored. |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `sxs.TimeSeries` | Interpolated complex waveform on ``new_time``. |

---

{% endraw %}
