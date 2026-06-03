# `WaveformModes` — Conceptual Guide

**Class:** `nrcatalogtools.WaveformModes`  
**Inherits from:** `sxs.WaveformModes` (ndarray subclass)

> **Auto-generated API reference**: For a complete method listing with signatures and
> docstrings, see [API Reference → waveform](api/waveform.md).

`WaveformModes` is the central data object returned by all catalog `get()` calls. It stores
complex gravitational-wave strain multipole modes $h_{\ell m}(t)$ in dimensionless NR units
and provides methods to convert to physical units, extract individual modes, sum to
polarizations, apply frame rotations, and compute mismatches.

---

## Construction

Three construction paths are available. In normal usage you will not call these directly —
use the catalog `get()` method instead.

### From HDF5 (RIT / MAYA waveform files)

```python
WaveformModes.load_from_h5(file_path, metadata=metadata_dict, verbosity=0)
```

Reads `amp_l{ell}_m{em}/X,Y` and `phase_l{ell}_m{em}/X,Y` datasets. Interpolates all modes
onto a common uniform time grid. Returns a `WaveformModes` with shape `(n_times, n_modes)`.

### From tar.gz (RIT psi4 files)

```python
WaveformModes.load_from_targz(file_path, metadata=metadata_dict, verbosity=0)
```

Reads ASCII `.asc` / `.dat` / `.txt` files (columns: `time`, `real`, `imag`) from inside
the archive. Missing modes are filled with zeros. Interpolates onto a uniform grid.

### Wrapping an `sxs.WaveformModes` (SXS catalog)

```python
WaveformModes(raw_sxs_obj.data, raw_sxs_obj.time, sim_metadata=metadata_dict, **meta)
```

Used internally by `SXSCatalog.get()`. The `sxs.WaveformModes.data` property may return a
memoryview; all arithmetic wraps it with `np.array(..., dtype=complex)`.

---

## Core usage patterns

### Getting a physically-scaled mode

`get_mode()` is the primary method for retrieving a scaled individual mode:

```python
mode22 = wfm.get_mode(
    2, 2,
    total_mass=60.0,          # M_sun
    distance=100.0,           # Mpc
    delta_t_seconds=1./4096,  # physical seconds (4096 Hz sampling)
)
# mode22 is a complex pycbc.types.TimeSeries
```

The amplitude is scaled by $G M_\text{tot} / (c^2 \cdot d_\text{Mpc} \cdot \text{Mpc})$ and
the time epoch is set so that $t=0$ coincides with the peak of the $(2,2)$ mode.

### Getting polarizations

`get_td_waveform()` sums over all modes with spin-weight $-2$ spherical harmonics:

$$H = h_+ + i h_\times = \sum_{\ell,m} {}^{-2}Y_{\ell m}(\iota, \phi_c) \, h_{\ell m}(t)$$

```python
pols = wfm.get_td_waveform(
    total_mass=40.,           # M_sun
    distance=100.,            # Mpc
    inclination=0.2,          # radians
    coa_phase=0.3,            # radians
    delta_t_seconds=1./4096,  # physical seconds
)
hp = pols.real()
hc = -1 * pols.imag()
```

> **Polarization convention:** Returns `conjugate(h)` so that `.real()` gives $h_+$ and
> `.imag()` gives $+h_\times$.  This differs from LAL convention where `imag() = -h_\times`.
> Pass `lal_convention=True` to get LAL-compatible output.

### Starting-frequency helpers

```python
# GW frequency at waveform start, normalized to 1 M_sun — divide by total_mass to get Hz
f_start_1msun = wfm.f_lower_at_1Msun()
f_start_hz = f_start_1msun / total_mass_msun

# GW frequency at the relaxation epoch
f_relax_hz = wfm.f_lower_at_relaxation(total_mass=60.0)
```

### Trimming to relaxation time

```python
wfm_trimmed = wfm.trim_to_relaxation_time(total_mass=60.0)
```

Reads the relaxation time from metadata (tries keys `'relaxed-time'`, `'relaxation_time'`,
`'reference_time'`).

---

## Sampling interval parameters

`get_mode()` and `get_td_waveform()` accept two explicit sampling parameters:

| Parameter | Interpretation |
|-----------|----------------|
| `delta_t_seconds` | **Physical seconds** (e.g. `1/4096` for 4096 Hz sampling) |
| `delta_t_Msun` | **Dimensionless M units** (NR native, e.g. `0.5` means $0.5\,GM/c^3$) |

The returned `TimeSeries.delta_t` is **always in physical seconds**.

> **Deprecated:** The old `delta_t` keyword argument (which inferred units from the magnitude)
> is deprecated. Use `delta_t_seconds` or `delta_t_Msun` explicitly.

---

## Frame rotation

### Rotating a `WaveformModes` object

```python
import quaternionic
R = quaternionic.array.from_euler_angles(alpha, beta, gamma)
wfm_rotated = wfm.rotated(R)
```

Applies Wigner D-matrix rotation: $h^{\text{rot}}_{\ell m}(t) = \sum_{m'} h_{\ell m'}(t) \, D^\ell_{m'm}(R)$

Uses the `spherical` package for Wigner D-matrix computation.

### Rotating a mode dict (surrogate / PyCBC format)

```python
from nrcatalogtools import apply_wigner_rotation_to_mode_dict

rotated_modes = apply_wigner_rotation_to_mode_dict(
    mode_dict,      # {(ell, m): pycbc.TimeSeries, ...}
    R,              # quaternionic rotation
    ell_max=4,
)
```

Useful for rotating surrogate model outputs into the NR frame for direct comparison.

---

## Mismatch methods

### `match_sphere_averaged(other, psd, f_lower, delta_t, ...)`

This method calculates the match (noise-weighted overlap) between two waveforms, integrated over the entire sphere of possible observer directions (sky-averaged) and maximized over coordinate time shift $t_c$, coalescence phase $\phi_c$, and active/passive $SO(3)$ rotation $R$.

#### Spherical Integration & Mode-by-Mode Equivalence
The overlap between the full multidimensional strain fields $h_1(t, \Omega)$ and $h_2(t, \Omega)$ over the sphere $S^2$ is defined as:
$$\mathcal{O}_{\text{sphere}}(h_1, h_2) = \frac{\int_{S^2} \langle h_1(t, \Omega) \mid h_2(t, \Omega) \rangle_t \, d\Omega}{\sqrt{\left[ \int_{S^2} \langle h_1(t, \Omega) \mid h_1(t, \Omega) \rangle_t \, d\Omega \right] \left[ \int_{S^2} \langle h_2(t, \Omega) \mid h_2(t, \Omega) \rangle_t \, d\Omega \right]}}$$

By expanding the strain fields in spin-weighted spherical harmonics ${}^{-2}Y_{\ell m}(\theta, \phi)$ and invoking their orthonormality:
$$\int_{S^2} {}^{-2}Y_{\ell m}^*(\Omega) \, {}^{-2}Y_{\ell' m'}(\Omega) \, d\Omega = \delta_{\ell \ell'} \, \delta_{m m'}$$
all cross-terms between different modes decouple and vanish. The integrated overlap simplifies to the sum of the mode-by-mode overlaps:
$$\int_{S^2} \langle h_1(t, \Omega) \mid h_2(t, \Omega) \rangle_t \, d\Omega = \sum_{\ell, m} \langle h_{1, \ell m} \mid h_{2, \ell m} \rangle_t$$

#### Extrinsic Optimization
To align the coordinate frames, the second waveform is transformed via:
$$h_{2, \ell m}^{\text{rot, shifted}}(t) = e^{-i m \phi_c} \sum_{m'=-\ell}^{\ell} h_{2, \ell m'}(t - t_c) \, D^{\ell}_{m' m}(R)$$
where $D^\ell_{m'm}(R)$ is the Wigner D-matrix. The returned mismatch $\mathcal{M}$ is defined as:
$$\mathcal{M} = 1 - \max_{t_c,\, \phi_c,\, R \in SO(3)} \left[ \frac{\sum_{\ell, m} \langle h_{1, \ell m} \mid h_{2, \ell m}^{\text{rot, shifted}}(t_c, \phi_c, R) \rangle_t}{\sqrt{\left( \sum_{\ell, m} \langle h_{1, \ell m} \mid h_{1, \ell m} \rangle_t \right) \left( \sum_{\ell, m} \langle h_{2, \ell m} \mid h_{2, \ell m} \rangle_t \right)}} \right]$$

```python
mismatch = wfm_a.match_sphere_averaged(
    wfm_b,
    psd=my_psd,
    f_lower=20.0,            # Hz
    delta_t_seconds=1./4096,
)
```

### `match_sphere_averaged_bms_maximized(other, psd, f_lower, j_max, ...)`

This extended method performs the optimization over the infinite-dimensional **BMS supertranslation** group at null infinity $\mathcal{I}^+$ in addition to the standard rigid rotations and shifts. 

Supertranslations correspond to direction-dependent retarded-time shifts:
$$u' = u - \alpha(\theta, \phi)$$
where the supertranslation field $\alpha(\theta, \phi)$ is decomposed into scalar spherical harmonics:
$$\alpha(\theta, \phi) = \sum_{j=0}^{j_{\text{max}}} \sum_{k=-j}^{j} \alpha_{jk} \, Y_{jk}(\theta, \phi)$$
For small supertranslations, the waveform modes undergo first-order mixing:
$$h'_{\ell m}(u) = h_{\ell m}(u) - \sum_{j=0}^{j_{\text{max}}} \sum_{k=-j}^{j} \sum_{p, q} \alpha_{jk} \, \mathcal{G}^{\ell m}_{jk,pq} \, \dot{h}_{pq}(u)$$
where $\mathcal{G}^{\ell m}_{jk,pq}$ are the spin-weighted Gaunt coefficients:
$$\mathcal{G}^{\ell m}_{jk,pq} = \int_{S^2} {}^{-2}Y_{\ell m}^*(\Omega) \, Y_{jk}(\Omega) \, {}^{-2}Y_{pq}(\Omega) \, d\Omega$$

This method optimizes $t_c$, $\phi_c$, $R \in SO(3)$, and the supertranslation coefficients $\alpha_{jk}$ (where $j=1$ shifts the coordinate origin to correct for center-of-mass drift, and $j \ge 2$ modes correct for proper supertranslations) using a Nelder-Mead simplex solver to find the global minimum mismatch.

> [!NOTE]
> The `scri` package is required to compute the spin-weighted Gaunt coupling coefficients $\mathcal{G}^{\ell m}_{jk,pq}$.

---

## Unit conventions

All data stored internally in `WaveformModes` uses **geometrized, mass-scaled dimensionless
units**:

| Quantity | Dimensionless unit |
|----------|--------------------|
| Time | $GM_\text{tot}/c^3$ |
| Amplitude ($r \, h_{\ell m}$) | $GM_\text{tot}/c^2$ |

Physical conversion factors (from [`nrcatalogtools.utils`](api/utils.md)):

```python
import lal

# Time: 1 NR M unit → seconds
m_secs = total_mass_msun * lal.MTSUN_SI

# Amplitude: dimensionless → strain at distance d_Mpc
amp_scale = (lal.G_SI * total_mass_msun * lal.MSUN_SI
             / (lal.C_SI**2 * d_Mpc * 1e6 * lal.PC_SI))
```

---

## Design notes

### Per-instance `_filepath` attribute

`_filepath` is extracted and stored as a per-instance attribute **before** the parent
`sxs.WaveformModes.__new__()` call. This prevents class-level sharing where loading a
second simulation would overwrite the first simulation's cached file path.

### Non-writable memoryview wrapping

`sxs.WaveformModes.data` may return a non-writable memoryview. All arithmetic in
`WaveformModes` wraps the data with `np.array(..., dtype=complex)` before any in-place
operations to avoid `ValueError: assignment destination is read-only`.
