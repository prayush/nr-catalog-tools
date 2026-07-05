---
title: nrcatalogtools.comparisons
parent: API Reference
nav_order: 11
---

<!-- GENERATED FILE — DO NOT EDIT. Regenerate with `python bin/generate_api_docs.py`. -->
{% raw %}

# `nrcatalogtools.comparisons`

NR vs NRSur7dq4 comparison utilities.

``compare_sim_vs_surrogate`` runs the full per-mode match pipeline for a
single simulation: loads the catalog waveform, generates surrogate modes,
computes noise-weighted matches and phase-drift metrics, writes a CSV, and
saves a figure.


## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Constants

| Name | Value |
|---|---|
| `DELTA_T` | `1.0 / 4096` |
| `DISTANCE` | `1.0` |

---

### `compare_sim_vs_surrogate`

```python
compare_sim_vs_surrogate(catalog_name: str, sim_name: str, total_mass: float = 40.0, psd_name: str = 'aLIGOZeroDetHighPower', outdir: str | None = None, figsdir: str | None = None, delta_t: float = DELTA_T, rotate: bool = False) -> dict
```

Run the full NR vs NRSur7dq4 comparison for one simulation.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `catalog_name` | `str` | One of ``'SXS'``, ``'RIT'``, ``'MAYA'``. |
| `sim_name` | `str` | Simulation identifier (e.g. ``'SXS:BBH:0001'``). |
| `total_mass` | `float` | Total mass in solar masses (default 40). |
| `psd_name` | `str` | PyCBC analytic PSD name (default ``'aLIGOZeroDetHighPower'``). |
| `outdir` | `str` | Directory for the output CSV (default ``'results'`` under cwd). |
| `figsdir` | `str` | Directory for the output figure (default ``'figs'`` under cwd). |
| `delta_t` | `float` | Sample spacing in physical seconds (default 1/4096). |
| `rotate` | `bool` | Also compute the SO(3)-rotation-optimized match (slow). |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `dict` | ``{(ell, em): {'match', 'f_lower_mode', 'phase_diff_per_cycle', 'n_cycles', 'match_rotated'}}`` |

---

{% endraw %}
