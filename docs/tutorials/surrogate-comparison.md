---
title: NR vs surrogate comparison
parent: Tutorials
nav_order: 3
---

# Tutorial: NR vs Surrogate (NRSur7dq4) Comparison
{: .no_toc }

This tutorial validates the NRSur7dq4 surrogate model against a catalog simulation by
computing per-mode noise-weighted matches — the workflow used for waveform-model
calibration and accuracy studies.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Prerequisites

```bash
pip install nr-catalog-tools gwsurrogate matplotlib
```

{: .note }
> `gwsurrogate` is an optional dependency of `nr-catalog-tools`. The NRSur7dq4 model
> data is downloaded automatically on first use. NRSur7dq4 is valid for mass ratios
> $$q \le 4$$ and spin magnitudes $$|\chi_{1,2}| \le 0.8$$; use
> [`check_surrogate_prior()`](../api/surrogate.md) before comparing.

## The one-call pipeline

[`compare_sim_vs_surrogate()`](../api/comparisons.md) runs the full comparison for a
single simulation — loads the NR waveform, generates matching surrogate modes, computes
per-mode matches and phase-drift metrics, and writes a CSV plus a diagnostic figure:

```python
from nrcatalogtools.comparisons import compare_sim_vs_surrogate

results = compare_sim_vs_surrogate(
    catalog_name="SXS",
    sim_name="SXS:BBH:0001",
    total_mass=40.0,                     # M_sun
    psd_name="aLIGOZeroDetHighPower",    # any PyCBC analytic PSD
)

for (ell, em), r in results.items():
    print(f"({ell},{em}): match={r['match']:.5f}  "
          f"f_lower={r['f_lower_mode']:.1f} Hz  "
          f"dphi/cycle={r['phase_diff_per_cycle']}")
```

Outputs land in `./results/<sim_id>_mode_matches.csv` and `./figs/` (override with the
`outdir` / `figsdir` arguments). For **precessing** simulations, the pipeline
automatically enables the $$SO(3)$$-rotation-optimized match (`rotate=True`), which
aligns the surrogate frame with the NR frame before matching; you can also request it
explicitly for aligned-spin runs.

The rest of this tutorial builds the same comparison step by step, so you can customize
each stage.

## Step 1: Load the NR waveform and its parameters

```python
import nrcatalogtools as nrcat

cat = nrcat.load_catalog("SXS")
sim_name = "SXS:BBH:0001"

wfm = cat.get(sim_name)
total_mass = 40.0
params = cat.get_parameters(sim_name, total_mass=total_mass)
print(params)
# {'mass1': ..., 'mass2': ..., 'spin1x': ..., ..., 'f_lower': ...}
```

`get_parameters()` returns a PyCBC-compatible dict — exactly the format the surrogate
helpers consume.

## Step 2: Check the surrogate prior

```python
from nrcatalogtools.surrogate import check_surrogate_prior

assert check_surrogate_prior(params), \
    f"{sim_name} lies outside the NRSur7dq4 validity volume"
```

To pre-filter an entire catalog instead, use
`nrcat.filter_by_surrogate_prior` (see the [API index](../api/index.md)).

## Step 3: Generate surrogate modes

```python
from nrcatalogtools.surrogate import generate_surrogate_modes

delta_t = 1.0 / 4096  # seconds

h_sur, f_low_eff = generate_surrogate_modes(
    params,
    total_mass=total_mass,
    distance=1.0,           # amplitude-irrelevant for matches
    delta_t_seconds=delta_t,
    sim_name=sim_name,      # enables epoch-aligned spins for precessing SXS runs
    catalog=cat,
)
# h_sur: {(ell, em): pycbc.types.TimeSeries}, epoch at the (2,2) amplitude peak
```

## Step 4: Extract the matching NR modes

```python
from nrcatalogtools.surrogate import SURROGATE_MODES

h_nr = {
    (ell, em): wfm.get_mode(ell, em,
                            total_mass=total_mass,
                            distance=1.0,
                            delta_t_seconds=delta_t)
    for (ell, em) in SURROGATE_MODES
}
```

Both dictionaries now hold complex `pycbc.types.TimeSeries` in physical units with a
common epoch convention (peak of the (2,2) amplitude at $$t = 0$$).

## Step 5: Compute per-mode matches

Each mode $$(\ell, m)$$ oscillates at roughly $$m/2$$ times the orbital frequency, so
its lower frequency cutoff must be scaled accordingly before matching:

```python
from nrcatalogtools.waveform.matching import (
    compute_mode_match,
    compute_phase_diff_per_cycle,
    mode_f_lower,
)

f_lower = params["f_lower"]

for (ell, em) in SURROGATE_MODES:
    f_low_mode = mode_f_lower(f_lower, em)
    match = compute_mode_match(
        h_nr[(ell, em)], h_sur[(ell, em)],
        f_lower_mode=f_low_mode,
        psd_name="aLIGOZeroDetHighPower",
    )
    dphi, n_cycles = compute_phase_diff_per_cycle(h_nr[(ell, em)],
                                                  h_sur[(ell, em)])
    print(f"({ell},{em}): match={match:.5f}, "
          f"phase drift/cycle={dphi:.2e} rad over {n_cycles:.0f} cycles")
```

`compute_mode_match()` maximizes over relative time and phase shifts and weights by the
chosen detector noise PSD; a match of 1 means the surrogate reproduces the NR mode
perfectly up to time/phase/amplitude conventions.

## Step 6 (optional): Optimize over frame rotations

For precessing systems the NR and surrogate frames generally differ by an unknown
$$SO(3)$$ rotation. [`apply_wigner_rotation_to_mode_dict()`](../api/waveform_matching.md)
rotates a whole mode dictionary via Wigner D-matrices:

```python
import quaternionic
from nrcatalogtools.waveform.matching import apply_wigner_rotation_to_mode_dict

R = quaternionic.array.from_euler_angles(alpha, beta, gamma)
h_sur_rot = apply_wigner_rotation_to_mode_dict(h_sur, R, ell_max=4)
```

Maximizing the sphere-averaged match over $$(\alpha, \beta, \gamma)$$ is what
`compare_sim_vs_surrogate(..., rotate=True)` does internally (and enables automatically
for precessing runs); see [Scientific Goal](../goal.md) for the formalism.

## Interpreting the results

- **Match $$\gtrsim 0.999$$** for the dominant (2,2) mode is typical inside the
  surrogate's calibration region.
- **Higher modes** — (3,3), (4,4), and especially (2,1) — carry less power and show
  lower matches; their `f_lower_mode` scaling also exposes them to more noise-weighted
  low-frequency content.
- **Phase drift per cycle** separates secular dephasing (modeling error) from
  overall time/phase offsets, which the match already maximizes over.
- The **total mass matters**: it sets where the waveform sits in the detector band, so
  quote matches together with `total_mass` and the PSD used.
