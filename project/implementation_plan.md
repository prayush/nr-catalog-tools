# NR Catalog vs NRSur7dq4 — Cross-Catalog Comparison Plan

## Background

The science goal is to compare every simulation in the SXS, RIT, and MAYA/GT
catalogs against the `NRSur7dq4` surrogate waveform at *identical* intrinsic
parameters.  Because the surrogate interpolates across the SXS training set, we
expect SXS simulations to agree very well.  Deviations in RIT and MAYA will
reveal the combined effect of code-dependent numerical errors and frame
conventions — the key scientific result.

The comparison must account for three classes of nuisance parameters that differ
between catalogs and the surrogate (as derived in `docs/goal.md` and
`notebooks/updates_v2.ipynb`):

1. **Time and phase shift** — trivially handled by `pycbc.filter.match()`.
2. **SO(3) rotation** — different NR codes define the z-axis differently;
   handled by optimizing over a unit quaternion using
   `WaveformModes.rotated(R)` / `apply_wigner_rotation_to_mode_dict()`.
3. **BMS supertranslation** — handled by
   `WaveformModes.match_sphere_averaged_bms_maximized()`.

The plan below is phased so that **Step 1 is immediately runnable** and each
subsequent step builds on it.

---

## Existing infrastructure (already in `nrcatalogtools`)

| Component | Location | Key method/function |
|---|---|---|
| Catalog loading | `nrcatalogtools/{sxs,rit,maya}.py` | `XCatalog.load()` |
| Parameter extraction | `nrcatalogtools/metadata.py` | `get_source_parameters_from_metadata()` |
| Single-mode extraction | `waveform/modes.py` | `WaveformModes.get_mode()` |
| Per-mode match | `waveform/modes.py` | `WaveformModes.match_single_mode()` |
| SO(3)-optimized overlap | `waveform/modes.py` | `WaveformModes.match_sphere_averaged()` |
| BMS-optimized overlap | `waveform/modes.py` | `WaveformModes.match_sphere_averaged_bms_maximized()` |
| Wigner rotation of dict | `waveform/matching.py` | `apply_wigner_rotation_to_mode_dict()` |
| Time/amp scaling | `nrcatalogtools/utils.py` | `time_to_physical()`, `amp_to_physical()` |

There is also a relevant existing script:
[compare_catalogs_q1_nospin.py](file:///home/prayush/src/nr-catalog-tools/bin/compare_catalogs_q1_nospin.py)
— this does cross-catalog amplitude/phase comparisons *without* the surrogate.

---

## Step 1 — Single-simulation NR vs NRSur7dq4 comparison script

### Goal
A self-contained, runnable script `scripts/compare_one_sim_vs_surrogate.py`
that:
1. Accepts a catalog name + simulation ID (or picks a default).
2. Loads the NR waveform and extracts intrinsic parameters via
   `catalog.get_parameters()`.
3. Generates `NRSur7dq4` modes at the same parameters, rescaled to 40 Mₛᵤₙ.
4. Computes `pycbc.filter.match()` for the 7 requested modes individually.
5. Writes a CSV result and a comparison figure.

### Design

#### 3.1 Parameter flow

```
NR simulation
  └─ catalog.get_parameters(sim_name, total_mass=40.)
       └─ {'mass1', 'mass2', 'spin1x/y/z', 'spin2x/y/z', 'f_lower'}
  └─ catalog.get(sim_name)
       └─ WaveformModes (dimensionless NR data)
            └─ wfm.get_mode(ell, em, total_mass=40., distance=1.,
                            delta_t_seconds=1./4096)
                 └─ pycbc.TimeSeries (complex, physical units)

NRSur7dq4
  └─ sur(q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=f_low_sur)
       └─ dict {(ell,em): h_lm (complex array)} in dimensionless units
  └─ rescale t_sur → physical:  t_sur * M * lal.MTSUN_SI
  └─ rescale h_sur → physical:  h_sur * G*M*Msun / (c²*D*Mpc)
  └─ wrap each mode as pycbc.TimeSeries
```

**Critical subtleties:**

- `NRSur7dq4` uses the SXS convention for `q = m1/m2 ≥ 1` and
  dimensionless spin components at the *reference time* (default:
  `-4500 M` before merger).  Use `f_ref = f_lower` for the NR simulation's
  starting orbital frequency.
- **Inertial frame output**: When called with `inclination=None`, the `gwsurrogate` evaluation of `NRSur7dq4` internally generates the modes in the co-orbital/co-precessing frame, and then rotates them to the **inertial (source) frame** before returning the mode dictionary. Therefore, the returned modes are already in the inertial frame of the surrogate, NOT the coprecessing frame. A direct mode-by-mode match (Step 1) evaluates these inertial modes without additional frame rotation.
- `pycbc.filter.match()` maximises over time and phase simultaneously.  For
  the (ℓ,m) mode we use `f_lower_mode = f_lower * |m| / 2` (GW frequency
  scales as `|m| * f_orbital`).
- The PSD should default to `aLIGOZeroDetHighPower` (as in the notebook),
  selectable via CLI.

#### 3.2 Modes

```python
MODES = [(2,2), (2,1), (3,3), (4,4), (5,5), (3,2), (4,3)]
```

For NR modes with negative m: `h_{ℓ,-m} = (-1)^ℓ h*_{ℓ,m}` for
non-precessing systems; for precessing systems all modes are independent and
the surrogate provides them.

#### 3.3 Key code outline

```python
# --- Load NR ---
cat = nrcat.SXSCatalog.load(download=False)  # or RIT/MAYA
wfm = cat.get(sim_name)
params = wfm.get_parameters(total_mass=40.)   # PyCBC-native dict

# --- Generate surrogate ---
import gwsurrogate as gws
sur = gws.LoadSurrogate('NRSur7dq4')
q   = params['mass1'] / params['mass2']       # m1 >= m2
chiA = [params['spin1x'], params['spin1y'], params['spin1z']]
chiB = [params['spin2x'], params['spin2y'], params['spin2z']]
t_sur, h_sur, dyn_sur = sur(q, chiA, chiB,
                             M=40., dist_mpc=1.,
                             dt=1./4096,
                             f_low=params['f_lower'],
                             mode_list=MODES)

# --- PSD ---
from pycbc.psd import from_string
psd = from_string('aLIGOZeroDetHighPower',
                  length=length_f, delta_f=delta_f,
                  low_freq_cutoff=f_lower)

# --- Per-mode match ---
results = {}
for (ell, em) in MODES:
    h_nr  = wfm.get_mode(ell, em, total_mass=40.,
                          distance=1., delta_t_seconds=1./4096)
    h_sur_lm = # wrap h_sur[(ell,em)] as pycbc.TimeSeries
    
    # align lengths
    N = max(len(h_nr), len(h_sur_lm))
    h_nr.resize(N); h_sur_lm.resize(N)
    psd_copy = psd.copy(); psd_copy.resize(len(h_nr.to_frequencyseries()))
    
    f_low_mode = params['f_lower'] * abs(em) / 2.
    mtc, _ = pycbc.filter.match(h_nr.real(), h_sur_lm.real(),
                                 psd=psd_copy,
                                 low_frequency_cutoff=f_low_mode)
    results[(ell, em)] = mtc
```

#### 3.4 Output

- **Console**: tabular summary of all 7 matches.
- **CSV**: `results/<sim_id>_mode_matches.csv` with columns
  `sim_id, catalog, ell, em, match, f_lower_mode`.
- **Figure**: 3-panel plot (amplitude, phase, match table) per simulation.

### Files to create

| File | Purpose |
|---|---|
| `scripts/compare_one_sim_vs_surrogate.py` | Main runnable script |
| `scripts/surrogate_utils.py` | Helper: load surrogate, wrap modes as pycbc TimeSeries, rescale |
| `scripts/match_utils.py` | Helper: per-mode match with correct f_lower, PSD loading |
| `scripts/catalog_utils.py` | Helper: load any catalog, filter by parameter cuts |

---

## Step 2 — SO(3) frame-rotation optimized match

### Goal
Extend Step 1 to maximize the match not just over time/phase but over the full
SO(3) source-frame rotation, implementing the formalism in
`docs/goal.md §Source frame ambiguity`.

### Design

The `WaveformModes.match_sphere_averaged()` method already implements this
using Nelder-Mead over `(t_c, φ_c, α, β, γ)`.  However, it operates on
two `WaveformModes` objects.

**Plan:**
1. **Inertial Frame Alignment**: Although the surrogate modes are in an inertial frame, that frame is defined by the binary's initial orbital angular momentum and separation vector at the reference epoch. This reference frame generally differs from the NR simulation's inertial frame by an SO(3) rotation. Therefore, optimizing over the rotation is required to align the two inertial frames.
2. Add a helper in `surrogate_utils.py` that wraps the surrogate mode dict as
   a `WaveformModes` object (respecting the same `ell_min/ell_max` and time
   grid as the NR data).
3. Call `wfm_nr.match_sphere_averaged(wfm_sur, psd, f_lower)`.
4. Record the optimal rotation quaternion `R*` alongside the match.
5. Study `R*` as a function of simulation properties (eccentricity, spin
   magnitudes, simulation length) across the catalog.

### New files

| File | Purpose |
|---|---|
| `scripts/compare_one_sim_rotation_optimized.py` | Step 2 script |
| `nrcatalogtools/waveform/surrogate_wrapper.py` | Wrap surrogate dict → WaveformModes |

---

## Step 3 — Batch processing over full catalogs

### Goal
Run Steps 1 and 2 over **all** simulations in each catalog that fall within
the `NRSur7dq4` prior volume:

- `q ∈ [1, 4]`
- `|χ₁|, |χ₂| ≤ 0.8`
- Total length `≥ some minimum`, orbital eccentricity `< 0.05`

### Design

```
for catalog in [SXS, RIT, MAYA]:
    for sim in catalog.simulations_list:
        params = catalog.get_parameters(sim, total_mass=40.)
        if not in_surrogate_prior(params): continue
        result = compute_mode_matches(catalog, sim, params)
        append_to_hdf5(results_file, sim, catalog, result)
```

Use `multiprocessing.Pool` for embarrassingly parallel execution.

### New files

| File | Purpose |
|---|---|
| `scripts/batch_compare_catalog.py` | Batch loop with multiprocessing |
| `scripts/results_io.py` | HDF5 I/O for match results |

---

## Step 4 — BMS supertranslation correction and analysis

### Goal
For the subset of simulations where the SO(3)-optimized match is below a
threshold (e.g. < 0.99), apply the BMS supertranslation optimization from
`WaveformModes.match_sphere_averaged_bms_maximized()` and measure the
improvement.  This isolates gauge artifacts from genuine numerical errors.

### New files

| File | Purpose |
|---|---|
| `scripts/compare_one_sim_bms.py` | BMS-optimized match for one sim |
| `scripts/batch_compare_bms.py` | Batch BMS pass over flagged sims |

---

## Step 5 — Visualization and scientific analysis

### Goal
Produce the figures for the paper:

1. **Match vs parameter space** — scatter plots of match vs q, χ_eff, χ_p for
   each catalog, color-coded by mode.
2. **Optimal rotation distribution** — histogram/scatter of Euler angles of
   R* across the catalog; look for systematic offsets between RIT/MAYA vs SXS.
3. **Mismatch comparison** — CDF of mismatches across catalogs.
4. **Mode-by-mode breakdown** — radar/bar plot of per-mode matches.

### New files

| File | Purpose |
|---|---|
| `notebooks/catalog_comparison_results.ipynb` | Analysis notebook |
| `scripts/plot_results.py` | Publication-quality figures |

---

## Verification Plan

### Step 1 verification (immediate)
- Run on `SXS:BBH:0001` (q≈1, no spin) — expect match > 0.999 for (2,2).
- Run on `RIT:BBH:0001-n100-id3` — expect match > 0.99 for (2,2).
- Sanity check: match of NR against itself = 1.0.
- Check that `f_lower_mode = f_lower * |m| / 2` gives physically sensible
  cutoffs for each mode.

### Step 2 verification
- For a q=1 non-spinning system, the optimal rotation R* should be very close
  to identity (or a known symmetry rotation π about z).
- Compare `match_sphere_averaged` result with Step 1 result; rotation
  optimization should improve or maintain the match.

### Step 3 verification
- Spot-check a few simulations from Step 1 to confirm batch results agree.
- Check that the number of simulations processed matches the expected prior
  volume.

---

## Open Questions

> [!NOTE]
> **Surrogate mode convention**: Both `gwsurrogate` and `WaveformModes.get_mode()` return modes as complex arrays `h_lm = h_+ - i h_×` (spin-weight -2 spherical harmonics) in dimensionless `r*h/M` units. This has been verified, so no complex convention conversion is needed—just amplitude scaling and time conversion.

> [!NOTE]
> **Reference frequency alignment**: The `NRSur7dq4` evaluator defaults `f_ref` to `f_low`. The `nrcatalogtools` library correctly computes `f_lower` as the instantaneous frequency of the (2,2) mode at the NR simulation's relaxation time. Since this frequency is passed as `f_low` to the surrogate, the surrogate's reference epoch for spin definition perfectly aligns with the NR spin extraction epoch.

> [!NOTE]
> **NRSur7dq4 coverage**: The surrogate covers q ∈ [1,4], |χ| ≤ 0.8,
> and simulation lengths ≥ ~4350 M.  SXS simulations shorter than this will
> be excluded from the comparison.  RIT and MAYA simulations at higher q or
> spin will also be excluded.

> [!NOTE]
> **Negative-m modes**: `gwsurrogate` returns both positive and negative m
> modes.  The NR catalogs may store only positive m (using the symmetry
> `h_{ℓ,-m} = (-1)^ℓ h*_{ℓ,m}` for non-precessing systems) or may store all
> modes.  The `WaveformModes` class has `get_mode(ell, em)` for any m.
> Confirm availability of (2,1), (3,2), (4,3) modes in each catalog.

---

## Immediate Next Action

→ **Implement Step 1**: `scripts/compare_one_sim_vs_surrogate.py`

Please confirm:
1. Which catalog/simulation to use as the default test case?
   (Suggestion: `SXS:BBH:0001` for non-spinning, or a spinning RIT sim)
2. Should the script also include an optional SO(3) rotation flag `--rotate`
   to run Step 2 inline?
3. PSD: use `aLIGOZeroDetHighPower` or `flat` for an initial code-correctness
   check?
