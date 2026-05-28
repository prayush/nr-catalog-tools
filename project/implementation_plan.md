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

## Simulation Classification

All simulations in the three catalogs are pre-classified into six mutually
exclusive categories based on spin morphology and orbital eccentricity.
The classification is implemented in
`nrcatalogtools/classification.py` (`NRCatalogClassifier`) and the results
are persisted as JSON files in `catalog_organization/`.

### Six categories

| Key | Full name | Spin condition | Eccentricity condition |
|---|---|---|---|
| **a** | non-spinning eccentric | $\|\boldsymbol{\chi}_{1,2}\| < \epsilon_\chi$ | $e \geq \epsilon_e$ |
| **b** | non-spinning non-eccentric | $\|\boldsymbol{\chi}_{1,2}\| < \epsilon_\chi$ | $e < \epsilon_e$ |
| **c** | aligned-spin eccentric | $\chi_\perp < \epsilon_\chi$, $|\chi_z| \geq \epsilon_\chi$ | $e \geq \epsilon_e$ |
| **d** | aligned-spin non-eccentric | $\chi_\perp < \epsilon_\chi$, $|\chi_z| \geq \epsilon_\chi$ | $e < \epsilon_e$ |
| **e** | precessing-spin eccentric | $\chi_\perp \geq \epsilon_\chi$ | $e \geq \epsilon_e$ |
| **f** | precessing-spin non-eccentric | $\chi_\perp \geq \epsilon_\chi$ | $e < \epsilon_e$ |

where $\chi_\perp^2 = \chi_x^2 + \chi_y^2$ is the in-plane spin magnitude for
either body.

### Thresholds

```python
spin_threshold  ε_χ = 0.001   # |component| below this → treated as zero
ecc_threshold   ε_e = 0.005   # eccentricity below this → treated as circular
```

These are the defaults in `NRCatalogClassifier.__init__()`.  They can be
overridden at instantiation time to explore threshold sensitivity.

### Metadata field mapping per catalog

The classifier reads different metadata keys depending on the catalog,
since each code stores spin and eccentricity information under distinct names:

| Quantity | SXS key | RIT key (preferred / fallback) | MAYA key |
|---|---|---|---|
| Eccentricity | `reference_eccentricity` | `eccentricity` | `eccentricity` |
| $\chi_{1x}$ | `reference_dimensionless_spin1[0]` | `relaxed-chi1x` / `initial-bh-chi1x` | `a1x` |
| $\chi_{1y}$ | `reference_dimensionless_spin1[1]` | `relaxed-chi1y` / `initial-bh-chi1y` | `a1y` |
| $\chi_{1z}$ | `reference_dimensionless_spin1[2]` | `relaxed-chi1z` / `initial-bh-chi1z` | `a1z` |
| $\chi_{2x}$ | `reference_dimensionless_spin2[0]` | `relaxed-chi2x` / `initial-bh-chi2x` | `a2x` |
| $\chi_{2y}$ | `reference_dimensionless_spin2[1]` | `relaxed-chi2y` / `initial-bh-chi2y` | `a2y` |
| $\chi_{2z}$ | `reference_dimensionless_spin2[2]` | `relaxed-chi2z` / `initial-bh-chi2z` | `a2z` |

For SXS, spin vectors are extracted at the **relaxation time** (the same epoch
as `f_lower`).  For RIT, the `relaxed-*` fields are preferred over
`initial-bh-*` fields because they correspond to a post-junk-radiation epoch.
Eccentricity strings of the form `"< 0.002"` or `"~0.01"` are cleaned by
stripping the `<` and `~` characters before conversion to float; `None` or
`NaN` values are treated as zero.

### Catalog counts by category

| Category | SXS | RIT | MAYA | **All** |
|---|---:|---:|---:|---:|
| (a) non-spinning eccentric | 206 | 499 | 74 | **779** |
| (b) non-spinning non-eccentric | 177 | 54 | 34 | **265** |
| (c) aligned-spin eccentric | 21 | 231 | 117 | **369** |
| (d) aligned-spin non-eccentric | 687 | 541 | 40 | **1,268** |
| (e) precessing-spin eccentric | 30 | 117 | 303 | **450** |
| (f) precessing-spin non-eccentric | 3,043 | 437 | 67 | **3,547** |
| **Total** | **4,164** | **1,879** | **635** | **6,678** |

Notable catalog-level differences:
- **SXS** is dominated by precessing quasi-circular (f): 73% of all SXS sims.
  Only 6% are eccentric.
- **RIT** has a large eccentric population: 45% eccentric overall, with
  non-spinning eccentric (a) as its single largest category (499 sims).
- **MAYA** is the most eccentricity-focused: 78% eccentric.  Its largest
  single category is precessing-spin eccentric (e) at 303 sims (48% of MAYA).

### NRSur7dq4 calibration sub-classification (SXS only)

The surrogate was trained on 1,731 SXS simulations, all quasi-circular
(categories b, d, f).  The JSON file `catalog_organization/sxs_classification.json`
stores an additional `nrsur7dq4_calibration` boolean per simulation and
provides `nrsur_calibration_count` per category:

| SXS category | Total | NRSur7dq4 training |
|---|---:|---:|
| (b) non-spinning non-eccentric | 177 | 60 |
| (d) aligned-spin non-eccentric | 687 | 282 |
| (f) precessing-spin non-eccentric | 3,043 | 1,389 |
| **Total** | **3,907** | **1,731** |

Eccentric categories (a, c, e) have zero calibration sims by construction —
NRSur7dq4 does not model eccentricity.

The `NRCatalogClassifier.get_simulations()` method accepts an
`only_nrsur_calibration=True` flag to restrict to training-set sims (SXS only).

### Accessing classifications in code

```python
from nrcatalogtools.classification import NRCatalogClassifier

clf = NRCatalogClassifier(spin_threshold=0.001, ecc_threshold=0.005)

# All aligned-spin quasi-circular SXS sims
sims_d_sxs = clf.get_simulations('SXS', 'd')
# or equivalently:
sims_d_sxs = clf.get_simulations('SXS', 'aligned-spin non-eccentric')

# Only SXS calibration sims in category (d)
sims_d_cal = clf.get_simulations('SXS', 'd', only_nrsur_calibration=True)

# All RIT non-spinning (both eccentric and quasi-circular)
sims_a_rit = clf.get_simulations('RIT', 'a')
sims_b_rit = clf.get_simulations('RIT', 'b')
```

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

- **Handling Precessing Spin Epochs**: The spins in `params` are extracted at the NR relaxation time (corresponding to `f_lower`). Because the spins precess over time, we must carefully align the evaluation epochs based on the relative lengths of the waveforms:
  1. **NR is shorter than Surrogate (`f_lower_NR > f_min_surrogate`)**:
     We can generate the full surrogate waveform by letting `f_low` default to its minimum (`f_low=None`), but we *must* explicitly set `f_ref = f_lower_NR`. The surrogate will correctly backward-evolve the spins from the NR epoch to its starting frequency.
  2. **NR is longer than Surrogate (`f_lower_NR < f_min_surrogate`)**:
     The surrogate cannot extrapolate backward before its minimum frequency (`~0.0171` dimensionless). We cannot pass the `params` spins to the surrogate because it will crash (`f_ref < f_low` error). Instead, we must set `f_low = None` (surrogate default min) and `f_ref = f_low_surrogate`, and we must **extract the instantaneous spins from the NR dynamics (e.g., `Horizons.h5` or `Spin*.asc`) at the time when the NR waveform reaches `f_min_surrogate`**. We cannot use the static `metadata.json` spins.
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

# --- Generate surrogate (Handle Epoch Mismatch) ---
import gwsurrogate as gws
sur = gws.LoadSurrogate('NRSur7dq4')
q   = params['mass1'] / params['mass2']       # m1 >= m2

# Determine surrogate's minimum valid frequency (approx 0.0171 for q=2)
# Here we represent the conditional logic for epochs
if params['f_lower'] >= 0.018:  # NR is shorter (approximate bound)
    chiA = [params['spin1x'], params['spin1y'], params['spin1z']]
    chiB = [params['spin2x'], params['spin2y'], params['spin2z']]
    f_low_sur = None            # Use surrogate's full length
    f_ref_sur = params['f_lower']
else:                           # NR is longer
    # Extract spins at f_min_surrogate from NR dynamics!
    chiA, chiB = extract_nr_spins_at_frequency(wfm, target_f=0.018)
    f_low_sur = None
    f_ref_sur = 0.018

t_sur, h_sur, dyn_sur = sur(q, chiA, chiB, ellMax=4, 
                            dt=dt_dimless, f_low=f_low_sur, f_ref=f_ref_sur)

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

## Step 2 — Batch comparison: non-spinning and aligned-spin catalogs (categories a–d)

### Goal

Apply the Step 1 infrastructure in batch mode to **all** simulations in the
SXS, RIT, and MAYA catalogs belonging to the non-spinning (a, b) and
aligned-spin (c, d) categories (both quasi-circular and eccentric variants),
filtered to lie within the NRSur7dq4 prior volume ($q \in [1,4]$,
$|\chi_{1,2}| \leq 0.8$).

### Rationale — why categories a–d, and why no SO(3) optimization yet

For non-spinning and aligned-spin systems the spin vectors are always parallel
(or anti-parallel, or zero) to the orbital angular momentum $\hat{L}$; they
carry no in-plane ($\chi_\perp = 0$) components and do not precess.
Consequently:

1. **z-axis alignment is guaranteed.** Both the NR simulation and NRSur7dq4
   define their source-frame z-axis as $\hat{L}$ evaluated at the reference
   epoch. For aligned/non-spinning binaries this axis is fixed throughout the
   evolution, so the two frames agree up to a rotation by some angle $\phi_0$
   about $z$.
2. **Rotation about $z$ ≡ overall phase shift.** A rotation $R_z(\phi_0)$
   transforms mode $h_{\ell m} \to e^{-im\phi_0} h_{\ell m}$, which is
   exactly the coalescence-phase degree of freedom already maximized over by
   `pycbc.filter.match()`. No additional SO(3) optimization step is needed.
3. **Clean, unambiguous matches.** Any residual mismatch after time+phase
   maximization reflects genuine waveform differences — numerical errors in the
   NR codes or surrogate interpolation error — free from frame-convention
   artifacts. This makes categories a–d the ideal first batch to establish
   cross-catalog accuracy baselines.
4. **Eccentric sims as a negative control.** Categories (a) and (c) contain
   eccentric orbits. NRSur7dq4 is trained on quasi-circular orbits, so matches
   for these are expected to be substantially lower. Including them in the same
   batch provides a built-in negative control: the pipeline should automatically
   return low matches for eccentric configurations, confirming it is sensitive
   to physical differences and not just returning 1 by default.

### Scope

From the classification data (`catalog_organization/`):

| Category | SXS | RIT | MAYA | Total |
|---|---:|---:|---:|---:|
| (a) non-spinning eccentric | 206 | 499 | 74 | 779 |
| (b) non-spinning non-eccentric | 177 | 54 | 34 | 265 |
| (c) aligned-spin eccentric | 21 | 231 | 117 | 369 |
| (d) aligned-spin non-eccentric | 687 | 541 | 40 | 1,268 |
| **Total (a–d)** | **1,091** | **1,325** | **265** | **2,681** |

After filtering by the NRSur7dq4 prior cuts and excluding simulations where
metadata retrieval fails, we expect roughly **400–700** processable
simulations.

### Design

#### 2.1 Simulation selection

Use `NRCatalogClassifier` to enumerate the target simulations, then filter
by the surrogate prior:

```python
from nrcatalogtools.classification import NRCatalogClassifier
from project.scripts.catalog_utils import load_catalog
from project.scripts.surrogate_utils import check_surrogate_prior

TARGET_CATEGORIES = [
    'non-spinning eccentric',       # a
    'non-spinning non-eccentric',   # b
    'aligned-spin eccentric',       # c
    'aligned-spin non-eccentric',   # d
]

clf = NRCatalogClassifier(spin_threshold=0.001, ecc_threshold=0.005)

sims_to_run = {}  # {catalog_name: [sim_name, ...]}
for catalog_name in ['SXS', 'RIT', 'MAYA']:
    clf.classify_all(catalog_name)
    candidates = []
    for cat in TARGET_CATEGORIES:
        candidates.extend(clf.get_simulations(catalog_name, cat))

    cat_obj = load_catalog(catalog_name)
    passing = []
    for sim in candidates:
        try:
            params = cat_obj.get_parameters(sim, total_mass=40.)
        except Exception:
            continue
        if check_surrogate_prior(params):
            passing.append(sim)
    sims_to_run[catalog_name] = passing
```

#### 2.2 Processing loop

For each `(catalog_name, sim_name)` pair call the existing
`compare_sim_vs_surrogate()` function from `compare_one_sim_vs_surrogate.py`
with `rotate=False` at the reference total mass $M = 40\,M_\odot$.
Accumulate results into per-catalog CSV files and a merged all-catalogs CSV.
An optional `--mass-scan` flag triggers `mass_scan.py`'s per-simulation mass
grid for a random 10% subsample of category (b) and (d) simulations.

#### 2.3 Parallelization

Use `multiprocessing.Pool` with one worker per simulation. Each worker
independently instantiates its own catalog and surrogate (gwsurrogate loads
a model-file singleton per process). Pool size defaults to
`min(os.cpu_count(), 8)`. A `tqdm`-based progress bar tracks completion.

#### 2.4 Output

- `results/batch_aligned_sxs.csv` — per-mode matches for all processed SXS
  simulations in categories a–d
- `results/batch_aligned_rit.csv` — same for RIT
- `results/batch_aligned_maya.csv` — same for MAYA
- `results/batch_aligned_all.csv` — merged cross-catalog table
- **Figure A**: 2×3 panel grid of $\log_{10}(1 - \mathcal{F}_{22})$ vs $q$
  and $\chi_{\rm eff}$ per catalog, points color-coded by sub-category (a/b/c/d)
- **Figure B**: Per-mode match CDFs for categories (b) and (d) across catalogs,
  to compare quasi-circular performance between SXS, RIT, and MAYA

### New files

| File | Purpose |
|---|---|
| `scripts/batch_aligned_catalogs.py` | Main batch loop for categories a–d |

---

## Step 3 — SO(3) frame-rotation optimized match

### Goal
Extend the per-simulation comparison to maximize the match not just over
time/phase but over the full SO(3) source-frame rotation, implementing the
formalism in `docs/goal.md §Source frame ambiguity`. This step is necessary
for precessing-spin systems (categories e and f) where the in-plane spin
components introduce a non-trivial rotation between the NR simulation frame
and the surrogate frame.

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
| `scripts/compare_one_sim_rotation_optimized.py` | Step 3 script |
| `nrcatalogtools/waveform/surrogate_wrapper.py` | Wrap surrogate dict → WaveformModes |

---

## Step 4 — Batch processing over full catalogs (all categories)

### Goal
Run Steps 1/2 and 3 over **all** simulations in each catalog that fall within
the `NRSur7dq4` prior volume, including the precessing-spin categories (e, f):

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
        if result.category in ['e', 'f']:
            result_rot = compute_rotation_optimized_match(catalog, sim, params)
        append_to_hdf5(results_file, sim, catalog, result, result_rot)
```

Use `multiprocessing.Pool` for embarrassingly parallel execution.

### New files

| File | Purpose |
|---|---|
| `scripts/batch_compare_catalog.py` | Batch loop with multiprocessing |
| `scripts/results_io.py` | HDF5 I/O for match results |

---

## Step 5 — BMS supertranslation correction and analysis

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

## Step 6 — Visualization and scientific analysis

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
- SXS category (b) non-spinning non-eccentric, NRSur7dq4 calibration sims
  (60 simulations): expect $\mathcal{F}_{22} > 0.99$ for the large majority.
- RIT/MAYA category (b): expect $\mathcal{F}_{22} > 0.99$ if codes are
  accurate, with any systematic offset revealing catalog-level numerical errors.
- Category (a) and (c) eccentric sims: expect noticeably lower matches
  ($< 0.95$), confirming the pipeline correctly penalizes eccentricity mismatch.
- Consistency check: Step 2 results for the four pilot SXS simulations should
  reproduce the Step 1 results within floating-point precision.

### Step 3 verification
- For a q=1 non-spinning system, the optimal rotation R* should be very close
  to identity (or a known symmetry rotation π about z).
- Compare `match_sphere_averaged` result with Step 1 result; rotation
  optimization should improve or maintain the match.

### Step 4 verification
- Spot-check a few simulations from Steps 1 and 2 to confirm batch results agree.
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

→ **Step 1 is complete.** Pilot results exist for four SXS simulations
  (SXS:BBH:0001, 0005, 0169, 0162) at $M = 40\,M_\odot$.

→ **Implement Step 2**: `scripts/batch_aligned_catalogs.py`
  — batch comparison of all categories (a–d) from SXS, RIT, and MAYA
  against NRSur7dq4, reusing the Step 1 infrastructure without SO(3)
  optimization.
