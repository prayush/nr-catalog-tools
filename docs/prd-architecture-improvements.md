# Product Requirements Document: Architecture Improvements

**Project:** nr-catalog-tools  
**Author:** Prayush Kumar  
**Date:** 2026-05-25  
**Status:** Draft — Pending Approval  

---

## 1. Background and Motivation

`nr-catalog-tools` provides a unified Python interface to three numerical-relativity (NR) waveform catalogs: RIT, SXS, and MAYA/GT. The library is used for cross-catalog validation, parameter extraction, and waveform injection studies.

An architectural review identified correctness risks, a 60 KB monolithic module, heavy coupling to the `sxs` package at the inheritance level, and several latent bugs. This PRD captures the full scope of proposed improvements so that work can be planned, reviewed, and executed in prioritised phases.

---

## 2. Goals

- Eliminate known correctness bugs and silent failure modes.
- Make the API explicit and self-documenting where it is currently heuristic-driven.
- Reduce coupling to the `sxs` third-party package to only what is necessary.
- Make the codebase easier to extend (new catalogs, new formats) without editing core files.
- Modernise the package build and versioning infrastructure.

## 3. Non-Goals

- Changing the scientific algorithms (mismatch computation, Wigner rotation, mode summation).
- Supporting new NR catalogs beyond RIT, SXS, and MAYA in this work (the registry architecture enables that as follow-on).
- Changing the public-facing API surface in Phase 1 or Phase 2.

---

## 4. Requirements

Requirements are grouped into three phases. Each has a priority (P0 = must-fix, P1 = should-fix, P2 = nice-to-have), an effort estimate (S/M/L/XL), and notes on backward compatibility.

---

### Phase 1 — Correctness and Reliability

These changes fix real bugs and brittle behaviours with no intentional API breaks.

---

#### REQ-1.1 Fix `scipy.stats.mode` deprecation

| Field | Value |
|---|---|
| **Priority** | P0 |
| **Effort** | S |
| **Files** | `nrcatalogtools/waveform.py` ~L154, ~L295 |
| **BC** | Non-breaking |

**Current behaviour:** `stat_mode(np.diff(amp_time), keepdims=True)[0][0]` is called to find the most common timestep when building the uniform interpolation grid. The `keepdims` parameter and the return-value shape of `scipy.stats.mode` changed between scipy 1.9 and 1.11, causing a `TypeError` or an off-by-one index error on newer scipy installs.

**Required behaviour:** Replace with `np.median(np.diff(time_array))`, or — if a mode (most common value) is genuinely needed — compute it via `np.unique(np.round(np.diff(time_array), 10), return_counts=True)` and pick the argmax. Add a unit test asserting that a non-uniform input time array produces a correct uniform grid.

---

#### REQ-1.2 Replace the `delta_t` dual-convention magic number

| Field | Value |
|---|---|
| **Priority** | P0 |
| **Effort** | M |
| **Files** | `nrcatalogtools/waveform.py` ~L562–573, ~L490 |
| **BC** | Deprecated old signature; new preferred signature |

**Current behaviour:** In `get_mode()` and `get_td_waveform()`, `delta_t > 1/128` is silently interpreted as dimensionless M units; `delta_t ≤ 1/128` as physical seconds. A caller passing `delta_t=0.001` meaning "one millisecond" gets NR-units behaviour with no warning.

**Required behaviour:**

```python
# New preferred API
get_mode(ell, em, total_mass=1.0, distance=1.0,
         delta_t_seconds=None, delta_t_Msun=None, ...)
get_td_waveform(..., delta_t_seconds=None, delta_t_Msun=None, ...)
```

- Exactly one of `delta_t_seconds` or `delta_t_Msun` must be provided; raise `ValueError` if both or neither are given.
- The old positional `delta_t` parameter is retained with a `DeprecationWarning` that instructs callers to switch to the explicit form. The old heuristic is used for the deprecated path so existing call-sites continue to work unchanged.
- A unit test verifies that `delta_t_seconds=1/4096` and `delta_t_Msun=0.5` produce differently-sampled outputs, and that passing the old `delta_t=1/4096` triggers the warning.

---

#### REQ-1.3 Inject `catalog_type` and remove sentinel-key detection

| Field | Value |
|---|---|
| **Priority** | P0 |
| **Effort** | S |
| **Files** | `nrcatalogtools/catalog.py`, `nrcatalogtools/rit.py`, `nrcatalogtools/sxs.py`, `nrcatalogtools/maya.py`, `nrcatalogtools/metadata.py` |
| **BC** | Non-breaking |

**Current behaviour:** `get_source_parameters_from_metadata()` infers which catalog a metadata dict belongs to by checking for sentinel keys (`"relaxed_mass1"` → RIT, `"GTID"` → MAYA, else SXS). This breaks if a catalog adds or renames those keys.

**Required behaviour:** Each catalog's `get_metadata()` method injects `metadata["catalog_type"] = "RIT"` (or `"SXS"`, `"MAYA"`) before returning. `get_source_parameters_from_metadata()` reads `metadata.get("catalog_type")` and raises a clear `ValueError` if the key is absent or unrecognised. The sentinel-key fallback is removed.

---

#### REQ-1.4 Fix `lru_cache` stale-result bug on `RITCatalog.load`

| Field | Value |
|---|---|
| **Priority** | P0 |
| **Effort** | S |
| **Files** | `nrcatalogtools/rit.py` L33–100 |
| **BC** | Non-breaking (behaviour fix) |

**Current behaviour:** `@functools.lru_cache()` on `RITCatalog.load` means a call with `download=True` after a prior call with `download=False` returns the cached (potentially incomplete) catalog object. The user expects a fresh download but silently gets the stale result.

**Required behaviour:** Remove `lru_cache`. Replace with a module-level `_rit_catalog_singleton: Optional[RITCatalog] = None` variable. `load()` populates it on the first call and returns it on subsequent calls unless `download=True` is passed, in which case it re-downloads and updates the singleton. Expose `RITCatalog.reload()` as the explicit "force refresh" path, consistent with the docstring that already mentions it.

---

#### REQ-1.5 Reduce download retry count and add exponential backoff

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | S |
| **Files** | `nrcatalogtools/utils.py` L73–160 |
| **BC** | Non-breaking |

**Current behaviour:** `url_exists(link, num_retries=100)` and `download_file(...)` retry up to 100 times with no delay. A transient or permanent failure causes code to hang for many minutes.

**Required behaviour:**
- Default `num_retries=5`.
- Use exponential backoff: wait `2^attempt` seconds between retries, capped at 30 s.
- Print a progress message at each retry if `verbosity > 0`.
- Raise `ConnectionError` (not silently return `False`) when all retries are exhausted in contexts where the file is required.

---

#### REQ-1.6 Warn on silent zero-padding of missing modes in `load_from_targz`

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | S |
| **Files** | `nrcatalogtools/waveform.py` L303–311 |
| **BC** | Non-breaking |

**Current behaviour:** When iterating modes `(ℓ, m)` to build the ordered mode array, modes absent from the tar archive are silently filled with zeros. Downstream analysis code that assumes "non-zero amplitude means mode is present" will silently produce wrong integrals, mismatches, or mode-by-mode comparisons.

**Required behaviour:** Keep the zero-padding (it is required by the `sxs.WaveformModes` data layout), but emit a `UserWarning` listing each zero-padded `(ℓ, m)`. Additionally, store the set of truly-present modes in `wfm._present_modes: set[tuple[int,int]]` so that callers can check `(ell, em) in wfm._present_modes` before trusting the data.

---

### Phase 2 — Code Organisation

Non-breaking internal refactors that improve maintainability and testability.

---

#### REQ-2.1 Split `waveform.py` into a sub-package

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | M |
| **Files** | `nrcatalogtools/waveform.py` → `nrcatalogtools/waveform/` |
| **BC** | Fully backward-compatible via re-export shim |

**Current behaviour:** `waveform.py` is ~1 000 lines / 60 KB and mixes file I/O, physical unit scaling, signal processing, and frame transformations.

**Required behaviour:** Create `nrcatalogtools/waveform/` sub-package with:

| New file | Responsibility |
|---|---|
| `modes.py` | `WaveformModes` class definition, `__new__`, `_load` |
| `loaders.py` | `load_from_h5`, `load_from_targz` |
| `units.py` | `get_mode`, `f_lower_at_1Msun`, `amp_to_physical`, `time_to_physical` |
| `matching.py` | `match_sphere_averaged`, `apply_wigner_rotation_to_mode_dict` |
| `__init__.py` | Re-exports everything currently in `waveform.py` |

The top-level `nrcatalogtools/waveform.py` is replaced by a one-line import shim:
```python
from nrcatalogtools.waveform import *  # noqa: F401,F403
```
so no import paths in user code change.

---

#### REQ-2.2 Move metadata key mappings to per-catalog YAML schemas

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | M |
| **Files** | `nrcatalogtools/metadata.py`, new `nrcatalogtools/schemas/` |
| **BC** | Non-breaking; the existing Python dicts remain as the loaded representation |

**Current behaviour:** `metadata.py` contains three manually-maintained dicts of 89, 163, and 195 entries. `CANONICAL_TO_CATALOG` is assembled from them. There is no schema validation.

**Required behaviour:**
- Create `nrcatalogtools/schemas/rit_keys.yaml`, `sxs_keys.yaml`, `maya_keys.yaml`.
- Each YAML file maps canonical names to catalog-specific keys (and notes for derived quantities).
- `metadata.py` loads these at import time: `RIT_KEYS = _load_schema("rit_keys.yaml")`. The loaded dicts have the same names and types as today; no external API changes.
- Include the YAML files in `MANIFEST.in` and `pyproject.toml` package data.

---

#### REQ-2.3 Make `RITCatalog` helper delegation explicit

| Field | Value |
|---|---|
| **Priority** | P2 |
| **Effort** | S |
| **Files** | `nrcatalogtools/rit.py` L28–30 |
| **BC** | Non-breaking |

**Current behaviour:**
```python
self.refresh_metadata_df_on_disk = self._helper.refresh_metadata_df_on_disk
self.download_data_for_catalog   = self._helper.download_data_for_catalog
self.write_metadata_df_to_disk   = self._helper.write_metadata_df_to_disk
```
This leaks `RITCatalogHelper`'s internals as public attributes and bypasses future interface boundaries.

**Required behaviour:** Replace with explicit methods:
```python
def refresh_metadata_df_on_disk(self, **kwargs):
    return self._helper.refresh_metadata_df_on_disk(**kwargs)
```
This makes the interface inspectable, docstringable, and mockable in tests.

---

#### REQ-2.4 Migrate build system to `pyproject.toml`

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | M |
| **Files** | `setup.py`, `requirements.txt` → `pyproject.toml` |
| **BC** | No user-facing change |

**Current behaviour:** `setup.py` uses a custom git-hash version scheme that does not work for installed wheels (the `.git/` directory is absent after `pip install`). `requirements.txt` is the dependency source, which is not the modern standard.

**Required behaviour:**
- Create `pyproject.toml` with `[build-system]` using `setuptools >= 64`.
- Use `setuptools_scm` for automatic version derivation from git tags. Remove the hand-written version logic in `setup.py`.
- Declare all runtime dependencies in `[project.dependencies]`.
- Declare optional dependencies for heavy extras (`lalsuite`, `mayawaves`) under `[project.optional-dependencies]` so lightweight installs are possible.
- Delete `setup.py` and `requirements.txt` after migration; keep `environment.yml` for conda users.

---

### Phase 3 — Architecture (Major Version)

Breaking changes that require a version bump (e.g. `v0.1.0 → v1.0.0`). These should be gated on Phase 1 and Phase 2 completing first.

---

#### REQ-3.1 Introduce `NRCatalogBase` independent of `sxs.Catalog`; migrate SXS backend to `sxs.Simulations`

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | L |
| **Files** | `nrcatalogtools/catalog.py`, all catalog modules |
| **BC** | Breaking — removes `sxs.Catalog` from MRO; `sxs >= 2024.0.0` required |

**Current behaviour:** `CatalogBase` inherits from `sxs.Catalog` (deprecated as of `sxs` 2024.0.0 with the warning *"THIS CAN ONLY PROVIDE OLD DATA!"*). The SXS backend loads via `sxs.load("catalog")`, which reads the old `catalog.zip` format. RIT and MAYA catalogs must also carry `_dict["records"]`, `_dict["modified"]`, and other SXS-internal keys that are meaningless for them. A breaking change in `sxs` can silently corrupt RIT/MAYA behaviour.

**Required behaviour — composition over inheritance:**

Replace `class CatalogBase(CatalogABC, sxs_Catalog)` with `class CatalogBase(CatalogABC)`. The class owns `self._simulations: dict[str, dict]` directly. The `sxs.Catalog` constructor call is removed; `RITCatalog` and `MayaCatalog` no longer carry any SXS-internal keys.

**Required behaviour — SXS backend migration:**

`SXSCatalog.load()` switches from the deprecated `sxs.load("catalog")` to `sxs.load("simulations")`, which returns an `sxs.Simulations` object (an `OrderedDict` subclass backed by the new `sxscatalog` package). The live `sxs.Simulations` object is stored on the singleton as `self._sxs_simulations` for use by `to_sxs()` and forwarding properties.

```python
# Before
sxs_catalog = sxs.load("catalog", download=download, **kwargs)
_sxs_catalog_singleton = cls(catalog=sxs_catalog._dict, verbosity=verbosity)

# After
sxs_sims = sxs.load("simulations", download=download, **kwargs)
simulations_dict = {k: dict(v) for k, v in sxs_sims.items()}
_sxs_catalog_singleton = cls(simulations_dict=simulations_dict, verbosity=verbosity)
_sxs_catalog_singleton._sxs_simulations = sxs_sims
```

The `extrapolation_order` fallback in `SXSCatalog.get()` must also be updated: the new `sxs.Simulation` API uses `extrapolation="N{n}"` (string), not `extrapolation_order=n` (int):

```python
# Before (old sxs path)
raw_obj = sxs.load(f"{sim_name}/rhOverM", extrapolation_order=extrapolation_order, ...)

# After (new Simulation factory)
sim_obj = sxs.load(sim_name, extrapolation=f"N{extrapolation_order}", download=download)
raw_obj = sim_obj.strain
```

**Required behaviour — forwarding properties on `CatalogBase`:**

`sxs.Catalog` exposed several properties that callers may rely on. The table below classifies each one after the migration to `sxs.Simulations`:

| `sxs.Catalog` property | Action | Implementation |
|---|---|---|
| `simulations` | **Forward** | `return self._simulations` (plain dict on base; `self._sxs_simulations` on SXSCatalog) |
| `simulations_list` | **Forward** | `return list(self._simulations)` (already exists; change backing) |
| `simulations_dataframe` / `table` | **Forward (SXS only)** | `return self._sxs_simulations.dataframe` — returns a `SimulationsDataFrame` with sub-properties `.BBH`, `.noneccentric`, `.deprecated`, `.undeprecated`, etc. RIT/MAYA use their own DataFrame construction as today. Note: column schema differs from the old `sxs.Catalog.simulations_dataframe`. |
| `tag` | **Forward (SXS only)** | `return self._sxs_simulations.tag` — git tag of the catalog snapshot, e.g. `"v3.0.28"` |
| `published_at` | **Forward (SXS only)** | `return self._sxs_simulations.published_at` — ISO timestamp from GitHub Releases |
| `modified` | **Approximate** | Map to `self.published_at`; document the rename in `CHANGELOG.md` |
| `reload` | **Forward** | Call `sxs.Simulations.reload()` then rebuild the singleton |
| `files` | **Drop** | Not available in `sxs.Simulations`; the new API fetches file info per-simulation on demand via `sxs.load(sim_id).files`. No global flat map exists. Document in `CHANGELOG.md`. |
| `select(path_pattern)` | **Drop** | Depended on `self.files`; no replacement |
| `select_files(path_pattern)` | **Drop** | Same |
| `records` | **Drop** | Zenodo/CaltechDATA record objects not present in `simulations.json` |
| `description` | **Drop** | `catalog_file_description` key absent from `simulations.json` |
| `open_access` | **Drop** | No equivalent; use `simulations.dataframe.undeprecated` as a partial substitute |
| `split_and_write(directory)` | **Drop** | Depended on `records` and `open_access` |
| `save(file)` | **Partial** | Can write `json.dump(dict(self._simulations), ...)` but schema differs from old `catalog.json` |

**Required behaviour — `to_sxs()` method (Option B):**

Add `to_sxs()` to `CatalogBase` as a best-effort constructor that works for all three catalog types. `SXSCatalog` overrides it to return the live object.

```python
# CatalogBase (catalog.py)
def to_sxs(self) -> "sxs.Simulations":
    """Return an sxs.Simulations view of this catalog's simulation metadata.

    For SXSCatalog, returns the live sxs.Simulations object loaded from
    sxs.load("simulations"), fully populated with .dataframe, .tag, etc.

    For RIT and MAYA catalogs, constructs a Simulations object from the
    local metadata dict.  The resulting object is a valid sxs.Simulations
    instance and can be passed to sxs tooling (e.g. closest_simulation()),
    but sxs-specific columns in .dataframe will be NaN because RIT/MAYA
    metadata keys do not match the SXS schema.
    """
    import sxs
    return sxs.Simulations(self._simulations)
```

```python
# SXSCatalog (sxs.py) — override
def to_sxs(self) -> "sxs.Simulations":
    """Return the live sxs.Simulations object backing this catalog."""
    return self._sxs_simulations
```

This lets RIT/MAYA callers pass their catalog into any sxs API that accepts a `Simulations` object (including `Simulation.closest_simulation()` and dataframe-based filtering), while `SXSCatalog.to_sxs()` gives back the full, authoritative object.

**Migration notes for downstream users:**

- `isinstance(catalog, sxs.Catalog)` is now `False` for all three catalog types. Use `isinstance(catalog, nrcatalogtools.CatalogBase)` instead, or `catalog.to_sxs()` to get back an `sxs`-native object.
- `catalog.files`, `catalog.select()`, `catalog.select_files()` are removed. Per-simulation file info is available via `sxs.load(sim_id).files`.
- `catalog.simulations_dataframe` for `SXSCatalog` now returns a `SimulationsDataFrame` (richer column schema; `.BBH`, `.noneccentric`, etc. sub-properties available).
- Provide a migration guide in `CHANGELOG.md`.

---

#### REQ-3.2 Standardize WaveformModes Attribute Propagation (Inheritance Retained)

| Field | Value |
|---|---|
| **Priority** | P1 |
| **Effort** | S |
| **Files** | `nrcatalogtools/waveform/modes.py` |
| **BC** | Non-breaking |

**Current behaviour:** `WaveformModes` inherits from `sxs.WaveformModes` (itself a `numpy.ndarray` subclass). During array operations like slicing (`wfm[100:500]`) or copying, standard `numpy` mechanisms do not automatically propagate custom instance-level attributes like `filepath`, `_present_modes`, and cached `peak_time_22` unless manually set, leading to potential data loss or incorrect states in downstream computations.

**Required behaviour:**
- Retain inheritance from `sxs.WaveformModes` to preserve direct NumPy array interoperability, zero-copy math performance in mismatch routines, and compatibility with downstream numerical-relativity packages (e.g., `pycbc`, `scri`).
- Implement robust attribute propagation inside `__array_finalize__` and standard copy operations (`__copy__` and `__deepcopy__`) so that `filepath`, `_present_modes`, `peak_time_22`, and `sim_metadata` are automatically and safely preserved across all slicing, view-casting, and copying boundaries.
- Retain and clean up the `__new__` constructor logic to safely initialize custom metadata arguments without leaking class-level variables or causing `np.ndarray` allocation warnings.

---

#### REQ-3.3 Add catalog plugin registry

| Field | Value |
|---|---|
| **Priority** | P2 |
| **Effort** | M |
| **Files** | New `nrcatalogtools/registry.py`, `nrcatalogtools/__init__.py` |
| **BC** | Additive — no existing API changes |

**Current behaviour:** Adding a new catalog requires editing `__init__.py`, `metadata.py`, and `CANONICAL_TO_CATALOG`.

**Required behaviour:**
- `registry.py` exposes `@register_catalog("TAG")` decorator and `get_catalog("TAG")` lookup.
- `RITCatalog`, `SXSCatalog`, `MayaCatalog` each apply the decorator.
- `__init__.py` imports `registry` and re-exports the three known catalogs; it does not need to be edited for a new catalog.
- A fourth-party catalog can be registered at runtime: `nrcatalogtools.registry.register_catalog("LVCNR")(LVCNRCatalog)`.

---

## 5. Testing Requirements

Each requirement above has an associated test obligation:

| Req | Test obligation |
|---|---|
| REQ-1.1 | Unit test: non-uniform input → correct uniform grid dt |
| REQ-1.2 | Unit test: `delta_t_seconds` and `delta_t_Msun` produce different outputs; old `delta_t` triggers `DeprecationWarning` |
| REQ-1.3 | Unit test: metadata dict without `catalog_type` raises `ValueError` in `get_source_parameters_from_metadata` |
| REQ-1.4 | Unit test: second `load(download=True)` call triggers re-download; second `load(download=False)` returns same object |
| REQ-1.5 | Unit test (mock): 5 failures raises `ConnectionError`; 1 failure + 1 success returns result |
| REQ-1.6 | Unit test: loading a file missing one mode emits `UserWarning`; `_present_modes` excludes the padded mode |
| REQ-2.1 | Existing tests pass unchanged after the refactor |
| REQ-2.4 | `pip install -e .` succeeds from a clean clone with no `setup.py` |
| REQ-3.1 | Cross-catalog tests pass; `isinstance(RITCatalog.load(), sxs.Catalog)` is now `False`; `SXSCatalog.load().to_sxs()` returns the live `sxs.Simulations` object; `RITCatalog.load().to_sxs()` returns a valid `sxs.Simulations` instance; `SXSCatalog.load().simulations_dataframe` exposes `.BBH` sub-property |
| REQ-3.2 | Slicing and copy operations automatically propagate custom attributes (filepath, _present_modes, peak_time_22); all waveform-level tests pass |

---

## 6. Dependency Constraints

| Concern | Constraint |
|---|---|
| `sxs` | Still a required runtime dependency (SXS catalog, Wigner rotations). Pin `sxs >= 2024.0.0` (first release where `sxs.Simulations` replaces the deprecated `sxs.Catalog`). Requires Python ≥ 3.10. |
| `lalsuite` | Move to optional extras (`pip install nrcatalogtools[lal]`). Required only for `lvc.py` and `f_lower` computation via LAL. |
| `mayawaves` | Move to optional extras (`pip install nrcatalogtools[maya]`). Required only for MAYA format reading. |
| `scipy` | Pin `scipy >= 1.11` after REQ-1.1 is resolved, to allow use of the new `stats.mode` API. |
| `setuptools_scm` | New build-time dependency (REQ-2.4). |

---

## 7. Phased Delivery Schedule

| Phase | Requirements | Merge target |
|---|---|---|
| Phase 1 | REQ-1.1 through REQ-1.6 | `main` (patch releases) |
| Phase 2 | REQ-2.1 through REQ-2.4 | `main` (minor release) |
| Phase 3 | REQ-3.1 through REQ-3.3 | `main` (major release, v1.0.0) |

Phases 1 and 2 can proceed in parallel on separate branches. Phase 3 must be branched from a merged Phase 2 state.

---

## 8. Open Questions

1. **REQ-1.2 deprecation horizon** — How long should the old `delta_t` positional argument be supported before removal? (Suggested: two minor releases.)
2. **REQ-2.2 YAML vs TOML** — Should key-mapping schemas use YAML (human-editable, supports comments) or TOML (consistent with `pyproject.toml`)? YAML is recommended unless there is a strong preference.
3. **REQ-3.2 attribute propagation scope** — Which custom attributes must propagate during array finalize or slicing? (Currently: `filepath`, `_present_modes`, `sim_metadata`, and `peak_time_22` are the target set.)
4. **REQ-2.4 optional extras** — Should `lalsuite` and `mayawaves` be hard or optional dependencies? This affects whether a `pip install nrcatalogtools` on a machine without LAL works at all.
