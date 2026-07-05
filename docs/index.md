---
title: Home
nav_order: 1
permalink: /
---

# nrcats
{: .fs-9 }

A stable, unified Python interface to public numerical-relativity (NR) binary black-hole
waveform catalogs — for LVK analyses, waveform modeling, and cross-catalog comparison.
{: .fs-6 .fw-300 }

[Install](installation.md){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Quick start](#quick-start){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/gwnrtools/nrcats){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## What it does

`nrcats` serves three overlapping communities with one API:

- **LIGO-Virgo-KAGRA analyses** — reliable, PyCBC-compatible waveform time series and
  source-parameter dicts for parameter estimation, injection studies, and template bank
  construction, with consistent physical units (masses in M☉, distances in Mpc, epoch at
  the (2,2) peak) across all backends
- **Waveform modeling** — consistent loading, physical scaling, and frame-alignment tools
  (Wigner D-matrix rotations, `f_lower` extraction, time/phase alignment) for calibrating
  and validating EOB, phenomenological, and surrogate models against any NR catalog
- **Cross-catalog studies** — noise-weighted mismatch computation maximized over time and
  phase shifts, $$SO(3)$$ rotations, and BMS supertranslations, for quantifying NR waveform
  accuracy across codes

All catalog backends expose an **identical interface** (defined by
[`CatalogBase`](api/catalog.md)), so analysis code written against one catalog works
against all others without modification.

## Supported catalogs

| Catalog | Code | Example simulation name |
|---------|------|------------------------|
| [SXS](https://data.black-holes.org/waveforms/catalog.html) | SpEC | `SXS:BBH:0001` |
| [RIT](https://ccrg.rit.edu/content/data/rit-waveform-catalog) | LazEv | `RIT:BBH:0001-n100-id3` |
| [MAYA / GT](https://einstein.gatech.edu/catalog/) | MayaKranc | `GT0001` |

New catalogs can be added without touching core code via the
[plugin registry](api/registry.md).

## Quick start

```python
import nrcats as nrcat

# Load catalogs (explicit class methods)
ritcat  = nrcat.RITCatalog.load()
sxscat  = nrcat.SXSCatalog.load(download=False)
mayacat = nrcat.MayaCatalog.load()

# ...or use the unified helper
ritcat = nrcat.load_catalog("RIT")

# Browse simulations
print(ritcat.simulations_dataframe.index)
# Index(['RIT:BBH:0001-n100-id3', 'RIT:BBH:0002-n100-id0', ...], length=1879)

# Load a waveform
wfm = ritcat.get("RIT:BBH:0003-n100-id0")
print(wfm.LM)     # available (ell, m) mode pairs

# Extract the (2,2) mode in physical units
mode22 = wfm.get_mode(2, 2, total_mass=60.0, distance=100.0,
                      delta_t_seconds=1./4096)

# Polarizations
pols = wfm.get_td_waveform(total_mass=40., distance=100.,
                           inclination=0.2, coa_phase=0.3)
hp, hc = pols.real(), -1 * pols.imag()

# PyCBC-compatible source parameters
params = ritcat.get_parameters("RIT:BBH:0001-n100-id3", total_mass=60.0)
# {'mass1': 30.0, 'mass2': 30.0, 'spin1x': 0.0, ..., 'f_lower': 23.4}
```

For a step-by-step walkthrough, start with the
[Loading a waveform](tutorials/load-waveform.md) tutorial.

## Documentation map

| Section | What you'll find |
|---------|------------------|
| [Installation](installation.md) | pip/conda install, dependencies, cache configuration |
| [Tutorials](tutorials/index.md) | Worked examples: loading waveforms, cross-catalog mismatches, surrogate comparison |
| [Catalog Reference](catalogs.md) | Per-catalog loading, metadata keys, file formats, cache layout |
| [Catalog Organization](catalog_organization.md) | Spin/eccentricity classification of the three catalogs |
| [WaveformModes Guide](waveform.md) | Conceptual guide to the central waveform object and its unit conventions |
| [Architecture](architecture.md) | Class hierarchy, data flows, key design decisions |
| [API Reference](api/index.md) | Per-module reference generated from the source docstrings |
| [Package Internals](package.md) | Detailed module descriptions, unit conventions, usage patterns, gotchas |
| [Scientific Goal](goal.md) | Motivation, source-frame ambiguity, BMS supertranslations, mismatch formalism |
| [Contributing](contributing.md) | Development setup, tests, docs builds, adding a new catalog |
| [Changelog](changelog.md) | Release history |

## Module structure

```
nrcats/
├── __init__.py        # Public API: RITCatalog, SXSCatalog, MayaCatalog,
│                      #   WaveformModes, load_catalog, registry, key maps
├── catalog.py         # Abstract base CatalogABC + shared CatalogBase
├── rit.py             # RITCatalog + RITCatalogHelper
├── sxs.py             # SXSCatalog
├── maya.py            # MayaCatalog
├── registry.py        # @register_catalog plugin mechanism
├── metadata.py        # Cross-catalog key maps + get_source_parameters_from_metadata
├── classification.py  # Spin/eccentricity classification of catalog simulations
├── surrogate.py       # NRSur7dq4 loading, evaluation, prior check
├── comparisons.py     # End-to-end NR vs surrogate comparison pipeline
├── waveform/          # WaveformModes sub-package
│   ├── modes.py       #   WaveformModes class
│   ├── loaders.py     #   load_from_h5, load_from_targz
│   ├── matching.py    #   mode matching, Wigner rotation, PSD helpers
│   └── units.py       #   waveform-level constants
├── lvc.py             # Frame-rotation helpers and LVCNR format utilities
├── utils.py           # Cache paths, download helpers, unit conversions
└── schemas/           # YAML key-mapping tables (rit/sxs/maya)
```
