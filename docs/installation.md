---
title: Installation
nav_order: 2
---

# Installation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Requirements

- Python ≥ 3.8 (CI tests 3.10–3.12)
- On Linux, `lalsuite` needs the GSL, FFTW3, and HDF5 system libraries. On
  Debian/Ubuntu:

```bash
sudo apt-get install libgsl-dev libfftw3-dev libhdf5-dev
```

## Install with pip

```bash
pip install nrcats
```

Or, for the latest development version:

```bash
pip install git+https://github.com/gwnrtools/nrcats.git
```

## Install with conda

The repository ships an [`environment.yml`](https://github.com/gwnrtools/nrcats/blob/master/environment.yml)
that creates a complete environment (including `lalsuite` from conda-forge):

```bash
git clone https://github.com/gwnrtools/nrcats.git
cd nrcats
conda env create -f environment.yml
conda activate nrcat
```

## Dependencies

Installed automatically by pip:

| Package | Version | Role |
|---------|---------|------|
| `sxs` | ≥ 2025.0.0 | SXS simulations access; base class `sxs.WaveformModes` |
| `pycbc` | any | `TimeSeries`, `match()`, `get_td_waveform_modes()`, `pnutils` |
| `lalsuite` | any | Physical constants (`MTSUN_SI`, `MSUN_SI`, `G_SI`, `C_SI`, `PC_SI`) |
| `h5py` | any | HDF5 reading (RIT waveform files) |
| `quaternionic` | any | Quaternion $$SO(3)$$ rotation representation |
| `spherical` | any | Wigner D-matrix computation |
| `scipy` | any | `InterpolatedUnivariateSpline` for mode resampling |
| `numpy`, `pandas`, `PyYAML` | any | Arrays, metadata DataFrames, key-map schemas |
| `mayawaves` | any | MAYA coalescence loading |
| `waveformtools` | any | Waveform post-processing helpers |

### Optional dependencies

| Package | Needed for |
|---------|-----------|
| `gwsurrogate` | NRSur7dq4 evaluation ([surrogate](api/surrogate.md), [comparisons](api/comparisons.md)). Importing those modules succeeds without it; the `ImportError` is deferred to the first surrogate call. |
| `scri` | Spin-weighted Gaunt coefficients for BMS supertranslation optimization |
| `matplotlib` | Plotting in the tutorials and the comparison pipeline's figures |

## Verify the installation

```python
import nrcatalogtools as nrcat
print(nrcat.list_catalogs())   # {'RIT', 'SXS', 'MAYA'}
```

The first `load()` of each catalog downloads its metadata; see the cache layout below
for where files land.

## Cache configuration

All downloaded catalog data lives under a single cache root, controlled by the
`NR_CATALOG_CACHE` environment variable (default: `~/.cache/`):

```bash
export NR_CATALOG_CACHE=/scratch/$USER/nrcache
```

```
$NR_CATALOG_CACHE/
├── RIT/
│   ├── metadata/
│   │   ├── metadata.csv                           # aggregated DataFrame
│   │   └── RIT:BBH:0001-n100-id3_Metadata.txt     # per-simulation files
│   └── data/
│       ├── ExtrapStrain_RIT-BBH-0001-n100.h5
│       └── ExtrapPsi4_RIT-BBH-0001-n100-id3.tar.gz
├── MAYA/
│   ├── metadata/
│   └── data/
│       └── catalog.zip                            # zipped MAYAmetadata.pkl
└── SXS/
    └── (managed by the sxs package; typically ~/.cache/sxs/)
```

{: .note }
> SXS downloads are managed by the `sxs` package itself and honor its own
> configuration (`sxs.write_config(download=..., cache=...)`); only the RIT and
> MAYA backends use the `NR_CATALOG_CACHE` layout directly.

Waveform data files are downloaded **on demand** — loading a catalog fetches only
metadata; the first `get(sim_name)` call for a simulation downloads its waveform file.

## Troubleshooting

**`ImportError` from `lal` / `lalsimulation`**
: Install the system libraries listed under [Requirements](#requirements), or prefer
  the conda route where `lalsuite` ships prebuilt.

**Catalog load fails with a network error**
: Pass `download=False` to `load()` to use only cached metadata, or retry —
  downloads use capped exponential-backoff retries (5 attempts).

**`gwsurrogate` missing when calling surrogate functions**
: `pip install gwsurrogate`. The NRSur7dq4 model data is downloaded automatically on
  the first `gwsurrogate.LoadSurrogate("NRSur7dq4")` call (wrapped by
  [`load_nrsur7dq4()`](api/surrogate.md)).
