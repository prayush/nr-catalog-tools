# API Reference

Auto-generated from source docstrings via [mkdocstrings](https://mkdocstrings.github.io/).

| Module | Description |
|--------|-------------|
| [catalog](catalog.md) | Abstract base classes `CatalogABC` and `CatalogBase` |
| [rit](rit.md) | `RITCatalog` and `RITCatalogHelper` for the RIT catalog |
| [sxs](sxs.md) | `SXSCatalog` wrapping the `sxs` package |
| [maya](maya.md) | `MayaCatalog` for the Georgia Tech MAYA catalog |
| [waveform](waveform.md) | `WaveformModes` and matching/rotation helpers |
| [surrogate](surrogate.md) | NRSur7dq4 loading, evaluation, and prior check |
| [comparisons](comparisons.md) | End-to-end NR vs surrogate comparison pipeline |
| [metadata](metadata.md) | Cross-catalog key mappings and `get_source_parameters_from_metadata` |
| [registry](registry.md) | Catalog plugin registry |
| [lvc](lvc.md) | Frame-rotation helpers and LVCNR format utilities |
| [utils](utils.md) | Cache paths, download helpers, unit conversions |

## Quick navigation

All public classes and functions are also importable directly from the top-level package:

```python
import nrcatalogtools as nrcat

# Catalogs
nrcat.RITCatalog
nrcat.SXSCatalog
nrcat.MayaCatalog

# Catalog helpers
nrcat.load_catalog          # load by name tag: "SXS", "RIT", or "MAYA"
nrcat.filter_by_surrogate_prior  # filter sims by NRSur7dq4 validity volume

# Waveform
nrcat.WaveformModes
nrcat.apply_wigner_rotation_to_mode_dict

# Registry
nrcat.register_catalog
nrcat.get_catalog
nrcat.list_catalogs

# Metadata key maps
nrcat.RIT_KEYS
nrcat.SXS_KEYS
nrcat.MAYA_KEYS
nrcat.CANONICAL_TO_CATALOG
nrcat.CANONICAL_TO_PYCBC
nrcat.PYCBC_KEYS
```

### Surrogate and matching (sub-module imports)

```python
from nrcatalogtools.surrogate import (
    load_nrsur7dq4,
    generate_surrogate_modes,
    check_surrogate_prior,
    SURROGATE_MODES,
    NR_MODES,
)

from nrcatalogtools.waveform.matching import (
    load_psd,
    compute_mode_match,
    compute_phase_diff_per_cycle,
    mode_f_lower,
)

from nrcatalogtools.comparisons import compare_sim_vs_surrogate, DELTA_T
```
