---
title: nrcatalogtools.catalog
parent: API Reference
nav_order: 2
---

<!-- GENERATED FILE — DO NOT EDIT. Regenerate with `python bin/generate_api_docs.py`. -->
{% raw %}

# `nrcatalogtools.catalog`

Abstract base classes and shared implementation for NR waveform catalogs.

### Classes

| Name | Description |
|---|---|
| `CatalogABC` | Pure abstract interface that every catalog must implement. Declares the filename/filepath/URL/download contract for waveform, psi4, and metadata data products. |
| `CatalogBase` | Concrete base class that combines ``CatalogABC`` with a plain-dict simulation registry and provides the shared ``get()``, ``get_metadata()``, ``get_parameters()``, and ``set_attribute_in_waveform_data_file()`` implementations used by all three catalog back-ends (RIT, SXS, MAYA). Subclasses must set ``CATALOG_TYPE`` (e.g. ``"RIT"``) and implement all abstract methods declared in ``CatalogABC``. |


## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## *class* `CatalogABC`

Bases: `ABC`

Pure abstract interface that every catalog back-end must implement.

Subclasses supply the catalog-specific logic for resolving file names,
local paths, remote URLs, and triggering downloads for the three data
products: *waveform strain* (HDF5), *psi4* (HDF5 or tar.gz), and
*per-simulation metadata* (text or JSON).

All methods take a ``sim_name`` string as their first positional
argument.  The naming convention for ``sim_name`` is catalog-specific
(e.g. ``"RIT:BBH:0001-n100-id3"``, ``"SXS:BBH:0001"``,
``"GT0001"``).


### `waveform_filename_from_simname`

```python
waveform_filename_from_simname(sim_name: str) -> str
```

Return the bare filename (no directory) for the waveform HDF5 file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Filename, e.g. ``"ExtrapStrain_RIT-BBH-0001-n100.h5"``. |

---

### `waveform_filepath_from_simname`

```python
waveform_filepath_from_simname(sim_name: str) -> str
```

Return the absolute local path for the waveform HDF5 file.

The file may not yet exist on disk if it has not been downloaded.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Absolute path under the catalog cache directory. |

---

### `waveform_url_from_simname`

```python
waveform_url_from_simname(sim_name: str) -> str
```

Return the remote URL for the waveform HDF5 file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Full HTTP(S) URL. |

---

### `download_waveform_data`

```python
download_waveform_data(sim_name: str)
```

Download the waveform data file for *sim_name* into the local cache.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

---

### `psi4_filename_from_simname`

```python
psi4_filename_from_simname(sim_name: str) -> str
```

Return the bare filename for the psi4 data file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Filename, e.g. ``"ExtrapPsi4_RIT-BBH-0001-n100-id3.tar.gz"``. |

#### Raises

| Exception | Condition |
|---|---|
| `NotImplementedError` | If the catalog does not distribute psi4 data. |

---

### `psi4_filepath_from_simname`

```python
psi4_filepath_from_simname(sim_name: str) -> str
```

Return the absolute local path for the psi4 data file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Absolute path under the catalog cache directory. |

#### Raises

| Exception | Condition |
|---|---|
| `NotImplementedError` | If the catalog does not distribute psi4 data. |

---

### `psi4_url_from_simname`

```python
psi4_url_from_simname(sim_name: str) -> str
```

Return the remote URL for the psi4 data file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Full HTTP(S) URL. |

#### Raises

| Exception | Condition |
|---|---|
| `NotImplementedError` | If the catalog does not distribute psi4 data. |

---

### `download_psi4_data`

```python
download_psi4_data(sim_name: str)
```

Download the psi4 data file for *sim_name* into the local cache.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Raises

| Exception | Condition |
|---|---|
| `NotImplementedError` | If the catalog does not distribute psi4 data. |

---

### `metadata_filename_from_simname`

```python
metadata_filename_from_simname(sim_name: str) -> str
```

Return the bare filename for the per-simulation metadata file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Filename, e.g. ``"RIT:BBH:0001-n100-id3_Metadata.txt"``. |

---

### `metadata_filepath_from_simname`

```python
metadata_filepath_from_simname(sim_name: str) -> str
```

Return the absolute local path for the per-simulation metadata file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Absolute path under the catalog cache directory. |

---

### `metadata_url_from_simname`

```python
metadata_url_from_simname(sim_name: str) -> str
```

Return the remote URL for the per-simulation metadata file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Simulation name tag. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `str` | `str` | Full HTTP(S) URL. |

---

## *class* `CatalogBase`

Bases: `CatalogABC`

Shared implementation base for all NR catalog back-ends.

Owns ``self._simulations: dict[str, dict]`` (simulation name →
metadata dict) directly, with no dependency on ``sxs.Catalog``.
Provides the default ``get()``, ``get_metadata()``,
``get_parameters()``, ``to_sxs()``, and
``set_attribute_in_waveform_data_file()`` implementations shared
by ``RITCatalog``, ``SXSCatalog``, and ``MayaCatalog``.

Subclasses must:
- Set ``CATALOG_TYPE`` to ``"RIT"``, ``"SXS"``, or ``"MAYA"``.
- Implement all ``@abstractmethod`` methods from ``CatalogABC``
  (filename, filepath, URL, and download helpers for waveform,
  psi4, and metadata products).
- Override ``get()`` if catalog-specific download logic is needed
  (``SXSCatalog`` overrides because SXS data is managed by the
  ``sxs`` package, not by local HDF5 files).


### *property* `simulations`

Mapping of simulation name → metadata dict for all simulations.

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `dict` | dict[str, dict]: The full simulation registry. |

---

### *property* `simulations_list`

List of all simulation name tags in the catalog.

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `list` | list[str]: Simulation names in insertion order. |

---

### `to_sxs`

```python
to_sxs() -> 'sxs.Simulations'
```

Return an ``sxs.Simulations`` view of this catalog's metadata.

For ``SXSCatalog``, the live ``sxs.Simulations`` object (fully
populated with ``.dataframe``, ``.tag``, etc.) is returned via an
override.  For RIT and MAYA catalogs this constructs an
``sxs.Simulations`` object from the local metadata dict; sxs-specific
columns in ``.dataframe`` will be NaN because RIT/MAYA keys do not
match the SXS schema.

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `'sxs.Simulations'` | sxs.Simulations: An sxs-native catalog object. |

---

### `save`

```python
save(file: str)
```

Save this catalog's metadata to a JSON file.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `file` | `str` | Path to the output JSON file. |

---

### `get`

```python
get(sim_name: str, quantity: str = 'waveform', **kwargs) -> waveform.WaveformModes
```

Retrieve specific quantities for one simulation

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Name of simulation in catalog |
| `quantity` | `str` | Name of quantity to fetch. Options: {waveform, psi4} |

#### Raises

| Exception | Condition |
|---|---|
| `IOError` | If `sim_name` not found in the catalog |
| `IOError` | If `quantity` is not one of the options above |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `waveform.WaveformModes` | nrcatalogtools.waveform.WaveformModes: Waveform modes |

---

### `get_metadata`

```python
get_metadata(sim_name: str) -> dict
```

Get Metadata for this simulation

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Name of simulation in catalog |

#### Raises

| Exception | Condition |
|---|---|
| `IOError` | If `sim_name` is not found in the catalog |

#### Returns

| Name | Type | Description |
|---|---|---|
| `dict` | `dict` | Metadata dictionary. Always contains a ``"catalog_type"`` |
|  | `dict` | key (value: ``"RIT"``, ``"SXS"``, or ``"MAYA"``) so that |
|  | `dict` | downstream code can dispatch without fragile sentinel-key checks. |

---

### `set_attribute_in_waveform_data_file`

```python
set_attribute_in_waveform_data_file(sim_name: str, attr_name: str, attr_value: Any)
```

Set attributes in the HDF5 file holding waveform data for a given
simulation

#### Parameters

| Name | Type | Description |
|---|---|---|
| `sim_name` | `str` | Name/Tag of the simulation |
| `attr_name` | `str` | Name of the attribute to set |
| `attr_value` | `any / serializable` | Value of the attribute |

---

### `get_parameters`

```python
get_parameters(sim_name: str, total_mass: float = 1.0) -> dict
```

Return the initial physical parameters for the simulation. Only for
quasicircular simulations are supported, orbital eccentricity is ignored

#### Parameters

| Name | Type | Description |
|---|---|---|
| `total_mass` | `float` | Total Mass of Binary (solar masses). Defaults to 1.0. |

#### Returns

| Name | Type | Description |
|---|---|---|
| `dict` | `dict` | Initial binary parameters with names compatible with PyCBC. |

---

{% endraw %}
