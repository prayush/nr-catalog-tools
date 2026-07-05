---
title: nrcatalogtools.registry
parent: API Reference
nav_order: 14
---

<!-- GENERATED FILE — DO NOT EDIT. Regenerate with `python bin/generate_api_docs.py`. -->
{% raw %}

# `nrcatalogtools.registry`

Catalog plugin registry.

Provides a lightweight decorator + lookup mechanism so that new catalogs
can be registered without editing ``__init__.py`` or any core module.

Built-in registration
---------------------
``RITCatalog``, ``SXSCatalog``, and ``MayaCatalog`` are all registered
automatically when the ``nrcatalogtools`` package is imported.

Third-party registration
------------------------
A downstream package (or a user in an interactive session) can register
an additional catalog at runtime::

    from nrcatalogtools.registry import register_catalog
    from nrcatalogtools.catalog import CatalogBase

    @register_catalog("LVCNR")
    class LVCNRCatalog(CatalogBase):
        CATALOG_TYPE = "LVCNR"
        ...

Lookup
------
::

    from nrcatalogtools.registry import get_catalog
    cls = get_catalog("RIT")   # → RITCatalog
    obj = cls.load()


## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### `register_catalog`

```python
register_catalog(tag: str) -> Callable[[Type], Type]
```

Class decorator that registers a catalog under *tag*.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `tag` | `str` | Short uppercase identifier (e.g. ``"RIT"``). Must be unique within the registry. |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `callable` | A class decorator; the class itself is returned unchanged so the decorator can be stacked with other decorators. |

#### Raises

| Exception | Condition |
|---|---|
| `ValueError` | If *tag* is already registered, to prevent silent overwrites. |

#### Examples

```python
>>> @register_catalog("LVCNR")
... class LVCNRCatalog(CatalogBase):
...     CATALOG_TYPE = "LVCNR"
```

---

### `get_catalog`

```python
get_catalog(tag: str) -> Type
```

Return the catalog class registered under *tag*.

#### Parameters

| Name | Type | Description |
|---|---|---|
| `tag` | `str` | Short uppercase identifier (e.g. ``"RIT"``). |

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `type` | The registered catalog class. |

#### Raises

| Exception | Condition |
|---|---|
| `KeyError` | If *tag* is not in the registry. |

#### Examples

```python
>>> cls = get_catalog("SXS")
>>> catalog = cls.load()
```

---

### `list_catalogs`

```python
list_catalogs() -> set
```

Return the set of all registered catalog tags.

#### Returns

| Name | Type | Description |
|---|---|---|
|  | `set[str]` | Copy of the current tag set; modifying it has no effect on the registry. |

#### Examples

```python
>>> list_catalogs()
{'MAYA', 'RIT', 'SXS'}
```

---

{% endraw %}
