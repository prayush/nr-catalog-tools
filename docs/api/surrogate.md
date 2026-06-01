# `nrcatalogtools.surrogate`

NRSur7dq4 loading, evaluation, and prior check utilities.

> **`gwsurrogate` is an optional dependency.** Importing this module succeeds even when
> `gwsurrogate` is not installed; the `ImportError` is raised only on the first call to
> `load_nrsur7dq4()` or `generate_surrogate_modes()`.

> **Conceptual background**: See [Package Internals](../package.md) § 11 for notes on the
> `f_ref` clipping strategy and spin-epoch conventions.

---

## Functions

::: nrcatalogtools.surrogate
    options:
      members:
        - load_nrsur7dq4
        - generate_surrogate_modes
        - check_surrogate_prior

---

## Constants

| Name | Description |
|------|-------------|
| `SURROGATE_MODES` | Mode list output by NRSur7dq4: `(2,1)`, `(2,2)`, `(3,2)`, `(3,3)`, `(4,3)`, `(4,4)` |
| `NR_MODES` | Superset of `SURROGATE_MODES` — all modes present in the NR catalogs, including `(5,5)` |
| `_SURROGATE_F_MIN_CYCLES_PER_M` | Minimum $M f_{GW}$ (cycles/M) for NRSur7dq4 training (`0.0165`) |
