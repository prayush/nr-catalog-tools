# `nrcatalogtools.comparisons`

End-to-end NR vs NRSur7dq4 comparison pipeline.

> **Conceptual background**: See [Package Internals](../package.md) § 12 for a description
> of the pipeline steps and output format.

---

## Functions

::: nrcatalogtools.comparisons
    options:
      members:
        - compare_sim_vs_surrogate

---

## Constants

| Name | Value | Description |
|------|-------|-------------|
| `DELTA_T` | `1/4096` s | Default waveform sample interval |
| `DISTANCE` | `1.0` Mpc | Reference distance used internally (amplitude-irrelevant for match) |
