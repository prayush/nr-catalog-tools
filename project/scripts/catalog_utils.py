"""Catalog loading and simulation filtering helpers.

Provides a unified ``load_catalog(name)`` entry point for SXS, RIT, and MAYA
catalogs, and a ``filter_by_surrogate_prior()`` helper that selects simulations
within the NRSur7dq4 validity domain.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import nrcatalogtools as nrcat


def load_catalog(name: str, **kwargs):
    """Load a catalog by name tag.

    Parameters
    ----------
    name : str
        One of ``'SXS'``, ``'RIT'``, or ``'MAYA'`` (case-insensitive).
    **kwargs
        Forwarded to the catalog's ``load()`` class method.

    Returns
    -------
    CatalogBase subclass instance

    Raises
    ------
    ValueError
        If *name* is not one of the three supported catalogs.
    """
    tag = name.upper()
    if tag == "SXS":
        return nrcat.SXSCatalog.load(download=False, **kwargs)
    elif tag == "RIT":
        return nrcat.RITCatalog.load(**kwargs)
    elif tag == "MAYA":
        return nrcat.MayaCatalog.load(**kwargs)
    else:
        raise ValueError(f"Unknown catalog '{name}'. Supported: 'SXS', 'RIT', 'MAYA'.")


def filter_by_surrogate_prior(
    catalog,
    total_mass: float = 40.0,
    q_max: float = 4.0,
    chi_max: float = 0.8,
    verbose: bool = False,
) -> list[str]:
    """Return the subset of simulation names within the NRSur7dq4 prior volume.

    NRSur7dq4 is valid for q ∈ [1, 4] and |χ₁|, |χ₂| ≤ 0.8.  Simulations
    whose parameters cannot be retrieved (metadata errors) are silently skipped.

    Parameters
    ----------
    catalog : CatalogBase
        Loaded catalog object.
    total_mass : float, optional
        Total mass for parameter extraction (default 40 M☉).
    q_max : float, optional
        Maximum allowed mass ratio (default 4).
    chi_max : float, optional
        Maximum allowed spin magnitude (default 0.8).
    verbose : bool, optional
        Print progress to stdout.

    Returns
    -------
    list[str]
        Simulation name tags that pass the prior cuts.
    """
    from surrogate_utils import check_surrogate_prior

    passing = []
    sims = catalog.simulations_list
    for i, sim_name in enumerate(sims):
        if verbose and i % 50 == 0:
            print(f"  Checking {i}/{len(sims)}: {sim_name}")
        try:
            params = catalog.get_parameters(sim_name, total_mass=total_mass)
        except Exception as exc:
            if verbose:
                print(f"  Skipping {sim_name}: {exc}")
            continue
        if check_surrogate_prior(params, q_max=q_max, chi_max=chi_max):
            passing.append(sim_name)
    return passing
