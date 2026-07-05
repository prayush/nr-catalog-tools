"""Unit tests for load_catalog() and filter_by_surrogate_prior() in nrcats."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

import nrcats as nrcat


# ── load_catalog ──────────────────────────────────────────────────────────────


def test_load_catalog_sxs_calls_sxs_load():
    mock_cat = MagicMock()
    with patch.object(nrcat.SXSCatalog, "load", return_value=mock_cat) as mock_load:
        result = nrcat.load_catalog("SXS")
    mock_load.assert_called_once_with(download=False)
    assert result is mock_cat


def test_load_catalog_rit_calls_rit_load():
    mock_cat = MagicMock()
    with patch.object(nrcat.RITCatalog, "load", return_value=mock_cat) as mock_load:
        result = nrcat.load_catalog("RIT")
    mock_load.assert_called_once_with()
    assert result is mock_cat


def test_load_catalog_maya_calls_maya_load():
    mock_cat = MagicMock()
    with patch.object(nrcat.MayaCatalog, "load", return_value=mock_cat) as mock_load:
        result = nrcat.load_catalog("MAYA")
    mock_load.assert_called_once_with()
    assert result is mock_cat


@pytest.mark.parametrize("name", ["sxs", "Sxs", "SXS", "sXs"])
def test_load_catalog_sxs_case_insensitive(name):
    with patch.object(nrcat.SXSCatalog, "load", return_value=MagicMock()) as mock_load:
        nrcat.load_catalog(name)
    mock_load.assert_called_once()


@pytest.mark.parametrize("name", ["rit", "Rit", "RIT"])
def test_load_catalog_rit_case_insensitive(name):
    with patch.object(nrcat.RITCatalog, "load", return_value=MagicMock()) as mock_load:
        nrcat.load_catalog(name)
    mock_load.assert_called_once()


def test_load_catalog_unknown_raises_value_error():
    with pytest.raises(ValueError, match="Unknown catalog"):
        nrcat.load_catalog("LIGO")


def test_load_catalog_unknown_message_lists_supported():
    with pytest.raises(ValueError, match="SXS.*RIT.*MAYA"):
        nrcat.load_catalog("ETDetector")


def test_load_catalog_kwargs_forwarded_to_rit():
    with patch.object(nrcat.RITCatalog, "load", return_value=MagicMock()) as mock_load:
        nrcat.load_catalog("RIT", verbosity=2)
    mock_load.assert_called_once_with(verbosity=2)


# ── filter_by_surrogate_prior ─────────────────────────────────────────────────


def _aligned_params(q=1.5, chi1z=0.3, chi2z=0.0):
    """Minimal PyCBC params dict within the NRSur7dq4 prior."""
    return {
        "mass1": q * 20.0,
        "mass2": 20.0,
        "spin1x": 0.0,
        "spin1y": 0.0,
        "spin1z": chi1z,
        "spin2x": 0.0,
        "spin2y": 0.0,
        "spin2z": chi2z,
        "f_lower": 20.0,
    }


def _make_catalog(sim_names, params_per_sim=None):
    """Build a mock catalog returning fixed params for each sim."""
    mock = MagicMock()
    mock.simulations_list = sim_names
    if params_per_sim is None:
        params_per_sim = {n: _aligned_params() for n in sim_names}
    mock.get_parameters.side_effect = lambda name, **kw: params_per_sim[name]
    return mock


def test_filter_all_pass():
    cat = _make_catalog(["SIM1", "SIM2", "SIM3"])
    result = nrcat.filter_by_surrogate_prior(cat)
    assert result == ["SIM1", "SIM2", "SIM3"]


def test_filter_some_fail_high_q():
    params = {
        "SIM1": _aligned_params(q=1.5),
        "SIM2": _aligned_params(q=5.0),  # outside prior
        "SIM3": _aligned_params(q=2.0),
    }
    cat = _make_catalog(list(params), params)
    result = nrcat.filter_by_surrogate_prior(cat)
    assert result == ["SIM1", "SIM3"]


def test_filter_some_fail_high_spin():
    params = {
        "SIM1": _aligned_params(chi1z=0.3),
        "SIM2": _aligned_params(chi1z=0.9),  # outside prior
    }
    cat = _make_catalog(list(params), params)
    result = nrcat.filter_by_surrogate_prior(cat)
    assert result == ["SIM1"]


def test_filter_metadata_error_silently_skipped():
    """Simulations that raise on get_parameters should be skipped, not crash."""
    mock = MagicMock()
    mock.simulations_list = ["GOOD", "BAD", "ALSO_GOOD"]

    def _get_params(name, **kw):
        if name == "BAD":
            raise KeyError("metadata missing")
        return _aligned_params()

    mock.get_parameters.side_effect = _get_params
    result = nrcat.filter_by_surrogate_prior(mock)
    assert result == ["GOOD", "ALSO_GOOD"]


def test_filter_empty_catalog():
    cat = _make_catalog([])
    assert nrcat.filter_by_surrogate_prior(cat) == []


def test_filter_respects_custom_q_max():
    params = {
        "SIM1": _aligned_params(q=3.0),
        "SIM2": _aligned_params(q=4.0),
    }
    cat = _make_catalog(list(params), params)
    result = nrcat.filter_by_surrogate_prior(cat, q_max=3.5)
    assert result == ["SIM1"]


def test_filter_respects_custom_chi_max():
    params = {
        "SIM1": _aligned_params(chi1z=0.5),
        "SIM2": _aligned_params(chi1z=0.7),
    }
    cat = _make_catalog(list(params), params)
    result = nrcat.filter_by_surrogate_prior(cat, chi_max=0.6)
    assert result == ["SIM1"]


def test_filter_verbose_prints_progress(capsys):
    """verbose=True should print at least one progress line."""
    cat = _make_catalog([f"SIM{i}" for i in range(60)])
    nrcat.filter_by_surrogate_prior(cat, verbose=True)
    out = capsys.readouterr().out
    assert "Checking" in out


def test_filter_uses_total_mass_kwarg():
    """total_mass should be forwarded to get_parameters."""
    cat = _make_catalog(["SIM1"])
    nrcat.filter_by_surrogate_prior(cat, total_mass=80.0)
    cat.get_parameters.assert_called_once_with("SIM1", total_mass=80.0)
