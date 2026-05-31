#!/usr/bin/env python
"""Publication-quality figures for the Step 2 batch comparison results.

Reads ``batch_aligned_all.csv`` (or individual per-catalog CSVs) and produces:

Figure 1a/1b — Mismatch vs source parameters (log y)
  1a: all categories (eccentricity on log x)
  1b: quasi-circular only (categories b+d)

Figure 2a/2b — Phase-difference-per-cycle vs source parameters
  2a: all categories (eccentricity on log x)
  2b: quasi-circular only, log y

Figure 3a — Mismatch CDFs per catalog × category (log x)
Figure 3b — Phase-diff CDFs per catalog × category

Figure 4  — Mismatch 1−ℱ heatmap in (q, χ_eff) space (log colour)

Figure 5  — SXS calibration-split mismatch CDFs (log x)

Usage
-----
    python plot_batch_results.py [--indir results] [--outdir results/figs]
        [--modes 22 21 33 44 32 43] [--cats-only b d]
        [--format png] [--dpi 150]
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPTS_DIR)
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

# ── Aesthetics ────────────────────────────────────────────────────────────────

CATALOG_COLOR = {"SXS": "#2166ac", "RIT": "#d6604d", "MAYA": "#4dac26"}
CATALOG_ZORDER = {"SXS": 3, "RIT": 2, "MAYA": 1}
CATEGORY_MARKER = {"a": "x", "b": "o", "c": "s", "d": "^"}
CATEGORY_LABEL = {
    "a": "ns eccentric",
    "b": "ns quasi-circ",
    "c": "aligned eccentric",
    "d": "aligned quasi-circ",
}
CATEGORY_ALPHA = {"a": 0.55, "b": 0.80, "c": 0.55, "d": 0.80}
MODE_LABEL = {
    "22": r"$(2,2)$",
    "21": r"$(2,1)$",
    "33": r"$(3,3)$",
    "44": r"$(4,4)$",
    "32": r"$(3,2)$",
    "43": r"$(4,3)$",
}
PARAM_LABEL = {
    "q": r"Mass ratio $q$",
    "chi_eff": r"$\chi_{\rm eff}$",
    "spin1z": r"$\chi_{1z}$",
    "spin2z": r"$\chi_{2z}$",
    "eccentricity": r"Eccentricity $e$",
}

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
    }
)


# ── Data helpers ──────────────────────────────────────────────────────────────


def load_results(indir: str) -> pd.DataFrame:
    merged = os.path.join(indir, "batch_aligned_all.csv")
    if os.path.exists(merged):
        df = pd.read_csv(merged)
    else:
        parts = []
        for fname in (
            "batch_aligned_sxs.csv",
            "batch_aligned_rit.csv",
            "batch_aligned_maya.csv",
        ):
            path = os.path.join(indir, fname)
            if os.path.exists(path):
                parts.append(pd.read_csv(path))
        if not parts:
            raise FileNotFoundError(f"No batch result CSVs found in {indir}")
        df = pd.concat(parts, ignore_index=True)

    # Drop failed rows
    if "error" in df.columns:
        df = df[df["error"].isna() | (df["error"] == "")]

    # Derived columns
    if "mismatch_22" not in df.columns and "match_22" in df.columns:
        for m in ("22", "21", "33", "44", "32", "43"):
            col = f"match_{m}"
            if col in df.columns:
                df[f"mismatch_{m}"] = 1.0 - df[col]

    return df


def _good(series: pd.Series) -> pd.Series:
    """Return finite, non-NaN values."""
    return series[np.isfinite(series)]


# ── Per-simulation figures ────────────────────────────────────────────────────

_MODE_ORDER = ["22", "21", "33", "44", "32", "43"]
_MODE_COLORS = ["#2166ac", "#74add1", "#abd9e9", "#fee090", "#f46d43", "#d73027"]

_DELTA_T = 1.0 / 4096
_DISTANCE = 1.0  # Mpc


def plot_individual_sim(
    row: "pd.Series | dict",
    outdir: str,
    fmt: str = "png",
    dpi: int = 150,
    wfm=None,
    h_sur: dict | None = None,
    total_mass: float | None = None,
    delta_t: float = _DELTA_T,
) -> str | None:
    """Four-panel (or two-panel fallback) per-simulation figure.

    Panels when waveforms are supplied:
      1. (2,2) amplitude: NR and NRSur7dq4 overlaid
      2. Re[h_{2,2}] near merger, phase-aligned
      3. Per-mode match bars
      4. Per-mode phase-diff bars

    Falls back to 2-panel (match bars + phase-diff bars) when waveforms are
    absent.  Returns the saved path or None when no mode data is available.
    """
    sim_id = str(row.get("sim_id", "unknown"))
    catalog = str(row.get("catalog", ""))

    if total_mass is None:
        total_mass = float(row.get("total_mass", 40.0))

    safe_id = sim_id.replace(":", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"{safe_id}_mode_matches.{fmt}")

    matches = [float(row.get(f"match_{m}", np.nan)) for m in _MODE_ORDER]
    phase_diffs = [float(row.get(f"phase_diff_{m}", np.nan)) for m in _MODE_ORDER]

    if all(np.isnan(v) for v in matches):
        return None

    mode_labels = [MODE_LABEL.get(m, m) for m in _MODE_ORDER]
    x = np.arange(len(_MODE_ORDER))
    has_wfm = wfm is not None and h_sur is not None

    if has_wfm:
        fig, axes = plt.subplots(
            4,
            1,
            figsize=(9, 12),
            gridspec_kw={"height_ratios": [1.8, 1.8, 1.0, 1.0]},
        )
        ax_amp, ax_re, ax_match, ax_phase = axes
    else:
        fig, (ax_match, ax_phase) = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True)

    # ── Panel 1: (2,2) amplitude ──────────────────────────────────────────
    if has_wfm:
        h22_nr = None
        h22_sur = h_sur.get((2, 2)) if isinstance(h_sur, dict) else None
        try:
            h22_nr = wfm.get_mode(
                2,
                2,
                total_mass=total_mass,
                distance=_DISTANCE,
                delta_t_seconds=delta_t,
            )
        except Exception:
            ax_amp.text(
                0.5,
                0.5,
                "NR (2,2) unavailable",
                transform=ax_amp.transAxes,
                ha="center",
                va="center",
                fontsize=7,
                color="gray",
            )
            h22_nr = None

        if h22_nr is not None and h22_sur is not None:
            t_start = max(float(h22_nr.start_time), float(h22_sur.start_time))
            t_end = min(float(h22_nr.end_time), float(h22_sur.end_time))
            if t_end > t_start:
                h22_nr = h22_nr.time_slice(t_start, t_end)
                h22_sur = h22_sur.time_slice(t_start, t_end)

        if h22_nr is not None:
            t_nr = np.array(h22_nr.sample_times)
            amp_nr = np.abs(np.array(h22_nr))
            ax_amp.plot(t_nr, amp_nr, color="C0", lw=1.5, label="NR", alpha=0.9)

        if h22_sur is not None:
            t_sur = np.array(h22_sur.sample_times)
            amp_sur = np.abs(np.array(h22_sur))
            ax_amp.plot(
                t_sur,
                amp_sur,
                color="C1",
                lw=1.5,
                ls="--",
                label="NRSur7dq4",
                alpha=0.9,
            )

        ax_amp.set_ylabel(r"$|h_{22}|$ [strain]", fontsize=9)
        ax_amp.set_title(r"$(2,2)$ mode amplitude", fontsize=9)
        ax_amp.legend(fontsize=8)
        if h22_nr is not None:
            ax_amp.set_xlim(left=max(-2.0, float(t_nr[0])))
        ax_amp.grid(True, lw=0.4, alpha=0.4)
        ax_amp.tick_params(labelsize=7)

        # ── Panel 2: Re[h_{2,2}] near merger ─────────────────────────────
        if h22_nr is not None and h22_sur is not None:
            try:
                t_nr_full = np.array(h22_nr.sample_times)
                re_nr = np.array(h22_nr.real())
                mask_nr = t_nr_full >= -0.25
                ax_re.plot(
                    t_nr_full[mask_nr], re_nr[mask_nr], color="C0", lw=1.2, label="NR"
                )

                # Phase-align surrogate at peak
                idx_nr_peak = int(np.argmax(np.abs(np.array(h22_nr))))
                idx_sur_peak = int(np.argmax(np.abs(np.array(h22_sur))))
                phase_nr_at_peak = np.angle(np.array(h22_nr)[idx_nr_peak])
                phase_sur_at_peak = np.angle(np.array(h22_sur)[idx_sur_peak])
                h22_sur_aligned = np.array(h22_sur) * np.exp(
                    1j * (phase_nr_at_peak - phase_sur_at_peak)
                )

                t_sur_full = np.array(h22_sur.sample_times)
                mask_sur = t_sur_full >= -0.25
                ax_re.plot(
                    t_sur_full[mask_sur],
                    h22_sur_aligned[mask_sur].real,
                    color="C1",
                    lw=1.2,
                    ls="--",
                    label="NRSur7dq4 (phase-aligned)",
                )
            except Exception:
                ax_re.text(
                    0.5,
                    0.5,
                    "alignment failed",
                    transform=ax_re.transAxes,
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="gray",
                )
        else:
            ax_re.text(
                0.5,
                0.5,
                "waveform data unavailable",
                transform=ax_re.transAxes,
                ha="center",
                va="center",
                fontsize=7,
                color="gray",
            )

        ax_re.set_ylabel(r"$\mathrm{Re}[h_{22}]$ [strain]", fontsize=9)
        ax_re.set_title(r"$(2,2)$ real part — merger region", fontsize=9)
        ax_re.legend(fontsize=8)
        ax_re.grid(True, lw=0.4, alpha=0.4)
        ax_re.tick_params(labelsize=7)

    # ── Panel 3 (or 1): match bars ────────────────────────────────────────
    bars1 = ax_match.bar(
        x, matches, color=_MODE_COLORS, edgecolor="k", linewidth=0.5, width=0.65
    )
    ax_match.axhline(1.0, color="k", lw=0.5, ls="--", alpha=0.5)
    ax_match.set_ylim(0, 1.08)
    ax_match.set_ylabel(r"Match $\mathcal{F}$", fontsize=9)
    ax_match.grid(True, axis="y", lw=0.4, alpha=0.4)
    ax_match.tick_params(labelsize=7)
    if not has_wfm:
        ax_match.set_xticks([])
    for bar, v in zip(bars1, matches):
        if np.isfinite(v):
            ax_match.text(
                bar.get_x() + bar.get_width() / 2,
                min(v + 0.015, 1.04),
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
            )

    # ── Panel 4 (or 2): phase-diff bars ──────────────────────────────────
    bars2 = ax_phase.bar(
        x, phase_diffs, color=_MODE_COLORS, edgecolor="k", linewidth=0.5, width=0.65
    )
    ax_phase.set_ylabel(r"$\Delta\Phi/{\rm cycle}$ [rad]", fontsize=9)
    ax_phase.set_xticks(x)
    ax_phase.set_xticklabels(mode_labels, fontsize=8.5)
    ax_phase.grid(True, axis="y", lw=0.4, alpha=0.4)
    ax_phase.tick_params(labelsize=7)
    max_pd = max((v for v in phase_diffs if np.isfinite(v)), default=0)
    for bar, v in zip(bars2, phase_diffs):
        if np.isfinite(v) and v > 0:
            ax_phase.text(
                bar.get_x() + bar.get_width() / 2,
                v + 0.005 * max(max_pd, 0.05),
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
            )

    # ── Title ─────────────────────────────────────────────────────────────
    q = row.get("q", np.nan)
    s1z = row.get("spin1z", np.nan)
    s2z = row.get("spin2z", np.nan)
    ecc = row.get("eccentricity", np.nan)
    M = row.get("total_mass", np.nan)
    is_cal = row.get("nrsur7dq4_calibration", False)
    cal_tag = " [NRSur cal.]" if is_cal else ""
    fig.suptitle(
        f"{catalog}:{sim_id}{cal_tag}  |  "
        f"$q={float(q):.2f}$, $\\chi_{{1z}}={float(s1z):.2f}$, "
        f"$\\chi_{{2z}}={float(s2z):.2f}$, $e={float(ecc):.3f}$, "
        f"$M={float(M):.0f}\\,M_\\odot$",
        fontsize=8.5,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return outpath


def make_individual_sim_figures(
    df: pd.DataFrame,
    outdir: str,
    fmt: str = "png",
    dpi: int = 120,
    skip_existing: bool = False,
) -> None:
    """Generate one 4-panel figure per simulation row from a batch-results DataFrame.

    Loads NR waveforms and generates surrogate modes on-the-fly from the
    parameters stored in the CSV.  Falls back to 2-panel (bar charts only) if
    waveform loading fails for a given simulation.
    """
    os.makedirs(outdir, exist_ok=True)
    n_saved = n_skipped = n_failed = 0

    # Pre-load surrogate and catalog utilities once
    _wfm_support = False
    _generate_sur = None
    _load_catalog = None
    try:
        from surrogate_utils import generate_surrogate_modes, load_nrsur7dq4
        from catalog_utils import load_catalog

        print("  Pre-loading NRSur7dq4 surrogate model ...")
        load_nrsur7dq4()
        _generate_sur = generate_surrogate_modes
        _load_catalog = load_catalog
        _wfm_support = True
        print("  Surrogate loaded.")
    except Exception as _e:
        print(
            f"  Warning: waveform utilities unavailable ({_e}); "
            "generating 2-panel figures from CSV data only."
        )

    cat_cache: dict = {}

    for idx, row in df.iterrows():
        err = row.get("error", "")
        if err and str(err) not in ("", "nan"):
            continue

        sim_id = str(row.get("sim_id", ""))
        safe_id = sim_id.replace(":", "_").replace("/", "_")
        outpath = os.path.join(outdir, f"{safe_id}_mode_matches.{fmt}")
        if skip_existing and os.path.exists(outpath):
            n_skipped += 1
            continue

        total_mass = float(row.get("total_mass", 40.0))
        wfm = None
        h_sur = None

        if _wfm_support:
            try:
                catalog_name = str(row.get("catalog", ""))
                if catalog_name not in cat_cache:
                    cat_cache[catalog_name] = _load_catalog(catalog_name)
                cat = cat_cache[catalog_name]
                wfm = cat.get(sim_id)

                # Reconstruct params from CSV (aligned-spin: perp components = 0)
                params = {
                    "mass1": float(row["mass1"]),
                    "mass2": float(row["mass2"]),
                    "spin1x": 0.0,
                    "spin1y": 0.0,
                    "spin1z": float(row.get("spin1z", 0.0)),
                    "spin2x": 0.0,
                    "spin2y": 0.0,
                    "spin2z": float(row.get("spin2z", 0.0)),
                    "f_lower": float(row.get("f_lower_nr", 20.0)),
                }
                h_sur, _ = _generate_sur(
                    params,
                    total_mass=total_mass,
                    distance=_DISTANCE,
                    delta_t_seconds=_DELTA_T,
                )
            except Exception:
                wfm = None
                h_sur = None

        try:
            result = plot_individual_sim(
                row,
                outdir,
                fmt=fmt,
                dpi=dpi,
                wfm=wfm,
                h_sur=h_sur,
                total_mass=total_mass,
            )
            if result:
                n_saved += 1
            else:
                n_failed += 1
        except Exception:
            n_failed += 1

    print(
        f"  Individual figures: {n_saved} new, {n_skipped} already exist, "
        f"{n_failed} failed  →  {outdir}"
    )


# ── Figure 1 & 2: scatter plots ───────────────────────────────────────────────


def _scatter_grid(
    df: pd.DataFrame,
    y_col_template: str,
    y_label_template: str,
    modes: list[str],
    params: list[str],
    outpath: str,
    y_lim: tuple | None = None,
    log_y: bool = False,
    log_x_cols: list[str] | None = None,
    categories: list[str] | None = None,
):
    cats_to_plot = categories or ["a", "b", "c", "d"]
    log_x_cols = log_x_cols or []
    n_modes = len(modes)
    n_params = len(params)

    fig, axes = plt.subplots(
        n_modes,
        n_params,
        figsize=(2.6 * n_params, 2.2 * n_modes),
        squeeze=False,
        sharex="col",
    )

    for col_idx, param in enumerate(params):
        use_log_x = param in log_x_cols
        for row_idx, mode in enumerate(modes):
            ax = axes[row_idx][col_idx]
            y_col = y_col_template.format(mode=mode)

            if y_col not in df.columns:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="gray",
                )
                continue

            for catname in ["MAYA", "RIT", "SXS"]:
                if catname not in df["catalog"].values:
                    continue
                sub = df[df["catalog"] == catname]
                for cat_key in cats_to_plot:
                    mask = sub["category"] == cat_key
                    if mask.sum() == 0:
                        continue
                    xv = sub.loc[mask, param]
                    yv = sub.loc[mask, y_col]
                    valid = np.isfinite(xv) & np.isfinite(yv)
                    if use_log_x:
                        valid &= xv > 0
                    if log_y:
                        valid &= yv > 0
                    if valid.sum() == 0:
                        continue
                    ax.scatter(
                        xv[valid],
                        yv[valid],
                        c=CATALOG_COLOR[catname],
                        marker=CATEGORY_MARKER[cat_key],
                        s=12,
                        linewidths=0.4,
                        alpha=CATEGORY_ALPHA[cat_key],
                        zorder=CATALOG_ZORDER[catname],
                        edgecolors="none"
                        if CATEGORY_MARKER[cat_key] == "o"
                        else CATALOG_COLOR[catname],
                    )

            if use_log_x:
                ax.set_xscale("log")
            if log_y:
                ax.set_yscale("log")
            if y_lim is not None:
                ax.set_ylim(y_lim)
            if row_idx == n_modes - 1:
                ax.set_xlabel(PARAM_LABEL.get(param, param), fontsize=9)
            if col_idx == 0:
                ax.set_ylabel(
                    y_label_template + f"\n{MODE_LABEL.get(mode, mode)}", fontsize=8
                )
            ax.grid(True, lw=0.4, alpha=0.4)
            ax.tick_params(labelsize=7)

    # shared legend
    legend_handles = []
    for catname, color in CATALOG_COLOR.items():
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker="o",
                markersize=5,
                linestyle="none",
                label=catname,
            )
        )
    for cat_key in cats_to_plot:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color="gray",
                marker=CATEGORY_MARKER[cat_key],
                markersize=5,
                linestyle="none",
                label=CATEGORY_LABEL[cat_key],
            )
        )
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(legend_handles),
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ── Figure 3: CDFs ────────────────────────────────────────────────────────────


def _cdf_figure(
    df: pd.DataFrame,
    modes: list[str],
    metric_col: str,
    x_label: str,
    outpath: str,
    categories: list[str] | None = None,
    x_lim: tuple | None = None,
    log_x: bool = False,
):
    cats_to_plot = categories or ["b", "d"]
    n = len(modes)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.5 * ncols, 3.0 * nrows), squeeze=False
    )

    cat_linestyle = {"b": "-", "d": "--", "a": ":", "c": "-."}

    for idx, mode in enumerate(modes):
        ax = axes[idx // ncols][idx % ncols]
        col = metric_col.format(mode=mode)

        if col not in df.columns:
            ax.set_visible(False)
            continue

        plotted = False
        for catname in ["SXS", "RIT", "MAYA"]:
            if catname not in df["catalog"].values:
                continue
            for cat_key in cats_to_plot:
                mask = (df["catalog"] == catname) & (df["category"] == cat_key)
                vals = _good(df.loc[mask, col])
                if len(vals) < 3:
                    continue
                if log_x:
                    vals = vals[vals > 0]
                if len(vals) < 3:
                    continue
                vals_sorted = np.sort(vals)
                cdf = np.arange(1, len(vals_sorted) + 1) / len(vals_sorted)
                ax.plot(
                    vals_sorted,
                    cdf,
                    color=CATALOG_COLOR[catname],
                    linestyle=cat_linestyle[cat_key],
                    lw=1.5,
                    alpha=0.85,
                    label=f"{catname}/{CATEGORY_LABEL[cat_key]}",
                )
                plotted = True

        ax.set_title(MODE_LABEL.get(mode, mode), fontsize=9)
        ax.set_ylabel("CDF", fontsize=8)
        ax.set_xlabel(x_label, fontsize=8)
        if log_x:
            ax.set_xscale("log")
        if x_lim is not None:
            ax.set_xlim(x_lim)
        ax.set_ylim(0, 1.05)
        ax.grid(True, lw=0.4, alpha=0.4)
        ax.tick_params(labelsize=7)
        if plotted and idx == 0:
            ax.legend(fontsize=6.5, loc="upper left", ncol=2, framealpha=0.7)

    for idx in range(len(modes), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ── Figure 4: (q, χ_eff) heatmaps ─────────────────────────────────────────────


def _heatmap_figure(
    df: pd.DataFrame,
    mode: str,
    metric_col: str,
    title: str,
    outpath: str,
    vmin: float = 0.0,
    vmax: float = 0.1,
    categories: list[str] | None = None,
    log_color: bool = False,
):
    cats = categories or ["b", "d"]
    sub = df[df["category"].isin(cats)].copy()
    col = metric_col.format(mode=mode)
    if col not in sub.columns:
        return

    sub = sub[
        np.isfinite(sub["q"]) & np.isfinite(sub["chi_eff"]) & np.isfinite(sub[col])
    ]
    if log_color:
        sub = sub[sub[col] > 0]

    catalogs = [c for c in ["SXS", "RIT", "MAYA"] if c in sub["catalog"].values]
    fig, axes = plt.subplots(
        1, len(catalogs), figsize=(4.5 * len(catalogs), 3.8), squeeze=False
    )

    # For mismatch: low = good = green; use reversed colormap
    cmap = plt.cm.RdYlGn_r if log_color else plt.cm.RdYlGn
    if log_color:
        norm = mcolors.LogNorm(vmin=max(vmin, 1e-6), vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    for ax, catname in zip(axes[0], catalogs):
        s = sub[sub["catalog"] == catname]
        sc = ax.scatter(
            s["q"],
            s["chi_eff"],
            c=s[col],
            cmap=cmap,
            norm=norm,
            s=25,
            edgecolors="k",
            linewidths=0.3,
            alpha=0.85,
        )
        ax.set_xlabel(r"Mass ratio $q$", fontsize=10)
        ax.set_ylabel(r"$\chi_{\rm eff}$", fontsize=10)
        ax.set_title(catname, fontsize=10)
        ax.grid(True, lw=0.4, alpha=0.4)
        plt.colorbar(sc, ax=ax, label=title)

    fig.suptitle(
        f"{MODE_LABEL.get(mode, mode)} — {title}  (categories: {', '.join(cats)})",
        fontsize=10,
        y=1.02,
    )
    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ── Summary stats helper ──────────────────────────────────────────────────────


def _print_summary(df: pd.DataFrame, modes: list[str]):
    print("\n=== Mismatch summary (median / 90th pctile) for quasi-circular (b+d) ===")
    qc = df[df["category"].isin(["b", "d"])]
    hdr = f"  {'Mode':<8}" + "".join(
        f"  {'':>6}{cat:<12}" for cat in ["SXS", "RIT", "MAYA"]
    )
    print(hdr)
    for mode in modes:
        col = f"mismatch_{mode}"
        if col not in qc.columns:
            col = f"match_{mode}"
            if col not in qc.columns:
                continue
        row_str = f"  {MODE_LABEL.get(mode, mode):<8}"
        for catname in ["SXS", "RIT", "MAYA"]:
            vals = _good(qc[qc["catalog"] == catname][col])
            if len(vals) == 0:
                row_str += f"  {'—':>18}"
            else:
                med = np.median(vals)
                p90 = np.percentile(vals, 90)
                row_str += f"  {med:.4f} / {p90:.4f}"
        print(row_str)


# ── SXS calibration-split CDF ─────────────────────────────────────────────────


def _cdf_calibration_split_figure(
    df: pd.DataFrame,
    modes: list[str],
    metric_col: str,
    x_label: str,
    outpath: str,
    categories: list[str] | None = None,
    x_lim: tuple | None = None,
    log_x: bool = False,
):
    """Mismatch CDFs where SXS is split into NRSur calibration vs non-calibration."""
    cats_to_plot = categories or ["b", "d"]
    has_cal_col = "nrsur7dq4_calibration" in df.columns

    n = len(modes)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.5 * ncols, 3.0 * nrows), squeeze=False
    )

    for idx, mode in enumerate(modes):
        ax = axes[idx // ncols][idx % ncols]
        col = metric_col.format(mode=mode)
        if col not in df.columns:
            ax.set_visible(False)
            continue

        plotted = False

        # SXS — split by calibration flag
        sxs = df[(df["catalog"] == "SXS") & df["category"].isin(cats_to_plot)]
        if has_cal_col:
            for is_cal, ls, suffix in [
                (True, "-", " (NRSur cal.)"),
                (False, "--", " (non-cal.)"),
            ]:
                mask = sxs["nrsur7dq4_calibration"].astype(bool) == is_cal
                vals = _good(sxs.loc[mask, col])
                if log_x:
                    vals = vals[vals > 0]
                if len(vals) >= 3:
                    vs = np.sort(vals)
                    cdf = np.arange(1, len(vs) + 1) / len(vs)
                    ax.plot(
                        vs,
                        cdf,
                        color=CATALOG_COLOR["SXS"],
                        ls=ls,
                        lw=1.5,
                        alpha=0.9,
                        label=f"SXS{suffix} (N={len(vs)})",
                    )
                    plotted = True
        else:
            vals = _good(sxs[col])
            if log_x:
                vals = vals[vals > 0]
            if len(vals) >= 3:
                vs = np.sort(vals)
                cdf = np.arange(1, len(vs) + 1) / len(vs)
                ax.plot(vs, cdf, color=CATALOG_COLOR["SXS"], lw=1.5, label="SXS")
                plotted = True

        for catname in ["RIT", "MAYA"]:
            mask = (df["catalog"] == catname) & df["category"].isin(cats_to_plot)
            vals = _good(df.loc[mask, col])
            if log_x:
                vals = vals[vals > 0]
            if len(vals) >= 3:
                vs = np.sort(vals)
                cdf = np.arange(1, len(vs) + 1) / len(vs)
                ax.plot(
                    vs,
                    cdf,
                    color=CATALOG_COLOR[catname],
                    ls="-",
                    lw=1.5,
                    alpha=0.9,
                    label=f"{catname} (N={len(vs)})",
                )
                plotted = True

        ax.set_title(MODE_LABEL.get(mode, mode), fontsize=9)
        ax.set_ylabel("CDF", fontsize=8)
        ax.set_xlabel(x_label, fontsize=8)
        if log_x:
            ax.set_xscale("log")
        if x_lim is not None:
            ax.set_xlim(x_lim)
        ax.set_ylim(0, 1.05)
        ax.grid(True, lw=0.4, alpha=0.4)
        ax.tick_params(labelsize=7)
        if plotted and idx == 0:
            ax.legend(fontsize=6.0, loc="upper left", framealpha=0.75)

    for idx in range(len(modes), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ── Main ──────────────────────────────────────────────────────────────────────


def make_all_figures(
    indir: str,
    outdir: str,
    modes: list[str],
    categories: list[str],
    fmt: str = "png",
    dpi: int = 150,
):
    os.makedirs(outdir, exist_ok=True)
    plt.rcParams["figure.dpi"] = dpi

    print(f"Loading results from {indir} ...")
    df = load_results(indir)
    print(f"  {len(df)} rows loaded from {df['catalog'].nunique()} catalog(s).")
    if df.empty:
        print("  No data found — check that batch_aligned_all.csv exists.")
        return

    _print_summary(df, modes)

    params_all = ["q", "chi_eff", "spin1z", "spin2z", "eccentricity"]
    params_qc = ["q", "chi_eff", "spin1z", "spin2z"]
    qc_cats = [c for c in categories if c in ("b", "d")]

    # ── Figure 1a: mismatch vs all params (log y, eccentricity log x) ────
    print("\nGenerating Figure 1a — mismatch vs source params (all categories)...")
    _scatter_grid(
        df,
        y_col_template="mismatch_{mode}",
        y_label_template=r"Mismatch $1{-}\mathcal{F}$",
        modes=modes,
        params=params_all,
        outpath=os.path.join(outdir, f"fig1a_mismatch_vs_params.{fmt}"),
        y_lim=(1e-4, 1.0),
        log_y=True,
        log_x_cols=["eccentricity"],
        categories=categories,
    )

    # ── Figure 1b: mismatch vs qc params (log y) ─────────────────────────
    if qc_cats:
        print(
            "Generating Figure 1b — mismatch vs source params (quasi-circular b+d)..."
        )
        _scatter_grid(
            df[df["category"].isin(qc_cats)],
            y_col_template="mismatch_{mode}",
            y_label_template=r"Mismatch $1{-}\mathcal{F}$",
            modes=modes,
            params=params_qc,
            outpath=os.path.join(outdir, f"fig1b_mismatch_vs_params_qc.{fmt}"),
            y_lim=(1e-4, 1.0),
            log_y=True,
            categories=qc_cats,
        )

    # ── Figure 2a: phase diff vs all params (eccentricity log x) ─────────
    print("Generating Figure 2a — phase diff vs source params...")
    _scatter_grid(
        df,
        y_col_template="phase_diff_{mode}",
        y_label_template=r"$\Delta\Phi/\mathrm{cycle}$ [rad]",
        modes=modes,
        params=params_all,
        outpath=os.path.join(outdir, f"fig2a_phasediff_vs_params.{fmt}"),
        y_lim=(0, None),
        log_y=False,
        log_x_cols=["eccentricity"],
        categories=categories,
    )

    # ── Figure 2b: log phase diff vs qc params ────────────────────────────
    if qc_cats:
        print("Generating Figure 2b — log phase diff vs source params (qc)...")
        _scatter_grid(
            df[df["category"].isin(qc_cats)],
            y_col_template="phase_diff_{mode}",
            y_label_template=r"$\Delta\Phi/\mathrm{cycle}$ [rad]",
            modes=modes,
            params=params_qc,
            outpath=os.path.join(outdir, f"fig2b_phasediff_vs_params_qc.{fmt}"),
            y_lim=(1e-4, 10),
            log_y=True,
            categories=qc_cats,
        )

    # ── Figure 3a: mismatch CDFs (log x) ─────────────────────────────────
    print("Generating Figure 3a — mismatch CDFs (log x)...")
    _cdf_figure(
        df,
        modes=modes,
        metric_col="mismatch_{mode}",
        x_label=r"Mismatch $1 - \mathcal{F}$",
        outpath=os.path.join(outdir, f"fig3a_mismatch_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(1e-4, 1.0),
        log_x=True,
    )

    # ── Figure 3b: phase-diff CDFs ────────────────────────────────────────
    print("Generating Figure 3b — phase diff CDFs per mode...")
    _cdf_figure(
        df,
        modes=modes,
        metric_col="phase_diff_{mode}",
        x_label=r"$\Delta\Phi/\mathrm{cycle}$ [rad]",
        outpath=os.path.join(outdir, f"fig3b_phasediff_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(0, 2.0),
    )

    # ── Figure 4: (q, χ_eff) mismatch heatmaps ───────────────────────────
    print("Generating Figure 4 — (q, χ_eff) mismatch heatmaps (qc)...")
    for mode in modes:
        _heatmap_figure(
            df,
            mode=mode,
            metric_col="mismatch_{mode}",
            title=rf"Mismatch $1{{-}}\mathcal{{F}}_{{{mode}}}$",
            outpath=os.path.join(outdir, f"fig4_mismatch{mode}_heatmap_qc.{fmt}"),
            vmin=1e-4,
            vmax=1.0,
            log_color=True,
            categories=qc_cats or categories,
        )

    # ── Figure 1c: mismatch vs eccentricity ─────────────────────────────
    print("Generating Figure 1c — mismatch vs eccentricity...")
    _scatter_grid(
        df,
        y_col_template="mismatch_{mode}",
        y_label_template=r"Mismatch $1{-}\mathcal{F}$",
        modes=modes,
        params=["eccentricity"],
        outpath=os.path.join(outdir, f"fig1c_mismatch_vs_eccentricity.{fmt}"),
        y_lim=(1e-4, 1.0),
        log_y=True,
        log_x_cols=["eccentricity"],
        categories=categories,
    )

    # ── Figure 5: SXS calibration-split mismatch CDFs (log x) ────────────
    print("Generating Figure 5 — SXS cal-split mismatch CDFs (log)...")
    _cdf_calibration_split_figure(
        df,
        modes=modes,
        metric_col="mismatch_{mode}",
        x_label=r"Mismatch $1 - \mathcal{F}$",
        outpath=os.path.join(outdir, f"fig5_sxs_cal_mismatch_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(1e-4, 1.0),
        log_x=True,
    )

    # ── Individual simulation figures ─────────────────────────────────────
    indiv_dir = os.path.join(outdir, "individual_sims")
    print(f"\nGenerating individual simulation figures → {indiv_dir} ...")
    make_individual_sim_figures(
        df, outdir=indiv_dir, fmt=fmt, dpi=dpi, skip_existing=False
    )

    print(f"\nAll figures saved to {outdir}/")


def _build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--indir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "results")),
        help="Directory containing batch CSV files (default: project/results/)",
    )
    p.add_argument(
        "--outdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "figs")),
        help="Output directory for figures (default: project/figs/)",
    )
    p.add_argument(
        "--modes",
        nargs="+",
        default=["22", "21", "33", "44", "32", "43"],
        help="Mode tags to include (default: 22 21 33 44 32 43)",
    )
    p.add_argument(
        "--categories",
        nargs="+",
        default=["a", "b", "c", "d"],
        help="Category keys to include in scatter plots",
    )
    p.add_argument(
        "--format",
        default="png",
        choices=["png", "pdf", "svg"],
        help="Output figure format (default: png)",
    )
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument(
        "--no-indiv",
        action="store_true",
        help="Skip per-simulation figure generation (fast summary only)",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    if args.no_indiv:
        # Patch make_all_figures to skip individual figures
        _orig = make_individual_sim_figures

        def _noop(*a, **kw):
            print("  --no-indiv: skipping per-simulation figures.")

        make_individual_sim_figures = _noop  # noqa: F811

    make_all_figures(
        indir=args.indir,
        outdir=args.outdir,
        modes=args.modes,
        categories=args.categories,
        fmt=args.format,
        dpi=args.dpi,
    )
