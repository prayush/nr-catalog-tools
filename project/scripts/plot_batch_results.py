#!/usr/bin/env python
"""Publication-quality figures for the Step 2 batch comparison results.

Reads ``batch_aligned_all.csv`` (or individual per-catalog CSVs) and produces:

Figure 1 — Match vs source parameters  (one column per parameter)
  Rows: dominant (2,2) mode and each sub-dominant mode
  Columns: q, χ_eff, χ₁z, χ₂z, eccentricity e
  Color-coded by catalog; category shown via marker shape

Figure 2 — Phase-difference-per-cycle vs source parameters
  Same layout as Figure 1

Figure 3 — Match CDFs per catalog × category
  One panel per mode; CDF lines colored by catalog, linestyle by category

Figure 4 — Mismatch 1−ℱ heatmap in (q, χ_eff) space
  One sub-panel per catalog for the (2,2) mode; only quasi-circular (b+d)

Usage
-----
    python plot_batch_results.py [--indir results] [--outdir results/figs]
        [--modes 22 21 33 44 32 43] [--cats-only b d]
        [--format png] [--dpi 150]
"""

from __future__ import annotations

import argparse
import os
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

    # Drop failed rows (error column non-empty)
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


def plot_individual_sim(
    row: "pd.Series | dict",
    outdir: str,
    fmt: str = "png",
    dpi: int = 150,
) -> str | None:
    """Two-panel per-simulation figure: match bar chart + phase-diff bar chart.

    Returns the saved file path, or None when no mode data is available.
    """
    sim_id = str(row.get("sim_id", "unknown"))
    catalog = str(row.get("catalog", ""))

    safe_id = sim_id.replace(":", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"{safe_id}_mode_matches.{fmt}")

    matches = [float(row.get(f"match_{m}", np.nan)) for m in _MODE_ORDER]
    phase_diffs = [float(row.get(f"phase_diff_{m}", np.nan)) for m in _MODE_ORDER]

    if all(np.isnan(v) for v in matches):
        return None

    mode_labels = [MODE_LABEL.get(m, m) for m in _MODE_ORDER]
    x = np.arange(len(_MODE_ORDER))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True)

    # ── match bars ────────────────────────────────────────────────────────────
    bars1 = ax1.bar(
        x, matches, color=_MODE_COLORS, edgecolor="k", linewidth=0.5, width=0.65
    )
    ax1.axhline(1.0, color="k", lw=0.5, ls="--", alpha=0.5)
    ax1.set_ylim(0, 1.08)
    ax1.set_ylabel(r"Match $\mathcal{F}$", fontsize=9)
    ax1.grid(True, axis="y", lw=0.4, alpha=0.4)
    for bar, v in zip(bars1, matches):
        if np.isfinite(v):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                min(v + 0.015, 1.04),
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
            )

    # ── phase-diff bars ───────────────────────────────────────────────────────
    bars2 = ax2.bar(
        x, phase_diffs, color=_MODE_COLORS, edgecolor="k", linewidth=0.5, width=0.65
    )
    ax2.set_ylabel(r"$\Delta\Phi/{\rm cycle}$ [rad]", fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(mode_labels, fontsize=8.5)
    ax2.grid(True, axis="y", lw=0.4, alpha=0.4)
    max_pd = max((v for v in phase_diffs if np.isfinite(v)), default=0)
    for bar, v in zip(bars2, phase_diffs):
        if np.isfinite(v) and v > 0:
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                v + 0.005 * max(max_pd, 0.05),
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
            )

    # ── title ─────────────────────────────────────────────────────────────────
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

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return outpath


def make_individual_sim_figures(
    df: pd.DataFrame,
    outdir: str,
    fmt: str = "png",
    dpi: int = 120,
    skip_existing: bool = True,
) -> None:
    """Generate one figure per simulation row from a batch-results DataFrame."""
    os.makedirs(outdir, exist_ok=True)
    n_saved = n_skipped = n_failed = 0

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
        try:
            result = plot_individual_sim(row, outdir, fmt=fmt, dpi=dpi)
            if result:
                n_saved += 1
        except Exception:
            n_failed += 1

    print(
        f"  Individual figures: {n_saved} new, {n_skipped} already exist, "
        f"{n_failed} failed  →  {outdir}"
    )


# ── Figure 1 & 2: scatter plots ───────────────────────────────────────────────


def _scatter_grid(
    df: pd.DataFrame,
    y_col_template: str,  # e.g. "match_{mode}" or "phase_diff_{mode}"
    y_label_template: str,  # e.g. "Match $\\mathcal{{F}}$" or "ΔΦ/cycle [rad]"
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

            for catname in ["MAYA", "RIT", "SXS"]:  # back-to-front z-order
                if catname not in df["catalog"].values:
                    continue
                sub = df[df["catalog"] == catname]
                for cat_key in cats_to_plot:
                    mask = sub["category"] == cat_key
                    if mask.sum() == 0:
                        continue
                    xv = sub.loc[mask, param]
                    yv = sub.loc[mask, y_col]
                    # For log-x columns filter non-positive x values
                    if use_log_x:
                        valid = np.isfinite(xv) & np.isfinite(yv) & (xv > 0)
                    else:
                        valid = np.isfinite(xv) & np.isfinite(yv)
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
                    y_label_template + f"\n{MODE_LABEL.get(mode, mode)}",
                    fontsize=8,
                )
            ax.grid(True, lw=0.4, alpha=0.4)
            ax.tick_params(labelsize=7)

    # ── shared legend ────────────────────────────────────────────────────────
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

    # hide unused panels
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
    vmin: float = 0.95,
    vmax: float = 1.0,
    categories: list[str] | None = None,
):
    cats = categories or ["b", "d"]
    sub = df[df["category"].isin(cats)].copy()
    col = metric_col.format(mode=mode)
    if col not in sub.columns:
        return

    sub = sub[
        np.isfinite(sub["q"]) & np.isfinite(sub["chi_eff"]) & np.isfinite(sub[col])
    ]

    catalogs = [c for c in ["SXS", "RIT", "MAYA"] if c in sub["catalog"].values]
    fig, axes = plt.subplots(
        1, len(catalogs), figsize=(4.5 * len(catalogs), 3.8), squeeze=False
    )

    cmap = plt.cm.RdYlGn
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
    print("\n=== Match summary (median / 10th pctile) for quasi-circular (b+d) ===")
    qc = df[df["category"].isin(["b", "d"])]
    hdr = f"  {'Mode':<8}" + "".join(
        f"  {'':>6}{cat:<12}" for cat in ["SXS", "RIT", "MAYA"]
    )
    print(hdr)
    for mode in modes:
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
                p10 = np.percentile(vals, 10)
                row_str += f"  {med:.4f} / {p10:.4f}"
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
    """CDFs where SXS is split into NRSur calibration vs non-calibration subsets.

    SXS calibration sims are shown as solid blue; SXS non-calibration as dashed
    blue; RIT and MAYA appear as solid lines in their respective catalog colors.
    """
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
            if len(vals) >= 3:
                vs = np.sort(vals)
                cdf = np.arange(1, len(vs) + 1) / len(vs)
                ax.plot(vs, cdf, color=CATALOG_COLOR["SXS"], lw=1.5, label="SXS")
                plotted = True

        # RIT and MAYA (all sims together)
        for catname in ["RIT", "MAYA"]:
            mask = (df["catalog"] == catname) & df["category"].isin(cats_to_plot)
            vals = _good(df.loc[mask, col])
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

    # ── Figure 1a: match vs all params, eccentricity on log x ────────────────
    print("\nGenerating Figure 1a — match vs source params (all categories)...")
    _scatter_grid(
        df,
        y_col_template="match_{mode}",
        y_label_template=r"Match $\mathcal{F}$",
        modes=modes,
        params=params_all,
        outpath=os.path.join(outdir, f"fig1a_match_vs_params.{fmt}"),
        y_lim=(0.0, 1.02),
        log_x_cols=["eccentricity"],
        categories=categories,
    )

    # ── Figure 1b: match vs params, quasi-circular only (b+d) ────────────────
    if qc_cats:
        print("Generating Figure 1b — match vs source params (quasi-circular b+d)...")
        _scatter_grid(
            df[df["category"].isin(qc_cats)],
            y_col_template="match_{mode}",
            y_label_template=r"Match $\mathcal{F}$",
            modes=modes,
            params=params_qc,
            outpath=os.path.join(outdir, f"fig1b_match_vs_params_qc.{fmt}"),
            y_lim=(0.7, 1.005),
            categories=qc_cats,
        )

    # ── Figure 2a: phase diff vs all params, eccentricity on log x ───────────
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

    # ── Figure 2b: log phase diff vs params, quasi-circular ──────────────────
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

    # ── Figure 3a: match CDFs ─────────────────────────────────────────────────
    print("Generating Figure 3a — match CDFs per mode...")
    _cdf_figure(
        df,
        modes=modes,
        metric_col="match_{mode}",
        x_label=r"Match $\mathcal{F}$",
        outpath=os.path.join(outdir, f"fig3a_match_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(0.5, 1.0),
    )

    # ── Figure 3b: mismatch CDFs (log scale) ─────────────────────────────────
    print("Generating Figure 3b — mismatch CDFs (log)...")
    _cdf_figure(
        df,
        modes=modes,
        metric_col="mismatch_{mode}",
        x_label=r"Mismatch $1 - \mathcal{F}$",
        outpath=os.path.join(outdir, f"fig3b_mismatch_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(1e-4, 1.0),
        log_x=True,
    )

    # ── Figure 3c: phase-diff CDFs ────────────────────────────────────────────
    print("Generating Figure 3c — phase diff CDFs per mode...")
    _cdf_figure(
        df,
        modes=modes,
        metric_col="phase_diff_{mode}",
        x_label=r"$\Delta\Phi/\mathrm{cycle}$ [rad]",
        outpath=os.path.join(outdir, f"fig3c_phasediff_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(0, 2.0),
    )

    # ── Figure 4: (q, χ_eff) heatmaps for (2,2) match ────────────────────────
    print("Generating Figure 4 — (q, χ_eff) match heatmaps (qc)...")
    _heatmap_figure(
        df,
        mode="22",
        metric_col="match_{mode}",
        title=r"Match $\mathcal{F}_{22}$",
        outpath=os.path.join(outdir, f"fig4_match22_heatmap_qc.{fmt}"),
        vmin=0.90,
        vmax=1.0,
        categories=qc_cats or categories,
    )

    # ── Figure 5a: mismatch vs all params (log y), eccentricity log x ─────────
    print("Generating Figure 5a — mismatch vs source params (log scale)...")
    _scatter_grid(
        df,
        y_col_template="mismatch_{mode}",
        y_label_template=r"Mismatch $1{-}\mathcal{F}$",
        modes=modes,
        params=params_all,
        outpath=os.path.join(outdir, f"fig5a_mismatch_vs_params.{fmt}"),
        y_lim=(1e-4, 1.0),
        log_y=True,
        log_x_cols=["eccentricity"],
        categories=categories,
    )

    # ── Figure 5b: mismatch vs qc params (log y) ──────────────────────────────
    if qc_cats:
        print("Generating Figure 5b — mismatch vs source params, qc (log scale)...")
        _scatter_grid(
            df[df["category"].isin(qc_cats)],
            y_col_template="mismatch_{mode}",
            y_label_template=r"Mismatch $1{-}\mathcal{F}$",
            modes=modes,
            params=params_qc,
            outpath=os.path.join(outdir, f"fig5b_mismatch_vs_params_qc.{fmt}"),
            y_lim=(1e-4, 1.0),
            log_y=True,
            categories=qc_cats,
        )

    # ── Figure 6a: SXS calibration-split match CDF ────────────────────────────
    print("Generating Figure 6a — SXS cal-split match CDFs...")
    _cdf_calibration_split_figure(
        df,
        modes=modes,
        metric_col="match_{mode}",
        x_label=r"Match $\mathcal{F}$",
        outpath=os.path.join(outdir, f"fig6a_sxs_cal_match_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(0.5, 1.0),
    )

    # ── Figure 6b: SXS calibration-split mismatch CDF (log) ──────────────────
    print("Generating Figure 6b — SXS cal-split mismatch CDFs (log)...")
    _cdf_calibration_split_figure(
        df,
        modes=modes,
        metric_col="mismatch_{mode}",
        x_label=r"Mismatch $1 - \mathcal{F}$",
        outpath=os.path.join(outdir, f"fig6b_sxs_cal_mismatch_cdf.{fmt}"),
        categories=qc_cats or categories,
        x_lim=(1e-4, 1.0),
        log_x=True,
    )

    # ── Individual simulation figures ─────────────────────────────────────────
    indiv_dir = os.path.join(outdir, "individual_sims")
    print(f"Generating individual simulation figures → {indiv_dir} ...")
    make_individual_sim_figures(df, outdir=indiv_dir, fmt=fmt, dpi=dpi)

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
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    make_all_figures(
        indir=args.indir,
        outdir=args.outdir,
        modes=args.modes,
        categories=args.categories,
        fmt=args.format,
        dpi=args.dpi,
    )
