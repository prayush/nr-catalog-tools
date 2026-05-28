#!/usr/bin/env python
"""Batch NR vs NRSur7dq4 comparison — non-spinning and aligned-spin (categories a–d).

Processes all simulations in the SXS, RIT, and MAYA catalogs that belong to
the non-spinning / aligned-spin sub-space (categories a–d) and lie within the
NRSur7dq4 prior volume (q ≤ 4, |χ₁,₂| ≤ 0.8).  For these systems χ_⊥ = 0 so
the source-frame z-axis is identical between the NR simulation and the
surrogate; the phase maximisation in pycbc.filter.match() absorbs the residual
rotation about z.  No SO(3) optimisation is needed.

Intermediate results are written to CSV after every simulation.  The script is
fully restartable: already-processed (catalog, sim_id) pairs are detected from
the existing CSV and skipped.

Usage
-----
    python batch_aligned_catalogs.py \\
        [--catalogs SXS RIT MAYA] \\
        [--categories a b c d] \\
        [--total-mass 40] \\
        [--workers 4] \\
        [--outdir results] \\
        [--psd aLIGOZeroDetHighPower] \\
        [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import warnings
import multiprocessing as mp

warnings.filterwarnings("ignore")

import numpy as np

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPTS_DIR)
sys.path.insert(0, os.path.join(_SCRIPTS_DIR, "..", ".."))

from surrogate_utils import (
    generate_surrogate_modes,
    check_surrogate_prior,
    SURROGATE_MODES,
    load_nrsur7dq4,
)
from match_utils import compute_mode_match, compute_phase_diff_per_cycle, mode_f_lower
from catalog_utils import load_catalog

# ── Constants ────────────────────────────────────────────────────────────────

DELTA_T = 1.0 / 4096
DISTANCE = 1.0  # Mpc

CATEGORY_FULL_NAME = {
    "a": "non-spinning eccentric",
    "b": "non-spinning non-eccentric",
    "c": "aligned-spin eccentric",
    "d": "aligned-spin non-eccentric",
}

_CLASSIFICATION_DIR = os.path.join(_SCRIPTS_DIR, "..", "..", "catalog_organization")

# CSV columns
_MODE_COLS = [
    f"{metric}_{ell}{abs(em)}"
    for (ell, em) in SURROGATE_MODES
    for metric in ("match", "phase_diff", "n_cycles", "f_lower_mode")
]
CSV_FIELDNAMES = [
    "catalog",
    "sim_id",
    "category",
    "nrsur7dq4_calibration",
    "total_mass",
    "q",
    "mass1",
    "mass2",
    "spin1z",
    "spin2z",
    "chi_eff",
    "eccentricity",
    "f_lower_nr",
    "f_lower_sur",
    "f_lower_match",
    "error",
] + _MODE_COLS


# ── Simulation enumeration ────────────────────────────────────────────────────


def _load_classification(catalog_name: str) -> dict:
    path = os.path.join(
        _CLASSIFICATION_DIR, f"{catalog_name.lower()}_classification.json"
    )
    with open(path) as fh:
        return json.load(fh)


def enumerate_sims(
    catalog_names: list[str], category_keys: list[str]
) -> list[tuple[str, str, str]]:
    """Return (catalog_name, sim_id, category_key) for all target sims."""
    rows = []
    for catname in catalog_names:
        data = _load_classification(catname)
        for key in category_keys:
            full_name = CATEGORY_FULL_NAME[key]
            if full_name not in data:
                continue
            sims = data[full_name]["simulations"]
            for s in sims:
                if isinstance(s, dict):
                    sim_id = s["id"]
                    is_cal = bool(s.get("nrsur7dq4_calibration", False))
                else:
                    sim_id = s
                    is_cal = False
                rows.append((catname, sim_id, key, is_cal))
    return rows


# ── Metadata phase: get params + eccentricity for each sim ───────────────────


def _extract_ecc(catalog_name: str, meta: dict) -> float:
    tag = catalog_name.upper()
    key = "reference_eccentricity" if tag == "SXS" else "eccentricity"
    ecc = meta.get(key, 0.0)
    if isinstance(ecc, str):
        ecc = ecc.replace("<", "").replace("~", "").replace(">", "").strip()
        try:
            ecc = float(ecc)
        except ValueError:
            ecc = 0.0
    if ecc is None or (isinstance(ecc, float) and np.isnan(ecc)):
        ecc = 0.0
    return float(ecc)


def collect_metadata(
    sim_list: list[tuple[str, str, str]],
    total_mass: float,
) -> list[dict]:
    """Fetch params + eccentricity for all sims; filter by surrogate prior.

    Returns a list of dicts ready to pass to workers.  Sims that fail
    parameter extraction or fail the prior cut are silently skipped.
    """
    cats: dict = {}  # cached catalog objects
    results = []

    print(f"\nCollecting metadata for {len(sim_list)} candidate simulations...")
    for i, (catname, sim_id, cat_key, is_cal) in enumerate(sim_list):
        if i > 0 and i % 200 == 0:
            print(f"  {i}/{len(sim_list)} ...", flush=True)

        if catname not in cats:
            cats[catname] = load_catalog(catname)
        cat = cats[catname]

        try:
            params = cat.get_parameters(sim_id, total_mass=total_mass)
            q = params["mass1"] / params["mass2"]
            if np.isnan(q) or np.isnan(params["spin1z"]):
                continue
        except Exception:
            continue

        if not check_surrogate_prior(params):
            continue

        chi1z = float(params["spin1z"])
        chi2z = float(params["spin2z"])
        m1, m2 = float(params["mass1"]), float(params["mass2"])
        chi_eff = (m1 * chi1z + m2 * chi2z) / (m1 + m2)

        try:
            meta = cat.get_metadata(sim_id)
            ecc = _extract_ecc(catname, meta)
        except Exception:
            ecc = float("nan")

        results.append(
            {
                "catalog": catname,
                "sim_id": sim_id,
                "category": cat_key,
                "nrsur7dq4_calibration": is_cal,
                "total_mass": total_mass,
                "q": q,
                "mass1": m1,
                "mass2": m2,
                "spin1z": chi1z,
                "spin2z": chi2z,
                "chi_eff": chi_eff,
                "eccentricity": ecc,
                "f_lower_nr": float(params["f_lower"]),
                "params": params,  # full dict for worker; not written to CSV
            }
        )

    print(f"  {len(results)} sims pass prior cuts.")
    return results


# ── Worker ────────────────────────────────────────────────────────────────────

_worker_cat_cache: dict = {}


def _init_worker():
    """Pre-load the surrogate once per worker process."""
    load_nrsur7dq4()


def _run_one_sim(job: dict, delta_t: float, psd_name: str, distance: float) -> dict:
    """Run the full NR vs NRSur7dq4 comparison for one simulation.

    ``job`` is a metadata dict as produced by ``collect_metadata()``.
    Returns a CSV-ready flat dict.
    """
    catname = job["catalog"]
    sim_id = job["sim_id"]
    params = job["params"]

    # Fields present in the job dict (metadata phase)
    _JOB_FIELDS = (
        "catalog",
        "sim_id",
        "category",
        "nrsur7dq4_calibration",
        "total_mass",
        "q",
        "mass1",
        "mass2",
        "spin1z",
        "spin2z",
        "chi_eff",
        "eccentricity",
        "f_lower_nr",
    )
    row = {k: job[k] for k in _JOB_FIELDS}
    # These are filled in during the worker run
    row["f_lower_sur"] = float("nan")
    row["f_lower_match"] = float("nan")
    row["error"] = ""
    for m in _MODE_COLS:
        row[m] = float("nan")

    # Declare waveform handles so they're accessible after the try block
    wfm = None
    h_sur = None

    try:
        # Load catalog (per-process singleton)
        if catname not in _worker_cat_cache:
            _worker_cat_cache[catname] = load_catalog(catname)
        cat = _worker_cat_cache[catname]

        # Generate surrogate modes
        h_sur, f_lower_sur = generate_surrogate_modes(
            params,
            total_mass=job["total_mass"],
            distance=distance,
            delta_t_seconds=delta_t,
        )
        row["f_lower_sur"] = f_lower_sur
        f_lower_match = max(job["f_lower_nr"], f_lower_sur)
        row["f_lower_match"] = f_lower_match

        # Fetch NR waveform
        wfm = cat.get(sim_id)

        # Per-mode matches
        for (ell, em) in SURROGATE_MODES:
            tag = f"{ell}{abs(em)}"
            try:
                h_nr_c = wfm.get_mode(
                    ell,
                    em,
                    total_mass=job["total_mass"],
                    distance=distance,
                    delta_t_seconds=delta_t,
                )
            except Exception:
                continue

            if (ell, em) not in h_sur:
                continue

            f_low_mode = mode_f_lower(f_lower_match, em)
            row[f"f_lower_mode_{tag}"] = f_low_mode

            mm = compute_mode_match(
                h_nr_c.real(),
                h_sur[(ell, em)].real(),
                f_low_mode,
                psd_name=psd_name,
            )
            dphase, n_cyc = compute_phase_diff_per_cycle(h_nr_c, h_sur[(ell, em)])

            row[f"match_{tag}"] = mm
            row[f"phase_diff_{tag}"] = dphase
            row[f"n_cycles_{tag}"] = n_cyc

    except Exception as exc:
        row["error"] = repr(exc)[:200]

    # Generate 4-panel individual figure while waveforms are in memory
    _indiv_dir = job.get("_indiv_dir")
    if _indiv_dir and wfm is not None and h_sur is not None and not row.get("error"):
        try:
            import matplotlib

            matplotlib.use("Agg")
            from plot_batch_results import plot_individual_sim as _plot_indiv

            _plot_indiv(
                row,
                _indiv_dir,
                wfm=wfm,
                h_sur=h_sur,
                total_mass=job["total_mass"],
                delta_t=delta_t,
            )
        except Exception:
            pass

    return row


def _run_one_sim_top(job: dict) -> dict:
    """Top-level wrapper required for multiprocessing pickling.

    Unpacks the run-time config stamped into the job dict by ``run_batch``.
    """
    return _run_one_sim(
        job,
        delta_t=job.get("_delta_t", DELTA_T),
        psd_name=job.get("_psd_name", "aLIGOZeroDetHighPower"),
        distance=job.get("_distance", DISTANCE),
    )


# ── CSV helpers ───────────────────────────────────────────────────────────────


def _load_done(csv_path: str) -> set[tuple[str, str]]:
    """Return the set of (catalog, sim_id) already written to csv_path."""
    if not os.path.exists(csv_path):
        return set()
    done = set()
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            done.add((row.get("catalog", ""), row.get("sim_id", "")))
    return done


def _ensure_csv(csv_path: str):
    """Write the CSV header if the file doesn't exist yet."""
    if not os.path.exists(csv_path):
        with open(csv_path, "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES).writeheader()


def _append_row(csv_path: str, row: dict):
    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writerow(row)


# ── CSV migration ─────────────────────────────────────────────────────────────


def _migrate_add_calibration(merged_csv: str, catalog_names: list[str]) -> bool:
    """Backfill nrsur7dq4_calibration column into an existing CSV if absent.

    Reads the classification JSONs to look up calibration status for each
    (catalog, sim_id) pair already in the CSV.  Returns True when a migration
    was performed, False when the column was already present.
    """
    import pandas as pd

    if not os.path.exists(merged_csv):
        return False

    # Cheap check: peek at the header only
    df_peek = pd.read_csv(merged_csv, nrows=0)
    if "nrsur7dq4_calibration" in df_peek.columns:
        return False

    # Build (catalog, sim_id) → is_calibration lookup from the JSON files
    cal_lookup: dict[tuple[str, str], bool] = {}
    for catname in catalog_names:
        try:
            data = _load_classification(catname)
        except FileNotFoundError:
            continue
        for cat_data in data.values():
            for s in cat_data.get("simulations", []):
                if isinstance(s, dict):
                    cal_lookup[(catname, s["id"])] = bool(
                        s.get("nrsur7dq4_calibration", False)
                    )
                else:
                    cal_lookup[(catname, s)] = False

    df = pd.read_csv(merged_csv)
    df["nrsur7dq4_calibration"] = df.apply(
        lambda r: cal_lookup.get((r["catalog"], r["sim_id"]), False), axis=1
    )
    # Re-order columns to match CSV_FIELDNAMES (extra cols appended at the end)
    ordered = [c for c in CSV_FIELDNAMES if c in df.columns]
    extra = [c for c in df.columns if c not in ordered]
    df[ordered + extra].to_csv(merged_csv, index=False)
    print(
        f"  Migrated {os.path.basename(merged_csv)}: "
        f"added nrsur7dq4_calibration to {len(df)} rows."
    )
    return True


# ── Main ──────────────────────────────────────────────────────────────────────


def run_batch(
    catalog_names: list[str],
    category_keys: list[str],
    total_mass: float = 40.0,
    workers: int = 4,
    psd_name: str = "aLIGOZeroDetHighPower",
    outdir: str = "results",
    dry_run: bool = False,
):
    os.makedirs(outdir, exist_ok=True)
    merged_csv = os.path.join(outdir, "batch_aligned_all.csv")

    # ── enumerate ────────────────────────────────────────────────────────────
    print("=" * 70)
    print("Step 2 — Batch comparison: non-spinning / aligned-spin (cats a–d)")
    print(f"Catalogs: {catalog_names}   Categories: {category_keys}")
    print(f"Total mass: {total_mass} M☉   Workers: {workers}   PSD: {psd_name}")
    print("=" * 70)

    all_sims = enumerate_sims(catalog_names, category_keys)
    print(f"\nFound {len(all_sims)} simulations in categories a–d (before prior cut).")

    if dry_run:
        from collections import Counter

        cnt = Counter((c, k) for c, _, k, _ic in all_sims)
        for (cat, key), n in sorted(cnt.items()):
            print(f"  {cat} category-{key}: {n} sims")
        return

    # ── migrate existing CSV if nrsur7dq4_calibration column is absent ────────
    _migrate_add_calibration(merged_csv, catalog_names)

    # ── metadata phase (sequential, fast) ────────────────────────────────────
    jobs = collect_metadata(all_sims, total_mass)

    # ── check already done ────────────────────────────────────────────────────
    done = _load_done(merged_csv)
    jobs_todo = [j for j in jobs if (j["catalog"], j["sim_id"]) not in done]
    print(f"\n  {len(done)} already processed, {len(jobs_todo)} remaining.")

    _ensure_csv(merged_csv)

    if not jobs_todo:
        print("  Nothing to do — all sims already processed.")

    # Individual figure output directory (created now so workers can write to it)
    _figs_base = os.path.join(_SCRIPTS_DIR, "..", "figs")
    _indiv_dir = os.path.abspath(os.path.join(_figs_base, "individual_sims"))
    os.makedirs(_indiv_dir, exist_ok=True)

    # Stamp run-time config into each job so the top-level worker can use it
    # (local closures cannot be pickled by multiprocessing).
    for j in jobs_todo:
        j["_delta_t"] = DELTA_T
        j["_psd_name"] = psd_name
        j["_distance"] = DISTANCE
        j["_indiv_dir"] = _indiv_dir

    # ── processing loop ───────────────────────────────────────────────────────
    t0 = time.time()
    n_done = 0
    n_fail = 0

    try:
        from tqdm import tqdm

        _tqdm = tqdm
    except ImportError:
        _tqdm = None

    def _progress_iter(it, total):
        if _tqdm is not None:
            return _tqdm(it, total=total, dynamic_ncols=True)
        return it

    if workers <= 1:
        # Sequential — simpler, easier to debug
        it = _progress_iter(iter(jobs_todo), total=len(jobs_todo))
        _init_worker()
        for job in it:
            row = _run_one_sim_top(job)
            _append_row(merged_csv, row)
            n_done += 1
            if row.get("error"):
                n_fail += 1
            if n_done % 10 == 0:
                elapsed = time.time() - t0
                rate = n_done / elapsed
                remaining = (len(jobs_todo) - n_done) / rate if rate > 0 else 0
                print(
                    f"  [{n_done}/{len(jobs_todo)}] "
                    f"rate={rate:.2f} sim/s  ETA={remaining/60:.1f} min  "
                    f"errors={n_fail}",
                    flush=True,
                )
    else:
        with mp.Pool(
            processes=workers,
            initializer=_init_worker,
        ) as pool:
            it = pool.imap_unordered(_run_one_sim_top, jobs_todo, chunksize=1)
            it = _progress_iter(it, total=len(jobs_todo))
            for row in it:
                _append_row(merged_csv, row)
                n_done += 1
                if row.get("error"):
                    n_fail += 1
                if n_done % 20 == 0:
                    elapsed = time.time() - t0
                    rate = n_done / elapsed
                    remaining = (len(jobs_todo) - n_done) / rate if rate > 0 else 0
                    print(
                        f"  [{n_done}/{len(jobs_todo)}] "
                        f"rate={rate:.2f} sim/s  ETA={remaining/60:.1f} min  "
                        f"errors={n_fail}",
                        flush=True,
                    )

    elapsed = time.time() - t0
    print(
        f"\nDone.  Processed {n_done} sims in {elapsed/60:.1f} min  "
        f"({n_fail} errors).  Results: {merged_csv}"
    )

    # ── write per-catalog CSVs ─────────────────────────────────────────────
    _split_by_catalog(merged_csv, outdir)


def _split_by_catalog(merged_csv: str, outdir: str):
    """Write per-catalog CSV files from the merged CSV."""
    import pandas as pd

    df = pd.read_csv(merged_csv)
    for catname in df["catalog"].unique():
        sub = df[df["catalog"] == catname]
        out_path = os.path.join(outdir, f"batch_aligned_{catname.lower()}.csv")
        sub.to_csv(out_path, index=False)
        print(f"  Wrote {len(sub)} rows → {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def _build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--catalogs",
        nargs="+",
        default=["SXS", "RIT", "MAYA"],
        help="Catalogs to process (default: SXS RIT MAYA)",
    )
    p.add_argument(
        "--categories",
        nargs="+",
        default=["a", "b", "c", "d"],
        help="Category keys to include (default: a b c d)",
    )
    p.add_argument(
        "--total-mass", type=float, default=40.0, help="Total mass in M☉ (default: 40)"
    )
    p.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4)",
    )
    p.add_argument(
        "--psd", default="aLIGOZeroDetHighPower", help="PyCBC analytic PSD name"
    )
    p.add_argument(
        "--outdir",
        default=os.path.abspath(os.path.join(_SCRIPTS_DIR, "..", "results")),
        help="Output directory (default: project/results/)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print sim counts, do not run comparisons",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    run_batch(
        catalog_names=args.catalogs,
        category_keys=args.categories,
        total_mass=args.total_mass,
        workers=args.workers,
        psd_name=args.psd,
        outdir=args.outdir,
        dry_run=args.dry_run,
    )
