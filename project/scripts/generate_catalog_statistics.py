import os
import nrcatalogtools as nrcat
import numpy as np


def parse_eccentricity(e_val):
    if e_val is None or e_val == "" or e_val == "None" or str(e_val).lower() == "nan":
        return 0.0
    if isinstance(e_val, str):
        # Handle cases like '<1e-4'
        e_val = e_val.replace("<", "").replace("~", "")
    try:
        return float(e_val)
    except Exception:
        return 0.0


def classify_sim(s1x, s1y, s1z, s2x, s2y, s2z, ecc):
    EPS = 1e-4
    mag_s1 = np.sqrt(s1x**2 + s1y**2 + s1z**2)
    mag_s2 = np.sqrt(s2x**2 + s2y**2 + s2z**2)

    is_non_spinning = (mag_s1 < EPS) and (mag_s2 < EPS)
    is_aligned = (
        not is_non_spinning
        and (abs(s1x) < EPS)
        and (abs(s1y) < EPS)
        and (abs(s2x) < EPS)
        and (abs(s2y) < EPS)
    )
    is_ecc = ecc >= 0.005

    if is_non_spinning:
        spin_cat = "Non-spinning"
    elif is_aligned:
        spin_cat = "Aligned-spinning"
    else:
        spin_cat = "Precessing-spinning"

    ecc_cat = "Eccentric" if is_ecc else "Non-eccentric"
    return spin_cat, ecc_cat


def analyze_catalog(name, catalog, get_spins_ecc):
    stats = {
        "Non-spinning": {"Non-eccentric": 0, "Eccentric": 0},
        "Aligned-spinning": {"Non-eccentric": 0, "Eccentric": 0},
        "Precessing-spinning": {"Non-eccentric": 0, "Eccentric": 0},
    }

    for sim in catalog.simulations_list:
        try:
            meta = catalog.get_metadata(sim)
            s1x, s1y, s1z, s2x, s2y, s2z, ecc = get_spins_ecc(meta)
            ecc_val = parse_eccentricity(ecc)
            spin_cat, ecc_cat = classify_sim(s1x, s1y, s1z, s2x, s2y, s2z, ecc_val)
            stats[spin_cat][ecc_cat] += 1
        except Exception:
            pass

    return stats


# Extractors
def sxs_extract(meta):
    s1 = meta.get("reference_dimensionless_spin1", [0, 0, 0])
    s2 = meta.get("reference_dimensionless_spin2", [0, 0, 0])
    ecc = meta.get("reference_eccentricity", 0.0)
    if not isinstance(s1, list) or len(s1) != 3:
        s1 = [0, 0, 0]
    if not isinstance(s2, list) or len(s2) != 3:
        s2 = [0, 0, 0]
    s1 = [0.0 if np.isnan(x) else x for x in s1]
    s2 = [0.0 if np.isnan(x) else x for x in s2]
    return s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], ecc


def rit_extract(meta):
    s1x = meta.get("relaxed-chi1x", 0.0)
    s1y = meta.get("relaxed-chi1y", 0.0)
    s1z = meta.get("relaxed-chi1z", 0.0)
    s2x = meta.get("relaxed-chi2x", 0.0)
    s2y = meta.get("relaxed-chi2y", 0.0)
    s2z = meta.get("relaxed-chi2z", 0.0)
    ecc = meta.get("eccentricity", 0.0)
    return s1x, s1y, s1z, s2x, s2y, s2z, ecc


def maya_extract(meta):
    s1x = meta.get("a1x", 0.0)
    s1y = meta.get("a1y", 0.0)
    s1z = meta.get("a1z", 0.0)
    s2x = meta.get("a2x", 0.0)
    s2y = meta.get("a2y", 0.0)
    s2z = meta.get("a2z", 0.0)
    ecc = meta.get("eccentricity", 0.0)
    return s1x, s1y, s1z, s2x, s2y, s2z, ecc


def main():
    print("Loading catalogs...")
    sxs = nrcat.SXSCatalog.load(download=False)
    rit = nrcat.RITCatalog.load(download=False)
    maya = nrcat.MayaCatalog.load(download=False)

    print("Analyzing SXS...")
    sxs_stats = analyze_catalog("SXS", sxs, sxs_extract)
    print("Analyzing RIT...")
    rit_stats = analyze_catalog("RIT", rit, rit_extract)
    print("Analyzing MAYA...")
    maya_stats = analyze_catalog("MAYA", maya, maya_extract)

    catalogs = {"SXS": sxs_stats, "RIT": rit_stats, "GT (MAYA)": maya_stats}

    # Prepare markdown table
    md_content = "# NR Catalog Statistics\n\n"
    md_content += (
        "The following table breaks down the simulations in the SXS, RIT, and "
        "GT (MAYA) catalogs by their spin configuration and initial eccentricity.\n\n"
    )
    md_content += "**Definitions:**\n"
    md_content += (
        "- **Non-spinning**: Magnitude of both dimensionless spins is $< 10^{-4}$.\n"
    )
    md_content += (
        "- **Aligned-spinning**: Not non-spinning, and all in-plane spin components "
        "($x$ and $y$) are $< 10^{-4}$.\n"
    )
    md_content += "- **Precessing-spinning**: All other spin configurations.\n"
    md_content += "- **Non-eccentric**: Initial eccentricity $< 0.005$.\n"
    md_content += "- **Eccentric**: Initial eccentricity $\\ge 0.005$.\n\n"

    md_content += "| Spin Category | Eccentricity | SXS | RIT | GT (MAYA) |\n"
    md_content += "|---------------|--------------|-----:|-----:|-----------:|\n"

    for spin_cat in ["Non-spinning", "Aligned-spinning", "Precessing-spinning"]:
        for ecc_cat in ["Non-eccentric", "Eccentric"]:
            sxs_val = catalogs["SXS"][spin_cat][ecc_cat]
            rit_val = catalogs["RIT"][spin_cat][ecc_cat]
            maya_val = catalogs["GT (MAYA)"][spin_cat][ecc_cat]
            md_content += f"| {spin_cat} | {ecc_cat} | {sxs_val:,} | {rit_val:,} | {maya_val:,} |\n"

    # Write to results directory
    # Script is in project/scripts/, results go to project/scripts/results/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results")

    os.makedirs(results_dir, exist_ok=True)
    out_file = os.path.join(results_dir, "catalog_statistics.md")

    with open(out_file, "w") as f:
        f.write(md_content)

    print(f"\nStatistics written successfully to: {out_file}")


if __name__ == "__main__":
    main()
