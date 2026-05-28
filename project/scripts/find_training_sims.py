import os
import numpy as np
import nrcatalogtools as nrcat


def main():
    print("Loading SXS catalog...")
    cat = nrcat.SXSCatalog.load(download=False)

    training_sims = []

    for sim in cat.simulations_list:
        try:
            num = int(sim.split(":")[-1])
            # The NRSur7dq4 model was completed in 2019 using simulations up to SXS:BBH:2082
            if num > 2082:
                continue

            meta = cat.get_metadata(sim)

            # Check eccentricity
            ecc = meta.get("reference_eccentricity", 0.0)
            if isinstance(ecc, str):
                ecc = ecc.replace("<", "").replace("~", "")
                try:
                    ecc = float(ecc)
                except ValueError:
                    ecc = 0.0
            elif ecc is None or np.isnan(ecc):
                ecc = 0.0

            if ecc >= 0.005:
                continue

            # Check mass ratio and spins
            params = cat.get_parameters(sim)
            q = params["mass1"] / params["mass2"]
            chi1 = np.sqrt(
                params["spin1x"] ** 2 + params["spin1y"] ** 2 + params["spin1z"] ** 2
            )
            chi2 = np.sqrt(
                params["spin2x"] ** 2 + params["spin2y"] ** 2 + params["spin2z"] ** 2
            )

            # NRSur7dq4 training limits: q <= 4.0, chi <= 0.8
            if q <= 4.0001 and chi1 <= 0.8001 and chi2 <= 0.8001:
                training_sims.append((num, sim))
        except Exception:
            pass

    # Sort by simulation number
    training_sims.sort()

    # Save to file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.abspath(os.path.join(script_dir, "..", "results"))
    os.makedirs(outdir, exist_ok=True)
    out_file = os.path.join(outdir, "nrsur7dq4_training_simulations.txt")

    with open(out_file, "w") as f:
        f.write(
            "# SXS simulation IDs used to train NRSur7dq4 (total count: {})\n".format(
                len(training_sims)
            )
        )
        for num, sim_id in training_sims:
            f.write(f"{sim_id}\n")

    print(f"Successfully identified {len(training_sims)} training simulations.")
    print(f"Written to: {out_file}")


if __name__ == "__main__":
    main()
