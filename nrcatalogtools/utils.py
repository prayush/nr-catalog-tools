import functools
import os
import pathlib
import shutil

import lal
import requests

import sxs

if os.getenv("NR_CATALOG_CACHE"):
    nrcatalog_cache_dir = (
        pathlib.Path(os.getenv("NR_CATALOG_CACHE")).expanduser().resolve()
    )
else:
    nrcatalog_cache_dir = pathlib.Path("~/.cache/").expanduser().resolve()
nr_group_tags = {}
nr_group_tags["SXS"] = "SXS"
nr_group_tags["RIT"] = "RIT"
nr_group_tags["MAYA"] = "MAYA"
nr_group_tags["UNKNOWN"] = "UNKNOWN"

rit_catalog_info = {}
rit_catalog_info["cache_dir"] = nrcatalog_cache_dir / "RIT"
rit_catalog_info["metadata_dir"] = rit_catalog_info["cache_dir"] / "metadata/"
rit_catalog_info["data_dir"] = rit_catalog_info["cache_dir"] / "data/"
rit_catalog_info["url"] = "https://ccrgpages.rit.edu/~RITCatalog/"
rit_catalog_info["metadata_url"] = rit_catalog_info["url"] + "/Metadata/"
rit_catalog_info["data_url"] = rit_catalog_info["url"] + "/Data/"
rit_catalog_info["possible_resolutions"] = [
    100,
    120,
    88,
    84,
    118,
    130,
    140,
    144,
    160,
    200,
]
rit_catalog_info["metadata_file_fmts"] = [
    "RIT:BBH:{:04d}-n{:03d}-id{:d}_Metadata.txt",
    "RIT:eBBH:{:04d}-n{:03d}-ecc_Metadata.txt",
]
rit_catalog_info["waveform_file_fmts"] = [
    "ExtrapStrain_RIT-BBH-{:04d}-n{:03d}.h5",
    "ExtrapStrain_RIT-eBBH-{:04d}-n{:03d}.h5",
]
rit_catalog_info["psi4_file_fmts"] = [
    "ExtrapPsi4_RIT-BBH-{:04d}-n{:03d}-id{:d}.tar.gz",
    "ExtrapPsi4_RIT-eBBH-{:04d}-n{:03d}-ecc.tar.gz",
]
rit_catalog_info["max_id_val"] = 6

maya_catalog_info = {
    "cache_dir": nrcatalog_cache_dir / "MAYA",
    "data_url": "https://cgpstorage.ph.utexas.edu/",
    "metadata_url": "https://cgpstorage.ph.utexas.edu/MAYAmetadata.pkl",
}
maya_catalog_info["data_dir"] = maya_catalog_info["cache_dir"] / "data/"
maya_catalog_info["metadata_dir"] = maya_catalog_info["cache_dir"] / "metadata"


def url_exists(link, num_retries=100):
    """Check if a given URL exists on the web.

    Args:
        link : complete web URL

    Returns:
        bool: True/False whether the URL could be found on WWW.
    """
    requests.packages.urllib3.disable_warnings()
    for n in range(num_retries):
        try:
            response = requests.head(link, verify=False)
            if response.status_code == requests.codes.ok:
                return True
            else:
                return False
        except Exception:
            continue
    return False


def download_file(url, path, progress=False, if_newer=True):
    if url_exists(url):
        try:
            return sxs.utilities.downloads.download_file(
                url, path, progress=progress, if_newer=if_newer
            )
        except Exception:
            requests.packages.urllib3.disable_warnings()
            for n in range(100):
                try:
                    r = requests.get(
                        url, verify=False, stream=True, allow_redirects=True
                    )
                    break
                except Exception:
                    continue
            if r.status_code != 200:
                print(f"An error occurred when trying to access <{url}>.")
                try:
                    print(r.json())
                except Exception:
                    pass
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
            r.raw.read = functools.partial(r.raw.read, decode_content=True)
            path = pathlib.Path(path).expanduser().resolve()
            with path.open("wb") as f:
                shutil.copyfileobj(r.raw, f)
    return path


def call_with_timeout(myfunc, args=(), kwargs={}, timeout=5):
    """
    This function calls user-provided `myfunc` with user-provided
    `args` and `kwargs` in a separate multiprocessing.Process.
    if the function evaluation takes more than `timeout` seconds, the
    `Process` is terminated and error raised. If it evalutes within
    `timeout` seconds, the results are fetched from the `Queue` and
    returned.
    """

    from multiprocessing import Process, Queue

    def funcwrapper(p, *args, **kwargs):
        """
        This thin wrapper calls the user-provided function, and puts
        its result into the multiprocessing `Queue`  so that it can be
        obtained via `Queue().get()`.
        """
        res = myfunc(*args, **kwargs)
        p.put(res)

    queue = Queue()
    task = Process(target=funcwrapper, args=(queue, *args))
    task.start()
    task.join(timeout=timeout)
    task.terminate()
    try:
        result = queue.get(timeout=0)
        return result
    except Exception:
        raise Exception("Timeout")


def time_to_physical(M):
    """Factor to convert time from dimensionless units to SI units

    parameters
    ----------
    M: mass of system in the units of solar mass

    Returns
    -------
    converting factor
    """

    return M * lal.MTSUN_SI


def amp_to_physical(M, D):
    """Factor to rescale strain to mass M and distance D convert from
    dimensionless units to SI units

    parameters
    ----------
    M: mass of the system in units of solar mass
    D: Luminosity distance in units of megaparsecs

    Returns
    -------
    Scaling factor
    """

    return lal.G_SI * M * lal.MSUN_SI / (lal.C_SI**2 * D * 1e6 * lal.PC_SI)
