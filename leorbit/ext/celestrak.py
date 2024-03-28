from requests import get as get_url
from requests import HTTPError
from typing import Callable
from tempfile import gettempdir
from json import loads, dumps
from json.decoder import JSONDecodeError
from pathlib import Path
from leorbit.math import Q_
from leorbit.time import Time
from os.path import getmtime

TEMPFILE_CELESTRAK_PREFIX = "python_earthorbit_celestrak_gpdata_"
MINIMAL_DURATION_UPDATE = Q_("1h")

def get_celestrak_gpdata_json(catnr: int, log: bool = False) -> dict[str, str | float | int]:
    """
    Retrieves GP data for object with given CATNR by fetching data on `celestrak.com.`
    If this GP was fetched recently, uses locally stored GP data instead of making a request.
    (https://celestrak.org/NORAD/documentation/gp-data-formats.php)

    Parameters
    ----------
    catnr : int
        NORAD catalog number of satellite

    Returns
    -------
    dict[str, str | float | int]
        GP data as a dictionary
    """
    if not (0 < catnr < int(1e10)):
        raise ValueError("NORAD Catalog ID must be a 1 to 9 digit number!")
    
    store_path = Path(gettempdir()) / f"{TEMPFILE_CELESTRAK_PREFIX}{catnr}.json"

    if store_path.exists():
        unixepoch_last_modified = int(getmtime(store_path))
        last_modified = Time(unixepoch_last_modified)
        if Time.now().delta(last_modified) < MINIMAL_DURATION_UPDATE:
            try:
                json_str: str = None
                with open(store_path, "r") as f:
                    json_str = f.read()
                gp = loads(json_str)
                print(f"Recent GP data for object n°{catnr} found in system.") if log else None
                return gp
            except: # GP data in storing file may be corrupted
                pass
                    
    print(f"Requesting GP data for object n°{catnr} on `celestrak.com...`") if log else None

    res = ""
    try:
        res = get_url(f"https://celestrak.com/NORAD/elements/gp.php?CATNR={catnr}&FORMAT=JSON").text
    except HTTPError as ex:
        msg = f"An error occured during TLE fetch on celestrak.org!\nHTTP code {ex.code}: "
        codes = {
            403: "Your IP has been blocked by celestrak.org... Sorry, not my fault"
        }
        msg += codes.get(ex.code, ex.msg)

        raise RuntimeError(msg)

    if res == "No GP data found":
        raise ValueError(f"ID {catnr} does not correspond to existing GP!")
    if any(error_kw in res for error_kw in ("403", "Forbidden", "denied", "Error")):
        raise HTTPError(
            "You have been temporarily blocked by Celestrak API!\n"
            "As you can see, this very robust and wonderful API cannot handle your 3 requests per minute,"
            "so you have been blocked for two hours.\n"
            "This should not happen since this library locally stores fetched data and reuses it "
            "if not too old, but hey, things happen!"
        )
    
    gp: dict = loads(res)
    if isinstance(gp, list):
        gp = gp[0]
    
    with open(store_path, "w") as f:
        json_str = dumps(gp)
        f.write(json_str)
    
    print(f"GP data successfully fetched and stored locally.") if log else None

    return gp