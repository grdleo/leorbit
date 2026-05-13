from genericpath import getmtime
import json
import os
from pathlib import Path
from tempfile import gettempdir
import numpy as np
from pydantic import BaseModel, Field
from requests import get
import requests
import urllib3

from leorbit.coordinates import OrbitalElements
from leorbit.mathematics import U, scalar
from leorbit.time import Timestamp

AngularAcceleration = U.radian / U.second ** 2
AngularJerk = U.radian / U.second ** 3
AngularVelocity = U.radian / U.second
InvLength = 1 / U.meter

class CelestrakDataGP(BaseModel):
    """Orbital elements as returned by Celestrak in JSON format

    The original data use uppercase field names (aliases) that differ from the
    Python attribute names.  Pydantic will **by default** only populate the
    model using the alias names, which means that if you dump the model using
    the regular field names and try to re‑create an instance with ``**`` the
    constructor will ignore all values.

    To make the object round‑trip nicely we enable ``populate_by_name`` in the
    configuration.  This lets us write e.g.::

        obj = CelestrakDataGP(**data)             # alias keys OK
        d = obj.model_dump()                      # lowercase names
        obj2 = CelestrakDataGP(**d)               # works because of config

    If you prefer to keep using aliases when dumping then pass ``by_alias=True``
    to ``model_dump``/``model_dump_json`` or use ``model_validate`` with the
    appropriate option.
    """
    # allow instantiation from either field names or aliases
    model_config = {"populate_by_name": True}

    epoch: str = Field(alias="EPOCH")
    eccentricity: float = Field(alias="ECCENTRICITY")
    inclination: float = Field(alias="INCLINATION")
    ra_of_asc_node: float = Field(alias="RA_OF_ASC_NODE")
    arg_of_pericenter: float = Field(alias="ARG_OF_PERICENTER")
    mean_motion: float = Field(alias="MEAN_MOTION")
    mean_anomaly: float = Field(alias="MEAN_ANOMALY")
    mean_motion_dot: float = Field(alias="MEAN_MOTION_DOT", default=0)
    mean_motion_ddot: float = Field(alias="MEAN_MOTION_DDOT", default=0)
    bstar: float = Field(alias="BSTAR", default=0)
    name: str = Field(alias="OBJECT_NAME", default="")
    norad_cat_id: int | None = Field(alias="NORAD_CAT_ID", default=None)

    def to_orbital_elements(self) -> "OrbitalElements":
        e = scalar(self.eccentricity).with_units(U.dimensionless)
        i = scalar(np.deg2rad(self.inclination)).with_units(U.radian)
        Ω = scalar(np.deg2rad(self.ra_of_asc_node)).with_units(U.radian)
        ω = scalar(np.deg2rad(self.arg_of_pericenter)).with_units(U.radian)
        n = scalar(self.mean_motion * 2 * np.pi / 86_400).with_units(U.radian / U.second)
        M = scalar(np.deg2rad(self.mean_anomaly)).with_units(U.radian)
        n_dot = scalar(self.mean_motion_dot * 2 * np.pi / 86_400 ** 2).with_units(U.radian / U.second ** 2)
        n_ddot = scalar(self.mean_motion_ddot * 2 * np.pi / 86_400 ** 3).with_units(U.radian / U.second ** 3)
        bstar = scalar(self.bstar / 6_378_137).with_units("1/meter")
        return OrbitalElements(
            epoch=Timestamp.fromisoformat(self.epoch),
            eccentricity=e,
            inclination=i,
            ra_of_asc_node=Ω,
            arg_of_pericenter=ω,
            mean_motion=n,
            mean_anomaly=M,
            mean_motion_dot=n_dot,
            mean_motion_ddot=n_ddot,
            bstar=bstar
        )

TEMPFILE_CELESTRAK_PREFIX = "python_leorbit_celestrak_gpdata_"
MINIMAL_DURATION_UPDATE_HOURS = 1

def get_celestrak_gpdata(catnr: int, log: bool = False) -> CelestrakDataGP:
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
    CelestrakDataGP
        GP data as a CelestrakDataGP instance
    """
    if not (0 < catnr <= 9_999_999_999):
        raise ValueError("NORAD Catalog ID must be a 1 to 9 digit number!")
    
    store_path = Path(gettempdir()) / f"{TEMPFILE_CELESTRAK_PREFIX}{catnr}.json"

    if store_path.exists():
        unixepoch_last_modified = int(getmtime(store_path))
        last_modified = Timestamp(unixepoch_last_modified)
        if Timestamp.now().delta(last_modified) < (MINIMAL_DURATION_UPDATE_HOURS * U.hour):
            try:
                return CelestrakDataGP(
                    **json.loads(store_path.read_text())
                )
            except: # GP data in storing file may be corrupted
                pass
                    
    print(f"Requesting GP data for object n°{catnr} on `celestrak.com...`") if log else None

    try:
        res = get(
            url="https://celestrak.com/NORAD/elements/gp.php", 
            params=dict(CATNR=catnr, FORMAT="JSON")
        )
    except requests.exceptions.SSLError:
        # Last-resort fallback for environments with broken/expired certificate chains.
        if os.environ.get("LEORBIT_ALLOW_INSECURE_SSL", "1") != "1":
            raise

        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        if log:
            print("Warning: SSL certificate validation failed; retrying Celestrak request with verify=False.")

        res = get(
            url="https://celestrak.com/NORAD/elements/gp.php",
            params=dict(CATNR=catnr, FORMAT="JSON"),
            verify=False,
        )

    if not res.ok:
        msg = f"An error occured during TLE fetch on celestrak.org!\nHTTP code {res.status_code}: "
        codes = {
            403: "Your IP has been blocked by celestrak.org... Sorry, not my fault"
        }
        msg += codes.get(res.status_code, res.reason)

        raise RuntimeError(msg)
    elif res.text == "No GP data found":
        raise ValueError(f"ID {catnr} does not correspond to existing GP!")\
    
    try:
        orbital_elements = CelestrakDataGP(
            **res.json()[0]
        )
    except Exception as e:
        raise RuntimeError(f"An error occured during parsing of GP data from celestrak.org! Original error message: {e}")
    
    store_path.write_text(orbital_elements.model_dump_json())
    
    print(f"GP data successfully fetched and stored locally.") if log else None

    return orbital_elements