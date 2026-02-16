from ast import Or
from genericpath import getmtime
import json
from pathlib import Path
from tempfile import gettempdir
from pydantic import BaseModel, Field
from requests import HTTPError, get

from leorbit2.coordinates import OrbitalElements
from leorbit2.m import D, Quantity, Scalar, cube, square
from leorbit2.time import Time

import pint

class CelestrakDataGP(BaseModel):
    """Orbital elements as returned by Celestrak in JSON format"""
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
        rad_per_second = (Quantity.rad / Quantity.second).cast(D.AngularVelocity)
        rad_per_second_squared = (Quantity.rad / square(Quantity.second)).cast(D.AngularAcc)
        rad_per_second_cubed = (Quantity.rad / cube(Quantity.second)).cast(D.AngularJerk)
        inv_meter = (1 / Quantity.meter).cast(D.InvLength)

        e = Scalar[D.Dimless](self.eccentricity)
        i = float(pint.Quantity(self.inclination, "degrees").m_as("radians")) * Quantity.rad
        Ω = float(pint.Quantity(self.ra_of_asc_node, "degrees").m_as("radians")) * Quantity.rad
        ω = float(pint.Quantity(self.arg_of_pericenter, "degrees").m_as("radians")) * Quantity.rad
        n = float(pint.Quantity(self.mean_motion, "turn/day").m_as("radians/second")) * rad_per_second
        M = float(pint.Quantity(self.mean_anomaly, "degrees").m_as("radians")) * Quantity.rad
        n_dot = float(pint.Quantity(self.mean_motion_dot, "turn/day^2").m_as("radians/second^2")) * rad_per_second_squared
        n_ddot = float(pint.Quantity(self.mean_motion_ddot, "turn/day^3").m_as("radians/second^3")) * rad_per_second_cubed
        bstar = float(pint.Quantity(self.bstar, "1/earthRadii").m_as("1/meter")) * inv_meter

        return OrbitalElements(
            epoch=Time.fromisoformat(self.epoch),
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
        last_modified = Time(unixepoch_last_modified)
        if Time.now().delta(last_modified) < (MINIMAL_DURATION_UPDATE_HOURS * Quantity.hour):
            try:
                return CelestrakDataGP(
                    **json.loads(store_path.read_text())
                )
            except: # GP data in storing file may be corrupted
                pass
                    
    print(f"Requesting GP data for object n°{catnr} on `celestrak.com...`") if log else None

    res = get(
        url="https://celestrak.com/NORAD/elements/gp.php", 
        params=dict(CATNR=catnr, FORMAT="JSON")
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