from ast import Or
from genericpath import getmtime
import json
from pathlib import Path
from tempfile import gettempdir
import numpy as np
from pydantic import BaseModel, Field
from requests import HTTPError, get

from leorbit.coordinates import OrbitalElements
from leorbit.m import D, Quantity, Scalar, cube, square
from leorbit.time import Time

from leorbit.utils import convert_quantity_units

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
        deg_to_rad = np.pi / 180
        turn_per_day_to_rad_per_second = 2 * np.pi / (24 * 3600)
        turn_per_day_sqr_to_rad_per_second_sqr = turn_per_day_to_rad_per_second / (24 * 3600)
        turn_per_day_cub_to_rad_per_second_cub = turn_per_day_sqr_to_rad_per_second_sqr / (24 * 3600)
        inv_radiiearth_to_inv_meter = float(1 / Quantity.radii_earth.magnitude("meter"))

        e = Scalar[D.Dimless](self.eccentricity)
        i = Scalar[D.Angle](self.inclination * deg_to_rad)
        Ω = Scalar[D.Angle](self.ra_of_asc_node * deg_to_rad)
        ω = Scalar[D.Angle](self.arg_of_pericenter * deg_to_rad)
        n = Scalar[D.AngularVelocity](self.mean_motion * turn_per_day_to_rad_per_second)
        M = Scalar[D.Angle](self.mean_anomaly * deg_to_rad)
        n_dot = Scalar[D.AngularAcc](
            self.mean_motion_dot * turn_per_day_sqr_to_rad_per_second_sqr
        )
        n_ddot = Scalar[D.AngularJerk](
            self.mean_motion_ddot * turn_per_day_cub_to_rad_per_second_cub
        )
        bstar = Scalar[D.InvLength](self.bstar * inv_radiiearth_to_inv_meter)

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