"""Special functions with special purposes. Should not be useful for the average user.
"""

import math
from math import log10, sin, cos, tan, atan2
from pint import Quantity

import numpy as np
from numpy.typing import NDArray
from mathematics.custom import HALF_REV

QtOrArray = TypeVar("FloatOrArray", Quantity, NDArray)

from typing import TYPE_CHECKING, TypeVar
if TYPE_CHECKING:
    from mathematics.vec3 import Vec3

from pint import Quantity
from mathematics.units import UREG

RADIIE_AA = (6_378_137 * UREG.meter)**2
RADIIE_A4 = RADIIE_AA**2
RADIIE_BB = (6_356_752 * UREG.meter)**2
RADIIE_B4 = RADIIE_BB**2

def geocentric_radius_earth(latitude: Quantity | float) -> Quantity:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    if isinstance(latitude, Quantity):
        assert latitude.check("rad")
    cc = cos(latitude)
    ss = sin(latitude)
    rr = (RADIIE_A4 * cc + RADIIE_B4 * ss) / (RADIIE_AA * cc + RADIIE_BB * ss)
    return rr**.5

def apparent_magnitude(sun: Vec3, sat: Vec3, observer: Vec3, std_mag: float) -> float:
    """Returns the apparent magnitude of a satellite from its standard magnitude,
    to a given observer.

    All positions must be given in the same orthogonal frame.
    
    Algorithm from: https://astronomy.stackexchange.com/questions/28744/calculating-the-apparent-magnitude-of-a-satellite"""
    dir_obs = observer - sat
    dir_sun = sun - sat
    dist_sat: Quantity = abs(dir_obs)
    phase = dir_sun.angle(dir_obs)

    phi_term = sin(phase) + (HALF_REV - phase).m_as("rad") * cos(phase)

    return std_mag + 5 * log10(dist_sat.m_as("megameter")) - 2.5 * log10(phi_term)

def humanize_duration(t: Quantity) -> str:
    """Make a `pint.Quantity` duration human-readable.
    
    Example:
    --------
    
    ```python
    humanize_duration(Quantity("723459.23min"))
    >>> '1 year 4 month 15 day 9 hour 39 min 13 s'
    ```
    """
    try:
        assert t.check("s")
    except:
        ValueError(f"{t} is not a `pint.Quantity` or does not have the 'duration' dimension.")

    stages = ["year", "month", "day", "hour", "min", "s"]
    components = {s: 0 for s in stages}

    for s in stages:
        stage_d = Quantity(f"1{s}")
        if t < stage_d:
            continue
        c = int(t.m_as(s))
        components[s] = c
        t -= stage_d * c
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)
    
def angle2dms(angle: Quantity) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle2convert = abs(angle.m_as(UREG.degrees))
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{deg}° {min}′ {sec}″"

def mean2true_anomaly(e: float, M: QtOrArray) -> QtOrArray:
    """ O(e**4)"""
    math_module = np if isinstance(M, np.ndarray) else math
    
    ee = e*e
    eee = e*ee
    return (
        M
        + (2*e - .25*eee) * math_module.sin(M)
        + 1.25*ee * math_module.sin(2*M)
        + (13/12)*eee * math_module.sin(3*M)
    )

def mean2eccentric_anomaly(e: float, M: QtOrArray) -> QtOrArray:
    math_module = np if isinstance(M, np.ndarray) else math

    E = M
    for _ in range(5):
        E = M + e * math_module.sin(E)
    return E

def eccentric2true_anomaly(e: float, E: QtOrArray) -> QtOrArray:
    """Returns the true anomaly from the eccentric anomaly and the excentricity.
    
    Parameters
    ----------
    e : float
        Excentricity of the orbit
    E : float
        Eccentric anomaly in radians
    
    Returns
    -------
    float
        True anomaly in radians"""
    math_module = np if isinstance(E, np.ndarray) else math

    υ = math_module.atan2(
        (1 - e*e)**.5 * math_module.sin(E), 
        math_module.cos(E) - e
    )

    return υ

def true2eccentric_anomaly(e: float, υ: QtOrArray) -> QtOrArray:
    """Returns the eccentric anomaly from the true anomaly and the excentricity.
    Parameters
    ----------
    e : float
        Excentricity of the orbit
    υ : float
        True anomaly in radians
    Returns
    
    -------
    float"""
    math_module = np if isinstance(υ, np.ndarray) else math

    E = math_module.atan2(
        (1 - e*e)**.5 * math_module.sin(υ),
        e + math_module.cos(υ)
    )
    return E