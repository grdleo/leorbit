"""Special functions with special purposes. Should not be useful for the average user.
"""

import math
from leorbit2.mathematics import D, Scalar, Quantity, Vector3
from leorbit2.mathematics.dimensions import DimEls
from leorbit2.mathematics.functions import atan2
import numpy as np
from numpy.typing import NDArray

from typing import TYPE_CHECKING, TypeVar, cast

from leorbit2.mathematics.scalar import scalar_class_factory

MU_EARTH = scalar_class_factory(
    DimEls(
        length=3, 
        time=-2
    )
).new(
    398_600_441_800_000
)
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = MU_EARTH.sqrt()
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

RADIIE_AA = (6_378_137 * Quantity.m).sqr()
RADIIE_A4 = RADIIE_AA.sqr()
RADIIE_BB = (6_356_752 * Quantity.m).sqr()
RADIIE_B4 = RADIIE_BB.sqr()

def mean_motion_to_semi_major_axis_earth(mean_motion: Scalar[D.AngularVelocity]) -> Scalar[D.Length]:
    return Scalar[D.Length].new(
        (MU_EARTH.base_unit_value / mean_motion.base_unit_value**2)**(1/3)
    )

def geocentric_radius_earth(latitude: Scalar[D.Angle]) -> Scalar[D.Length]:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    cc = latitude.cos()
    ss = latitude.sin()
    return ((RADIIE_A4 * cc + RADIIE_B4 * ss) / (RADIIE_AA * cc + RADIIE_BB * ss)).sqrt().cast(D.Length)

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

def humanize_duration(t: Scalar[Dim.time]) -> str:
    """Make a duration human-readable.
    
    Example:
    --------
    
    ```python
    humanize_duration(Quantity("723459.23min"))
    >>> '1 year 4 month 15 day 9 hour 39 min 13 s'
    ```
    """
    stages = ["year", "month", "day", "hour", "min", "s"]
    components = {s: 0 for s in stages}

    for s in stages:
        stage_d: Scalar[Dim.time] = quantity(f"1 {s}")
        if t < stage_d:
            continue
        c = int(t.magnitude(s))
        components[s] = c
        d = stage_d * c
        t = t - d
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)

def mean2true_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """ O(e**4)"""

    _e = e.cast(D.Angle)
    _ee = (_e*_e).cast(D.Angle)
    _eee = (_e*_ee).cast(D.Angle)
    return (
        M
        + (2 * _e - .25 * _eee) * M.sin()
        + 1.25 * _ee * (2 * M).sin()
        + (13 / 12) * _eee * (3 * M).sin()
    )

def mean2eccentric_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]:
    math_module = np if isinstance(M, np.ndarray) else math

    E = M
    for _ in range(5):
        E = M + (e * E.sin()).cast(D.Angle)
    return E

def eccentric2true_anomaly(e: Scalar[D.Dimless], E: Scalar[D.Angle]) -> Scalar[D.Angle]:
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

    return atan2(
        (1 - e*e)**.5 * E.sin(), 
        E.cos() - e
    )


def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """Returns the eccentric anomaly from the true anomaly and the excentricity.
    Parameters
    ----------
    e : float
        Excentricity of the orbit
    nu : float
        True anomaly in radians
    Returns
    
    -------
    float"""

    E = atan2(
        (1 - e*e)**.5 * nu.sin(),
        e + nu.cos()
    )
    return E