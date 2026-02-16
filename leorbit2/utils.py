"""Special functions with special purposes. Should not be useful for the average user.
"""

import math
from multiprocessing import Value
import numpy as np
from numpy.typing import NDArray

from typing import TYPE_CHECKING, Any, Type, TypeAlias, TypeVar, cast, overload

from leorbit2.m import N2, P1, P3, DimCoords, ProductDim, PowerDim, Scalar, ScalarArray, Tensor_S, Vector3, Quantity, Vector3Array, atan2, cos, ensure_tensor, D, Dim, Tensor_V3, Matrix33, sin, sqrt, square

GravParam: TypeAlias = ProductDim[
    PowerDim[D.Length, P3, P1],
    PowerDim[D.Time, N2, P1]
]

MU_EARTH: Scalar[GravParam] = cast(
    Scalar[GravParam],
    Scalar[DimCoords(length=3, time=-2).to_dimension()](398_600_441_800_000)
)
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = sqrt(MU_EARTH)
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

RADIIE_AA = square(6_378_137 * Quantity.meter)
RADIIE_A4 = square(RADIIE_AA)
RADIIE_BB = square(6_356_752 * Quantity.meter)
RADIIE_B4 = square(RADIIE_BB)

#####################################

def angle2dms(angle: Scalar[D.Angle]) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle_deg = float(angle.get_raw_array("deg"))
    angle2convert = abs(angle_deg)
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{deg}° {min}′ {sec}″"

@overload
def mean_motion_to_semi_major_axis_earth(mean_motion: Scalar[D.AngularVelocity]) -> Scalar[D.Length]: ...

@overload
def mean_motion_to_semi_major_axis_earth(mean_motion: ScalarArray[D.AngularVelocity]) -> ScalarArray[D.Length]: ...

def mean_motion_to_semi_major_axis_earth(mean_motion: Tensor_S[D.AngularVelocity]) -> Tensor_S[D.Length]:
    sma = np.cbrt(MU_EARTH._values / np.square(mean_motion._values))
    return Tensor_S[D.Length](sma)

def geocentric_radius_earth(latitude: Scalar[D.Angle]) -> Scalar[D.Length]:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    cc = cos(latitude)
    ss = sin(latitude)
    quo = (RADIIE_A4 * cc + RADIIE_B4 * ss) / (RADIIE_AA * cc + RADIIE_BB * ss)

    return sqrt(quo).cast(D.Length)

# def apparent_magnitude(sun: Vec3, sat: Vec3, observer: Vec3, std_mag: float) -> float:
#     """Returns the apparent magnitude of a satellite from its standard magnitude,
#     to a given observer.

#     All positions must be given in the same orthogonal frame.
    
#     Algorithm from: https://astronomy.stackexchange.com/questions/28744/calculating-the-apparent-magnitude-of-a-satellite"""
#     dir_obs = observer - sat
#     dir_sun = sun - sat
#     dist_sat: Quantity = abs(dir_obs)
#     phase = dir_sun.angle(dir_obs)

#     phi_term = sin(phase) + (HALF_REV - phase).m_as("rad") * cos(phase)

#     return std_mag + 5 * log10(dist_sat.m_as("megameter")) - 2.5 * log10(phi_term)

def humanize_duration(t: Scalar[D.Time]) -> str:
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
        stage_d = Quantity.get(f"1 {s}").cast(D.Time)
        if t < stage_d:
            continue
        c = int(t.magnitude(s))
        components[s] = c
        d = stage_d * c
        t = t - d
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)

@overload
def mean2true_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]: ...

@overload
def mean2true_anomaly(e: Scalar[D.Dimless], M: ScalarArray[D.Angle]) -> ScalarArray[D.Angle]: ...

@overload
def mean2true_anomaly(e: Tensor_S[D.Dimless], M: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]: ...

def mean2true_anomaly(e: Tensor_S[D.Dimless], M: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]:
    """ O(e**4)"""

    _e = e.cast(D.Angle)
    _ee = (_e*_e).cast(D.Angle)
    _eee = (_e*_ee).cast(D.Angle)
    return (
        M
        + (2 * _e - .25 * _eee) * sin(M)
        + 1.25 * _ee * sin(2 * M)
        + (13 / 12) * _eee * sin(3 * M)
    )

@overload
def mean2eccentric_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]: ...

@overload
def mean2eccentric_anomaly(e: Scalar[D.Dimless], M: ScalarArray[D.Angle]) -> ScalarArray[D.Angle]: ...

@overload
def mean2eccentric_anomaly(e: Tensor_S[D.Dimless], M: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]: ...

def mean2eccentric_anomaly(e: Tensor_S[D.Dimless], M: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]:
    E = M
    for _ in range(5):
        E = M + (e * sin(E)).cast(D.Angle)
    return E

@overload
def eccentric2true_anomaly(e: Scalar[D.Dimless], E: Scalar[D.Angle]) -> Scalar[D.Angle]: ...

@overload
def eccentric2true_anomaly(e: Scalar[D.Dimless], E: ScalarArray[D.Angle]) -> ScalarArray[D.Angle]: ...

@overload
def eccentric2true_anomaly(e: Tensor_S[D.Dimless], E: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]: ...

def eccentric2true_anomaly(e: Tensor_S[D.Dimless], E: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]:
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
        sqrt(1 - square(e)) * sin(E), 
        cos(E) - e
    )

@overload
def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: Scalar[D.Angle]) -> Scalar[D.Angle]: ...

@overload
def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: ScalarArray[D.Angle]) -> ScalarArray[D.Angle]: ...

@overload
def true2eccentric_anomaly(e: Tensor_S[D.Dimless], nu: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]: ...

def true2eccentric_anomaly(e: Tensor_S[D.Dimless], nu: Tensor_S[D.Angle]) -> Tensor_S[D.Angle]:
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
        sqrt(1 - square(e)) * sin(nu),
        cos(nu) + e
    )
    return E

### ELEMENTS TO ORTHOGONAL GCRF ###
### ########################### ###

@overload
def elements2orthogonal_gcrf(
    υ: Scalar[D.Angle], 
    e: Scalar[D.Dimless], 
    a: Scalar[D.Length], 
    Ω: Scalar[D.Angle], 
    ω: Scalar[D.Angle], 
    i: Scalar[D.Angle]
) -> tuple[
    Vector3[D.Length], 
    Vector3[D.Velocity]
]: ...

@overload
def elements2orthogonal_gcrf(
    υ: ScalarArray[D.Angle], 
    e: Scalar[D.Dimless], 
    a: Scalar[D.Length], 
    Ω: Scalar[D.Angle], 
    ω: Scalar[D.Angle], 
    i: Scalar[D.Angle]
) -> tuple[
    Vector3Array[D.Length], 
    Vector3Array[D.Velocity]
]: ...

def elements2orthogonal_gcrf(
    υ: Tensor_S[D.Angle], 
    e: Tensor_S[D.Dimless], 
    a: Tensor_S[D.Length], 
    Ω: Tensor_S[D.Angle], 
    ω: Tensor_S[D.Angle], 
    i: Tensor_S[D.Angle]
) -> tuple[
    Tensor_V3[D.Length], 
    Tensor_V3[D.Velocity]
]:
    """Returns the position and velocity of a satellite in GCRF coordinates (meters, meters/second)
    All parameters are in radians, except `e` dimensionless and `a` in meters."""
    
    # position of satellite in orbit plane (with z = 0)
    sqrt_mu_earth = cast(Tensor_S[GravParam], SQRT_MU_EARTH)
    
    one_ee = (1 - square(e))
    E = true2eccentric_anomaly(e, υ)
    esinE = e * sin(E)
    
    r = a * one_ee / (1 + e * cos(υ))
    rd = (sqrt_mu_earth * sqrt(a) * esinE / r).cast(D.Velocity)
    rυd = rd * one_ee / esinE
    
    c_raan, s_raan = cos(Ω), sin(Ω)
    c_i, s_i = cos(i), sin(i)
    υpω = υ + ω
    c_theta, s_theta = cos(υpω), sin(υpω)
     
    def vector_gcrf_factory(x: Tensor_S[D.Dimless], y: Tensor_S[D.Dimless]) -> Tensor_V3[D.Dimless]:
        return Tensor_V3[D.Dimless].from_components(
            x = c_raan * x - s_raan * c_i * y,
            y = s_raan * x + c_raan * c_i * y,
            z = s_i * y
        )

    ur = vector_gcrf_factory(c_theta, s_theta)
    ut = vector_gcrf_factory(-s_theta, c_theta)

    return (
        ur * r,
        ur * rd + ut * rυd
    )