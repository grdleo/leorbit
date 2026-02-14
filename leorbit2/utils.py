"""Special functions with special purposes. Should not be useful for the average user.
"""

import math
from leorbit2.mathematics import D, Scalar, Quantity, Vector3
from leorbit2.mathematics.tensor import Tensor
from leorbit2.mathematics.dimensions import P3, P2, P1, OO, N1, N2, N3, DimCoords, PowerDim
from leorbit2.mathematics.functions import atan2, cos, sin, square, sqrt
import numpy as np
from numpy.typing import NDArray

from typing import TYPE_CHECKING, Any, TypeAlias, TypeVar, cast, overload

from leorbit2.mathematics.matrix33 import ProductDim
from leorbit2.mathematics.scalar_array import ScalarArray

GravParam: TypeAlias = ProductDim[
    PowerDim[D.Length, P3, P1],
    PowerDim[D.Time, N2, P1]
]

MU_EARTH: Scalar[GravParam] = Scalar.dimensionalize(
    DimCoords(
        length=3, 
        time=-2
    )
).new(
    398_600_441_800_000
)
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = sqrt(MU_EARTH)
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

RADIIE_AA = square(6_378_137 * Quantity.m)
RADIIE_A4 = square(RADIIE_AA)
RADIIE_BB = square(6_356_752 * Quantity.m)
RADIIE_B4 = square(RADIIE_BB)

def mean_motion_to_semi_major_axis_earth(mean_motion: Scalar[D.AngularVelocity]) -> Scalar[D.Length]:
    return Scalar[D.Length].new(
        (MU_EARTH.base_unit_value / mean_motion.base_unit_value**2)**(1/3)
    )

def geocentric_radius_earth(latitude: Scalar[D.Angle]) -> Scalar[D.Length]:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    cc = cos(latitude)
    ss = sin(latitude)
    return sqrt((RADIIE_A4 * cc + RADIIE_B4 * ss) / (RADIIE_AA * cc + RADIIE_BB * ss)).cast(D.Length)

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

def mean2true_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]:
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

def mean2eccentric_anomaly(e: Scalar[D.Dimless], M: Scalar[D.Angle]) -> Scalar[D.Angle]:
    E = M
    for _ in range(5):
        E = M + (e * sin(E)).cast(D.Angle)
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
        sqrt(1 - square(e)) * sin(E), 
        cos(E) - e
    )

@overload
def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: Scalar[D.Angle]) -> Scalar[D.Angle]: ...

@overload
def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: ScalarArray[D.Angle]) -> ScalarArray[D.Angle]: ...

def true2eccentric_anomaly(e: Scalar[D.Dimless], nu: Scalar[D.Angle] | ScalarArray[D.Angle]) -> Scalar[D.Angle] | ScalarArray[D.Angle]:
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


STD_GRAV_PARAM_TERRA_FLOAT = SQRT_MU_EARTH_PINT.m_as("m**1.5/s")

def elements2orthogonal_gcrf(
    υ: Scalar[D.Angle] | ScalarArray[D.Angle], 
    e: Scalar[D.Dimless], 
    a: Scalar[D.Length], 
    Ω: Scalar[D.Angle], 
    ω: Scalar[D.Angle], 
    i: Scalar[D.Angle]
) -> Any:
	"""Returns the position and velocity of a satellite in GCRF coordinates (meters, meters/second)
	All parameters are in radians, except `e` dimensionless and `a` in meters."""

	# position of satellite in orbit plane (with z = 0)
	
	ee = e**2
	one_ee = (1 - ee)
	E = true2eccentric_anomaly(e, υ)
	esinE = e * sin(E)
	r = a * one_ee / (1 + e * cos(υ))
	rd = (STD_GRAV_PARAM_TERRA_FLOAT * sqrt(a) * esinE / r).cast(D.Velocity)
	rυd = rd * one_ee / esinE

	c_raan, s_raan = np.cos(Ω), np.sin(Ω)
	c_i, s_i = np.cos(i), np.sin(i)
	υpω = υ + ω
	c_theta, s_theta = np.cos(υpω), np.sin(υpω)

	def unitvec_gcrf(x: np.float64 | NDArray, y: np.float64 | NDArray) -> tuple[NDArray, NDArray, NDArray]:
		return (
			c_raan * x - s_raan * c_i * y, # X
			s_raan * x + c_raan * c_i * y, # Y
			s_i * y                        # Z
		)

	ur_x, ur_y, ur_z = unitvec_gcrf(c_theta, s_theta)
	ut_x, ut_y, ut_z = unitvec_gcrf(-s_theta, c_theta)

	return PosVelGCRF(
		ur_x * r, 
		ur_y * r, 
		ur_z * r,
		ur_x * rd + ut_x * rυd,
		ur_y * rd + ut_y * rυd,
		ur_z * rd + ut_z * rυd
	)