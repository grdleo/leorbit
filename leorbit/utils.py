"""Special functions with special purposes. Should not be useful for the average user.
"""

import numpy as np

from typing import Iterable, NamedTuple, TypeAlias, TypeVar, Union, cast, overload, TYPE_CHECKING

if TYPE_CHECKING:
    import pint

from leorbit.m import N2, P1, P3, DimCoords, Number, ProductDim, PowerDim, Scalar, ScalarArray, Tensor_S, Vector3, Quantity, Vector3Array, atan2, cos, abs, cube, ensure_tensor, D, Dim, Tensor_V3, Matrix33, normalize_angle, sin, sqrt, square

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
    
    angle_deg = angle.magnitude("deg")
    angle2convert = np.abs(angle_deg)
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{int(deg): 04}° {int(min):02}′ {int(sec):02}″"

@overload
def mean_motion_to_semi_major_axis_earth(mean_motion: Scalar[D.AngularVelocity]) -> Scalar[D.Length]: ...

@overload
def mean_motion_to_semi_major_axis_earth(mean_motion: ScalarArray[D.AngularVelocity]) -> ScalarArray[D.Length]: ...

def mean_motion_to_semi_major_axis_earth(mean_motion: Tensor_S[D.AngularVelocity]) -> Tensor_S[D.Length]:
    sma = np.cbrt(MU_EARTH._values / np.square(mean_motion._values))
    return Tensor_S[D.Length](sma)

@overload
def semi_major_axis_earth_to_mean_motion(sma: Scalar[D.Length]) -> Scalar[D.AngularVelocity]: ...

@overload
def semi_major_axis_earth_to_mean_motion(sma: ScalarArray[D.Length]) -> ScalarArray[D.AngularVelocity]: ...

def semi_major_axis_earth_to_mean_motion(sma: Tensor_S[D.Length]) -> Tensor_S[D.AngularVelocity]:
    mm = np.sqrt(MU_EARTH._values / sma._values ** 3)
    return Tensor_S[D.AngularVelocity](mm)

@overload
def geocentric_radius_earth(latitude: Scalar[D.Angle]) -> Scalar[D.Length]: ...

@overload
def geocentric_radius_earth(latitude: ScalarArray[D.Angle]) -> ScalarArray[D.Length]: ...

@overload
def geocentric_radius_earth(latitude: Tensor_S[D.Angle]) -> Tensor_S[D.Length]: ...

def geocentric_radius_earth(latitude: Tensor_S[D.Angle]) -> Tensor_S[D.Length]:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    cc = cos(latitude)
    ss = sin(latitude)
    cc2 = square(cc)
    ss2 = square(ss)
    quo = (cc2 * RADIIE_A4 + ss2 * RADIIE_B4) / (cc2 * RADIIE_AA + ss2 * RADIIE_BB)

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

ScalarType = TypeVar("ScalarType", bound=Scalar | ScalarArray)

class OrbitalElementsTuple(NamedTuple):
    eccentricity: Scalar[D.Dimless]
    inclination: Scalar[D.Angle]
    ra_of_asc_node: Scalar[D.Angle]
    arg_of_pericenter: Scalar[D.Angle]
    mean_motion: Scalar[D.AngularVelocity]
    mean_anomaly: Scalar[D.Angle]
    

def gcrf_state_vectors2elements(pos: Vector3[D.Length], vel: Vector3[D.Velocity]) -> OrbitalElementsTuple:
    north = Vector3[D.Dimless].from_components(0.0, 0.0, 1.0)

    kinetic = pos.cross(vel)
    kinetic_sq = kinetic.length_squared
    pos_dir = pos.normalized()
    ecc_vec = (vel.cross(kinetic) / MU_EARTH).cast(D.Dimless) - pos_dir
    descending_node = kinetic.cross(north).normalized() # descending node line
    asc = -descending_node

    # create a 2D frame on the ellipsis, x along ascending node line
    xaxis_asc = asc.normalized()
    yaxis_asc = kinetic.cross(asc).normalized()

    xe_asc = xaxis_asc.dot(ecc_vec)
    ye_asc = yaxis_asc.dot(ecc_vec)
    argp = atan2(ye_asc, xe_asc)

    xp_asc = xaxis_asc.dot(pos)
    yp_asc = yaxis_asc.dot(pos)
    nu = normalize_angle(atan2(yp_asc, xp_asc) - argp)

    e = ecc_vec.length
    ee = e * e
    eee = ee * e
    eeee = eee * e
    i = north.angle(kinetic.normalized())
    raan = normalize_angle(atan2(asc.y, asc.x))
    a = (kinetic_sq / (MU_EARTH * (1 - ee))).cast(D.Length)
    aaa = cube(a)
    n = sqrt(MU_EARTH / aaa).cast(D.AngularVelocity)
    M = (
        nu
        - 2 * e * sin(nu).cast(D.Angle)
        + (ee * .75 + eeee * .125) * sin(2 * nu).cast(D.Angle)
        - eee * sin(3 * nu).cast(D.Angle) / 3
        + eeee * sin(4 * nu).cast(D.Angle) * .15625
    )

    return OrbitalElementsTuple(
        e,
        i,
        raan,
        argp,
        n,
        M
    )

class TupleGPS(NamedTuple):
    latitude: Tensor_S[D.Angle]
    longitude: Tensor_S[D.Angle]
    altitude: Tensor_S[D.Length]

def itrf2gps(itrf_pos: Tensor_V3[D.Length]) -> TupleGPS:
    """Convert a position in ITRS coordinates to GPS coordinates, at a given time."""
    lon = itrf_pos.theta
    lat = itrf_pos.delta
    alt = itrf_pos.length - geocentric_radius_earth(lat)

    return TupleGPS(
        latitude=lat,
        longitude=lon,
        altitude=alt
    )

def convert_quantity_units(quantity: Number, units_from: Union[str, "pint.Unit"], units_to: Union[str, "pint.Unit"]) -> Number:
    """Converts a number `quantity` that is expressed in `units_from`, to `units_to` (uses Pint)"""
    import pint

    return cast(Number, pint.Quantity(quantity, units_from).m_as(units_to))

def truth_array_to_indices_intervals(
    truth_array: Iterable[bool] | np.ndarray,
    min_size_intervals: int = 1
) -> list[tuple[int, int]]:
    """Return contiguous ``True`` runs as inclusive index intervals.

    Parameters
    ----------
    truth_array:
        Boolean sequence (Python iterable or NumPy array). If not already a
        NumPy array, it is converted with ``np.asarray(..., dtype=bool)``.
        The input is flattened to 1-D.
    min_size_intervals:
        Minimum run length to keep. Runs shorter than this value are filtered
        out. Must be >= 1.

    Returns
    -------
    list[tuple[int, int]]
        Inclusive ``(start, stop)`` index pairs for each kept ``True`` run.

    Notes
    -----
    The implementation is fully vectorized (padding + ``np.diff``), avoiding
    Python loops over elements.

    Example
    -------
    >>> t = [False, True, True, True, False, False, True, True, False]
    >>> truth_array_to_indices_intervals(t)
    [(1, 3), (6, 7)]
    """
    if min_size_intervals < 1:
        raise ValueError("min_size_intervals must be >= 1")

    if isinstance(truth_array, np.ndarray):
        arr = np.asarray(truth_array, dtype=bool).reshape(-1)
    else:
        arr = np.asarray(list(truth_array), dtype=bool).reshape(-1)
    if arr.size == 0:
        return []

    padded = np.pad(arr.astype(np.int8), (1, 1), mode="constant")
    edges = np.diff(padded)

    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1) - 1

    lengths = stops - starts + 1
    keep = lengths >= int(min_size_intervals)

    return list(zip(starts[keep].tolist(), stops[keep].tolist()))