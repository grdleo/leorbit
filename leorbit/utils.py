"""Special functions with special purposes. Should not be useful for the average user.
"""

from fractions import Fraction
import numpy as np
import numpy.typing as npt

from typing import Annotated, Any, Iterable, NamedTuple, TypeAlias, TypeVar, Union, cast, overload, TYPE_CHECKING

from leorbit.transforms import TransformVector3Affine

if TYPE_CHECKING:
    from leorbit.frames import EarthLocalFrame
    import pint

from leorbit.mathematics import AngularVelocity, Dimless, Length, RealNumber, Angle, Tensor, Quantity, TensorBound, TensorKind, Time, atan2, cos, normalize_angle, sin, tan, vector3, tensor_check

MU_EARTH = 398_600_441_800_000 * (Quantity.meter ** 3 / Quantity.second ** 2)
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = MU_EARTH ** .5
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

RADIIE_AA = (6_378_137 * Quantity.meter) ** 2
RADIIE_A4 = RADIIE_AA ** 2
RADIIE_BB = (6_356_752 * Quantity.meter) ** 2
RADIIE_B4 = RADIIE_BB ** 2

TWOPI = 2 * np.pi
TWELF_PI = np.pi / 12

#####################################

def angle2dms(angle: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle_deg = angle.scalar.value("deg")
    angle2convert = np.abs(angle_deg)
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{int(deg): 04}° {int(min):02}′ {int(sec):02}″"

def unixepoch_to_j2000(unixepoch: npt.NDArray | RealNumber) -> npt.NDArray:
    """Representation of this `Time` object as "Julian year (J2000)", aka 
    the number of days since 2000/01/01T12:00:00."""
    return np.asarray(unixepoch) / 86_400 - 10_957.5

def j2000_to_stl0(j2000: npt.NDArray | RealNumber) -> npt.NDArray | RealNumber:
    """The 
    [Sideral Time](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
    (angle) of Latitude 0 at this `Time`, in radians.
    """
    return ((np.longdouble(18.697374558) + np.longdouble(24.06570982441908) * np.asarray(j2000)) * TWELF_PI) % TWOPI

def jd(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(dimension=Time, kind=TensorKind.SCALAR)]:
    """Representation of this `Timestamp` object as "Julian day (JD)", aka 
    the number of days since -4712/01/01."""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    days = (unixepoch / 86_400 + 2_440_587.5)
    return Tensor(days) * Quantity.day

def j2000(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(dimension=Time, kind=TensorKind.SCALAR)]:
    """Representation of this `Timestamp` object as "Julian year (J2000)", aka 
    the number of days since 2000/01/01T12:00:00."""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    return Tensor(unixepoch_to_j2000(unixepoch)) * Quantity.day

def from_mil(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(dimension=Time, kind=TensorKind.SCALAR)]:
    """Representation as a fraction of days since 1 january 2000 00:00.

    Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    days = (unixepoch / 86_400 - 10_957.5) - .5
    return Tensor(days) * Quantity.day

def stl0(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
    """Sidereal time (angle) of latitude 0 at this unixepoch."""
    j2k = j2000(unixepoch)
    sidereal = j2000_to_stl0(j2k.raw_data_array("day"))
    return Tensor(np.asarray(sidereal, dtype=np.float64)) * Quantity.radian

@tensor_check
def mean_motion_to_semi_major_axis_earth(
    mean_motion: Annotated[Tensor, TensorBound(dimension=Angle / Time, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
    return (MU_EARTH / mean_motion ** 2) ** Fraction(1, 3)

@tensor_check
def semi_major_axis_earth_to_mean_motion(
        sma: Annotated[Tensor,  TensorBound(dimension=Length, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Angle / Time, kind=TensorKind.SCALAR)]:
    return (MU_EARTH / sma ** 3) ** .5

@tensor_check
def geocentric_radius_earth(
    latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    cc = cos(latitude)
    ss = sin(latitude)
    cc2 = cc ** 2
    ss2 = ss ** 2
    quo = (cc2 * RADIIE_A4 + ss2 * RADIIE_B4) / (cc2 * RADIIE_AA + ss2 * RADIIE_BB)

    return quo ** .5

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

def humanize_duration(
    t: Annotated[Tensor, TensorBound(dimension=Time, kind=TensorKind.SCALAR, size=1)]
) -> str:
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
        stage_d = Quantity.get(s)
        if t < stage_d:
            continue
        c = int(t.scalar.value(s))
        components[s] = c
        d = stage_d * c
        t = t - d
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)

@tensor_check
def mean2true_anomaly(
    e: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)],
    M: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
    """ O(e**4)"""

    ee = e ** 2
    eee = e ** 3

    return (
        M
        + (2 * e - .25 * eee) * sin(M)
        + 1.25 * ee * sin(2 * M)
        + (13 / 12) * eee * sin(3 * M)
    )

@tensor_check
def mean2eccentric_anomaly(
    e: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)],
    M: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
    E = M
    for _ in range(5):
        E = M + e * sin(E)
    return E

@tensor_check
def eccentric2true_anomaly(
    e: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)], 
    E: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
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
        sin(E) * (1 - e ** 2) ** .5, 
        cos(E) - e
    )

@tensor_check
def true2eccentric_anomaly(
    e: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)], 
    nu: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
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
        sin(nu) * (1 - e ** 2) ** .5,
        cos(nu) + e
    )
    return E

### ELEMENTS TO ORTHOGONAL GCRF ###
### ########################### ###

@tensor_check
def elements2orthogonal_gcrf(
    υ: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)], 
    e: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)], 
    a: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)], 
    Ω: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)], 
    ω: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
    i: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
) -> tuple[
    Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)], 
    Annotated[Tensor, TensorBound(dimension=Length / Time, kind=TensorKind.VECTOR3)]
]:
    """Returns the position and velocity of a satellite in GCRF coordinates (meters, meters/second)
    All parameters are in radians, except `e` dimensionless and `a` in meters."""
    
    # position of satellite in orbit plane (with z = 0)
    one_ee = 1 - e ** 2
    E = true2eccentric_anomaly(e, υ)
    esinE = e * sin(E)
    
    r = a * one_ee / (1 + e * cos(υ))
    rd = (SQRT_MU_EARTH * a ** .5 * esinE / r).secure(dimension=Length / Time)
    rυd = rd * one_ee / esinE
    
    c_raan, s_raan = cos(Ω), sin(Ω)
    c_i, s_i = cos(i), sin(i)
    υpω = υ + ω
    c_theta, s_theta = cos(υpω), sin(υpω)
     
    def vector_gcrf_factory(
        x: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)],
        y: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.VECTOR3)]:
        vx = vector3(1, 0, 0)
        vy = vector3(0, 1, 0)
        vz = vector3(0, 0, 1)

        return (
            vx * (c_raan * x - s_raan * c_i * y)
            + vy * (s_raan * x + c_raan * c_i * y)
            + vz * (s_i * y)
        )

    ur = vector_gcrf_factory(c_theta, s_theta)
    ut = vector_gcrf_factory(-s_theta, c_theta)

    return (
        ur * r,
        ur * rd + ut * rυd
    )

class OrbitalElementsTuple(NamedTuple):
    eccentricity: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]
    inclination: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    ra_of_asc_node: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    arg_of_pericenter: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    mean_motion: Annotated[Tensor, TensorBound(dimension=Angle / Time, kind=TensorKind.SCALAR)]
    mean_anomaly: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]

@tensor_check
def gcrf_state_vectors2elements(
    pos: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)],
    vel: Annotated[Tensor, TensorBound(dimension=Length / Time, kind=TensorKind.VECTOR3)]
) -> OrbitalElementsTuple:
    north = vector3(0.0, 0.0, 1.0)

    kinetic = pos.vector3.cross(vel.vector3)
    kinetic_sq = kinetic.vector3.length_squared
    pos_dir = pos.vector3.normalized()
    ecc_vec = (vel.vector3.cross(kinetic) / MU_EARTH) - pos_dir
    descending_node = kinetic.vector3.cross(north.vector3).vector3.normalized() # descending node line
    asc = -descending_node

    # create a 2D frame on the ellipsis, x along ascending node line
    xaxis_asc = asc.vector3.normalized()
    yaxis_asc = kinetic.vector3.cross(asc.vector3).vector3.normalized()

    xe_asc = xaxis_asc.vector3.dot(ecc_vec.vector3)
    ye_asc = yaxis_asc.vector3.dot(ecc_vec.vector3)
    argp = atan2(ye_asc, xe_asc)

    xp_asc = xaxis_asc.vector3.dot(pos.vector3)
    yp_asc = yaxis_asc.vector3.dot(pos.vector3)
    nu = normalize_angle(atan2(yp_asc, xp_asc) - argp)

    e = ecc_vec.vector3.length
    ee = e * e
    eee = ee * e
    eeee = eee * e
    i = north.vector3.angle(kinetic.vector3.normalized())
    raan = normalize_angle(atan2(asc.vector3.y, asc.vector3.x))
    a = (kinetic_sq / (MU_EARTH * (1 - ee)))
    aaa = a ** 3
    n = (MU_EARTH / aaa) ** .5
    M = (
        nu
        - 2 * e * sin(nu)
        + (ee * .75 + eeee * .125) * sin(2 * nu)
        - eee * sin(3 * nu) / 3
        + eeee * sin(4 * nu) * .15625
    )

    return OrbitalElementsTuple(
        e,
        i,
        raan,
        argp,
        n,
        M
    )

class _TupleGPS(NamedTuple):
    latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    longitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]

@tensor_check
def itrf2gps(
    itrf_pos: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)]
) -> _TupleGPS:
    """Convert a position in ITRS coordinates to GPS coordinates, at a given time."""
    lon = itrf_pos.vector3.theta
    lat = itrf_pos.vector3.delta
    alt = itrf_pos.vector3.length - geocentric_radius_earth(lat)

    return _TupleGPS(
        latitude=lat,
        longitude=lon,
        altitude=alt
    )

class _TupleHorizontal(NamedTuple):
    azimuth: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    distance: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]

@tensor_check
def itrf2horizontal(
    itrf_pos: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)], 
    earth_local_frame: Any
) -> _TupleHorizontal:
    t = cast(
        TransformVector3Affine,
        earth_local_frame.transform
    )
    horizontal_pos = t.do(itrf_pos)

    return _TupleHorizontal(
        azimuth=horizontal_pos.vector3.theta,
        altitude=horizontal_pos.vector3.delta,
        distance=itrf_pos.vector3.length
    )

def convert_quantity_units(quantity: RealNumber, units_from: Union[str, "pint.Unit"], units_to: Union[str, "pint.Unit"]) -> RealNumber:
    """Converts a number `quantity` that is expressed in `units_from`, to `units_to` (uses Pint)"""
    import pint

    return cast(RealNumber, pint.Quantity(quantity, units_from).m_as(units_to))

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