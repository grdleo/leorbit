"""Special functions with special purposes. Should not be useful for the average user.
"""

from dataclasses import dataclass
from fractions import Fraction
from multiprocessing import Value
import numpy as np
import numpy.typing as npt

from typing import Annotated, Any, Generator, Iterable, NamedTuple, TypeAlias, TypeVar, Union, cast, overload, TYPE_CHECKING

from leorbit.transforms import TransformVector3Affine

if TYPE_CHECKING:
    from leorbit.time import TimeInterval
    from leorbit.frames import EarthLocalFrame
    import pint

from leorbit.mathematics import U, RealNumber, Tensor, TensorBound, TensorKind, atan2, cos, normalize_angle, scalar, sin, tan, vector3, tensor_check

MU_EARTH = scalar(398_600_441_800_000).with_units(U.meter ** 3 / U.second ** 2)
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = MU_EARTH ** .5
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

# Earth mean radius (WGS84) in meters
_radii_earth_unit = "radii_earth"
U.define(f"{_radii_earth_unit} = 6378137 * meter")

RADIIE_AA = scalar(6_378_137).with_units("meter") ** 2
RADIIE_A4 = RADIIE_AA ** 2
RADIIE_BB = scalar(6_356_752).with_units("meter") ** 2
RADIIE_B4 = RADIIE_BB ** 2

TWOPI = 2 * np.pi
TWELF_PI = np.pi / 12

#####################################

def angle2dms(angle: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle_deg = angle.scalar.value("deg")
    angle2convert = np.abs(angle_deg)
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{int(deg): 04}° {int(min):02}′ {int(sec):02}″"

def unixepoch_to_j2000(unixepoch: npt.NDArray | RealNumber) -> npt.NDArray:
    """Representation of this `U.second` object as "Julian year (J2000)", aka 
    the number of days since 2000/01/01T12:00:00."""
    return np.asarray(unixepoch) / 86_400 - 10_957.5

def j2000_to_stl0(j2000: npt.NDArray | RealNumber) -> npt.NDArray | RealNumber:
    """The 
    [Sideral U.second](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
    (angle) of Latitude 0 at this `U.second`, in radians.
    """
    return ((np.longdouble(18.697374558) + np.longdouble(24.06570982441908) * np.asarray(j2000)) * TWELF_PI) % TWOPI

def jd(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(units=U.second, kind=TensorKind.SCALAR)]:
    """Representation of this `Timestamp` object as "Julian day (JD)", aka 
    the number of days since -4712/01/01."""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    days = (unixepoch / 86_400 + 2_440_587.5)
    return Tensor(days).with_units(U.day)

def j2000(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(units=U.second, kind=TensorKind.SCALAR)]:
    """Representation of this `Timestamp` object as "Julian year (J2000)", aka 
    the number of days since 2000/01/01T12:00:00."""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    return Tensor(unixepoch_to_j2000(unixepoch)).with_units(U.day)

def from_mil(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(units=U.second, kind=TensorKind.SCALAR)]:
    """Representation as a fraction of days since 1 january 2000 00:00.

    Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
    unixepoch = np.asarray(unixepoch, dtype=np.float64)
    days = (unixepoch / 86_400 - 10_957.5) - .5
    return Tensor(days).with_units(U.day)

def stl0(
    unixepoch: int | float | npt.NDArray[np.float64]
) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
    """Sidereal time (angle) of latitude 0 at this unixepoch."""
    j2k = j2000(unixepoch)
    sidereal = j2000_to_stl0(j2k.raw_data_array("day"))
    return Tensor(np.asarray(sidereal, dtype=np.float64)).with_units(U.radian)

@tensor_check
def mean_motion_to_semi_major_axis_earth(
    mean_motion: Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
    """Convert Earth mean motion to semi-major axis using Kepler's third law."""
    return (MU_EARTH / mean_motion ** 2) ** Fraction(1, 3)

@tensor_check
def semi_major_axis_earth_to_mean_motion(
        sma: Annotated[Tensor,  TensorBound(units=U.meter, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)]:
    """Convert Earth semi-major axis to mean motion using Kepler's third law."""
    return (MU_EARTH / sma ** 3) ** .5

@tensor_check
def geocentric_radius_earth(
    latitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
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
#     dist_sat = abs(dir_obs)
#     phase = dir_sun.angle(dir_obs)

#     phi_term = sin(phase) + (HALF_REV - phase).m_as("rad") * cos(phase)

#     return std_mag + 5 * log10(dist_sat.m_as("megameter")) - 2.5 * log10(phi_term)

def humanize_duration(
    t: Annotated[Tensor, TensorBound(units=U.second, kind=TensorKind.SCALAR, size=1)]
) -> str:
    """Make a duration human-readable.
    
    Example:
    --------
    
    ```python
    humanize_duration(scalar(723459.23).with_units(U.minute))
    >>> '1 year 4 month 15 day 9 hour 39 min 13 s'
    ```
    """
    stages = ["year", "month", "day", "hour", "min", "s"]
    components = {s: 0 for s in stages}

    for s in stages:
        stage_d = scalar(1).with_units(U.parse_units(s))
        if t < stage_d:
            continue
        c = int(t.scalar.value(s))
        components[s] = c
        d = stage_d * c
        t = t - d
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)

@tensor_check
def mean2true_anomaly(
    e: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)],
    M: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
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
    e: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)],
    M: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
    """Solve Kepler's equation for eccentric anomaly by fixed-point iterations."""
    E = M
    for _ in range(5):
        E = M + e * sin(E)
    return E

@tensor_check
def eccentric2true_anomaly(
    e: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)], 
    E: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
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
    e: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)], 
    nu: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
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
    υ: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)], 
    e: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)], 
    a: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)], 
    Ω: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)], 
    ω: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
    i: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
) -> tuple[
    Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)], 
    Annotated[Tensor, TensorBound(units=U.meter / U.second, kind=TensorKind.VECTOR3)]
]:
    """Returns the position and velocity of a satellite in GCRF coordinates (meters, meters/second)
    All parameters are in radians, except `e` dimensionless and `a` in meters."""
    
    # position of satellite in orbit plane (with z = 0)
    one_ee = 1 - e ** 2
    E = true2eccentric_anomaly(e, υ)
    esinE = e * sin(E)
    
    r = a * one_ee / (1 + e * cos(υ))
    rd = (SQRT_MU_EARTH * a ** .5 * esinE / r).secure(units=U.meter / U.second)
    rυd = rd * one_ee / esinE
    
    c_raan, s_raan = cos(Ω), sin(Ω)
    c_i, s_i = cos(i), sin(i)
    υpω = υ + ω
    c_theta, s_theta = cos(υpω), sin(υpω)
     
    def vector_gcrf_factory(
        x: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)],
        y: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.VECTOR3)]:
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
    """Compact orbital-element tuple reconstructed from state vectors."""

    eccentricity: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]
    inclination: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    ra_of_asc_node: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    arg_of_pericenter: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    mean_motion: Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)]
    mean_anomaly: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]

@tensor_check
def gcrf_state_vectors2elements(
    pos: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)],
    vel: Annotated[Tensor, TensorBound(units=U.meter / U.second, kind=TensorKind.VECTOR3)]
) -> OrbitalElementsTuple:
    """Estimate Keplerian elements from GCRF position and velocity vectors."""
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

@dataclass
class GPSTrajectory:
    latitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    longitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]

    def to_csv(self, timeline: "TimeInterval") -> str:
        """Serialize the trajectory samples to a CSV string.

        The output uses one row per sample in ``timeline`` and the following
        header:

        ``timestamp,latitude,longitude,altitude``

        CSV columns:
        - ``timestamp``: ISO-8601 timestamp for the sample time
            (``U.second.isoformat``).
        - ``latitude``: geocentric latitude in degrees.
        - ``longitude``: geocentric longitude in degrees.
        - ``altitude``: altitude above the reference Earth spheroid in
            kilometers.

        Numeric columns are formatted with one decimal place.

        Parameters
        ----------
        timeline : TimeInterval
                Timeline used to pair each trajectory sample with its timestamp.
                Its number of steps must match the trajectory sample count.

        Returns
        -------
        str
                Complete CSV document including the header row.

        Raises
        ------
        ValueError
                If ``timeline.steps`` does not match the number of trajectory
                samples.
        """
        latitude_deg = self.latitude.scalar.values('deg').flatten()
        longitude_deg = self.longitude.scalar.values('deg').flatten()
        altitude_km = self.altitude.scalar.values('km').flatten()

        assert len(latitude_deg) == len(longitude_deg) == len(altitude_km)
        if len(latitude_deg) != timeline.steps:
            raise ValueError("Given timeline is not matching the size of trajectory")
        
        def _float2str(el: float) -> str:
            return f"{el:.1f}"
        
        line_gen = zip(
            map(_float2str, latitude_deg),
            map(_float2str, longitude_deg),
            map(_float2str, altitude_km),
            map(lambda t: t.isoformat, timeline)
        )

        lines = "\n".join(
            f"{ts},{lat},{lon},{alt}"
            for lat, lon, alt, ts in line_gen
        )
        
        return f"timestamp,latitude,longitude,altitude\n{lines}"

@tensor_check
def itrf2gps(
    itrf_pos: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)]
) -> GPSTrajectory:
    """Convert a position in ITRS coordinates to GPS coordinates, at a given time."""
    lon = itrf_pos.vector3.theta
    lat = itrf_pos.vector3.delta
    alt = itrf_pos.vector3.length - geocentric_radius_earth(lat)

    return GPSTrajectory(
        latitude=lat,
        longitude=lon,
        altitude=alt
    )

@dataclass
class HorizontalTrajectory:
    azimuth: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    distance: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]

    def to_csv(self, timeline: "TimeInterval") -> str:
        """Serialize the trajectory samples to a CSV string.

        The output uses one row per sample in ``timeline`` and the following
        header:

        ``timestamp,azimuth,elevation,range``

        CSV columns:
        - ``timestamp``: ISO-8601 timestamp for the sample time
          (``U.second.isoformat``).
        - ``azimuth``: horizontal azimuth angle in degrees.
        - ``elevation``: horizontal elevation/altitude angle in degrees.
        - ``range``: line-of-sight distance in kilometers.

        Numeric columns are formatted with one decimal place.

        Parameters
        ----------
        timeline : TimeInterval
            Timeline used to pair each trajectory sample with its timestamp.
            Its number of steps must match the trajectory sample count.

        Returns
        -------
        str
            Complete CSV document including the header row.

        Raises
        ------
        ValueError
            If ``timeline.steps`` does not match the number of trajectory
            samples.
        """
        azimuth_deg = self.azimuth.scalar.values('deg').flatten()
        elevation_deg = self.altitude.scalar.values('deg').flatten()
        _range_km = self.distance.scalar.values('km').flatten()

        assert len(azimuth_deg) == len(elevation_deg) == len(_range_km)
        if len(azimuth_deg) != timeline.steps:
            raise ValueError("Given timeline is not matching the size of trajectory")
        
        def _float2str(el: float) -> str:
            return f"{el:.1f}"
        
        line_gen = zip(
            map(_float2str, azimuth_deg),
            map(_float2str, elevation_deg),
            map(_float2str, _range_km),
            map(lambda t: t.isoformat, timeline)
        )

        lines = "\n".join(
            f"{ts},{az},{el},{_range}"
            for az, el, _range, ts in line_gen
        )
        
        return f"timestamp,azimuth,elevation,range\n{lines}"

@tensor_check
def itrf2horizontal(
    itrf_pos: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)], 
    earth_local_frame: Any
) -> HorizontalTrajectory:
    """Convert an ITRF position to horizontal coordinates for a local frame."""
    t = cast(
        TransformVector3Affine,
        earth_local_frame.transform
    )
    horizontal_pos = t.do(itrf_pos)

    return HorizontalTrajectory(
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