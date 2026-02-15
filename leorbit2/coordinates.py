from dataclasses import dataclass
from functools import cache, lru_cache
from typing import Literal, NamedTuple, Self

from enum import Enum
from functools import lru_cache
from typing import ParamSpec, Callable, TypeVar, cast

from leorbit2.algorithms import OrbitalElementsComputeTuple
from leorbit2.frames import AbsoluteFrame, EarthLocalFrame, frame_transform_factory, Frame
from leorbit2.m import D, Dim, Quantity, Scalar, Vector3
from leorbit2.time import Time, TimeInterval
from leorbit2.transforms import Transform, TransformVector3Affine
from leorbit2.utils import geocentric_radius_earth, mean2eccentric_anomaly, mean_motion_to_semi_major_axis_earth

PosVec = Vector3[D.Length]
VelVec = Vector3[D.Velocity]

SomeDim = TypeVar("SomeDim", bound=Dim)
DynamicVec = Vector3[D.Length] | Vector3[D.Velocity]
SomeDynamicVec = TypeVar("SomeDynamicVec", bound=DynamicVec)
DynamicD = D.Length | D.Velocity
SomeDynamicD = TypeVar("SomeDynamicD", bound=DynamicD)

class PosVel(NamedTuple):
	pos: PosVec
	vel: VelVec | None


class Coordinates:
    def __init__(self, epoch: Time, frame: Frame, pos: PosVec, vel: VelVec | None = None):
        self.epoch = epoch
        self.positions: dict[Frame, PosVel] = {frame: PosVel(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None
        
        self._already_computed_repr: dict[type[CoordinatesRepresentation], CoordinatesRepresentation] = {}

    def _compute_new_frame(self, frame: Frame):
        if frame in self.positions.keys():
            return
        
        p, v = self.positions[self.privileged_frame]
        transform = cast(Transform[DynamicVec, DynamicVec], frame_transform_factory(self.privileged_frame, frame)(self.epoch))
        p = cast(PosVec, transform.do(p))
        if v is not None:
            v = cast(VelVec, transform.do(v))

        self.positions[frame] = PosVel(p, v)

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame
    
    def get_pos(self, frame: Frame) -> PosVec:
        self._compute_new_frame(frame)
        return self.positions[frame].pos
    
    def get_vel(self, frame: Frame) -> VelVec | None:
        if not self.vel_available:
            return None
        self._compute_new_frame(frame)
        return self.positions[frame].vel
    
    @property
    def pos(self) -> PosVec:
        return self.get_pos(self.privileged_frame)
    
    @property
    def vel(self) -> VelVec | None:
        return self.get_vel(self.privileged_frame)
    
    ### GPS ###
    ### ### ###
    
    @staticmethod
    def from_gps(
        longitude: Scalar[D.Angle], 
        latitude: Scalar[D.Angle], 
        altitude: Scalar[D.Length] = 0 * Quantity.meter, 
        epoch: Time | None = None
    ) -> Coordinates:
        theta = longitude
        delta = latitude
        rho = geocentric_radius_earth(delta) + altitude
        pos = Vector3.from_spherical(theta, delta, rho)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, AbsoluteFrame.ITRF, pos)
    
    @lru_cache
    def gps(self) -> GPS:
        """Returns this coordinates as their GPS representation"""
        gps = self._already_computed_repr.get(GPS, None)
        if gps is not None:
            return cast(GPS, gps)
        
        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        lon = itrf_pos.theta
        lat = itrf_pos.delta
        alt = itrf_pos.length - geocentric_radius_earth(lat)

        return GPS(
            longitude=lon, 
            latitude=lat, 
            altitude=alt
        )
    
    ### ### ###

    ### HORIZONTAL ###
    ### ########## ###

    @staticmethod
    def from_horizontal(
        azimuth: Scalar[D.Angle], 
        altitude: Scalar[D.Angle], 
        distance: Scalar[D.Length], 
        frame: EarthLocalFrame, 
        epoch: Time | None = None
    ) -> Coordinates:
        assert isinstance(frame, EarthLocalFrame)
        
        pos_local = Vector3.from_spherical(azimuth, altitude, distance)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, frame, pos_local)
    
    @lru_cache
    def horizontal(self, local_frame: EarthLocalFrame) -> Horizontal:
        """Returns the representation of this coordinates as 'horizontal', in the given frame"""
        assert local_frame.reference_frame == AbsoluteFrame.ITRF # FIXME

        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        t = cast(
            TransformVector3Affine[D.Length],
            local_frame.transform
        )
        horizontal_pos = t.do(itrf_pos)

        return Horizontal(
            azimuth=horizontal_pos.theta,
            altitude=horizontal_pos.delta,
            distance=itrf_pos.length
        )
    
    ### ########## ###
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Coordinates):
            return False
        
        return (
            self.get_pos(AbsoluteFrame.ITRF) == o.get_pos(AbsoluteFrame.ITRF)
            and self.get_vel(AbsoluteFrame.ITRF) == o.get_vel(AbsoluteFrame.ITRF)
            and self.epoch == o.epoch
        )
    
    def __neq__(self, o: Self) -> bool:
        return not self.__eq__(o)
    
#################################
# REPRESENTATIONS

_HASH_VAL = "__hash_val"
    
class CoordinatesRepresentation:
    def __hash__(self) -> int:
        if not hasattr(self, _HASH_VAL):
            h = hash(
                "".join(
                    (
                        f"{attr_name}:{getattr(self, attr_name, None)}"
                        for attr_name in self.__annotations__.keys()
                    )
                )
            )

            setattr(self, _HASH_VAL, h)
        
        return getattr(self, _HASH_VAL)
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, CoordinatesRepresentation):
            return False
        
        return hash(self) == hash(o)

from dataclasses import dataclass
from functools import cached_property, lru_cache


class GPS(CoordinatesRepresentation):
    """Representation of a point on Earth (or around) using GPS standards.
    """
    longitude: Scalar[D.Angle]
    latitude: Scalar[D.Angle]
    altitude: Scalar[D.Length]

    def __init__(self,
        longitude: Scalar[D.Angle],
        latitude: Scalar[D.Angle],
        altitude: Scalar[D.Length],
    ):
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude
        
    @cached_property
    def dms(self) -> str:
        """Representation of the GPS coordinates in DSM notation (degrees, minutes, seconds)

        Example: `39° 17′ N, 76° 36′ O`"""
        lon = normalize_angle_symmetric(self.longitude)
        lat = normalize_angle_symmetric(self.latitude)

        return (
            angle2dms(lon) + ("E" if lon >= 0 else "O") 
            + ", "
            + angle2dms(lat) + ("N" if lat >= 0 else "S") 
        )
    
    def __repr__(self) -> str:
        return f"<GPS: {self.dms}>"
    
    @lru_cache
    def to_coordinates(self, epoch: Time) -> "Coordinates":
        coordinates = Coordinates.from_gps(
            longitude=self.longitude,
            latitude=self.latitude,
            altitude=self.altitude,
            epoch=epoch
        )
        coordinates._already_computed_repr[GPS] = self
        
        return coordinates
    
    @cached_property
    def earth_local_frame(self) -> EarthLocalFrame:
        return EarthLocalFrame(
            self.to_coordinates()
        )
    
class Horizontal(CoordinatesRepresentation):
    """Horizontal coordinates use a celestial sphere centered on the observer. 
    Azimuth is measured eastward from the north point of the horizon.
    Altitude is the angle above the horizon.

    If distance if ommited, the coordinates represents a direction in the sky, 
    and therefore cannot be converted to a 'real' coordinate.
    
    (see https://en.wikipedia.org/wiki/Horizontal_coordinate_system)"""

    azimuth: Scalar[D.Angle]
    altitude: Scalar[D.Angle]
    distance: Scalar[D.Length] | None

    def __init__(self,
        azimuth: Scalar[D.Angle],
        altitude: Scalar[D.Angle],
        distance: Scalar[D.Length] | None = None
    ):
        self.azimuth = azimuth
        self.altitude = altitude
        self.distance = distance

    @cached_property
    def dms(self) -> str:
        azi = normalize_angle(self.azimuth)
        alt = normalize_angle_symmetric(self.altitude)

        return (
            f"Azimuth: {angle2dms(azi)}, "
            f"Altitude: {'-' if alt < 0 else ''}{angle2dms(alt)}"
        )
    
    def __repr__(self) -> str:
        return f"<Horizontal: {self.dms}>"
    
    def to_coordinates(self, frame: EarthLocalFrame, epoch: Time) -> "Coordinates":
        if self.distance is None:
            raise ValueError("Cannot convert `Horizontal` representation with unset `distance` to `Coordinates`.")
        
        return Coordinates.from_horizontal(
            azimuth=self.azimuth,
            altitude=self.altitude,
            distance=self.distance,
            frame=frame,
            epoch=epoch
        )

"""Implementation of "orbital elements" of an object orbiting Earth."""

class OrbitalElements(CoordinatesRepresentation):
    """Dataclass holding orbital elements, at a given epoch, gathered from Celestrak.org 
    (also known as GP data)"""

    def __init__(self,
        epoch: Time,
        eccentricity: Scalar[D.Dimless],
        inclination: Scalar[D.Angle],
        ra_of_asc_node: Scalar[D.Angle],
        arg_of_pericenter: Scalar[D.Angle],
        mean_motion: Scalar[D.AngularVelocity],
        mean_anomaly: Scalar[D.Angle],
        mean_motion_dot: Scalar[D.AngularAcc] = Scalar[D.AngularAcc].new(0),
        mean_motion_ddot: Scalar[D.AngularJerk] = Scalar[D.AngularJerk].new(0),
        bstar: Scalar[D.InvLength] = Scalar[D.InvLength].new(0),
    ):
        deg_0 = 0 * Quantity.deg
        deg_180 = 180 * Quantity.deg
        deg_360 = 360 * Quantity.deg

        self.name: str = "No name"
        self.norad_cat_id: int | None = None

        self.epoch = epoch

        self.eccentricity = e = eccentricity
        if e < 0:
            raise ValueError()
        
        self.inclination = i = inclination
        if not (deg_0 <= i <= deg_180):
            raise ValueError()

        self.ra_of_asc_node = raan = ra_of_asc_node
        if not (deg_0 <= raan <= deg_360):
            raise ValueError()

        self.arg_of_pericenter = argp = arg_of_pericenter
        if not (deg_0 <= argp <= deg_360):
            raise ValueError()

        self.mean_motion = n = mean_motion
        if not n.check(D.AngularVelocity):
            raise ValueError()

        self.mean_anomaly = M = mean_anomaly
        if not (deg_0 <= M <= deg_360):
            raise ValueError()

        self.mean_motion_dot = mean_motion_dot
        if not mean_motion_dot.check(D.AngularAcc):
            raise ValueError()

        self.mean_motion_ddot = mean_motion_ddot
        if not mean_motion_ddot.check(D.AngularJerk):
            raise ValueError()

        self.bstar = bstar
        if not bstar.check(D.InvLength):
            raise ValueError()

        e = self.eccentricity
        M = self.mean_anomaly
        self.eccentric_anomaly = E = mean2eccentric_anomaly(e, M)

        tan_half_nu = sqrt((1 + e) / (1 - e)) * tan(.5 * E)
        self.true_anomaly = 2 * atan(tan_half_nu)

        self.semi_major_axis = mean_motion_to_semi_major_axis_earth(self.mean_motion)

        dt = (self.mean_anomaly / self.mean_motion).cast(D.Time)
        self.time_at_periaster = self.epoch - dt

        self.semi_minor_axis = self.semi_major_axis * sqrt(1 - square(e))

    @cached_property
    def compute_tuple(self) -> OrbitalElementsComputeTuple:
        n = self.mean_motion.base_unit_value # [rad/s]
        n = n * 60 # [rad/min]

        bstar = self.bstar.base_unit_value # [1/m]
        bstar = bstar * Quantity.radii_earth.base_unit_value # [1/Earth radii]

        return OrbitalElementsComputeTuple(
            n=n,
            i=self.inclination.base_unit_value,
            e=self.eccentricity.base_unit_value,
            argp=self.arg_of_pericenter.base_unit_value,
            raan=self.ra_of_asc_node.base_unit_value,
            M=self.mean_anomaly.base_unit_value,
            bstar=bstar
        )
    
    @property
    def period(self) -> Scalar[D.Time]:
        """The period of a full revolution."""
        from math import tau

        return ((tau * Quantity.rad) / self.mean_motion).cast(D.Time)
    
    def to_coordinates(self) -> Coordinates:
        """Returns current orbital elements, at given epoch, 
        as a `Coordinates` object."""
        # position of satellite in orbit plane (with z = 0)
        pos_vel_gcrf = elements2orthogonal_gcrf(
            self.true_anomaly.m_as("rad"),
            self.eccentricity.m_as("1"),
            self.semi_major_axis.m_as("m"),
            self.ra_of_asc_node.m_as("rad"),
            self.arg_of_pericenter.m_as("rad"),
            self.inclination.m_as("rad")
        )

        assert not pos_vel_gcrf.array_values

        return Coordinates(
            AbsoluteFrame.GCRF, 
            Vec3(pos_vel_gcrf.x, pos_vel_gcrf.y, pos_vel_gcrf.z, UREG.meter),
            Vec3(pos_vel_gcrf.vx, pos_vel_gcrf.vy, pos_vel_gcrf.vz, UREG.meter / UREG.second),
            self.epoch
        )

    @staticmethod
    def from_state_vectors(coordinate: Coordinates) -> "OrbitalElements":
        """
        From the state vectors of a given satellite (position and velocity, both condensed in a `Coordinates` object),
        returns a `OrbitalElements` object corresponding to its orbit.

        .. caution::
            Creating orbital elements this way discard any information about drag (`BSTAR` quantity in NORAD elements),
            and any other perturbations.
            Propagating those elements for more than a few hours will lead in huge errors.
        """
        raise NotImplementedError("should be updated")

        if coordinate.epoch is None:
            raise ValueError(f"Epoch of given coordinate {coordinate} needs to be specified.")
        if not coordinate.vel_specified:
            raise ValueError(f"Velocity of given coordinate {coordinate} needs to be specified.")
        
        pos = coordinate.to_gcrf()
        vel = coordinate.to_gcrf_vel()
        north = Vec3.zaxis()

        kinetic = pos.cross(vel)
        kinetic_sq = kinetic.sqr()
        pos_dir = pos.normalize()
        ecc_vec: Vec3 = vel.cross(kinetic) / MU_EARTH - pos_dir
        descending_node = kinetic.cross(north).normalize() # descending node line
        asc = -descending_node

        # create a 2D frame on the ellipsis, x along ascending node line
        xaxis_asc = asc.normalize()
        yaxis_asc = kinetic.cross(asc).normalize()

        xe_asc = xaxis_asc.dot(ecc_vec).m
        ye_asc = yaxis_asc.dot(ecc_vec).m
        argp = atan2(ye_asc, xe_asc) * UREG("rad")

        xp_asc = xaxis_asc.dot(pos).m_as("m")
        yp_asc = yaxis_asc.dot(pos).m_as("m")
        nu = (atan2(yp_asc, xp_asc) * UREG("rad") - argp) % tau

        e = abs(ecc_vec) # [1]
        ee = e * e
        eee = ee * e
        eeee = eee * e
        i = north.angle(kinetic) # [rad]
        raan = atan2(asc.y, asc.x) % tau # [rad]
        raan = raan * UREG("rad") # convert to Quantity
        a = kinetic_sq / (MU_EARTH * (1 - ee)) # [m]
        aaa = a**3
        n = ((MU_EARTH / aaa)**0.5).to("rad/s")  # [rad/s]
        M = (
            nu
            - 2 * e * sin(nu)
            + (ee * .75 + eeee * .125) * sin(2 * nu)
            - eee * sin(3 * nu) / 3
            + eeee * sin(4 * nu) * .15625
        )

        return OrbitalElements(
            coordinate.epoch,
            e,
            i,
            raan,
            argp,
            n,
            M
        )

    @staticmethod
    def from_celestrak_json(celestrak_json: str | dict[str, str | float | int]) -> "OrbitalElements":
        """Creates an instance of `OrbitalElements` object from a Celestrak query in JSON format"""
        query = celestrak_json
        if isinstance(query, str):
            query = json.loads(celestrak_json)
            if isinstance(query, list):
                query = query[0]
        
        if not isinstance(query, dict):
            raise ValueError(f"Given argument {celestrak_json} is not an acceptable JSON Celestrak query")

        try:
            return OrbitalElements(
                epoch=Time.fromisoformat(query["EPOCH"]),
                eccentricity=query["ECCENTRICITY"] * UREG("dimensionless"),
                inclination=query["INCLINATION"] * UREG("°"),
                ra_of_asc_node=query["RA_OF_ASC_NODE"] * UREG("°"),
                arg_of_pericenter=query["ARG_OF_PERICENTER"] * UREG("°"),
                mean_motion=query["MEAN_MOTION"] * UREG("turn/day"),
                mean_anomaly=query["MEAN_ANOMALY"] * UREG("°"),
                mean_motion_dot=query["MEAN_MOTION_DOT"] * UREG("turn/day^2") * 2, # NOTE: factor is cancelled when loading directly from Celestrack
                mean_motion_ddot=query["MEAN_MOTION_DDOT"] * UREG("turn/day^3") * 6, # NOTE: factor is cancelled when loading directly from Celestrack
                bstar=query["BSTAR"] * UREG("1/earthRadii"),
                name=query["OBJECT_NAME"],
                norad_cat_id=int(query["NORAD_CAT_ID"])
            )
        except KeyError as ex:
            raise ValueError(f"Given argument {celestrak_json} is not an acceptable JSON Celestrak query")

    @staticmethod
    def from_celestrak_norad_cat_id(catnr: int, log: bool = False) -> "OrbitalElements":
        """
        Fetches lastest GP data on `celestrak.org` corresponding to given NORAD catalog ID, 
        and creates an instance of `OrbitalElements` using those.

        If last fetch on Celestrak is recent enough, uses cached GP data.
        """
        gp_dict = get_celestrak_gpdata_json(catnr, log)
        return OrbitalElements.from_celestrak_json(gp_dict)
    





class Interpolation(Enum):
    SNAP = "snap"
    """Snaps to closest"""

    LINEAR = "linear"
    """Linear interpolation between two closest"""

class PosVelArray(NamedTuple):
    pos: "PosVecArray"
    vel: "VelVecArray"

class Trajectory:
    def __init__(self, interval: TimeInterval, frame: Frame, pos: "PosVecArray", vel: "VelVecArray" | None = None):
        self.interval = interval
        self.positions: dict[Frame, PosVel] = {frame: PosVelArray(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None
        
        self._already_computed_repr: dict[type[CoordinatesRepresentation], CoordinatesRepresentation] = {}

    def coordinates_at(self, epoch: Time, interpolation: Interpolation = Interpolation.SNAP) -> Coordinates:
        raise NotImplementedError()

    def get_pos(self, epoch: Time, frame: Frame, interpolation: Interpolation = Interpolation.SNAP) -> PosVec: 
        raise NotImplementedError()
    
    def get_vel(self, epoch: Time, frame: Frame, interpolation: Interpolation = Interpolation.SNAP) -> VelVec: 
        raise NotImplementedError()
    
    def gps(self) -> dict[Time, GPS]:
        raise NotImplementedError()
    
    def horizontal(self) -> dict[Time, Horizontal]:
        raise NotImplementedError()