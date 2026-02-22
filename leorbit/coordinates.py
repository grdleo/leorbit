from dataclasses import dataclass
from functools import cache, lru_cache
from statistics import mean
from typing import Literal, NamedTuple, Self

from enum import Enum
from functools import lru_cache
from typing import ParamSpec, Callable, TypeVar, cast

from leorbit.algorithms import OrbitalElementsComputeTuple
from leorbit.frames import AbsoluteFrame, EarthLocalFrame, frame_transform_factory, Frame
from leorbit.m import D, Dim, Quantity, Scalar, Vector3, Vector3Array, atan, normalize_angle, normalize_angle_symmetric, sqrt, square, tan
from leorbit.time import Time, TimeInterval
from leorbit.transforms import Transform, TransformVector3Affine
from leorbit.utils import angle2dms, eccentric2true_anomaly, elements2orthogonal_gcrf, gcrf_state_vectors2elements, geocentric_radius_earth, mean2eccentric_anomaly, mean_motion_to_semi_major_axis_earth

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

    def __repr__(self) -> str:
        km = "kilo_meter"
        kmph = "kilo_meter_per_hour"
        pos_repr = f"(x={self.pos.x.magnitude(km):.2f}, y={self.pos.y.magnitude(km):.2f}, z={self.pos.z.magnitude(km):.2f})"
        vel_repr = None if self.vel is None else f"(x={self.vel.x.magnitude(kmph):.2f}, y={self.vel.y.magnitude(kmph):.2f}, z={self.vel.z.magnitude(kmph):.2f})"
        return (
            f"<Coordinates: epoch={self.epoch}, "
            f"frame={self.privileged_frame}, "
            f"pos={pos_repr}, "
            f"vel={vel_repr}>" if vel_repr is not None else ""
            ">"
        )
    
    def __hash__(self) -> int:
        pos, vel = self.positions[self.privileged_frame]
        hash_str = f"{hash(self.epoch)}${hash(pos)}${hash(vel)}${self.name}"
        return hash(hash_str)

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
        mean_motion_dot: Scalar[D.AngularAcc] = Scalar[D.AngularAcc](0),
        mean_motion_ddot: Scalar[D.AngularJerk] = Scalar[D.AngularJerk](0),
        bstar: Scalar[D.InvLength] = Scalar[D.InvLength](0),
    ):
        deg_0 = 0 * Quantity.deg
        deg_180 = 180 * Quantity.deg
        deg_360 = 360 * Quantity.deg

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
        self.true_anomaly = eccentric2true_anomaly(e, E)
        self.semi_major_axis = mean_motion_to_semi_major_axis_earth(self.mean_motion)
        self.semi_minor_axis = self.semi_major_axis * sqrt(1 - square(e))
        dt = (self.mean_anomaly / self.mean_motion).cast(D.Time)
        self.time_at_periaster = self.epoch - dt

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
        pos_gcrf, vel_gcrf = elements2orthogonal_gcrf(
            self.true_anomaly,
            self.eccentricity,
            self.semi_major_axis,
            self.ra_of_asc_node,
            self.arg_of_pericenter,
            self.inclination
        )

        return Coordinates(
            self.epoch,
            AbsoluteFrame.GCRF, 
            pos_gcrf,
            vel_gcrf
        )

    @staticmethod
    def from_state_vectors(epoch: Time, pos: PosVec, vel: VelVec) -> "OrbitalElements":
        """
        From the state vectors of a given satellite (position and velocity, both condensed in a `Coordinates` object),
        returns a `OrbitalElements` object corresponding to its orbit.

        .. caution::
            Creating orbital elements this way discard any information about drag (`BSTAR` quantity in NORAD elements),
            and any other perturbations.
            Propagating those elements for more than a few hours will lead in huge errors.
        """
        els = gcrf_state_vectors2elements(pos, vel)

        return OrbitalElements(
            epoch=epoch,
            eccentricity=els.eccentricity,
            inclination=els.inclination,
            ra_of_asc_node=els.ra_of_asc_node,
            arg_of_pericenter=els.arg_of_pericenter,
            mean_motion=els.mean_motion,
            mean_anomaly=els.mean_anomaly,
            mean_motion_dot=Scalar[D.AngularAcc](0),
            mean_motion_ddot=Scalar[D.AngularJerk](0),
            bstar=Scalar[D.InvLength](0)
        )

    @staticmethod
    def from_celestrak_norad_cat_id(catnr: int, log: bool = False) -> "OrbitalElements":
        """
        Fetches lastest GP data on `celestrak.org` corresponding to given NORAD catalog ID, 
        and creates an instance of `OrbitalElements` using those.

        If last fetch on Celestrak is recent enough, uses cached GP data.
        """
        from leorbit.ext import get_celestrak_gpdata
        return get_celestrak_gpdata(catnr, log).to_orbital_elements()

class Interpolation(Enum):
    CONSTANT = "constant"
    """Snaps to closest"""

    LINEAR = "linear"
    """Linear interpolation between two closest"""

PosVecArray = Vector3Array[D.Length]
VelVecArray = Vector3Array[D.Velocity]

class PosVelArray(NamedTuple):
    pos: PosVecArray
    vel: VelVecArray | None

class Trajectory:
    def __init__(self, interval: TimeInterval, frame: Frame, pos: PosVecArray, vel: VelVecArray | None = None):
        if not (interval.steps == pos.size and (vel is None or interval.steps == vel.size)):
            raise ValueError("Size of position and velocity arrays must match the number of steps in the given interval.")
        
        self.interval = interval
        self.positions: dict[Frame, PosVelArray] = {frame: PosVelArray(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None
        
        self._already_computed_repr: dict[type[CoordinatesRepresentation], CoordinatesRepresentation] = {}

    def _compute_new_frame(self, frame: Frame):
        if frame in self.positions.keys():
            return
        
        positions, velocities = self.positions[self.privileged_frame]
        transformed_positions: list[PosVec] = []
        transformed_velocities: list[VelVec] = []

        for idx, epoch in enumerate(self.interval):
            p = positions[idx]
            v = velocities[idx] if velocities is not None else None

            transform = cast(
                Transform[DynamicVec, DynamicVec], 
                frame_transform_factory(self.privileged_frame, frame)(epoch)
            )

            transformed_positions.append(
                cast(PosVec, transform.do(p))
            )
            if v is not None:
                transformed_velocities.append(
                    cast(VelVec, transform.do(v))
                )
        
        self.positions[frame] = PosVelArray(
            Vector3Array[D.Length].from_vectors(*transformed_positions),
            (
                Vector3Array[D.Velocity].from_vectors(*transformed_velocities)
                if transformed_velocities 
                else None
            )
        )

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame

    @lru_cache(4096)
    def coordinates_at(self, epoch: Time, interpolation: Interpolation = Interpolation.CONSTANT) -> Coordinates:
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")
        
        idx = self.interval._time2idx(epoch) # check if epoch is in interval
        pos, vel = self.positions[self.privileged_frame]

        return Coordinates(
            epoch,
            self.privileged_frame,
            pos[idx],
            vel[idx] if vel is not None else None
        )

    @lru_cache(4096)
    def get_pos(self, epoch: Time, frame: Frame, interpolation: Interpolation = Interpolation.CONSTANT) -> PosVec: 
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")
        
        self._compute_new_frame(frame)
        
        idx = self.interval._time2idx(epoch) # check if epoch is in interval
        pos, vel = self.positions[frame]

        return pos[idx]
    
    @lru_cache(4096)
    def get_vel(self, epoch: Time, frame: Frame, interpolation: Interpolation = Interpolation.CONSTANT) -> VelVec:
        if self.vel_available is False:
            raise ValueError("Velocity data is not available for this trajectory.")
         
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")
        
        self._compute_new_frame(frame)
        
        idx = self.interval._time2idx(epoch) # check if epoch is in interval
        pos, vel = self.positions[frame]

        return cast(VelVecArray, vel)[idx]
    
    @lru_cache(4096)
    def gps_at(self, epoch: Time) -> GPS:
        return self.coordinates_at(epoch).gps()
    
    @lru_cache(4096)
    def horizontal_at(self, epoch: Time) -> Horizontal:
        return self.coordinates_at(epoch).horizontal()