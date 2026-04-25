from enum import Enum
from functools import cached_property, lru_cache
from typing import Annotated, NamedTuple, Self, cast

from leorbit.algorithms import OrbitalElementsComputeTuple
from leorbit.frames import AbsoluteFrame, EarthLocalFrame, Frame, frame_transform_factory
from leorbit.mathematics import (
    Angle,
    AngularAcceleration,
    AngularJerk,
    AngularVelocity,
    Dimless,
    InvLength,
    Length,
    Quantity,
    Tensor,
    TensorAsVector3,
    TensorBound,
    TensorKind,
    Time,
    normalize_angle,
    normalize_angle_symmetric,
)
from leorbit.time import TimeInterval, Timestamp
from leorbit.transforms import Transform
from leorbit.utils import (
    angle2dms,
    eccentric2true_anomaly,
    elements2orthogonal_gcrf,
    gcrf_state_vectors2elements,
    geocentric_radius_earth,
    itrf2gps,
    itrf2horizontal,
    mean2eccentric_anomaly,
    mean_motion_to_semi_major_axis_earth,
)

PosVec = Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)]
VelVec = Annotated[Tensor, TensorBound(dimension=Length / Time, kind=TensorKind.VECTOR3)]
PosVecArray = Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)]
VelVecArray = Annotated[Tensor, TensorBound(dimension=Length / Time, kind=TensorKind.VECTOR3)]


class PosVel(NamedTuple):
    pos: PosVec
    vel: VelVec | None


class Coordinates:
    def __init__(self, epoch: Timestamp, frame: Frame, pos: PosVec, vel: VelVec | None = None):
        self.epoch = epoch
        self.positions: dict[Frame, PosVel] = {frame: PosVel(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None

    def __repr__(self) -> str:
        km = "kilo_meter"
        kmph = "kilo_meter / hour"
        p = self.pos.vector3
        pos_repr = (
            f"(x={p.x.scalar.value(km):.2f}, y={p.y.scalar.value(km):.2f}, z={p.z.scalar.value(km):.2f})"
        )
        vel_repr = "None"
        if self.vel is not None:
            v = self.vel.vector3
            vel_repr = (
                f"(x={v.x.scalar.value(kmph):.2f}, y={v.y.scalar.value(kmph):.2f}, z={v.z.scalar.value(kmph):.2f})"
            )
        return f"<Coordinates: epoch={self.epoch}, frame={self.privileged_frame}, pos={pos_repr}, vel={vel_repr}>"

    def __hash__(self) -> int:
        pos, vel = self.positions[self.privileged_frame]
        hash_str = f"{self.__class__.__name__}${hash(self.epoch)}${hash(pos)}${hash(vel)}${self.name}"
        return hash(hash_str)

    def _compute_new_frame(self, frame: Frame):
        if frame in self.positions:
            return

        p, v = self.positions[self.privileged_frame]
        transform = cast(Transform, frame_transform_factory(self.privileged_frame, frame)(self.epoch))
        p = transform.do(p)
        if v is not None:
            v = transform.do(v)

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

    @staticmethod
    def from_gps(
        longitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)] = 0 * Quantity.meter,
        epoch: Timestamp | None = None,
    ) -> "Coordinates":
        theta = longitude
        delta = latitude
        rho = geocentric_radius_earth(delta) + altitude
        pos = TensorAsVector3.from_spherical(theta, delta, rho)
        return Coordinates(Timestamp.now() if epoch is None else epoch, AbsoluteFrame.ITRF, pos)

    @lru_cache
    def gps(self) -> "GPS":
        tuple_gps = itrf2gps(self.get_pos(AbsoluteFrame.ITRF))
        return GPS(
            longitude=tuple_gps.longitude,
            latitude=tuple_gps.latitude,
            altitude=tuple_gps.altitude,
            epoch=self.epoch,
        )

    @staticmethod
    def from_horizontal(
        azimuth: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        distance: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)],
        frame: EarthLocalFrame,
        epoch: Timestamp | None = None,
    ) -> "Coordinates":
        pos_local = TensorAsVector3.from_spherical(azimuth, altitude, distance)
        return Coordinates(Timestamp.now() if epoch is None else epoch, frame, pos_local)

    @lru_cache
    def horizontal(self, local_frame: EarthLocalFrame) -> "Horizontal":
        tuple_hor = itrf2horizontal(self.get_pos(AbsoluteFrame.ITRF), local_frame)
        return Horizontal(
            azimuth=tuple_hor.azimuth,
            altitude=tuple_hor.altitude,
            distance=tuple_hor.distance,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Coordinates):
            return False
        return (
            self.get_pos(AbsoluteFrame.ITRF) == other.get_pos(AbsoluteFrame.ITRF)
            and self.get_vel(AbsoluteFrame.ITRF) == other.get_vel(AbsoluteFrame.ITRF)
            and self.epoch == other.epoch
        )

    def __neq__(self, other: Self) -> bool:
        return not self.__eq__(other)


_HASH_VAL = "__hash_val"


class CoordinatesRepresentation:
    def __hash__(self) -> int:
        if not hasattr(self, _HASH_VAL):
            annotations = getattr(type(self), "__annotations__", {})
            h = hash(
                f"${type(self).__name__}$".join(
                    f"{attr_name}:{getattr(self, attr_name, None)}" for attr_name in annotations.keys()
                )
            )
            setattr(self, _HASH_VAL, h)
        return getattr(self, _HASH_VAL)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CoordinatesRepresentation):
            return False
        return hash(self) == hash(other)


class GPS(CoordinatesRepresentation):
    longitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]

    def __init__(
        self,
        longitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)],
        epoch: Timestamp | None = None,
    ):
        self.longitude = normalize_angle_symmetric(longitude)
        self.latitude = normalize_angle_symmetric(latitude)
        self.altitude = altitude
        self.epoch = epoch

    @cached_property
    def dms(self) -> str:
        lon = self.longitude
        lat = self.latitude
        zero_ang = 0 * Quantity.radian
        return (
            angle2dms(abs(lon))
            + ("E" if lon >= zero_ang else "O")
            + ", "
            + angle2dms(abs(lat))
            + ("N" if lat >= zero_ang else "S")
        )

    def __repr__(self) -> str:
        return f"<GPS: {self.dms}>"

    @lru_cache(4096)
    def to_coordinates(self, epoch: Timestamp | None = None) -> Coordinates:
        effective_epoch = epoch if epoch is not None else (self.epoch if self.epoch is not None else Timestamp.now())
        return Coordinates.from_gps(
            longitude=self.longitude,
            latitude=self.latitude,
            altitude=self.altitude,
            epoch=effective_epoch,
        )

    @cached_property
    def earth_local_frame(self) -> EarthLocalFrame:
        return EarthLocalFrame(self.to_coordinates())


class Horizontal(CoordinatesRepresentation):
    azimuth: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
    distance: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)] | None

    def __init__(
        self,
        azimuth: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        distance: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)] | None = None,
    ):
        self.azimuth = azimuth
        self.altitude = altitude
        self.distance = distance

    @cached_property
    def dms(self) -> str:
        azi = normalize_angle(self.azimuth)
        alt = normalize_angle_symmetric(self.altitude)
        zero_ang = 0 * Quantity.radian
        return f"Azimuth: {angle2dms(azi)}, Altitude: {'-' if alt < zero_ang else ''}{angle2dms(abs(alt))}"

    def __repr__(self) -> str:
        return f"<Horizontal: {self.dms}>"

    def to_coordinates(self, frame: EarthLocalFrame, epoch: Timestamp) -> Coordinates:
        if self.distance is None:
            raise ValueError("Cannot convert `Horizontal` representation with unset `distance` to `Coordinates`.")
        return Coordinates.from_horizontal(
            azimuth=self.azimuth,
            altitude=self.altitude,
            distance=self.distance,
            frame=frame,
            epoch=epoch,
        )


class OrbitalElements(CoordinatesRepresentation):
    def __init__(
        self,
        epoch: Timestamp,
        eccentricity: Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)],
        inclination: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        ra_of_asc_node: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        arg_of_pericenter: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        mean_motion: Annotated[Tensor, TensorBound(dimension=AngularVelocity, kind=TensorKind.SCALAR)],
        mean_anomaly: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)],
        mean_motion_dot: Annotated[Tensor, TensorBound(dimension=AngularAcceleration, kind=TensorKind.SCALAR)] = 0
        * (Quantity.radian / Quantity.second ** 2),
        mean_motion_ddot: Annotated[Tensor, TensorBound(dimension=AngularJerk, kind=TensorKind.SCALAR)] = 0
        * (Quantity.radian / Quantity.second ** 3),
        bstar: Annotated[Tensor, TensorBound(dimension=InvLength, kind=TensorKind.SCALAR)] = 0 * (1 / Quantity.radii_earth),
    ):
        deg_0 = 0 * Quantity.degree
        deg_180 = 180 * Quantity.degree
        deg_360 = 360 * Quantity.degree

        if eccentricity < 0 * Quantity.dimensionless:
            raise ValueError()
        if not (deg_0 <= inclination <= deg_180):
            raise ValueError()
        if not (deg_0 <= ra_of_asc_node <= deg_360):
            raise ValueError()
        if not (deg_0 <= arg_of_pericenter <= deg_360):
            raise ValueError()
        if not mean_motion.check(AngularVelocity, TensorKind.SCALAR):
            raise ValueError()
        if not (deg_0 <= mean_anomaly <= deg_360):
            raise ValueError()
        if not mean_motion_dot.check(AngularAcceleration, TensorKind.SCALAR):
            raise ValueError()
        if not mean_motion_ddot.check(AngularJerk, TensorKind.SCALAR):
            raise ValueError()
        if not bstar.check(InvLength, TensorKind.SCALAR):
            raise ValueError()

        self.epoch = epoch
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.ra_of_asc_node = ra_of_asc_node
        self.arg_of_pericenter = arg_of_pericenter
        self.mean_motion = mean_motion
        self.mean_anomaly = mean_anomaly
        self.mean_motion_dot = mean_motion_dot
        self.mean_motion_ddot = mean_motion_ddot
        self.bstar = bstar

        e = self.eccentricity
        M = self.mean_anomaly
        self.eccentric_anomaly = mean2eccentric_anomaly(e, M)
        self.true_anomaly = eccentric2true_anomaly(e, self.eccentric_anomaly)
        self.semi_major_axis = mean_motion_to_semi_major_axis_earth(self.mean_motion)
        self.semi_minor_axis = self.semi_major_axis * (1 - e * e) ** 0.5
        dt = (self.mean_anomaly / self.mean_motion).secure(Time, TensorKind.SCALAR)
        self.time_at_periaster = self.epoch - dt

    @cached_property
    def compute_tuple(self) -> OrbitalElementsComputeTuple:
        return OrbitalElementsComputeTuple(
            n=self.mean_motion.scalar.value(Quantity.radian / Quantity.minute),
            i=self.inclination.scalar.value(Quantity.radian),
            e=self.eccentricity.scalar.value(Quantity.dimensionless),
            argp=self.arg_of_pericenter.scalar.value(Quantity.radian),
            raan=self.ra_of_asc_node.scalar.value(Quantity.radian),
            M=self.mean_anomaly.scalar.value(Quantity.radian),
            bstar=self.bstar.scalar.value(1 / Quantity.radii_earth),
        )

    @property
    def period(self) -> Annotated[Tensor, TensorBound(dimension=Time, kind=TensorKind.SCALAR)]:
        from math import tau
        return ((tau * Quantity.radian) / self.mean_motion).secure(Time, TensorKind.SCALAR)

    def to_coordinates(self) -> Coordinates:
        pos_gcrf, vel_gcrf = elements2orthogonal_gcrf(
            self.true_anomaly,
            self.eccentricity,
            self.semi_major_axis,
            self.ra_of_asc_node,
            self.arg_of_pericenter,
            self.inclination,
        )
        return Coordinates(self.epoch, AbsoluteFrame.GCRF, pos_gcrf, vel_gcrf)

    @staticmethod
    def from_state_vectors(epoch: Timestamp, pos: PosVec, vel: VelVec) -> "OrbitalElements":
        els = gcrf_state_vectors2elements(pos, vel)
        return OrbitalElements(
            epoch=epoch,
            eccentricity=els.eccentricity,
            inclination=els.inclination,
            ra_of_asc_node=els.ra_of_asc_node,
            arg_of_pericenter=els.arg_of_pericenter,
            mean_motion=els.mean_motion,
            mean_anomaly=els.mean_anomaly,
            mean_motion_dot=0 * (Quantity.radian / Quantity.second ** 2),
            mean_motion_ddot=0 * (Quantity.radian / Quantity.second ** 3),
            bstar=0 * (1 / Quantity.radii_earth),
        )

    @staticmethod
    def from_celestrak_norad_cat_id(catnr: int, log: bool = False) -> "OrbitalElements":
        from leorbit.ext import get_celestrak_gpdata

        return get_celestrak_gpdata(catnr, log).to_orbital_elements()


class Interpolation(Enum):
    CONSTANT = "constant"
    LINEAR = "linear"


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

    def __hash__(self) -> int:
        pos, vel = self.positions[self.privileged_frame]
        hash_str = f"{self.__class__.__name__}${hash(self.interval)}${hash(pos)}${hash(vel)}${self.name}"
        return hash(hash_str)

    def _compute_new_frame(self, frame: Frame):
        if frame in self.positions:
            return

        p, v = self.positions[self.privileged_frame]
        transform = cast(Transform, frame_transform_factory(self.privileged_frame, frame)(self.interval))

        p = transform.do(p)
        if v is not None:
            v = transform.do(v)

        self.positions[frame] = PosVelArray(p, v)

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame

    @lru_cache(4096)
    def coordinates_at(self, epoch: Timestamp, interpolation: Interpolation = Interpolation.CONSTANT) -> Coordinates:
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        idx = self.interval._time2idx(epoch)
        pos, vel = self.positions[self.privileged_frame]
        return Coordinates(epoch, self.privileged_frame, pos[idx], vel[idx] if vel is not None else None)

    @lru_cache(4096)
    def get_pos(self, epoch: Timestamp, frame: Frame, interpolation: Interpolation = Interpolation.CONSTANT) -> PosVec:
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        self._compute_new_frame(frame)
        idx = self.interval._time2idx(epoch)
        pos, _ = self.positions[frame]
        return pos[idx]

    @lru_cache(4096)
    def get_vel(self, epoch: Timestamp, frame: Frame, interpolation: Interpolation = Interpolation.CONSTANT) -> VelVec:
        if not self.vel_available:
            raise ValueError("Velocity data is not available for this trajectory.")
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        self._compute_new_frame(frame)
        idx = self.interval._time2idx(epoch)
        _, vel = self.positions[frame]
        return cast(Tensor, vel)[idx]

    @lru_cache(16)
    def trajectory_pos(self, frame: Frame) -> PosVecArray:
        self._compute_new_frame(frame)
        pos, _ = self.positions[frame]
        return pos

    @lru_cache(16)
    def trajectory_vel(self, frame: Frame) -> VelVecArray:
        if not self.vel_available:
            raise ValueError("Velocity data is not available for this trajectory.")

        self._compute_new_frame(frame)
        _, vel = self.positions[frame]
        return cast(Tensor, vel)

    class _GPSArray(NamedTuple):
        longitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
        latitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
        altitude: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]

    @lru_cache(16)
    def gps(self) -> _GPSArray:
        tuple_gps = itrf2gps(self.trajectory_pos(AbsoluteFrame.ITRF))
        return cast(Trajectory._GPSArray, tuple_gps)

    @lru_cache(4096)
    def gps_at(self, epoch: Timestamp, interpolation: Interpolation = Interpolation.CONSTANT) -> GPS:
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        gps_array = self.gps()
        idx = self.interval._time2idx(epoch)
        return GPS(
            longitude=gps_array.longitude[idx],
            latitude=gps_array.latitude[idx],
            altitude=gps_array.altitude[idx],
            epoch=epoch,
        )

    class _HorizontalArray(NamedTuple):
        azimuth: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
        altitude: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]
        distance: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]

    @lru_cache(16)
    def horizontal(self, local_frame: EarthLocalFrame) -> _HorizontalArray:
        tuple_hor = itrf2horizontal(self.trajectory_pos(AbsoluteFrame.ITRF), local_frame)
        return cast(Trajectory._HorizontalArray, tuple_hor)

    @lru_cache(4096)
    def horizontal_at(
        self,
        epoch: Timestamp,
        local_frame: EarthLocalFrame,
        interpolation: Interpolation = Interpolation.CONSTANT,
    ) -> Horizontal:
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        horizontal_array = self.horizontal(local_frame)
        idx = self.interval._time2idx(epoch)
        return Horizontal(
            azimuth=horizontal_array.azimuth[idx],
            altitude=horizontal_array.altitude[idx],
            distance=horizontal_array.distance[idx],
        )
