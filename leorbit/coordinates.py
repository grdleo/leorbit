from enum import Enum
from functools import cached_property, lru_cache
from typing import Annotated, NamedTuple, Self, cast

from leorbit.algorithms import OrbitalElementsComputeTuple
from leorbit.frames import AbsoluteFrame, EarthLocalFrame, Frame, frame_transform_factory
from leorbit.mathematics import U, Tensor, TensorAsVector3, TensorBound, TensorKind, normalize_angle, normalize_angle_symmetric, scalar
from leorbit.time import TimeInterval, Timestamp
from leorbit.transforms import Transform, _secure_position_transform, _secure_velocity_transform
from leorbit.utils import (
    GPSTrajectory,
    HorizontalTrajectory,
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
PosVec = Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)]
VelVec = Annotated[Tensor, TensorBound(units=U.meter / U.second, kind=TensorKind.VECTOR3)]
PosVecArray = Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.VECTOR3)]
VelVecArray = Annotated[Tensor, TensorBound(units=U.meter / U.second, kind=TensorKind.VECTOR3)]


class PosVel(NamedTuple):
    """Pair of position vector and optional velocity vector at one epoch."""

    pos: PosVec
    vel: VelVec | None


class Coordinates:
    """State vector at a single epoch with lazy frame conversions."""

    def __init__(self, epoch: Timestamp, frame: Frame, pos: PosVec, vel: VelVec | None = None):
        """Create coordinates at ``epoch`` in ``frame`` from position/velocity."""
        self.epoch = epoch
        self.positions: dict[Frame, PosVel] = {frame: PosVel(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None

    def __repr__(self) -> str:
        km = U.kilometer
        kmph = U.kilometer / U.hour
        p = self.pos.vector3
        pos_repr = (
            f"(x={p.x.scalar.value(km):.0f}km, y={p.y.scalar.value(km):.0f}km, z={p.z.scalar.value(km):.0f}km)"
        )
        vel_repr = "None"
        if self.vel is not None:
            v = self.vel.vector3
            vel_repr = (
                f"(x={v.x.scalar.value(kmph):.0f} km/h, y={v.y.scalar.value(kmph):.0f} km/h, z={v.z.scalar.value(kmph):.0f} km/h)"
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
        transform = frame_transform_factory(self.privileged_frame, frame)(self.epoch)
        p = _secure_position_transform(transform).do(p)
        v = None if v is None else _secure_velocity_transform(transform).do(v)

        self.positions[frame] = PosVel(p, v)

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame

    def get_pos(self, frame: Frame) -> PosVec:
        """Return the position expressed in ``frame``."""
        self._compute_new_frame(frame)
        return self.positions[frame].pos

    def get_vel(self, frame: Frame) -> VelVec | None:
        """Return the velocity expressed in ``frame`` when available."""
        if not self.vel_available:
            return None
        self._compute_new_frame(frame)
        return self.positions[frame].vel

    @property
    def pos(self) -> PosVec:
        """Position in the current privileged frame."""
        return self.get_pos(self.privileged_frame)

    @property
    def vel(self) -> VelVec | None:
        """Velocity in the current privileged frame, if available."""
        return self.get_vel(self.privileged_frame)

    @staticmethod
    def from_gps(
        longitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        latitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)] = scalar(0).with_units(U.meter),
        epoch: Timestamp | None = None,
    ) -> "Coordinates":
        """Create ITRF coordinates from geodetic longitude/latitude/altitude."""
        theta = longitude
        delta = latitude
        rho = geocentric_radius_earth(delta) + altitude
        pos = TensorAsVector3.from_spherical(theta, delta, rho)
        return Coordinates(Timestamp.now() if epoch is None else epoch, AbsoluteFrame.ITRF, pos)

    @lru_cache
    def gps(self) -> "GPS":
        """Convert coordinates to a geodetic ``GPS`` representation."""
        tuple_gps = itrf2gps(self.get_pos(AbsoluteFrame.ITRF))
        return GPS(
            longitude=tuple_gps.longitude,
            latitude=tuple_gps.latitude,
            altitude=tuple_gps.altitude,
            epoch=self.epoch,
        )

    @staticmethod
    def from_horizontal(
        azimuth: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        distance: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)],
        frame: EarthLocalFrame,
        epoch: Timestamp | None = None,
    ) -> "Coordinates":
        """Create coordinates from local horizontal angles and range."""
        pos_local = TensorAsVector3.from_spherical(azimuth, altitude, distance)
        return Coordinates(Timestamp.now() if epoch is None else epoch, frame, pos_local)

    @lru_cache
    def horizontal(self, local_frame: EarthLocalFrame) -> "Horizontal":
        """Convert coordinates to horizontal coordinates in ``local_frame``."""
        tuple_hor = itrf2horizontal(self.get_pos(AbsoluteFrame.ITRF), local_frame)
        return Horizontal(
            azimuth=tuple_hor.azimuth,
            altitude=tuple_hor.altitude,
            distance=tuple_hor.distance,
        )
    
    def gcrf(self) -> "TensorAsVector3":
        """Return position in GCRF frame."""
        return self.get_pos(AbsoluteFrame.GCRF).vector3
    
    def itrf(self) -> "TensorAsVector3":
        """Return position in ITRF frame."""
        return self.get_pos(AbsoluteFrame.ITRF).vector3

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
    """Hash/equality helper mixin for immutable coordinate-like value objects."""

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
    """Geodetic coordinates (longitude, latitude, altitude) on Earth."""

    longitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    latitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]

    def __init__(
        self,
        longitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        latitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)],
        epoch: Timestamp | None = None,
    ):
        """Create a GPS coordinate, normalizing angular components."""
        self.longitude = normalize_angle_symmetric(longitude)
        self.latitude = normalize_angle_symmetric(latitude)
        self.altitude = altitude
        self.epoch = epoch

    @cached_property
    def dms(self) -> str:
        """Return a human-readable DMS representation with cardinal letters."""
        lon = self.longitude
        lat = self.latitude
        zero_ang = scalar(0).with_units(U.radian)
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
        """Convert this geodetic point to ``Coordinates`` in the ITRF frame."""
        effective_epoch = epoch if epoch is not None else (self.epoch if self.epoch is not None else Timestamp.now())
        return Coordinates.from_gps(
            longitude=self.longitude,
            latitude=self.latitude,
            altitude=self.altitude,
            epoch=effective_epoch,
        )

    @cached_property
    def earth_local_frame(self) -> EarthLocalFrame:
        """Return the local topocentric frame centered on this GPS position."""
        return EarthLocalFrame(self.to_coordinates())


class Horizontal(CoordinatesRepresentation):
    """Horizontal coordinates (azimuth, altitude, optional distance)."""

    azimuth: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    altitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]
    distance: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)] | None

    def __init__(
        self,
        azimuth: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        altitude: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        distance: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)] | None = None,
    ):
        self.azimuth = azimuth
        self.altitude = altitude
        self.distance = distance

    @cached_property
    def dms(self) -> str:
        """Return a human-readable azimuth/altitude string."""
        azi = normalize_angle(self.azimuth)
        alt = normalize_angle_symmetric(self.altitude)
        zero_ang = scalar(0).with_units(U.radian)
        return f"Azimuth: {angle2dms(azi)}, Altitude: {'-' if alt < zero_ang else ''}{angle2dms(abs(alt))}"

    def __repr__(self) -> str:
        return f"<Horizontal: {self.dms}>"

    def to_coordinates(self, frame: EarthLocalFrame, epoch: Timestamp) -> Coordinates:
        """Convert this horizontal observation to 3D coordinates.

        Requires ``distance`` to be defined.
        """
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
    """Keplerian elements describing an Earth-centered osculating orbit."""

    def __init__(
        self,
        epoch: Timestamp,
        eccentricity: Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)],
        inclination: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        ra_of_asc_node: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        arg_of_pericenter: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        mean_motion: Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)],
        mean_anomaly: Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)],
        mean_motion_dot: Annotated[Tensor, TensorBound(units=U.radian / U.second ** 2, kind=TensorKind.SCALAR)] = scalar(0).with_units(U.radian / U.second ** 2),
        mean_motion_ddot: Annotated[Tensor, TensorBound(units=U.radian / U.second ** 3, kind=TensorKind.SCALAR)] = scalar(0).with_units(U.radian / U.second ** 3),
        bstar: Annotated[Tensor, TensorBound(units=1 / U.meter, kind=TensorKind.SCALAR)] = scalar(0).with_units(1 / U.meter),
    ):
        deg_0 = scalar(0).with_units(U.degree)
        deg_180 = scalar(180).with_units(U.degree)
        deg_360 = scalar(360).with_units(U.degree)

        if eccentricity < scalar(0).with_units(U.dimensionless):
            raise ValueError()
        if not (deg_0 <= inclination <= deg_180):
            raise ValueError()
        if not (deg_0 <= ra_of_asc_node <= deg_360):
            raise ValueError()
        if not (deg_0 <= arg_of_pericenter <= deg_360):
            raise ValueError()
        if not mean_motion.check(U.radian / U.second, TensorKind.SCALAR):
            raise ValueError()
        if not (deg_0 <= mean_anomaly <= deg_360):
            raise ValueError()
        if not mean_motion_dot.check(U.radian / U.second ** 2, TensorKind.SCALAR):
            raise ValueError()
        if not mean_motion_ddot.check(U.radian / U.second ** 3, TensorKind.SCALAR):
            raise ValueError()
        if not bstar.check(1 / U.meter, TensorKind.SCALAR):
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
        dt = (self.mean_anomaly / self.mean_motion).secure(U.second, TensorKind.SCALAR)
        self.time_at_periaster = self.epoch - dt

    @cached_property
    def compute_tuple(self) -> OrbitalElementsComputeTuple:
        """Return orbital elements in scalar units expected by SGP4 kernels."""
        return OrbitalElementsComputeTuple(
            n=self.mean_motion.scalar.value("radian/minute"),
            i=self.inclination.scalar.value("radian"),
            e=self.eccentricity.scalar.value("dimensionless"),
            argp=self.arg_of_pericenter.scalar.value("radian"),
            raan=self.ra_of_asc_node.scalar.value("radian"),
            M=self.mean_anomaly.scalar.value("radian"),
            bstar=self.bstar.scalar.value("1/meter") * 6_378_137,
        )

    @property
    def period(self) -> Annotated[Tensor, TensorBound(units=U.second, kind=TensorKind.SCALAR)]:
        """Orbital period derived from the mean motion."""
        from math import tau
        return ((scalar(tau).with_units(U.radian)) / self.mean_motion).secure(U.second, TensorKind.SCALAR)

    def to_coordinates(self) -> Coordinates:
        """Convert orbital elements to GCRF position and velocity vectors."""
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
        """Estimate orbital elements from GCRF state vectors at ``epoch``."""
        els = gcrf_state_vectors2elements(pos, vel)
        return OrbitalElements(
            epoch=epoch,
            eccentricity=els.eccentricity,
            inclination=els.inclination,
            ra_of_asc_node=els.ra_of_asc_node,
            arg_of_pericenter=els.arg_of_pericenter,
            mean_motion=els.mean_motion,
            mean_anomaly=els.mean_anomaly,
            mean_motion_dot=scalar(0).with_units(U.radian / U.second ** 2),
            mean_motion_ddot=scalar(0).with_units(U.radian / U.second ** 3),
            bstar=scalar(0).with_units(1 / U.meter),
        )

    @staticmethod
    def from_celestrak_norad_cat_id(catnr: int, log: bool = False) -> "OrbitalElements":
        """Fetch latest GP data from Celestrak and convert to orbital elements."""
        from leorbit.ext import get_celestrak_gpdata

        return get_celestrak_gpdata(catnr, log).to_orbital_elements()


class Interpolation(Enum):
    """Interpolation modes used when sampling trajectories."""

    CONSTANT = "constant"
    LINEAR = "linear"


class PosVelArray(NamedTuple):
    """Array pair of position vectors and optional velocity vectors."""

    pos: PosVecArray
    vel: VelVecArray | None


class Trajectory:
    """U.second-sampled coordinates over an interval with lazy frame conversions."""

    def __init__(self, interval: TimeInterval, frame: Frame, pos: PosVecArray, vel: VelVecArray | None = None):
        """Create a trajectory sampled on ``interval`` in ``frame``."""
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
            # Relative frame transforms are affine (include translation), which
            # applies to positions but not to velocities.
            if isinstance(self.privileged_frame, AbsoluteFrame) and isinstance(frame, AbsoluteFrame):
                v = transform.do(v)
            else:
                v = None

        self.positions[frame] = PosVelArray(p, v)

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame

    @lru_cache(4096)
    def coordinates_at(self, epoch: Timestamp, interpolation: Interpolation = Interpolation.CONSTANT) -> Coordinates:
        """Return coordinates snapshot at ``epoch`` from the sampled trajectory."""
        if interpolation != Interpolation.CONSTANT:
            raise NotImplementedError("Only `CONSTANT` interpolation is implemented for now.")
        if epoch not in self.interval:
            raise ValueError("Epoch is out of bounds of this trajectory.")

        idx = self.interval._time2idx(epoch)
        pos, vel = self.positions[self.privileged_frame]
        return Coordinates(epoch, self.privileged_frame, pos[idx], vel[idx] if vel is not None else None)

    @lru_cache(4096)
    def get_pos(self, epoch: Timestamp, frame: Frame, interpolation: Interpolation = Interpolation.CONSTANT) -> PosVec:
        """Return trajectory position at ``epoch`` expressed in ``frame``."""
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
        """Return trajectory velocity at ``epoch`` expressed in ``frame``."""
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
        """Return the full sampled position array in ``frame``."""
        self._compute_new_frame(frame)
        pos, _ = self.positions[frame]
        return pos

    @lru_cache(16)
    def trajectory_vel(self, frame: Frame) -> VelVecArray:
        """Return the full sampled velocity array in ``frame``."""
        if not self.vel_available:
            raise ValueError("Velocity data is not available for this trajectory.")

        self._compute_new_frame(frame)
        _, vel = self.positions[frame]
        return cast(Tensor, vel)

    @lru_cache(16)
    def gps(self) -> GPSTrajectory:
        """Return geodetic arrays (lon/lat/alt) for the full trajectory."""
        return itrf2gps(self.trajectory_pos(AbsoluteFrame.ITRF))

    @lru_cache(4096)
    def gps_at(self, epoch: Timestamp, interpolation: Interpolation = Interpolation.CONSTANT) -> GPS:
        """Return geodetic coordinates at ``epoch``."""
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

    @lru_cache(16)
    def horizontal(self, local_frame: EarthLocalFrame) -> HorizontalTrajectory:
        """Return horizontal azimuth/altitude/range arrays for ``local_frame``."""
        return itrf2horizontal(self.trajectory_pos(AbsoluteFrame.ITRF), local_frame)

    @lru_cache(4096)
    def horizontal_at(
        self,
        epoch: Timestamp,
        local_frame: EarthLocalFrame,
        interpolation: Interpolation = Interpolation.CONSTANT,
    ) -> Horizontal:
        """Return horizontal coordinates at ``epoch`` for ``local_frame``."""
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
