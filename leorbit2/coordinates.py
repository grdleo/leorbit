from dataclasses import dataclass
from functools import lru_cache
from typing import NamedTuple, Self

from enum import Enum
from functools import lru_cache
from typing import ParamSpec, Callable, TypeVar, cast

from leorbit2.frames import AbsoluteFrame, frame_transform_factory
from leorbit2.mathematics import D, Dim, Scalar, TransformChain, Vector3, Transform, TransformVector3RotationZ, TransformIdentify, D.Velocity
from leorbit2.time import Time
from leorbit.frames import Frame

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
        if frame in self.positions:
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
    def from_gps(longitude: Q_, latitude: Q_, altitude: Q_ = 0 * POS_UNIT, epoch: Time = None) -> Self:
        assert longitude.check(UREG.radians)
        assert latitude.check(UREG.radians)
        assert altitude.check(UREG.meters)

        theta = longitude
        delta = latitude
        rho = geocentric_radius_earth(delta) + altitude
        pos = Vec3.from_spherical(theta, delta, rho)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, AbsoluteFrame.ITRF, pos)
    
    @lru_cache
    def gps(self) -> GPS:
        """Returns this coordinates as their GPS representation"""
        gps = self._already_computed_repr.get(GPS, None)
        if gps is not None:
            return gps
        
        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        lon = itrf_pos.theta
        lat = itrf_pos.delta
        alt = itrf_pos.rho - geocentric_radius_earth(lat)

        return GPS(
            longitude=lon, 
            latitude=lat, 
            altitude=alt
        )
    
    ### ### ###

    ### HORIZONTAL ###
    ### ########## ###

    @staticmethod
    def from_horizontal(azimuth: Q_, altitude: Q_, distance: Q_, frame: EarthLocalFrame, epoch: Time = None) -> Self:
        assert isinstance(frame, EarthLocalFrame)
        
        pos_local = Vec3.from_spherical(azimuth, altitude, distance)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, frame, pos_local)
    
    @lru_cache
    def horizontal(self, local_frame: EarthLocalFrame) -> Horizontal:
        """Returns the representation of this coordinates as 'horizontal', in the given frame"""
        assert local_frame.reference_frame == AbsoluteFrame.ITRF # FIXME

        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        horizontal_pos = local_frame.transform.apply(itrf_pos)

        return Horizontal(
            azimuth=horizontal_pos.theta,
            altitude=horizontal_pos.delta,
            distance=itrf_pos.rho
        )
    
    ### ########## ###
    
    def __eq__(self, o: Self) -> bool:
        common_computed_frames = set(self.positions.keys()).intersection(o.positions.keys())

        if common_computed_frames:
            f: Frame
            for f in common_computed_frames: break # Get any element from set
            p, v = self.positions[f]
            po, vo = o.positions[f]

            return (p == po) and (v == vo)
        
        # We get a frame already computed in `self` and compute it for the other
        o._compute_new_frame(self.privileged_frame)
        # Now they have a common frame. We can call `__eq__` again
        return self.__eq__(o)
    
    def __neq__(self, o: Self) -> bool:
        return not self.__eq__(o)
    
#################################
# REPRESENTATIONS

    
class CoordinatesRepresentation:
    pass

"""Implementation of "orbital elements" of an object orbiting Earth."""

from collections import namedtuple
import json
from math import atan, cos, sin, sqrt, tan, atan2, tau
from typing import Optional, Self
from dataclasses import dataclass, field
from typing import TypedDict, NamedTuple

import numpy as np

def angle_between(angle: Quantity, mini_deg: float, maxi_deg: float) -> bool:
    return (mini_deg * UREG.degree) <= angle <= (maxi_deg * UREG.degree)

class OrbitalElementsComputeTuple(NamedTuple):
    n: float # [rad/min]
    i: float # [rad]
    e: float # [1]
    argp: float # [rad]
    raan: float # [rad]
    M: float # [rad]
    bstar: float # [1/earthRadii]

@dataclass(frozen=True)
class OrbitalElements(CoordinatesRepresentation):
    """Dataclass holding orbital elements, at a given epoch, gathered from Celestrak.org 
    (also known as GP data)"""
    epoch: Time
    eccentricity: Scalar[D.Dimless] # [1]
    inclination: Scalar[D.Angle] # [rad]
    ra_of_asc_node: Scalar[D.Angle] # [rad]
    arg_of_pericenter: Scalar[D.Angle] # [rad]
    mean_motion: Scalar[D.Angle] # [rad]
    mean_anomaly: Scalar[D.Angle] # [rad]
    mean_motion_dot: Scalar[D.AngularAcc] = field(init=True, default_factory=lambda: Scalar[D.AngularAcc].new(0)) # [rad/s²]
    mean_motion_ddot: Scalar[D.AngularJerk] = field(init=True, default_factory=lambda: Scalar[D.AngularJerk].new(0)) # [rad/s3]
    bstar: Scalar[D.InvLength] = field(init=True, default_factory=lambda: Scalar[D.InvLength].new(0)) # [1/m]

    name: str = "No name"
    norad_cat_id: Optional[int] = None

    # \/ CACHED PROPERTIES \/
    eccentric_anomaly: Scalar[D.Length] = field(init=False) # [rad]
    true_anomaly: Scalar[D.Angle] = field(init=False) # [rad]
    semi_major_axis: Scalar[D.Length] = field(init=False) # [m]
    semi_minor_axis: Scalar[D.Length] = field(init=False) # [m]
    time_at_periaster: Time = field(init=False)

    def __post_init__(self):
        assert self.eccentricity.check("1")
        if not self.eccentricity >= 0:
            raise ValueError()

        assert self.inclination.check("°")
        if not angle_between(self.inclination, 0, 180):
            raise ValueError()

        assert self.ra_of_asc_node.check("°")
        if not angle_between(self.ra_of_asc_node, 0, 360):
            raise ValueError()

        assert self.arg_of_pericenter.check("°")
        if not angle_between(self.arg_of_pericenter, 0, 360):
            raise ValueError()

        assert self.mean_motion.check("rad/s")

        self.mean_anomaly.check("°")
        if not angle_between(self.mean_anomaly, 0, 360):
            raise ValueError()
        
        assert self.mean_motion_dot.check("rad/s**2")
        assert self.mean_motion_ddot.check("rad/s**3")
        assert self.bstar.check("1/m")

        e = self.eccentricity
        M = self.mean_anomaly
        E = mean2eccentric_anomaly(e, M)
        object.__setattr__(self, "eccentric_anomaly", E) # NOTE: to avoid `FrozenInstanceError`...

        tan_half_nu = sqrt((1 + e) / (1 - e)) * tan(.5 * E)
        true_anomaly = 2 * atan(tan_half_nu)
        object.__setattr__(self, "true_anomaly", true_anomaly)

        semi_major_axis = ((MU_EARTH / self.mean_motion**2)**(1/3)).to("m")
        object.__setattr__(self, "semi_major_axis", semi_major_axis)

        time_at_periaster = self.epoch - self.mean_anomaly / self.mean_motion
        object.__setattr__(self, "time_at_periaster", time_at_periaster)

        semi_minor_axis = self.semi_major_axis * sqrt(1 - e**2)
        object.__setattr__(self, "semi_minor_axis", semi_minor_axis)

        _els_as_float_tuple = OrbitalElementsComputeTuple(
            n=self.mean_motion.m_as("rad/min"),
            i=self.inclination.m_as("rad"),
            e=self.eccentricity.m,
            argp=self.arg_of_pericenter.m_as("rad"),
            raan=self.ra_of_asc_node.m_as("rad"),
            M=self.mean_anomaly.m_as("rad"),
            bstar=self.bstar.m_as("1/earthRadii")
        )
        object.__setattr__(self, "_els_as_float_tuple", _els_as_float_tuple)
        self._els_as_float_tuple: OrbitalElementsComputeTuple
    
    @property
    def period(self) -> Quantity:
        """The period of a full revolution."""
        return (tau * UREG("rad")) / self.mean_motion
    
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