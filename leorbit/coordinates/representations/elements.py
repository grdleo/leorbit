"""Implementation of "orbital elements" of an object orbiting Earth."""

from collections import namedtuple
import json
from math import atan, cos, sin, sqrt, tan, atan2, tau
from typing import Optional, Self
from dataclasses import dataclass, field
from typing import TypedDict, NamedTuple

import numpy as np
from pint import Quantity

from coordinates.representations import CoordinatesRepresentation

from coordinates.coordinates import Coordinates
from frames.absolute_frame import AbsoluteFrame
from leorbit.algorithms.elements2orthogonal_gcrf import elements2orthogonal_gcrf
from leorbit.algorithms.elements2orthogonal_gcrf import orb_tels2posvel_gcrf, orb_tels2posvel_gcrf_numpy
from leorbit.algorithms.utils import mean2eccentric_anomaly
from mathematics.vec3 import Vec3
from mathematics.units import UREG
from ext.celestrak import get_celestrak_gpdata_json
from physics.constants import MU_EARTH, SQRT_MU_EARTH
from physics.time import Time

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
    eccentricity: Quantity # [1]
    inclination: Quantity # [rad]
    ra_of_asc_node: Quantity # [rad]
    arg_of_pericenter: Quantity # [rad]
    mean_motion: Quantity # [rad]
    mean_anomaly: Quantity # [rad]
    mean_motion_dot: Quantity = field(init=True, default_factory=lambda: UREG("0 rad/s**2")) # [rad/s²]
    mean_motion_ddot: Quantity = field(init=True, default_factory=lambda: UREG("0 rad/s**3")) # [rad/s3]
    bstar: Quantity = field(init=True, default_factory=lambda: UREG("0 1/m"))

    name: str = "No name"
    norad_cat_id: Optional[int] = None

    eccentric_anomaly: Quantity = field(init=False) # [rad]
    true_anomaly: Quantity = field(init=False) # [rad]
    semi_major_axis: Quantity = field(init=False) # [m]
    semi_minor_axis: Quantity = field(init=False) # [m]
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

    @classmethod
    def from_state_vectors(cls: "OrbitalElements", coordinate: Coordinates) -> "OrbitalElements":
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

    @classmethod
    def from_celestrak_json(cls: "OrbitalElements", celestrak_json: str | dict[str, str | float | int]) -> "OrbitalElements":
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

    @classmethod
    def from_celestrak_norad_cat_id(cls, catnr: int, log: bool = False) -> "OrbitalElements":
        """
        Fetches lastest GP data on `celestrak.org` corresponding to given NORAD catalog ID, 
        and creates an instance of `OrbitalElements` using those.

        If last fetch on Celestrak is recent enough, uses cached GP data.
        """
        gp_dict = get_celestrak_gpdata_json(catnr, log)
        return OrbitalElements.from_celestrak_json(gp_dict)