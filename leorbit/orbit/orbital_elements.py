"""Implementation of "orbital elements" of an object orbiting Earth."""

import json
from math import atan, cos, sin, sqrt, tan, atan2, tau
from typing import Optional, Self
from dataclasses import dataclass, field
from typing import TypedDict

from leorbit.math.coordinate import Coordinates, AbsoluteFrame
from leorbit.math.vector import Vec3
from leorbit.ext.celestrak import get_celestrak_gpdata_json
from leorbit.orbit.constants import MU_EARTH, SQRT_MU_EARTH, Q_
from leorbit.time import Time

class _FastComputeDictOrbitalElements(TypedDict):
    """Storing orbital with proper units, ready to be computed with SGP algo"""
    n: float # [rad/min]
    i: float # [rad]
    e: float # [1]
    argp: float # [rad]
    raan: float # [rad]
    M: float # [rad]
    bstar: float # [1/earthRadii]

@dataclass(frozen=True)
class OrbitalElements:
    """Dataclass holding orbital elements, at a given epoch, gathered from Celestrak.org 
    (also known as GP data)"""
    epoch: Time
    eccentricity: Q_ # [1]
    inclination: Q_ # [rad]
    ra_of_asc_node: Q_ # [rad]
    arg_of_pericenter: Q_ # [rad]
    mean_motion: Q_ # [rad]
    mean_anomaly: Q_ # [rad]
    mean_motion_dot: Q_ = field(init=True, default_factory=lambda: Q_("0 rad/s**2")) # [rad/s²]
    mean_motion_ddot: Q_ = field(init=True, default_factory=lambda: Q_("0 rad/s**3")) # [rad/s3]
    bstar: Q_ = field(init=True, default_factory=lambda: Q_("0 1/m"))

    name: str = "No name"
    norad_cat_id: Optional[int] = None

    eccentric_anomaly: Q_ = field(init=False) # [rad]
    true_anomaly: Q_ = field(init=False) # [rad]
    semi_major_axis: Q_ = field(init=False) # [m]
    semi_minor_axis: Q_ = field(init=False) # [m]
    time_at_periaster: Time = field(init=False)

    def __post_init__(self):
        assert self.eccentricity.check("1")
        if not self.eccentricity >= 0:
            raise ValueError()

        assert self.inclination.check("°")
        if not Q_(0, "°") <= self.inclination <= Q_(180, "°"):
            raise ValueError()

        assert self.ra_of_asc_node.check("°")
        if not Q_(0, "°") <= self.ra_of_asc_node <= Q_(360, "°"):
            raise ValueError()

        assert self.arg_of_pericenter.check("°")
        if not Q_(0, "°") <= self.arg_of_pericenter <= Q_(360, "°"):
            raise ValueError()

        assert self.mean_motion.check("rad/s")

        self.mean_anomaly.check("°")
        if not Q_(0, "°") <= self.mean_anomaly <= Q_(360, "°"):
            raise ValueError()
        
        assert self.mean_motion_dot.check("rad/s**2")
        assert self.mean_motion_ddot.check("rad/s**3")
        assert self.bstar.check("1/m")

        e = self.eccentricity
        M = self.mean_anomaly
        E = M
        for _ in range(5):
            E = M + e * sin(E)
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

        _els_ready2compute = _FastComputeDictOrbitalElements(
            n=self.mean_motion.m_as("rad/min"),
            i=self.inclination.m_as("rad"),
            e=self.eccentricity.m,
            argp=self.arg_of_pericenter.m_as("rad"),
            raan=self.ra_of_asc_node.m_as("rad"),
            M=self.mean_anomaly.m_as("rad"),
            bstar=self.bstar.m_as("1/earthRadii")
        )
        object.__setattr__(self, "_els_ready2compute", _els_ready2compute)
        self._els_ready2compute: _FastComputeDictOrbitalElements
    
    @property
    def period(self) -> Q_:
        """The period of a full revolution."""
        return Q_(tau, "rad") / self.mean_motion
    
    def to_coordinates(self) -> Coordinates:
        """Returns current orbital elements, at given epoch, 
        as a `Coordinates` object."""
        # position of satellite in orbit plane (with z = 0)
        υ = self.true_anomaly
        e = self.eccentricity
        a = self.semi_major_axis
        Ω = self.ra_of_asc_node
        ω = self.arg_of_pericenter
        i = self.inclination
        E = self.eccentric_anomaly

        ee = e**2
        one_ee = (1 - ee)
        esinE = e * sin(E)
        r = a * one_ee / (1 + e * cos(υ))
        rd = SQRT_MU_EARTH * a**.5 * esinE / r
        rυd = rd * one_ee / esinE

        c_raan, s_raan = cos(Ω), sin(Ω)
        c_i, s_i = cos(i), sin(i)
        υpω = υ + ω
        c_theta, s_theta = cos(υpω), sin(υpω)

        def unitvec_gcrf(x: Q_, y: Q_) -> Vec3:
            return Vec3(
                c_raan * x - s_raan * c_i * y,
                s_raan * x + c_raan * c_i * y,
                s_i * y
            )

        ur = unitvec_gcrf(c_theta, s_theta)
        ut = unitvec_gcrf(-s_theta, c_theta)
        pos_gcrf = ur * r
        vel_gcrf = ur * rd + ut * rυd

        return Coordinates(AbsoluteFrame.GCRF, pos_gcrf, vel_gcrf, self.epoch)

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
        argp = Q_(atan2(ye_asc, xe_asc), "rad")

        xp_asc = xaxis_asc.dot(pos).m_as("m")
        yp_asc = yaxis_asc.dot(pos).m_as("m")
        nu = (Q_(atan2(yp_asc, xp_asc), "rad") - argp) % tau

        e = abs(ecc_vec) # [1]
        ee = e * e
        eee = ee * e
        eeee = eee * e
        i = north.angle(kinetic) # [rad]
        raan = atan2(asc.y, asc.x) % tau # [rad]
        raan = Q_(raan, "rad") # convert to Quantity
        a = kinetic_sq / (MU_EARTH * (1 - ee)) # [m]
        aaa = a**3
        n = ((MU_EARTH / aaa)**0.5).to("rad/s")  # [rad/s]
        M = (
            nu
            - 2 * e * sin(nu)
            + (3 / 4 * ee + 1 / 8 * eeee) * sin(2 * nu)
            - 1 / 3 * eee * sin(3 * nu)
            + 5 / 32 * eeee * sin(4 * nu)
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
                eccentricity=Q_(query["ECCENTRICITY"], "dimensionless"),
                inclination=Q_(query["INCLINATION"], "°"),
                ra_of_asc_node=Q_(query["RA_OF_ASC_NODE"], "°"),
                arg_of_pericenter=Q_(query["ARG_OF_PERICENTER"], "°"),
                mean_motion=Q_(query["MEAN_MOTION"], "turn/day"),
                mean_anomaly=Q_(query["MEAN_ANOMALY"], "°"),
                mean_motion_dot=Q_(query["MEAN_MOTION_DOT"], "turn/day^2") * 2, # NOTE: factor is cancelled when loading directly from Celestrack
                mean_motion_ddot=Q_(query["MEAN_MOTION_DDOT"], "turn/day^3") * 6, # NOTE: factor is cancelled when loading directly from Celestrack
                bstar=Q_(query["BSTAR"], "1/earthRadii"),
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
    
    def _to_tle( # TODO
        self,
        obj_id="0000-000X",
        classi="U",
        rev_at_epoch="99999",
        name="NONAME",
        norad_cat_id="00000",
        el_set_no="999",
    ) -> str:
        """Returns orbital elements as a TLE"""
        n_dot_tle = ("-" if self.mean_motion_dot < 0 else " ") + str((self.mean_motion_dot.to("turn/day^2") / 2).m)[
            1:
        ].ljust(9, "0")
        n_dot_tle = n_dot_tle[:10]

        n_ddot_tle = ("-" if self.mean_motion_ddot < 0 else " ") + "{:.5e}".format(
            abs((self.mean_motion_ddot.to("turn/day^3") / 6).m)
        ).replace("e", "")[2:]
        n_ddot_tle = n_ddot_tle[:8]

        bstar_tle = ("-" if self.bstar < 0 else " ") + "{:.5e}".format(abs(self.bstar.to("1/earthRadii"))).replace("e", "")[2:]
        bstar_tle = bstar_tle[:8]

        def checksum(st: str) -> int:
            somme = 0
            for s in st.replace("-", "1"):
                try:
                    somme += int(s)
                except:
                    pass
            return str(somme % 10)

        def nb_format(nb: float, lsize: int, rsize: int) -> str:
            l = int(nb)
            r = str(nb - l)[1:]
            return (
                f"{str(l)[:lsize].rjust(lsize, '0')}.{str(r)[:rsize].ljust(rsize, '0')}"
            )

        one = f"1 {norad_cat_id[:5]}{classi[:1]} {obj_id[2:].replace('-', '')[:6]}   {self.epoch.year_day} {n_dot_tle} {n_ddot_tle} {bstar_tle} 0  {el_set_no[:3]}"
        one += checksum(one)

        i_tle = nb_format(self.inclination.to("deg").m, 3, 4)
        raan_tle = nb_format(self.ra_of_asc_node.to("deg").m, 3, 4)
        argp_tle = nb_format(self.arg_of_pericenter.to("deg").m, 3, 4)
        M0_tle = nb_format(self.mean_anomaly.to("deg").m, 3, 4)
        n0_tle = nb_format(self.mean_motion.to("deg").m, 2, 8)
        e_tle = str(self.eccentricity.m)[2:][:7].ljust(7)

        two = f"2 {norad_cat_id[:5]} {i_tle} {raan_tle} {e_tle} {argp_tle} {M0_tle} {n0_tle}{rev_at_epoch[5:]}"
        two += checksum(two)

        return f"{name.upper()}\n{one}\n{two}"
    
    @classmethod
    def _from_tle(cls: "OrbitalElements", tle: str) -> Self: # TODO
        """Creates an instance of Orbit object, using TLE (Two-Line Elements) string.

        String may contain a title line. If not, this is not a problem, a name will be added automatically.
        Conventions used can be found here: https://en.wikipedia.org/wiki/Two-line_element_set

        Parameters
        ----------
        tle : str
            TLE for the orbit, at given time. Three lines separated with '\n' (slash n)

        """
        raise NotImplementedError("Needs to be rewritten")
    
        if tle.endswith("\n"):  # if ends with break line, remove the last one
            tle = tle[0 : len(tle) - 2]

        tle = tle.replace("\r", "")
        splited = tle.split("\n")

        name = "noname"
        first = ""
        second = ""

        l = 1 if len(splited) > 2 else 0 # first line may be name
        first = splited[l]
        second = splited[l+1]

        y = int(first[18:20])
        dec_day = float(first[20:32])

        epoch = Time.fromisoformat(f"20{y}-01-01T00:00:00") + Q_(dec_day - 1, "day")

        inclination = Q_(float(second[8:15]), "°")
        raan = Q_(float(second[17:24]), "°")
        eccentricity = Q_(float(f"0.{second[26:33]}"), "dimensionless")
        argp = Q_(float(second[34:42]), "°")
        mean_anomaly = Q_(float(second[43:51]), "°")
        mean_motion = Q_(float(second[52:65]), "turn/day")
        mean_motion_dot = Q_(2 * float(first[33:42]), "turn/day^2") # NOTE: factor is cancelled when loading directly from Celestrack
        mmtay = first[45:52]
        str_mean_motion_taylor = f"{mmtay[:5]}e{mmtay[5:]}"
        mean_motion_ddot = Q_(6 * float(str_mean_motion_taylor), "turn/day^3") # NOTE: factor is cancelled when loading directly from Celestrack
        bstar_tle = first[54:61]
        bstar_str = f"{bstar_tle[:6]}e{bstar_tle[6:]}"
        bstar = Q_(float(bstar_str), "1/earthRadii")

        return OrbitalElements(
            epoch,
            eccentricity,
            inclination,
            raan,
            argp,
            mean_motion,
            mean_anomaly,
            mean_motion_dot,
            mean_motion_ddot,
            bstar
        )