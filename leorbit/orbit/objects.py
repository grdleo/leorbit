"""Astronomical objects definitions (Satellite, Sun, Moon, etc.)"""

__pdoc__ = dict()

from abc import ABC, abstractmethod
from functools import lru_cache
from math import atan, cos, sin, sqrt, tan
from typing import Optional, overload
from leorbit.ext.heavens_above import get_standard_magnitude
from leorbit.math import FULL_REV
from leorbit.math.coordinate import Coordinates, PosVelGCRF
from leorbit.orbit.constants import MU_EARTH, Q_, RADII_EARTH
from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.orbit.propagator import Propagator
from leorbit.time import Time
from leorbit.time.timeline import Timeline

MAX_LRU_CACHE = 10_000_000
        
class AstroObject(ABC):
    """Astronomical object whose position can be computed at any epoch."""
    def __init__(self, name: str, radius: Q_, mass: Q_):
        try:
            assert radius > Q_("0m")
            assert mass > Q_("0g")
            assert isinstance(name, str)
        except Exception as ex:
            raise ex
        
        self.radius = radius
        self.mass = mass
        self.name = name

    def coordinates(self, epoch: Time) -> Coordinates:
        """Returns the position and velocity as a `Coordinates` object, 
        at given epoch."""
        return Coordinates.from_gcrf_pos_vel_tuple(
            self._gcrf_pos_vel_tuple(epoch),
            epoch
        )

    def simulation(self, timeline: Timeline) -> list[Coordinates]:
        """Returns the positions and velocities as a list of `Coordinates`,
        for every instant in the given `Timeline`"""
        return [
            Coordinates.from_gcrf_pos_vel_tuple(
                self._gcrf_pos_vel_tuple(e),
                e
            )
            for e in timeline
        ]

    def earth_osculating_orbit(self, epoch: Time) -> OrbitalElements:
        """Returns the [osculating orbit](https://en.wikipedia.org/wiki/Osculating_orbit) 
        of this `AstroObject` at given epoch as a `OrbitalElements` class."""
        raise NotImplementedError()

    def _gcrf_pos_vel_tuple(self, epoch: Time) -> PosVelGCRF:
        raise NotImplementedError()

class Satellite(AstroObject):
    """Astronomical object orbiting Earth. Its dynamic is computed using orbital elements,
    whose are subject to evolution due to perturbation (drag, gravitational, etc).
    Coordinates at given epoch are computed with the given propagator. Osculating orbit can be computed thanks to the coordinates."""
    def __init__(self, propagator: Propagator, name: str = "", catnr: Optional[int] = None, radius: Q_ = Q_("1m"), mass: Q_ = Q_("1000kg")):
        super().__init__(name, radius, mass)
        self.catnr = catnr
        self.propagator = propagator
    
    def __repr__(self) -> str:
        name = f"name='{self.name}' " if self.name else ""
        return f"<Satellite {name}t₀='{self.propagator._start_elements.epoch.human}' [{self.propagator.__class__.__name__}]>"
    
    def earth_osculating_orbit(self, epoch: Time) -> OrbitalElements:
        return OrbitalElements.from_state_vectors(
            self.coordinates(epoch)
        )
    
    @lru_cache(MAX_LRU_CACHE)
    def _gcrf_pos_vel_tuple(self, epoch: Time) -> PosVelGCRF:
        return self.propagator._fast_propagate_gcrf(epoch)
    
    @property
    def std_mag(self) -> float | None:
        """Standard magnitude of this satellite."""
        if not hasattr(self, "_std_mag"):
            self._std_mag = get_standard_magnitude(self.catnr) if isinstance(self.catnr, int) else None
        return self._std_mag
    
    @std_mag.setter
    def std_mag(self, value: float | None):
        if isinstance(value, float) or value is None:
            self._std_mag = value
            return
        raise ValueError()
    
    def magnitude(self, epoch: Time) -> float:
        """.. todo::
            Not implemented yet.
        
        Apparent magnitude of this satellite at given epoch"""
        raise NotImplementedError()
        
class Star(AstroObject, ABC):
    """Astronomical object (planets, moons, stars, etc) evolving in the sky. Those objects do not have to be revolving Earth per se.
    
    Yet their osculating orbit elements are known at every epoch, and thanks to them one can compute its coordinates at a given epoch.
    
    So in order to create a `Star` one has to manually implement the value of their orbital elements.
    """
    def __init__(self, name: str = "No name star", radius: Q_ = Q_("1e7 m"), mass: Q_ = Q_("1e25 kg")):
        super().__init__(name, radius, mass)

    def earth_osculating_orbit(self, epoch: Time) -> OrbitalElements:
        return OrbitalElements(
            epoch,
            self.eccentricity(epoch),
            self.inclination(epoch),
            self.ra_of_asc_node(epoch),
            self.arg_of_pericenter(epoch),
            self.mean_motion(epoch),
            self.mean_anomaly(epoch)
        )
    
    @lru_cache(MAX_LRU_CACHE)
    def _gcrf_pos_vel_tuple(self, epoch: Time) -> PosVelGCRF:
        oe = self.earth_osculating_orbit(epoch)
        return oe.to_coordinates().to_gcrf_pos_vel_tuple()

    def semi_major_axis(self, epoch: Time) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def a(self, epoch: Time) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return self.semi_major_axis(epoch)
    
    def eccentricity(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def e(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        return self.eccentricity(epoch)
    
    def inclination(self, epoch: Time) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def i(self, epoch: Time) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return self.inclination(epoch)
    
    def ra_of_asc_node(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def Ω(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return self.ra_of_asc_node(epoch)
    
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def ω(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        return self.arg_of_pericenter(epoch)
    
    def mean_anomaly(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def M(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        return self.mean_anomaly(epoch)
    
    def mean_motion(self, epoch: Time) -> Q_:
        """Returns the mean motion (n) of this star at given epoch."""
        return (MU_EARTH / self.semi_major_axis(epoch)**3)**0.5
    def n(self, epoch: Time) -> Q_:
        """Returns the mean motion (n) of this star at given epoch."""
        return self.mean_motion(epoch)
    
    def radius(self, epoch: Time) -> Q_:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        a = self.semi_major_axis(epoch)
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        return a * (1 - e * cos(E))
    def r(self, epoch: Time) -> Q_:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        return self.radius(epoch)

    def true_anomaly(self, epoch: Time) -> Q_:
        """Returns the true anomaly (υ) of this star at given epoch."""
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        tan_half_nu = sqrt((1 + e) / (1 - e)) * tan(.5 * E)
        nu = 2 * atan(tan_half_nu)
        return nu
    def υ(self, epoch: Time) -> Q_:
        """Returns the true anomaly (υ) of this star at given epoch."""
        return self.true_anomaly(epoch)
    
    def eccentric_anomaly(self, epoch: Time) -> Q_:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        M = self.mean_anomaly(epoch)
        e = self.eccentricity(epoch)
        E = Q_(M)
        for _ in range(5):  # compute excentric anomaly
            E = e * sin(E) + M
        return E
    def E(self, epoch: Time) -> Q_:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        return self.eccentric_anomaly(epoch)

class Sun(Star):
    """`Star` object representing the Sun.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """

    def __init__(self):
        super().__init__("Sun", Q_(696_340, "km"), Q_(1.988e30, "kg"))
    
    def semi_major_axis(self, epoch: Time = None) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return Q_(149_597_870_700, "m")
    
    def inclination(self, epoch: Time = None) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return Q_(0, "rad")
    
    def eccentricity(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        e = 0.016709 - 1.151e-9 * t.m_as("day")
        return Q_(e, "dimensionless")
    
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        argp = 282.9404 + 4.70935e-5 * t.m_as("day")
        return Q_(argp, "°") % FULL_REV
    
    def ra_of_asc_node(self, epoch: Time = None) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return Q_(0, "rad")
    
    def mean_anomaly(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        M = 356.0470 + 0.9856002585 * t.m_as("day")
        return Q_(M, "°") % FULL_REV

class Moon(Star):
    """`Star` object representing the Moon.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """
    def __init__(self):
        super().__init__("Moon", Q_(1737.4, "km"), Q_(7.342e22, "kg"))

    def semi_major_axis(self, epoch: Time = None) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return 60.2666 * RADII_EARTH
    
    def inclination(self, epoch: Time = None) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return Q_(0.08980417133211624, "rad")
    
    def eccentricity(self, epoch: Time = None) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        return Q_(0.054900, "dimensionless")
    
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        argp = 318.0634 + 0.1643573223 * t.m_as("day")
        return Q_(argp, "°") % FULL_REV
    
    def ra_of_asc_node(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        raan = 125.1228 - 0.0529538083 * t.m_as("day")
        return Q_(raan, "°") % FULL_REV
    
    def mean_anomaly(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        M = 115.3654 + 13.0649929509 * t.m_as("day")
        return Q_(M, "°") % FULL_REV