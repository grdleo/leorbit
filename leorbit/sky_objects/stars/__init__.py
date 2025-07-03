from abc import ABC, abstractmethod
from functools import lru_cache
from math import atan, cos, sin, sqrt, tan
from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from events.timeline import CoordinatesTimeline
from leorbit.coordinates.trajectory import Trajectory
from leorbit.propagators.no import NoPropagator
from physics.constants import MU_EARTH
from physics.time import Time

from pint import Quantity as Q_

from physics.time_interval import TimeInterval
from sky_objects import SkyObject

class Star(SkyObject, ABC):
    """Astronomical object (planets, moons, stars, etc) evolving in the sky. Those objects do not have to be revolving Earth per se.
    
    Yet their osculating orbit elements are known at every epoch, and thanks to them one can compute its coordinates at a given epoch.
    
    So in order to create a `Star` one has to manually implement the value of their orbital elements.
    """
    def __init__(self, name: str = "No name star", radius: Q_ = Q_("1e7 m"), mass: Q_ = Q_("1e25 kg")):
        super().__init__(name, radius, mass)
        
    def coordinates(self, at: Time) -> Coordinates:
        return self.earth_osculating_orbit(at).to_coordinates()
    
    def trajectory(self, during: TimeInterval) -> Trajectory:
        els = self.earth_osculating_orbit(during.start)
        propagator = NoPropagator(els)
        return propagator.propagate_timeline(during)

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

    @abstractmethod
    def semi_major_axis(self, epoch: Time) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def a(self, epoch: Time) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return self.semi_major_axis(epoch)
    
    @abstractmethod
    def eccentricity(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def e(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        return self.eccentricity(epoch)
    
    @abstractmethod
    def inclination(self, epoch: Time) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def i(self, epoch: Time) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return self.inclination(epoch)
    
    @abstractmethod
    def ra_of_asc_node(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def Ω(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return self.ra_of_asc_node(epoch)
    
    @abstractmethod
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def ω(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        return self.arg_of_pericenter(epoch)
    
    @abstractmethod
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