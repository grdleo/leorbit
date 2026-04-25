from abc import ABC, abstractmethod
from functools import lru_cache
from leorbit.coordinates import Coordinates, OrbitalElements, Trajectory
from typing import Annotated

from leorbit.mathematics import Angle, AngularVelocity, Dimless, Length, Mass, Quantity, Tensor, TensorBound, TensorKind, atan, cos, normalize_angle, scalar, sin, tan
from leorbit.propagator import NoPropagator, Propagator
from leorbit.time import TimeInterval, Timestamp
from leorbit.utils import semi_major_axis_earth_to_mean_motion


class SkyObject(ABC):
    @abstractmethod
    def coordinates(self, at: Timestamp) -> Coordinates:
        ...
    
    @abstractmethod
    def trajectory(self, during: TimeInterval) -> Trajectory:
        ...

class Satellite(SkyObject):
    """Small `SkyObject` with GP data to be propagated by a `Propagator`"""
    def __init__(self, 
        name: str, 
        orbital_elements: OrbitalElements, 
        propagator: type[Propagator] = NoPropagator
    ):
        self.name = name
        self.propagator = propagator(orbital_elements)
    
    def coordinates(self, at: Timestamp) -> Coordinates:
        return self.propagator.propagate(at)
    
    def trajectory(self, during: TimeInterval) -> Trajectory:
        return self.propagator.propagate(during)
    
class Body(SkyObject, ABC):
    """Astronomical object (planets, moons, stars, etc) evolving in the sky. Those objects do not have to be revolving Earth per se.
    
    Yet their osculating orbit elements are known at every epoch, and thanks to them one can compute its coordinates at a given epoch.
    
    So in order to create a `Star` one has to manually implement the value of their orbital elements.
    """
    def __init__(self, 
        name: str, 
        radius: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)], 
        mass: Annotated[Tensor, TensorBound(dimension=Mass, kind=TensorKind.SCALAR)]
    ):
        self.name = name
        self.body_radius = radius
        self.body_mass = mass
    
    @lru_cache(512)
    def coordinates(self, at: Timestamp) -> Coordinates:
        return self.earth_osculating_orbit(at).to_coordinates()
    
    @lru_cache(16)
    def trajectory(self, during: TimeInterval) -> Trajectory:
        els = self.earth_osculating_orbit(during.start)
        propagator = NoPropagator(els)
        return propagator.propagate(during)

    def earth_osculating_orbit(self, epoch: Timestamp) -> OrbitalElements:
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
    def semi_major_axis(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def a(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        return self.semi_major_axis(epoch)
    
    @abstractmethod
    def eccentricity(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def e(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        return self.eccentricity(epoch)
    
    @abstractmethod
    def inclination(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def i(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        return self.inclination(epoch)
    
    @abstractmethod
    def ra_of_asc_node(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def Ω(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return self.ra_of_asc_node(epoch)
    
    @abstractmethod
    def arg_of_pericenter(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def ω(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        return self.arg_of_pericenter(epoch)
    
    @abstractmethod
    def mean_anomaly(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def M(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        return self.mean_anomaly(epoch)
    
    def mean_motion(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=AngularVelocity, kind=TensorKind.SCALAR)]:
        """Returns the mean motion (n) of this star at given epoch."""
        return semi_major_axis_earth_to_mean_motion(self.semi_major_axis(epoch))
    def n(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=AngularVelocity, kind=TensorKind.SCALAR)]:
        """Returns the mean motion (n) of this star at given epoch."""
        return self.mean_motion(epoch)
    
    def radius(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        a = self.semi_major_axis(epoch)
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        return a * (1 - e * cos(E))
    def r(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        return self.radius(epoch)

    def true_anomaly(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the true anomaly (υ) of this star at given epoch."""
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        tan_half_nu = ((1 + e) / (1 - e)) ** 0.5 * tan(.5 * E)
        nu = 2 * atan(tan_half_nu)
        return nu
    def υ(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the true anomaly (υ) of this star at given epoch."""
        return self.true_anomaly(epoch)
    
    def eccentric_anomaly(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        M = self.mean_anomaly(epoch)
        e = self.eccentricity(epoch)
        E = M.copy()
        for _ in range(5):  # compute excentric anomaly
            E = e * sin(E) + M
        return E
    def E(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        return self.eccentric_anomaly(epoch)
    

class Moon(Body):
    """`Body` object representing the Moon.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """
    def __init__(self):
        super().__init__("Moon", 1737.4 * Quantity.kilo_meter, 7.342e22 * Quantity.kilo_gram)

    def semi_major_axis(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        return 60.2666 * Quantity.radii_earth
    
    def inclination(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        return 0.08980417133211624 * Quantity.radian
    
    def eccentricity(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        return scalar(0.054900)
    
    def arg_of_pericenter(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t = epoch.from_mil
        argp = 318.0634 + 0.1643573223 * t.scalar.value("day")
        return normalize_angle(argp * Quantity.degree)
    
    def ra_of_asc_node(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        t = epoch.from_mil
        raan = 125.1228 - 0.0529538083 * t.scalar.value("day")
        return normalize_angle(raan * Quantity.degree)
    
    def mean_anomaly(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t = epoch.from_mil
        M = 115.3654 + 13.0649929509 * t.scalar.value("day")
        return normalize_angle(M * Quantity.degree)
    
class Sun(Body):
    """`Body` object representing the Sun.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """

    def __init__(self):
        super().__init__("Sun", 696_340 * Quantity.kilo_meter, 1.988e30 * Quantity.kilo_gram)
    
    def semi_major_axis(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        return 149_597_870_700 * Quantity.meter
    
    def inclination(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        return 0 * Quantity.radian
    
    def eccentricity(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Dimless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        t = epoch.from_mil
        e = 0.016709 - 1.151e-9 * t.scalar.value("day")
        return scalar(e)
    
    def arg_of_pericenter(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t = epoch.from_mil
        argp = 282.9404 + 4.70935e-5 * t.scalar.value("day")
        return normalize_angle(argp * Quantity.degree)
    
    def ra_of_asc_node(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return 0 * Quantity.radian
    
    def mean_anomaly(self, epoch: Timestamp) -> Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t = epoch.from_mil
        M = 356.0470 + 0.9856002585 * t.scalar.value("day")
        return normalize_angle(M * Quantity.degree)