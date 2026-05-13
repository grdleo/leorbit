from abc import ABC, abstractmethod
from functools import lru_cache

import numpy as np
from leorbit.coordinates import Coordinates, OrbitalElements, Trajectory
from typing import Annotated

from leorbit.mathematics import U, Tensor, TensorBound, TensorKind, atan, cos, normalize_angle, scalar, sin, tan
from leorbit.propagator import NoPropagator, Propagator
from leorbit.time import TimeInterval, Timestamp
from leorbit.utils import from_mil, semi_major_axis_earth_to_mean_motion


def __one(epoch: Timestamp | TimeInterval) -> Tensor:
    """Transforms epoch into dimentionless scalar tensor with value 1, 
    keeping the shape of the input epoch (scalar or array) in the output."""
    shape = epoch.to_unixepoch().shape
    return Tensor(
        data=np.full(shape, 1.0, dtype=np.float64), 
        units=U.dimensionless
    )

class SkyObject(ABC):
    """Abstract object that can provide coordinates and trajectories."""

    @abstractmethod
    def coordinates(self, at: Timestamp) -> Coordinates:
        """Return object coordinates at instant ``at``."""
        ...
    
    @abstractmethod
    def trajectory(self, during: TimeInterval) -> Trajectory:
        """Return object trajectory sampled over interval ``during``."""
        ...

class Satellite(SkyObject):
    """Small `SkyObject` with GP data to be propagated by a `Propagator`"""
    def __init__(self, 
        name: str, 
        orbital_elements: OrbitalElements, 
        propagator: type[Propagator] = NoPropagator
    ):
        """Create a satellite from orbital elements and a propagation strategy."""
        self.name = name
        self.propagator = propagator(orbital_elements)

    def __repr__(self) -> str:
        epoch_iso = self.propagator.elements.epoch.isoformat
        propagator_name = type(self.propagator).__name__
        return f"<Satellite name='{self.name}' t₀='{epoch_iso}' [{propagator_name}]>"
    
    def coordinates(self, at: Timestamp) -> Coordinates:
        """Propagate and return coordinates at instant ``at``."""
        return self.propagator.propagate(at)
    
    def trajectory(self, during: TimeInterval) -> Trajectory:
        """Propagate and return trajectory over ``during``."""
        return self.propagator.propagate(during)
    
class Body(SkyObject, ABC):
    """Astronomical object (planets, moons, stars, etc) evolving in the sky. Those objects do not have to be revolving Earth per se.
    
    Yet their osculating orbit elements are known at every epoch, and thanks to them one can compute its coordinates at a given epoch.
    
    So in order to create a `Star` one has to manually implement the value of their orbital elements.
    """
    def __init__(self, 
        name: str, 
        radius: Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)], 
        mass: Annotated[Tensor, TensorBound(units=U.kilogram, kind=TensorKind.SCALAR)]
    ):
        """Create a celestial body with physical radius and mass."""
        self.name = name
        self.body_radius = radius.secure(units=U.meter, kind=TensorKind.SCALAR, size=1)
        self.body_mass = mass.secure(units=U.kilogram, kind=TensorKind.SCALAR, size=1)

    def __repr__(self) -> str:
        radius_km = self.body_radius.scalar.value("km")
        mass_kg = self.body_mass.scalar.value("kilo_gram")
        return (
            f"<{self.__class__.__name__} name='{self.name}' "
            f"radius={radius_km:.3f} km mass={mass_kg:.3e} kg>"
        )
    
    @lru_cache(512)
    def coordinates(self, at: Timestamp) -> Coordinates:
        return self.earth_osculating_orbit(at).to_coordinates()
    
    @lru_cache(16)
    def trajectory(self, during: TimeInterval) -> Trajectory:
        els = self.earth_osculating_orbit(during.start)
        propagator = NoPropagator(els)
        return propagator.propagate(during)

    def earth_osculating_orbit(self, epoch: Timestamp) -> OrbitalElements:
        """Build Earth-centered osculating elements at ``epoch`` from body models."""
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
    def semi_major_axis(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def a(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        return self.semi_major_axis(epoch)
    
    @abstractmethod
    def eccentricity(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def e(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        return self.eccentricity(epoch)
    
    @abstractmethod
    def inclination(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def i(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        return self.inclination(epoch)
    
    @abstractmethod
    def ra_of_asc_node(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def Ω(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return self.ra_of_asc_node(epoch)
    
    @abstractmethod
    def arg_of_pericenter(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def ω(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        return self.arg_of_pericenter(epoch)
    
    @abstractmethod
    def mean_anomaly(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        raise NotImplementedError("A `Star` object has to implement all orbital elements: `a, e, i, Ω, ω, M`")
    def M(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        return self.mean_anomaly(epoch)
    
    def mean_motion(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)]:
        """Returns the mean motion (n) of this star at given epoch."""
        return semi_major_axis_earth_to_mean_motion(self.semi_major_axis(epoch))
    def n(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian / U.second, kind=TensorKind.SCALAR)]:
        """Returns the mean motion (n) of this star at given epoch."""
        return self.mean_motion(epoch)
    
    def radius(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        a = self.semi_major_axis(epoch)
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        return a * (1 - e * cos(E))
    def r(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the radius from the focal of the ellipse's orbit 
        of this star at given epoch."""
        return self.radius(epoch)

    def true_anomaly(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the true anomaly (υ) of this star at given epoch."""
        e = self.eccentricity(epoch)
        E = self.eccentric_anomaly(epoch)
        tan_half_nu = ((1 + e) / (1 - e)) ** 0.5 * tan(.5 * E)
        nu = 2 * atan(tan_half_nu)
        return nu
    def υ(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the true anomaly (υ) of this star at given epoch."""
        return self.true_anomaly(epoch)
    
    def eccentric_anomaly(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        M = self.mean_anomaly(epoch)
        e = self.eccentricity(epoch)
        E = M.copy()
        for _ in range(5):  # compute excentric anomaly
            E = e * sin(E) + M
        return E
    def E(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the eccentric anomaly (E) of this star at given epoch."""
        return self.eccentric_anomaly(epoch)
    

class Moon(Body):
    """`Body` object representing the Moon.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """
    def __init__(self):
        super().__init__("Moon", scalar(1737.4).with_units(U.kilometer), scalar(7.342e22).with_units(U.kilogram))

    def semi_major_axis(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * scalar(60.2666 * 6_378_137).with_units(U.meter)
    
    def inclination(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * (0.08980417133211624 * U.radian)
    
    def eccentricity(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * 0.054900
    
    def arg_of_pericenter(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        argp = Tensor(
            318.0634 + 0.1643573223 * t.scalar.values("day")
        ) * U.degree
        return normalize_angle(argp)
    
    def ra_of_asc_node(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        raan = Tensor(
            125.1228 - 0.0529538083 * t.scalar.values("day")
        ) * U.degree
        return normalize_angle(raan)
    
    def mean_anomaly(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        M = Tensor(
            115.3654 + 13.0649929509 * t.scalar.values("day")
        ) * U.degree
        return normalize_angle(M)
    
class Sun(Body):
    """`Body` object representing the Sun.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """

    def __init__(self):
        super().__init__("Sun", scalar(696_340).with_units(U.kilometer), scalar(1.988e30).with_units(U.kilogram))
    
    def semi_major_axis(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.meter, kind=TensorKind.SCALAR)]:
        """Returns the semi major axis (a) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * (149_597_870_700 * U.meter)
    
    def inclination(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the inclination (i) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * (0 * U.radian)
    
    def eccentricity(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.dimensionless, kind=TensorKind.SCALAR)]:
        """Returns the eccentricity (e) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        e = Tensor(
            0.016709 - 1.151e-9 * t.scalar.values("day")
        )
        return e
    
    def arg_of_pericenter(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        argp = Tensor(
            282.9404 + 4.70935e-5 * t.scalar.values("day")
        ) * U.degree
        return normalize_angle(argp)
    
    def ra_of_asc_node(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        # to keep the shape of the input epoch (scalar or array) in the output
        one = __one(epoch)
        return one * (0 * U.radian)
    
    def mean_anomaly(self, epoch: Timestamp | TimeInterval) -> Annotated[Tensor, TensorBound(units=U.radian, kind=TensorKind.SCALAR)]:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t = from_mil(epoch.to_unixepoch())
        M = Tensor(
            356.0470 + 0.9856002585 * t.scalar.values("day")
        ) * U.degree
        return normalize_angle(M)