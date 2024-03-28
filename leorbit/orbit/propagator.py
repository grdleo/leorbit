"""Implementation of orbital elements propagation."""

from abc import ABC, abstractmethod
from functools import lru_cache
from leorbit.math.vector import Vec3
from leorbit.math import Q_, FULL_REV
from leorbit.math.coordinate import Coordinates, AbsoluteFrame, PosVelGCRF
from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.time import Time
from math import tau

from leorbit.algorithms.sgp4 import sgp4

from time import perf_counter_ns

from leorbit.time import Timeline

class Propagator(ABC):
    """Object that propagates given `OrbitalElements` in time.
    
    The principal element that is propagated is the mean anomaly (M), 
    which corresponds to the position of the satellite on its orbit.

    Perturbations might affect the actual orbit, so all the elements 
    might change.
    """
    start_elements: OrbitalElements

    def __init__(self, start_elements: OrbitalElements):
        self._start_elements = start_elements
    
    @property
    def start_elements(self) -> OrbitalElements:
        """The `OrbitalElements` to be propagated."""
        return self._start_elements
    
    def propagate(self, epoch: Time) -> Coordinates:
        """Returns the position and velocity as a `Coordinates` object 
        at given epoch."""
        gcrf_x, gcrf_y, gcrf_z, gcrf_vx, gcrf_vy, gcrf_vz = self._fast_propagate_gcrf(epoch)

        gcrf_pos = Vec3(
            Q_(gcrf_x, "m"),
            Q_(gcrf_y, "m"),
            Q_(gcrf_z, "m")
        )
        gcrf_vel = Vec3(
            Q_(gcrf_vx, "m/s"),
            Q_(gcrf_vy, "m/s"),
            Q_(gcrf_vz, "m/s")
        )

        return Coordinates(AbsoluteFrame.GCRF, gcrf_pos, gcrf_vel, epoch)
    
    def propagate_timeline(self, timeline: Timeline) -> list[Coordinates]:
        """Returns the positions and velocities as a list of `Coordinates`,
        for every instant in the given `Timeline`"""
        return [self.propagate(t) for t in timeline]

    def _fast_propagate_gcrf(self, epoch: Time) -> PosVelGCRF:
        """Propagates coordinates at given epoch 
        and returns the (pos, vel) GCRF vector as a tuple.
        
        Holds computed value in cache"""
        raise NotImplementedError()

class NoPropagator(Propagator):
    """Propagator that does not implement any perturbations.

    Result of the propagation will be the same ellipse with 
    mean anomaly (M) shifted as a regular Keplerian orbit."""
    def __init__(self, start_elements: OrbitalElements):
        super().__init__(start_elements)
        
    @lru_cache(maxsize=None)
    def _fast_propagate_gcrf(self, epoch: Time) -> PosVelGCRF:
        els = self.start_elements
        M = (els.mean_anomaly + els.mean_motion * epoch.delta(els.epoch)) % FULL_REV
        new_els = OrbitalElements(
            epoch=epoch,
            eccentricity=els.eccentricity,
            inclination=els.inclination,
            ra_of_asc_node=els.ra_of_asc_node,
            arg_of_pericenter=els.arg_of_pericenter,
            mean_motion=els.mean_motion,
            mean_anomaly=M,
            mean_motion_dot=els.mean_motion_dot,
            mean_motion_ddot=els.mean_motion_ddot,
            bstar=els.bstar
        )
        return new_els.to_coordinates().to_gcrf_pos_vel_tuple()


class SGP4(Propagator):
    """
        Propagator that implements SGP4 propagation.
        Is relevant for small objects (satellites) orbiting Earth.

        It takes osculting orbital elements of the object at a given epoch,
        and propagates its position to a later epoch, by using laws of physics
        (aka using Kepler equations but taking count perturbations such as 
        astmospheric drag and the geopotential model of the Earth)
    """
    def __init__(self, els: OrbitalElements):
        super().__init__(els)
        self._els_as_list = list(els._els_ready2compute.values()) # Careful, order is important for SGP4 algorithm arguments!
        
    @lru_cache(maxsize=None)
    def _fast_propagate_gcrf(self, epoch: Time, print_delay: bool = False) -> PosVelGCRF:
        tsince = epoch.delta(self.start_elements.epoch).m_as("min")
        t0 = perf_counter_ns()
        pos_vel_gcrf = sgp4(*self._els_as_list, tsince)
        t1 = perf_counter_ns()
        if print_delay:
            delay_mis = (t1 - t0) * 1e-3
            print(f"--- SGP4 algorithm took {int(delay_mis)}µs to run")
        return pos_vel_gcrf
        