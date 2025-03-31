from functools import lru_cache

from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from events.timeline import CoordinatesTimeline
from mathematics.custom import FULL_REV
from physics.time import Time
from physics.time_interval import TimeInterval
from propagators import Propagator

class NoPropagator(Propagator):
    def __init__(self, elements: OrbitalElements):
        super().__init__(elements)
    
    @lru_cache
    def propagate(self, at: Time) -> Coordinates:
        els = self.elements
        shifted_mean_anomaly = (
            els.mean_anomaly 
            + els.mean_motion * at.delta(els.epoch)
        ) % FULL_REV
        
        return OrbitalElements(
            epoch=at,
            eccentricity=els.eccentricity,
            inclination=els.inclination,
            ra_of_asc_node=els.ra_of_asc_node,
            arg_of_pericenter=els.arg_of_pericenter,
            mean_motion=els.mean_motion,
            mean_anomaly=shifted_mean_anomaly,
            mean_motion_dot=els.mean_motion_dot,
            mean_motion_ddot=els.mean_motion_ddot,
            bstar=els.bstar
        ).to_coordinates()
        
    def propagate_timeline(self, on: TimeInterval) -> CoordinatesTimeline:
        coordinates_computer = self.propagate
        return CoordinatesTimeline(on, coordinates_computer)