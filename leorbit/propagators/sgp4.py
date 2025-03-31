from functools import lru_cache
from algorithms.sgp4 import PosVelGCRF, sgp4
from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from events.timeline import CoordinatesTimeline
from physics.time import Time
from physics.time_interval import TimeInterval
from propagators import Propagator
from mathematics.units import UREG
    

class SGP4(Propagator):
    def __init__(self, elements: OrbitalElements):
        super().__init__(elements)
    
    @lru_cache
    def propagate(self, to: Time):
        return sgp4(
            self.elements._els_as_float_tuple,
            to.delta(self.elements.epoch).m_as(UREG.min)
        ).to_coordinates(to)
    
    def propagate_timeline(self, on: TimeInterval) -> CoordinatesTimeline:
        tsince: list[float] = [t.delta(self.elements.epoch).m_as(UREG.min) for t in on]
        gcrf_positions_arrays = sgp4(self.elements._els_as_float_tuple, tsince)
        
        def coordinates_computer(at: Time) -> Coordinates:
            i = on._time2idx(at)
            return gcrf_positions_arrays.get_values_at(i).to_coordinates(at)
        
        return CoordinatesTimeline(on, coordinates_computer)
    
    