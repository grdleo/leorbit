from functools import lru_cache
from algorithms.sgp4 import PosVelGCRF, sgp4
from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from leorbit.coordinates.trajectory import Trajectory
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
    
    def propagate_timeline(self, on: TimeInterval) -> Trajectory:
        tsince_min = (on.to_time_stamps() - on.start.unixepoch) / 60.
        pos_vel = sgp4(
            self.elements._els_as_float_tuple,
            tsince_min
        )

        return Trajectory(pos_vel, on)