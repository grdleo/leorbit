from functools import lru_cache

from pint import Quantity
from numpy.typing import NDArray

from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from events.timeline import CoordinatesTimeline
from leorbit.algorithms.elements2orthogonal_gcrf import elements2orthogonal_gcrf
from leorbit.algorithms.utils import mean2true_anomaly
from leorbit.coordinates.pos_vel_tuple.gcrf import PosVelGCRF
from leorbit.coordinates.trajectory import Trajectory
from leorbit.mathematics.units import UREG
from mathematics.custom import FULL_REV
from physics.time import Time
from physics.time_interval import TimeInterval
from propagators import Propagator

class NoPropagator(Propagator):
    def __init__(self, elements: OrbitalElements):
        super().__init__(elements)

    def _compute_posvel_gcrf(self, time: Time | TimeInterval) -> PosVelGCRF:
        els = self.elements
        shifted_mean_ano: Quantity | NDArray

        if isinstance(time, Time):
            shifted_mean_ano: Quantity = els.mean_anomaly + els.mean_motion * time.delta(els.epoch)
            shifted_mean_ano %= FULL_REV
            shifted_mean_ano = shifted_mean_ano.m_as(UREG.radians)
        elif isinstance(time, TimeInterval):
            time_line = time.to_time_stamps() - time.start.unixepoch
            shifted_mean_ano: NDArray = els.mean_anomaly.m_as(UREG.radians) + els.mean_motion.m_as(UREG.radians / UREG.second) * time_line
            shifted_mean_ano %= FULL_REV.m_as(UREG.radians)
        else:
            raise TypeError()

        nu: NDArray | Quantity = mean2true_anomaly(els.eccentricity.magnitude, shifted_mean_ano)
        
        return elements2orthogonal_gcrf(
            nu,
            els.eccentricity.m,
            els.semi_major_axis.m_as(UREG.meter),
            els.ra_of_asc_node.m_as(UREG.radians),
            els.arg_of_pericenter.m_as(UREG.radians),
            els.inclination.m_as(UREG.radians)
        )
    
    @lru_cache
    def propagate(self, at: Time) -> Coordinates:
        pos_vel = self._compute_posvel_gcrf(at)
        return pos_vel.to_coordinates(at)
        
    def propagate_timeline(self, on: TimeInterval) -> Trajectory:
        pos_vel = self._compute_posvel_gcrf(on)
        return Trajectory(pos_vel, on)