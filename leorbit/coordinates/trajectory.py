from functools import lru_cache
from leorbit.coordinates.coordinates import Coordinates
from leorbit.coordinates.pos_vel_tuple.gcrf import PosVelGCRF
from leorbit.frames.absolute_frame import AbsoluteFrame
from leorbit.physics.time import Time
from leorbit.physics.time_interval import TimeInterval


class Trajectory:
	def __init__(self, positions_velocities: PosVelGCRF, interval: TimeInterval):
		assert positions_velocities.count == interval.steps
		
		self._pos_vel_list = positions_velocities
		self._interval = interval

	@lru_cache
	def _get_posvel_at_index(self, idx: int) -> Coordinates:
		pos_vel = self._pos_vel_list.get_values_at(idx)
		assert not pos_vel.array_values
		return pos_vel.to_coordinates(self._interval._idx2time(idx))

	def __getitem__(self, time: Time) -> Coordinates:
		idx: int
		try:
			idx = self._interval._time2idx(time)
		except ValueError:
			raise ValueError(f"Given time {time} lies outside trajectory's time interval")
		
		return self._get_posvel_at_index(idx)