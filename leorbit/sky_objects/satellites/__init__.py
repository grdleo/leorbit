from abc import ABC
from leorbit.coordinates.coordinates import Coordinates
from leorbit.coordinates.representations.elements import OrbitalElements
from leorbit.coordinates.trajectory import Trajectory
from leorbit.ext.celestrak import get_celestrak_gpdata_json
from leorbit.mathematics.units import Q_
from leorbit.physics.time import Time
from leorbit.physics.time_interval import TimeInterval
from leorbit.propagators import Propagator
from leorbit.propagators.no import NoPropagator
from leorbit.sky_objects import SkyObject


class Satellite(SkyObject):
	def __init__(self, name: str, propagator: Propagator):
		super().__init__(name, Q_("1m"), Q_("1000kg"))
		self.propagator = propagator

	@staticmethod
	def from_celestrak_norad_cat_id(norad_cat_id: int, propagator_cls: type[Propagator] = NoPropagator) -> "Satellite":
		els = OrbitalElements.from_celestrak_norad_cat_id(norad_cat_id)
		gp = get_celestrak_gpdata_json(norad_cat_id)

		return Satellite(
			gp["OBJECT_NAME"],
			propagator_cls(els)
		)
	
	def coordinates(self, at: Time) -> Coordinates:
		return self.propagator.propagate(at)
	
	def trajectory(self, during: TimeInterval) -> Trajectory:
		return self.propagator.propagate_timeline(during)