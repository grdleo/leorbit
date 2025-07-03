
from math import acos, pi

from pint import Quantity
from leorbit.events import Event
from leorbit.mathematics.units import Q_, UREG
from leorbit.mathematics.vec3 import Vec3
from leorbit.physics.constants import RADII_EARTH
from leorbit.physics.time_interval import TimeInterval
from leorbit.sky_objects.satellites import Satellite
from leorbit.sky_objects.stars.sun import Sun

HALF_PI_RAD = Q_(pi / 2, UREG.radians)

def sun_visi(sat_gcrf: Vec3, sun_gcrf: Vec3) -> bool:
	shadow_angle: Quantity = acos(RADII_EARTH / sat_gcrf.rho)
	return sat_gcrf.angle(sun_gcrf) <= HALF_PI_RAD + shadow_angle

class Daytime(Event):
	"""When given satellite experiences daytime (i.e. visible by the Sun)"""

	def __init__(self, satellite: Satellite):
		self._sat = satellite
		self._sun = Sun()

	def compute(self, on: TimeInterval):
		traj_sat = self._sat.trajectory(on)
		traj_sun = self._sun.trajectory(on)

		predicate_values = [
			sun_visi(
				traj_sat._pos_vel_list.get_values_at(i).pos,
				traj_sun._pos_vel_list.get_values_at(i).pos
			)
			for i in range(on.steps)
		]

		self._build_by_predicate(predicate_values, on.start, on.dt)
