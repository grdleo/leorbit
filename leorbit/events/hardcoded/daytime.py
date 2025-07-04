
from math import acos, pi

import numpy as np
from pint import Quantity
from numpy.typing import NDArray
from leorbit.coordinates.trajectory import Trajectory
from leorbit.events import Event
from leorbit.mathematics.units import Q_, UREG
from leorbit.mathematics.vec3 import Vec3
from leorbit.physics.constants import RADII_EARTH
from leorbit.physics.time_interval import TimeInterval
from leorbit.sky_objects.satellites import Satellite
from leorbit.sky_objects.stars.sun import Sun


def sun_visi(sat: Trajectory, sun: Trajectory) -> NDArray:
	assert sat._interval == sun._interval

	sat_x, sat_y, sat_z = sat._pos_vel_list.x, sat._pos_vel_list.y, sat._pos_vel_list.z
	sun_x, sun_y, sun_z = sun._pos_vel_list.x, sun._pos_vel_list.y, sun._pos_vel_list.z

	sat_rho: NDArray = (sat_x**2 + sat_y**2 + sat_z**2)**.5
	sun_rho: NDArray = (sun_x**2 + sun_y**2 + sun_z**2)**.5
	sat_sun_dp: NDArray = sat_x * sun_x + sat_y * sun_y + sat_z * sun_z
	sat_sun_angle: NDArray = np.arccos(sat_sun_dp / (sat_rho * sun_rho))
	shadow_angle = np.acos(RADII_EARTH.m_as(UREG.meter) / sat_rho)

	return sat_sun_angle <= pi*.5 + shadow_angle

class Daytime(Event):
	"""When given satellite experiences daytime (i.e. visible by the Sun)"""

	def __init__(self, satellite: Satellite):
		self._sat = satellite
		self._sun = Sun()

	def compute(self, on: TimeInterval):
		traj_sat = self._sat.trajectory(on)
		traj_sun = self._sun.trajectory(on)

		predicate_values: list[bool] = list(sun_visi(traj_sat, traj_sun))

		self._build_by_predicate(predicate_values, on.start, on.dt)
