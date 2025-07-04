
from math import acos, pi

import numpy as np
from pint import Quantity
from numpy.typing import NDArray
from leorbit.coordinates.trajectory import Trajectory
from leorbit.events import Event
from leorbit.mathematics.units import Q_, UREG
from leorbit.mathematics.vec3 import Vec3
from leorbit.physics.constants import RADII_EARTH, RADII_SUN
from leorbit.physics.time_interval import TimeInterval
from leorbit.sky_objects.satellites import Satellite
from leorbit.sky_objects.stars.sun import Sun

R_EARTH_SI = RADII_EARTH.m_as(UREG.meter)
R_SUN_SI = RADII_SUN.m_as(UREG.meter)
HALF_PI = pi/2

def sun_visi(sat: Trajectory, sun: Trajectory) -> NDArray:
	assert sat._interval == sun._interval

	sat_x, sat_y, sat_z = sat._pos_vel_list.x, sat._pos_vel_list.y, sat._pos_vel_list.z
	sun_x, sun_y, sun_z = sun._pos_vel_list.x, sun._pos_vel_list.y, sun._pos_vel_list.z

	sat_rho: NDArray = (sat_x**2 + sat_y**2 + sat_z**2)**.5
	sun_rho: NDArray = (sun_x**2 + sun_y**2 + sun_z**2)**.5
	sat_sun_dp: NDArray = sat_x * sun_x + sat_y * sun_y + sat_z * sun_z
	sat_sun_angle: NDArray = np.arccos(sat_sun_dp / (sat_rho * sun_rho))
	theta: NDArray = np.asin((R_SUN_SI - R_EARTH_SI) / sun_rho)
	phi: NDArray = np.asin((R_SUN_SI + R_EARTH_SI) / sun_rho)

	if sat_sun_angle >= HALF_PI - phi:
		pass # full!

	b_fact = R_EARTH_SI / (R_SUN_SI + R_EARTH_SI)
	bx, by, bz = sun_x * b_fact, sun_y * b_fact, sun_z * b_fact
	bs_x, bs_y, bs_z = sat_x - bx, sat_y - by, sat_z - bz
	bs_rho = (bs_x**2 + bs_y**2 + bs_z**2)**.5
	bs_angle = np.acos((bs_x * sun_x + bs_y * sun_y + bs_z * sun_z) / (bs_rho * sun_rho))

	a_fact = -R_EARTH_SI / (R_SUN_SI - R_EARTH_SI)
	ax, ay, az = sun_x * a_fact, sun_y * a_fact, sun_z * a_fact
	as_x, as_y, as_z = sat_x - ax, sat_y - ay, sat_z - az
	as_rho = (as_x**2 + as_y**2 + as_z**2)**.5
	as_angle = np.acos((as_x * sun_x + as_y * sun_y + as_z * sun_z) / (as_rho * sun_rho))

	raise NotImplementedError()
	


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
