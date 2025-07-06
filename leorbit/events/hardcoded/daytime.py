
from enum import Enum
from math import acos, pi

import numpy as np
from pint import Quantity
from numpy.typing import NDArray
from leorbit.coordinates.trajectory import Trajectory
from leorbit.events import Event
from leorbit.mathematics.units import Q_, UREG
from leorbit.mathematics.vec3 import Vec3
from leorbit.mathematics.vec3_numpy_array import Vec3NumpyArray
from leorbit.physics.constants import RADII_EARTH, RADII_SUN
from leorbit.physics.time_interval import TimeInterval
from leorbit.sky_objects.satellites import Satellite
from leorbit.sky_objects.stars.sun import Sun

RT: float = RADII_EARTH.m_as(UREG.meter)
RS: float = RADII_SUN.m_as(UREG.meter)
HALF_PI: float = pi/2

class SunShadowState(Enum):
	# https://en.wikipedia.org/wiki/Umbra,_penumbra_and_antumbra
	ENLIGHTEN = 1
	UMBRA = 2
	PENUMBRA = 3
	ANTUMBRA = 4

def sun_visi(traj_sat: Trajectory, traj_sun: Trajectory) -> NDArray:
	# See `leorbit\assets\sun-visibility-schema.ggb``
	assert traj_sat._interval == traj_sun._interval

	sat = Vec3NumpyArray(traj_sat._pos_vel_list.x, traj_sat._pos_vel_list.y, traj_sat._pos_vel_list.z)
	sun = Vec3NumpyArray(traj_sun._pos_vel_list.x, traj_sun._pos_vel_list.y, traj_sun._pos_vel_list.z)

	u = sun.normalized
	w = u.cross(sat.normalized)
	v = u.cross(w)

	x = sat.dot(u)
	y = sat.dot(v)

	alpha: NDArray = np.asin((RS - RT) / sun.rho)
	beta: NDArray = np.asin((RS + RT) / sun.rho)

	b = sun.rho * (RT / (RS + RT))
	a = -sun.rho * (RT / (RS - RT))
	z0 = RT * np.cos(HALF_PI - beta)
	z1 = -RT * np.cos(HALF_PI - alpha)

	_below_bline = np.abs(b - x) * np.tan(beta) >= np.abs(y)
	_below_aline = np.abs(a - x) * np.tan(alpha) >= np.abs(y)

	state = np.zeros(sat.shape, dtype=np.uint8)
	state[x >= z0] = SunShadowState.ENLIGHTEN.value

	within_z0_z1 = z0 > x & x >= z1
	state[within_z0_z1 & _below_bline] = SunShadowState.PENUMBRA.value
	state[within_z0_z1 & ~_below_bline] = SunShadowState.ENLIGHTEN.value

	umbra_zone = z1 > x & x >= a
	state[umbra_zone & ~_below_bline] = SunShadowState.ENLIGHTEN.value
	state[umbra_zone & _below_aline] = SunShadowState.UMBRA.value
	state[umbra_zone & _below_bline & ~_below_aline] = SunShadowState.PENUMBRA.value

	atumbra_zone = a > x
	state[atumbra_zone & ~_below_bline] = SunShadowState.ENLIGHTEN.value
	state[atumbra_zone & _below_aline] = SunShadowState.ANTUMBRA.value
	state[atumbra_zone & _below_bline & ~_below_aline] = SunShadowState.PENUMBRA.value

	raise state

class Daytime(Event):
	"""When given satellite experiences daytime (i.e. visible by the Sun)"""

	def __init__(self, satellite: Satellite):
		self._sat = satellite
		self._sun = Sun()
 
	def compute(self, on: TimeInterval):
		traj_sat = self._sat.trajectory(on)
		traj_sun = self._sun.trajectory(on)

		predicate_values: list[bool] = [SunShadowState.ENLIGHTEN.value == s for s in sun_visi(traj_sat, traj_sun).flatten()]

		self._build_by_predicate(predicate_values, on.start, on.dt)
