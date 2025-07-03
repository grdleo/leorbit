import numpy as np
from pint import Quantity
from leorbit.algorithms.utils import true2eccentric_anomaly, true2eccentric_anomaly_numpy
from leorbit.coordinates.coordinates import Coordinates
from math import sin, cos
from types import ModuleType
from numpy.typing import NDArray

from leorbit.coordinates.pos_vel_tuple.gcrf import PosVelGCRF
from leorbit.frames.absolute_frame import AbsoluteFrame
from leorbit.mathematics.vec3 import Vec3
from leorbit.physics.constants import SQRT_MU_EARTH as SQRT_MU_EARTH_PINT

STD_GRAV_PARAM_TERRA_FLOAT = SQRT_MU_EARTH_PINT.m_as("m**1.5/s")

def elements2orthogonal_gcrf(υ: NDArray | np.float64, e: float, a: float, Ω: float, ω: float, i: float) -> PosVelGCRF:
	"""Returns the position and velocity of a satellite in GCRF coordinates (meters, meters/second)
	All parameters are in radians, except `e` dimensionless and `a` in meters."""

	# position of satellite in orbit plane (with z = 0)
	
	ee = e**2
	one_ee = (1 - ee)
	E = true2eccentric_anomaly_numpy(e, υ)
	esinE = e * np.sin(E)
	r = a * one_ee / (1 + e * cos(υ))
	rd = STD_GRAV_PARAM_TERRA_FLOAT * a**.5 * esinE / r
	rυd = rd * one_ee / esinE

	c_raan, s_raan = np.cos(Ω), np.sin(Ω)
	c_i, s_i = np.cos(i), np.sin(i)
	υpω = υ + ω
	c_theta, s_theta = np.cos(υpω), np.sin(υpω)

	def unitvec_gcrf(x: np.float64 | NDArray, y: np.float64 | NDArray) -> tuple[NDArray, NDArray, NDArray]:
		return (
			c_raan * x - s_raan * c_i * y, # X
			s_raan * x + c_raan * c_i * y, # Y
			s_i * y                        # Z
		)

	ur_x, ur_y, ur_z = unitvec_gcrf(c_theta, s_theta)
	ut_x, ut_y, ut_z = unitvec_gcrf(-s_theta, c_theta)

	return PosVelGCRF(
		ur_x * r, 
		ur_y * r, 
		ur_z * r,
		ur_x * rd + ut_x * rυd,
		ur_y * rd + ut_y * rυd,
		ur_z * rd + ut_z * rυd
	)