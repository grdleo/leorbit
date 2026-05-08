"""Public end-user API for LEOrbit.

This module centralizes the most useful classes and helpers for typical
satellite tracking workflows.
"""

from typing import TYPE_CHECKING

from leorbit.coordinates import Coordinates, GPS, Horizontal, OrbitalElements, Trajectory
from leorbit.events import Event, TimeMap, VisibleFromEarthLocationEvent
from leorbit.ext import CelestrakDataGP
from leorbit.frames import AbsoluteFrame, EarthLocalFrame
from leorbit.mathematics import Quantity, Tensor, matrix33, scalar, vector3
from leorbit.propagator import NoPropagator, Propagator, SGP4
from leorbit.sky_object import Moon, Satellite, SkyObject, Sun
from leorbit.time import TimeInterval, TimeIntervalSet, Timestamp, Timeline, get_intersections_timelines

if TYPE_CHECKING:
	from leorbit.mathematics import Tensor as _Tensor


def get_satellite(
	norad_cat_id: int,
	propagator: type[Propagator] = SGP4,
	log: bool = True,
) -> Satellite:
	"""Build a ``Satellite`` from latest Celestrak GP data.

	Parameters
	----------
	norad_cat_id:
		NORAD catalog identifier of the satellite.
	propagator:
		Propagator class to use. Defaults to ``SGP4``.
	log:
		Whether to print fetch logs.
	"""
	from leorbit.ext import get_celestrak_gpdata

	gp_data = get_celestrak_gpdata(norad_cat_id, log)

	return Satellite(
		gp_data.name,
		gp_data.to_orbital_elements(),
		propagator,
	)

def get_passes(
    satellite: Satellite,
	during: TimeInterval,
	gps_observer: GPS,
	altitude_angle_min_degrees: float = 0.
) -> TimeIntervalSet:
	return VisibleFromEarthLocationEvent(
		satellite.trajectory(during),
		gps_observer,
		altitude_angle_min_degrees * Quantity.degree
    ).visible_intervals
	

def Q_(unit_name: str) -> Tensor:
	"""Return a unit quantity by its registered name.

	Examples
	--------
	``Q_("day")``
	``5 * Q_("second")``
	"""
	return Quantity.get(unit_name)


__all__ = [
	"AbsoluteFrame",
	"CelestrakDataGP",
	"Coordinates",
	"EarthLocalFrame",
	"Event",
	"GPS",
	"Horizontal",
	"Moon",
	"NoPropagator",
	"OrbitalElements",
	"Propagator",
	"Q_",
	"Quantity",
	"SGP4",
	"Satellite",
	"SkyObject",
	"Sun",
	"Tensor",
	"TimeInterval",
	"TimeMap",
	"Timeline",
	"Timestamp",
	"Trajectory",
	"VisibleFromEarthLocationEvent",
	"get_intersections_timelines",
	"get_satellite",
	"matrix33",
	"scalar",
	"vector3",
]
