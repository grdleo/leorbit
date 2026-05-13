"""Public end-user API for LEOrbit.

This module centralizes the most useful classes and helpers for typical
satellite tracking workflows.
"""

from typing import TYPE_CHECKING

from leorbit.mathematics import U as UnitRegistry
from leorbit.mathematics import Tensor, matrix33, scalar, vector3
from leorbit.coordinates import Coordinates, GPS, Horizontal, OrbitalElements, Trajectory
from leorbit.events import Event, TimeMap, VisibleFromEarthLocationEvent
from leorbit.ext import CelestrakDataGP
from leorbit.frames import AbsoluteFrame, EarthLocalFrame

from leorbit.propagator import NoPropagator, Propagator, SGP4
from leorbit.sky_object import Moon, Satellite, SkyObject, Sun
from leorbit.time import TimeInterval, TimeIntervalSet, Timestamp, Timeline, get_intersections_timelines


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
	"""Compute visibility intervals of a satellite from a ground location.

	A pass is considered visible when the satellite altitude (elevation)
	above the local horizon is greater than or equal to
	``altitude_angle_min_degrees``.

	Parameters
	----------
	satellite:
		Satellite to evaluate.
	during:
		Time interval over which visibility is searched.
	gps_observer:
		Observer geodetic location.
	altitude_angle_min_degrees:
		Minimum elevation angle in degrees for visibility. Defaults to ``0.0``
		(geometric horizon).

	Returns
	-------
	TimeIntervalSet
		Set of visible pass intervals for ``satellite`` during ``during`` from
		``gps_observer``.
	"""
	return VisibleFromEarthLocationEvent(
		satellite.trajectory(during),
		gps_observer,
		altitude_angle_min_degrees * UnitRegistry.degree
    ).visible_intervals


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
	"UnitRegistry",
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
