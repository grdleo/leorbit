"""Public end-user API for LEOrbit.

This module centralizes the most useful classes and helpers for typical
satellite tracking workflows.
"""

from typing import TYPE_CHECKING, ClassVar, Literal

from leorbit.mathematics import U as UnitRegistry, RealNumber
from leorbit.mathematics import Tensor, matrix33, scalar, vector3
from leorbit.coordinates import Coordinates, GPS, Horizontal, OrbitalElements, Trajectory
from leorbit.events import Event, TimeMap, VisibleFromEarthLocationEvent
from leorbit.ext import CelestrakDataGP
from leorbit.frames import AbsoluteFrame, EarthLocalFrame
from leorbit.utils import _radii_earth_unit

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
		U.second interval over which visibility is searched.
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
		scalar(altitude_angle_min_degrees).with_units(UnitRegistry.degree)
    ).visible_intervals

class Qty:
	"""Convenience class that exposes scalar tensor of value 1,
	for all the most used units
	
	Example usage:
	```
	from leorbit.api import Qty
	altitude = 400 * Qty.km
	```
	"""

	deg: ClassVar[Tensor] = scalar(1).with_units("degree")
	"""degree"""

	degree: ClassVar[Tensor] = deg
	"""degree"""

	radian: ClassVar[Tensor] = scalar(1).with_units("radian")
	"""radian"""

	rad: ClassVar[Tensor] = radian
	"""radian"""

	m: ClassVar[Tensor] = scalar(1).with_units("meter")
	"""meter"""

	meter: ClassVar[Tensor] = m
	"""meter"""

	km: ClassVar[Tensor] = scalar(1).with_units("kilometer")
	"""kilometer"""

	kilometer: ClassVar[Tensor] = km
	"""kilometer"""

	radii_earth: ClassVar[Tensor] = scalar(1).with_units(_radii_earth_unit)
	"""Earth radius"""

	second: ClassVar[Tensor] = scalar(1).with_units("second")
	"""second"""

	s: ClassVar[Tensor] = scalar(1).with_units("second")
	"""second"""

	minute: ClassVar[Tensor] = scalar(1).with_units("minute")
	"""minute"""

	min: ClassVar[Tensor] = scalar(1).with_units("minute")
	"""minute"""

	h: ClassVar[Tensor] = scalar(1).with_units("hour")
	"""hour"""

	hour: ClassVar[Tensor] = scalar(1).with_units("hour")
	"""hour"""

	day: ClassVar[Tensor] = scalar(1).with_units("day")
	"""day"""

	year: ClassVar[Tensor] = scalar(1).with_units("year")
	"""year"""

	kg: ClassVar[Tensor] = scalar(1).with_units("kilogram")
	"""kilogram"""

	g: ClassVar[Tensor] = scalar(1).with_units("g")
	"""gram"""

class _ConvenienceBuilderGPS:
	"""Convenience builder for GPS coordinates on the cardinal axes."""

	def __init__(self, 
		coords_type: Literal["north", "south", "east", "west"]
	):
		self._lon = scalar(
			1 if coords_type == "east" else (
				-1 if coords_type == "west" else 0
			)
		).with_units("deg")

		self._lat = scalar(
			1 if coords_type == "north" else (
				-1 if coords_type == "south" else 0
			)
		).with_units("deg")

	def __mul__(self, other: object) -> GPS:
		if not isinstance(other, RealNumber):
			raise TypeError(f"Multiplication with type {type(other)} not supported.")
		return GPS(
			longitude=self._lon * other,
			latitude=self._lat * other,
			altitude=scalar(0).with_units("m")
		)
	
	def __rmul__(self, other: object) -> GPS:
		if not isinstance(other, RealNumber):
			raise TypeError(f"Multiplication with type {type(other)} not supported.")
		return GPS(
			longitude=self._lon * other,
			latitude=self._lat * other,
			altitude=scalar(0).with_units("m")
		)


N = _ConvenienceBuilderGPS("north")
"""GPS builder for a point on the northern cardinal axis (latitude = +1°)."""

S = _ConvenienceBuilderGPS("south")
"""GPS builder for a point on the southern cardinal axis (latitude = -1°)."""

E = _ConvenienceBuilderGPS("east")
"""GPS builder for a point on the eastern cardinal axis (longitude = +1°)."""

W = _ConvenienceBuilderGPS("west")
"""GPS builder for a point on the western cardinal axis (longitude = -1°)."""

__all__ = [
	"Qty",
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
	"N",
	"S",
	"E",
	"W",
]
