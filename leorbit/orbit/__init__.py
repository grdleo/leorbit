"""Utilitaries for objets orbiting Earth and its propagation."""

from leorbit.orbit.objects import Satellite, Star, Sun, Moon
from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.orbit.propagator import Propagator, NoPropagator, SGP4

SUN = Sun()
"""Global instance of the `leorbit.orbit.objects.Sun` object."""


MOON = Moon()
"""Global instance of the `leorbit.orbit.objects.Moon` object."""