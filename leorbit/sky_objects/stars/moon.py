from abc import ABC
from functools import lru_cache
from math import atan, cos, sin, sqrt, tan
from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from mathematics.custom import FULL_REV
from physics.constants import MU_EARTH, RADII_EARTH
from physics.time import Time

from pint import Quantity as Q_

from sky_objects.stars import Star

class Moon(Star):
    """`Star` object representing the Moon.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """
    def __init__(self):
        super().__init__("Moon", Q_(1737.4, "km"), Q_(7.342e22, "kg"))

    def semi_major_axis(self, epoch: Time = None) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return 60.2666 * RADII_EARTH
    
    def inclination(self, epoch: Time = None) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return Q_(0.08980417133211624, "rad")
    
    def eccentricity(self, epoch: Time = None) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        return Q_(0.054900, "dimensionless")
    
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        argp = 318.0634 + 0.1643573223 * t.m_as("day")
        return Q_(argp, "°") % FULL_REV
    
    def ra_of_asc_node(self, epoch: Time) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        raan = 125.1228 - 0.0529538083 * t.m_as("day")
        return Q_(raan, "°") % FULL_REV
    
    def mean_anomaly(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        M = 115.3654 + 13.0649929509 * t.m_as("day")
        return Q_(M, "°") % FULL_REV