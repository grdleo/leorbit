from abc import ABC
from functools import lru_cache
from math import atan, cos, sin, sqrt, tan
from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from mathematics.custom import FULL_REV
from physics.constants import MU_EARTH
from physics.time import Time

from pint import Quantity as Q_

from sky_objects.stars import Star

class Sun(Star):
    """`Star` object representing the Sun.
    
    Orbital elements for this object taken from here: https://stjarnhimlen.se/comp/ppcomp.html#4
    """

    def __init__(self):
        super().__init__("Sun", Q_(696_340, "km"), Q_(1.988e30, "kg"))
    
    def semi_major_axis(self, epoch: Time = None) -> Q_:
        """Returns the semi major axis (a) of this star at given epoch."""
        return Q_(149_597_870_700, "m")
    
    def inclination(self, epoch: Time = None) -> Q_:
        """Returns the inclination (i) of this star at given epoch."""
        return Q_(0, "rad")
    
    def eccentricity(self, epoch: Time) -> Q_:
        """Returns the eccentricity (e) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        e = 0.016709 - 1.151e-9 * t.m_as("day")
        return Q_(e, "dimensionless")
    
    def arg_of_pericenter(self, epoch: Time) -> Q_:
        """Returns the argument of pericenter (ω) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        argp = 282.9404 + 4.70935e-5 * t.m_as("day")
        return Q_(argp, "°") % FULL_REV
    
    def ra_of_asc_node(self, epoch: Time = None) -> Q_:
        """Returns the right ascension of the ascending node (Ω) of this star at given epoch."""
        return Q_(0, "rad")
    
    def mean_anomaly(self, epoch: Time) -> Q_:
        """Returns the mean anomaly (M) of this star at given epoch."""
        t: Q_ = epoch.from_mil
        M = 356.0470 + 0.9856002585 * t.m_as("day")
        return Q_(M, "°") % FULL_REV