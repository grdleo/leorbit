from typing import NamedTuple

from leorbit2.mathematics.dimensions import Number


class OrbitalElementsComputeTuple(NamedTuple):
    n: Number
    """Mean motion [rad/min]"""

    i: Number
    """Inclination [rad]"""

    e: Number
    """Eccentricity [1]"""

    argp: Number
    """Argument of pericenter [rad]"""

    raan: Number
    """Right ascension of ascending node [rad]"""

    M: Number
    """Mean anomaly [rad]"""

    bstar: Number
    """BSTAR drag term [1/earthRadii]"""