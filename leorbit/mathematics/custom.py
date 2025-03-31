from math import pi
from mathematics.units import UREG

from pint import Quantity as Q_

Number = float | int
AnyNumber = float | int | Q_

TWELF_PI = pi / 12
TWOPI = 2 * pi

FULL_REV = TWOPI * UREG.radians
HALF_REV = FULL_REV / 2

def normalize_angle(angle: AnyNumber) -> Q_:
    """Returns the given angle in its [0, 2π] range."""
    return angle % FULL_REV

def normalize_angle_symmetric(angle: AnyNumber) -> Q_:
    """Returns the given angle in its [-π, π] range."""
    normalized = normalize_angle(angle)
    return normalized if normalize_angle <= HALF_REV else normalized - FULL_REV