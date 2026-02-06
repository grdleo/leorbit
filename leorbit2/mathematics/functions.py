from math import pi
from leorbit2.mathematics.dimensions import Number, D
from leorbit2.mathematics.scalar import Scalar

TWELF_PI = pi / 12
TWOPI = 2 * pi

FULL_REV = TWOPI * UREG.radians
HALF_REV = FULL_REV / 2

def normalize_angle(angle: Number | Scalar[D.Angle]) -> Q_:
    """Returns the given angle in its [0, 2π] range."""
    return angle % FULL_REV

def normalize_angle_symmetric(angle: AnyNumber) -> Q_:
    """Returns the given angle in its [-π, π] range."""
    normalized = normalize_angle(angle)
    return normalized if normalize_angle <= HALF_REV else normalized - FULL_REV