from math import pi
from leorbit2.mathematics.dimensions import D, Number
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar

TWELF_PI = pi / 12
TWOPI = 2 * pi

FULL_REV = (2 * pi) * Quantity.rad
HALF_REV = FULL_REV / 2

def normalize_angle(angle: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """Returns the given angle in its [0, 2π] range."""
    return angle % FULL_REV

def normalize_angle_symmetric(angle: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """Returns the given angle in its [-π, π] range."""
    normalized = normalize_angle(angle)
    return normalized if normalized <= HALF_REV else normalized - FULL_REV

def angle2dms(angle: Scalar[D.Angle]) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle2convert = abs(angle.magnitude("deg"))
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{deg}° {min}′ {sec}″"