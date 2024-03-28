"""Math tools such as vectors, quaternions, etc. """

__pdoc__ = dict()

from math import pi
from pint import UnitRegistry

UREG = UnitRegistry()
PintQuantity = UREG.Quantity
Q_ = PintQuantity

PI = pi
HALFPI = 0.5 * pi
TWELF_PI = pi / 12
TWOPI = 2 * pi
SQRT_2 = 1.4142135623730951
SQRT_3 = 1.7320508075688772
RIGHT_ANGLE = Q_("90°")
FULL_REV = Q_("360°")
HALF_REV = Q_("180°")

from leorbit.math.vector import Vec3
from leorbit.math.coordinate import Coordinates, AbsoluteFrame, RelativeFrame, GPSCoordinates, HorizontalCoordinates, PosVelGCRF

# DOC: hide from documentation
__pdoc__["leorbit.math.orientation"] = False
__pdoc__["leorbit.math.quaternion"] = False
__pdoc__["leorbit.math.visibility"] = False