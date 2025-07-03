from pint import UnitRegistry

UREG = UnitRegistry()
UREG.define("earthRadii = 6378135 * m")

POS_UNIT = UREG.meter
VEL_UNIT = UREG.meter / UREG.second
ACC_UNIT = UREG.meter / UREG.second**2

Q_ = UREG.Quantity