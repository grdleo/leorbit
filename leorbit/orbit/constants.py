"""Astronomy related constants.

As already mentionned, `leorbit` uses the `pint.Quantity` class to handle units.
This library also allows to implement custom units.

Therefore, `leorbit` implements these custom units:

- `earthRadii`, which equals 6378135 meters
"""

from leorbit.math import UREG, Q_

UREG.define("earthRadii = 6378135 * m")

MU_EARTH = Q_(398_600_441_800_000.0, "m**3/s**2")
"""Gravitational parameter for planet Earth (µ🜨) 

`µ🜨 = 3.986e14 m**3/s**2`
"""

SQRT_MU_EARTH = MU_EARTH**0.5
"""
`√µ🜨 = 1.996e7 m**1.5/s`
"""

RADII_EARTH = Q_("1 earthRadii")
"""Mean radius of planet Earth (R🜨). 

`R🜨 = 6378135 m`
"""