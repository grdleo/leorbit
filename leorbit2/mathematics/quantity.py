
from math import pi
from leorbit2.mathematics.dimensions import D, registered_units
from leorbit2.mathematics.scalar import Scalar


class QuantityMeta(type):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """

    def __getattr__(cls, name: str) -> Scalar:
        """Resolve a unit name into its corresponding scalar quantity."""
        try:
            dim, factor = registered_units()[name]
            return Scalar[dim].new(factor)
        except KeyError:
            raise ValueError(f"No unit named '{name}'")   

class Quantity(metaclass=QuantityMeta):
    """Convenience namespace for creating unit-scaled scalar values."""

    @classmethod
    def get(cls, value: str) -> Scalar:
        """Return the scalar unit associated with ``value``."""
        return cls.__getattr__(value)

    # ANGLES

    rad: Scalar[D.Angle]
    """radians"""

    deg: Scalar[D.Angle]
    """degrees"""

    # DISTANCES

    meter: Scalar[D.Length]
    """meter"""

    kilo_meter: Scalar[D.Length]
    """kilometer"""

    radii_earth: Scalar[D.Length]
    """Mean radius of planet Earth (R🜨). 

    `R🜨 = 6378135 m`
    """

    radii_sun: Scalar[D.Length]
    """Mean radius of Sun (R☉). 

    `R☉ = 6.957e8 m`
    """

    # DURATIONS

    second: Scalar[D.Time]
    """second"""

    minute: Scalar[D.Time]
    """minute"""

    hour: Scalar[D.Time]
    """hour"""

    day: Scalar[D.Time]
    """day"""

    month: Scalar[D.Time]
    """month (30 days)"""

    year: Scalar[D.Time]
    """year (365 days)"""