
from math import pi
from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.scalar import Scalar


class QuantityMeta(type):
    def __getattr__(cls, name: str) -> Scalar:
        if name == "rad":
            return Scalar[D.Angle].new(1)
        elif name == "deg":
            return Scalar[D.Angle].new(pi / 180)
        elif name == "m":
            return Scalar[D.Length].new(D.Length.meter)
        elif name == "km":
            return Scalar[D.Length].new(D.Length.kilo_meter)
        elif name == "radii_earth":
            return Scalar[D.Length].new(6378135)
        elif name == "radii_sun":
            return Scalar[D.Length].new(6.957e8)
        elif name == "s":
            return Scalar[D.Time].new(D.Time.second)
        elif name == "min":
            return Scalar[D.Time].new(D.Time.minute)
        elif name == "hour":
            return Scalar[D.Time].new(D.Time.hour)
        elif name == "day":
            return Scalar[D.Time].new(D.Time.day)
        elif name == "month":
            return Scalar[D.Time].new(D.Time.month)
        elif name == "year":
            return Scalar[D.Time].new(D.Time.year)
        
        raise ValueError(f"No unit named '{name}'")
        

class Quantity(metaclass=QuantityMeta):
    # ANGLES

    rad: Scalar[D.Angle]
    """radians"""

    deg: Scalar[D.Angle]
    """degrees"""

    # DISTANCES

    m: Scalar[D.Length]
    """meter"""

    km: Scalar[D.Length]
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

    s: Scalar[D.Time]
    """second"""

    min: Scalar[D.Time]
    """minute"""

    hour: Scalar[D.Time]
    """hour"""

    day: Scalar[D.Time]
    """day"""

    month: Scalar[D.Time]
    """month (30 days)"""

    year: Scalar[D.Time]
    """year (365 days)"""