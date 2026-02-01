"""Time handling"""

from datetime import datetime, timezone
from typing import cast

from leorbit2.mathematics import Dimensionless, TimeDim, Scalar, Number, DimensionRegister, Quantity

import numpy as np
TWOPI = 2 * np.pi
TWELF_PI = np.pi / 12

class Time:
    """Class representing a time instant."""
    def __init__(self, unixepoch: Number):
        """
        Parameters:
        -----------
        - `unixepoch : float` Unixepoch (timestamp) corresponding to this instant. 
        Units: seconds"""
        self._unixepoch = float(unixepoch)
    
    @property
    def unixepoch(self) -> Number:
        """Unixepoch (timestamp) corresponding to this instant. 
        Units: seconds"""
        return self._unixepoch
    
    def __repr__(self) -> str:
        return f"Time(unixepoch={self._unixepoch})"
    
    def __hash__(self) -> int:
        return hash(self.__repr__())

    @staticmethod
    def now() -> "Time":
        """Returns a `Time` object corresponding to when 
        this function was executed (aka: now)"""
        return Time(
            datetime.now().timestamp()
        )

    @classmethod
    def fromisoformat(cls, iso_date: str) -> "Time":
        """Creates a `Time` object from a date given in ISO format as a string

        Parameters:
        -----------
        - `iso_date : str` Date in ISO format. Note: if timezone not precised, GTM is assumed
        """
        if "+" not in iso_date:
            iso_date += "+00:00"
        return cls(datetime.fromisoformat(iso_date).timestamp())

    def __eq__(self: "Time", other: object) -> bool:
        if not isinstance(other, Time):
            raise TypeError()
        
        return self._unixepoch == other._unixepoch
    
    def __contains__(self: "Time", other: "Time") -> bool:
        return self.__eq__(other)

    def __ne__(self: "Time", other: object) -> bool:
        if not isinstance(other, Time):
            raise TypeError()
        
        return self._unixepoch != other._unixepoch

    def __lt__(self: "Time", other: "Time") -> bool:
        return self._unixepoch < other._unixepoch

    def __le__(self: "Time", other: "Time") -> bool:
        return self._unixepoch <= other._unixepoch

    def __gt__(self: "Time", other: "Time") -> bool:
        return self._unixepoch > other._unixepoch

    def __ge__(self: "Time", other: "Time") -> bool:
        return self._unixepoch >= other._unixepoch

    def copy(self: "Time") -> "Time":
        """Returns a copy of this `Time` object."""
        return self.__class__(self._unixepoch)
    
    def __copy__(self, *args, **kwargs) -> "Time":
        return self.copy()
    
    def __deepcopy__(self, *args, **kwargs) -> "Time":
        return self.copy()

    def __add__(self: "Time", other: Scalar[TimeDim]) -> "Time":
        try:
            assert other._dimension == DimensionRegister.time
            delta_seconds = other.base_units_value
            return self.__class__(self._unixepoch + delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __iadd__(self: "Time", other: Scalar[TimeDim]) -> None:
        try:
            assert other._dimension == DimensionRegister.time
            delta_seconds = other.base_units_value
            self._unixepoch += float(delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex
        
    def __sub__(self: "Time", other: Scalar[TimeDim]) -> "Time":
        try:
            assert other._dimension == DimensionRegister.time
            delta_seconds = other.base_units_value
            return self.__class__(self._unixepoch - delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a duration"
            ) from ex

    def __isub__(self: "Time", other: Scalar[TimeDim]) -> None:
        try:
            assert other._dimension == DimensionRegister.time
            delta_seconds = other.base_units_value
            self._unixepoch -= float(delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def delta(self: "Time", other: "Time") -> Scalar[TimeDim]:
        """Return the duration between two given `Time` objects (i.e `self - other`), as a `pint.Quantity`.

        If `other > self`, the returned duration will be negative. 
        """
        return Scalar[TimeDim](self._unixepoch - other._unixepoch)

    @property
    def isoformat(self: "Time") -> str:
        """Representation of this `Time` object in 
        [ISO format.](https://en.wikipedia.org/wiki/ISO_8601)"""
        return datetime.fromtimestamp(self._unixepoch, timezone.utc).isoformat()
    
    @property
    def human(self) -> str:
        """Representation of this `Time` object in human readable format."""
        date = datetime.fromtimestamp(self._unixepoch, timezone.utc)
        return date.strftime("%Y-%m-%d at %H:%M:%S")

    @property
    def jd(self: "Time") -> Scalar[TimeDim]:
        """Representation of this `Time` object as "Julian day (JD)", aka 
        the number of days since -4712/01/01."""
        days = (self._unixepoch / 86_400 + 2_440_587.5)
        return cast(Scalar[TimeDim], days * Quantity.day)

    @property
    def j2000(self: "Time") -> Scalar[TimeDim]:
        """Representation of this `Time` object as "Julian year (J2000)", aka 
        the number of days since 2000/01/01T12:00:00."""
        days = (self._unixepoch / 86_400 - 10_957.5)
        return cast(Scalar[TimeDim], days * Quantity.day)

    @property
    def from_mil(self: "Time") -> Scalar[TimeDim]:
        """Representation of this `Time` object as a fraction of days since 1 january 2000 00:00.

        Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
        days = (self._unixepoch / 86_400 - 10_957.5) - .5
        return cast(Scalar[TimeDim], days * Quantity.day)

    @property
    def year_day(self: "Time") -> str:
        """Representation of this `Time` object as a `yyddd.dddddddd` string  
        where `yy` is the last two digits of the year and 
        `ddd.dddddddd` is the fractionnal day of the year."""
        iso = self.isoformat
        full_y = iso[0:4]
        newyear = Time.fromisoformat(f"{full_y}-01-01T00:00:00")
        from_newyear = self.delta(newyear)
        days = from_newyear.base_units_value / 86_400
        return f"{full_y[2:4]}{days:012.8f}"

    @property
    def stl0(self: "Time") -> Scalar[Dimensionless]: # FIXME: better algorithm on the Wiki page
        """The 
        [Sideral Time](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
        (angle) of Latitude 0 at this `Time`.
        """
        d = self.j2000.base_units_value / 86_400
        angle_rad = ((np.float128(18.697374558) + np.float128(24.06570982441908) * d) * TWELF_PI) % TWOPI
        return cast(Scalar[Dimensionless], angle_rad * Quantity.rad)
    
if Time.now() >= Time.fromisoformat("2100-01-01T00:00:00"):
    raise RuntimeError(f"Nobody will ever see this but considering you "
                       f"are living in the 22th century, parts of this "
                       f"code will no longer work properly. Please check "
                       f"and correct with Python 7.32")