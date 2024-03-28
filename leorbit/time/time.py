"""Time handling"""

from datetime import datetime, timezone
from typing import Self

from leorbit.math import TWELF_PI, TWOPI
from leorbit.math import Q_

class Time:
    """Class representing a time instant."""
    def __init__(self, unixepoch: float):
        """
        Parameters:
        -----------
        - `unixepoch : float` Unixepoch (timestamp) corresponding to this instant. 
        Units: seconds"""
        self._unixepoch = unixepoch
    
    @property
    def unixepoch(self) -> float:
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
    def fromisoformat(cls: "Time", iso_date: str) -> "Time":
        """Creates a `Time` object from a date given in ISO format as a string

        Parameters:
        -----------
        - `iso_date : str` Date in ISO format. Note: if timezone not precised, GTM is assumed
        """
        if "+" not in iso_date:
            iso_date += "+00:00"
        return cls(datetime.fromisoformat(iso_date).timestamp())

    def __eq__(self: Self, other: Self) -> bool:
        return self.unixepoch == other.unixepoch

    def __ne__(self: Self, other: Self) -> bool:
        return self.unixepoch != other.unixepoch

    def __lt__(self: Self, other: Self) -> bool:
        return self.unixepoch < other.unixepoch

    def __le__(self: Self, other: Self) -> bool:
        return self.unixepoch <= other.unixepoch

    def __gt__(self: Self, other: Self) -> bool:
        return self.unixepoch > other.unixepoch

    def __ge__(self: Self, other: Self) -> bool:
        return self.unixepoch >= other.unixepoch

    def __repr__(self: Self) -> str:
        return f"<Time unixepoch={self.unixepoch}>"

    def __hash__(self: Self) -> int:
        return hash(self.unixepoch)

    def copy(self: Self) -> Self:
        """Returns a copy of this `Time` object."""
        return type(self)(self.unixepoch)

    def __add__(self: Self, other: Q_) -> Self:
        try:
            return type(self)(self.unixepoch + other.m_as("s"))
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __iadd__(self: Self, other: Q_) -> None:
        try:
            self._unixepoch += other.m_as("s")
            return self
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __isub__(self: Self, other: Q_) -> None:
        try:
            self._unixepoch -= other.m_as("s")
            return self
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __sub__(self: Self, other: Q_) -> Self:
        try:
            return type(self)(self.unixepoch - other.m_as("s"))
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a duration"
            ) from ex

    def delta(self: Self, other: Self) -> Q_:
        """Return the duration between two given `Time` objects, as a `pint.Quantity`.
        The operator `~` overloads this function.

        If `other > self`, the returned duration will be negative. 
        """
        return Q_(self.unixepoch - other.unixepoch, "s")
    
    def __invert__(self: Self, other: Self) -> Q_: # ~
        """Operator for `delta` method"""
        return self.delta(other)

    @property
    def isoformat(self: Self) -> str:
        """Representation of this `Time` object in 
        [ISO format.](https://en.wikipedia.org/wiki/ISO_8601)"""
        return datetime.fromtimestamp(self.unixepoch, timezone.utc).isoformat()
    
    @property
    def human(self) -> str:
        """Representation of this `Time` object in human readable format."""
        date = datetime.fromtimestamp(self.unixepoch, timezone.utc)
        return date.strftime("%Y-%m-%d at %H:%M:%S")

    @property
    def jd(self: Self) -> Q_:
        """Representation of this `Time` object as "Julian day (JD)", aka 
        the number of days since -4712/01/01."""
        return Q_(self.unixepoch / 86400 + 2440587.5, "day")

    @property
    def j2000(self: Self) -> Q_:
        """Representation of this `Time` object as "Julian year (J2000)", aka 
        the number of days since 2000/01/01T12:00:00."""
        return Q_(self.unixepoch / 86400 - 10957.5, "day")

    @property
    def from_mil(self: Self) -> Q_:
        """Representation of this `Time` object as a fraction of days since 1 january 2000 00:00.

        Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
        return self.j2000 - Q_(0.5, "day")

    @property
    def year_day(self: Self) -> str:
        """Representation of this `Time` object as a `yyddd.ddddddd` string  
        where `yy` is the last two digits of the year and 
        `ddd.ddddddd` is the fractionnal day of the year."""
        iso = self.isoformat
        y = iso[0:4]
        newyear = Time.fromisoformat(f"{y}-01-01T00:00:00")
        from_newyear = self.delta(newyear)
        d_full = round(from_newyear.to("day").m, 8)
        d_int = int(d_full)
        d_dec = str(d_full - d_int)[2:]
        d = f"{str(d_int).rjust(3, '0')}.{d_dec.ljust(8, '0')}"
        return f"{y[2:4]}{d}"

    @property
    def stl0(self: Self) -> Q_: # FIXME: better algorithm on the Wiki page
        """The 
        [Sideral Time](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
        (angle) of Latitude 0 at this `Time`.
        """
        return Q_(
            ((18.697374558 + 24.06570982441908 * self.j2000.m) * TWELF_PI) % TWOPI, "rad"
        )
    
if Time.now() >= Time.fromisoformat("2100-01-01T00:00:00"):
    raise RuntimeError(f"Nobody will ever see this but considering you "
                       f"are living in the 22th century, parts of this "
                       f"code will no longer work properly. Please check "
                       f"and correct with Python 7.32")
