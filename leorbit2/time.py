"""Time handling"""

from datetime import datetime, timezone
from typing import cast

from leorbit2.mathematics import Dim, Scalar, Number, Quantity

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

    def __add__(self: "Time", other: Scalar[Dim.time]) -> "Time":
        try:
            assert other._dimension == Dim.Time._dim
            delta_seconds = other.base_units_value
            return self.__class__(self._unixepoch + delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __iadd__(self: "Time", other: Scalar[Dim.time]) -> None:
        try:
            assert other._dimension == Dim.Time._dim
            delta_seconds = other.base_units_value
            self._unixepoch += float(delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex
        
    def __sub__(self: "Time", other: Scalar[Dim.time]) -> "Time":
        try:
            assert other._dimension == Dim.Time._dim
            delta_seconds = other.base_units_value
            return self.__class__(self._unixepoch - delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a duration"
            ) from ex

    def __isub__(self: "Time", other: Scalar[Dim.time]) -> None:
        try:
            assert other._dimension == Dim.Time._dim
            delta_seconds = other.base_units_value
            self._unixepoch -= float(delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def delta(self: "Time", other: "Time") -> Scalar[Dim.time]:
        """Return the duration between two given `Time` objects (i.e `self - other`), as a `pint.Quantity`.

        If `other > self`, the returned duration will be negative. 
        """
        return Scalar[Dim.time](self._unixepoch - other._unixepoch)

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
    def jd(self: "Time") -> Scalar[Dim.time]:
        """Representation of this `Time` object as "Julian day (JD)", aka 
        the number of days since -4712/01/01."""
        days = (self._unixepoch / 86_400 + 2_440_587.5)
        return cast(Scalar[Dim.time], days * Quantity.day)

    @property
    def j2000(self: "Time") -> Scalar[Dim.time]:
        """Representation of this `Time` object as "Julian year (J2000)", aka 
        the number of days since 2000/01/01T12:00:00."""
        days = (self._unixepoch / 86_400 - 10_957.5)
        return cast(Scalar[Dim.time], days * Quantity.day)

    @property
    def from_mil(self: "Time") -> Scalar[Dim.time]:
        """Representation of this `Time` object as a fraction of days since 1 january 2000 00:00.

        Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
        days = (self._unixepoch / 86_400 - 10_957.5) - .5
        return cast(Scalar[Dim.time], days * Quantity.day)

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
    def stl0(self: "Time") -> Scalar[Dim.dimensionless]: # FIXME: better algorithm on the Wiki page
        """The 
        [Sideral Time](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
        (angle) of Latitude 0 at this `Time`.
        """
        d = self.j2000.magnitude("day")
        angle_rad = ((np.float128(18.697374558) + np.float128(24.06570982441908) * d) * TWELF_PI) % TWOPI
        return cast(Scalar[Dim.dimensionless], angle_rad * Quantity.rad)
    
from typing import Iterable, Self, Iterator, Optional
from physics.time import Time
from mathematics.units import Q_
from math import ceil
from algorithms.utils import humanize_duration

import numpy as np
from numpy.typing import NDArray

MIN_DURATION = Q_("1ns")

class TimeInterval:
    """A time interval between two `Time` objects. """
    def __init__(self, start: Time, stop: Time, dt=Q_("1s")):
        if not stop > start:
            raise ValueError()
        if start + dt > stop:
            raise ValueError()
        if dt < MIN_DURATION:
            raise ValueError(f"Time delta cannot be lower than minimal duration {MIN_DURATION}")
        
        self.start = start
        self.stop = stop
        self.dt = dt
        self.duration = stop.delta(start)
        steps = ceil(self.duration / dt)
        self.steps = int(steps)
    
    def __eq__(self, other: "TimeInterval") -> bool:
        if not isinstance(other, type(self)):
            return False
        return (
            self.start == other.start
            and self.stop == other.stop
            and self.dt == other.dt
        )
    
    def __repr__(self) -> str:
        return f"TimeInterval(start={self.start}, stop={self.stop}, dt={self.dt})"
    
    def __hash__(self) -> int:
        return hash(self.__repr__())

    def __iter__(self) -> Iterator[Time]:
        for i in range(self.steps):
            yield self.start + i * self.dt
    
    def duplicate(self, dt: Q_ | None = None) -> "TimeInterval":
        """Duplicates this `TimeInterval`.
        A new `dt` can be passed."""
        dt = dt if dt is not None else self.dt
        return TimeInterval(self.start, self.stop, dt)
    
    def _idx2time(self, idx: int) -> Time:
        if not isinstance(idx, int):
            raise TypeError()
        t = self.start + self.dt * idx
        if self.start <= t <= self.stop:
            return t
        raise ValueError()
    
    def _time2idx(self, time: Time) -> int:
        if not (self.start <= time <= self.stop):
            raise ValueError()
        i = (time.unixepoch - self.start.unixepoch) / self.dt.m_as("s")
        return round(i)
    
    def snap_to_discretization(self, time: Time) -> Time:
        """Returns the time closest to given time, that would be part of the interval's discretization"""
        if time <= self.start:
            return self.start
        elif time <= self.stop:
            return self.stop
        i = self._time2idx(time)
        return self._idx2time(i)
    
    def __contains__(self, t: Time | Self) -> bool:
        """Returns `True` if given `Time` or `TimeInterval` is **fully** contained in this `TimeInterval.`"""
        if isinstance(t, Time):
            return (self.start <= t <= self.stop)
        elif isinstance(t, type(self)):
            return self.start <= t.start <= t.stop <= self.stop
        raise TypeError()
    
    def intersects(self, t: Time | Self) -> bool:
        if isinstance(t, Time):
            return (self.start <= t <= self.stop)
        
        return (
            self.start <= t.start <= self.stop
            or self.start <= t.stop <= self.stop
        )
    
    def intersection(self: "TimeInterval", timeline: "TimeInterval", dt: Q_ | None = None) -> Optional["TimeInterval"]:
        """Returns the intersection of current timeline with given timeline"""
        if self.stop <= timeline.start or self.start >= timeline.stop:
            return None
        smallest, biggest = (self, timeline) if self.duration <= timeline.duration else (timeline, self)
        start_in = smallest.start in biggest
        stop_in = smallest.stop in biggest

        dt = dt if dt is not None else self.dt
        if start_in and stop_in:
            return smallest.duplicate(dt)
        elif start_in:
            return TimeInterval(smallest.start, biggest.stop, dt)
        elif stop_in:
            return TimeInterval(biggest.start, smallest.stop, dt)
        return None
    
    def progress(self, t: Time) -> float | None:
        """Returns the proportion of given time over the current timeline"""
        p = (t.delta(self.start) / self.duration).m
        if not (0 <= p <= 1):
            return None
        return p
    
    def divide(self, nb_segments: int, dt: Q_ | None = None) -> list["TimeInterval"]:
        """Divides the current TimeInterval in a given number of segments"""
        dt = dt if dt is not None else self.dt
        dur: Q_ = self.duration / nb_segments
        start = self.start
        return [TimeInterval(start + i * dur, start + (i + 1) * dur, dt) for i in range(nb_segments)]
    
    @property
    def human(self) -> str:
        """Representation of this `TimeInterval` as a human-friendly string.

        Example: `"TimeInterval: from '2024-04-11 at 17:12:49' to '2024-04-12 at 05:33:33', duration=12 hour 20 min 43 s"`
        """
        return f"TimeInterval: from '{self.start.human}' to '{self.stop.human}', duration={humanize_duration(self.duration)}"
    
    @property
    def ponctual(self) -> bool:
        """Returns `True` if this `TimeInterval` is ponctual (start == stop + MIN_DURATION)"""
        return self.start == self.stop + MIN_DURATION
    
    def to_time_stamps(self) -> NDArray:
        return np.linspace(self.start.unixepoch, self.stop.unixepoch, self.steps)
    
    @staticmethod
    def make_ponctual(time: Time) -> "TimeInterval":
        """Creates a ponctual `TimeInterval` at given `Time`"""
        return TimeInterval(time, time + MIN_DURATION, MIN_DURATION)
    
def get_intersections_timelines(first_set: Iterable[TimeInterval], second_set: Iterable[TimeInterval]) -> list[TimeInterval]:
    """Returns the intersections of the two sets of timelines"""
    all_pairs: dict[set[TimeInterval, TimeInterval], TimeInterval | None] = {}
    for t in first_set:
        for tt in second_set:
            k = t, tt
            kk = tt, t
            if k in all_pairs.keys() or kk in all_pairs.keys():
                continue
            all_pairs[k] = t.intersection(tt)
    return list(tl for tl in all_pairs.values() if tl is not None)
    
if Time.now() >= Time.fromisoformat("2100-01-01T00:00:00"):
    raise RuntimeError(f"Nobody will ever see this but considering you "
                       f"are living in the 22th century, parts of this "
                       f"code will no longer work properly. Please check "
                       f"and correct with Python 7.32")