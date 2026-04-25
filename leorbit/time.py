"""Time handling"""

from datetime import datetime, timezone, timedelta
from typing import Annotated, Iterable, Self, Iterator, Optional
from math import ceil

import numpy as np
import numpy.typing as npt

from leorbit.mathematics import Angle, Quantity, Tensor, TensorBound, TensorKind, Time
from leorbit.utils import from_mil, humanize_duration, j2000, j2000_to_stl0, jd, stl0, unixepoch_to_j2000

MIN_DURATION = 1e-9 * Quantity.second

class Timestamp:
    """Class representing a time instant."""
    def __init__(self, unixepoch: float | int):
        """
        Parameters:
        -----------
        - `unixepoch : float` Unixepoch (timestamp) corresponding to this instant. 
        Units: seconds"""
        self._unixepoch = float(unixepoch)
    
    @property
    def unixepoch(self) -> float:
        """Unixepoch (timestamp) corresponding to this instant. 
        Units: seconds"""
        return self._unixepoch
    
    def __repr__(self) -> str:
        return f"Timestamp(unixepoch={self._unixepoch})"
    
    def __hash__(self) -> int:
        return hash(self.__repr__())

    @staticmethod
    def now() -> "Timestamp":
        """Returns a `Timestamp` object corresponding to when 
        this function was executed (aka: now)"""
        return Timestamp(
            datetime.now().timestamp()
        )

    @classmethod
    def fromisoformat(cls, iso_date: str) -> "Timestamp":
        """Creates a `Timestamp` object from a date given in ISO format as a string

        Parameters:
        -----------
        - `iso_date : str` Date in ISO format. Note: if timezone not precised, GTM is assumed
        """
        if "+" not in iso_date:
            iso_date += "+00:00"
        return cls(datetime.fromisoformat(iso_date).timestamp())

    def __eq__(self: "Timestamp", other: object) -> bool:
        if not isinstance(other, Timestamp):
            raise TypeError()
        
        return self._unixepoch == other._unixepoch
    
    def __contains__(self: "Timestamp", other: "Timestamp") -> bool:
        return self.__eq__(other)

    def __ne__(self: "Timestamp", other: object) -> bool:
        if not isinstance(other, Timestamp):
            raise TypeError()
        
        return self._unixepoch != other._unixepoch

    def __lt__(self: "Timestamp", other: "Timestamp") -> bool:
        return self._unixepoch < other._unixepoch

    def __le__(self: "Timestamp", other: "Timestamp") -> bool:
        return self._unixepoch <= other._unixepoch

    def __gt__(self: "Timestamp", other: "Timestamp") -> bool:
        return self._unixepoch > other._unixepoch

    def __ge__(self: "Timestamp", other: "Timestamp") -> bool:
        return self._unixepoch >= other._unixepoch

    def copy(self: "Timestamp") -> "Timestamp":
        """Returns a copy of this `Timestamp` object."""
        return self.__class__(self._unixepoch)
    
    def __copy__(self, *args, **kwargs) -> "Timestamp":
        return self.copy()
    
    def __deepcopy__(self, *args, **kwargs) -> "Timestamp":
        return self.copy()

    @staticmethod
    def _duration_seconds(other: Tensor | timedelta) -> float:
        if isinstance(other, timedelta):
            return other.total_seconds()

        if not isinstance(other, Tensor):
            raise TypeError("Unsupported duration type")

        if not other.check(dimension=Time, kind=TensorKind.SCALAR):
            raise TypeError("Unsupported duration type")

        delta_seconds = other.scalar.value("second")
        return float(np.asarray(delta_seconds).reshape(-1)[0])

    def __add__(self: "Timestamp", other: Tensor | timedelta) -> "Timestamp":
        try:
            delta_seconds = self._duration_seconds(other)
            return self.__class__(self._unixepoch + delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def __iadd__(self: "Timestamp", other: Tensor | timedelta) -> "Timestamp":
        try:
            self._unixepoch += self._duration_seconds(other)
            return self
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex
        
    def __sub__(self: "Timestamp", other: Tensor | timedelta) -> "Timestamp":
        try:
            delta_seconds = self._duration_seconds(other)
            return self.__class__(self._unixepoch - delta_seconds)
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a duration"
            ) from ex

    def __isub__(self: "Timestamp", other: Tensor | timedelta) -> "Timestamp":
        try:
            self._unixepoch -= self._duration_seconds(other)
            return self
        except Exception as ex:
            raise ValueError(
                f"Could not do operation with {other} and {self} since it is not a time"
            ) from ex

    def delta(self: "Timestamp", other: "Timestamp") -> Tensor:
        """Return the duration between two given `Timestamp` objects (i.e `self - other`), as a `pint.Quantity`.

        If `other > self`, the returned duration will be negative. 
        """
        return (self._unixepoch - other._unixepoch) * Quantity.second

    @property
    def isoformat(self: "Timestamp") -> str:
        """Representation of this `Timestamp` object in 
        [ISO format.](https://en.wikipedia.org/wiki/ISO_8601)"""
        return datetime.fromtimestamp(self._unixepoch, timezone.utc).isoformat()
    
    @property
    def human(self) -> str:
        """Representation of this `Timestamp` object in human readable format."""
        date = datetime.fromtimestamp(self._unixepoch, timezone.utc)
        return date.strftime("%Y-%m-%d at %H:%M:%S")

    @property
    def jd(self: "Timestamp") -> Tensor:
        """Representation of this `Timestamp` object as "Julian day (JD)", aka 
        the number of days since -4712/01/01."""
        return jd(self._unixepoch)

    @property
    def j2000(self: "Timestamp") -> Tensor:
        """Representation of this `Timestamp` object as "Julian year (J2000)", aka 
        the number of days since 2000/01/01T12:00:00."""
        return j2000(self._unixepoch)

    @property
    def from_mil(self: "Timestamp") -> Tensor:
        """Representation of this `Timestamp` object as a fraction of days since 1 january 2000 00:00.

        Taken from: https://stjarnhimlen.se/comp/ppcomp.html#3"""
        return from_mil(self._unixepoch)

    @property
    def year_day(self: "Timestamp") -> str:
        """Representation of this `Timestamp` object as a `yyddd.dddddddd` string  
        where `yy` is the last two digits of the year and 
        `ddd.dddddddd` is the fractionnal day of the year."""
        iso = self.isoformat
        full_y = iso[0:4]
        newyear = Timestamp.fromisoformat(f"{full_y}-01-01T00:00:00")
        from_newyear = self.delta(newyear)
        days = from_newyear.scalar.value("second") / 86_400
        return f"{full_y[2:4]}{days:012.8f}"

    @property
    def stl0(self: "Timestamp") -> Tensor: # FIXME: better algorithm on the Wiki page
        """The 
        [Sideral Time](https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral#Calcul_de_l'heure_sid%C3%A9rale) 
        (angle) of Latitude 0 at this `Timestamp`.
        """
        return stl0(self._unixepoch)
    
    def to_unixepoch(self) -> npt.NDArray[np.float64]:
        """Representation of this `Timestamp` object as a numpy array of unixepoch (timestamp) in seconds."""
        return np.asarray(self._unixepoch, dtype=np.float64)

class TimeInterval:
    """A time interval between two `Timestamp` objects. """
    def __init__(self, start: Timestamp, stop: Timestamp, dt=(1 * Quantity.second)):
        dt_seconds = Timestamp._duration_seconds(dt)
        dt = dt_seconds * Quantity.second

        if not stop > start:
            raise ValueError()
        if start + dt > stop:
            raise ValueError()
        if dt_seconds < float(MIN_DURATION.scalar.value("second")):
            raise ValueError(f"Time delta cannot be lower than minimal duration {MIN_DURATION}")
        
        self.start = start
        self.stop = stop
        self.dt = dt
        self.duration = stop.delta(start)
        steps = ceil((self.duration / dt).scalar.value())
        self.steps = int(steps)
    
    def __eq__(self, other: object) -> bool:
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

    def __iter__(self) -> Iterator[Timestamp]:
        for i in range(self.steps):
            yield self.start + i * self.dt
    
    def duplicate(self, dt: Tensor | None = None) -> "TimeInterval":
        """Duplicates this `TimeInterval`.
        A new `dt` can be passed."""
        dt = dt if dt is not None else self.dt
        return TimeInterval(self.start, self.stop, dt)
    
    def _idx2time(self, idx: int) -> Timestamp:
        if not isinstance(idx, int):
            raise TypeError()
        t = self.start + self.dt * idx
        if self.start <= t <= self.stop:
            return t
        raise ValueError()
    
    def _time2idx(self, time: Timestamp) -> int:
        if not (self.start <= time <= self.stop):
            raise ValueError()
        i = (time.unixepoch - self.start.unixepoch) / self.dt.scalar.value("second")
        return round(i)
    
    def snap_to_discretization(self, time: Timestamp) -> Timestamp:
        """Returns the time closest to given time, that would be part of the interval's discretization"""
        if time <= self.start:
            return self.start
        elif time <= self.stop:
            return self.stop
        i = self._time2idx(time)
        return self._idx2time(i)
    
    def __contains__(self, t: Timestamp | Self) -> bool:
        """Returns `True` if given `Timestamp` or `TimeInterval` is **fully** contained in this `TimeInterval.`"""
        if isinstance(t, Timestamp):
            return (self.start <= t <= self.stop)
        elif isinstance(t, type(self)):
            return self.start <= t.start <= t.stop <= self.stop
        raise TypeError()
    
    def intersects(self, t: Timestamp | Self) -> bool:
        if isinstance(t, Timestamp):
            return (self.start <= t <= self.stop)
        
        return (
            self.start <= t.start <= self.stop
            or self.start <= t.stop <= self.stop
        )
    
    def intersection(self: "TimeInterval", timeline: "TimeInterval", dt: Tensor | None = None) -> Optional["TimeInterval"]:
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
    
    def progress(self, t: Timestamp) -> float | None:
        """Returns the proportion of given time over the current timeline"""
        p = float((t.delta(self.start) / self.duration).scalar.value())
        if not (0 <= p <= 1):
            return None
        return p
    
    def divide(self, nb_segments: int, dt: Tensor | None = None) -> list["TimeInterval"]:
        """Divides the current TimeInterval in a given number of segments"""
        dt = dt if dt is not None else self.dt
        dur = self.duration / nb_segments
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
    
    def to_time_stamps(self) -> Tensor:
        return Tensor(self.to_unixepoch(), Time)
    
    @staticmethod
    def make_ponctual(time: Timestamp) -> "TimeInterval":
        """Creates a ponctual `TimeInterval` at given `Timestamp`"""
        return TimeInterval(time, time + MIN_DURATION, MIN_DURATION)
    
    def to_unixepoch(self) -> npt.NDArray[np.float64]:
        """Representation of this `TimeInterval` object as a numpy array of unixepoch (timestamp) in seconds."""
        return np.linspace(self.start.unixepoch, self.stop.unixepoch, self.steps, dtype=np.float64)

Timeline = TimeInterval
    
def get_intersections_timelines(first_set: Iterable[TimeInterval], second_set: Iterable[TimeInterval]) -> list[TimeInterval]:
    """Returns the intersections of the two sets of timelines"""
    all_pairs: dict[tuple[TimeInterval, TimeInterval], TimeInterval | None] = {}
    for t in first_set:
        for tt in second_set:
            k = t, tt
            kk = tt, t
            if k in all_pairs.keys() or kk in all_pairs.keys():
                continue
            all_pairs[k] = t.intersection(tt)
    return list(tl for tl in all_pairs.values() if tl is not None)
    
if Timestamp.now() >= Timestamp.fromisoformat("2100-01-01T00:00:00"):
    raise RuntimeError(f"Nobody will ever see this but considering you "
                       f"are living in the 22th century, parts of this "
                       f"code will no longer work properly. Please check "
                       f"and correct with Python 7.32")