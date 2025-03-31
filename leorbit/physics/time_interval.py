from typing import Iterable, Self, Iterator, Optional
from physics.time import Time
from mathematics import Q_
from math import ceil
from algorithms.utils import humanize_duration

class TimeInterval:
    """A time interval between two `Time` objects. """
    def __init__(self, start: Time, stop: Time, dt=Q_("1s")):
        if not stop > start:
            raise ValueError()
        if start + dt > stop:
            raise ValueError()
        if dt <= 0:
            raise ValueError()
        
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
        for i in range(self.steps + 1):
            yield self.start + i * self.dt
    
    def duplicate(self, dt: Q_ | None = None) -> "TimeInterval":
        """Duplicates this `TimeInterval`.
        A new `dt` can be passed."""
        dt = dt if dt is not None else self.dt
        return TimeInterval(self.start, self.stop, dt)
    
    def _idx2time(self, idx: int) -> Time:
        if not isinstance(idx, int):
            raise ValueError()
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
    
    def intersection(self, timeline: "TimeInterval", dt: Q_ | None = None) -> Optional["TimeInterval"]:
        """Returns the intersection of current timeline with given timeline"""
        if self.stop <= timeline.start or self.start >= timeline.stop:
            return None
        smallest, biggest = (self, timeline) if self.duration <= timeline.duration else (timeline, self)
        start_in = biggest.contains(smallest.start)
        stop_in = biggest.contains(smallest.stop)

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