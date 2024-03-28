from typing import Iterable, Self, Iterator, Optional
from leorbit.time import Time
from leorbit.math import Q_
from math import ceil
from leorbit.algorithms.utils import humanize_duration

class Timeline:
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
    
    def __eq__(self, other: "Timeline") -> bool:
        if not isinstance(other, type(self)):
            return False
        return (
            self.start == other.start
            and self.stop == other.stop
            and self.dt == other.dt
        )
    
    def __repr__(self) -> str:
        return f"Timeline(start={self.start}, stop={self.stop}, dt={self.dt})"
    
    def __hash__(self) -> int:
        return hash(self.__repr__())

    def __iter__(self) -> Iterator[Time]:
        for i in range(self.steps + 1):
            yield self.start + i * self.dt
    
    def duplicate(self, dt: Q_ | None = None) -> "Timeline":
        """Duplicates this `Timeline`.
        A new `dt` can be passed."""
        dt = dt if dt is not None else self.dt
        return Timeline(self.start, self.stop, dt)
    
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
    
    def contains(self, t: Time | Self) -> bool:
        """Returns `True` if given `Time` or `Timeline` is **fully** contained in this `Timeline.`"""
        if isinstance(t, Time):
            return (self.start <= t <= self.stop)
        elif isinstance(t, type(self)):
            return self.start <= t.start <= t.stop <= self.stop
        raise TypeError()
    
    def intersection(self, timeline: "Timeline", dt: Q_ | None = None) -> Optional["Timeline"]:
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
            return Timeline(smallest.start, biggest.stop, dt)
        elif stop_in:
            return Timeline(biggest.start, smallest.stop, dt)
        return None
    
    def progress(self, t: Time) -> float | None:
        """Returns the proportion of given time over the current timeline"""
        p = (t.delta(self.start) / self.duration).m
        if not (0 <= p <= 1):
            return None
        return p
    
    def divide(self, nb_segments: int, dt: Q_ | None = None) -> list["Timeline"]:
        """Divides the current Timeline in a given number of segments"""
        dt = dt if dt is not None else self.dt
        dur: Q_ = self.duration / nb_segments
        start = self.start
        return [Timeline(start + i * dur, start + (i + 1) * dur, dt) for i in range(nb_segments)]
    
    @property
    def human(self) -> str:
        """Representation of this `Timeline` as a human-friendly string.

        Example: `"Timeline: from '2024-04-11 at 17:12:49' to '2024-04-12 at 05:33:33', duration=12 hour 20 min 43 s"`
        """
        return f"Timeline: from '{self.start.human}' to '{self.stop.human}', duration={humanize_duration(self.duration)}"

def get_intersections_timelines(first_set: Iterable[Timeline], second_set: Iterable[Timeline]) -> list[Timeline]:
    """Returns the intersections of the two sets of timelines"""
    all_pairs: dict[set[Timeline, Timeline], Timeline | None] = {}
    for t in first_set:
        for tt in second_set:
            k = t, tt
            kk = tt, t
            if k in all_pairs.keys() or kk in all_pairs.keys():
                continue
            all_pairs[k] = t.intersection(tt)
    return list(tl for tl in all_pairs.values() if tl is not None)