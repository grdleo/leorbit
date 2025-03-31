from typing import Callable, TypeVar, NamedTuple, TYPE_CHECKING

from algorithms.sgp4 import PosVelGCRF
from coordinates.coordinates import Coordinates
from physics.time import Time
from physics.time_interval import TimeInterval
from mathematics.units import UREG

from pint import Quantity

T = TypeVar("T")

class TimelineKey(NamedTuple):
    start_ts: float
    stop_ts: float
    
    @staticmethod
    def create(t: Time | TimeInterval) -> "TimelineKey":
        if isinstance(t, Time):
            ts = t.unixepoch
            return TimelineKey(ts, ts)
        elif isinstance(t, TimeInterval):
            return TimelineKey(t.start.unixepoch, t.stop.unixepoch)
    
    @property
    def ponctual(self) -> bool:
        return self.start_ts == self.stop_ts
    
    def intersects(self, key: "TimelineKey") -> bool:
        if self.ponctual:
            return key.start_ts <= self.start_ts <= key.stop_ts
        if key.ponctual:
            return self.start_ts <= key.start_ts <= self.stop_ts
        
        return (
            self.start_ts <= key.start_ts < self.stop_ts
            or self.start_ts < key.stop_ts <= self.stop_ts
        )
        
    def time_in(self, time: Time) -> bool:
        if self.ponctual:
            return self.start_ts == time.unixepoch
        return self.start_ts <= time.unixepoch < self.stop_ts

class Timeline[T]:
    """Assotiation of certain types of events with time.
    Events on a timeline must be exclusive (aka two events cannot happend at the same time)"""
    
    def __init__(self):
        self._pieces: dict[TimelineKey, T] = {}
        
    def append(self, event: T, epoch: Time, duration: Quantity = None):
        key = TimelineKey.create(
            TimeInterval(epoch, epoch + duration, duration*.1) if duration else epoch
        )
        for k in self._pieces.keys():
            if k.intersects(key):
                raise ValueError()
            
        self._pieces[key] = event
        
    def get(self, at: Time) -> T:
        for k in self._pieces.keys():
            if k.time_in(at):
                return self._pieces[k]
        raise KeyError()

class CoordinatesTimeline(Timeline):
    def __init__(self, 
        time_interval: TimeInterval, 
        coordinates_computer: Callable[[Time, ], Coordinates]
    ):
        self.time_interval = time_interval
        self.coordinates_computer = coordinates_computer
        self._coordinates: dict[Time, Coordinates]
        
    def append(self, *args, **kwargs):
        raise NotImplementedError("Cannot append manually to a `CoordinatesTimeline`")
    
    def get(self, at: Time) -> Coordinates:
        snapped: Time
        try:
            snapped = self.time_interval.snap_to_discretization(at)
        except ValueError:
            raise ValueError(f"Given time {at} is not in this time interval")

        coords = self._coordinates.get(snapped)
        if coords is not None:
            return coords

        coords = self.coordinates_computer(snapped)
        self._coordinates[snapped] = coords
        return coords