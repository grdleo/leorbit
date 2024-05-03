"""Handling for astro events: compute when specific conditions are met on a satellite"""

from abc import abstractmethod
from leorbit.algorithms.utils import apparent_magnitude
from leorbit.math.coordinate import GPSCoordinates
from leorbit.math.vector import Vec3
from leorbit.math import Q_
from leorbit.orbit.objects import AstroObject, Satellite, Sun
from leorbit.simulation.utils import compute_magnitude
from leorbit.time import Time
from leorbit.time import Timeline
from leorbit.orbit import SUN, MOON

from typing import Self

class AstroEvent:
    pass

class AstroEvent:
    """Object representing an astronomical event, 
    aka time periods (`Timeline`) where a certain condition is met."""
    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def predicate(self, epoch: Time) -> bool:
        """The condition of this event.
        Must return `True` if condition is met at given epoch."""
        raise NotImplementedError()
    
    def compute_intervals(self, over: Timeline, minimal_duration: Q_ = Q_("0s"), print_progress: bool = False) -> list[Timeline]:
        """Returns all the `Timeline` objects when predicate of this event is `True`
        over the given timeline.

        Parameters:
        -----------
        - `over: Timeline` The timeline to compute over
        - `minimal_duration: Q_` Intervals that validates the event's predicate have to be at least this long. 
            If not, the interval will not be considered.
        - `print_progress: bool` If `True`, prints progress to the console.
        """
        start = None
        intervals: list[Timeline] = []
        last_percent = None
        for t in over:
            p = self.predicate(t)

            if print_progress:
                percent = int(over.progress(t) * 100)
                if percent != last_percent:
                    last_percent = percent
                    print(f"Progress: {percent}%")
            
            if start is None and not p:
                continue
            elif start is not None and p:
                continue
            elif start is None and p:
                start = t
            elif start is not None and not p:
                intervals.append(
                    Timeline(start, t, over.dt)
                )
                start = None
            
        return [tl for tl in intervals if tl.duration >= minimal_duration]
    
    def intersection(event1, event2: Self, name: str = None) -> Self:
        """Creates a new `AstroEvent` that corresponds to the intersection of both events.

        Boolean operator `&` overloads this function.
        """
        if not issubclass(event2.__class__, AstroEvent):
            raise TypeError()
        new_event = AstroEvent(f"({event1.name} and {event2.name})" if not isinstance(name, str) else name)
        def new_event_predicate(epoch: Time) -> bool:
            return event1.predicate(epoch) and event2.predicate(epoch)
        new_event.predicate = new_event_predicate
        return new_event
    
    def __and__(self, other: Self) -> Self:
        """Creates a new `AstroEvent` that corresponds to the intersection of both events"""
        return self.intersection(other)
    
    def union(event1, event2: Self, name: str = None) -> Self:
        """Creates a new `AstroEvent` that corresponds to the union of both events.
        
        Boolean operator `|` overloads this function.
        """
        if not issubclass(event2.__class__, AstroEvent):
            raise TypeError()
        new_event = AstroEvent(f"({event1.name} or {event2.name})" if not isinstance(name, str) else name)
        def new_event_predicate(epoch: Time) -> bool:
            return event1.predicate(epoch) or event2.predicate(epoch)
        new_event.predicate = new_event_predicate
        return new_event
    
    def __or__(self, other: Self) -> Self:
        """Creates a new `AstroEvent` that corresponds to the union of both events"""
        return self.union(other)

class InNorthenAuroraArea(AstroEvent):
    """This event corresponds to when the given satellite is located
    in a zone around the North pole where it is likely to encounter polar lights"""
    def __init__(self, sat: AstroObject):
        super().__init__("In northen hemisphere")
        self.sat = sat
    
    def predicate(self, epoch: Time) -> bool:
        return 65 <= self.sat.coordinates(epoch).to_gps().latitude_deg <= 75

class VisibleFromLocation(AstroEvent):
    """This event corresponds to when the given satellite is visible
    from a specified location on Earth"""
    def __init__(self, sat: AstroObject, gps_location: GPSCoordinates, min_altitude_visi: Q_ = Q_("10°")):
        assert isinstance(min_altitude_visi, Q_)
        assert min_altitude_visi.check("rad")

        super().__init__(f"Visibility from {gps_location.custom_repr()}")

        self.sat = sat
        self.gps = gps_location
        self.local_frame = gps_location.earth_local_frame()
        self.min_alti = min_altitude_visi.m_as("°") # [°]
    
    def predicate(self, epoch: Time) -> bool:
        sat_coord = self.sat.coordinates(epoch)
        hor = self.local_frame.get_horizontal_coordinates(sat_coord)
        return hor.altitude_deg >= self.min_alti

class NightTime(AstroEvent):
    """This event corresponds to night time at given location."""
    def __init__(self, gps_location: GPSCoordinates, sun: Sun = SUN):
        super().__init__(f"Night time at {gps_location.custom_repr()}")

        self.sun = sun
        self.gps = gps_location
        self.local_frame = gps_location.earth_local_frame()
    
    def predicate(self, epoch: Time) -> bool:
        sun_coord = self.sun.coordinates(epoch)
        hor = self.local_frame.get_horizontal_coordinates(sun_coord)
        return hor.altitude_deg <= 0

class BrightnessBelow(AstroEvent):
    """This event corresponds to when sat has brightness (magnitude) below given value"""
    def __init__(self, sat: Satellite, observer: GPSCoordinates, below: float, sun: Sun = SUN):
        super().__init__(f"Brightness below {below}")

        self.sat_std_mag = sat.std_mag
        if self.sat_std_mag is None:
            raise RuntimeError()
        self.sun = sun
        self.sat = sat
        self.local_frame = observer.earth_local_frame()
        self.below = below
    
    def predicate(self, epoch: Time) -> bool:
        mag = compute_magnitude(epoch, self.local_frame, self.sat)
        return mag <= self.below