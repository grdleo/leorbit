from coordinates.coordinates import Coordinates
from events.timeline import CoordinatesTimeline
from physics.time import Time
from physics.time_interval import TimeInterval


class SkyObject:
    def coordinates(self, at: Time) -> Coordinates:
        raise NotImplementedError()
    
    def coordinates_timeline(self, on: TimeInterval) -> CoordinatesTimeline:
        raise NotImplementedError()