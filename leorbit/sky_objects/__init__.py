from abc import ABC, abstractmethod
from coordinates.coordinates import Coordinates
from leorbit.coordinates.trajectory import Trajectory
from physics.time import Time
from physics.time_interval import TimeInterval


class SkyObject(ABC):
    @abstractmethod
    def coordinates(self, at: Time) -> Coordinates:
        raise NotImplementedError()
    
    @abstractmethod
    def trajectory(self, during: TimeInterval) -> Trajectory:
        raise NotImplementedError()