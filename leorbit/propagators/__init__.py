from abc import ABC, abstractmethod

from coordinates.coordinates import Coordinates
from coordinates.representations.elements import OrbitalElements
from leorbit.coordinates.trajectory import Trajectory
from physics.time import Time
from physics.time_interval import TimeInterval


class Propagator(ABC):
    """Algorithm to propagate given orbital elements at given time"""
    
    def __init__(self, elements: OrbitalElements):
        self.elements = elements
    
    @abstractmethod
    def propagate(self, to: Time) -> Coordinates:
        ...
        
    @abstractmethod
    def propagate_timeline(self, on: TimeInterval) -> Trajectory:
        ...