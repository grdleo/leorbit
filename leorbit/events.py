from abc import ABC, abstractmethod
from typing import Any, cast

import numpy as np

from leorbit.coordinates import GPS
from leorbit.mathematics import Angle, Dimless, Quantity, Scalar, ScalarArray, acos
from leorbit.propagator import Trajectory
from leorbit.time import TimeInterval, Timestamp

import numpy.typing as npt

from leorbit.utils import truth_array_to_indices_intervals

def _truth_array_to_time_intervals(truth_array: npt.NDArray, timeline: TimeInterval) -> list[TimeInterval]:
    return [
        TimeInterval(
            timeline._idx2time(i0),
            timeline._idx2time(i1),
            timeline.dt
        )
        for i0, i1 in truth_array_to_indices_intervals(truth_array)
    ]


class TimeMap:
    def __init__(self, interval: TimeInterval, **values: npt.NDArray):
        for arr in values.values():
            if len(arr) != interval.steps:
                raise ValueError(f"All value arrays must have the same length as the interval, expected {interval.steps}, got {len(arr)}")
        
        self.interval = interval
        self._mapped_values = values

    def get_values(self, value_name: str) -> npt.NDArray:
        values = self._mapped_values.get(value_name)
        if values is None:
            raise ValueError(f"Value name {value_name} not found in TimeMap")
        
        return values

    def get_value(self, epoch: Timestamp, value_name: str) -> Any:
        if epoch not in self.interval:
            raise ValueError(f"Epoch {epoch} not contained in TimeMap interval {self.interval}")
        
        values = self.get_values(value_name)
        idx = self.interval._time2idx(epoch)
        return values[idx]

class Event(ABC):
    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory
        self._time_map = self.compute()

    @property
    def timeline(self) -> TimeInterval:
        return self.trajectory.interval
    
    @abstractmethod
    def compute(self) -> TimeMap:
        """Compute the event for the given trajectory and return a TimeMap of the event values"""
        ...

class VisibleFromEarthLocationEvent(Event):
    def __init__(self, 
        trajectory: Trajectory, 
        gps_observer: GPS,
        altitude_angle_min: Scalar[Angle] = 0 * Quantity.radian
    ):
        self.gps_observer = gps_observer
        self.altitude_angle_min = altitude_angle_min

        if gps_observer.altitude < 0 or gps_observer.altitude > 10 * Quantity.kilo_meter:
            print(
                f"Warning: Observer GPS altitude {gps_observer.altitude} is out of typical range for Earth's surface! "
                "This may lead to inaccurate results."
            )

        super().__init__(trajectory)
    
    def compute(self) -> TimeMap:
        local_frame = self.gps_observer.earth_local_frame
        local_pos = self.trajectory.trajectory_pos(local_frame)
    
        visible = local_pos.z.get_raw_array("meter") > 0 # visible if satellite is above the horizon

        return TimeMap(
            self.timeline,
            visible=visible
        )

    @property
    def visible_intervals(self) -> list[TimeInterval]:
        return _truth_array_to_time_intervals(
            self._time_map.get_values("visible"),
            self.timeline
        )
    
    @property
    def not_visible_intervals(self) -> list[TimeInterval]:
        return _truth_array_to_time_intervals(
            ~self._time_map.get_values("visible"),
            self.timeline
        )
    
