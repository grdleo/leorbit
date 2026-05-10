from abc import ABC, abstractmethod
from typing import Annotated, Any

import numpy as np

from leorbit.coordinates import GPS, Trajectory
from leorbit.mathematics import Angle, Quantity, Tensor, TensorBound, TensorKind, sin
from leorbit.time import TimeInterval, TimeIntervalSet, Timestamp

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
    """Time-indexed container for arrays computed on a fixed interval."""

    def __init__(self, interval: TimeInterval, **values: npt.NDArray):
        """Store named arrays sampled over ``interval``.

        All arrays must have exactly ``interval.steps`` elements.
        """
        for arr in values.values():
            if len(arr) != interval.steps:
                raise ValueError(f"All value arrays must have the same length as the interval, expected {interval.steps}, got {len(arr)}")
        
        self.interval = interval
        self._mapped_values = values

    def get_values(self, value_name: str) -> npt.NDArray:
        """Return the full sampled array registered under ``value_name``."""
        values = self._mapped_values.get(value_name)
        if values is None:
            raise ValueError(f"Value name {value_name} not found in TimeMap")
        
        return values

    def get_value(self, epoch: Timestamp, value_name: str) -> Any:
        """Return a single mapped value at ``epoch`` for ``value_name``."""
        if epoch not in self.interval:
            raise ValueError(f"Epoch {epoch} not contained in TimeMap interval {self.interval}")
        
        values = self.get_values(value_name)
        idx = self.interval._time2idx(epoch)
        return values[idx]

class Event(ABC):
    def __init__(self, trajectory: Trajectory):
        """Initialize an event evaluator for a trajectory."""
        self.trajectory = trajectory
        self._time_map = self._compute()

    @property
    def timeline(self) -> TimeInterval:
        """Sampling timeline used to compute this event."""
        return self.trajectory.interval
    
    @abstractmethod
    def _compute(self) -> TimeMap:
        """Compute the event for the given trajectory and return a TimeMap of the event values"""
        ...

class VisibleFromEarthLocationEvent(Event):
    """Visibility event of a trajectory from a fixed Earth observer location."""

    _KW_VISIBLE = "visible"

    def __init__(self, 
        trajectory: Trajectory, 
        gps_observer: GPS,
        altitude_angle_min: Annotated[Tensor, TensorBound(dimension=Angle, kind=TensorKind.SCALAR)] = 0 * Quantity.radian
    ):
        """Build a visibility event from observer GPS coordinates.

        ``altitude_angle_min`` is kept for API compatibility and future
        visibility thresholds.
        """
        self.gps_observer = gps_observer
        self.altitude_angle_min = altitude_angle_min

        if gps_observer.altitude < 0 * Quantity.meter or gps_observer.altitude > 10 * Quantity.kilo_meter:
            print(
                f"Warning: Observer GPS altitude {gps_observer.altitude} is out of typical range for Earth's surface! "
                "This may lead to inaccurate results."
            )

        super().__init__(trajectory)
    
    def _compute(self) -> TimeMap:
        """Compute per-step visibility booleans for the underlying trajectory."""
        local_frame = self.gps_observer.earth_local_frame
        local_pos = self.trajectory.trajectory_pos(local_frame)

        visible = (
            local_pos.vector3.normalized().vector3.z.raw_data_array()
            > sin(self.altitude_angle_min).scalar.raw_data_array()
        )

        return TimeMap(
            self.timeline,
            **{
                self._KW_VISIBLE: visible
            }
        )

    @property
    def visible_intervals(self) -> TimeIntervalSet:
        """Contiguous intervals where the object is visible from the observer."""
        intervals = _truth_array_to_time_intervals(
            self._time_map.get_values(self._KW_VISIBLE),
            self.timeline
        )
        return TimeIntervalSet(intervals)
    
    @property
    def not_visible_intervals(self) -> TimeIntervalSet:
        """Contiguous intervals where the object is below the local horizon."""
        intervals = _truth_array_to_time_intervals(
            ~self._time_map.get_values(self._KW_VISIBLE),
            self.timeline
        )
        return TimeIntervalSet(intervals)
    