import numpy as np
import pytest

from leorbit.coordinates import GPS, Trajectory
from leorbit.events import TimeMap, VisibleFromEarthLocationEvent
from leorbit.frames import AbsoluteFrame
from leorbit.mathematics import U, Tensor, scalar
from leorbit.time import TimeInterval, Timestamp



def _make_interval() -> TimeInterval:
    t0 = Timestamp.fromisoformat("2024-01-01T00:00:00")
    return TimeInterval(t0, t0 + 6 * U.second, 1 * U.second)


def test_timemap_rejects_wrong_array_lengths():
    interval = _make_interval()
    with pytest.raises(ValueError):
        TimeMap(interval, visible=np.asarray([True, False]))


def test_timemap_get_values_and_value_lookup():
    interval = _make_interval()
    visible = np.asarray([False, False, True, True, False, False], dtype=bool)
    time_map = TimeMap(interval, visible=visible)

    assert np.array_equal(time_map.get_values("visible"), visible)
    assert bool(time_map.get_value(interval._idx2time(2), "visible")) is True
    assert bool(time_map.get_value(interval._idx2time(0), "visible")) is False


def test_visible_from_earth_location_event_full():
    interval = _make_interval()
    gps_paris = GPS(
        longitude=scalar(2.333333).with_units(U.degree), 
        latitude=scalar(48.866667).with_units(U.degree), 
        altitude=scalar(0).with_units(U.meter),
        epoch=interval.start,
    )

    local_frame = gps_paris.earth_local_frame

    # Visibility in this event is based on local-frame z > 0.
    z = np.asarray([-1.0, -0.2, 3.0, 4.0, -2.0, -1.0], dtype=np.float64)
    pos_local = Tensor(
        np.asarray(
            [
                np.zeros_like(z),
                np.zeros_like(z),
                z,
            ],
            dtype=np.float64,
        ),
        U.meter,
    )

    trajectory = Trajectory(interval, local_frame, pos_local)
    event = VisibleFromEarthLocationEvent(trajectory, gps_paris)

    visible = event._time_map.get_values("visible")
    assert np.array_equal(visible, np.asarray([False, False, True, True, False, False], dtype=bool))

    visible_intervals = event.visible_intervals
    assert len(visible_intervals) == 1
    assert visible_intervals[0].start == interval._idx2time(2)
    assert visible_intervals[0].stop == interval._idx2time(3)

    not_visible_intervals = event.not_visible_intervals
    assert len(not_visible_intervals) == 2
    assert not_visible_intervals[0].start == interval._idx2time(0)
    assert not_visible_intervals[0].stop == interval._idx2time(1)
    assert not_visible_intervals[1].start == interval._idx2time(4)
    assert not_visible_intervals[1].stop == interval._idx2time(5)


def test_visible_from_earth_location_event_with_velocity_transform():
    interval = _make_interval()
    gps_paris = GPS(
        longitude=scalar(2.333333).with_units(U.degree),
        latitude=scalar(48.866667).with_units(U.degree),
        altitude=scalar(0).with_units(U.meter),
        epoch=interval.start,
    )
    local_frame = gps_paris.earth_local_frame

    # Build a local-frame trajectory with known visibility, then map it to ITRF.
    z = np.asarray([-1.0, -0.2, 3.0, 4.0, -2.0, -1.0], dtype=np.float64)
    pos_local = Tensor(
        np.asarray(
            [
                np.zeros_like(z),
                np.zeros_like(z),
                z,
            ],
            dtype=np.float64,
        ),
        U.meter,
    )
    pos_itrf = local_frame.transform.undo(pos_local)

    vel_itrf = Tensor(
        np.zeros((3, interval.steps), dtype=np.float64),
        U.meter / U.second,
    )

    trajectory = Trajectory(interval, AbsoluteFrame.ITRF, pos_itrf, vel_itrf)
    event = VisibleFromEarthLocationEvent(trajectory, gps_paris)

    visible = event._time_map.get_values("visible")
    assert np.array_equal(visible, np.asarray([False, False, True, True, False, False], dtype=bool))
