import pytest

from leorbit.coordinates import Coordinates
from leorbit.ext import CelestrakDataGP
from leorbit.frames import EarthLocalFrame
from leorbit.m import Quantity, normalize_angle_symmetric
from leorbit.propagator import SGP4
from leorbit.sky_object import Satellite
from leorbit.time import Timestamp, TimeInterval


def _find_visibility_windows(sat: Satellite, local_frame: EarthLocalFrame, timeline: TimeInterval, min_altitude):
    windows: list[tuple[Timestamp, Timestamp]] = []

    start = None
    stop = None

    for t in timeline:
        hor = sat.coordinates(t).horizontal(local_frame)
        visible = hor.altitude >= min_altitude

        if visible:
            if start is None:
                start = t
            stop = t
        elif start is not None and stop is not None:
            windows.append((start, stop))
            start = None
            stop = None

    if start is not None and stop is not None:
        windows.append((start, stop))

    return windows


def test_event_visibility_integration():
    oe_09fev24 = {
        "OBJECT_NAME": "ISS (ZARYA)",
        "OBJECT_ID": "1998-067A",
        "EPOCH": "2024-02-08T00:14:27.133728",
        "MEAN_MOTION": 15.49623054,
        "ECCENTRICITY": 0.0002067,
        "INCLINATION": 51.6398,
        "RA_OF_ASC_NODE": 240.4555,
        "ARG_OF_PERICENTER": 215.9229,
        "MEAN_ANOMALY": 294.1058,
        "EPHEMERIS_TYPE": 0,
        "CLASSIFICATION_TYPE": "U",
        "NORAD_CAT_ID": 25544,
        "ELEMENT_SET_NO": 999,
        "REV_AT_EPOCH": 43836,
        "BSTAR": 0.00031432,
        "MEAN_MOTION_DOT": 0.00017286,
        "MEAN_MOTION_DDOT": 0,
    }

    oe_iss = CelestrakDataGP(**oe_09fev24).to_orbital_elements()
    iss = Satellite("ISS", oe_iss, SGP4)

    t0 = Timestamp.fromisoformat("2024-02-09T08:15:00")
    t1 = Timestamp.fromisoformat("2024-02-09T08:30:00")
    timeline = TimeInterval(t0, t1, 1 * Quantity.second)

    gre_coords = Coordinates.from_gps(
        longitude=5.71667 * Quantity.deg,
        latitude=45.166672 * Quantity.deg,
        altitude=0 * Quantity.meter,
        epoch=t0,
    )
    local_frame = EarthLocalFrame(gre_coords)

    min_altitude = 10 * Quantity.deg
    windows = _find_visibility_windows(iss, local_frame, timeline, min_altitude)
    assert windows

    event_start, event_stop = windows[0]

    assert event_start.unixepoch == pytest.approx(Timestamp.fromisoformat("2024-02-09T08:17:15").unixepoch, abs=2)
    assert event_stop.unixepoch == pytest.approx(Timestamp.fromisoformat("2024-02-09T08:22:57").unixepoch, abs=2)

    c0 = iss.coordinates(event_start)
    hor0 = c0.horizontal(local_frame)
    azi0 = normalize_angle_symmetric(hor0.azimuth).magnitude("deg")
    alt0 = hor0.altitude.magnitude("deg")
    assert azi0 == pytest.approx(-164, abs=2)
    assert alt0 == pytest.approx(10, abs=1)

    c1 = iss.coordinates(event_stop)
    hor1 = c1.horizontal(local_frame)
    azi1 = normalize_angle_symmetric(hor1.azimuth).magnitude("deg")
    alt1 = hor1.altitude.magnitude("deg")
    assert azi1 == pytest.approx(81, abs=2)
    assert alt1 == pytest.approx(10, abs=1)

    tmid = event_start + 206 * Quantity.second
    cmid = iss.coordinates(tmid)
    hormid = cmid.horizontal(local_frame)
    azimid = normalize_angle_symmetric(hormid.azimuth).magnitude("deg")
    altmid = hormid.altitude.magnitude("deg")
    assert azimid == pytest.approx(120, abs=2)
    assert altmid == pytest.approx(24, abs=1)
