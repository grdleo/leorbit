from typing import NamedTuple

import pytest
import numpy as np

from leorbit.api import get_passes
from leorbit.coordinates import GPS, Coordinates
from leorbit.ext import CelestrakDataGP
from leorbit.frames import EarthLocalFrame
from leorbit.mathematics import U, normalize_angle_symmetric, scalar
from leorbit.propagator import SGP4
from leorbit.sky_object import Satellite
from leorbit.time import Timeline, Timestamp, TimeInterval


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

def test_pass_melbourne():
    """
    2026-05-14

    Integration test by computing a ISS pass over a location in Melbourne, Australia.
    Coordinates compared using [`in-the-sky.org`](https://in-the-sky.org/satpasseschart.php?utc1=1778787426&utc3=1778787788&satid=25544)
    Accuracy well within margin errors (usually always less than 1 degree when comparing altitude and azimuth).
    Therefore the values are considered true and hardcoded as reference values
    """

    oe14may26 = CelestrakDataGP(
        **{
            "OBJECT_NAME": "ISS (ZARYA)",
            "OBJECT_ID": "1998-067A",
            "EPOCH": "2026-05-14T04:45:57.957408",
            "MEAN_MOTION": 15.49211692,
            "ECCENTRICITY": 0.00075358,
            "INCLINATION": 51.6312,
            "RA_OF_ASC_NODE": 108.3512,
            "ARG_OF_PERICENTER": 56.9254,
            "MEAN_ANOMALY": 303.2457,
            "EPHEMERIS_TYPE": 0,
            "CLASSIFICATION_TYPE": "U",
            "NORAD_CAT_ID": 25544,
            "ELEMENT_SET_NO": 999,
            "REV_AT_EPOCH": 56648,
            "BSTAR": 0.00010032304,
            "MEAN_MOTION_DOT": 5.122e-5,
            "MEAN_MOTION_DDOT": 0
        }
    )

    iss = Satellite(
        "ISS", 
        oe14may26.to_orbital_elements(),
        SGP4
    )

    gps_melbourne_port = GPS(
        longitude=scalar(144.93084).with_units("deg"),
        latitude=scalar(-37.84614).with_units("deg"),
        altitude=scalar(0).with_units("meter"),
    )
    melbourne_frame = gps_melbourne_port.earth_local_frame

    now = Timestamp.fromisoformat("2026-05-14T20:46:00+02:00")
    timeline = Timeline(
        start=now,
        stop=now + scalar(24).with_units("hour"),
        dt=scalar(5).with_units("second")
    )

    next_pass, *others = get_passes(iss, timeline, gps_melbourne_port)
    next_pass: TimeInterval

    assert next_pass.start.unixepoch == pytest.approx(1778787285., abs=1)
    assert next_pass.stop.unixepoch == pytest.approx(1778787920., abs=1)

    class EpochAzimutAltitude(NamedTuple):
        epoch: Timestamp
        azimuth_deg: float
        altitude_deg: float

    test_epochs = [
        EpochAzimutAltitude(
            Timestamp(unixepoch=1778787480.0),
            -6.802126292549803,
            17.090051058943523
        ),
        EpochAzimutAltitude(
            Timestamp(unixepoch=1778787540.0), 
            12.122262102801265, 
            25.608692561518314
        ),
        EpochAzimutAltitude(
            Timestamp(unixepoch=1778787600.0), 
            45.120883892470964, 
            31.159225799278193
        ),
        EpochAzimutAltitude(
            Timestamp(unixepoch=1778787660.0), 
            79.02958595520127, 
            26.30998255767017
        ),
        EpochAzimutAltitude(
            Timestamp(unixepoch=1778787720.0), 
            98.9348867597868, 
            17.758579003225403
        )
    ]

    for (epoch, az, alt) in test_epochs:
        coords = iss.coordinates(epoch).horizontal(melbourne_frame)
        assert coords.azimuth.scalar.value("deg") == pytest.approx(az, rel=1e-6)
        assert coords.altitude.scalar.value("deg") == pytest.approx(alt, rel=1e-6)


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
    timeline = TimeInterval(t0, t1, scalar(1).with_units(U.second))

    gre_coords = Coordinates.from_gps(
        longitude=scalar(np.deg2rad(5.71667)).with_units(U.radian),
        latitude=scalar(np.deg2rad(45.166672)).with_units(U.radian),
        altitude=scalar(0).with_units(U.meter),
        epoch=t0,
    )
    local_frame = EarthLocalFrame(gre_coords)

    min_altitude = scalar(np.deg2rad(10)).with_units(U.radian)
    windows = _find_visibility_windows(iss, local_frame, timeline, min_altitude)
    assert windows

    event_start, event_stop = windows[0]

    assert event_start.unixepoch == pytest.approx(Timestamp.fromisoformat("2024-02-09T08:17:15").unixepoch, abs=2)
    assert event_stop.unixepoch == pytest.approx(Timestamp.fromisoformat("2024-02-09T08:22:57").unixepoch, abs=2)

    c0 = iss.coordinates(event_start)
    hor0 = c0.horizontal(local_frame)
    azi0 = normalize_angle_symmetric(hor0.azimuth).scalar.value("deg")
    alt0 = hor0.altitude.scalar.value("deg")
    assert azi0 == pytest.approx(-164, abs=2)
    assert alt0 == pytest.approx(10, abs=1)

    c1 = iss.coordinates(event_stop)
    hor1 = c1.horizontal(local_frame)
    azi1 = normalize_angle_symmetric(hor1.azimuth).scalar.value("deg")
    alt1 = hor1.altitude.scalar.value("deg")
    assert azi1 == pytest.approx(81, abs=2)
    assert alt1 == pytest.approx(10, abs=1)

    tmid = event_start + scalar(206).with_units(U.second)
    cmid = iss.coordinates(tmid)
    hormid = cmid.horizontal(local_frame)
    azimid = normalize_angle_symmetric(hormid.azimuth).scalar.value("deg")
    altmid = hormid.altitude.scalar.value("deg")
    assert azimid == pytest.approx(120, abs=2)
    assert altmid == pytest.approx(24, abs=1)
