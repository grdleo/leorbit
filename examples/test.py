from typing import NamedTuple

from leorbit import TimeInterval
from leorbit.api import scalar, get_passes, GPS
from leorbit.ext import CelestrakDataGP
from leorbit.propagator import SGP4
from leorbit.sky_object import Satellite
from leorbit.time import Timeline, Timestamp

import pytest

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