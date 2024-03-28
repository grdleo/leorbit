from leorbit.time import Time
from leorbit.math import Q_

import pytest


@pytest.mark.parametrize(
    "unix, iso",
    (
        (1652780572, "2022-05-17T09:42:52"),
        (946684800, "2000-01-01T00:00:00"),
        (0, "1970-01-01T00:00:00"),
        (-10000000000, "1653-02-10T06:13:20"),
    ),
)
def test_instance(unix: float, iso: str):
    assert Time(unix) == Time.fromisoformat(iso)


@pytest.mark.parametrize(
    "unix, shift",
    (
        (1652780572, Q_(123, "s")),
        (1652780572, Q_(-2.5, "week")),
        (-10000000000, Q_(623, "year")),
    ),
)
def test_shift(unix: float, shift: Q_):
    assert (Time(unix) + shift).unixepoch == (unix + shift.m_as("s"))


@pytest.mark.parametrize(
    "inp, outp",
    (
        ("2022-01-01T00:00:00", "22000.00000000"),
        ("2022-01-07T00:00:00", "22006.00000000"),
        ("2022-01-01T12:00:00", "22000.50000000"),
    ),
)
def test_yearday(inp: str, outp: str):
    t: Time = Time.fromisoformat(inp)
    assert t.year_day == outp

@pytest.mark.parametrize(
    "iso, stl0",
    (
        ("2017-01-01T00:00:00", Q_("100.83793932378164°")),
    )
)
def test_stl0(iso: str, stl0: Q_):
    t: Time = Time.fromisoformat(iso)
    assert t.stl0 == stl0