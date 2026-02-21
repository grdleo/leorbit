import pytest

from leorbit2.m import Quantity
from leorbit2.time import Time, TimeInterval, get_intersections_timelines


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
        (1652780572, 123 * Quantity.second),
        (1652780572, -2.5 * 7 * Quantity.day),
        (-10000000000, 623 * Quantity.year),
    ),
)
def test_shift(unix: float, shift):
    assert (Time(unix) + shift).unixepoch == pytest.approx(unix + shift.magnitude("second"))


@pytest.mark.parametrize(
    "inp, outp",
    (
        ("2022-01-01T00:00:00", "22000.00000000"),
        ("2022-01-07T00:00:00", "22006.00000000"),
        ("2022-01-01T12:00:00", "22000.50000000"),
    ),
)
def test_yearday(inp: str, outp: str):
    t = Time.fromisoformat(inp)
    assert t.year_day == outp


@pytest.mark.parametrize(
    "iso, stl0_deg",
    (
        ("2017-01-01T00:00:00", 100.83793932378164),
    ),
)
def test_stl0(iso: str, stl0_deg: float):
    t = Time.fromisoformat(iso)
    assert t.stl0.magnitude("deg") == pytest.approx(stl0_deg)


@pytest.mark.parametrize(
    "start, stop, dt",
    [
        (
            Time.fromisoformat("2024-02-11T18:00:00"),
            Time.fromisoformat("2024-02-11T18:05:00"),
            1 * Quantity.minute,
        ),
        (
            Time.fromisoformat("2023-11-29T03:23:59"),
            Time.fromisoformat("2023-11-29T04:05:05"),
            42.3 * Quantity.second,
        ),
        (
            Time.fromisoformat("2012-01-13T14:44:09"),
            Time.fromisoformat("2018-12-29T01:02:03"),
            3.44 * 7 * Quantity.day,
        ),
    ],
)
def test_iter(start: Time, stop: Time, dt):
    timeline = TimeInterval(start, stop, dt)
    for i, t in enumerate(timeline):
        assert t.unixepoch == pytest.approx((start + dt * i).unixepoch)


@pytest.mark.parametrize(
    "tl, iters",
    [
        (
            TimeInterval(
                Time.fromisoformat("2024-02-11T18:00:12.34"),
                Time.fromisoformat("2024-02-11T18:00:16.44"),
                1 * Quantity.second,
            ),
            [
                Time.fromisoformat("2024-02-11T18:00:12.34"),
                Time.fromisoformat("2024-02-11T18:00:13.34"),
                Time.fromisoformat("2024-02-11T18:00:14.34"),
                Time.fromisoformat("2024-02-11T18:00:15.34"),
                Time.fromisoformat("2024-02-11T18:00:16.34"),
            ],
        ),
        (
            TimeInterval(
                Time.fromisoformat("2024-02-11T19:45:25.92"),
                Time.fromisoformat("2024-02-11T19:45:29.11"),
                1 * Quantity.second,
            ),
            [
                Time.fromisoformat("2024-02-11T19:45:25.92"),
                Time.fromisoformat("2024-02-11T19:45:26.92"),
                Time.fromisoformat("2024-02-11T19:45:27.92"),
                Time.fromisoformat("2024-02-11T19:45:28.92"),
            ],
        ),
    ],
)
def test_iter2(tl: TimeInterval, iters: list[Time]):
    assert tl.steps == len(iters)
    for t, tt in zip(tl, iters):
        assert t.unixepoch == pytest.approx(tt.unixepoch)


@pytest.mark.parametrize(
    "tl1, tl2, intersect",
    [
        (
            TimeInterval(
                Time.fromisoformat("2024-02-11T18:01:23"),
                Time.fromisoformat("2024-02-11T18:14:44"),
                1 * Quantity.second,
            ),
            TimeInterval(
                Time.fromisoformat("2024-02-11T18:08:08"),
                Time.fromisoformat("2024-02-11T18:35:22"),
                2.5 * Quantity.second,
            ),
            TimeInterval(
                Time.fromisoformat("2024-02-11T18:08:08"),
                Time.fromisoformat("2024-02-11T18:14:44"),
                1 * Quantity.second,
            ),
        ),
        (
            TimeInterval(
                Time.fromisoformat("2024-02-09T14:00:00"),
                Time.fromisoformat("2024-02-10T15:00:00"),
                30 * Quantity.second,
            ),
            TimeInterval(
                Time.fromisoformat("2024-02-10T18:00:30"),
                Time.fromisoformat("2024-02-12T09:00:00"),
                4 * Quantity.second,
            ),
            None,
        ),
    ],
)
def test_intersection(tl1: TimeInterval, tl2: TimeInterval, intersect: TimeInterval | None):
    dt = 1 * Quantity.second
    assert intersect == tl1.intersection(tl2, dt) == tl2.intersection(tl1, dt)


A_DATE = Time.fromisoformat("2024-02-11T18:01:23")
S = 1 * Quantity.second


@pytest.mark.parametrize(
    "tls1, tls2, inters",
    [
        (
            [
                TimeInterval(A_DATE + 3 * S, A_DATE + 16 * S),
                TimeInterval(A_DATE + 18 * S, A_DATE + 22 * S),
                TimeInterval(A_DATE + 26 * S, A_DATE + 30 * S),
            ],
            [
                TimeInterval(A_DATE + 7 * S, A_DATE + 19 * S),
                TimeInterval(A_DATE + 23 * S, A_DATE + 28 * S),
            ],
            [
                TimeInterval(A_DATE + 7 * S, A_DATE + 16 * S),
                TimeInterval(A_DATE + 18 * S, A_DATE + 19 * S),
                TimeInterval(A_DATE + 26 * S, A_DATE + 28 * S),
            ],
        )
    ],
)
def test_intersections(tls1: list[TimeInterval], tls2: list[TimeInterval], inters: list[TimeInterval]):
    inters_computed = get_intersections_timelines(tls1, tls2)
    assert set(inters_computed) == set(inters)
