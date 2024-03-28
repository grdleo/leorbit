import pytest
from leorbit.time import Time
from leorbit.time import Timeline, get_intersections_timelines
from leorbit.math import Q_

@pytest.mark.parametrize(
    "start, stop, dt",
    [
        (Time.fromisoformat("2024-02-11T18:00:00"), Time.fromisoformat("2024-02-11T18:05:00"), Q_("1min")),
        (Time.fromisoformat("2023-11-29T03:23:59"), Time.fromisoformat("2023-11-29T04:05:05"), Q_("42.3s")),
        (Time.fromisoformat("2012-01-13T14:44:09"), Time.fromisoformat("2018-12-29T01:02:03"), Q_("3.44 week")),
    ] 
)
def test_iter(start: Time, stop: Time, dt: Q_):
    timeline = Timeline(start, stop, dt)
    for i, t in enumerate(timeline):
        assert t.unixepoch == pytest.approx((start + dt * i).unixepoch)

@pytest.mark.parametrize(
    "tl, iters",
    [
        (
            Timeline(Time.fromisoformat("2024-02-11T18:00:12.34"), Time.fromisoformat("2024-02-11T18:00:16.44"), Q_("1s")),
            [
                Time.fromisoformat("2024-02-11T18:00:12.34"),
                Time.fromisoformat("2024-02-11T18:00:13.34"),
                Time.fromisoformat("2024-02-11T18:00:14.34"),
                Time.fromisoformat("2024-02-11T18:00:15.34"),
                Time.fromisoformat("2024-02-11T18:00:16.34")
            ]
        ),
        (
            Timeline(Time.fromisoformat("2024-02-11T19:45:25.92"), Time.fromisoformat("2024-02-11T19:45:29.11"), Q_("1s")),
            [
                Time.fromisoformat("2024-02-11T19:45:25.92"),
                Time.fromisoformat("2024-02-11T19:45:26.92"),
                Time.fromisoformat("2024-02-11T19:45:27.92"),
                Time.fromisoformat("2024-02-11T19:45:28.92")
            ]
        ),
    ]
)
def test_iter2(tl: Timeline, iters: list[Time]):
    assert tl.steps == len(iters)
    for t, tt in zip(tl, iters):
        assert t == tt

@pytest.mark.parametrize(
    "tl1, tl2, intersect",
    [
        (
            Timeline(Time.fromisoformat("2024-02-11T18:01:23"), Time.fromisoformat("2024-02-11T18:14:44"), Q_("1s")),
            Timeline(Time.fromisoformat("2024-02-11T18:08:08"), Time.fromisoformat("2024-02-11T18:35:22"), Q_("2.5s")),
            Timeline(Time.fromisoformat("2024-02-11T18:08:08"), Time.fromisoformat("2024-02-11T18:14:44"), Q_("1s"))
        ),
        (
            Timeline(Time.fromisoformat("2024-02-09T14:00:00"), Time.fromisoformat("2024-02-10T15:00:00"), Q_("30s")),
            Timeline(Time.fromisoformat("2024-02-10T18:00:30"), Time.fromisoformat("2024-02-12T09:00:00"), Q_("4s")),
            None
        ),
    ] 
)
def test_intersection(tl1: Timeline, tl2: Timeline, intersect: Timeline | None):
    dt = Q_("1s")
    assert intersect == tl1.intersection(tl2, dt) == tl2.intersection(tl1, dt)

A_DATE = Time.fromisoformat("2024-02-11T18:01:23")
S = Q_("1s")
@pytest.mark.parametrize(
    "tls1, tls2, inters",
    [
        (
            [
                Timeline(A_DATE + 3*S, A_DATE + 16*S),
                Timeline(A_DATE + 18*S, A_DATE + 22*S),
                Timeline(A_DATE + 26*S, A_DATE + 30*S),
            ],
            [
                Timeline(A_DATE + 7*S, A_DATE + 19*S),
                Timeline(A_DATE + 23*S, A_DATE + 28*S)
            ],
            [
                Timeline(A_DATE + 7*S, A_DATE + 16*S),
                Timeline(A_DATE + 18*S, A_DATE + 19*S),
                Timeline(A_DATE + 26*S, A_DATE + 28*S),
            ]
        )
    ]
)
def test_intersections(tls1: list[Timeline], tls2: list[Timeline], inters: list[Timeline]):
    inters_computed = get_intersections_timelines(tls1, tls2)
    assert set(inters_computed) == set(inters)