"""
    Quick benchmark for LEOrbit propagation performance.  Similar to
    ``iss_accuracy_test.py`` but instead of querying a remote API we just
    propagate the ISS state repeatedly and measure how long it takes.

    Run with ``python examples2/iss_benchmark.py`` and you should see
    per-propagation timings printed.
"""

from datetime import timedelta
from time import perf_counter

from leorbit2.m import Quantity
from leorbit2.time import Time, TimeInterval
from leorbit2 import get_satellite


def main():
    # build a satellite object once
    iss = get_satellite(25544, log=True)

    # pick a fixed epoch to propagate from
    t0 = Time.now()
    timeline = TimeInterval(
        t0,
        t0 + timedelta(hours=1),
    )

    _time_flag = perf_counter()
    trajectory = iss.trajectory(timeline)
    print(f"Computed trajectory for 1 hour with 1s step in {perf_counter() - _time_flag:.4f}s")
    print(trajectory.name)

if __name__ == "__main__":
    main()
