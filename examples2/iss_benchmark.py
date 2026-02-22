"""LEOrbit ISS benchmark against Open Notify live position.

This script:
1) builds a 4-hour trajectory for the ISS using LEOrbit,
2) polls Open Notify for the live ISS GPS position,
3) compares LEOrbit GPS at the same epoch,
4) prints a readable summary with latitude/longitude errors.

Run:
    python examples2/iss_benchmark.py
"""

from datetime import timedelta
from time import perf_counter, sleep

import numpy as np
from pydantic import BaseModel
import requests

from leorbit.coordinates import GPS
from leorbit.m import Quantity
from leorbit.time import Time, TimeInterval
from leorbit import get_satellite

class OpenNotifyIssResponse(BaseModel):
    """Response payload for ``http://api.open-notify.org/iss-now.json``."""

    message: str
    timestamp: int
    iss_position: "OpenNotifyIssResponse.OpenNotifyPosition"

    class OpenNotifyPosition(BaseModel):
        """Raw latitude/longitude values returned by the API, in degrees."""

        latitude: str
        longitude: str

    @staticmethod
    def retrieve() -> OpenNotifyIssResponse:
        """Fetch and validate one Open Notify sample."""
        response = requests.get("http://api.open-notify.org/iss-now.json")
        response.raise_for_status()
        return OpenNotifyIssResponse(**response.json())
    
    @property
    def epoch(self) -> Time:
        """UTC epoch associated with this sample."""
        return Time(self.timestamp)

    @property
    def gps(self) -> GPS:
        """Sample converted to LEOrbit ``GPS`` representation.

        Open Notify does not provide altitude; we use a rough LEO value
        to build a full GPS point for display/comparison purposes.
        """
        return GPS(
            latitude=float(self.iss_position.latitude) * Quantity.deg,
            longitude=float(self.iss_position.longitude) * Quantity.deg,
            altitude=400 * Quantity.kilo_meter
        )


def _delta_deg(a, b) -> float:
    """Return absolute angular delta in degrees for two angle scalars."""
    return abs(float(np.asarray((a - b).magnitude("deg")).reshape(-1)[0]))


def _print_sample(sample_idx: int, epoch: Time, open_notify_gps: GPS, leorbit_gps: GPS) -> None:
    """Pretty-print one comparison sample."""
    dlat = _delta_deg(leorbit_gps.latitude, open_notify_gps.latitude)
    dlon = _delta_deg(leorbit_gps.longitude, open_notify_gps.longitude)

    print(f"\n── Sample #{sample_idx:02d} @ {epoch.isoformat}")
    print(f"   Open Notify : {open_notify_gps.dms}")
    print(f"   LEOrbit     : {leorbit_gps.dms}")
    print(f"   Error       : Δlat={dlat:7.4f}°, Δlon={dlon:7.4f}°")

def main():
    """Run the benchmark loop and print live comparison samples."""

    print("\n══════════════════════════════════════════════════════════════")
    print("  LEOrbit ISS benchmark vs Open Notify")
    print("══════════════════════════════════════════════════════════════")

    iss = get_satellite(25544, log=True)

    t0 = Time.now() - timedelta(hours=2)
    timeline = TimeInterval(
        t0,
        t0 + timedelta(hours=4),
    )

    _time_flag = perf_counter()
    trajectory = iss.trajectory(timeline)
    dt = perf_counter() - _time_flag
    print(f"\nTrajectory precompute: 4h @ 1s step in {dt:.4f}s")

    for i in range(1, 17):
        open_notify_coords = OpenNotifyIssResponse.retrieve()
        open_notify_gps = open_notify_coords.gps
        epoch = open_notify_coords.epoch
        leorbit_gps = trajectory.gps_at(epoch)

        _print_sample(i, epoch, open_notify_gps, leorbit_gps)

        sleep(5)

    print("\nDone.")

if __name__ == "__main__":
    main()
