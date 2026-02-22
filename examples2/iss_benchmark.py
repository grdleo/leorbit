"""
    Quick benchmark for LEOrbit propagation performance.  Similar to
    ``iss_accuracy_test.py`` but instead of querying a remote API we just
    propagate the ISS state repeatedly and measure how long it takes.

    Run with ``python examples2/iss_benchmark.py`` and you should see
    per-propagation timings printed.
"""

from datetime import timedelta
from time import perf_counter, sleep
from typing import Any

import numpy as np
from pydantic import BaseModel
import requests

from leorbit.coordinates import GPS
from leorbit.frames import AbsoluteFrame
from leorbit.m import Quantity
from leorbit.time import Time, TimeInterval
from leorbit import get_satellite
from leorbit.m import Vector3

class OpenNotifyPosition(BaseModel):
    latitude: str # [degrees]
    longitude: str # [degrees]

class OpenNotifyIssResponse(BaseModel):
    """http://api.open-notify.org/iss-now.json"""
    message: str
    timestamp: int
    iss_position: OpenNotifyPosition

    @staticmethod
    def retrieve() -> OpenNotifyIssResponse:
        response = requests.get("http://api.open-notify.org/iss-now.json")
        return OpenNotifyIssResponse(**response.json())
    
    @property
    def epoch(self) -> Time:
        return Time(self.timestamp)

    @property
    def gps(self) -> GPS:
        return GPS(
            latitude=float(self.iss_position.latitude) * Quantity.deg,
            longitude=float(self.iss_position.longitude) * Quantity.deg,
            altitude=400 * Quantity.kilo_meter
        )
    
def _repr_vec_pos(v: Vector3[Any]) -> str:
    vals = v._values.flatten()
    return f"Position: ({vals[0]:.0f} m, {vals[1]:.0f} m, {vals[2]:.0f} m)"

def main():
    # build a satellite object once
    iss = get_satellite(25544, log=True)

    # pick a fixed epoch to propagate from and include margin for API timestamp drift
    t0 = Time.now() - timedelta(hours=2)
    timeline = TimeInterval(
        t0,
        t0 + timedelta(hours=4),
    )

    _time_flag = perf_counter()
    trajectory = iss.trajectory(timeline)
    print(f"Computed trajectory for 4 hours with 1s step in {perf_counter() - _time_flag:.4f}s")

    for _ in range(16):
        open_notify_coords = OpenNotifyIssResponse.retrieve()
        open_notify_gps = open_notify_coords.gps
        epoch = open_notify_coords.epoch
        gps = trajectory.gps_at(epoch)

        on_pos_itrf = open_notify_coords.gps.to_coordinates(epoch).get_pos(AbsoluteFrame.ITRF)
        leorbit_pos_itrf = trajectory.get_pos(epoch, AbsoluteFrame.ITRF)

        # display results
        print(f"At epoch {epoch.isoformat}: ")
        print(f"Open notify position: {open_notify_gps.dms}")
        print(f"Leorbit position:     {gps.dms}")

        # compute simple latitude/longitude differences in degrees
        dlat = float(np.asarray((gps.latitude - open_notify_gps.latitude).magnitude("deg")).reshape(-1)[0])
        dlon = float(np.asarray((gps.longitude - open_notify_gps.longitude).magnitude("deg")).reshape(-1)[0])
        dlat = abs(dlat)
        dlon = abs(dlon)
        print(f"Difference: Δlat={dlat:.4f}°, Δlon={dlon:.4f}°")

        sleep(5)

if __name__ == "__main__":
    main()
