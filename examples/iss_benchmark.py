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
import os
import ssl
from time import perf_counter, sleep

import certifi
from pydantic import BaseModel
import requests
import urllib3

from leorbit.coordinates import GPS
from leorbit.mathematics import Quantity
from leorbit.time import Timestamp, TimeInterval
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
    def epoch(self) -> Timestamp:
        """UTC epoch associated with this sample."""
        return Timestamp(self.timestamp)

    @property
    def gps(self) -> GPS:
        """Sample converted to LEOrbit ``GPS`` representation.

        Open Notify does not provide altitude; we use a rough LEO value
        to build a full GPS point for display/comparison purposes.
        """
        return GPS(
            latitude=float(self.iss_position.latitude) * Quantity.degree,
            longitude=float(self.iss_position.longitude) * Quantity.degree,
            altitude=400 * Quantity.kilo_meter
        )


def _delta_deg(a, b) -> float:
    """Return absolute angular delta in degrees for two angle scalars."""
    return abs((a - b).scalar.value("deg"))


def _configure_ssl_for_https() -> None:
    """Point Python/requests TLS validation to certifi CA bundle."""
    ca_bundle = certifi.where()
    os.environ.setdefault("SSL_CERT_FILE", ca_bundle)
    os.environ.setdefault("REQUESTS_CA_BUNDLE", ca_bundle)
    os.environ.setdefault("CURL_CA_BUNDLE", ca_bundle)
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=ca_bundle)


def _configure_insecure_ssl_fallback() -> None:
    """Disable TLS certificate validation as a last-resort fallback."""
    os.environ["PYTHONHTTPSVERIFY"] = "0"
    os.environ["CURL_CA_BUNDLE"] = ""
    os.environ["REQUESTS_CA_BUNDLE"] = ""
    ssl._create_default_https_context = ssl._create_unverified_context
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def _print_sample(sample_idx: int, epoch: Timestamp, open_notify_gps: GPS, leorbit_gps: GPS) -> None:
    """Pretty-print one comparison sample."""
    dlat = _delta_deg(leorbit_gps.latitude, open_notify_gps.latitude)
    dlon = _delta_deg(leorbit_gps.longitude, open_notify_gps.longitude)

    print(f"\n── Sample #{sample_idx:02d} @ {epoch.isoformat}")
    print(f"   Open Notify : {open_notify_gps.dms}")
    print(f"   LEOrbit     : {leorbit_gps.dms}")
    print(f"   Error       : Δlat={dlat:7.4f}°, Δlon={dlon:7.4f}°")

def main():
    """Run the benchmark loop and print live comparison samples."""

    _configure_ssl_for_https()

    print("\n══════════════════════════════════════════════════════════════")
    print("  LEOrbit ISS benchmark vs Open Notify")
    print("══════════════════════════════════════════════════════════════")

    try:
        iss = get_satellite(25544, log=True)
        print(f"   ISS Epoch   : {iss.propagator.elements.epoch.isoformat}")
    except requests.exceptions.SSLError:
        print("Warning: TLS certificate validation failed; retrying with insecure SSL fallback.")
        _configure_insecure_ssl_fallback()
        iss = get_satellite(25544, log=True)

    t0 = Timestamp.now() - timedelta(hours=2)
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
