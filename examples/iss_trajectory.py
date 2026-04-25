import os
import ssl
from time import perf_counter
from typing import cast
from urllib.error import URLError

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import certifi
from cartopy.mpl.geoaxes import GeoAxes

from leorbit import get_satellite
from leorbit.mathematics import Quantity
from leorbit.time import Timestamp, TimeInterval


def _configure_ssl_for_cartopy() -> None:
    ca_bundle = certifi.where()
    os.environ.setdefault("SSL_CERT_FILE", ca_bundle)
    os.environ.setdefault("REQUESTS_CA_BUNDLE", ca_bundle)
    os.environ.setdefault("CURL_CA_BUNDLE", ca_bundle)
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=ca_bundle)


def main() -> None:
    _configure_ssl_for_cartopy()

    iss = get_satellite(25544, log=True)

    start = Timestamp.now()
    timeline = TimeInterval(
        start,
        start + 72 * Quantity.hour,
        dt=5 * Quantity.second,
    )

    t0 = perf_counter()
    trajectory = iss.trajectory(timeline)
    print(f"Trajectory precompute: 72h @ 5s step in {perf_counter() - t0:.4f}s")

    t0 = perf_counter()
    gps = trajectory.gps()
    longitudes: list[float] = gps.longitude.raw_data_array("deg").tolist()
    latitudes: list[float] = gps.latitude.raw_data_array("deg").tolist()
    print(f"GPS extraction: 72h @ 5s step in {perf_counter() - t0:.4f}s")

    fig = plt.figure(figsize=(12, 6))
    ax = cast(GeoAxes, plt.axes(projection=ccrs.Robinson(central_longitude=90)))
    ax.set_global()
    ax.set_facecolor("0.9")
    ax.gridlines(draw_labels=False, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
    try:
        ax.add_feature(cfeature.LAND, facecolor="0.85", edgecolor="none")
        ax.add_feature(cfeature.OCEAN, facecolor="0.75", edgecolor="none")
        ax.coastlines(linewidth=0.6)
    except URLError as exc:
        raise RuntimeError(
            "Could not download Cartopy Natural Earth map data. "
            "Check internet/certificate configuration, or pre-download Cartopy data."
        ) from exc

    ax.plot(
        longitudes,
        latitudes,
        color="red",
        linewidth=0.8,
        transform=ccrs.Geodetic(),
    )

    ax.set_title("ISS trajectory for the next 72 hours (sampled every 5s)")
    plt.show()


if __name__ == "__main__":
    main()
