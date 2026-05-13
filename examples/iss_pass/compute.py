from __future__ import annotations

import csv
from datetime import datetime, timedelta
from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
from matplotlib.projections.polar import PolarAxes
import numpy as np

from leorbit import get_satellite
from leorbit.api import scalar
from leorbit.coordinates import Coordinates
from leorbit.frames import EarthLocalFrame
from leorbit.time import Timestamp, TimeInterval


NORAD_ID = 25544 # ISS
EPOCH_DT = datetime.fromisoformat("2026-03-01T06:00:00+01:00")

# Grenoble (FR)
OBS_LAT_DEG = 45.17
OBS_LON_DEG = 5.72


def _float_array_deg(tensor) -> np.ndarray:
    return np.asarray(tensor.raw_data_array("deg"), dtype=float).reshape(-1)


def _equatorial_to_horizontal(
    eq_vectors: np.ndarray,
    latitude_deg: float,
    longitude_deg: float,
    epoch: Timestamp,
) -> tuple[np.ndarray, np.ndarray]:
    lat = np.deg2rad(latitude_deg)
    lon = np.deg2rad(longitude_deg)
    stl0 = epoch.stl0.scalar.value("radian")
    local_sidereal = stl0 + lon

    x = eq_vectors[0, :]
    y = eq_vectors[1, :]
    z = eq_vectors[2, :]

    sin_t = np.sin(local_sidereal)
    cos_t = np.cos(local_sidereal)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    east = -sin_t * x + cos_t * y
    north = -sin_lat * cos_t * x - sin_lat * sin_t * y + cos_lat * z
    up = cos_lat * cos_t * x + cos_lat * sin_t * y + sin_lat * z

    up = np.clip(up, -1.0, 1.0)
    alt = np.arcsin(up)
    az = (np.arctan2(east, north) + 2 * np.pi) % (2 * np.pi)

    return np.rad2deg(az), np.rad2deg(alt)


def _plane_horizontal_curve(
    plane: str,
    latitude_deg: float,
    longitude_deg: float,
    epoch: Timestamp,
    n_samples: int = 721,
) -> tuple[np.ndarray, np.ndarray]:
    angle = np.linspace(0.0, 2 * np.pi, n_samples)

    if plane == "ecliptic":
        obliquity = np.deg2rad(23.439291)
        eq_vec = np.vstack(
            [
                np.cos(angle),
                np.sin(angle) * np.cos(obliquity),
                np.sin(angle) * np.sin(obliquity),
            ]
        )
    elif plane == "galactic":
        r_eq_to_gal = np.array(
            [
                [-0.0548755604, -0.8734370902, -0.4838350155],
                [0.4941094279, -0.4448296300, 0.7469822445],
                [-0.8676661490, -0.1980763734, 0.4559837762],
            ]
        )
        r_gal_to_eq = r_eq_to_gal.T
        gal_vec = np.vstack(
            [
                np.cos(angle),
                np.sin(angle),
                np.zeros_like(angle),
            ]
        )
        eq_vec = r_gal_to_eq @ gal_vec
    else:
        raise ValueError(f"Unknown plane '{plane}'")

    return _equatorial_to_horizontal(eq_vec, latitude_deg, longitude_deg, epoch)


def _plot_wrapped_polar(
    ax: PolarAxes,
    az_deg: np.ndarray,
    el_deg: np.ndarray,
    *,
    label: str,
    color: str = "black",
    linestyle: str = "-",
    linewidth: float = 1.0,
) -> None:
    visible = el_deg >= 0.0
    if not np.any(visible):
        return

    labeled = False
    indices = np.where(visible)[0]
    splits = np.where(np.diff(indices) > 1)[0]
    start_k = 0

    for split_k in list(splits) + [len(indices) - 1]:
        seg_indices = indices[start_k:split_k + 1]
        seg_az = az_deg[seg_indices]
        seg_el = el_deg[seg_indices]

        if seg_indices.size < 2:
            start_k = split_k + 1
            continue

        theta = np.deg2rad(seg_az % 360.0)
        r = seg_el

        wrap_breaks = np.where(np.abs(np.diff(seg_az)) > 180.0)[0]
        wstart = 0
        for idx in list(wrap_breaks) + [len(seg_az) - 1]:
            ax.plot(
                theta[wstart:idx + 1],
                r[wstart:idx + 1],
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label if not labeled else None,
            )
            labeled = True
            wstart = idx + 1

        start_k = split_k + 1


def plot_horizontal_trajectory(
    az_deg: np.ndarray,
    el_deg: np.ndarray,
    timestamp_indices: np.ndarray,
    timestamp_labels: list[str],
    norad_id: int,
    start_dt: datetime,
    stop_dt: datetime,
) -> None:
    r = el_deg
    theta = np.deg2rad(az_deg)

    fig = plt.figure(figsize=(8, 8))
    ax = cast(PolarAxes, fig.add_subplot(111, projection="polar"))

    # ecl_az_deg, ecl_el_deg = _plane_horizontal_curve(
    #     "ecliptic",
    #     latitude_deg=OBS_LAT_DEG,
    #     longitude_deg=OBS_LON_DEG,
    #     epoch=Timestamp.fromisoformat(EPOCH_DT.isoformat()),
    # )
    # gal_az_deg, gal_el_deg = _plane_horizontal_curve(
    #     "galactic",
    #     latitude_deg=OBS_LAT_DEG,
    #     longitude_deg=OBS_LON_DEG,
    #     epoch=Timestamp.fromisoformat(EPOCH_DT.isoformat()),
    # )

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.set_rlim(90, 0)

    ax.set_thetagrids(np.arange(0, 360, 15))
    ax.set_rticks(np.arange(0, 91, 5))
    ax.set_rlabel_position(225)

    ax.grid(True, color="0.85", linewidth=0.6)
    # _plot_wrapped_polar(ax, ecl_az_deg, ecl_el_deg, label="Ecliptic plane", color="black", linestyle="-", linewidth=1.1)
    # _plot_wrapped_polar(ax, gal_az_deg, gal_el_deg, label="Galactic plane", color="black", linestyle="--", linewidth=1.1)
    _plot_wrapped_polar(ax, az_deg, el_deg, label="ISS trajectory", color="tab:red", linestyle="-", linewidth=1.8)

    for i, label in zip(timestamp_indices, timestamp_labels):
        if el_deg[int(i)] < 0:
            continue
        theta_i = theta[int(i)]
        r_i = r[int(i)]
        ax.plot([theta_i], [r_i], marker="o", markersize=2.0, color="tab:red")
        ax.text(theta_i, min(90.0, r_i + 1.2), label, fontsize=6, color="black", ha="center", va="bottom")

    ax.legend(loc="lower left", framealpha=0.85)

    ax.set_title(
        f"ISS pass ({norad_id}) — {start_dt.isoformat()} to {stop_dt.isoformat()}\\n"
        "Polar sky plot (top = South)"
    )

    fig.tight_layout()
    plt.show()


def export_trajectory_csv(
    file_path: Path,
    epochs: list[datetime],
    altitude_deg: np.ndarray,
    azimuth_deg: np.ndarray,
) -> None:
    file_path.parent.mkdir(parents=True, exist_ok=True)

    with file_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["epoch", "altitude", "azimuth"])
        for epoch, alt, az in zip(epochs, altitude_deg, azimuth_deg):
            writer.writerow([epoch.isoformat(), f"{float(alt):.8f}", f"{float(az):.8f}"])


def main() -> None:
    start_dt = EPOCH_DT - timedelta(minutes=2)
    stop_dt = EPOCH_DT + timedelta(minutes=2)

    tmax = Timestamp.fromisoformat(EPOCH_DT.isoformat())
    timeline = TimeInterval(
        tmax - scalar("2 minute"),
        tmax + scalar("2 minute"),
        dt=scalar("1 second")
    )

    satellite = get_satellite(NORAD_ID, log=True)
    trajectory = satellite.trajectory(timeline)

    observer = Coordinates.from_gps(
        longitude=scalar(OBS_LON_DEG).with_units("degree"),
        latitude=scalar(OBS_LAT_DEG).with_units("degree"),
        altitude=scalar("0 meter"),
        epoch=timeline.start,
    )
    local_frame = EarthLocalFrame(observer)

    horizontal = trajectory.horizontal(local_frame)
    az_deg = _float_array_deg(horizontal.azimuth) % 360
    el_deg = _float_array_deg(horizontal.altitude)
    epochs = [start_dt + timedelta(seconds=i) for i in range(az_deg.size)]

    csv_path = Path(__file__).resolve().with_name("trajectory.csv")
    export_trajectory_csv(
        file_path=csv_path,
        epochs=epochs,
        altitude_deg=el_deg,
        azimuth_deg=az_deg,
    )
    print(f"CSV exported: {csv_path}")

    timestamp_indices = np.arange(0, az_deg.size, 10, dtype=int)
    timestamp_labels = [(start_dt + timedelta(seconds=int(i))).strftime("%H:%M") for i in timestamp_indices]

    plot_horizontal_trajectory(
        az_deg=az_deg,
        el_deg=el_deg,
        timestamp_indices=timestamp_indices,
        timestamp_labels=timestamp_labels,
        norad_id=NORAD_ID,
        start_dt=start_dt,
        stop_dt=stop_dt,
    )


if __name__ == "__main__":
    main()
