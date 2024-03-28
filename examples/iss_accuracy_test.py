"""
    This example is more a proof of concept: ensure LEOrbit's propagation's accuracy.
    Using an online service, we get the true GPS coordinates of the ISS, and we compare it with LEOrbit's propagation.
    Long story short: it *is* accurate.
"""

import json

from leorbit import get_sat, Time, OrbitalElements

import requests
from time import sleep

# Create a `Satellite` ready to be propagated
iss = get_sat(25544, log=True)

for _ in range(60):
    # Get ISS GPS coordinates from an online API 
    # (therefore considered as true position, to be compared to)
    iss_loc_raw = requests.get("http://api.open-notify.org/iss-now.json").text
    iss_loc = json.loads(iss_loc_raw)

    # Ensure API returned correct data
    assert iss_loc["message"] == "success"

    # Get timestamp and coordinates from request (coordinates given in decimal degrees °)
    ts, lon, lat = ( 
        float(iss_loc["timestamp"]),
        float(iss_loc["iss_position"]["longitude"]),
        float(iss_loc["iss_position"]["latitude"]),
    )

    # Get ISS coordinates using LEOrbit
    epoch = Time(ts)
    gps = iss.coordinates(epoch).to_gps()

    # ... And print difference between true position and our computed position.
    # The more those values are close to zero, the better LEOrbit propagation is!
    print(
        f"∆Longitude = {(lon - gps.longitude_deg):.3f}°",
        f"∆Latitude = {(lat - gps.latitude_deg):.3f}°",
        "-------------------"
    )

    sleep(1)
