"""Tests the accuracy of the `Orbit` propagation using ISS as a reference,
since it is easy to find its coordinates online"""

import json

from leorbit.math import Q_
from leorbit.orbit.objects import Satellite
from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.orbit.propagator import SGP4
from leorbit.simulation.utils import compute_magnitude
from leorbit.time import Time
from leorbit.time import Timeline
from leorbit.simulation.event import VisibleFromLocation, NightTime, AstroEvent
from leorbit.math.coordinate import GPSCoordinates

from time import perf_counter_ns

GP_ISS = {
    "OBJECT_NAME": "ISS (ZARYA)",
    "OBJECT_ID": "1998-067A",
    "EPOCH": "2024-03-17T10:36:39.376224",
    "MEAN_MOTION": 15.49037779,
    "ECCENTRICITY": 0.000446,
    "INCLINATION": 51.6409,
    "RA_OF_ASC_NODE": 50.1325,
    "ARG_OF_PERICENTER": 336.4542,
    "MEAN_ANOMALY": 166.8282,
    "EPHEMERIS_TYPE": 0,
    "CLASSIFICATION_TYPE": "U",
    "NORAD_CAT_ID": 25544,
    "ELEMENT_SET_NO": 999,
    "REV_AT_EPOCH": 44431,
    "BSTAR": 0.00023745,
    "MEAN_MOTION_DOT": 0.00012672,
    "MEAN_MOTION_DDOT": 0
}

if __name__ == "__main__":
    # Create `Satellite` object from already existant GP data
    oe_iss: OrbitalElements = OrbitalElements.from_celestrak_json(GP_ISS)
    iss = Satellite(SGP4(oe_iss), name="ISS", catnr=25544)

    # Create GPS coordinates of "Grenoble (France)" location
    gre_gps = GPSCoordinates(
        longitude_deg=5.71667,
        latitude_deg=45.166672,
        location_name="Grenoble (FR)"
    )

    # Build an event corresponding to when ISS is visible from Grenoble AND also when it is night time 
    iss_visible_in_gre = VisibleFromLocation(iss, gre_gps, min_altitude_visi=Q_("10°"))
    night_time_in_gre = NightTime(gre_gps)
    night_visibility_in_gre: AstroEvent = iss_visible_in_gre & night_time_in_gre

    # Create a timeline covering the next 4 days (dt=5s)
    epoch = oe_iss.epoch
    timeline = Timeline(oe_iss.epoch, oe_iss.epoch + Q_("4 day"), Q_("5s"))

    # Compute the time intervals where our event occurs
    t0 = perf_counter_ns()
    time_intervals = night_visibility_in_gre.compute_intervals(timeline, minimal_duration=Q_("30s"))
    t1 = perf_counter_ns()
    print(f"Took {(t1 - t0) * 1e-9:.1f}s to compute")

    # Print the computed intervals
    for tl in time_intervals:
        horis = [iss.coordinates(t).to_horizontal_coordinates(gre_gps) for t in tl]
        mags = [compute_magnitude(t, gre_gps.earth_local_frame(), iss) for t in tl] # Compute the magnitude values of ISS over the timeline
        max_alti = max(h.altitude_deg for h in horis)
        print(f"Visibility in {gre_gps.location_name} from {tl.start.human} to {tl.stop.human}. "
              f"Max altitude={max_alti:.0f}°. Max brightness={min(mags):.2f}. ")
    
    # Get the first interval that matches our event,
    # and make it more precise with a new dt
    custom_event_tl = time_intervals[0].duplicate(dt=Q_("1s"))

    # Recompute coordinates over new more precise timeline
    coordinates_custom_event = iss.simulation(custom_event_tl)

    # Create a text export file of horizontal coordinates of selected observation interval
    text_export = (
        "ISS Observation in Grenoble (FR) - Angles in degrees [°]\n"
        "--------------------------------------------------------\n"
        "              Time in ISO format	Azimuth	Altitude\n"
    )
    for c in coordinates_custom_event:
        hori = c.to_horizontal_coordinates(gre_gps)
        text_export += f"{c.epoch.isoformat}	{hori.azimuth_deg:03.2f}	{hori.altitude_deg:03.2f}\n"
    
    print(text_export)
    