from leorbit.math.coordinate import GPSCoordinates
from leorbit.math import Q_
from leorbit.orbit.objects import Satellite
from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.orbit.propagator import SGP4
from leorbit.time import Timeline, Time
from leorbit.simulation.event import VisibleFromLocation

from pytest import approx

def test_event_visibility():
    oe_09fev24 = {
        "OBJECT_NAME":"ISS (ZARYA)",
        "OBJECT_ID":"1998-067A",
        "EPOCH":"2024-02-08T00:14:27.133728",
        "MEAN_MOTION":15.49623054,
        "ECCENTRICITY":0.0002067,
        "INCLINATION":51.6398,
        "RA_OF_ASC_NODE":240.4555,
        "ARG_OF_PERICENTER":215.9229,
        "MEAN_ANOMALY":294.1058,
        "EPHEMERIS_TYPE":0,
        "CLASSIFICATION_TYPE":"U",
        "NORAD_CAT_ID":25544,
        "ELEMENT_SET_NO":999,
        "REV_AT_EPOCH":43836,
        "BSTAR":0.00031432,
        "MEAN_MOTION_DOT":0.00017286,
        "MEAN_MOTION_DDOT":0
    }
    oe_iss: OrbitalElements = OrbitalElements.from_celestrak_json(oe_09fev24)
    iss = Satellite(SGP4(oe_iss), "ISS")
    t0 = Time.fromisoformat("2024-02-09T08:15:00")
    t1: Time = Time.fromisoformat("2024-02-09T08:30:00")
    tl = Timeline(t0, t1, Q_("1s"))
    gre_gps = GPSCoordinates(
        longitude_deg=5.71667,
        latitude_deg=45.166672,
        location_name="Grenoble (FR)"
    )

    # https://in-the-sky.org/satpasses.php?town=3014728
    # FEB 09 2024
    # ######	Time    	Dir	Alt	Mag
    # START:	08:17:15	SSW	10°	0.4
    # MAX  :	08:20:37	ESE	24°	4.8	
    # END  :	08:22:58	E	10°	1.6

    iss_visible_gre = VisibleFromLocation(iss, gre_gps, Q_("10°"))
    event, *_ = iss_visible_gre.compute_intervals(tl)

    assert event.start == Time.fromisoformat("2024-02-09T08:17:15")
    assert event.stop == Time.fromisoformat("2024-02-09T08:22:57")

    c0 = iss.coordinates(event.start)
    hor0 = iss_visible_gre.local_frame.get_horizontal_coordinates(c0)
    assert hor0.azimuth_deg == approx(-164, abs=2)
    assert hor0.altitude_deg == approx(10, abs=1)

    c1 = iss.coordinates(event.stop)
    hor1 = iss_visible_gre.local_frame.get_horizontal_coordinates(c1)
    assert hor1.azimuth_deg == approx(81, abs=2)
    assert hor1.altitude_deg == approx(10, abs=1)

    tmid = event._idx2time(206)
    cmid = iss.coordinates(tmid)
    hormid = iss_visible_gre.local_frame.get_horizontal_coordinates(cmid)
    assert hormid.azimuth_deg == approx(120, abs=2)
    assert hormid.altitude_deg == approx(24, abs=1)