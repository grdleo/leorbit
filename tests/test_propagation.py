from leorbit.orbit.orbital_elements import OrbitalElements
from leorbit.orbit.propagator import NoPropagator, SGP4
from leorbit.orbit.objects import Satellite
from leorbit.math import Q_

import pytest
from math import tau

from tests import CELESTRAK_JSON, qt_almost_equals

def test_propagation():
    """Testing the NoPropagator algorithm and the creation of OrbitalElements with `from_state_vectors`"""
    full_angle = Q_(tau, "rad")
    oe: OrbitalElements = OrbitalElements.from_celestrak_json(CELESTRAK_JSON)
    period = oe.period
    propa = NoPropagator(oe)
    sgp4_propa = SGP4(oe)
    sat = Satellite(propa, "idk")
    sgp4_sat = Satellite(sgp4_propa, "sgp4-idk")

    coe = oe.to_coordinates().to_gps()
    c0 = sgp4_sat.coordinates(oe.epoch).to_gps()
    assert coe.longitude_deg  == pytest.approx(c0.longitude_deg, abs=0.1)
    assert coe.latitude_deg  == pytest.approx(c0.latitude_deg, abs=0.1)
    # SGP4 corrects altitude in a weird way...
    assert coe.altitude_from_sea_level_meter == pytest.approx(c0.altitude_from_sea_level_meter, abs=10_000)

    maxk = 40
    for i in range(0, 50):
        k = i / maxk
        prop_period = period * k
        prop_angle = full_angle * k

        c = sat.coordinates(oe.epoch + prop_period)
        oee: OrbitalElements = OrbitalElements.from_state_vectors(c)

        assert qt_almost_equals(oe.eccentricity, oee.eccentricity, abso=1e-5)
        assert qt_almost_equals(oe.inclination, oee.inclination, abso=1e-5)
        assert qt_almost_equals(oe.mean_motion, oee.mean_motion, abso=1e-5)
        assert qt_almost_equals(oe.arg_of_pericenter, oee.arg_of_pericenter, abso=1e-2)
        assert qt_almost_equals(oe.ra_of_asc_node, oee.ra_of_asc_node, abso=1e-5)
        assert qt_almost_equals((oe.mean_anomaly + prop_angle) % full_angle, oee.mean_anomaly, abso=1e-2)