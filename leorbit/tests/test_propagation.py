from math import tau

import numpy as np
import pytest

from leorbit.coordinates import OrbitalElements
from leorbit.ext import CelestrakDataGP
from leorbit.frames import AbsoluteFrame
from leorbit.m import Quantity, normalize_angle
from leorbit.propagator import NoPropagator, SGP4
from leorbit.sky_object import Satellite


CELESTRAK_JSON = {
    "OBJECT_NAME": "ISS (ZARYA)",
    "OBJECT_ID": "1998-067A",
    "EPOCH": "2024-01-30T13:06:52.283808",
    "MEAN_MOTION": 15.49383526,
    "ECCENTRICITY": 0.0002899,
    "INCLINATION": 51.641,
    "RA_OF_ASC_NODE": 282.3593,
    "ARG_OF_PERICENTER": 174.8339,
    "MEAN_ANOMALY": 252.8219,
    "EPHEMERIS_TYPE": 0,
    "CLASSIFICATION_TYPE": "U",
    "NORAD_CAT_ID": 25544,
    "ELEMENT_SET_NO": 999,
    "REV_AT_EPOCH": 43705,
    "BSTAR": 0.00053423,
    "MEAN_MOTION_DOT": 0.00029493,
    "MEAN_MOTION_DDOT": 0,
}


def _scalar_close(a, b, abs_tol: float) -> bool:
    return a.base_unit_value == pytest.approx(b.base_unit_value, abs=abs_tol)


def test_propagation():
    full_angle = tau * Quantity.rad
    oe: OrbitalElements = CelestrakDataGP(**CELESTRAK_JSON).to_orbital_elements()

    period = oe.period
    sat = Satellite("idk", oe, NoPropagator)
    sgp4_sat = Satellite("sgp4-idk", oe, SGP4)

    c_ref = oe.to_coordinates()
    c0 = sgp4_sat.coordinates(oe.epoch)

    p_ref = c_ref.get_pos(AbsoluteFrame.GCRF)
    p0 = c0.get_pos(AbsoluteFrame.GCRF)
    v_ref = c_ref.get_vel(AbsoluteFrame.GCRF)
    v0 = c0.get_vel(AbsoluteFrame.GCRF)

    np.testing.assert_allclose(p_ref._values.flatten(), p0._values.flatten(), atol=10_000)
    assert v_ref is not None and v0 is not None
    np.testing.assert_allclose(v_ref._values.flatten(), v0._values.flatten(), atol=100)

    maxk = 40
    for i in range(0, 50):
        k = i / maxk
        prop_period = period * k
        prop_angle = full_angle * k

        c = sat.coordinates(oe.epoch + prop_period)
        pos = c.get_pos(AbsoluteFrame.GCRF)
        vel = c.get_vel(AbsoluteFrame.GCRF)
        assert vel is not None

        oee = OrbitalElements.from_state_vectors(c.epoch, pos, vel)

        assert _scalar_close(oe.eccentricity, oee.eccentricity, abs_tol=1e-5)
        assert _scalar_close(oe.inclination, oee.inclination, abs_tol=1e-5)
        assert _scalar_close(oe.mean_motion, oee.mean_motion, abs_tol=1e-5)
        assert _scalar_close(oe.arg_of_pericenter, oee.arg_of_pericenter, abs_tol=1e-2)
        assert _scalar_close(oe.ra_of_asc_node, oee.ra_of_asc_node, abs_tol=1e-5)

        expected_m = normalize_angle(oe.mean_anomaly + prop_angle)
        assert _scalar_close(expected_m, oee.mean_anomaly, abs_tol=1e-2)
