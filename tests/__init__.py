from pytest import approx
from leorbit.math import Q_

def qt_almost_equals(a: Q_, b: Q_, rela: float = None, abso: float = None):
    am = a.to_base_units().m
    bm = b.to_base_units().m
    return am == approx(bm, rela, abso)

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
    "MEAN_MOTION_DDOT": 0
}