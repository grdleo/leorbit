import json

import pytest

from leorbit.ext import CelestrakDataGP


def make_sample():
    # minimal example using alias keys like the API returns
    return {
        "EPOCH": "2024-01-30T13:06:52.283808",
        "ECCENTRICITY": 0.0002899,
        "INCLINATION": 51.641,
        "RA_OF_ASC_NODE": 282.3593,
        "ARG_OF_PERICENTER": 174.8339,
        "MEAN_MOTION": 15.49383526,
        "MEAN_ANOMALY": 252.8219,
        "BSTAR": 0.00053423,
        "NORAD_CAT_ID": 25544,
    }


def test_roundtrip_dict():
    data = make_sample()
    obj = CelestrakDataGP(**data)

    # produce a dump with *field names* (default behaviour)
    dumped = obj.model_dump()
    # keys should be lower-case attribute names, not the aliases
    assert "epoch" in dumped and "EPOCH" not in dumped

    # re-create using the dumped dict; this works because of populate_by_name
    obj2 = CelestrakDataGP(**dumped)
    assert obj2 == obj


def test_roundtrip_json():
    data = make_sample()
    obj = CelestrakDataGP(**data)

    # if we want aliases in the JSON, supply by_alias
    j = obj.model_dump_json(by_alias=True)
    # parse back using the pydantic helper
    obj3 = CelestrakDataGP.model_validate_json(j)
    assert obj3 == obj

    # and the plain JSON (no alias) also round-trips thanks to the config
    j2 = obj.model_dump_json()
    obj4 = CelestrakDataGP.model_validate_json(j2)
    assert obj4 == obj
