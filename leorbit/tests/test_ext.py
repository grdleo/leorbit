import json
import os
from pathlib import Path
import time

import pytest

from leorbit.ext import CelestrakDataGP, get_celestrak_gpdata


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


def test_get_celestrak_gpdata_uses_fresh_cache(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setenv("LEORBIT_CACHE_DIR", str(tmp_path))
    monkeypatch.setenv("LEORBIT_CELESTRAK_CACHE_HOURS", "12")

    catnr = 25544
    cache_file = tmp_path / f"python_leorbit_celestrak_gpdata_{catnr}.json"
    cache_file.write_text(CelestrakDataGP(**make_sample()).model_dump_json())

    now = time.time()
    os.utime(cache_file, (now, now))

    def _should_not_call_network(*args, **kwargs):
        raise AssertionError("Network fetch should not be called when cache is fresh")

    monkeypatch.setattr("leorbit.ext.get", _should_not_call_network)

    obj = get_celestrak_gpdata(catnr)
    assert obj.norad_cat_id == 25544


def test_get_celestrak_gpdata_refreshes_stale_cache(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setenv("LEORBIT_CACHE_DIR", str(tmp_path))
    monkeypatch.setenv("LEORBIT_CELESTRAK_CACHE_HOURS", "1")

    catnr = 25544
    cache_file = tmp_path / f"python_leorbit_celestrak_gpdata_{catnr}.json"
    cache_file.write_text(CelestrakDataGP(**make_sample()).model_dump_json())

    stale = time.time() - 2 * 3600
    os.utime(cache_file, (stale, stale))

    sample = make_sample()
    sample["OBJECT_NAME"] = "ISS (ZARYA)"

    class _Response:
        ok = True
        status_code = 200
        reason = "OK"
        text = "ok"

        @staticmethod
        def json():
            return [sample]

    called = {"value": False}

    def _fake_get(*args, **kwargs):
        called["value"] = True
        return _Response()

    monkeypatch.setattr("leorbit.ext.get", _fake_get)

    obj = get_celestrak_gpdata(catnr)
    assert called["value"] is True
    assert obj.name == "ISS (ZARYA)"
