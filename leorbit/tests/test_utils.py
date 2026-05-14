import numpy as np
import pytest

from leorbit.mathematics import U, TensorKind, scalar, vector3
import leorbit.utils as u



@pytest.mark.parametrize(
    "truth_array, expected",
    [
        ([False, True, True, True, False, False, True, True, False], [(1, 3), (6, 7)]),
        ([False, False, False], []),
        ([True, True, True], [(0, 2)]),
        ([True], [(0, 0)]),
        ([False], []),
        ([False, True, False], [(1, 1)]),
        ([True, False, True], [(0, 0), (2, 2)]),
        ([False, True, True, False, True, False, True, True, True, False], [(1, 2), (4, 4), (6, 8)]),
    ],
)
def test_truth_array_to_indices_intervals_various_patterns(truth_array, expected):
    assert u.truth_array_to_indices_intervals(truth_array) == expected


@pytest.mark.parametrize(
    "truth_array, min_size, expected",
    [
        ([False, True, True, True, False, False, True, True, False], 1, [(1, 3), (6, 7)]),
        ([False, True, True, True, False, False, True, True, False], 2, [(1, 3), (6, 7)]),
        ([False, True, True, True, False, False, True, True, False], 3, [(1, 3)]),
        ([False, True, True, True, False, False, True, True, False], 4, []),
        ([True, True, False, True, True, True, False, True], 2, [(0, 1), (3, 5)]),
        ([True, True, False, True, True, True, False, True], 3, [(3, 5)]),
        ([True, True, False, True, True, True, False, True], 10, []),
    ],
)
def test_truth_array_to_indices_intervals_min_size_filtering(truth_array, min_size, expected):
    assert u.truth_array_to_indices_intervals(truth_array, min_size_intervals=min_size) == expected


@pytest.mark.parametrize(
    "arr",
    [
        np.array([0, 1, 1, 1, 0, 1, 1, 0, 1], dtype=bool),
        np.array([0, 1, 1, 1, 0, 1, 1, 0, 1], dtype=np.int8),
        np.array([0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=np.float64),
        np.array([[False, True, True], [False, True, False]], dtype=bool),
    ],
)
def test_truth_array_to_indices_intervals_numpy_inputs(arr):
    got = u.truth_array_to_indices_intervals(arr, min_size_intervals=2)
    ref = u.truth_array_to_indices_intervals(np.asarray(arr, dtype=bool).reshape(-1), min_size_intervals=2)
    assert got == ref


@pytest.mark.parametrize(
    "iterable_input, expected",
    [
        ((x for x in [False, True, True, False]), [(1, 2)]),
        ((x for x in [True, False, True, True]), [(0, 0), (2, 3)]),
        (tuple([False, False, True, True, True]), [(2, 4)]),
    ],
)
def test_truth_array_to_indices_intervals_iterable_inputs(iterable_input, expected):
    assert u.truth_array_to_indices_intervals(iterable_input) == expected


def test_truth_array_to_indices_intervals_empty_and_invalid_min_size():
    assert u.truth_array_to_indices_intervals([]) == []

    assert u.truth_array_to_indices_intervals(np.array([], dtype=bool)) == []

    with pytest.raises(ValueError):
        u.truth_array_to_indices_intervals([True, False], min_size_intervals=0)

    with pytest.raises(ValueError):
        u.truth_array_to_indices_intervals([True, False], min_size_intervals=-3)


def test_angle_and_time_helpers():
    assert u.angle2dms(scalar(39.5).with_units(U.degree)) == " 039° 30′ 00″"

    out = u.unixepoch_to_j2000(np.array([0.0, 86_400.0]))
    np.testing.assert_allclose(out, np.array([-10_957.5, -10_956.5]))

    stl0 = u.j2000_to_stl0(0.0)
    assert np.isfinite(np.asarray(stl0)).all()


def test_dms2angle_converts_correctly():
    angle = u.dms2angle(39.0, 30.0, 0.0)
    expected = np.deg2rad(39.5)

    assert angle.check(units=U.radian, kind=TensorKind.SCALAR)
    assert angle.scalar.value("radian") == pytest.approx(expected)


@pytest.mark.parametrize(
    "unixepoch",
    [
        0.0,
        86_400.0,
        np.array([0.0, 86_400.0, 172_800.0], dtype=np.float64),
    ],
)
def test_jd_parametrized(unixepoch):
    got = u.jd(unixepoch)
    expected = np.asarray(unixepoch, dtype=np.float64) / 86_400 + 2_440_587.5

    assert got.check(units=U.second, kind=TensorKind.SCALAR)
    np.testing.assert_allclose(got.raw_data_array("day"), expected)


@pytest.mark.parametrize(
    "unixepoch",
    [
        0.0,
        86_400.0,
        np.array([0.0, 86_400.0, 172_800.0], dtype=np.float64),
    ],
)
def test_j2000_parametrized(unixepoch):
    got = u.j2000(unixepoch)
    expected = np.asarray(unixepoch, dtype=np.float64) / 86_400 - 10_957.5

    assert got.check(units=U.second, kind=TensorKind.SCALAR)
    np.testing.assert_allclose(got.raw_data_array("day"), expected)


@pytest.mark.parametrize(
    "unixepoch",
    [
        0.0,
        86_400.0,
        np.array([0.0, 86_400.0, 172_800.0], dtype=np.float64),
    ],
)
def test_from_mil_parametrized(unixepoch):
    got = u.from_mil(unixepoch)
    expected = np.asarray(unixepoch, dtype=np.float64) / 86_400 - 10_958.0

    assert got.check(units=U.second, kind=TensorKind.SCALAR)
    np.testing.assert_allclose(got.raw_data_array("day"), expected)


@pytest.mark.parametrize(
    "unixepoch",
    [
        0.0,
        86_400.0,
        np.array([0.0, 86_400.0, 172_800.0], dtype=np.float64),
    ],
)
def test_stl0_parametrized(unixepoch):
    got = u.stl0(unixepoch)
    j2k = np.asarray(unixepoch, dtype=np.float64) / 86_400 - 10_957.5
    expected = np.asarray(u.j2000_to_stl0(j2k), dtype=np.float64)

    assert got.check(units=U.radian, kind=TensorKind.SCALAR)
    assert np.isfinite(got.raw_data_array("radian")).all()
    np.testing.assert_allclose(got.raw_data_array("radian"), expected)


def test_orbital_scalar_conversions_and_contracts():
    sma = scalar(7_000_000).with_units(U.meter)
    n = u.semi_major_axis_earth_to_mean_motion(sma)
    assert n.check(units=U.radian / U.second, kind=TensorKind.SCALAR)

    sma_rt = u.mean_motion_to_semi_major_axis_earth(n)
    assert sma_rt.check(units=U.meter, kind=TensorKind.SCALAR)

    with pytest.raises((TypeError, ValueError)):
        u.mean_motion_to_semi_major_axis_earth(scalar(1.0).with_units(U.meter))


def test_anomaly_helpers_and_contracts():
    e = scalar(0.01).with_units(U.dimensionless)
    m = scalar(0.2).with_units(U.radian)

    nu = u.mean2true_anomaly(e, m)
    assert nu.check(units=U.radian, kind=TensorKind.SCALAR)

    ecc = u.mean2eccentric_anomaly(e, m)
    assert ecc.check(units=U.radian, kind=TensorKind.SCALAR)

    nu2 = u.eccentric2true_anomaly(e, ecc)
    assert nu2.check(units=U.radian, kind=TensorKind.SCALAR)

    with pytest.raises((TypeError, ValueError)):
        u.mean2true_anomaly(vector3(0.01, 0.0, 0.0), m)


def test_itrf_to_gps_contract():
    pos = vector3(6_378_135.0, 0.0, 0.0) * U.meter
    gps = u.itrf2gps(pos)

    assert gps.latitude.check(units=U.radian, kind=TensorKind.SCALAR)
    assert gps.longitude.check(units=U.radian, kind=TensorKind.SCALAR)
    assert gps.altitude.check(units=U.meter, kind=TensorKind.SCALAR)

    with pytest.raises((TypeError, ValueError)):
        u.itrf2gps(vector3(1.0, 0.0, 0.0))
