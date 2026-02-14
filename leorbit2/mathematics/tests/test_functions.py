import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.functions import (
    HALF_REV,
    FULL_REV,
    acos,
    angle2dms,
    asin,
    atan,
    atan2,
    cos,
    normalize_angle,
    normalize_angle_symmetric,
    sin,
    sqrt,
    square,
    tan,
)
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.scalar_array import ScalarArray


def test_square_and_sqrt_scalar_and_scalar_array():
    s = 3 * Quantity.meter
    ss = square(s)
    assert ss.dim_coords == (D.Length._d ** 2)
    assert ss.magnitude() == pytest.approx(9)
    assert sqrt(ss).magnitude("meter") == pytest.approx(3)

    arr = ScalarArray[D.Length].new([1, 4, 9])
    arr_sq = square(arr)
    assert arr_sq.dim_coords == (D.Length._d ** 2)
    np.testing.assert_allclose(arr_sq.magnitude(), [1, 16, 81])
    np.testing.assert_allclose(sqrt(arr_sq).magnitude("meter"), [1, 4, 9])


def test_trigonometric_functions_scalar():
    angle = (np.pi / 2) * Quantity.rad
    assert sin(angle).magnitude() == pytest.approx(1)
    assert cos(angle).magnitude() == pytest.approx(0, abs=1e-12)
    assert tan(0 * Quantity.rad).magnitude() == pytest.approx(0)

    u = Scalar[D.Dimless].new(0.5)
    assert asin(u).check(D.Angle)
    assert acos(u).check(D.Angle)
    assert atan(u).check(D.Angle)


def test_trigonometric_functions_scalar_array_and_atan2():
    angles = ScalarArray[D.Angle].new([0, np.pi / 2])
    np.testing.assert_allclose(sin(angles).magnitude(), [0, 1], atol=1e-12)
    np.testing.assert_allclose(cos(angles).magnitude(), [1, 0], atol=1e-12)

    y = ScalarArray[D.Length].new([0, 1])
    x = ScalarArray[D.Length].new([1, 0])
    a = atan2(y, x)
    assert a.check(D.Angle)
    np.testing.assert_allclose(a.magnitude("rad"), [0, np.pi / 2], atol=1e-12)


def test_angle_normalization_helpers_and_constants():
    assert FULL_REV.magnitude("rad") == pytest.approx(2 * np.pi)
    assert HALF_REV.magnitude("rad") == pytest.approx(np.pi)

    angle = (5 * np.pi) * Quantity.rad
    normalized = normalize_angle(angle)
    assert normalized.magnitude("rad") == pytest.approx(np.pi)

    sym = normalize_angle_symmetric((3 * np.pi / 2) * Quantity.rad)
    assert sym.magnitude("rad") == pytest.approx(-np.pi / 2)


def test_angle2dms_format():
    out = angle2dms((1.5) * Quantity.deg)
    assert "°" in out
    assert "′" in out
    assert "″" in out
