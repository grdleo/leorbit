from fractions import Fraction

import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.scalar_array import ScalarArray


def test_scalar_array_new_size_getitem():
    arr = ScalarArray[D.Length].new([1, 2, 3])
    assert arr.size == 3
    assert arr[0].magnitude("meter") == pytest.approx(1)
    assert arr[2].magnitude("meter") == pytest.approx(3)
    with pytest.raises(KeyError):
        _ = arr[3]


def test_scalar_array_add_sub_mul_div():
    a = ScalarArray[D.Length].new([1, 2, 3])
    b = ScalarArray[D.Length].new([4, 5, 6])

    np.testing.assert_allclose((a + b).magnitude(), [5, 7, 9])
    np.testing.assert_allclose((b - a).magnitude(), [3, 3, 3])

    scaled = a * Scalar[D.Dimless].new(2)
    np.testing.assert_allclose(scaled.magnitude(), [2, 4, 6])

    speed = a / Scalar[D.Time].new(2)
    assert speed.dim_coords == (D.Length._d / D.Time._d)
    np.testing.assert_allclose(speed.magnitude(), [0.5, 1.0, 1.5])


def test_scalar_array_modulo_and_reverse_modulo_dimless():
    a = ScalarArray[D.Dimless].new([7, 8, 9])
    b = ScalarArray[D.Dimless].new([2, 3, 4])

    np.testing.assert_allclose((a % b).magnitude(), [1, 2, 1])
    np.testing.assert_allclose((10 % b).magnitude(), [0, 1, 2])


def test_scalar_array_comparisons():
    a = ScalarArray[D.Dimless].new([1, 2, 3])
    b = ScalarArray[D.Dimless].new([2, 2, 1])

    np.testing.assert_array_equal(a < b, [True, False, False])
    np.testing.assert_array_equal(a <= b, [True, True, False])
    np.testing.assert_array_equal(a > b, [False, False, True])
    np.testing.assert_array_equal(a >= b, [False, True, True])
    np.testing.assert_array_equal(a == a, [True, True, True])
    np.testing.assert_array_equal(a.__neq__(b), [True, False, True])


def test_scalar_array_pow_updates_dimension_and_values():
    arr = ScalarArray[D.Length].new([1, 4, 9])

    squared = arr ** 2
    assert squared.dim_coords == (D.Length._d ** 2)
    np.testing.assert_allclose(squared.magnitude(), [1, 16, 81])

    rooted = arr ** Fraction(1, 2)
    assert rooted.dim_coords == (D.Length._d ** Fraction(1, 2))
    np.testing.assert_allclose(rooted.magnitude(), [1, 2, 3])


def test_scalar_array_matmul_not_supported():
    arr = ScalarArray[D.Dimless].new([1, 2])
    with pytest.raises(RuntimeError):
        _ = arr @ arr


def test_scalar_array_dimension_mismatch_raises():
    arr = ScalarArray[D.Length].new([1, 2])
    with pytest.raises(RuntimeError):
        _ = arr + Scalar[D.Time].new(1)
    with pytest.raises(RuntimeError):
        _ = arr + Quantity.second
