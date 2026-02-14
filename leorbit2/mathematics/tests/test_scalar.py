from fractions import Fraction

import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar


def test_scalar_new_and_magnitude_and_cast():
    s = Scalar[D.Length].new(2000)
    assert s.magnitude("meter") == pytest.approx(2000)
    assert s.magnitude("kilo_meter") == pytest.approx(2)

    same = s.cast(D.Length)
    assert same is s

    with pytest.raises(RuntimeError):
        s.cast(D.Time)


def test_scalar_add_sub_same_dimension():
    a = 2 * Quantity.meter
    b = 500 * Quantity.meter
    assert (a + b).magnitude("meter") == pytest.approx(502)
    assert (a - b).magnitude("meter") == pytest.approx(-498)
    assert (b - a).magnitude("meter") == pytest.approx(498)


def test_scalar_add_raises_on_dimension_mismatch():
    with pytest.raises(RuntimeError):
        _ = Quantity.meter + Quantity.second


def test_scalar_mul_and_div_update_dimensions():
    distance = 3 * Quantity.meter
    duration = 2 * Quantity.second

    speed = distance / duration
    assert speed.dim_coords == (D.Length._d / D.Time._d)
    assert speed.magnitude() == pytest.approx(1.5)

    prod = distance * duration
    assert prod.dim_coords == (D.Length._d * D.Time._d)
    assert prod.magnitude() == pytest.approx(6.0)


def test_scalar_modulo_and_reverse_modulo_dimless():
    a = Scalar[D.Dimless].new(7)
    b = Scalar[D.Dimless].new(3)

    assert (a % b).magnitude() == pytest.approx(1)
    assert (7 % b).magnitude() == pytest.approx(1)


def test_scalar_modulo_raises_on_dimension_mismatch():
    with pytest.raises(RuntimeError):
        _ = Quantity.meter % Quantity.second


def test_scalar_comparisons():
    a = 2 * Quantity.second
    b = 3 * Quantity.second

    assert a < b
    assert a <= b
    assert b > a
    assert b >= a
    assert a == 2 * Quantity.second
    assert a != b


def test_scalar_pow_updates_dimension_and_values():
    length = 3 * Quantity.meter
    squared = length ** 2
    assert squared.dim_coords == (D.Length._d ** 2)
    assert squared.magnitude() == pytest.approx(9)

    sqrt_length = length ** Fraction(1, 2)
    assert sqrt_length.dim_coords == (D.Length._d ** Fraction(1, 2))
    assert sqrt_length.magnitude() == pytest.approx(3 ** 0.5)


def test_scalar_matmul_not_supported():
    s = Scalar[D.Dimless].new(1)
    with pytest.raises(RuntimeError):
        _ = s @ s
