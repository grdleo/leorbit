import numpy as np
import pytest

from leorbit.mathematics import (
    Dimless,
    Length,
    Quantity,
    Tensor,
    TensorAsMatrix33,
    TensorKind,
    ensure_same_dimensions,
    mat33,
    scalar,
    vector3,
)


def test_scalar_creation_and_kind_inference():
    a = scalar(123)
    assert a.kind == TensorKind.SCALAR
    assert a.phy_dimension == Dimless
    assert a.shape == (1,)
    assert a.size == 1
    assert a.scalar.value() == pytest.approx(123.0)


def test_dimension_composition_with_units():
    distance = 123.0 * Quantity.meter
    duration = 123.0 * Quantity.second
    speed = (distance / duration) * 10.0

    assert speed.kind == TensorKind.SCALAR
    assert speed.phy_dimension == Length / Quantity.second.phy_dimension
    assert speed.scalar.value(Quantity.kilo_meter_per_hour) == pytest.approx(36.0)


def test_tensor_concatenation_with_or_operator():
    a = vector3(1.0, 2.0, 3.0)
    b = vector3(4.0, 5.0, 6.0)
    c = a | b

    assert c.kind == TensorKind.VECTOR3
    assert c.size == 2
    np.testing.assert_allclose(c.raw_data_array(), np.array([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]]))


def test_tensor_getitem_preserves_kind_and_size_one():
    v = vector3(10.0, 20.0, 30.0) | vector3(40.0, 50.0, 60.0)
    one = v[1]

    assert one.kind == TensorKind.VECTOR3
    assert one.size == 1
    np.testing.assert_allclose(one.raw_data_array().reshape(3), [40.0, 50.0, 60.0])


def test_vector_normalized_is_dimensionless():
    v = vector3(3.0, 4.0, 0.0) * Quantity.meter
    n = v.vector3.normalized()

    assert n.phy_dimension == Dimless
    assert n.kind == TensorKind.VECTOR3
    assert n.shape == (3, 1)
    np.testing.assert_allclose(n.raw_data_array().reshape(3), [0.6, 0.8, 0.0], atol=1e-12)


def test_matrix_from_elements_and_inverse_roundtrip():
    m = TensorAsMatrix33.from_elements(
        2.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 5.0,
    )
    inv = m.matrix33.inverse()

    np.testing.assert_allclose(m.raw_data_array().reshape(3, 3), np.diag([2.0, 3.0, 5.0]))
    np.testing.assert_allclose(inv.raw_data_array().reshape(3, 3), np.diag([0.5, 1 / 3, 0.2]), atol=1e-12)


def test_raw_data_array_with_string_and_tensor_units():
    d = 1500.0 * Quantity.meter
    assert d.scalar.value("kilo_meter") == pytest.approx(1.5)
    assert d.scalar.value(Quantity.kilo_meter) == pytest.approx(1.5)


def test_ensure_same_dimensions_short_circuit_for_single_arg():
    t = 42.0 * Quantity.second
    assert ensure_same_dimensions(t)


def test_invalid_unit_lookup_raises_value_error():
    with pytest.raises(ValueError):
        Quantity.get("this_unit_does_not_exist")


def test_matrix_factory_kind_and_shape():
    m = mat33(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
    )
    assert m.kind == TensorKind.MATRIX33
    assert m.shape == (3, 3, 1)

