from typing import Annotated

import numpy as np
import pytest

from leorbit.mathematics import (
    U,
    Tensor,
    TensorAsVector3,
    TensorBound,
    TensorKind,
    acos,
    asin,
    atan,
    atan2,
    cos,
    ensure_same_dimensions,
    ensure_tensor,
    interpolate,
    matrix33,
    normalize_angle,
    normalize_angle_symmetric,
    scalar,
    sin,
    tan,
    tensor_check,
    vector3,
)

Angle = U.radian
Dimless = U.dimensionless
Length = U.meter


def test_quantity_registry_and_synonyms():
    assert Tensor(1.0, U.parse_units("meter")).scalar.value("meter") == pytest.approx(1.0)
    assert Tensor(1.0, U.parse_units("m")).scalar.value("meter") == pytest.approx(1.0)
    assert Tensor(1.0, U.parse_units("kilometer")).scalar.value("meter") == pytest.approx(1000.0)
    assert Tensor(1.0, U.parse_units("km")).scalar.value("meter") == pytest.approx(1000.0)

def test_matmul():
    m = matrix33(
        0, 0, 1,
        1, 0, 0,
        0, 1, 0
    )
    v = vector3(1, 2, 3)
    assert m.matrix33 @ v.vector3 == vector3(3, 1, 2)

    m = matrix33(
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    ) * U.meter
    v = vector3(11, 22, 33) / (U.second ** 2)

    assert m.matrix33 @ v.vector3 == vector3(154, 352, 550) * (U.meter / U.second ** 2)


def test_matmul_matrix_size4_with_vector3_size4():
    m1 = matrix33(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    )
    m2 = matrix33(
        2, 0, 0,
        0, 3, 0,
        0, 0, 4,
    )
    m3 = matrix33(
        1, 2, 3,
        0, 1, 0,
        0, 0, 1,
    )
    m4 = matrix33(
        0, 1, 0,
        -1, 0, 0,
        0, 0, 1,
    )

    v1 = vector3(1, 2, 3)
    v2 = vector3(1, 2, 3)
    v3 = vector3(1, 2, 3)
    v4 = vector3(7, 8, 9)

    m_batch = m1 | m2 | m3 | m4
    v_batch = v1 | v2 | v3 | v4

    out = m_batch.matrix33 @ v_batch.vector3

    assert out.kind == TensorKind.VECTOR3
    assert out.size == 4
    assert out.shape == (3, 4)

    np.testing.assert_allclose(
        out.raw_data_array(out.units),
        np.array([
            [1.0, 2.0, 14.0, 8.0],
            [2.0, 6.0, 2.0, -7.0],
            [3.0, 12.0, 3.0, 9.0],
        ]),
    )


def test_matmul_matrix_by_matrix_size1_and_size5():
    left_1 = matrix33(
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
    )
    right_1 = matrix33(
        9, 8, 7,
        6, 5, 4,
        3, 2, 1,
    )

    out_1 = left_1.matrix33 @ right_1.matrix33

    assert out_1.kind == TensorKind.MATRIX33
    assert out_1.size == 1
    assert out_1.shape == (3, 3, 1)
    np.testing.assert_allclose(
        out_1.raw_data_array(out_1.units).reshape(3, 3),
        np.array([
            [30.0, 24.0, 18.0],
            [84.0, 69.0, 54.0],
            [138.0, 114.0, 90.0],
        ]),
    )

    left_batch = (
        matrix33(1, 0, 0, 0, 1, 0, 0, 0, 1)
        | matrix33(2, 0, 0, 0, 2, 0, 0, 0, 2)
        | matrix33(1, 2, 0, 0, 1, 0, 0, 0, 1)
        | matrix33(0, -1, 0, 1, 0, 0, 0, 0, 1)
        | matrix33(3, 1, 0, 0, 3, 0, 0, 0, 3)
    )
    right_batch = (
        matrix33(1, 2, 3, 4, 5, 6, 7, 8, 9)
        | matrix33(1, 0, 0, 0, 1, 0, 0, 0, 1)
        | matrix33(2, 0, 0, 0, 2, 0, 0, 0, 2)
        | matrix33(1, 0, 0, 0, 1, 0, 0, 0, 1)
        | matrix33(0, 1, 0, 1, 0, 0, 0, 0, 1)
    )

    out_5 = left_batch.matrix33 @ right_batch.matrix33

    assert out_5.kind == TensorKind.MATRIX33
    assert out_5.size == 5
    assert out_5.shape == (3, 3, 5)
    np.testing.assert_allclose(
        out_5.raw_data_array(out_5.units),
        np.stack(
            [
                np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]),
                np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
                np.array([[2.0, 4.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
                np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
                np.array([[1.0, 3.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 3.0]]),
            ],
            axis=2,
        ),
    )

def test_scalar_vector_matrix_factories_and_kinds():
    s = scalar(7)
    v = vector3(1.0, 2.0, 3.0)
    m = matrix33(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )

    assert s.kind == TensorKind.SCALAR
    assert s.shape == (1,)
    assert s.size == 1

    assert v.kind == TensorKind.VECTOR3
    assert v.shape == (3, 1)
    assert v.size == 1

    assert m.kind == TensorKind.MATRIX33
    assert m.shape == (3, 3, 1)
    assert m.size == 1


def test_scalar_value_and_unit_conversion():
    speed = scalar(36.0).with_units(U.kilometer / U.hour)
    assert speed.scalar.value("kilometer/hour") == pytest.approx(36.0)
    assert speed.scalar.value("meter/second") == pytest.approx(10.0)


def test_tensor_check_and_secure():
    d = scalar(5.0).with_units(U.meter)
    assert d.check(units=Length)
    assert d.check(kind=TensorKind.SCALAR)
    assert d.secure(units=Length, kind=TensorKind.SCALAR) is d

    with pytest.raises(ValueError):
        d.secure(units=Angle)
    assert d.check() is True


def test_ensure_tensor_and_dimension_guard():
    from_float = ensure_tensor(3.5)
    assert from_float.check(units=Dimless, kind=TensorKind.SCALAR)
    assert from_float.scalar.value("dimensionless") == pytest.approx(3.5)

    t1 = scalar(1.0).with_units(U.meter)
    t2 = scalar(2.0).with_units(U.meter)
    assert ensure_same_dimensions(t1, t2)

    with pytest.raises(RuntimeError):
        ensure_same_dimensions(t1, scalar(1.0).with_units(U.second))


def test_vector_components_norm_dot_cross_and_angle():
    v = vector3(3.0, 4.0, 0.0) * U.meter
    w = vector3(0.0, 4.0, 3.0) * U.meter

    assert v.vector3.x.scalar.value("meter") == pytest.approx(3.0)
    assert v.vector3.y.scalar.value("meter") == pytest.approx(4.0)
    assert v.vector3.z.scalar.value("meter") == pytest.approx(0.0)
    assert v.vector3.length.scalar.values("meter") == pytest.approx([5.0])

    dot = v.vector3.dot(w)
    assert dot.units.is_compatible_with(Length * Length)
    assert dot.scalar.value(dot.units) == pytest.approx(16.0)

    cross = v.vector3.cross(w)
    np.testing.assert_allclose(cross.raw_data_array("meter**2").reshape(3), [12.0, -9.0, 12.0])

    angle = v.vector3.angle(v)
    assert angle.check(units=Angle, kind=TensorKind.SCALAR)
    assert angle.scalar.value("rad") == pytest.approx(0.0)


def test_vector_spherical_roundtrip():
    theta = scalar(np.pi / 4.0).with_units(U.radian)
    delta = scalar(0.0).with_units(U.radian)
    radius = scalar(np.sqrt(2.0)).with_units(U.meter)

    v = TensorAsVector3.from_spherical(theta=theta, delta=delta, radius=radius)
    np.testing.assert_allclose(v.raw_data_array("meter").reshape(3), [1.0, 1.0, 0.0], atol=1e-12)


def test_matrix_inverse_and_product():
    m = matrix33(
        2.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 4.0,
    )
    inv = m.matrix33.inverse()

    np.testing.assert_allclose(
        inv.raw_data_array(inv.units).reshape(3, 3),
        np.diag([0.5, 1.0 / 3.0, 0.25]),
        atol=1e-12,
    )


def test_trigonometric_helpers_and_domains():
    right_angle = scalar(np.pi / 2.0).with_units(U.radian)

    assert sin(right_angle).scalar.value("dimensionless") == pytest.approx(1.0)
    assert cos(right_angle).scalar.value("dimensionless") == pytest.approx(0.0, abs=1e-12)
    assert tan(scalar(0.0).with_units(U.radian)).scalar.value("dimensionless") == pytest.approx(0.0)

    u = scalar(0.5).with_units(U.dimensionless)
    assert asin(u).check(units=Angle)
    assert acos(u).check(units=Angle)
    assert atan(u).check(units=Angle)

    with pytest.raises(ValueError):
        sin(scalar(1.0).with_units(U.meter))


def test_atan2_normalization_and_interpolation():
    y = scalar(1.0).with_units(U.meter)
    x = scalar(1.0).with_units(U.meter)
    assert atan2(y, x).scalar.value("rad") == pytest.approx(np.pi / 4.0)

    wrapped = normalize_angle(scalar(5.0 * np.pi).with_units(U.radian))
    assert wrapped.scalar.value("rad") == pytest.approx(np.pi)

    symmetric = normalize_angle_symmetric(scalar(1.5 * np.pi).with_units(U.radian))
    assert symmetric.scalar.value("rad") == pytest.approx(-np.pi / 2.0)

    a = scalar(10.0).with_units(U.meter)
    b = scalar(14.0).with_units(U.meter)
    p = interpolate(a, b, 0.25)
    assert p.scalar.value("meter") == pytest.approx(11.0)


def test_tensor_check_decorator_rejects_bad_argument_and_return():
    @tensor_check
    def scale_length_vector(
        v: Annotated[Tensor, TensorBound(units=Length, kind=TensorKind.VECTOR3)],
        k: float,
    ) -> Annotated[Tensor, TensorBound(units=Length, kind=TensorKind.VECTOR3)]:
        return v * k

    good = vector3(1.0, -2.0, 3.0) * U.meter
    out = scale_length_vector(good, 2.0)
    np.testing.assert_allclose(out.raw_data_array("meter").reshape(3), [2.0, -4.0, 6.0])

    with pytest.raises(ValueError):
        scale_length_vector(scalar(1.0).with_units(U.meter), 2.0)

    @tensor_check
    def broken_return(
        v: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)],
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return v

    with pytest.raises(ValueError):
        broken_return(vector3(1.0, 2.0, 3.0))
