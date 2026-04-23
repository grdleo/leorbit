from typing import Annotated

import numpy as np
import pytest

from leorbit.mathematics import (
    Angle,
    Dimless,
    Length,
    Quantity,
    Tensor,
    TensorBound,
    TensorKind,
    acos,
    asin,
    atan,
    atan2,
    cbrt,
    cos,
    cube,
    ensure_same_dimensions,
    ensure_tensor,
    interpolate,
    mat33,
    normalize_angle,
    normalize_angle_symmetric,
    scalar,
    scalar_array,
    sin,
    sqrt,
    square,
    tan,
    tensor_check,
    vec3,
    vec3_array,
)


def test_quantity_and_unit_conversion():
    assert Quantity.meter.check(Length)
    assert Quantity.second.check(dimension=Quantity.second.phy_dimension)
    assert Quantity.get("kilo_meter").magnitude("meter") == pytest.approx(1000)
    assert (2 * Quantity.kilo_meter).magnitude("meter") == pytest.approx(2000)


def test_tensor_factories_and_kinds():
    s = scalar(3)
    sa = scalar_array([1.0, 2.0, 3.0])
    v = vec3(x=1.0, y=2.0, z=3.0)
    va = vec3_array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    m = mat33(
        a11=1.0, a21=0.0, a31=0.0,
        a12=0.0, a22=1.0, a32=0.0,
        a13=0.0, a23=0.0, a33=1.0,
    )

    assert s.kind == TensorKind.SCALAR and s.size == 1
    assert sa.kind == TensorKind.SCALAR and sa.size == 3
    assert v.kind == TensorKind.VECTOR3 and v.size == 1
    assert va.kind == TensorKind.VECTOR3 and va.size == 2
    assert m.kind == TensorKind.MATRIX33 and m.size == 1


def test_dimension_composition_and_cast():
    length = 2 * Quantity.meter
    time = 4 * Quantity.second
    speed = length / time

    assert speed.phy_dimension.triplet() == (Length / Quantity.second.phy_dimension).triplet()

    with pytest.raises(RuntimeError):
        length.cast(Angle)


def test_basic_vector_algebra():
    v = vec3(x=1.0, y=2.0, z=3.0) * Quantity.meter
    w = vec3(x=3.0, y=2.0, z=1.0) * Quantity.meter

    np.testing.assert_allclose((v + w)._values.flatten(), [4.0, 4.0, 4.0])
    np.testing.assert_allclose(v.cross(w)._values.flatten(), [-4.0, 8.0, -4.0])
    assert v.dot(w).magnitude("meter") == pytest.approx(10.0)


def test_vector_angles_and_spherical_helpers():
    v = vec3(x=1.0, y=1.0, z=0.0) * Quantity.meter
    assert v.theta.check(Angle)
    assert v.theta.magnitude("rad") == pytest.approx(np.pi / 4)
    assert v.delta.magnitude("rad") == pytest.approx(0.0)

    u = Tensor.from_spherical(
        theta=np.pi / 4 * Quantity.radian,
        delta=0 * Quantity.radian,
        radius=np.sqrt(2) * Quantity.meter,
    )
    np.testing.assert_allclose(u._values.flatten(), [1.0, 1.0, 0.0], atol=1e-12)


def test_matrix_products_with_vector_arrays():
    vectors = vec3_array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    m = Tensor.from_elements(
        a11=2.0, a21=0.0, a31=0.0,
        a12=0.0, a22=3.0, a32=0.0,
        a13=0.0, a23=0.0, a33=4.0,
    )

    out = m @ vectors
    np.testing.assert_allclose(out._values, np.array([[2.0, 8.0], [6.0, 15.0], [12.0, 24.0]]))


def test_ensure_helpers():
    t = ensure_tensor(3.5)
    assert t.check(Dimless)

    a = 1 * Quantity.meter
    b = 2 * Quantity.meter
    assert ensure_same_dimensions(a, b)

    with pytest.raises(RuntimeError):
        ensure_same_dimensions(1 * Quantity.meter, 1 * Quantity.second)


def test_trig_and_power_helpers():
    ang = np.pi / 2 * Quantity.radian

    assert sin(ang).magnitude() == pytest.approx(1.0)
    assert cos(ang).magnitude() == pytest.approx(0.0, abs=1e-12)
    assert tan(0 * Quantity.radian).magnitude() == pytest.approx(0.0)

    unit = 0.5 * Quantity.dimensionless
    assert asin(unit).check(Angle)
    assert acos(unit).check(Angle)
    assert atan(unit).check(Angle)

    y = 1 * Quantity.meter
    x = 1 * Quantity.meter
    assert atan2(y, x).magnitude("rad") == pytest.approx(np.pi / 4)

    sq = square(3 * Quantity.meter)
    assert sq.magnitude() == pytest.approx(9.0)
    assert sqrt(sq).magnitude("meter") == pytest.approx(3.0)
    assert cube(4 * Quantity.meter).magnitude() == pytest.approx(64.0)
    assert cbrt(64 * (Quantity.meter ** 3)).magnitude("meter") == pytest.approx(4.0)


def test_angle_normalization_and_interpolation():
    n = normalize_angle(5 * np.pi * Quantity.radian)
    ns = normalize_angle_symmetric(3 * np.pi / 2 * Quantity.radian)

    assert n.magnitude("rad") == pytest.approx(np.pi)
    assert ns.magnitude("rad") == pytest.approx(-np.pi / 2)

    a = 10 * Quantity.meter
    b = 14 * Quantity.meter
    assert interpolate(a, b, 0.25).magnitude("meter") == pytest.approx(11.0)


def test_tensor_check_with_runtime_bounds():
    @tensor_check
    def scale_vector(
        value: Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)],
        gain: float,
    ) -> Annotated[Tensor, TensorBound(dimension=Length, kind=TensorKind.VECTOR3)]:
        return value * gain

    v = vec3(x=1.0, y=-2.0, z=3.0) * Quantity.meter
    out = scale_vector(v, 2.0)
    np.testing.assert_allclose(out._values.flatten(), [2.0, -4.0, 6.0])

    with pytest.raises(ValueError):
        scale_vector(1 * Quantity.meter, 2.0)
