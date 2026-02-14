import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.vector3 import Vector3


def test_vector3_named_constants():
    np.testing.assert_allclose(Vector3.O._values.flatten(), [0, 0, 0])
    np.testing.assert_allclose(Vector3.X._values.flatten(), [1, 0, 0])
    np.testing.assert_allclose(Vector3.Y._values.flatten(), [0, 1, 0])
    np.testing.assert_allclose(Vector3.Z._values.flatten(), [0, 0, 1])
    np.testing.assert_allclose(Vector3.ONE._values.flatten(), [1, 1, 1])


def test_vector3_new_properties_and_getters():
    v = Vector3[D.Length].new(1, 2, 2)
    assert v.x.magnitude("m") == pytest.approx(1)
    assert v.y.magnitude("m") == pytest.approx(2)
    assert v.z.magnitude("m") == pytest.approx(2)
    assert v.length.magnitude("m") == pytest.approx(3)


def test_vector3_new_from_scalars_requires_same_dim():
    with pytest.raises(RuntimeError):
        _ = Vector3.new(Quantity.m, Quantity.s, Quantity.m)


def test_vector3_add_sub_mul_div_and_dimension_changes():
    v = Vector3[D.Length].new(1, 2, 3)
    w = Vector3[D.Length].new(3, 2, 1)

    np.testing.assert_allclose((v + w)._values.flatten(), [4, 4, 4])
    np.testing.assert_allclose((v - w)._values.flatten(), [-2, 0, 2])

    scaled = v * Scalar[D.Dimless].new(2)
    np.testing.assert_allclose(scaled._values.flatten(), [2, 4, 6])

    speed = v / Scalar[D.Time].new(2)
    assert speed.dim_coords == (D.Length._d / D.Time._d)
    np.testing.assert_allclose(speed._values.flatten(), [0.5, 1.0, 1.5])


def test_vector3_dot_cross_matmul_and_normalized():
    x = Vector3.X
    y = Vector3.Y

    dot = x @ y
    assert dot.magnitude() == pytest.approx(0)

    cross = x.cross(y)
    np.testing.assert_allclose(cross._values.flatten(), [0, 0, 1])

    v = Vector3[D.Dimless].new(3, 0, 4)
    n = v.normalized()
    np.testing.assert_allclose(n._values.flatten(), [0.6, 0.0, 0.8])


def test_vector3_angle_and_spherical_conversion():
    x = Vector3.X
    y = Vector3.Y
    angle = x.angle(y)
    assert angle.magnitude("rad") == pytest.approx(np.pi / 2)

    theta = Scalar[D.Angle].new(0)
    delta = Scalar[D.Angle].new(0)
    rho = Scalar[D.Length].new(2)
    v = Vector3.from_spherical(theta, delta, rho)
    np.testing.assert_allclose(v._values.flatten(), [2, 0, 0], atol=1e-12)


def test_vector3_comparison_operators_raise():
    v = Vector3.X
    with pytest.raises(RuntimeError):
        _ = v < v
    with pytest.raises(RuntimeError):
        _ = v <= v
    with pytest.raises(RuntimeError):
        _ = v > v
    with pytest.raises(RuntimeError):
        _ = v >= v
