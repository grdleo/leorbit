import numpy as np
import pytest

from leorbit2.m import (
    D,
    Matrix33,
    Quantity,
    Scalar,
    ScalarArray,
    Tensor_V3,
    Vector3,
    Vector3Array,
    abs,
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
    normalize_angle,
    normalize_angle_symmetric,
    sin,
    sqrt,
    square,
    tan,
)


def test_quantity_units_and_getter():
    assert Quantity.meter.check(D.Length)
    assert Quantity.second.check(D.Time)
    assert Quantity.rad.check(D.Angle)
    assert Quantity.get("kilo_meter").magnitude("meter") == pytest.approx(1000)


def test_scalar_arithmetic_dimensions_and_cast():
    a = Scalar[D.Length](2)
    b = Scalar[D.Length](500)

    assert (a + b).magnitude("meter") == pytest.approx(502)
    assert (a - b).magnitude("meter") == pytest.approx(-498)

    speed = a / Scalar[D.Time](2)
    assert speed.dim_coords == (D.Length._d / D.Time._d)

    area = a * Scalar[D.Length](3)
    assert area.dim_coords == (D.Length._d * D.Length._d)

    same = a.cast(D.Length)
    assert same is a
    with pytest.raises(RuntimeError):
        a.cast(D.Time)


def test_scalar_modulo_and_comparisons():
    x = Scalar[D.Dimless](7)
    y = Scalar[D.Dimless](3)

    assert (x % y).magnitude() == pytest.approx(1)
    assert (Scalar[D.Dimless](7) % y).magnitude() == pytest.approx(1)

    t1 = Scalar[D.Time](2)
    t2 = Scalar[D.Time](3)
    assert t1 < t2
    assert t1 <= t2
    assert t2 > t1
    assert t2 >= t1


def test_ensure_tensor_and_dimension_helpers():
    t = ensure_tensor(3.5)
    assert t.check(D.Dimless)
    assert float(np.asarray(t._values)) == pytest.approx(3.5)

    s1 = Scalar[D.Time](1)
    s2 = Scalar[D.Time](2)
    assert ensure_same_dimensions(s1, s2) is True

    with pytest.raises(RuntimeError):
        ensure_same_dimensions(Quantity.second, Quantity.meter)


def test_vector3_from_components_and_basic_vector_ops():
    v = Vector3[D.Length].from_components(1.0, 2.0, 3.0)
    w = Vector3[D.Length].from_components(3.0, 2.0, 1.0)

    vsum = v + w
    np.testing.assert_allclose(vsum._values.flatten(), [4, 4, 4])

    vscaled = v * Scalar[D.Dimless](2)
    np.testing.assert_allclose(vscaled._values.flatten(), [2, 4, 6])

    cross = v.cross(w)
    np.testing.assert_allclose(cross._values.flatten(), [-4, 8, -4])


def test_vector3_from_components_rejects_mixed_dims():
    with pytest.raises(RuntimeError):
        Tensor_V3[D.Length].from_components(Quantity.meter, Quantity.second, Quantity.meter)  # type: ignore[arg-type]


def test_vector3array_from_components_and_matrix_products():
    arr = Vector3Array[D.Dimless].from_components(
        np.array([1.0, 0.0]),
        np.array([0.0, 1.0]),
        np.array([0.0, 0.0]),
    )
    assert arr.size == 2

    m = Matrix33[D.Dimless].from_elements(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    out = m @ arr
    np.testing.assert_allclose(out._values, np.array([[2, 0], [0, 2], [0, 0]], dtype=float))


def test_matrix33_arithmetic_inverse_and_products():
    i = Matrix33[D.Dimless].from_elements(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    )
    x2 = Matrix33[D.Dimless].from_elements(
        2, 0, 0,
        0, 2, 0,
        0, 0, 2,
    )

    np.testing.assert_allclose((i + x2)._values, np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=float))
    np.testing.assert_allclose((x2 - i)._values, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float))

    mm = i @ x2
    np.testing.assert_allclose(mm._values, x2._values)

    inv = x2.inverse()
    np.testing.assert_allclose(inv._values, np.linalg.inv(x2._values))

    v = Vector3[D.Length].from_components(1.0, 2.0, 3.0)
    mv = x2 @ v
    np.testing.assert_allclose(mv._values.flatten(), [2, 4, 6])


def test_square_sqrt_and_trig_functions():
    s = Scalar[D.Length](3)
    sq = square(s)
    assert sq.dim_coords == (D.Length._d ** 2)
    assert sq.magnitude() == pytest.approx(9)

    root = sqrt(sq)
    assert root.magnitude("meter") == pytest.approx(3)

    ang = Scalar[D.Angle](np.pi / 2)
    assert sin(ang).magnitude() == pytest.approx(1)
    assert cos(ang).magnitude() == pytest.approx(0, abs=1e-12)
    assert tan(Scalar[D.Angle](0)).magnitude() == pytest.approx(0)

    u = Scalar[D.Dimless](0.5)
    assert asin(u).check(D.Angle)
    assert acos(u).check(D.Angle)
    assert atan(u).check(D.Angle)


def test_atan2_and_angle_normalization_helpers():
    y = Scalar[D.Length](1)
    x = Scalar[D.Length](1)
    a = atan2(y, x)
    assert a.check(D.Angle)
    assert a.magnitude("rad") == pytest.approx(np.pi / 4)

    n = normalize_angle(Scalar[D.Angle](5 * np.pi))
    assert n.magnitude("rad") == pytest.approx(np.pi)

    with pytest.raises(TypeError):
        _ = normalize_angle_symmetric(Scalar[D.Angle](3 * np.pi / 2))


def test_tensor_v3_direct_cross_operation():
    x = Tensor_V3[D.Dimless].from_components(1.0, 0.0, 0.0)
    y = Tensor_V3[D.Dimless].from_components(0.0, 1.0, 0.0)
    cross = x.cross(y)
    np.testing.assert_allclose(cross._values.flatten(), [0.0, 0.0, 1.0])


def test_abs_helper_on_scalar_vector_and_matrix():
    s = Scalar[D.Length](-12)
    assert abs(s).magnitude("meter") == pytest.approx(12)

    v = Tensor_V3[D.Length].from_components(-1.0, 2.0, -3.0)
    av = abs(v)
    np.testing.assert_allclose(av._values.flatten(), [1.0, 2.0, 3.0])
    assert av.dim_coords == D.Length._d

    m = Matrix33[D.Time].from_elements(
        -1, 2, -3,
        4, -5, 6,
        -7, 8, -9,
    )
    am = abs(m)
    np.testing.assert_allclose(
        am._values,
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float),
    )
    assert am.dim_coords == D.Time._d


def test_cube_and_cbrt_roundtrip_and_dimensions():
    length = Scalar[D.Length](4)
    volume = cube(length)
    assert volume.dim_coords == (D.Length._d ** 3)
    assert volume.magnitude() == pytest.approx(64)

    restored = cbrt(volume)
    assert restored.dim_coords == D.Length._d
    assert restored.magnitude("meter") == pytest.approx(4)


def test_interpolate_scalar_vector_and_matrix():
    a = Scalar[D.Length](10)
    b = Scalar[D.Length](14)
    mid = interpolate(a, b, 0.25)
    assert mid.magnitude("meter") == pytest.approx(11)

    v1 = Vector3[D.Dimless].from_components(0.0, 0.0, 0.0)
    v2 = Vector3[D.Dimless].from_components(2.0, 4.0, 6.0)
    vm = interpolate(v1, v2, 0.5)
    np.testing.assert_allclose(vm._values.flatten(), [1.0, 2.0, 3.0])

    m1 = Matrix33[D.Dimless].from_elements(
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
    )
    m2 = Matrix33[D.Dimless].from_elements(
        2, 2, 2,
        2, 2, 2,
        2, 2, 2,
    )
    mm = interpolate(m1, m2, 0.5)
    np.testing.assert_allclose(mm._values, np.ones((3, 3)))


def test_dot_for_single_vector_returns_scalar():
    v = Tensor_V3[D.Dimless].from_components(3.0, 4.0, 12.0)
    dot = v.dot(v)
    assert dot.check(D.Dimless)
    assert dot.magnitude() == pytest.approx(169)


def test_normalize_angle_wrap_for_negative_scalar():
    a = Scalar[D.Angle](-np.pi / 2)
    wrapped = normalize_angle(a)
    assert wrapped.magnitude("rad") == pytest.approx(3 * np.pi / 2)
