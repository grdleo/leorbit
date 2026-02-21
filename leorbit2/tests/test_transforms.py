import numpy as np
import pytest

from leorbit2.m import D, Matrix33, Scalar, Vector3, Vector3Array
from leorbit2.transforms import (
    TransformChain,
    TransformIdentify,
    TransformVector3Affine,
    TransformVector3Linear,
    TransformVector3RotationZ,
)


def test_transform_identity_do_undo_and_copy():
    t = TransformIdentify[Vector3[D.Length]]()
    v = Vector3[D.Length].from_components(
        Scalar[D.Length](1.0),
        Scalar[D.Length](2.0),
        Scalar[D.Length](3.0),
    )

    out = t.do(v)
    back = t.undo(out)
    t_copy = t.copy()

    np.testing.assert_allclose(out._values.flatten(), [1.0, 2.0, 3.0])
    np.testing.assert_allclose(back._values.flatten(), [1.0, 2.0, 3.0])
    assert type(t_copy) is type(t)


def test_transform_vector3_linear_do_and_undo_on_vector():
    scale2 = Matrix33[D.Dimless].from_elements(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    t = TransformVector3Linear[D.Length](scale2)
    v = Vector3[D.Length].from_components(
        Scalar[D.Length](1.0),
        Scalar[D.Length](-2.0),
        Scalar[D.Length](3.0),
    )

    transformed = t.do(v)
    restored = t.undo(transformed)

    np.testing.assert_allclose(transformed._values.flatten(), [2.0, -4.0, 6.0])
    np.testing.assert_allclose(restored._values.flatten(), [1.0, -2.0, 3.0])


def test_transform_vector3_linear_do_on_vector_array():
    scale3 = Matrix33[D.Dimless].from_elements(
        3.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 3.0,
    )
    t = TransformVector3Linear[D.Dimless](scale3)
    arr = Vector3Array[D.Dimless].from_components(
        np.array([1.0, -1.0]),
        np.array([2.0, -2.0]),
        np.array([3.0, -3.0]),
    )

    out = t.do(arr)

    np.testing.assert_allclose(
        out._values,
        np.array([[3.0, -3.0], [6.0, -6.0], [9.0, -9.0]], dtype=float),
    )


def test_transform_vector3_affine_do_and_undo():
    ident = Matrix33[D.Dimless].from_elements(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    translation = Vector3[D.Length].from_components(
        Scalar[D.Length](10.0),
        Scalar[D.Length](-5.0),
        Scalar[D.Length](2.0),
    )
    t = TransformVector3Affine[D.Length](ident, translation)
    v = Vector3[D.Length].from_components(
        Scalar[D.Length](1.0),
        Scalar[D.Length](2.0),
        Scalar[D.Length](3.0),
    )

    transformed = t.do(v)
    restored = t.undo(transformed)

    np.testing.assert_allclose(transformed._values.flatten(), [11.0, -3.0, 5.0])
    np.testing.assert_allclose(restored._values.flatten(), [1.0, 2.0, 3.0])


def test_transform_rotation_z_constructor_current_behavior():
    with pytest.raises(ValueError):
        _ = TransformVector3RotationZ[D.Dimless](Scalar[D.Angle](np.pi / 2))


def test_transform_reverse_swaps_do_and_undo():
    scale2 = Matrix33[D.Dimless].from_elements(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    t = TransformVector3Linear[D.Dimless](scale2)
    r = t.reverse()

    v = Vector3[D.Dimless].from_components(2.0, 4.0, 6.0)

    np.testing.assert_allclose(r.do(v)._values.flatten(), [1.0, 2.0, 3.0])
    np.testing.assert_allclose(r.undo(v)._values.flatten(), [4.0, 8.0, 12.0])


def test_transform_chain_do_and_undo():
    scale2 = Matrix33[D.Dimless].from_elements(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    linear = TransformVector3Linear[D.Length](scale2)

    ident = Matrix33[D.Dimless].from_elements(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    translation = Vector3[D.Length].from_components(
        Scalar[D.Length](5.0),
        Scalar[D.Length](0.0),
        Scalar[D.Length](-1.0),
    )
    affine = TransformVector3Affine[D.Length](ident, translation)

    chain = TransformChain[Vector3[D.Length], Vector3[D.Length]](linear, affine)
    v = Vector3[D.Length].from_components(
        Scalar[D.Length](1.0),
        Scalar[D.Length](2.0),
        Scalar[D.Length](3.0),
    )

    out = chain.do(v)
    back = chain.undo(out)

    np.testing.assert_allclose(out._values.flatten(), [7.0, 4.0, 5.0])
    np.testing.assert_allclose(back._values.flatten(), [1.0, 2.0, 3.0])


def test_transform_copy_current_behavior_for_linear_and_chain():
    mat = Matrix33[D.Dimless].from_elements(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    linear = TransformVector3Linear[D.Dimless](mat)
    linear_copy = linear.copy()

    linear.matrix._values[0, 0] = 10.0

    v = Vector3[D.Dimless].from_components(1.0, 0.0, 0.0)
    np.testing.assert_allclose(linear.do(v)._values.flatten(), [10.0, 0.0, 0.0])
    np.testing.assert_allclose(linear_copy.do(v)._values.flatten(), [10.0, 0.0, 0.0])

    chain = TransformChain[Vector3[D.Dimless], Vector3[D.Dimless]](linear)
    chain_copy = chain.copy()

    linear.matrix._values[1, 1] = 20.0

    y = Vector3[D.Dimless].from_components(0.0, 1.0, 0.0)
    np.testing.assert_allclose(chain.do(y)._values.flatten(), [0.0, 20.0, 0.0])
    np.testing.assert_allclose(chain_copy.do(y)._values.flatten(), [0.0, 20.0, 0.0])
