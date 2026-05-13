import numpy as np
import pytest

from leorbit.mathematics import U, Tensor, matrix33, scalar, vector3
from leorbit.transforms import (
    TransformChain,
    TransformIdentity,
    TransformVector3Affine,
    TransformVector3Linear,
    TransformVector3RotationZ,
)

Angle = U.radian
Dimless = U.dimensionless
Length = U.meter


def test_transform_identity_do_undo_and_copy():
    t = TransformIdentity()
    v = vector3(1.0, 2.0, 3.0) * U.meter

    out = t.do(v)
    back = t.undo(out)
    t_copy = t.copy()

    np.testing.assert_allclose(out.raw_data_array("meter").reshape(3), [1.0, 2.0, 3.0])
    np.testing.assert_allclose(back.raw_data_array("meter").reshape(3), [1.0, 2.0, 3.0])
    assert type(t_copy) is type(t)


def test_transform_vector3_linear_do_and_undo_on_vector():
    scale2 = matrix33(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    t = TransformVector3Linear(scale2)
    v = vector3(1.0, -2.0, 3.0) * U.meter

    transformed = t.do(v)
    restored = t.undo(transformed)

    np.testing.assert_allclose(transformed.raw_data_array("meter").reshape(3), [2.0, -4.0, 6.0])
    np.testing.assert_allclose(restored.raw_data_array("meter").reshape(3), [1.0, -2.0, 3.0])


def test_transform_vector3_linear_do_on_vector_array():
    scale3 = matrix33(
        3.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 3.0,
    )
    t = TransformVector3Linear(scale3)
    arr = vector3(1.0, 2.0, 3.0) | vector3(-1.0, -2.0, -3.0)

    out = t.do(arr)

    np.testing.assert_allclose(
        out.raw_data_array(out.units),
        np.array([[3.0, -3.0], [6.0, -6.0], [9.0, -9.0]], dtype=float),
    )


def test_transform_vector3_affine_do_and_undo():
    ident = matrix33(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    translation = vector3(10.0, -5.0, 2.0) * U.meter
    t = TransformVector3Affine(ident, translation)
    v = vector3(1.0, 2.0, 3.0) * U.meter

    transformed = t.do(v)
    restored = t.undo(transformed)

    np.testing.assert_allclose(transformed.raw_data_array("meter").reshape(3), [11.0, -3.0, 5.0])
    np.testing.assert_allclose(restored.raw_data_array("meter").reshape(3), [1.0, 2.0, 3.0])


def test_transform_rotation_z_quarter_turn():
    t = TransformVector3RotationZ(scalar(np.pi / 2).with_units(U.radian))
    x = vector3(1.0, 0.0, 0.0)

    y = t.do(x)
    x_back = t.undo(y)

    np.testing.assert_allclose(y.raw_data_array(y.units).reshape(3), [0.0, 1.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(x_back.raw_data_array(x_back.units).reshape(3), [1.0, 0.0, 0.0], atol=1e-12)


def test_transform_reverse_swaps_do_and_undo():
    scale2 = matrix33(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    t = TransformVector3Linear(scale2)
    r = t.reverse()

    v = vector3(2.0, 4.0, 6.0)

    rv_do = r.do(v)
    rv_undo = r.undo(v)
    np.testing.assert_allclose(rv_do.raw_data_array(rv_do.units).reshape(3), [1.0, 2.0, 3.0])
    np.testing.assert_allclose(rv_undo.raw_data_array(rv_undo.units).reshape(3), [4.0, 8.0, 12.0])


def test_transform_chain_do_and_undo():
    scale2 = matrix33(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )
    linear = TransformVector3Linear(scale2)

    ident = matrix33(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    translation = vector3(5.0, 0.0, -1.0) * U.meter
    affine = TransformVector3Affine(ident, translation)

    chain = TransformChain(linear, affine)
    v = vector3(1.0, 2.0, 3.0) * U.meter

    out = chain.do(v)
    back = chain.undo(out)

    np.testing.assert_allclose(out.raw_data_array("meter").reshape(3), [7.0, 4.0, 5.0])
    np.testing.assert_allclose(back.raw_data_array("meter").reshape(3), [1.0, 2.0, 3.0])


def test_transform_copy_current_behavior_for_linear_and_chain():
    mat = matrix33(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    linear = TransformVector3Linear(mat)
    linear_copy = linear.copy()

    linear.matrix._data[0, 0, 0] = 10.0

    v = vector3(1.0, 0.0, 0.0)
    lv = linear.do(v)
    lv_copy = linear_copy.do(v)
    np.testing.assert_allclose(lv.raw_data_array(lv.units).reshape(3), [10.0, 0.0, 0.0])
    np.testing.assert_allclose(lv_copy.raw_data_array(lv_copy.units).reshape(3), [1.0, 0.0, 0.0])

    chain = TransformChain(linear)
    chain_copy = chain.copy()

    linear.matrix._data[1, 1, 0] = 20.0

    y = vector3(0.0, 1.0, 0.0)
    cy = chain.do(y)
    cy_copy = chain_copy.do(y)
    np.testing.assert_allclose(cy.raw_data_array(cy.units).reshape(3), [0.0, 20.0, 0.0])
    np.testing.assert_allclose(cy_copy.raw_data_array(cy_copy.units).reshape(3), [0.0, 1.0, 0.0])
