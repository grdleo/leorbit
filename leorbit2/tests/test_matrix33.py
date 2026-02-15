import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.matrix33 import Matrix33
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.vector3 import Vector3
from leorbit2.mathematics.vector3 import Vector3Array


def test_matrix33_new_from_numbers_and_scalars():
    m = Matrix33[D.Dimless].new(
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
    )
    np.testing.assert_allclose(m._values, np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))

    ms = Matrix33[D.Length].new(
        Scalar[D.Length].new(1), Scalar[D.Length].new(0), Scalar[D.Length].new(0),
        Scalar[D.Length].new(0), Scalar[D.Length].new(1), Scalar[D.Length].new(0),
        Scalar[D.Length].new(0), Scalar[D.Length].new(0), Scalar[D.Length].new(1),
    )
    assert ms.check(D.Length)


def test_matrix33_new_rejects_mixed_scalar_dimensions():
    with pytest.raises(RuntimeError):
        _ = Matrix33.new(
            Scalar[D.Length].new(1), Scalar[D.Length].new(0), Scalar[D.Length].new(0),
            Scalar[D.Length].new(0), Scalar[D.Time].new(1), Scalar[D.Length].new(0),
            Scalar[D.Length].new(0), Scalar[D.Length].new(0), Scalar[D.Length].new(1),
        )


def test_matrix33_add_sub_mul_div():
    a = Matrix33[D.Dimless].new(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    )
    b = Matrix33[D.Dimless].new(
        2, 2, 2,
        2, 2, 2,
        2, 2, 2,
    )

    np.testing.assert_allclose((a + b)._values, np.array([[3, 2, 2], [2, 3, 2], [2, 2, 3]]))
    np.testing.assert_allclose((b - a)._values, np.array([[1, 2, 2], [2, 1, 2], [2, 2, 1]]))

    two = Scalar[D.Dimless].new(2)
    np.testing.assert_allclose((a * two)._values, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
    np.testing.assert_allclose((a / two)._values, np.array([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]))


def test_matrix33_matmul_with_matrix_vector_and_vector_array():
    i = Matrix33[D.Dimless].new(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    )
    x2 = Matrix33[D.Dimless].new(
        2, 0, 0,
        0, 2, 0,
        0, 0, 2,
    )

    mm = i @ x2
    np.testing.assert_allclose(mm._values, x2._values)

    v = Vector3[D.Length].new(1, 2, 3)
    mv = x2 @ v
    np.testing.assert_allclose(mv._values.flatten(), [2, 4, 6])
    assert mv.dim_coords == v.dim_coords

    va = Vector3Array[D.Length](np.array([[1, 2], [0, 1], [3, 4]], dtype=float))
    mva = x2 @ va
    np.testing.assert_allclose(mva._values, np.array([[2, 4], [0, 2], [6, 8]], dtype=float))


def test_matrix33_inverse_computes_inverse():
    m = Matrix33[D.Dimless].new(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    )
    inv = m.inverse()
    # inverse of identity is identity
    np.testing.assert_allclose(inv._values, np.linalg.inv(m._values))
    assert inv.check(D.Dimless)
