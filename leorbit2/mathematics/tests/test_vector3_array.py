import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.scalar_array import ScalarArray
from leorbit2.mathematics.vector3 import Vector3
from leorbit2.mathematics.vector3_array import Vector3Array


def test_vector3_array_builders_and_indexing():
    v1 = Vector3[D.Length].new(1, 0, 0)
    v2 = Vector3[D.Length].new(0, 1, 0)
    arr = Vector3Array.new_from_vectors(v1, v2)

    assert arr.size == 2
    np.testing.assert_allclose(arr[0]._values.flatten(), [1, 0, 0])
    np.testing.assert_allclose(arr[1]._values.flatten(), [0, 1, 0])
    with pytest.raises(KeyError):
        _ = arr[2]

    repeated = Vector3Array.new_from_single_vector(v1, 3)
    assert repeated.size == 3
    np.testing.assert_allclose(repeated._values, np.array([[1, 1, 1], [0, 0, 0], [0, 0, 0]]))


def test_vector3_array_builders_reject_invalid_inputs():
    with pytest.raises(ValueError):
        _ = Vector3Array.new_from_vectors()

    with pytest.raises(ValueError):
        _ = Vector3Array.new_from_vectors(Vector3[D.Length].new(1, 0, 0), Vector3[D.Time].new(1, 0, 0))


def test_vector3_array_properties_and_ops():
    arr = Vector3Array[D.Dimless](np.array([[3, 0], [0, 4], [4, 0]], dtype=float))

    np.testing.assert_allclose(arr.length.magnitude(), [5, 4])

    n = arr.normalized()
    np.testing.assert_allclose(n._values[:, 0], [0.6, 0.0, 0.8])
    np.testing.assert_allclose(n._values[:, 1], [0.0, 1.0, 0.0])

    plus = arr + ScalarArray[D.Dimless].new([1, 2])
    np.testing.assert_allclose(plus._values, np.array([[4, 2], [1, 6], [5, 2]], dtype=float))


def test_vector3_array_dot_cross_and_matmul():
    arr = Vector3Array[D.Dimless](np.array([[1, 0], [0, 1], [0, 0]], dtype=float))
    x = Vector3.X

    dot = arr @ x
    np.testing.assert_allclose(dot.magnitude(), [1, 0])

    cross = arr.cross(x)
    np.testing.assert_allclose(cross._values, np.array([[0, 0], [0, 0], [0, -1]], dtype=float))


def test_vector3_array_comparisons_and_invalid_comparators():
    a = Vector3Array[D.Dimless](np.array([[1, 2], [0, 1], [0, 1]], dtype=float))
    b = Vector3Array[D.Dimless](np.array([[1, 3], [0, 1], [0, 0]], dtype=float))

    np.testing.assert_array_equal(a == b, [True, False])
    np.testing.assert_array_equal(a != b, [False, True])

    with pytest.raises(RuntimeError):
        _ = a < b
    with pytest.raises(RuntimeError):
        _ = a <= b
    with pytest.raises(RuntimeError):
        _ = a > b
    with pytest.raises(RuntimeError):
        _ = a >= b
