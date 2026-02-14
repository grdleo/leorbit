import numpy as np
import pytest

from leorbit2.mathematics.dimensions import D
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.tensor import ensure_same_dimensions, ensure_tensor


def test_ensure_tensor_wraps_numbers_as_dimless_tensor():
    t = ensure_tensor(3.5)
    assert t.check(D.Dimless)
    assert np.float64(t._values) == pytest.approx(3.5)


def test_ensure_tensor_keeps_tensor_instances():
    s = Scalar[D.Length].new(2.0)
    wrapped = ensure_tensor(s)
    assert wrapped is s


def test_ensure_same_dimensions_accepts_compatible_tensors():
    a = Scalar[D.Time].new(1.0)
    b = Scalar[D.Time].new(2.0)
    assert ensure_same_dimensions(a, b) is True


def test_ensure_same_dimensions_raises_for_incompatible_tensors():
    a = Scalar[D.Time].new(1.0)
    b = Scalar[D.Length].new(1.0)
    with pytest.raises(RuntimeError):
        ensure_same_dimensions(a, b)


def test_tensor_unary_and_copy_behavior():
    s = Scalar[D.Length].new(3.0)
    assert (+s).magnitude("m") == pytest.approx(3.0)
    assert (-s).magnitude("m") == pytest.approx(-3.0)

    c = s.copy()
    assert c is not s
    assert c.magnitude("m") == pytest.approx(s.magnitude("m"))
