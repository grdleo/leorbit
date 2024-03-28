import pytest
from leorbit.math.vector import Vec3
from leorbit.math import HALFPI, SQRT_3, SQRT_2
from leorbit.math import Q_


@pytest.fixture
def vectors():
    return (Vec3.zero(), Vec3.xaxis(), Vec3.one(), Vec3(1 / SQRT_3))


def test_operators():
    zero = Vec3(0.0, 0.0, 0.0)
    one = Vec3.one()
    x = Vec3(1.0, 0.0, 0.0)
    y = Vec3(0.0, 1.0, 0.0)

    assert zero.unit == "dimensionless"

    with pytest.raises(Exception):
        one += 1.0
    with pytest.raises(Exception):
        one - Q_(0.0, "m")

    onee = zero + one
    1 * onee * 1
    normed = one / SQRT_3
    assert abs(onee) == abs(one)
    assert abs(normed) == 1.0

    a: Vec3 = one - x
    assert a.x == 0.0
    assert a.y == 1.0
    assert a.z == 1.0
    assert abs(a) == SQRT_2
    assert a @ one == 2.0

    a += x
    assert a.x == 1.0
    a *= 2.0
    a /= 2.0

    assert abs(a) == SQRT_3

    assert one.dot(zero) == 0.0

    a = Vec3(1, 2, 3)
    b = Vec3(-1, 1, -1)

    assert a.dot(b) == Vec3.dot(a, b) == -2

    assert x.cross(x) == zero
    assert x & y == Vec3.zaxis()

def test_dimensions():
    d = Vec3(Q_("1m"), Q_("0m"), Q_("0m"))
    t = Q_(1.0, "s")
    v: Vec3 = d / t
    assert v.x.check("m/s")

    a = v / t
    assert a.x.check("m/s²")

    f = Vec3.one(unit="N")
    m = Q_(2.0, "kg")

    fa = a * m
    Q_(1, f.unit).to(fa.unit)

def test_rotation():
    x = Vec3.xaxis()
    yy = x.rt_z(HALFPI)
    assert yy.x == pytest.approx(0)
    assert yy.y == pytest.approx(1)
    assert yy.z == pytest.approx(0)
