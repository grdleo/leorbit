from leorbit.math.coordinate import AbsoluteFrame, RelativeFrame
from leorbit.math.vector import Vec3

from leorbit.math import Q_

import pytest

@pytest.mark.parametrize(
        "pos, nx, ny, nz, tr",
        [
            (
                Vec3(Q_("0m"), Q_("3m"), Q_("-3m")),
                Vec3(-1, 0, -1),
                Vec3(0, 1, -2),
                Vec3(1, 2, 0),
                Vec3(Q_("1m"), Q_("1.5m"), Q_("2m"))
            ),
            (
                Vec3(Q_("1m"), Q_("2m"), Q_("3m")),
                Vec3(1.5, 2, -1),
                Vec3(-.5, 1, 2),
                Vec3(-2, -5, -6),
                Vec3(Q_("0m"), Q_("10m"), Q_("1m"))
            ),
        ]
)
def test_change_coordinates(pos: Vec3, nx: Vec3, ny: Vec3, nz: Vec3, tr: Vec3):
    frame = RelativeFrame(AbsoluteFrame.ITRF, nx, ny, nz, tr)
    pp = frame.new_pos(pos)
    ppp = frame.old_pos(pp)

    assert pos.x.to("m").m == pytest.approx(ppp.x.to("m").m)
    assert pos.y.to("m").m == pytest.approx(ppp.y.to("m").m)
    assert pos.z.to("m").m == pytest.approx(ppp.z.to("m").m)

def test_wrong_coordinates_system():
    # Given axises do not generate the whole 3D space: matrix not inversible
    with pytest.raises(ZeroDivisionError):
        RelativeFrame(
            AbsoluteFrame.GCRF,
            Vec3(1, 0, 0),
            Vec3(0, 1, 0),
            Vec3(1, 1, 0),
            Vec3.zero("m"),
        )