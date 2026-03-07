import numpy as np
import pytest
from typing import cast

from leorbit.frames import (
    AbsoluteFrame,
    RelativeFrame,
    absolute_frame_transform_factory,
    frame_transform_factory,
)
from leorbit.mathematics import Dimless, Length, Matrix33, Scalar, Vector3, Velocity
from leorbit.time import Timestamp
from leorbit.transforms import TransformVector3Affine
from leorbit.transforms import Transform


def _vec_len(x: float, y: float, z: float) -> Vector3[Length]:
    return Vector3[Length].from_components(
        Scalar[Length](x),
        Scalar[Length](y),
        Scalar[Length](z),
    )


def _vec_vel(x: float, y: float, z: float) -> Vector3[Velocity]:
    return Vector3[Velocity].from_components(
        Scalar[Velocity](x),
        Scalar[Velocity](y),
        Scalar[Velocity](z),
    )


@pytest.mark.parametrize(
    "pos, matrix, tr",
    [
        (
            _vec_len(0.0, 3.0, -3.0),
            Matrix33[Dimless].from_elements(
                -1.0, 0.0, 1.0,
                0.0, 1.0, 2.0,
                -1.0, -2.0, 0.0,
            ),
            _vec_len(1.0, 1.5, 2.0),
        ),
        (
            _vec_len(1.0, 2.0, 3.0),
            Matrix33[Dimless].from_elements(
                1.5, -0.5, -2.0,
                2.0, 1.0, -5.0,
                -1.0, 2.0, -6.0,
            ),
            _vec_len(0.0, 10.0, 1.0),
        ),
    ],
)
def test_relative_frame_roundtrip(pos: Vector3[Length], matrix: Matrix33[Dimless], tr: Vector3[Length]):
    rel = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine[Length](matrix, tr)),
    )
    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")

    to_rel = frame_transform_factory(AbsoluteFrame.ITRF, rel)(epoch)
    to_abs = frame_transform_factory(rel, AbsoluteFrame.ITRF)(epoch)

    p_rel = to_rel.do(pos)
    p_abs = to_abs.do(p_rel)

    np.testing.assert_allclose(p_abs._values.flatten(), pos._values.flatten())


def test_wrong_relative_frame_matrix_not_invertible():
    singular = Matrix33[Dimless].from_elements(
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0,
        0.0, 0.0, 0.0,
    )
    rel = RelativeFrame(
        AbsoluteFrame.GCRF,
        cast(
            Transform,
            TransformVector3Affine[Length](
                singular,
                _vec_len(0.0, 0.0, 0.0),
            ),
        ),
    )

    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")
    to_abs = frame_transform_factory(rel, AbsoluteFrame.GCRF)(epoch)

    with pytest.raises(np.linalg.LinAlgError):
        _ = to_abs.do(_vec_len(1.0, 2.0, 3.0))


def test_absolute_frame_transform_factory_roundtrip_position_and_velocity():
    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")
    itrf_to_gcrf = absolute_frame_transform_factory(AbsoluteFrame.ITRF, AbsoluteFrame.GCRF)(epoch)
    gcrf_to_itrf = absolute_frame_transform_factory(AbsoluteFrame.GCRF, AbsoluteFrame.ITRF)(epoch)

    p = _vec_len(10.0, -2.0, 7.5)
    v = _vec_vel(1.0, 2.0, -3.0)

    p2 = gcrf_to_itrf.do(itrf_to_gcrf.do(p))
    v2 = gcrf_to_itrf.do(itrf_to_gcrf.do(v))

    np.testing.assert_allclose(p2._values.flatten(), p._values.flatten(), atol=1e-9)
    np.testing.assert_allclose(v2._values.flatten(), v._values.flatten(), atol=1e-9)


def test_frame_transform_factory_between_relative_frames_chain():
    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")

    m1 = Matrix33[Dimless].from_elements(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    m2 = Matrix33[Dimless].from_elements(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )

    r1 = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine[Length](m1, _vec_len(1.0, 0.0, 0.0))),
    )
    r2 = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine[Length](m2, _vec_len(0.0, 2.0, 0.0))),
    )

    direct = frame_transform_factory(r1, r2)(epoch)
    back = frame_transform_factory(r2, r1)(epoch)

    p_r1 = _vec_len(3.0, 4.0, 5.0)
    p_r2 = direct.do(p_r1)
    p_r1_back = back.do(p_r2)

    np.testing.assert_allclose(p_r1_back._values.flatten(), p_r1._values.flatten(), atol=1e-9)
