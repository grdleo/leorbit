import numpy as np
import pytest
from typing import cast

from leorbit.frames import (
    AbsoluteFrame,
    RelativeFrame,
    absolute_frame_transform_factory,
    frame_transform_factory,
)
from leorbit.mathematics import Quantity, Tensor, matrix33, vector3
from leorbit.time import Timestamp
from leorbit.transforms import TransformVector3Affine
from leorbit.transforms import Transform


def _vec_len(x: float, y: float, z: float) -> Tensor:
    return vector3(x, y, z) * Quantity.meter


def _vec_vel(x: float, y: float, z: float) -> Tensor:
    return vector3(x, y, z) * Quantity.meter / Quantity.second


@pytest.mark.parametrize(
    "pos, matrix, tr",
    [
        (
            _vec_len(0.0, 3.0, -3.0),
            matrix33(
                -1.0, 0.0, 1.0,
                0.0, 1.0, 2.0,
                -1.0, -2.0, 0.0,
            ),
            _vec_len(1.0, 1.5, 2.0),
        ),
        (
            _vec_len(1.0, 2.0, 3.0),
            matrix33(
                1.5, -0.5, -2.0,
                2.0, 1.0, -5.0,
                -1.0, 2.0, -6.0,
            ),
            _vec_len(0.0, 10.0, 1.0),
        ),
    ],
)
def test_relative_frame_roundtrip(pos: Tensor, matrix: Tensor, tr: Tensor):
    rel = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine(matrix, tr)),
    )
    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")

    to_rel = frame_transform_factory(AbsoluteFrame.ITRF, rel)(epoch)
    to_abs = frame_transform_factory(rel, AbsoluteFrame.ITRF)(epoch)

    p_rel = to_rel.do(pos)
    p_abs = to_abs.do(p_rel)

    np.testing.assert_allclose(p_abs.raw_data_array("meter").reshape(3), pos.raw_data_array("meter").reshape(3))


def test_wrong_relative_frame_matrix_not_invertible():
    singular = matrix33(
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0,
        0.0, 0.0, 0.0,
    )
    rel = RelativeFrame(
        AbsoluteFrame.GCRF,
        cast(
            Transform,
            TransformVector3Affine(
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

    np.testing.assert_allclose(p2.raw_data_array("meter").reshape(3), p.raw_data_array("meter").reshape(3), atol=1e-9)
    np.testing.assert_allclose(v2.raw_data_array().reshape(3), v.raw_data_array().reshape(3), atol=1e-9)


def test_frame_transform_factory_between_relative_frames_chain():
    epoch = Timestamp.fromisoformat("2024-02-11T18:00:00")

    m1 = matrix33(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    m2 = matrix33(
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    )

    r1 = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine(m1, _vec_len(1.0, 0.0, 0.0))),
    )
    r2 = RelativeFrame(
        AbsoluteFrame.ITRF,
        cast(Transform, TransformVector3Affine(m2, _vec_len(0.0, 2.0, 0.0))),
    )

    direct = frame_transform_factory(r1, r2)(epoch)
    back = frame_transform_factory(r2, r1)(epoch)

    p_r1 = _vec_len(3.0, 4.0, 5.0)
    p_r2 = direct.do(p_r1)
    p_r1_back = back.do(p_r2)

    np.testing.assert_allclose(p_r1_back.raw_data_array("meter").reshape(3), p_r1.raw_data_array("meter").reshape(3), atol=1e-9)
