from enum import Enum
from functools import lru_cache
from typing import ParamSpec, Callable, TypeVar, cast

from leorbit2.mathematics import LengthDim, SomeDim, Vector3, VelocityDim, Transform, TransformVector3RotationZ
from leorbit2.time import Time

PosVec = Vector3[LengthDim]
VelVec = Vector3[VelocityDim]

DynamicVec = PosVec | VelVec
SomeDynamicVec = TypeVar("SomeDynamicVec", bound=DynamicVec)
DynamicDim = LengthDim | VelocityDim
SomeDynamicDim = TypeVar("SomeDynamicDim", bound=DynamicDim)

class Frame:
    """Base class for frames"""

class AbsoluteFrame(Frame, Enum):
    """Absolute frame"""

    GCRF = "GCRF"
    ITRF = "ITRF"

AbsoluteFrameTransformFactory = Callable[[Time], Transform[DynamicVec, DynamicVec]]

@lru_cache(4096)
def itrf2gcrf(epoch: Time) -> Transform[SomeDynamicVec, SomeDynamicVec]:
    """NOTE: This rotation can transform any position or velocity"""

    t = TransformVector3RotationZ[DynamicDim](epoch.stl0)
    return cast(Transform[SomeDynamicVec, SomeDynamicVec], t)

ABS_FRAME_TRANSFORMS: dict[tuple[AbsoluteFrame, AbsoluteFrame], AbsoluteFrameTransformFactory] = {
    (AbsoluteFrame.GCRF, AbsoluteFrame.ITRF): itrf2gcrf
}
"""Transformations between every absolute frames"""

@lru_cache
def absolute_frame_transform_factory(from_frame: AbsoluteFrame, to_frame: AbsoluteFrame) -> AbsoluteFrameTransformFactory:
    """Returns a function to be called with epoch as parameter. 
    That function returns a transformation that, applied to a vector `v` (whose coordinates are expressed in `from_frame`), 
    returns the same vector but whose coordinates are expressed in `to_frame`"""
    frames = from_frame, to_frame
    factory = ABS_FRAME_TRANSFORMS.get(frames, None)
    if factory is not None:
        return factory
    
    reverse_frames = to_frame, from_frame
    reverse_factory = ABS_FRAME_TRANSFORMS.get(reverse_frames, None)
    if reverse_factory is not None:
        return lambda epoch: reverse_factory(epoch).reverse()
    
    raise NotImplementedError("No algorithm to compute composed transformations")