from enum import Enum
from functools import lru_cache
from typing import Callable
from mathematics.transformation import Identity, Transform, RotationZ, ChainTransform
from physics.time_custom import Time
from mathematics.vec3 import Vec3
from physics.time_custom import Time
from frames import Frame

class AbsoluteFrame(Frame, Enum):
    GCRF = "GCRF"
    ITRF = "ITRF"

@lru_cache(2**16)
def itrf2gcrf(epoch: Time) -> Transform:
    return RotationZ(epoch.stl0)

ABS_FRAME_TRANSFORMS = {
    (AbsoluteFrame.GCRF, AbsoluteFrame.ITRF): itrf2gcrf
}

TransformFactory = Callable[[Time, ], Transform]
@lru_cache
def absolute_frame_transform_factory(from_frame: AbsoluteFrame, to_frame: AbsoluteFrame) -> TransformFactory:
    """Returns a function to be called with epoch as parameter. 
    That function returns a transformation that, applied to a vector `v` (whose coordinates are expressed in `from_frame`), 
    returns the same vector but whose coordinates are expressed in `to_frame`"""
    frames = from_frame, to_frame
    factory = ABS_FRAME_TRANSFORMS.get(frames, None)
    if factory is not None:
        return factory
    
    reverse_frames = to_frame, from_frame
    reverse_factory = ABS_FRAME_TRANSFORMS.get(reverse_frames, None)
    if reverse_frames is not None:
        return lambda epoch: reverse_factory(epoch).reverse()
    
    raise NotImplementedError("No algorithm to compute composed transformations")