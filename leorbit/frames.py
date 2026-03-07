from enum import Enum
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, TypeVar, cast

from leorbit.mathematics import Angle, Dim, Dimless, Length, Quantity, Scalar, Tensor_M33, Tensor_V3, Vector3, Velocity, cos
from leorbit.transforms import Transform, TransformChain, TransformIdentify, TransformVector3Affine, TransformVector3Linear
from leorbit.time import Timestamp, TimeInterval

import numpy as np
import numpy.typing as npt

from leorbit.utils import unixepoch_to_j2000, j2000_to_stl0

if TYPE_CHECKING:
    from leorbit.coordinates import Coordinates

PosVec = Vector3[Length]
VelVec = Vector3[Velocity]

DynamicVec = Vector3[Length] | Vector3[Velocity]
SomeDynamicVec = TypeVar("SomeDynamicVec", bound=DynamicVec)
DynamicD = Length | Velocity
SomeDynamicD = TypeVar("SomeDynamicD", bound=Dim)

class AbsoluteFrame(Enum):
    """Absolute frame"""

    GCRF = "GCRF"
    ITRF = "ITRF"

FrameTransformFactory = Callable[[Timestamp | TimeInterval], Transform[SomeDynamicVec, SomeDynamicVec]]

@lru_cache(4096)
def itrf2gcrf(epoch: Timestamp | TimeInterval) -> Transform[SomeDynamicVec, SomeDynamicVec]:
    """NOTE: This rotation can transform any position or velocity"""

    unixepoch: npt.NDArray
    if isinstance(epoch, Timestamp):
        unixepoch = np.array(epoch._unixepoch)
    elif isinstance(epoch, TimeInterval):
        unixepoch = epoch.to_time_stamps()._values
    else:
        raise ValueError(...)
    
    stl0 = j2000_to_stl0(unixepoch_to_j2000(unixepoch))
    cos_stl0 = np.cos(stl0)
    sin_stl0 = np.sin(stl0)
    one = np.ones_like(stl0)
    zero = np.zeros_like(stl0)

    if stl0.ndim == 0:
        rot_mat = np.array(
            [
                [cos_stl0.item(), -sin_stl0.item(), zero.item()],
                [sin_stl0.item(), cos_stl0.item(),  zero.item()],
                [zero.item(),     zero.item(),      one.item() ],
            ]
        )
        return cast(
            Transform[SomeDynamicVec, SomeDynamicVec],
            TransformVector3Linear(Tensor_M33[Dimless](rot_mat))
        )
    elif stl0.ndim == 1:
        rot_mat = np.stack(
            [
                np.stack([cos_stl0, -sin_stl0, zero], axis=0),
                np.stack([sin_stl0, cos_stl0,  zero], axis=0),
                np.stack([zero,     zero,      one ], axis=0),
            ], 
            axis=0
        )
        return cast(
            Transform[SomeDynamicVec, SomeDynamicVec],
            TransformVector3Linear(Tensor_M33[Dimless](rot_mat))
        )
    else:
        raise ValueError("Unsupported sidereal angle shape")

    

ABS_FRAME_TRANSFORMS: dict[tuple[AbsoluteFrame, AbsoluteFrame], FrameTransformFactory] = {
    (AbsoluteFrame.ITRF, AbsoluteFrame.GCRF): itrf2gcrf
}
"""Transformations between every absolute frames"""

@lru_cache
def absolute_frame_transform_factory(from_frame: AbsoluteFrame, to_frame: AbsoluteFrame) -> FrameTransformFactory:
    """Returns a function to be called with epoch as parameter. 
    That function returns a transformation that, applied to a vector `v` (whose coordinates are expressed in `from_frame`), 
    returns the same vector but whose coordinates are expressed in `to_frame`"""
    if from_frame == to_frame:
        return lambda epoch: cast(Transform[DynamicVec, DynamicVec], TransformIdentify[DynamicVec]())

    frames = from_frame, to_frame
    factory = ABS_FRAME_TRANSFORMS.get(frames, None)
    if factory is not None:
        return factory
    
    reverse_frames = to_frame, from_frame
    reverse_factory = ABS_FRAME_TRANSFORMS.get(reverse_frames, None)
    if reverse_factory is not None:
        return lambda epoch: reverse_factory(epoch).reverse()
    
    raise NotImplementedError("No algorithm to compute composed transformations")

### RELATIVE FRAMES

_AbsPos = TypeVar("_AbsPos")
_RelPos = TypeVar("_RelPos")

class RelativeFrame:
    transform: Transform
    """The transformation that takes a vector expressed in the `reference_frame` and returns the same vector expressed in this `RelativeFrame`"""

    reference_frame: AbsoluteFrame
    """The absolute frame to which this frame is relative"""

    def __init__(self, reference_frame: AbsoluteFrame, transform: Transform[_AbsPos, _RelPos]):
        """A frame relative to a reference frame. """
        self.reference_frame = reference_frame
        self.transform = transform

Frame = AbsoluteFrame | RelativeFrame

@lru_cache
def frame_transform_factory(from_frame: Frame, to_frame: Frame) -> FrameTransformFactory:
    """Finds and returns the transform factory for any frames.
    The idea is to transform to absolute frame reference of the source, 
    then transform to the absolute frame of the target, and then transform to source frame.
    Always works as long as transforms between absolute frames are defined properly"""

    first: Transform[DynamicVec, DynamicVec]
    last: Transform[DynamicVec, DynamicVec]
    first = last = TransformIdentify[DynamicVec]()

    abs_frame_from: AbsoluteFrame
    if isinstance(from_frame, AbsoluteFrame):
        abs_frame_from = from_frame
    elif isinstance(from_frame, RelativeFrame):
        first = cast(Transform[DynamicVec, DynamicVec], from_frame.transform.reverse())
        abs_frame_from = from_frame.reference_frame
    else:
        raise RuntimeError("Unreachable?")

    abs_frame_to: AbsoluteFrame
    if isinstance(to_frame, AbsoluteFrame):
        abs_frame_to = to_frame
    elif isinstance(to_frame, RelativeFrame):
        last = cast(Transform[DynamicVec, DynamicVec], to_frame.transform)
        abs_frame_to = to_frame.reference_frame
    else:
        raise RuntimeError("Unreachable?")
    
    abs_transform = absolute_frame_transform_factory(abs_frame_from, abs_frame_to)

    if isinstance(first, TransformIdentify) and isinstance(last, TransformIdentify):
        return abs_transform
    
    def _factory(epoch: Timestamp | TimeInterval) -> TransformChain:
        return TransformChain(
            cast(Transform, first),
            cast(Transform, abs_transform(epoch)), 
            cast(Transform, last)
        )

    return _factory

########################################################################

class EarthLocalFrame(RelativeFrame):
    """Coordinates frame relative to a given location on Earth.
    
        - Origin: Given location
        - `z:` towards zenith (aka. towards the sky, perpendicular to ground)
        - `y:` towards "East"
        - `x × y = -z`
    """
    location: "Coordinates"
    transform: Transform[Tensor_V3[Length], Tensor_V3[Length]]

    def __init__(self, location: "Coordinates"):
        itrf = location.get_pos(AbsoluteFrame.ITRF)

        z = itrf.normalized()
        north = Vector3.Z
        ang = z.angle(north)

        half_turn = 180 * Quantity.degree
        quart_turn = 90 * Quantity.degree
        x: Vector3[Dimless]

        if ang % half_turn == 0: # FIXME
            raise ValueError("Cannot create `EarthLocalFrame` in Earth's poles!")
        elif ang == quart_turn: # FIXME
            x = north
        else:
            cos_ang = cos(ang)
            x = (north / cos_ang - z).normalized()
            if ang > quart_turn:
                x = -x
        
        y = x.cross(z) # towards "east"

        mat = Tensor_M33[Dimless].from_elements(
            x.x, y.x, z.x,
            x.y, y.y, z.y,
            x.z, y.z, z.z
        )

        # Transform : Local @ v -> ITRF @ v
        transform_local2itrf = TransformVector3Affine[Length](
            mat,
            itrf,
        )

        # FIXME: please check that this is correct... And could optimize
        # Transform : ITRF @ v -> Local @ v
        transform_itrf2local = transform_local2itrf.reverse()

        super().__init__(
            AbsoluteFrame.ITRF, 
            cast(Transform, transform_itrf2local)
        )

        self.location = location
    
    def __repr__(self) -> str:
        gps = self.location.gps()
        return f"<EarthLocalFrame at GPS location {gps.dms}>"