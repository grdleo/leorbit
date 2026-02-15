from enum import Enum
from functools import lru_cache
from typing import TYPE_CHECKING, ParamSpec, Callable, TypeVar, cast

from leorbit2.m import D, Vector3
from leorbit2.transforms import Transform, TransformChain, TransformIdentify, TransformVector3Affine, TransformVector3RotationZ
from leorbit2.time import Time

if TYPE_CHECKING:
    from leorbit2.coordinates import Coordinates

PosVec = Vector3[D.Length]
VelVec = Vector3[D.Velocity]

DynamicVec = Vector3[D.Length] | Vector3[D.Velocity]
SomeDynamicVec = TypeVar("SomeDynamicVec", bound=DynamicVec)
DynamicD = D.Length | D.Velocity
SomeDynamicD = TypeVar("SomeDynamicD", bound=DynamicD)

class AbsoluteFrame(Enum):
    """Absolute frame"""

    GCRF = "GCRF"
    ITRF = "ITRF"

FrameTransformFactory = Callable[[Time], Transform[SomeDynamicVec, SomeDynamicVec]]

@lru_cache(4096)
def itrf2gcrf(epoch: Time) -> Transform[SomeDynamicVec, SomeDynamicVec]:
    """NOTE: This rotation can transform any position or velocity"""

    t = TransformVector3RotationZ[DynamicD](epoch.stl0)
    return cast(Transform[SomeDynamicVec, SomeDynamicVec], t)

ABS_FRAME_TRANSFORMS: dict[tuple[AbsoluteFrame, AbsoluteFrame], FrameTransformFactory] = {
    (AbsoluteFrame.GCRF, AbsoluteFrame.ITRF): itrf2gcrf
}
"""Transformations between every absolute frames"""

@lru_cache
def absolute_frame_transform_factory(from_frame: AbsoluteFrame, to_frame: AbsoluteFrame) -> FrameTransformFactory:
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

### RELATIVE FRAMES

_AbsPos = TypeVar("_AbsPos", bound=PosVec)
_RelPos = TypeVar("_RelPos", bound=PosVec)

class RelativeFrame:
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
    
    def _factory(epoch: Time) -> TransformChain:
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
    transform: TransformVector3Affine[DynamicD]

    def __init__(self, location: "Coordinates"):
        itrf = location.get_pos(AbsoluteFrame.ITRF)
        
        z = itrf.normalized() # towards zenith
        north = Vector3.Z
        ang = z.angle(north)

        half_turn = 180 * Quantity.deg
        quart_turn = 90 * Quantity.deg

        if ang % half_turn == 0: # FIXME
            raise ValueError("Cannot create `EarthLocalFrame` in Earth's poles!")
        elif ang == quart_turn: # FIXME
            x = north.copy()
        else:
            x = (north / cos(ang) - z).normalized()
            if ang > quart_turn:
                x = -x
        
        y = x.cross(z) # towards "east"

        mat = Matrix33[D.Dimless].new(
            x.x, y.x, z.x,
            x.y, y.y, z.y,
            x.z, y.z, z.z
        )

        transform = TransformVector3Affine(mat, itrf)

        super().__init__(AbsoluteFrame.ITRF, transform)

        self.location = location
    
    def __repr__(self) -> str:
        gps = self.location.gps()
        return f"<EarthLocalFrame at GPS location {gps.dms}>"