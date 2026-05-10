from enum import Enum
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, cast

from leorbit.mathematics import Dimless, Length, Quantity, Tensor, cos, matrix33, vector3
from leorbit.transforms import Transform, TransformChain, TransformIdentity, TransformVector3Affine, TransformVector3Linear
from leorbit.time import Timestamp, TimeInterval

import numpy as np
import numpy.typing as npt

from leorbit.utils import unixepoch_to_j2000, j2000_to_stl0

if TYPE_CHECKING:
    from leorbit.coordinates import Coordinates

class AbsoluteFrame(Enum):
    """Absolute frame"""

    GCRF = "GCRF"
    ITRF = "ITRF"

FrameTransformFactory = Callable[[Timestamp | TimeInterval], Transform]

@lru_cache(4096)
def itrf2gcrf(epoch: Timestamp | TimeInterval) -> Transform:
    """Build the Earth rotation transform from ITRF to GCRF at ``epoch``.

    The returned transform applies equally to position and velocity vectors.
    """

    unixepoch: npt.NDArray[np.float64]
    if isinstance(epoch, Timestamp):
        unixepoch = np.asarray(epoch.unixepoch, dtype=np.float64)
    elif isinstance(epoch, TimeInterval):
        unixepoch = np.asarray(epoch.to_time_stamps().raw_data_array("second"), dtype=np.float64)
    else:
        raise ValueError(...)
    
    stl0 = np.asarray(j2000_to_stl0(unixepoch_to_j2000(unixepoch)), dtype=np.float64)
    cos_stl0 = np.cos(stl0)
    sin_stl0 = np.sin(stl0)
    one = np.ones_like(stl0)
    zero = np.zeros_like(stl0)

    if stl0.ndim == 0:
        rot_mat = matrix33(
            cos_stl0.item(), -sin_stl0.item(), zero.item(),
            sin_stl0.item(), cos_stl0.item(), zero.item(),
            zero.item(), zero.item(), one.item(),
        )
        return TransformVector3Linear(rot_mat)
    elif stl0.ndim == 1:
        rot_mat = np.stack(
            [
                np.stack([cos_stl0, -sin_stl0, zero], axis=0),
                np.stack([sin_stl0, cos_stl0,  zero], axis=0),
                np.stack([zero,     zero,      one ], axis=0),
            ], 
            axis=0
        )
        return TransformVector3Linear(Tensor(rot_mat, Dimless))
    else:
        raise ValueError("Unsupported sidereal angle shape")

    

ABS_FRAME_TRANSFORMS: dict[tuple[AbsoluteFrame, AbsoluteFrame], FrameTransformFactory] = {
    (AbsoluteFrame.ITRF, AbsoluteFrame.GCRF): itrf2gcrf,
    (AbsoluteFrame.GCRF, AbsoluteFrame.ITRF): lambda epoch: itrf2gcrf(epoch).reverse(),
}
"""Transformations between every absolute frames"""

@lru_cache
def absolute_frame_transform_factory(from_frame: AbsoluteFrame, to_frame: AbsoluteFrame) -> FrameTransformFactory:
    """Returns a function to be called with epoch as parameter. 
    That function returns a transformation that, applied to a vector `v` (whose coordinates are expressed in `from_frame`), 
    returns the same vector but whose coordinates are expressed in `to_frame`"""
    if from_frame == to_frame:
        return lambda epoch: TransformIdentity()

    frames = from_frame, to_frame
    factory = ABS_FRAME_TRANSFORMS.get(frames, None)
    if factory is None:
        raise NotImplementedError("No algorithm to compute composed transformations")
    
    return factory

### RELATIVE FRAMES

class RelativeFrame:
    """Frame defined by a transform relative to an absolute reference frame."""

    transform: Transform
    """The transformation that takes a vector expressed in the `reference_frame` and returns the same vector expressed in this `RelativeFrame`"""

    reference_frame: AbsoluteFrame
    """The absolute frame to which this frame is relative"""

    def __init__(self, reference_frame: AbsoluteFrame, transform: Transform):
        """Create a frame expressed from ``reference_frame`` via ``transform``."""
        self.reference_frame = reference_frame
        self.transform = transform

Frame = AbsoluteFrame | RelativeFrame

@lru_cache
def frame_transform_factory(from_frame: Frame, to_frame: Frame) -> FrameTransformFactory:
    """Finds and returns the transform factory for any frames.
    The idea is to transform to absolute frame reference of the source, 
    then transform to the absolute frame of the target, and then transform to source frame.
    Always works as long as transforms between absolute frames are defined properly"""

    first = TransformIdentity()
    last = TransformIdentity()

    abs_frame_from: AbsoluteFrame
    if isinstance(from_frame, AbsoluteFrame):
        abs_frame_from = from_frame
    elif isinstance(from_frame, RelativeFrame):
        first = from_frame.transform.reverse()
        abs_frame_from = from_frame.reference_frame
    else:
        raise RuntimeError("Unreachable?")

    abs_frame_to: AbsoluteFrame
    if isinstance(to_frame, AbsoluteFrame):
        abs_frame_to = to_frame
    elif isinstance(to_frame, RelativeFrame):
        last = to_frame.transform
        abs_frame_to = to_frame.reference_frame
    else:
        raise RuntimeError("Unreachable?")
    
    abs_transform = absolute_frame_transform_factory(abs_frame_from, abs_frame_to)

    if isinstance(first, TransformIdentity) and isinstance(last, TransformIdentity):
        return abs_transform
    
    def _factory(epoch: Timestamp | TimeInterval) -> TransformChain:
        return TransformChain(
            first,
            abs_transform(epoch),
            last,
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

    def __init__(self, location: "Coordinates"):
        """Create a local topocentric frame centered at ``location``."""
        itrf = location.get_pos(AbsoluteFrame.ITRF)

        z = itrf.vector3.normalized()
        north = vector3(0.0, 0.0, 1.0)
        ang = z.vector3.angle(north)

        half_turn = 180 * Quantity.degree
        quart_turn = 90 * Quantity.degree
        x: Tensor

        if ang % half_turn == 0: # FIXME
            raise ValueError("Cannot create `EarthLocalFrame` in Earth's poles!")
        elif ang == quart_turn: # FIXME
            x = north
        else:
            cos_ang = cos(ang)
            x = (north / cos_ang - z).vector3.normalized()
            if ang > quart_turn:
                x = -x
        
        y = x.vector3.cross(z)  # towards "east"

        mat = matrix33(
            x.vector3.x.scalar.value(), y.vector3.x.scalar.value(), z.vector3.x.scalar.value(),
            x.vector3.y.scalar.value(), y.vector3.y.scalar.value(), z.vector3.y.scalar.value(),
            x.vector3.z.scalar.value(), y.vector3.z.scalar.value(), z.vector3.z.scalar.value(),
        )

        # Transform : Local @ v -> ITRF @ v
        transform_local2itrf = TransformVector3Affine(
            mat,
            itrf,
        )

        # FIXME: please check that this is correct... And could optimize
        # Transform : ITRF @ v -> Local @ v
        transform_itrf2local = transform_local2itrf.reverse()

        super().__init__(
            AbsoluteFrame.ITRF, 
            transform_itrf2local,
        )

        self.location = location
    
    def __repr__(self) -> str:
        """Return a compact textual representation of the local frame location."""
        gps = self.location.gps()
        return f"<EarthLocalFrame at GPS location {gps.dms}>"