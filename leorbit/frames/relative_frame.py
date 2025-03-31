from functools import lru_cache
from frames import Frame
from frames.absolute_frame import AbsoluteFrame, TransformFactory, absolute_frame_transform_factory
from mathematics.transformation import ChainTransform, Identity, Transform


class RelativeFrame(Frame):
    def __init__(self, reference_frame: AbsoluteFrame, transform: Transform):
        """A frame relative to a reference frame. """
        self.reference_frame = reference_frame
        self.transform = transform

@lru_cache
def frame_transform_factory(from_frame: Frame, to_frame: Frame) -> TransformFactory:
    first = last = Identity()

    abs_frame_from: AbsoluteFrame = from_frame
    if isinstance(from_frame, RelativeFrame):
        first = from_frame.transform.reverse()
        abs_frame_from = from_frame.reference_frame

    abs_frame_to: AbsoluteFrame = to_frame
    if isinstance(to_frame, RelativeFrame):
        last = to_frame.transform
        abs_frame_to = to_frame.reference_frame
    
    abs_transform = absolute_frame_transform_factory(abs_frame_from, abs_frame_to)

    if isinstance(first, Identity) and isinstance(last, Identity):
        return abs_transform

    return lambda epoch: ChainTransform(first, abs_transform(epoch), last)