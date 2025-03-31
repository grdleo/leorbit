from dataclasses import dataclass
from functools import cached_property
from pint import Quantity as Q_

from algorithms.utils import angle2dms
from coordinates.coordinates import Coordinates
from coordinates.representations import CoordinatesRepresentation
from frames.earth_local_frame import EarthLocalFrame
from mathematics.custom import normalize_angle, normalize_angle_symmetric
from physics.time import Time

@dataclass(frozen=True)
class Horizontal(CoordinatesRepresentation):
    """Horizontal coordinates use a celestial sphere centered on the observer. 
    Azimuth is measured eastward from the north point of the horizon.
    Altitude is the angle above the horizon.

    If distance if ommited, the coordinates represents a direction in the sky, 
    and therefore cannot be converted to a 'real' coordinate.
    
    (see https://en.wikipedia.org/wiki/Horizontal_coordinate_system)"""

    azimuth: Q_
    altitude: Q_
    distance: Q_ | None = None

    @cached_property
    def dms(self) -> str:
        azi = normalize_angle(self.azimuth)
        alt = normalize_angle_symmetric(self.altitude)

        return (
            f"Azimuth: {angle2dms(azi)}, "
            f"Altitude: {'-' if alt < 0 else ''}{angle2dms(alt)}"
        )
    
    def __repr__(self) -> str:
        return f"<Horizontal: {self.dms}>"
    
    def to_coordinates(self, frame: EarthLocalFrame, epoch: Time = None) -> "Coordinates":
        """If no `epoch` is provided, uses the time of this functions execution."""

        if self.distance is None:
            raise ValueError("Cannot convert `Horizontal` representation with unset `distance` to `Coordinates`.")
        
        return Coordinates.from_horizontal(
            azimuth=self.azimuth,
            altitude=self.altitude,
            distance=self.distance,
            frame=frame,
            epoch=Time.now() if epoch is None else epoch
        )