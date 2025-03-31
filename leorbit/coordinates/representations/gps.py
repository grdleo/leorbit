from dataclasses import dataclass
from functools import cached_property, lru_cache

from algorithms.utils import angle2dms
from coordinates.coordinates import Coordinates
from coordinates.representations import CoordinatesRepresentation

from pint import Quantity as Q_

from frames.earth_local_frame import EarthLocalFrame
from mathematics.custom import normalize_angle_symmetric
from mathematics.units import POS_UNIT
from physics.time import Time


@dataclass(frozen=True)
class GPS(CoordinatesRepresentation):
    """Representation of a point on Earth (or around) using GPS standards.
    """
    longitude: Q_
    latitude: Q_
    altitude: Q_ = 0 * POS_UNIT
        
    @cached_property
    def dms(self) -> str:
        """Representation of the GPS coordinates in DSM notation (degrees, minutes, seconds)

        Example: `39° 17′ N, 76° 36′ O`"""
        lon = normalize_angle_symmetric(self.longitude)
        lat = normalize_angle_symmetric(self.latitude)

        return (
            angle2dms(lon) + ("E" if lon >= 0 else "O") 
            + ", "
            + angle2dms(lat) + ("N" if lat >= 0 else "S") 
        )
    
    def __repr__(self) -> str:
        return f"<GPS: {self.dms}>"
    
    @lru_cache
    def to_coordinates(self, epoch: Time = None) -> "Coordinates":
        """If no `epoch` is provided, uses the time of this functions execution."""

        coordinates = Coordinates.from_gps(
            longitude=self.longitude,
            latitude=self.latitude,
            altitude=self.altitude,
            epoch=Time.now() if epoch is None else epoch
        )
        coordinates._already_computer_repr[GPS] = self
        
        return coordinates
    
    @cached_property
    def earth_local_frame(self) -> EarthLocalFrame:
        return EarthLocalFrame(
            self.to_coordinates()
        )