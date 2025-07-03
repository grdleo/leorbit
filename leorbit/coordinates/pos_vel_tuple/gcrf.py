
from functools import cache, cached_property, lru_cache
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from numpy.typing import NDArray


from leorbit.frames.absolute_frame import AbsoluteFrame
from leorbit.mathematics.units import Q_, UREG
from leorbit.mathematics.vec3 import Vec3
from leorbit.physics.time import Time

if TYPE_CHECKING:
	from leorbit.coordinates.coordinates import Coordinates


class PosVelGCRF(NamedTuple):
    """A tuple to hold position and velocity in GCRF coordinates.
	All values are in SI units (meters, meters/second).
    Can hold coordinates as NumPy arrays or floats.
	If NumPy arrays are used, all values must have the same length.
    """
    x: np.float64 | NDArray # [m]
    y: np.float64 | NDArray # [m]
    z: np.float64 | NDArray # [m]
    vx: np.float64 | NDArray # [m/s]
    vy: np.float64 | NDArray # [m/s]
    vz: np.float64 | NDArray # [m/s]

    @property
    def count(self) -> None | int:
        if self.array_values:
            return len(self.x)
        return None
    
    @property
    def array_values(self) -> bool:
        """Is `True` if the tuple is made from NumPy arrays, `False` if it is made from floats."""
        return isinstance(self.x, np.ndarray)
    
    def get_values_at(self, i: int) -> "PosVelGCRF":
        """If current tuple is made from NumPy arrays instead of floats,
        retrieve the values at given index and returns a float `PosVelGCRF`"""
        if not self.array_values:
            raise RuntimeError()
        
        PosVelGCRF(*(el[i] for el in self))
    
    def to_coordinates(self, epoch: Time) -> Coordinates:
        if self.array_values:
            raise ValueError()
        raise Coordinates(epoch, AbsoluteFrame.GCRF, self.pos, self.vel)
    
    @cached_property
    def pos(self) -> Vec3:
        assert not self.array_values
        return Vec3(self.x, self.y, self.z, UREG.meter)
    
    @cached_property
    def vel(self) -> Vec3:
        assert not self.array_values
        return Vec3(self.vx, self.vy, self.vz, UREG.meter / UREG.second)