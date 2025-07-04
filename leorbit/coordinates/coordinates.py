from functools import lru_cache
from typing import NamedTuple, Self
from algorithms.utils import geocentric_radius_earth
from coordinates.representations import CoordinatesRepresentation
from coordinates.representations.gps import GPS
from coordinates.representations.horizontal import Horizontal
from frames.absolute_frame import AbsoluteFrame, Frame
from frames.earth_local_frame import EarthLocalFrame
from frames.relative_frame import frame_transform_factory
from physics.time import Time
from mathematics.vec3 import Vec3, UREG, POS_UNIT

from pint import Quantity as Q_

class PosVel(NamedTuple):
	"""A tuple to hold position and velocity.
	All values are in SI units (meters, meters/second).
	Can hold coordinates as NumPy arrays or floats.
	If NumPy arrays are used, all values must have the same length.
	"""
	pos: Vec3  # Position in meters
	vel: Vec3  # Velocity in meters/second


class Coordinates:
    def __init__(self, epoch: Time, frame: Frame, pos: Vec3, vel: Vec3 = None):
        self.epoch = epoch
        self.positions: dict[Frame, PosVel] = {frame: PosVel(pos, vel)}
        self.privileged_frame = frame
        self.vel_available = vel is not None
        self.name = None
        
        self._already_computed_repr: dict[type[CoordinatesRepresentation], CoordinatesRepresentation] = {}

    def _compute_new_frame(self, frame: Frame):
        if frame in self.positions:
            return
        
        p, v = self.positions[self.privileged_frame]
        transform = frame_transform_factory(self.privileged_frame, frame)(self.epoch)
        p = transform.apply(p)
        if v is not None:
            v = transform.apply(v)

        self.positions[frame] = PosVel(p, v)

        if isinstance(frame, AbsoluteFrame) and not isinstance(self.privileged_frame, AbsoluteFrame):
            self.privileged_frame = frame
    
    def get_pos(self, frame: Frame) -> Vec3:
        self._compute_new_frame(frame)
        return self.positions[frame].pos
    
    def get_vel(self, frame: Frame) -> Vec3 | None:
        if not self.vel_available:
            return None
        self._compute_new_frame(frame)
        return self.positions[frame].vel
    
    @property
    def pos(self) -> Vec3:
        return self.get_pos(self.privileged_frame)
    
    @property
    def vel(self) -> Vec3 | None:
        return self.get_vel(self.privileged_frame)
    
    ### GPS ###
    ### ### ###
    
    @staticmethod
    def from_gps(longitude: Q_, latitude: Q_, altitude: Q_ = 0 * POS_UNIT, epoch: Time = None) -> Self:
        assert longitude.check(UREG.radians)
        assert latitude.check(UREG.radians)
        assert altitude.check(UREG.meters)

        theta = longitude
        delta = latitude
        rho = geocentric_radius_earth(delta) + altitude
        pos = Vec3.from_spherical(theta, delta, rho)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, AbsoluteFrame.ITRF, pos)
    
    @lru_cache
    def gps(self) -> GPS:
        """Returns this coordinates as their GPS representation"""
        gps = self._already_computed_repr.get(GPS, None)
        if gps is not None:
            return gps
        
        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        lon = itrf_pos.theta
        lat = itrf_pos.delta
        alt = itrf_pos.rho - geocentric_radius_earth(lat)

        return GPS(
            longitude=lon, 
            latitude=lat, 
            altitude=alt
        )
    
    ### ### ###

    ### HORIZONTAL ###
    ### ########## ###

    @staticmethod
    def from_horizontal(azimuth: Q_, altitude: Q_, distance: Q_, frame: EarthLocalFrame, epoch: Time = None) -> Self:
        assert isinstance(frame, EarthLocalFrame)
        
        pos_local = Vec3.from_spherical(azimuth, altitude, distance)
        epoch = Time.now() if epoch is None else epoch

        return Coordinates(epoch, frame, pos_local)
    
    @lru_cache
    def horizontal(self, local_frame: EarthLocalFrame) -> Horizontal:
        """Returns the representation of this coordinates as 'horizontal', in the given frame"""
        assert local_frame.reference_frame == AbsoluteFrame.ITRF # FIXME

        itrf_pos = self.get_pos(AbsoluteFrame.ITRF)
        horizontal_pos = local_frame.transform.apply(itrf_pos)

        return Horizontal(
            azimuth=horizontal_pos.theta,
            altitude=horizontal_pos.delta,
            distance=itrf_pos.rho
        )
    
    ### ########## ###
    
    def __eq__(self, o: Self) -> bool:
        common_computed_frames = set(self.positions.keys()).intersection(o.positions.keys())

        if common_computed_frames:
            f: Frame
            for f in common_computed_frames: break # Get any element from set
            p, v = self.positions[f]
            po, vo = o.positions[f]

            return (p == po) and (v == vo)
        
        # We get a frame already computed in `self` and compute it for the other
        o._compute_new_frame(self.privileged_frame)
        # Now they have a common frame. We can call `__eq__` again
        return self.__eq__(o)
    
    def __neq__(self, o: Self) -> bool:
        return not self.__eq__(o)