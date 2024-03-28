"""Frames definition, conversions, coordinates handling, ..."""

from collections import namedtuple
from enum import Enum
from functools import lru_cache
from leorbit.algorithms.utils import geocentric_radius_earth
from leorbit.math.vector import Matrix33AsTuple, Vec3, inverse_mat3x3
from leorbit.time import Time
from leorbit.math import HALF_REV, Q_, RIGHT_ANGLE
from math import cos, sin
from dataclasses import dataclass, field
from typing import overload

__pdoc__ = dict()

PosVelGCRF = namedtuple("PosVelGCRF", ["x", "y", "z", "vx", "vy", "vz"])
__pdoc__["PosVelGCRF"] = "Tuple corresponding to a satellite's coordinates given in GCRF."
__pdoc__["PosVelGCRF.x"] = "`float` Coordinate's X projection of position (GCRF). Units: `meters`"
__pdoc__["PosVelGCRF.y"] = "`float` Coordinate's Y projection of position (GCRF). Units: `meters`"
__pdoc__["PosVelGCRF.z"] = "`float` Coordinate's Z projection of position (GCRF). Units: `meters`"
__pdoc__["PosVelGCRF.vx"] = "`float` Coordinate's X projection of velocity (GCRF). Units: `meters/second`"
__pdoc__["PosVelGCRF.vy"] = "`float` Coordinate's Y projection of velocity (GCRF). Units: `meters/second`"
__pdoc__["PosVelGCRF.vz"] = "`float` Coordinate's Z projection of velocity (GCRF). Units: `meters/second`"

class AbsoluteFrame(Enum):
    """Celestial coordinate frames well defined, considered as absolute.
    Conversions from any `AbsoluteFrame` to one another has to be defined and implemented.
    
    GCRF (Geocentric Celestial Reference Frame, also Equatorial coordinate system):

    - Origin: Earth's mass center
    - `x:` towards the Vernal point
    - `z:` towards Earth's true north
    - `x × y = z`
    
    ITRF (International Terrestrial Reference Frame):
        
    - Origin: Earth's mass center
    - `x:` towards the origin of GPS coordinates
    - `z:` towards Earth's true north
    - `x × y = z`
    """
    GCRF = "GCRF"
    ITRF = "ITRF"

    def __repr__(self) -> str:
        return f"<AbsoluteFrame: {self.value}>"

    @staticmethod
    def itrf2gcrf(itrf: Vec3, epoch: Time) -> Vec3:
        return itrf.rt_z(epoch.stl0)
    
    @staticmethod
    def gcrf2itrf(gcrf: Vec3, epoch: Time) -> Vec3:
        return gcrf.rt_z(-epoch.stl0)

class RelativeFrame:
    """New frame constructed relative to another `AbsoluteFrame`.
    One can convert a vector whose coordinates are given in the same `AbsoluteFrame`"""
    def __init__(self, frame: AbsoluteFrame, xaxis: Vec3, yaxis: Vec3, zaxis: Vec3, translation: Vec3):
        if not isinstance(frame, AbsoluteFrame):
            raise ValueError()
        
        try:
            assert xaxis.x.check("1")
            assert yaxis.x.check("1")
            assert zaxis.x.check("1")
            translation.x.check("m")
        except:
            raise ValueError()
        
        self.absolute_frame = frame
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.zaxis = zaxis
        self.translation = translation

        self._matrix_old_base = Matrix33AsTuple(
            self.xaxis._x.m, self.yaxis._x.m, self.zaxis._x.m, 
            self.xaxis._y.m, self.yaxis._y.m, self.zaxis._y.m, 
            self.xaxis._z.m, self.yaxis._z.m, self.zaxis._z.m, 
        )

        self._matrix_new_base = inverse_mat3x3(self._matrix_old_base)
    
    def new_pos(self, of: Vec3) -> Vec3:
        """Arguments:
        -------------
        `of: Vec3` Position vector to be transformed, given in the **same frame as this `RelativeFrame` was defined.**
        
        Returns:
        --------
        `Vec3` Transformed vector, given in this `RelativeFrame.`"""
        return (of - self.translation).apply_mat_list(self._matrix_new_base)
    
    def new_vel(self, of: Vec3) -> Vec3:
        """Arguments:
        -------------
        `of: Vec3` Velocity vector to be transformed, given in the **same frame as this `RelativeFrame` was defined.**
        
        Returns:
        --------
        `Vec3` Transformed vector, given in this `RelativeFrame.`"""
        return of.apply_mat_list(self._matrix_new_base)
    
    def old_pos(self, of: Vec3) -> Vec3:
        """Arguments:
        -------------
        `of: Vec3` Position vector to be transformed, given in **this `RelativeFrame.`**
        
        Returns:
        --------
        `Vec3` Transformed vector, given in the **same frame as this `RelativeFrame` was defined.**"""
        return of.apply_mat_list(self._matrix_old_base) + self.translation
    
    def old_vel(self, of: Vec3) -> Vec3:
        """Arguments:
        -------------
        `of: Vec3` Velocity vector to be transformed, given in **this `RelativeFrame.`**
        
        Returns:
        --------
        `Vec3` Transformed vector, given in the **same frame as this `RelativeFrame` was defined.**"""
        return of.apply_mat_list(self._matrix_old_base)

AnyFrame = AbsoluteFrame | RelativeFrame

@dataclass(frozen=True)
class HorizontalCoordinates:
    """Coordinates of an object given in a specific `EarthLocalFrame`.

    Both angle are given in degrees.
    
    [See on Wikipedia](https://en.wikipedia.org/wiki/Horizontal_coordinate_system)"""
    azimuth_deg: float
    altitude_deg: float

class EarthLocalFrame(RelativeFrame):
    """Coordinates frame relative to a given location on Earth.
    
        - Origin: Given location
        - `z:` towards zenith (aka. towards the sky, perpendicular to ground)
        - `y:` towards "East"
        - `x × y = -z`
    """
    def __init__(self, loc: "Coordinates"):
        itrf: Vec3 = loc._get_existant_pos(AbsoluteFrame.ITRF)
        if itrf is None:
            raise RuntimeError()
        
        z = itrf.normalize() # towards zenith
        north = Vec3.zaxis()
        ang = Vec3.angle(z, north)
        x: Vec3 # towards "true north"

        if ang in (0, HALF_REV):
            raise ValueError("Cannot create `EarthLocalFrame` in Earth's poles!")
        elif ang == RIGHT_ANGLE:
            x = north.copy()
        else:
            x = (north / cos(ang) - z).normalize()
            if ang > RIGHT_ANGLE:
                x *= -1
        
        y = x.cross(z) # towards "east"

        self._location = loc
        super().__init__(AbsoluteFrame.ITRF, x, y, z, itrf)
    
    def __repr__(self) -> str:
        gps = self._location.to_gps()
        return f"<EarthLocalFrame at GPS location {gps.dms}>"
    
    def get_horizontal_coordinates(self, of: "Coordinates") -> HorizontalCoordinates:
        """Transforms the given coordinates in the current `EarthLocalFrame`,
        and returns the new coordinates as a `HorizontalCoordinates` object."""
        itrf_pos = of.to_itrf()
        loc_pos = self.new_pos(itrf_pos)
        return HorizontalCoordinates(
            azimuth_deg=loc_pos.theta.m_as("°"),
            altitude_deg=loc_pos.delta.m_as("°")
        )

@dataclass(frozen=True)
class GPSCoordinates:
    """Coordinates of an object as GPS.
    Angles are given in degrees, distances in meters.
    
    Longitude ∈ [-180°, 180°]
    Latitude ∈ [-90°, 90°]
    """
    longitude_deg: float
    latitude_deg: float
    altitude_from_sea_level_meter: float = field(default=0., repr=False)
    location_name: str = field(default="", repr=False)

    def custom_repr(self) -> str:
        place = f"'{self.location_name}'" if self.location_name else ""
        return f"(lon={self.longitude_deg}°, lat={self.latitude_deg}° {place})"

    @lru_cache
    def to_itrf(self) -> Vec3:
        """Returns this coordinates as a vector in the ITRF frame."""
        theta = Q_(self.longitude_deg, "°")
        delta = Q_(self.latitude_deg, "°")
        rho = geocentric_radius_earth(delta) + Q_(self.altitude_from_sea_level_meter, "m")
        return Vec3.from_spherical(theta, delta, rho)
    
    @lru_cache
    def to_coordinates(self) -> "Coordinates":
        """Returns this GPS coordinates as a `Coordinates` object"""
        return Coordinates(AbsoluteFrame.ITRF, self.to_itrf())
    
    @lru_cache
    def earth_local_frame(self) -> EarthLocalFrame:
        """Returns the local frame corresponding to the GPS coordinates"""
        return EarthLocalFrame(
            self.to_coordinates()
        )
    
    @property
    def dms(self) -> str:
        """Representation of the GPS coordinates in DSM notation (degrees, minutes, seconds)

        Example: `39° 17′ N, 76° 36′ O`"""
        lon, lat = abs(self.longitude_deg), abs(self.latitude_deg)
        lon_dir = "E" if self.longitude_deg >= 0 else "O"
        lat_dir = "N" if self.latitude_deg >= 0 else "S"

        lon_deg, lon_deg_dec = divmod(lon, 1)
        lat_deg, lat_deg_dec = divmod(lat, 1)

        lon_min, lon_min_dec = divmod(lon_deg_dec * 60, 1)
        lat_min, lat_min_dec = divmod(lat_deg_dec * 60, 1)

        lon_sec, lon_sec_dec = divmod(lon_min_dec * 60, 1)
        lat_sec, lat_sec_dec = divmod(lat_min_dec * 60, 1)
        
        return (
            f"{lat_deg}° {lat_min}′ {lat_sec}″ {lat_dir}, "
            f"{lon_deg}° {lon_min}′ {lon_sec}″ {lon_dir}"
        )

class Coordinates:
    """This object corresponds to the 3D coordinates of an object.
    The coordinates can be expressed in any given frame.
    A velocity associated and an epoch can also be specified."""
    def __init__(self, frame: AnyFrame, pos: Vec3, vel: Vec3 | None = None, epoch: Time | None = None):
        self._positions: dict[AnyFrame, Vec3] = {frame: pos}
        self._velocities: dict[AnyFrame, Vec3] = {frame: vel} if vel is not None else {}
        self._epoch = epoch

        if isinstance(frame, RelativeFrame): # if given position is a relative frame, gather the old coordinates
            self._positions[frame.absolute_frame] = frame.old_pos(pos)
            if vel is not None:
                self._velocities[frame.absolute_frame] = frame.old_vel(vel)
    
    def __repr__(self) -> str:
        f, p = next((f, p) for f, p in self._positions.items())
        epoch_txt = f"epoch={self._epoch.human}, " if self._epoch is not None else ""

        return f"<Coordinates: {epoch_txt}{p} in {f}>"
            
    @property
    def epoch(self) -> Time:
        # FIXME: why raise ??
        if self._epoch is None:
            raise ValueError("Epoch was not set in coordinate instance")
        return self._epoch
    
    @property
    def vel_specified(self) -> bool:
        """Returns `True` if a velocity was specified during initialization"""
        return len(self._velocities) > 0
    
    def _get_existant_pos(self, frame: AnyFrame) -> Vec3 | None:
        return self._positions.get(frame, None)
    
    def _get_existant_vel(self, frame: AnyFrame) -> Vec3 | None:
        return self._velocities.get(frame, None)
    
    def pos_in_frame(self, frame: AnyFrame) -> Vec3:
        """Returns this coordinates as a vector in given frame."""
        p = self._get_existant_pos(frame)
        if p is not None:
            return p
        
        # FIXME: by construction, `Coordinates` object at least
        # has one position in an absolute frame, so it works.
        # But if we introduce another `AbsoluteFrame`, will
        # cause infinite recursive calls in some cases...
        if isinstance(frame, RelativeFrame):
            pos_absolute_frame = self.pos_in_frame(frame.absolute_frame)
            local = frame.new_pos(pos_absolute_frame)
            self._positions[frame] = local
            return local
        elif frame == AbsoluteFrame.ITRF:
            gcrf = self.pos_in_frame(AbsoluteFrame.GCRF)
            itrf = AbsoluteFrame.gcrf2itrf(gcrf, self.epoch)
            self._positions[AbsoluteFrame.ITRF] = itrf
            return itrf
        elif frame == AbsoluteFrame.GCRF:
            itrf = self.pos_in_frame(AbsoluteFrame.ITRF)
            gcrf = AbsoluteFrame.itrf2gcrf(itrf, self.epoch)
            self._positions[AbsoluteFrame.GCRF] = gcrf
            return gcrf
        
        raise RuntimeError(f"Given frame `{frame}` not recognized.")
    
    def to_gcrf(self) -> Vec3:
        """Returns this coordinates as a vector in GCRF frame."""
        return self.pos_in_frame(AbsoluteFrame.GCRF)
    
    def to_itrf(self) -> Vec3:
        """Returns this coordinates as a vector in ITRF frame."""
        return self.pos_in_frame(AbsoluteFrame.ITRF)
    
    @overload
    def to_horizontal_coordinates(self, earth_local_frame: EarthLocalFrame) -> HorizontalCoordinates:
        ...
    
    @overload
    def to_horizontal_coordinates(self, gps: GPSCoordinates) -> HorizontalCoordinates:
        ...
    
    def to_horizontal_coordinates(self, from_location: EarthLocalFrame | GPSCoordinates) -> HorizontalCoordinates:
        """Returns this coordinates as a `HorizontalCoordinates` object."""
        if isinstance(from_location, EarthLocalFrame):
            return from_location.get_horizontal_coordinates(self)
        elif isinstance(from_location, GPSCoordinates):
            return from_location.earth_local_frame().get_horizontal_coordinates(self)
        raise TypeError()

    def to_gps(self) -> GPSCoordinates:
        """Returns this coordinates as a `GPSCoordinates` object."""
        itrf = self.to_itrf()
        lon = itrf.theta
        lat = itrf.delta
        sea_level = geocentric_radius_earth(lat)
        alt = itrf.rho - sea_level
        return GPSCoordinates(lon.m_as("°"), lat.m_as("°"), alt.m_as("m"))
    
    def to_gcrf_vel(self) -> Vec3:
        """Returns this coordinates' associated velocity as a vector in GCRF frame."""
        if not self.vel_specified:
            raise RuntimeError("No velocity given in coordinate's instance.")
        if (vgcrf := self._get_existant_vel(AbsoluteFrame.GCRF)) is not None:
            return vgcrf
        elif (vitrf := self._get_existant_vel(AbsoluteFrame.ITRF)) is not None:
            vgcrf = AbsoluteFrame.to_itrf2gcrf(vitrf, self.epoch)
            self._velocities[AbsoluteFrame.GCRF] = vgcrf
            return vgcrf
        raise RuntimeError("Cannot get GCRF velocity from this coordinate. Object not instantied correctly?")
    
    def to_itrf_vel(self) -> Vec3:
        """Returns this coordinates' associated velocity as a vector in ITRF frame."""
        if len(self._velocities) == 0:
            raise RuntimeError("No velocity given in coordinate's instance.")
        if (vitrf := self._get_existant_vel(AbsoluteFrame.ITRF)) is not None:
            return vitrf
        if (vgcrf := self._get_existant_vel(AbsoluteFrame.GCRF)) is not None:
            vitrf = AbsoluteFrame.gcrf2itrf(vgcrf, self.epoch)
            self._velocities[AbsoluteFrame.ITRF] = vitrf
            return vitrf
        raise RuntimeError("Cannot get GCRF velocity from this coordinate. Object not instantied correctly?")
    
    def to_gcrf_pos_vel_tuple(self) -> PosVelGCRF:
        """Returns this coordinates' position and velocity as a tuple in GCRF frame.
        
        Distances are given in meters and velocities in meter/second.
        """
        pos = [i.m_as("m") for i in self.to_gcrf().xyz()]
        vel = [i.m_as("m/s") for i in self.to_gcrf_vel().xyz()]
        return (*pos, *vel)
    
    @staticmethod
    def from_gcrf_pos_vel_tuple(
        gcrf_pos_vel_tuple: PosVelGCRF, 
        epoch: Time = None
    ) -> "Coordinates":
        """Creates a `Coordinates` object from a tuple that corresponds to the object's position and velocity.
        Distances are given in meters and velocities in meter/second.
        The tuple's frame is considered as GCRF.
        """
        gcrf_x, gcrf_y, gcrf_z, gcrf_vx, gcrf_vy, gcrf_vz = gcrf_pos_vel_tuple

        gcrf_pos = Vec3(
            Q_(gcrf_x, "m"),
            Q_(gcrf_y, "m"),
            Q_(gcrf_z, "m")
        )
        gcrf_vel = Vec3(
            Q_(gcrf_vx, "m/s"),
            Q_(gcrf_vy, "m/s"),
            Q_(gcrf_vz, "m/s")
        )

        return Coordinates(AbsoluteFrame.GCRF, gcrf_pos, gcrf_vel, epoch)
