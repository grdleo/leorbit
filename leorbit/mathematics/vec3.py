from functools import cached_property
from math import acos, atan2, cos, pi, sin
from typing import Self
from pint import Quantity as Q_
from pint import Unit, DimensionalityError

from mathematics.custom import AnyNumber, Number
from mathematics.units import ACC_UNIT, POS_UNIT, UREG, VEL_UNIT

from mathematics.quaternions import Quaternion



def get_base_conversion(unit: Unit) -> tuple[Unit, float]:
    conv = (1.0 * unit).to_base_units()
    return conv.units, conv.magnitude if conv.units != unit else 1.

def get_conversion_factor(ufrom: Unit, uto: Unit) -> float:
    return (1.0 * ufrom).to(uto).magnitude

def pint_force_float(x: AnyNumber, unit: Unit) -> float:
    if isinstance(x, Q_):
        return x.to(unit).magnitude
    return x

class Vec3:
    def __init__(self, 
        x: Number, 
        y: Number, 
        z: Number, 
        unit: Unit = UREG.dimensionless
    ):
        if any((not isinstance(e, Number) for e in (x, y, z))):
            raise ValueError("Projections should be numeric, and the unit given as parameter. Use `Vec3.from_quantities` instead")
        
        bunit, fact = get_base_conversion(unit)

        self._x = x * fact
        self._y = y * fact
        self._z = z * fact
        self._unit = bunit
        self._ureg = self._unit._REGISTRY
    
    @staticmethod
    def from_quantities(x: AnyNumber, y: AnyNumber, z: AnyNumber) -> "Vec3":
        unit = x.to_base_units().units
        if any((not e.check(unit) for e in (x, y, z))):
            raise ValueError("Projections' units have to be from the same dimension.")
        
        return Vec3(
            x.m_as(unit),
            y.m_as(unit),
            z.m_as(unit),
            unit
        )
    
    @staticmethod
    def from_spherical(theta: AnyNumber, delta: AnyNumber, rho: AnyNumber) -> "Vec3":
        """
        Creates and returns a 3D vector from spherical coordinates. 

        Uses "radius-longitude-latitude" convention, [see in Wikipedia.](https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques#Convention_rayon-longitude-latitude))

        Arguments
        ---------
        - `theta:` Longitude angle (θ) from given convention. If is a `pint.Quantity`, must have angle dimension.
        - `delta:` Latitude angle (δ) from given convention. If is a `pint.Quantity`, must have angle dimension.
        - `rho:` Radius (ρ) from given convention.
        """
        θ = pint_force_float(theta, UREG.radians)
        δ = pint_force_float(delta, UREG.radians)
        ρ = pint_force_float(rho, UREG.meters)

        cos_delta = cos(δ)
        return Vec3(
            ρ * cos(θ) * cos_delta,
            ρ * sin(θ) * cos_delta,
            ρ * sin(δ),
            rho.unit if isinstance(rho, Q_) else UREG.dimensionless
        )
    
    def to(self, unit: Unit) -> tuple[float, float, float]:
        """Returns this vector as a `float` tuple, in the given unit"""
        fact = get_conversion_factor(self._unit, unit)
        return (
            self._x * fact,
            self._y * fact,
            self._z * fact
        )
    
    @cached_property
    def as_tuple(self) -> tuple[float, float, float]:
        return (
            self._x,
            self._y,
            self._z,
        )

    @cached_property
    def x(self) -> Q_:
        """X projection of this vector."""
        return self._x * self._unit

    @cached_property
    def y(self) -> Q_:
        """Y projection of this vector."""
        return self._y * self._unit

    @cached_property
    def z(self) -> Q_:
        """Z projection of this vector."""
        return self._z * self._unit
    
    @cached_property
    def xyz(self) -> tuple[Q_, Q_, Q_]:
        """The (x, y, z) projections as `pint.Quantity`."""
        return self.x, self.y, self.z
    
    @cached_property
    def unit(self) -> Unit:
        return self._unit
    
    @cached_property
    def theta(self) -> Q_:
        """Angle between the projection of `self` on the (xy) plane, and the x axis.
        In [-π, π] range
        """
        return atan2(self._y, self._x) * self._ureg.rad
    
    @cached_property
    def delta(self) -> Q_:
        """The complementary angle between `self` and z axis.
        In [-π/2, π/2] range
        """
        xy = (self._x**2 + self._y**2)**.5
        return atan2(self._z, xy) * self._ureg.rad
    
    @cached_property
    def _rho2_float(self) -> float:
        """The magnitude of this vector, *squared.* (as a float)"""
        return self._x**2 + self._y**2 + self._z**2
    
    @cached_property
    def rho2(self) -> Q_:
        """Magnitude of this vector, *squared.*"""
        return self._rho2_float * self._unit**2
    
    @cached_property
    def _rho_float(self) -> float:
        """The magnitude of this vector. (as a float)"""
        return self._rho2_float**.5
    
    @cached_property
    def rho(self) -> Q_:
        """Magnitude of this vector."""
        return self._rho_float * self._unit
        
    def __abs__(self) -> Q_:
        """Returns the magnitude of the vector."""
        return self.rho
    
    def __add__(self, other: "Vec3") -> "Vec3":
        if self.unit != other.unit:
            raise DimensionalityError(self.unit, other.unit)
        
        return Vec3(
            self._x + other._x,
            self._y + other._y,
            self._z + other._z,
            self._unit
        )
    
    def __sub__(self, other: "Vec3") -> "Vec3":
        if self.unit != other.unit:
            raise DimensionalityError(self.unit, other.unit)
        
        return Vec3(
            self._x - other._x,
            self._y - other._y,
            self._z - other._z,
            self._unit
        )
    
    def __mult__(self, other: AnyNumber) -> "Vec3":
        assert isinstance(other, AnyNumber)

        is_qt = isinstance(other, Q_)
        u = self._unit * (other.units if is_qt else 1)
        m = other.magnitude if is_qt else other

        return Vec3(
            self._x * m,
            self._y * m,
            self._z * m,
            u
        )
    
    def __truediv__(self, other: AnyNumber) -> Self:
        assert isinstance(other, AnyNumber)

        return self.__mul__(1/other)
    
    def dot(self, other: "Vec3") -> Q_:
        assert isinstance(other, Vec3)

        return (
            self._x * other._x
            + self._y * other._y
            + self._z * other._z
        ) * (self._unit * other._unit)
    
    def cross(self, other: "Vec3") -> "Vec3":
        assert isinstance(other, Vec3)

        return Vec3(
            self._y*other._z - other._y*self._z,
            other._x*self._z - self._x*other._z,
            self._x*other._y - other._x*self._y,
            self._unit * other._unit
        )
    
    def angle(a: Self, b: Self) -> Q_:
        """Returns the angle between the two given vectors.
        Returned angle is in `[0;π]` range.
        """
        if a == b:
            return 0 * UREG.radians
        cos_angle = Vec3.dot(a, b) / (a.rho * b.rho)
        if cos_angle >= 1:
            return 0 * UREG.radians
        elif cos_angle <= -1:
            return pi * UREG.radians
        return acos(cos_angle) * UREG.radians
    
    def normalized(self, keep_dim = False) -> Self:
        return Vec3(
            self._x / self._rho_float,
            self._y / self._rho_float,
            self._z / self._rho_float,
            self.unit if keep_dim else UREG.dimensionless
        ) 
    
    def __matmult__(self, other: int | float | Q_) -> "Vec3":
        return self.dot(other)
    
    def __xor__(self, other: "Vec3") -> "Vec3":
        return self.cross(other)

    def to_quaternion(self, scalar: float = 0.) -> Quaternion:
        return Quaternion(scalar, self._x, self._y, self._z)
    
    def rotate_quaternion(self, q: Quaternion) -> "Vec3":
        rotated = self.to_quaternion().conjugate(q)
        return Vec3(*rotated.bcd, self._unit)
    
    def rotate_axis_angle(self, axis: "Vec3", angle: AnyNumber) -> "Vec3":
        hangle: float = pint_force_float(angle, self._ureg.radians) / 2
        u: Vec3 = axis*sin(hangle)
        q = u.to_quaternion(cos(hangle))
        return self.rotate_quaternion(q)
    
    def rotate_zaxis(self, angle: AnyNumber) -> "Vec3":
        c, s = cos(angle), sin(angle)
        return Vec3(
            c * self._x + s * self._y,
            -s * self._x + c * self._y,
            self._z,
            self._unit
        )
    
    @cached_property
    def dynamic_vector(self) -> bool:
        """Returns `True` if vector is a 'dynamic' vector (aka. either position, velocity, acceleration)"""
        return self.unit in (POS_UNIT, VEL_UNIT, ACC_UNIT)
    
VEC_0 = ORIGIN = Vec3(0, 0, 0)
XAXIS = Vec3(1, 0, 0)
YAXIS = Vec3(0, 1, 0)
ZAXIS = Vec3(0, 0, 1)