"""Implementation of 3D vectors, with units handling. 
Supports most of vectors operations."""

__pdoc__ = dict()

from collections import namedtuple
from functools import cache
from typing import Any, TypeVar, Iterator
from types import UnionType
from copy import copy, deepcopy

import pint
from pint import DimensionalityError
from leorbit.math import HALFPI
from leorbit.math import Q_
from pint import Unit
from pint.errors import DimensionalityError

from math import sqrt, cos, sin, acos, atan2, pi, tau

Matrix33AsTuple = namedtuple("Matrix33AsTuple", list("abcdefghi"))
__pdoc__["Matrix33AsTuple"] = """3x3 Matrix as a tuple. 
Each element should be a real number.
Corresponds to the concatenation of each lines:

```
[ a, b, c,
  d, e, f,
  g, h, i ]
```
"""
for letter in "abcdefghi":
    __pdoc__[f"Matrix33AsTuple.{letter}"] = False

class WrongOperandTypeError(Exception):
    def __init__(self, wrong_operand: Any, expected_operand_types: type | UnionType | str):
        super().__init__(f"Operand {wrong_operand} does not have expected type {expected_operand_types}")

def _raise_if_wrong_operand(operand: Any, expected_operand_types: type | UnionType):
    if isinstance(operand, expected_operand_types):
        return
    raise WrongOperandTypeError(operand, expected_operand_types)

UREG = pint.get_application_registry()

PINT_DIMLESS = Unit("1")

class Vec3:
    """3D Vector with common linear algebra. Projections can have any physical dimension."""
    def __init__(
        self: "Vec3",
        x: float | Q_,
        y: float | Q_,
        z: float | Q_
    ):
        xyz = [x, y, z]
        if any(not isinstance(e, float | int | Q_) for e in xyz):
            raise TypeError("...")
        
        self._x: Q_ = Q_(x).to_base_units()
        self._y: Q_ = Q_(y).to_base_units()
        self._z: Q_ = Q_(z).to_base_units()

        try:
            self._x + self._y + self._z
        except DimensionalityError:
            raise ValueError(f"{self._x} {self._y} {self._z} do not have same dimension.")

        self._theta: Q_ = None # atan2(y, x)
        self._delta: Q_ = None # atan2(z, [x²+y²])
        self._rho2: Q_ = None # x²+y²+z²
        self._rho: Q_ = None # sqrt(rho)
        
        self._unit: Unit = self._x.units
    
    @staticmethod
    def from_spherical(theta: Q_ | float, delta: Q_ | float, rho: Q_ | float) -> "Vec3":
        """
        Creates and returns a 3D vector from spherical coordinates. 

        Uses "radius-longitude-latitude" convention, [see in Wikipedia.](https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques#Convention_rayon-longitude-latitude))

        Arguments
        ---------
        - `theta:` Longitude angle (θ) from given convention. If is a `pint.Quantity`, must have angle dimension.
        - `delta:` Latitude angle (δ) from given convention. If is a `pint.Quantity`, must have angle dimension.
        - `rho:` Radius (ρ) from given convention.
        """
        cos_delta = cos(delta)
        return Vec3(
            rho * cos(theta) * cos_delta,
            rho * sin(theta) * cos_delta,
            rho * sin(delta)
        )

    def __repr__(self: "Vec3") -> str:
        return f"Vec3(x={self._x}, y={self._y}, z={self._z})"
    
    @property
    def unit(self: "Vec3") -> str:
        """Returns the units of this vector as a string."""
        return str(self._unit)
    
    def check(self: "Vec3", dimension: str) -> bool:
        """Returns `True` if given dimension matches vector's units."""
        return self.x.check(dimension)

    @property
    def x(self) -> Q_:
        """X projection of this vector."""
        return self._x

    @property
    def y(self) -> Q_:
        """Y projection of this vector."""
        return self._y

    @property
    def z(self) -> Q_:
        """Z projection of this vector."""
        return self._z
    
    @property
    def theta(self) -> Q_:
        """Angle between the projection of `self` on the (xy) plane, and the x axis.
        In [-π, π] range
        """
        if self._theta is None:
            ang = atan2(self._y.m, self._x.m)
            self._theta = Q_(ang, "rad")
        return self._theta
    
    @property
    def delta(self) -> Q_:
        """The complementary angle between `self` and z axis.
        In [-π/2, π/2] range
        """
        if self._delta is None:
            xy = (self._x.m**2 + self._y.m**2)**.5
            ang = atan2(self._z.m, xy)
            self._delta = Q_(ang, "rad")
        return self._delta
    
    @property
    def rho(self) -> Q_:
        """The magnitude of this vector."""
        if self._rho is None:
            self._rho2 = self._x**2 + self._y**2 + self._z**2
            self._rho = self._rho2**.5
        return self._rho
    
    @property
    def rho2(self) -> Q_:
        """The magnitude of this vector, *squared.*"""
        if self._rho2 is None:
            self._rho2 = self._x**2 + self._y**2 + self._z**2
            self._rho = self._rho2**.5
        return self._rho2
    
    def sqr(self) -> Q_:
        """Returns the squared magnitude of this vector."""
        return self.rho2
        
    def __abs__(self) -> Q_:
        """Returns the magnitude of the vector."""
        return self.rho
    
    def xyz(self) -> Iterator[Q_]:
        """Returns an iterator that will iterate throught 
        the X, Y, Z coordinates of this vector."""
        for i in (self._x, self._y, self._z):
            yield i

    def __copy__(self: "Vec3") -> "Vec3":
        return Vec3(
            copy(self._x), 
            copy(self._y), 
            copy(self._z)
        )
    
    def __deepcopy__(self: "Vec3", *args, **kwargs) -> "Vec3":
        return Vec3(
            deepcopy(self._x), 
            deepcopy(self._y), 
            deepcopy(self._z)
        )

    def __eq__(self: "Vec3", other: "Vec3") -> bool:
        _raise_if_wrong_operand(other, Vec3)
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self: "Vec3", other: "Vec3") -> bool:
        return not self == other

    ### Dot product ###
    def dot(vec: "Vec3", other: "Vec3") -> Q_:
        """Returns the dot product of two given vectors.
        Can also use the `@` operator for this function.
        
        Example:
        --------
        ```python
        a = Vec3(1, 2, 3)
        b = Vec3(-1, 0, 1)
        c = Vec3(1, .5, .25)

        a.dot(b)
        >>> 2

        Vec3.dot(b, c)
        >>> 0.75

        a @ c
        >>> 2.75
        ```"""
        _raise_if_wrong_operand(other, Vec3)
        return Q_(vec.x * other.x + vec.y * other.y + vec.z * other.z)

    # v @ w
    def __matmul__(self: "Vec3", other: "Vec3") -> Q_:
        return self.dot(other)

    ### ########### ###

    ### Cross product ###
    def cross(vec: "Vec3", other: "Vec3") -> "Vec3":
        """Returns the cross product of two given vectors. Can also use the `&` operator for this function.
        
        Example:
        --------
        ```python
        a = Vec3(1, 2, 3)
        b = Vec3(-1, 0, 1)
        c = Vec3(1, .5, .25)

        a.cross(b)
        >>> Vec3(x=2, y=-4, z=2)

        Vec3.cross(b, c)
        >>> Vec3(x=-0.5, y=1.25, z=-0.5)

        a & c
        >>> Vec3(x=-1.0, y=2.75, z=-1.5)
        ```
        """
        _raise_if_wrong_operand(other, Vec3)
        return Vec3(
            vec.y * other.z - vec.z * other.y,
            -vec.x * other.z + vec.z * other.x,
            vec.x * other.y - vec.y * other.x
        )

    # v & w
    def __and__(self: "Vec3", other: "Vec3") -> "Vec3":
        return self.cross(other)

    ### ############# ###
    
    def normalize(self: "Vec3") -> "Vec3":
        """Returns normalized vector: it is scaled so that its magnitude is equal to 1. 
        
        Also if the vector has a physical dimension, it will be cancelled and become dimensionless.
        """
        norm = abs(self)
        return self * norm**-1 if norm != 1 else self

    def __add__(self: "Vec3", other: "Vec3") -> "Vec3":
        _raise_if_wrong_operand(other, Vec3)
        return Vec3(self._x + other._x, self._y + other._y, self._z + other._z)

    def __iadd__(self: "Vec3", other: "Vec3") -> None:
        return self + other

    def __neg__(self: "Vec3") -> "Vec3":
        return Vec3(-self._x, -self._y, -self._z)

    def __sub__(self: "Vec3", other: "Vec3") -> "Vec3":
        _raise_if_wrong_operand(other, Vec3)
        return Vec3(self._x - other._x, self._y - other._y, self._z - other._z)

    def __isub__(self: "Vec3", other: "Vec3") -> "Vec3":
        return self - other

    def __mul__(self: "Vec3", other: Q_ | int | float) -> "Vec3":
        q = Q_(other)
        return Vec3(self.x * q, self.y * q, self.z * q)

    def __rmul__(self: "Vec3", other: Q_ | int | float) -> "Vec3":
        return self.__mul__(other)

    def __imul__(self: "Vec3", other: Q_) -> "Vec3":
        return self.__mul__(other)

    def __truediv__(self: "Vec3", other: Q_) -> "Vec3":
        return self.__mul__(1 / other)

    def __itruediv__(self: "Vec3", other: Q_) -> "Vec3":
        return self.__imul__(1 / other)

    def angle(a: "Vec3", b: "Vec3") -> Q_:
        """Returns the angle between the two given vectors.
        Returned angle is in `[0;π]` range.
        """
        if a == b:
            return Q_(0, "rad")
        cos_angle = Vec3.dot(a, b) / (a.rho * b.rho)
        if cos_angle >= 1:
            return Q_(0, "rad")
        elif cos_angle <= -1:
            return Q_(180, "°")
        return Q_(acos(cos_angle), "rad")
    
    def rotate_axis_angle(self: "Vec3", axis: "Vec3", angle: float | Q_) -> "Vec3":
        """Rotates a given vector along a given rotation axis, by a given angle,
        and returns the result.
        
        Arguments
        ---------
        `self: Vec3` Vector to be rotated.
        `axis: Vec3` Rotation axis. Must be normalized (aka has no dimension and has magnitude equal to 1).
        `angle: float | Q_` Rotation angle. If is a `pint.Quantity`, must have angle dimension.
        """
        if axis.sqr() != 1:
            raise ValueError(f"axis vector should be normalized")

        c = cos(angle)
        s = sin(angle)
        if isinstance(angle, Q_):
            c = c.m_as("rad")
            s = s.m_as("rad")
        omc = 1 - c
        
        x = axis.x.m
        y = axis.y.m
        z = axis.z.m

        xs = x * s
        ys = y * s
        zs = z * s
        xy_omc = x * y * omc
        xz_omc = x * z * omc
        yz_omc = y * z * omc

        rot_mat = [
            c + omc * x*x, xy_omc - zs, xz_omc + ys,
            xy_omc + zs, c + omc * y*y, yz_omc - xs,
            xz_omc - ys, yz_omc + xs, c + omc * z*z
        ]

        return self.apply_mat_list(rot_mat)

    def to_array(self) -> list:
        """Returns this vector as an array of `pint.Quantity`"""
        return Q_([self.x, self.y, self.z])

    @staticmethod
    def zero(unit: str = "") -> "Vec3":
        """Creates and returns the `(0, 0, 0)` vector,
        with a specified unit, if given."""
        return Vec3(Q_(0, unit), Q_(0, unit), Q_(0, unit))

    @staticmethod
    def xaxis(unit: str = "") -> "Vec3":
        """Creates and returns the `(1, 0, 0)` vector,
        with a specified unit, if given."""
        return Vec3(Q_(1, unit), Q_(0, unit), Q_(0, unit))

    @staticmethod
    def yaxis(unit: str = "") -> "Vec3":
        """Creates and returns the `(0, 1, 0)` vector,
        with a specified unit, if given."""
        return Vec3(Q_(0, unit), Q_(1, unit), Q_(0, unit))

    @staticmethod
    def zaxis(unit: str = "") -> "Vec3":
        """Creates and returns the `(0, 0, 1)` vector,
        with a specified unit, if given."""
        return Vec3(Q_(0, unit), Q_(0, unit), Q_(1, unit))

    @staticmethod
    def one(unit: str = "") -> "Vec3":
        """Creates and returns the `(1, 1, 1)` vector,
        with a specified unit, if given."""
        return Vec3(Q_(1, unit), Q_(1, unit), Q_(1, unit))

    def rt_z(self: "Vec3", angle: float | Q_) -> "Vec3":
        """Rotates vector along Z axis of an angle given by user"""
        return Vec3(
            self._x * cos(angle) - self._y * sin(angle),
            self._x * sin(angle) + self._y * cos(angle),
            self._z
        )
    
    def apply_mat_list(self: "Vec3", mat: Matrix33AsTuple) -> "Vec3":
        """Applies a given matrix to given vector. Computation follows usual matrix product.

        Arguments
        ---------
        `self: Vec3` Vector to be transformed
        `mat: Matrix33AsTuple` 3x3 matrix to apply the given vector.

        Returns
        -------
        `Vec3` Vector transformed by given matrix.
        """
        try:
            return Vec3(
                mat[0] * self._x + mat[1] * self._y + mat[2] * self._z,
                mat[3] * self._x + mat[4] * self._y + mat[5] * self._z,
                mat[6] * self._x + mat[7] * self._y + mat[8] * self._z,
            )
        except IndexError:
            raise ValueError("Matrix given does not have 3x3 shape!")

def inverse_mat3x3(mat: Matrix33AsTuple) -> Matrix33AsTuple:
    """Computes and returns the inverse of the given matrix.
    
    Raises:
    -------
    `ZeroDivisionError:` if given matrix has no inverse.
    """
    a, b, c, d, e, f, g, h, i = mat
    det_mat = a*e*i + d*h*c + g*b*f - c*e*g - f*h*a - i*b*d
    if det_mat == 0:
        raise ZeroDivisionError(f"Given matrix {mat} has no inverse.")
    idm = det_mat**-1

    aa = (e*i - f*h)
    bb = -(b*i - c*h)
    cc = (b*f - c*e)
    dd = -(d*i - f*g)
    ee = (a*i - c*g)
    ff = -(a*f - c*d)
    gg = (d*h - e*g)
    hh = -(a*h - b*g)
    ii = (a*e - b*d)

    return Matrix33AsTuple(
        aa * idm, bb * idm, cc * idm,
        dd * idm, ee * idm, ff * idm,
        gg * idm, hh * idm, ii * idm
    )

