from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from math import pi
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics import sqrt
from leorbit2.mathematics.dimensions import Dim, DimCoords, D, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.functions import atan2, acos, cos, sin, square
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.tensor import Tensor, ensure_same_dimensions

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)
    
NumberOrScalarT = TypeVar("NumberOrScalarT", bound=Number | Scalar)

def all_simple_numbers(els: list[object]) -> TypeGuard[list[Number]]:
    return all(isinstance(el, Number) for el in els)

def all_scalars_numbers(els: list[object]) -> TypeGuard[list[Scalar]]:
    return all(isinstance(el, Scalar) for el in els)

class Vector3(Generic[SomeDim], Tensor[SomeDim]):
    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[Vector3]:
        return cast(
            type[Vector3],
            super().dimensionalize(dim_coords)
        )
    
    O: ClassVar[Vector3[D.Dimless]]
    X: ClassVar[Vector3[D.Dimless]]
    Y: ClassVar[Vector3[D.Dimless]]
    Z: ClassVar[Vector3[D.Dimless]]
    ONE: ClassVar[Vector3[D.Dimless]]

    @classmethod
    def new(cls, x: NumberOrScalarT, y: NumberOrScalarT, z: NumberOrScalarT):
        v: TensorData
        xyz = cast(list[object], [x, y, z])

        if all_simple_numbers(xyz):
            v = np.array(xyz).reshape((3,1))
        elif all_scalars_numbers(xyz):
            ensure_same_dimensions(*xyz)
            v = np.array([el.base_unit_value for el in xyz]).reshape((3,1))
        else:
            raise RuntimeError("...")

        return cls(v)

    def cast(self, dim: type[SomeOtherDim]) -> Vector3[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @cached_property
    def x(self) -> Scalar[SomeDim]:
        return Scalar[self.dim](self._values[0])
    
    @cached_property
    def y(self) -> Scalar[SomeDim]:
        return Scalar[self.dim](self._values[1])
    
    @cached_property
    def z(self) -> Scalar[SomeDim]:
        return Scalar[self.dim](self._values[2])
    
    @cached_property
    def length(self) -> Scalar[SomeDim]:
        return Scalar[self.dim](
            np.sum(self._values ** 2) ** .5
        )
    
    @cached_property
    def theta(self) -> Scalar[D.Angle]:
        """Angle between the projection of `self` on the (xy) plane, and the x axis.
        In [-π, π] range
        """

        return atan2(self.y, self.x)
    
    @cached_property
    def delta(self) -> Scalar[D.Angle]:
        """The complementary angle between `self` and z axis.
        In [-π/2, π/2] range
        """

        xy = sqrt(square(self.x) + square(self.y))
        return atan2(self.z, xy)
    
    def angle(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Scalar[D.Angle]:
        """Returns the angle between the two given vectors.
        Returned angle is in `[0;π]` range.
        """
        if self == o:
            return 0 * Quantity.rad
        cos_angle = cast(Scalar[D.Dimless], self.dot(o) / (self.length * o.length))
        if cos_angle >= 1:
            return 0 * Quantity.rad
        elif cos_angle <= -1:
            return pi * Quantity.rad
        return acos(cos_angle)
    
    def normalized(self) -> Vector3[D.Dimless]:
        l = self.length
        return self / l
    
    @staticmethod
    def from_spherical(theta: Scalar[D.Angle], delta: Scalar[D.Angle], rho: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]:
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
        sin_delta = sin(delta)
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        return Vector3[rho.dim].new(
            x = rho * cos_theta * cos_delta,
            y = rho * sin_theta * cos_delta,
            z = rho * sin_delta
        )
    
    @overload
    def dot(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Scalar[ProductDim[SomeDim, SomeDim]]: ...

    @overload
    def dot(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...
    
    def dot(self, o: object) -> Scalar[Any]:
        if not isinstance(o, Vector3):
            raise TypeError("dot product requires two Vector3 instances")

        return Scalar.dimensionalize(
            self.dim_coords * o.dim_coords
        ).new(
            np.sum(self._values * o._values) ** .5
        )
    
    @overload
    def cross(self: Vector3[D.Dimless], o: Vector3[D.Dimless]) -> Vector3[D.Dimless]: ...
    
    @overload
    def cross(self: Vector3[D.Dimless], o: Vector3[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...
    
    @overload
    def cross(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...
    
    def cross(self, o: object) -> Vector3[Any]:
        if not isinstance(o, Vector3):
            raise TypeError("cross product requires two Vector3 instances")

        return Vector3.dimensionalize(
            self.dim_coords * o.dim_coords
        )(
            np.cross(self._values, o._values, axis=0)
        )
    
    def __repr__(self) -> str:
        return (
            f"Vector3[D.{self.dim.__class__.__name__}]("
            f"x={self._values[0]}, "
            f"y={self._values[1]}, "
            f"z={self._values[2]}"
            ")"
        )

    ### + OPERATOR ###

    @overload
    def __add__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __add__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __add__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__add__(o))

    @overload
    def __radd__(self: Vector3[D.Dimless], o: Vector3[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __radd__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __radd__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __radd__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __sub__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __sub__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__sub__(o))

    @overload
    def __rsub__(self: Vector3[D.Dimless], o: Vector3[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __rsub__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __rsub__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __rsub__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __mul__(self: Vector3[D.Dimless], o: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[SomeDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Vector3[Any]: ...

    def __mul__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __rmul__(self: Vector3[D.Dimless], o: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __rmul__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[SomeDim]: ...

    @overload
    def __rmul__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __rmul__(self, o: Scalar) -> Vector3[Any]: ...

    def __rmul__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __truediv__(self: Vector3[D.Dimless], o: Scalar[SomeDim]) -> Vector3[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[SomeDim]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Vector3[Any]: ...

    def __truediv__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __rtruediv__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rtruediv__(o))
    
    ### @ OPERATOR (dot product) ###

    @overload
    def __matmul__(self: Vector3[D.Dimless], o: Vector3[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __matmul__(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> Scalar[Any]:
        if isinstance(o, Vector3):
            return self.dot(o)
        
        raise TypeError("...")
    
    ################

    def __eq__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> bool:
        return bool(np.all(super().__eq__(o)._values))
    
    def __neq__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> bool:
        return not self.__eq__(o)
    
    def __le__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3")
    
    def __lt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3")
    
    def __ge__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3")
    
    def __gt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3")
    


_vector3_dimless = Vector3[D.Dimless]
Vector3.O = _vector3_dimless.new(0, 0, 0)
Vector3.X = _vector3_dimless.new(1, 0, 0)
Vector3.Y = _vector3_dimless.new(0, 1, 0)
Vector3.Z = _vector3_dimless.new(0, 0, 1)
Vector3.ONE = _vector3_dimless.new(1, 1, 1)