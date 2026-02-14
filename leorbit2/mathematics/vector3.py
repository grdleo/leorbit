from __future__ import annotations

from functools import cached_property
from itertools import repeat
from math import pi
from typing import Any, ClassVar, Generic, Literal, Never, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import D, DimCoords, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.functions import acos, atan2, cos, sin, sqrt, square
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar, ScalarArray, TensorScalar
from leorbit2.mathematics.tensor import Tensor, ensure_same_dimensions

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]


def _is_simple_number(value: object) -> bool:
    return isinstance(value, (int, float, np.floating, np.integer))


def all_simple_numbers(els: list[object]) -> bool:
    return all(_is_simple_number(el) for el in els)


def all_scalars_numbers(els: list[object]) -> bool:
    return all(isinstance(el, TensorScalar) for el in els)


def _scalar_from_values(dim: type, values: np.ndarray) -> TensorScalar[Any]:
    flat = np.asarray(values).reshape(-1)
    if flat.size == 1:
        return Scalar[dim].new(flat.item())
    return ScalarArray[dim].new(flat)


class TensorVector3(Tensor[SomeDim], Generic[SomeDim]):
    """Core vector tensor type supporting one or many 3D vectors."""

    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[TensorVector3]:
        return cast(type[TensorVector3], super().dimensionalize(dim_coords))

    @classmethod
    def new_from_vectors(cls, *vectors: TensorVector3[SomeDim]) -> TensorVector3[SomeDim]:
        if not vectors:
            raise ValueError("...")
        v0 = vectors[0]
        try:
            for vector in vectors[1:]:
                ensure_same_dimensions(v0, vector)
        except RuntimeError as exc:
            raise ValueError("...") from exc
        return cls(np.concatenate([v._values for v in vectors], axis=1))

    @classmethod
    def new_from_components(
        cls,
        x: TensorScalar[SomeDim],
        y: TensorScalar[SomeDim],
        z: TensorScalar[SomeDim],
    ) -> TensorVector3[SomeDim]:
        ensure_same_dimensions(x, y, z)
        x_vals = np.asarray(x.base_unit_value).reshape(-1)
        y_vals = np.asarray(y.base_unit_value).reshape(-1)
        z_vals = np.asarray(z.base_unit_value).reshape(-1)
        if not (x_vals.size == y_vals.size == z_vals.size):
            raise ValueError("...")
        return cls(np.stack([x_vals, y_vals, z_vals]))

    @classmethod
    def new_from_single_vector(cls, vector: TensorVector3[SomeDim], size: int) -> TensorVector3[SomeDim]:
        return cls.new_from_vectors(*repeat(vector, size))

    @property
    def size(self) -> int:
        _, s = self._values.shape
        return int(s)

    def cast(self, dim: type[SomeOtherDim]) -> TensorVector3[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self  # type: ignore
        raise RuntimeError("Cannot cast")

    @cached_property
    def x(self) -> TensorScalar[SomeDim]:
        return _scalar_from_values(self.dim, self._values[0, :])

    @cached_property
    def y(self) -> TensorScalar[SomeDim]:
        return _scalar_from_values(self.dim, self._values[1, :])

    @cached_property
    def z(self) -> TensorScalar[SomeDim]:
        return _scalar_from_values(self.dim, self._values[2, :])

    @cached_property
    def length(self) -> TensorScalar[SomeDim]:
        return _scalar_from_values(self.dim, np.sum(self._values ** 2, axis=0) ** 0.5)

    def __getitem__(self, index: int) -> Vector3[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("...")
        return Vector3[self.dim].new(
            x=self._values[0, index],
            y=self._values[1, index],
            z=self._values[2, index],
        )

    @cached_property
    def theta(self) -> TensorScalar[D.Angle]:
        return cast(TensorScalar[D.Angle], atan2(self.y, self.x))

    @cached_property
    def delta(self) -> TensorScalar[D.Angle]:
        xy = sqrt(square(self.x) + square(self.y))
        return cast(TensorScalar[D.Angle], atan2(self.z, xy))

    def angle(self, o: TensorVector3[SomeDim]) -> TensorScalar[D.Angle]:
        cos_angle = cast(TensorScalar[D.Dimless], self.dot(o) / (self.length * o.length))
        if self.size == 1:
            if cast(Scalar[D.Dimless], cos_angle) >= 1:
                return cast(TensorScalar[D.Angle], 0 * Quantity.rad)
            if cast(Scalar[D.Dimless], cos_angle) <= -1:
                return cast(TensorScalar[D.Angle], pi * Quantity.rad)
        return cast(TensorScalar[D.Angle], acos(cos_angle))

    def normalized(self) -> TensorVector3[D.Dimless]:
        return cast(TensorVector3[D.Dimless], self / self.length)

    @staticmethod
    def from_spherical(
        theta: TensorScalar[D.Angle],
        delta: TensorScalar[D.Angle],
        rho: TensorScalar[SomeOtherDim],
    ) -> TensorVector3[SomeOtherDim]:
        cos_delta = cos(delta)
        sin_delta = sin(delta)
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        x = rho * cos_theta * cos_delta
        y = rho * sin_theta * cos_delta
        z = rho * sin_delta
        return TensorVector3[rho.dim].new_from_components(
            cast(TensorScalar[SomeOtherDim], x),
            cast(TensorScalar[SomeOtherDim], y),
            cast(TensorScalar[SomeOtherDim], z),
        )

    @overload
    def dot(self: TensorVector3[SomeDim], o: TensorVector3[SomeDim]) -> TensorScalar[ProductDim[SomeDim, SomeDim]]: ...

    @overload
    def dot(self: TensorVector3[SomeDim], o: TensorVector3[SomeOtherDim]) -> TensorScalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    def dot(self, o: object) -> TensorScalar[Any]:
        if not isinstance(o, TensorVector3):
            raise TypeError("dot product requires two TensorVector3 instances")
        vals = np.sum(self._values * o._values, axis=0)
        return _scalar_from_values((self.dim_coords * o.dim_coords).to_dimension(), vals)

    @overload
    def cross(self: TensorVector3[D.Dimless], o: TensorVector3[D.Dimless]) -> TensorVector3[D.Dimless]: ...

    @overload
    def cross(self: TensorVector3[D.Dimless], o: TensorVector3[SomeOtherDim]) -> TensorVector3[SomeOtherDim]: ...

    @overload
    def cross(self: TensorVector3[SomeDim], o: TensorVector3[SomeOtherDim]) -> TensorVector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    def cross(self, o: object) -> TensorVector3[Any]:
        if not isinstance(o, TensorVector3):
            raise TypeError("cross product requires two TensorVector3 instances")
        return TensorVector3.dimensionalize(self.dim_coords * o.dim_coords)(
            np.cross(self._values, o._values, axis=0)
        )

    @overload
    def __add__(self: TensorVector3[SomeDim], o: TensorVector3[SomeDim] | TensorScalar[SomeDim]) -> TensorVector3[SomeDim]: ...

    def __add__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__add__(o))

    @overload
    def __radd__(self: TensorVector3[D.Dimless], o: Number) -> TensorVector3[D.Dimless]: ...

    @overload
    def __radd__(self: TensorVector3[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__radd__(o))

    @overload
    def __sub__(self: TensorVector3[SomeDim], o: TensorVector3[SomeDim] | TensorScalar[SomeDim]) -> TensorVector3[SomeDim]: ...

    def __sub__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__sub__(o))

    @overload
    def __rsub__(self: TensorVector3[D.Dimless], o: Number) -> TensorVector3[D.Dimless]: ...

    @overload
    def __rsub__(self: TensorVector3[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__rsub__(o))

    @overload
    def __mul__(self: TensorVector3[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorVector3[D.Dimless]: ...

    @overload
    def __mul__(self: TensorVector3[D.Dimless], o: TensorScalar[SomeOtherDim]) -> TensorVector3[SomeOtherDim]: ...

    @overload
    def __mul__(self: TensorVector3[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorVector3[SomeDim]: ...

    @overload
    def __mul__(self: TensorVector3[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorVector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __mul__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__mul__(o))

    @overload
    def __rmul__(self: TensorVector3[D.Dimless], o: Number) -> TensorVector3[D.Dimless]: ...

    @overload
    def __rmul__(self: TensorVector3[SomeDim], o: Number) -> TensorVector3[SomeDim]: ...

    def __rmul__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__rmul__(o))

    @overload
    def __truediv__(self: TensorVector3[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorVector3[D.Dimless]: ...

    @overload
    def __truediv__(self: TensorVector3[D.Dimless], o: TensorScalar[SomeDim]) -> TensorVector3[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: TensorVector3[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorVector3[SomeDim]: ...

    @overload
    def __truediv__(self: TensorVector3[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorVector3[QuotientDim[SomeDim, SomeOtherDim]]: ...

    def __truediv__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: TensorVector3[D.Dimless], o: Number) -> TensorVector3[D.Dimless]: ...

    @overload
    def __rtruediv__(self: TensorVector3[SomeDim], o: Number) -> TensorVector3[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> TensorVector3[Any]:
        return cast(TensorVector3[Any], super().__rtruediv__(o))

    @overload
    def __matmul__(self: TensorVector3[D.Dimless], o: TensorVector3[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __matmul__(self: TensorVector3[SomeDim], o: TensorVector3[SomeOtherDim]) -> TensorScalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> TensorScalar[Any]:
        if not isinstance(o, TensorVector3):
            raise TypeError("...")
        return self.dot(o)

    def __eq__(self, o: object) -> np.ndarray:
        return np.all(super().__eq__(o)._values, axis=0)

    def __neq__(self, o: object) -> np.ndarray:
        return ~self.__eq__(o)

    def __le__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for vectors")

    def __lt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for vectors")

    def __ge__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for vectors")

    def __gt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for vectors")


class Vector3(TensorVector3[SomeDim], Generic[SomeDim]):
    """End-user convenience single-vector class."""

    O: ClassVar[Vector3[D.Dimless]]
    X: ClassVar[Vector3[D.Dimless]]
    Y: ClassVar[Vector3[D.Dimless]]
    Z: ClassVar[Vector3[D.Dimless]]
    ONE: ClassVar[Vector3[D.Dimless]]

    @classmethod
    def new(cls, x: Number | TensorScalar[Any], y: Number | TensorScalar[Any], z: Number | TensorScalar[Any]) -> Vector3:
        xyz = [x, y, z]
        if all(_is_simple_number(v) for v in xyz):
            return cls(np.array(xyz).reshape((3, 1)))
        if all(isinstance(v, TensorScalar) for v in xyz):
            ensure_same_dimensions(cast(TensorScalar[Any], x), cast(TensorScalar[Any], y), cast(TensorScalar[Any], z))
            vals = [cast(TensorScalar[Any], v).base_unit_value for v in xyz]
            vals = [np.asarray(v).item() for v in vals]
            return cls(np.array(vals).reshape((3, 1)))
        raise RuntimeError("...")


class Vector3Array(TensorVector3[SomeDim], Generic[SomeDim]):
    """End-user convenience vector-array class."""


_vector3_dimless = Vector3[D.Dimless]
Vector3.O = _vector3_dimless.new(0, 0, 0)
Vector3.X = _vector3_dimless.new(1, 0, 0)
Vector3.Y = _vector3_dimless.new(0, 1, 0)
Vector3.Z = _vector3_dimless.new(0, 0, 1)
Vector3.ONE = _vector3_dimless.new(1, 1, 1)
