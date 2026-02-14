from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from itertools import repeat
from multiprocessing import Value
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.functions import sqrt
from leorbit2.mathematics.dimensions import Dim, DimCoords, D, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.functions import atan2, acos, cos, sin, square
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.scalar_array import ScalarArray
from leorbit2.mathematics.tensor import Tensor, TensorType, ensure_same_dimensions
from leorbit2.mathematics.vector3 import Vector3

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)
    
NumberOrScalarT = TypeVar("NumberOrScalarT", bound=Number | Scalar)

def all_vector3_same_dim(els: list[object]) -> TypeGuard[list[Vector3]]:
    """Return whether all elements are vectors sharing the same dimension."""
    v0, *others = els
    if not isinstance(v0, Vector3):
        return False
    
    return all(
        isinstance(el, Vector3) and v0.ensure_compatible_dimensions(el)
        for el in others
    )

class Vector3Array(Tensor[SomeDim], Generic[SomeDim]):
    """Array of 3D vectors sharing the same dimension."""

    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[Vector3Array]:
        """Return a vector-array class bound to ``dim_coords``."""
        return cast(
            type[Vector3Array],
            super().dimensionalize(dim_coords)
        )
    
    @classmethod
    def new_from_vectors(cls, *vectors: Vector3[SomeDim]):
        """Create an array by concatenating vectors column-wise."""
        if not vectors:
            raise ValueError("...")
        list_vectors = cast(list[object], list(vectors))

        if not all_vector3_same_dim(list_vectors):
            raise ValueError("...")
        
        return cls(
            np.concatenate(
                [v._values for v in list_vectors],
                axis=1
            )
        )
    
    @classmethod
    def new_from_components(cls, x: ScalarArray[SomeDim], y: ScalarArray[SomeDim], z: ScalarArray[SomeDim]):
        """Create an array from x, y and z components."""
        if not (x.dim == y.dim == z.dim):
            raise ValueError("...")
        
        return cls(
            np.stack([x._values, y._values, z._values])
        )
    
    @classmethod
    def new_from_single_vector(cls, vector: Vector3[SomeDim], size: int):
        """Create an array by repeating one vector ``size`` times."""
        return cls.new_from_vectors(*repeat(vector, size))
    
    @property
    def tensor_type(self) -> TensorType:
        return TensorType.VECTOR3_ARRAY
        

    def cast(self, dim: type[SomeOtherDim]) -> Vector3[SomeOtherDim]:
        """Type-cast to another dimension when coordinates are identical."""
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @cached_property
    def x(self) -> ScalarArray[SomeDim]:
        """X components as a scalar array."""
        return ScalarArray[self.dim](self._values[0,:])
    
    @cached_property
    def y(self) -> ScalarArray[SomeDim]:
        """Y components as a scalar array."""
        return ScalarArray[self.dim](self._values[1,:])
    
    @cached_property
    def z(self) -> ScalarArray[SomeDim]:
        """Z components as a scalar array."""
        return ScalarArray[self.dim](self._values[2,:])
    
    @cached_property
    def length(self) -> ScalarArray[SomeDim]:
        """Euclidean norm for each vector in the array."""
        return ScalarArray[self.dim](
            np.sum(self._values ** 2, axis=0) ** .5
        )
    
    @property
    def size(self) -> int:
        """Number of vectors stored in the array."""
        _, s = self._values.shape
        return s
    
    def __getitem__(self, index: int) -> Vector3[SomeDim]:
        """Return vector at ``index``."""
        if index < 0 or index >= self.size:
            raise KeyError("...")
        
        return Vector3[self.dim].new(
            x=self._values[0, index],
            y=self._values[1, index],
            z=self._values[2, index],
        )
    
    @cached_property
    def theta(self) -> ScalarArray[D.Angle]:
        """Angle between the projection of `self` on the (xy) plane, and the x axis.
        In [-π, π] range
        """

        return atan2(self.y, self.x)

    @cached_property
    def delta(self) -> ScalarArray[D.Angle]:
        """The complementary angle between `self` and z axis.
        In [-π/2, π/2] range
        """

        xy = sqrt(square(self.x) + square(self.y))
        return atan2(self.z, xy)

    def angle(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> ScalarArray[D.Angle]:
        """Returns the angle between the two given vectors.
        Returned angle is in `[0;π]` range.
        """
        cos_angle = cast(ScalarArray[D.Dimless], self.dot(o) / (self.length * o.length))
        return acos(cos_angle)

    def normalized(self) -> Vector3Array[D.Dimless]:
        """Return vectors scaled to unit norm (dimensionless)."""
        l = self.length
        return self / l

    @staticmethod
    def from_spherical(theta: ScalarArray[D.Angle], delta: ScalarArray[D.Angle], rho: ScalarArray[SomeOtherDim]) -> Vector3Array[SomeOtherDim]:
        """
        Create a vector array from spherical coordinates.

        Uses "radius-longitude-latitude" convention, [see in Wikipedia.](https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques#Convention_rayon-longitude-latitude))

        Arguments
        ---------
        - `theta:` Longitude angles (θ), as ``ScalarArray[D.Angle]``.
        - `delta:` Latitude angles (δ), as ``ScalarArray[D.Angle]``.
        - `rho:` Radius (ρ) from given convention.
        """

        cos_delta = cos(delta)
        sin_delta = sin(delta)
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        x = rho * cos_theta * cos_delta
        y = rho * sin_theta * cos_delta
        z = rho * sin_delta

        return Vector3Array[rho.dim](
            np.stack([x.base_unit_value, y.base_unit_value, z.base_unit_value])
        )

    @overload
    def dot(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> ScalarArray[ProductDim[SomeDim, SomeDim]]: ...

    @overload
    def dot(self: Vector3Array[SomeDim], o: Vector3[SomeOtherDim] | Vector3Array[SomeOtherDim]) -> ScalarArray[ProductDim[SomeDim, SomeOtherDim]]: ...

    def dot(self, o: object) -> ScalarArray[Any]:
        """Return dot product with a vector or vector array."""
        if not isinstance(o, (Vector3, Vector3Array)):
            raise TypeError("dot product requires two Vector3/Vector3Array instances")

        return ScalarArray.dimensionalize(
            self.dim_coords * o.dim_coords
        )(
            np.sum(self._values * o._values, axis=0)
        )

    @overload
    def cross(self: Vector3Array[D.Dimless], o: Vector3[D.Dimless] | Vector3Array[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def cross(self: Vector3Array[D.Dimless], o: Vector3[SomeOtherDim] | Vector3Array[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def cross(self: Vector3Array[SomeDim], o: Vector3[SomeOtherDim] | Vector3Array[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...

    def cross(self, o: object) -> Vector3Array[Any]:
        """Return cross product with a vector or vector array."""
        if not isinstance(o, (Vector3, Vector3Array)):
            raise TypeError("cross product requires two Vector3/Vector3Array instances")

        return Vector3Array.dimensionalize(
            self.dim_coords * o.dim_coords,
        )(
            np.cross(self._values, o._values, axis=0)
        )

    def __repr__(self) -> str:
        return (
            f"Vector3Array[D.{self.dim.__class__.__name__}]("
            f"x={self._values[0]}, "
            f"y={self._values[1]}, "
            f"z={self._values[2]}"
            ")"
        )

    ### + OPERATOR ###

    @overload
    def __add__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> Vector3Array[SomeDim]: ...

    @overload
    def __add__(self: Vector3Array[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    def __add__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__add__(o))

    @overload
    def __radd__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> Vector3Array[SomeDim]: ...

    @overload
    def __radd__(self: Vector3Array[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    def __radd__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> Vector3Array[SomeDim]: ...

    @overload
    def __sub__(self: Vector3Array[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    def __sub__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__sub__(o))

    @overload
    def __rsub__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> Vector3Array[SomeDim]: ...

    @overload
    def __rsub__(self: Vector3Array[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    def __rsub__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __mul__(self: Vector3Array[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def __mul__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[SomeDim]: ...

    @overload
    def __mul__(self: Vector3Array[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar | ScalarArray) -> Vector3Array[Any]: ...

    def __mul__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __rmul__(self: Vector3Array[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def __rmul__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[SomeDim]: ...

    @overload
    def __rmul__(self: Vector3Array[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __rmul__(self, o: Scalar | ScalarArray) -> Vector3Array[Any]: ...

    def __rmul__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __truediv__(self: Vector3Array[D.Dimless], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[SomeDim]: ...

    @overload
    def __truediv__(self: Vector3Array[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: Vector3Array[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar | ScalarArray) -> Vector3Array[Any]: ...

    def __truediv__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __rtruediv__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Vector3Array[Any]:
        return cast(Vector3Array[Any], super().__rtruediv__(o))

    ### @ OPERATOR (dot product) ###

    @overload
    def __matmul__(self: Vector3Array[D.Dimless], o: Vector3[D.Dimless] | Vector3Array[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __matmul__(self: Vector3Array[SomeDim], o: Vector3[SomeOtherDim] | Vector3Array[SomeOtherDim]) -> ScalarArray[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> ScalarArray[Any]:
        if isinstance(o, (Vector3, Vector3Array)):
            return self.dot(o)

        raise TypeError("...")

    ################

    def __eq__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> np.ndarray:
        return np.all(super().__eq__(o)._values, axis=0)

    def __neq__(self: Vector3Array[SomeDim], o: Vector3[SomeDim] | Vector3Array[SomeDim]) -> np.ndarray:
        return ~self.__eq__(o)

    def __le__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3Array")

    def __lt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3Array")

    def __ge__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3Array")

    def __gt__(self, o: Never) -> Never:
        raise RuntimeError("Comparisons operations not defined for Vector3Array")