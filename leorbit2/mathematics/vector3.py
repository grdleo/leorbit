from abc import ABC, abstractmethod
from dataclasses import dataclass
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import Dim, DimEls, D, ProductDim, QuotientDim
from leorbit2.mathematics.scalar import Scalar, scalar_class_factory
from leorbit2.mathematics.tensor import Tensor, dimensional_tensor_class_factory, ensure_same_dimensions

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

def vector3_class_factory(dim: DimEls) -> type[Vector3]:
    return cast(
        type[Vector3], 
        dimensional_tensor_class_factory(
            dim,
            Vector3,
        )
    )
    
NumberOrScalarT = TypeVar("NumberOrScalarT", bound=Number | Scalar)

def all_simple_numbers(els: list[object]) -> TypeGuard[list[Number]]:
    return all(isinstance(el, Number) for el in els)

def all_scalars_numbers(els: list[object]) -> TypeGuard[list[Scalar]]:
    return all(isinstance(el, Scalar) for el in els)

class Vector3[SomeDim](Tensor[SomeDim]):
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
        dim_els = getattr(dim, "_d", None)
        if dim_els == self._dimension:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @property
    def x(self) -> Scalar[SomeDim]:
        return scalar_class_factory(self._dimension)(self._values[0])
    
    @property
    def y(self) -> Scalar[SomeDim]:
        return scalar_class_factory(self._dimension)(self._values[1])
    
    @property
    def z(self) -> Scalar[SomeDim]:
        return scalar_class_factory(self._dimension)(self._values[2])
    
    def dot(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]:
        if not isinstance(o, Vector3):
            raise TypeError("dot product requires two Vector3 instances")

        # Dot product: transpose first vector and matrix multiply
        result_values = self._values.transpose().dot(o._values)
        return scalar_class_factory(self._dimension)(result_values)
    
    def cross(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]:
        if not isinstance(o, Vector3):
            raise TypeError("cross product requires two Vector3 instances")

        # Dot product: transpose first vector and matrix multiply
        result_values = np.cross(self._values, o._values, axis=0)
        return vector3_class_factory(self._dimension)(result_values)
    
    def __repr__(self) -> str:
        return (
            f"Vector3[D.{self._dimension.__class__.__name__}]("
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
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[D.Dimless]: ...

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