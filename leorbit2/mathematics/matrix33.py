from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import Dim, D, DimCoords, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.scalar import Scalar, TensorScalar
from leorbit2.mathematics.tensor import Tensor, ensure_same_dimensions
from leorbit2.mathematics.vector3 import Vector3, Vector3Array, TensorVector3, all_scalars_numbers, all_simple_numbers

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)
NumberOrScalarT = TypeVar("NumberOrScalarT", bound=Number | Scalar)
    
class Matrix33(Tensor[SomeDim], Generic[SomeDim]):
    """3×3 matrix carrying a physical dimension."""

    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[Matrix33]:
        """Return a matrix class bound to ``dim_coords``."""
        return cast(
            type[Matrix33],
            super().dimensionalize(dim_coords)
        )
    
    @classmethod
    def new(cls,
        a: NumberOrScalarT, b: NumberOrScalarT, c: NumberOrScalarT,
        d: NumberOrScalarT, e: NumberOrScalarT, f: NumberOrScalarT,
        g: NumberOrScalarT, h: NumberOrScalarT, i: NumberOrScalarT,
    ):
        """Create a matrix from row-major coefficients."""
        mat: TensorData
        mat_els = cast(list[object], [a, b, c, d, e, f, g, h, i])

        if all_simple_numbers(mat_els):
            mat = np.array(mat_els).reshape((3,3))
        elif all_scalars_numbers(mat_els):
            ensure_same_dimensions(*mat_els)
            mat = np.array([el.base_unit_value for el in mat_els]).reshape((3,3))
        else:
            raise RuntimeError("...")

        return cls(mat)
    
    def cast(self, dim: type[SomeOtherDim]) -> Matrix33[SomeOtherDim]:
        """Type-cast to another dimension when coordinates are identical."""
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @cached_property
    def det(self) -> Number:
        return np.linalg.det(self._values)
    
    @overload
    def inverse(self: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def inverse(self: Matrix33[SomeDim]) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...
    
    def inverse(self) -> Matrix33:
        """Return matrix inverse.
        """
        return Matrix33.dimensionalize(
            self.dim_coords ** -1
        )(
            np.linalg.inv(self._values)
        )
    
    def __repr__(self) -> str:
        return f"Matrix33[D.{self.dim.__class__.__name__}]({self._values})"

    ### + OPERATOR ###

    @overload
    def __add__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __add__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __add__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __add__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__add__(o))

    @overload
    def __radd__(self: Matrix33[D.Dimless], o: Number) -> Matrix33[D.Dimless]: ...

    @overload
    def __radd__(self: Matrix33[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __sub__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __sub__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __sub__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__sub__(o))

    @overload
    def __rsub__(self: Matrix33[D.Dimless], o: Number) -> Matrix33[D.Dimless]: ...

    @overload
    def __rsub__(self: Matrix33[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rsub__(o))

    ### * OPERATOR (element-wise scaling) ###

    @overload
    def __mul__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __mul__(self: Matrix33[D.Dimless], o: Scalar[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[SomeDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Matrix33[Any]: ...

    def __mul__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Matrix33[D.Dimless], o: Number) -> Matrix33[D.Dimless]: ...

    @overload
    def __rmul__(self: Matrix33[SomeDim], o: Number) -> Matrix33[SomeDim]: ...

    def __rmul__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __truediv__(self: Matrix33[D.Dimless], o: Scalar[SomeDim]) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[SomeDim]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Matrix33[Any]: ...

    def __truediv__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Matrix33[D.Dimless], o: Number) -> Matrix33[D.Dimless]: ...

    @overload
    def __rtruediv__(self: Matrix33[SomeDim], o: Number) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rtruediv__(o))

    ### @ OPERATOR (matrix multiplication) ###

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Matrix33[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Matrix33[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Vector3[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3Array[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3Array[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Vector3Array[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> Matrix33[Any] | Vector3[Any] | Vector3Array[Any]:
        """Apply matrix product with matrix, vector, or vector array."""
        if isinstance(o, Matrix33):
            return Matrix33.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )
        elif isinstance(o, Vector3):
            return Vector3.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )
        elif isinstance(o, Vector3Array):
            return Vector3Array.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )

        raise TypeError("@ operation requires Matrix33 or Vector3 or Vector3Array")


TensorMatrix33 = Matrix33