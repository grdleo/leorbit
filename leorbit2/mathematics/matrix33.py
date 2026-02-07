from abc import ABC, abstractmethod
from dataclasses import dataclass
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import Dim, D, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.tensor import Tensor, ensure_same_dimensions
from leorbit2.mathematics.vector3 import Vector3, all_scalars_numbers, all_simple_numbers

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)
NumberOrScalarT = TypeVar("NumberOrScalarT", bound=Number | Scalar)
    
class Matrix33(Generic[SomeDim], Tensor[SomeDim]):
    @classmethod
    def new(cls,
        a: NumberOrScalarT, b: NumberOrScalarT, c: NumberOrScalarT,
        d: NumberOrScalarT, e: NumberOrScalarT, f: NumberOrScalarT,
        g: NumberOrScalarT, h: NumberOrScalarT, i: NumberOrScalarT,
    ):
        """order: by lines"""
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
        if dim._d == self._dimension:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @overload
    def inverse(self: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def inverse(self: Matrix33[SomeDim]) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...
    
    def inverse(self) -> Matrix33:
        raise NotImplementedError()
    
    def __repr__(self) -> str:
        return f"Matrix33[D.{self._dimension.__class__.__name__}]({self._values})"

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
    def __radd__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __radd__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __radd__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

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
    def __rsub__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __rsub__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __rsub__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

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
    def __rmul__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __rmul__(self: Matrix33[D.Dimless], o: Scalar[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __rmul__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[SomeDim]: ...

    @overload
    def __rmul__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __rmul__(self, o: Scalar) -> Matrix33[Any]: ...

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
    def __rtruediv__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __rtruediv__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...

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

    def __matmul__(self, o: object) -> Matrix33[Any] | Vector3[Any]:
        if isinstance(o, Matrix33):
            return cast(Matrix33[Any], super().__matmul__(o))
        elif isinstance(o, Vector3):
            return cast(Vector3[Any], super().__matmul__(o))

        raise TypeError("@ operation requires Matrix33 or Vector3")