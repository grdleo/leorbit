from functools import cached_property
from typing import Any, Generic, Never, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import D, DimCoords, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.scalar import TensorScalar
from leorbit2.mathematics.tensor import Tensor, ensure_same_dimensions
from leorbit2.mathematics.vector3 import TensorVector3, all_scalars_numbers, all_simple_numbers

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
NumberOrTensorScalarT = TypeVar("NumberOrTensorScalarT", bound=Number | TensorScalar[Any])


class TensorMatrix33(Tensor[SomeDim], Generic[SomeDim]):
    """3×3 matrix carrying a physical dimension."""

    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[TensorMatrix33]:
        """Return a matrix class bound to ``dim_coords``."""
        return cast(
            type[TensorMatrix33],
            super().dimensionalize(dim_coords)
        )
    
    @classmethod
    def new(cls,
        a: NumberOrTensorScalarT, b: NumberOrTensorScalarT, c: NumberOrTensorScalarT,
        d: NumberOrTensorScalarT, e: NumberOrTensorScalarT, f: NumberOrTensorScalarT,
        g: NumberOrTensorScalarT, h: NumberOrTensorScalarT, i: NumberOrTensorScalarT,
    ) -> TensorMatrix33[SomeDim]:
        """Create a matrix from row-major coefficients."""
        mat: TensorData
        mat_els = cast(list[object], [a, b, c, d, e, f, g, h, i])

        if all_simple_numbers(mat_els):
            mat = np.array(mat_els).reshape((3,3))
        elif all_scalars_numbers(mat_els):
            scalars = cast(list[TensorScalar[Any]], mat_els)
            ensure_same_dimensions(*scalars)
            mat = np.array([el.base_unit_value for el in scalars]).reshape((3,3))
        else:
            raise RuntimeError("...")

        return cls(mat)
    
    def cast(self, dim: type[SomeOtherDim]) -> TensorMatrix33[SomeOtherDim]:
        """Type-cast to another dimension when coordinates are identical."""
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @cached_property
    def det(self) -> Number:
        return np.linalg.det(self._values)
    
    @overload
    def inverse(self: TensorMatrix33[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def inverse(self: TensorMatrix33[SomeDim]) -> TensorMatrix33[QuotientDim[D.Dimless, SomeDim]]: ...
    
    def inverse(self) -> TensorMatrix33:
        """Return matrix inverse.
        """
        return TensorMatrix33.dimensionalize(
            self.dim_coords ** -1
        )(
            np.linalg.inv(self._values)
        )
    
    def __repr__(self) -> str:
        return f"TensorMatrix33[D.{self.dim.__class__.__name__}]({self._values})"

    ### + OPERATOR ###

    @overload
    def __add__(self: TensorMatrix33[D.Dimless], o: TensorMatrix33[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __add__(self: TensorMatrix33[SomeDim], o: TensorMatrix33[SomeDim]) -> TensorMatrix33[SomeDim]: ...

    @overload
    def __add__(self: TensorMatrix33[SomeDim], o: TensorScalar[SomeDim]) -> TensorMatrix33[SomeDim]: ...

    def __add__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__add__(o))

    @overload
    def __radd__(self: TensorMatrix33[D.Dimless], o: Number) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __radd__(self: TensorMatrix33[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: TensorMatrix33[D.Dimless], o: TensorMatrix33[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __sub__(self: TensorMatrix33[SomeDim], o: TensorMatrix33[SomeDim]) -> TensorMatrix33[SomeDim]: ...

    @overload
    def __sub__(self: TensorMatrix33[SomeDim], o: TensorScalar[SomeDim]) -> TensorMatrix33[SomeDim]: ...

    def __sub__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__sub__(o))

    @overload
    def __rsub__(self: TensorMatrix33[D.Dimless], o: Number) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __rsub__(self: TensorMatrix33[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__rsub__(o))

    ### * OPERATOR (element-wise scaling) ###

    @overload
    def __mul__(self: TensorMatrix33[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __mul__(self: TensorMatrix33[D.Dimless], o: TensorScalar[SomeOtherDim]) -> TensorMatrix33[SomeOtherDim]: ...

    @overload
    def __mul__(self: TensorMatrix33[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorMatrix33[SomeDim]: ...

    @overload
    def __mul__(self: TensorMatrix33[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorMatrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: TensorScalar) -> TensorMatrix33[Any]: ...

    def __mul__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__mul__(o))

    @overload
    def __rmul__(self: TensorMatrix33[D.Dimless], o: Number) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __rmul__(self: TensorMatrix33[SomeDim], o: Number) -> TensorMatrix33[SomeDim]: ...

    def __rmul__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: TensorMatrix33[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __truediv__(self: TensorMatrix33[D.Dimless], o: TensorScalar[SomeDim]) -> TensorMatrix33[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: TensorMatrix33[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorMatrix33[SomeDim]: ...

    @overload
    def __truediv__(self: TensorMatrix33[SomeDim], o: TensorScalar[SomeDim]) -> TensorMatrix33[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: TensorMatrix33[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorMatrix33[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: TensorScalar) -> TensorMatrix33[Any]: ...

    def __truediv__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: TensorMatrix33[D.Dimless], o: Number) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __rtruediv__(self: TensorMatrix33[SomeDim], o: Number) -> TensorMatrix33[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> TensorMatrix33[Any]:
        return cast(TensorMatrix33[Any], super().__rtruediv__(o))

    ### @ OPERATOR (matrix multiplication) ###

    @overload
    def __matmul__(self: TensorMatrix33[D.Dimless], o: TensorMatrix33[D.Dimless]) -> TensorMatrix33[D.Dimless]: ...

    @overload
    def __matmul__(self: TensorMatrix33[D.Dimless], o: TensorMatrix33[SomeOtherDim]) -> TensorMatrix33[SomeOtherDim]: ...

    @overload
    def __matmul__(self: TensorMatrix33[SomeDim], o: TensorMatrix33[SomeOtherDim]) -> TensorMatrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: TensorMatrix33[D.Dimless], o: TensorVector3[D.Dimless]) -> TensorVector3[D.Dimless]: ...

    @overload
    def __matmul__(self: TensorMatrix33[D.Dimless], o: TensorVector3[SomeOtherDim]) -> TensorVector3[SomeOtherDim]: ...

    @overload
    def __matmul__(self: TensorMatrix33[SomeDim], o: TensorVector3[SomeOtherDim]) -> TensorVector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> TensorMatrix33[Any] | TensorVector3[Any]:
        """Apply matrix product with matrix, vector, or vector array."""
        if isinstance(o, TensorMatrix33):
            return TensorMatrix33.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )
        elif isinstance(o, TensorVector3):
            return TensorVector3.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )
        raise TypeError("@ operation requires TensorMatrix33 or TensorVector3")


class Matrix33(TensorMatrix33[SomeDim], Generic[SomeDim]):
    """End-user convenience matrix class."""