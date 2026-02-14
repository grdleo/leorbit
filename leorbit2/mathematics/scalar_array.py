from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, IntEnum
from fractions import Fraction
from math import acos, asin, atan, cos, sin, tan
from pyclbr import Class
from typing import Annotated, Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np
from numpy._typing import _UFunc_Nin1_Nout1

from leorbit2.mathematics.dimensions import Dim, DimCoords, D, PowerDim, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.tensor import Tensor
from leorbit2.mathematics.scalar import Scalar

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

class ScalarArray(Generic[SomeDim], Tensor[SomeDim]):
    @classmethod
    def new(cls, values: list[Number] | TensorData):
        return cls(
            np.array(values).flatten()
        )

    def cast(self, dim: type[SomeOtherDim]) -> ScalarArray[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @property
    def base_unit_value(self) -> np.ndarray:
        return cast(
            np.ndarray,
            np.float64(self._values)
        )
    
    @property
    def size(self) -> int:
        s, = self._values.shape
        return s
    
    def __getitem__(self, index: int) -> Scalar[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("...")
        
        return Scalar[self.dim].new(
            self._values[index]
        )
    
    def magnitude(self, units: str = "1") -> list[Number]:
        raise NotImplementedError()
    
    def __repr__(self) -> str:
        return f"Scalar[D.{self.dim.__class__.__name__}]({self._values})"

    ### + OPERATOR ###

    @overload
    def __add__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __add__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    def __add__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__add__(o))

    @overload
    def __radd__(self: ScalarArray[D.Dimless], o: Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __radd__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __sub__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    def __sub__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__sub__(o))

    @overload
    def __rsub__(self: ScalarArray[D.Dimless], o: Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __rsub__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __mul__(self: ScalarArray[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[SomeOtherDim]: ...

    @overload
    def __mul__(self: ScalarArray[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[SomeDim]: ...

    @overload
    def __mul__(self: ScalarArray[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar | ScalarArray) -> ScalarArray[Any]: ...

    def __mul__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__mul__(o))

    @overload
    def __rmul__(self: ScalarArray[D.Dimless], o: Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __rmul__(self: ScalarArray[SomeDim], o: Number) -> ScalarArray[SomeDim]: ...

    def __rmul__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __truediv__(self: ScalarArray[D.Dimless], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: ScalarArray[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[SomeDim]: ...

    @overload
    def __truediv__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: ScalarArray[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar | ScalarArray) -> ScalarArray[Any]: ...

    def __truediv__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: ScalarArray[D.Dimless], o: Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __rtruediv__(self: ScalarArray[SomeDim], o: Number) -> ScalarArray[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__rtruediv__(o))

    ### % OPERATOR ###

    @overload
    def __mod__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __mod__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    @overload
    def __mod__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    def __mod__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__mod__(o))

    @overload
    def __rmod__(self: ScalarArray[D.Dimless], o: Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __rmod__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    def __rmod__(self, o: object) -> ScalarArray[Any]:
        return cast(ScalarArray[Any], super().__rmod__(o))

    ### @ OPERATOR ###

    def __matmul__(self, o: Never) -> Never:
        raise RuntimeError("@ operation not defined for scalar array")

    ### COMPARISONS

    @overload
    def __eq__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __eq__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> np.ndarray: ...

    def __eq__(self, o: object) -> np.ndarray:
        return super().__eq__(o)._values

    @overload
    def __neq__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __neq__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> np.ndarray: ...

    def __neq__(self, o: object) -> np.ndarray:
        return super().__neq__(o)._values

    @overload
    def __lt__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __lt__(self: ScalarArray[D.Dimless], o: Scalar[D.Dimless] | ScalarArray[D.Dimless] | Number) -> np.ndarray: ...

    @overload
    def __lt__(self: ScalarArray[SomeDim], o: Literal[0]) -> np.ndarray: ...

    def __lt__(self, o: object) -> np.ndarray:
        return super().__lt__(o)._values

    @overload
    def __le__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __le__(self: ScalarArray[D.Dimless], o: Scalar[D.Dimless] | ScalarArray[D.Dimless] | Number) -> np.ndarray: ...

    @overload
    def __le__(self: ScalarArray[SomeDim], o: Literal[0]) -> np.ndarray: ...

    def __le__(self, o: object) -> np.ndarray:
        return super().__le__(o)._values

    @overload
    def __gt__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __gt__(self: ScalarArray[D.Dimless], o: Scalar[D.Dimless] | ScalarArray[D.Dimless] | Number) -> np.ndarray: ...

    @overload
    def __gt__(self: ScalarArray[SomeDim], o: Literal[0]) -> np.ndarray: ...

    def __gt__(self, o: object) -> np.ndarray:
        return super().__gt__(o)._values

    @overload
    def __ge__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> np.ndarray: ...

    @overload
    def __ge__(self: ScalarArray[D.Dimless], o: Scalar[D.Dimless] | ScalarArray[D.Dimless] | Number) -> np.ndarray: ...

    @overload
    def __ge__(self: ScalarArray[SomeDim], o: Literal[0]) -> np.ndarray: ...

    def __ge__(self, o: object) -> np.ndarray:
        return super().__ge__(o)._values

    #############################################

    # ** OPERATOR

    @overload
    def __pow__(self: ScalarArray[D.Dimless], o: Fraction | Number) -> ScalarArray[D.Dimless]: ...

    @overload
    def __pow__(self: ScalarArray[D.Dimless], o: Scalar[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    def __pow__(self, o: object) -> Tensor[Any]:
        if isinstance(o, Scalar) and o.check(D.Dimless):
            o = o.base_unit_value

        if not isinstance(o, Fraction | Number):
            raise ValueError()

        return self.transform(
            self.dim_coords ** Fraction(o),
            lambda data: np.power(data, float(o))
        )

    #############################################