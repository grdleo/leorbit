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

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

class Scalar(Generic[SomeDim], Tensor[SomeDim]):
    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[Scalar]:
        return cast(
            type[Scalar],
            super().dimensionalize(dim_coords)
        )

    @classmethod
    def new(cls, value: Number | TensorData):
        if isinstance(value, Number):
            pass
        elif isinstance(value, np.ndarray):
            value = value.item()
        
        return cls(value)

    def cast(self, dim: type[SomeOtherDim]) -> Scalar[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @property
    def base_unit_value(self) -> Number:
        return np.float64(self._values)
    
    def magnitude(self, units: str = "1") -> np.float64:
        return np.float64(self.get_raw_array(units))
    
    def __repr__(self) -> str:
        return f"Scalar[D.{self.dim.__class__.__name__}]({self._values})"

    ### + OPERATOR ###

    @overload
    def __add__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __add__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__add__(o))
    
    @overload
    def __radd__(self: Scalar[D.Dimless], o: Number) -> Scalar[D.Dimless]: ...

    @overload
    def __radd__(self: Scalar[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__radd__(o))
    
    ### - OPERATOR ###
    
    @overload
    def __sub__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __sub__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__sub__(o))
    
    @overload
    def __rsub__(self: Scalar[D.Dimless], o: Number) -> Scalar[D.Dimless]: ...

    @overload
    def __rsub__(self: Scalar[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __mul__(self: Scalar[D.Dimless], o: Scalar[SomeOtherDim]) -> Scalar[SomeOtherDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Number | Scalar[D.Dimless]) -> Scalar[SomeDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Scalar[Any]: ...

    def __mul__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Scalar[D.Dimless], o: Number) -> Scalar[D.Dimless]: ...

    @overload
    def __rmul__(self: Scalar[SomeDim], o: Number) -> Scalar[SomeDim]: ...

    def __rmul__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __truediv__(self: Scalar[D.Dimless], o: Scalar[SomeDim]) -> Scalar[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Number | Scalar[D.Dimless]) -> Scalar[SomeDim]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[D.Dimless]: ... # type: ignore

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Scalar[Any]: ...

    def __truediv__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Scalar[D.Dimless], o: Number) -> Scalar[D.Dimless]: ...

    @overload
    def __rtruediv__(self: Scalar[SomeDim], o: Number) -> Scalar[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rtruediv__(o))
    
    ### % OPERATOR ###

    @overload
    def __mod__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __mod__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __mod__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __mod__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__mod__(o))
    
    @overload
    def __rmod__(self: Scalar[D.Dimless], o: Number) -> Scalar[D.Dimless]: ...

    @overload
    def __rmod__(self: Scalar[SomeDim], o: Number) -> Never: ...

    def __rmod__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__radd__(o))
    
    ### @ OPERATOR ###

    def __matmul__(self, o: Never) -> Never:
        raise RuntimeError("@ operation not defined for scalar")
    
    ### COMPARISONS

    def __eq__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool:
        return bool(super().__eq__(o)._values)
    
    def __neq__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool:
        return bool(super().__eq__(o)._values)

    @overload
    def __lt__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool: ...

    @overload
    def __lt__(self: Scalar[D.Dimless], o: Scalar[D.Dimless] | Number) -> bool: ...

    @overload
    def __lt__(self: Scalar[SomeDim], o: Literal[0]) -> bool: ...

    def __lt__(self, o: object) -> bool:
        return bool(super().__lt__(o)._values)

    @overload
    def __le__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool: ...

    @overload
    def __le__(self: Scalar[D.Dimless], o: Scalar[D.Dimless] | Number) -> bool: ...

    @overload
    def __le__(self: Scalar[SomeDim], o: Literal[0]) -> bool: ...

    def __le__(self, o: object) -> bool:
        return bool(super().__le__(o)._values)

    @overload
    def __gt__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool: ...

    @overload
    def __gt__(self: Scalar[D.Dimless], o: Scalar[D.Dimless] | Number) -> bool: ...

    @overload
    def __gt__(self: Scalar[SomeDim], o: Literal[0]) -> bool: ...

    def __gt__(self, o: object) -> bool:
        return bool(super().__gt__(o)._values)

    @overload
    def __ge__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> bool: ...

    @overload
    def __ge__(self: Scalar[D.Dimless], o: Scalar[D.Dimless] | Number) -> bool: ...

    @overload
    def __ge__(self: Scalar[SomeDim], o: Literal[0]) -> bool: ...

    def __ge__(self, o: object) -> bool:
        return bool(super().__ge__(o)._values)
    
    #############################################

    # ** OPERATOR

    @overload
    def __pow__(self: Scalar[D.Dimless], o: Fraction | Number) -> Scalar[D.Dimless]: ...

    @overload
    def __pow__(self: Scalar[D.Dimless], o: Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

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