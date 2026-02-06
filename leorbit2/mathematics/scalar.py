from abc import ABC, abstractmethod
from dataclasses import dataclass
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import Dim, DimEls, D, ProductDim, QuotientDim
from leorbit2.mathematics.tensor import Tensor, dimensional_tensor_class_factory

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

def scalar_class_factory(dim: DimEls) -> type[Scalar]:
    return cast(
        type[Scalar], 
        dimensional_tensor_class_factory(
            dim,
            Scalar,
        )
    )

class Scalar[SomeDim](Tensor[SomeDim]):
    @classmethod
    def new(cls, value: Number | TensorData):
        if isinstance(value, Number):
            pass
        elif isinstance(value, np.ndarray):
            value = value.item()
        
        return cls(value)

    def cast(self, dim: type[SomeOtherDim]) -> Scalar[SomeOtherDim]:
        dim_els = getattr(dim, "_d", None) 
        if dim_els == self._dimension:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @property
    def base_unit_value(self) -> Number:
        return np.floating(self._values)
    
    def magnitude(self, units: str = "1") -> Number:
        raise NotImplementedError()
    
    def __repr__(self) -> str:
        return f"Scalar[D.{self._dimension.__class__.__name__}]({self.base_unit_value})"

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
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[D.Dimless]: ...

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
    
    ### @ OPERATOR ###

    def __matmul__(self, o: object) -> Scalar[Any]:
        raise RuntimeError("@ operation not defined for scalar")