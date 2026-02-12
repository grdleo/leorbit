from abc import ABC, abstractmethod
from dataclasses import dataclass
from pyclbr import Class
from typing import Any, Callable, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

from leorbit2.mathematics.dimensions import Dim, DimEls, D

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

class Tensor[SomeDim = D.Dimless]():
    _dimension: DimEls

    def __init__(self, values: Number | TensorData):
        """CAREFUL!! User should never instantiante any tensor using `__init__`.
        Always use the classmethod `new`."""
        if not hasattr(self.__class__, "_dimension"):
            raise RuntimeError("!!")
        
        self._values = np.array(values)

    @property
    def dim(self) -> DimEls:
        return self._dimension

    def __repr__(self) -> str:
        return f"Tensor[D.{self._dimension.__class__.__name__}]({self._values})"

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        dim_els: DimEls = getattr(dim, "_d")

        return dimensional_tensor_class_factory(
            dim_els,
            cast(type[Tensor], cls)
        )
    
    def cast(self, dim: type[SomeOtherDim]) -> Tensor[SomeOtherDim]:
        if dim._d == self._dimension:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    def copy(self) -> Self:
        return self.__class__(self._values)

    def ensure_compatible_dimensions(self, o: Tensor) -> TypeIs[Tensor[SomeDim]]:
        return isinstance(o, Tensor) and o._dimension == self._dimension
    
    def check(self, dim: type[Dim]) -> bool:
        """Returns `True` if tensor is of dimension `dim`"""
        return self._dimension == dim._d
    
    def transform(self, new_dim: DimEls, function: Callable[[TensorData], TensorData]) -> Tensor[Any]: # type: ignore
        return dimensional_tensor_class_factory(
            new_dim,
            cast(type[Tensor], self.__class__)
        )(
            function(self._values)
        )
    
    def __add__(self, o: object) -> Tensor[SomeDim]: # self + o
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return self.__class__(self._values + o._values)
        except Exception as e:
            raise RuntimeError("...") from e
    
    def __radd__(self, o: object) -> Tensor[SomeDim]: # o + self
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return self.__class__(o._values + self._values)
        except:
            raise RuntimeError("...")
    
    def __sub__(self, o: object) -> Tensor[SomeDim]: # self - o
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return self.__class__(self._values - o._values)
        except:
            raise RuntimeError("...")

    def __rsub__(self, o: object) -> Tensor[SomeDim]: # o - self
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return self.__class__(o._values - self._values)
        except:
            raise RuntimeError("...")

    def __mul__(self, o: object) -> Tensor: # self * o
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                self._dimension * o._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(self._values * o._values)
        except:
            raise RuntimeError("...")

    def __rmul__(self, o: object) -> Tensor: # o * self
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                o._dimension * self._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(o._values * self._values)
        except:
            raise RuntimeError("...")

    def __truediv__(self, o: object) -> Tensor: # self / o
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                self._dimension / o._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(self._values / o._values)
        except:
            raise RuntimeError("...")

    def __rtruediv__(self, o: object) -> Tensor: # o / self
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                o._dimension / self._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(o._values / self._values)
        except:
            raise RuntimeError("...")

    def __matmul__(self, o: object) -> Tensor: # self @ o
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                self._dimension * o._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(self._values @ o._values)
        except:
            raise RuntimeError("...")

    def __rmatmul__(self, o: object) -> Tensor: # o @ self
        o = ensure_tensor(o)
        try:
            tensor_class = dimensional_tensor_class_factory(
                o._dimension * self._dimension, 
                cast(type[Tensor], self.__class__)
            )
            return tensor_class(o._values @ self._values)
        except:
            raise RuntimeError("...")
        
    def __mod__(self, o: object) -> Tensor: # self % o
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return self.__class__(self._values % o._values)
        except:
            raise RuntimeError("...")

    def __rmod__(self, o: object) -> Tensor: # o % self
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
         raise RuntimeError("Tensors dimensions are not compatible!")

        try:
            return self.__class__(o._values % self._values)
        except:
            raise RuntimeError("...")
        
    def __pos__(self) -> Self:
        return self.__class__(+self._values)
    
    def __neg__(self) -> Self:
        return self.__class__(-self._values)
    
    def __eq__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                np.equal
            )
        except:
            raise RuntimeError("...")
        
    def __neq__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                np.not_equal
            )
        except:
            raise RuntimeError("...")
    
    def __lt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                lambda data: data < o._values
            )
        except:
            raise RuntimeError("...")
    
    def __le__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                lambda data: data <= o._values
            )
        except:
            raise RuntimeError("...")
    
    def __gt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                lambda data: data > o._values
            )
        except:
            raise RuntimeError("...")
    
    def __ge__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                lambda data: data >= o._values
            )
        except:
            raise RuntimeError("...")
    
def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[SomeDim] | Tensor[D.Dimless]:
    """ensure tensor. if not a tensor object, creates a dimless tensor"""

    if isinstance(o, Tensor):
        return o
    elif isinstance(o, Number | np.ndarray):
        return Tensor[D.Dimless](o)

    raise RuntimeError("...")

def ensure_same_dimensions(*tensors: Tensor[Any]) -> Literal[True]:
    if len(tensors) == 0:
        return True
    
    t0, *_ = tensors
    if all(t._dimension == t0._dimension for t in tensors):
        return True
    
    raise RuntimeError("...")

def dimensional_tensor_class_factory(
    dim: DimEls, 
    parent_class: type[Tensor] = Tensor
) -> type[Tensor]:
    # FIXME
    return type(
        "DimTensor",
        (parent_class, ),
        dict(_dimension=dim)
    )