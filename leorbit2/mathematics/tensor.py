from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import singledispatchmethod
from pyclbr import Class
from typing import Any, Callable, ClassVar, Generic, Literal, Never, Self, TypeAlias, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np
import numpy.typing as npt

from leorbit2.mathematics.dimensions import Dim, DimCoords, D, registered_units

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
TensorData = np.typing.NDArray[np.floating[Any]]
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

TensorDataTransformer: TypeAlias = Callable[[TensorData], TensorData]

class Tensor[SomeDim = D.Dimless]():
    _dim: type[Dim]
    _base_tensor_class: type[Tensor]

    def __init__(self, values: Number | TensorData):
        """CAREFUL!! User should never instantiante any tensor using `__init__`.
        Always use the classmethod `new`."""
        if not hasattr(self.__class__, "_dim"):
            raise RuntimeError("...")
        
        self._values = np.array(values)

    @property
    def dim_coords(self) -> DimCoords:
        return self._dim._d
    
    @property
    def dim(self) -> type[SomeDim]:
        return cast(
            type[SomeDim],
            self._dim
        )
    
    @classmethod
    def _is_base_tensor_class(cls) -> bool:
        return not hasattr(cls, "_base_tensor_class")

    def __repr__(self) -> str:
        return f"Tensor[D.{self.dim.__class__.__name__}]({self._values})"

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        if not issubclass(dim, Dim):
            raise TypeError("...")

        return _dimensional_tensor_class_factory(
            dim,
            cast(type[Tensor], cls)
        )
    
    def cast(self, dim: type[SomeOtherDim]) -> Tensor[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    def copy(self) -> Self:
        return self.__class__(self._values)

    def ensure_compatible_dimensions(self: Tensor[SomeDim], o: Tensor[SomeOtherDim]) -> TypeIs[Tensor[SomeOtherDim]]:
        return (
            isinstance(o, Tensor) 
            and o.dim_coords == self.dim_coords
        )
    
    def check(self, dim: type[Dim]) -> bool:
        """Returns `True` if tensor is of dimension `dim`"""
        return self.dim_coords == dim._d
    
    def get_raw_array(self, units: str = "1") -> npt.NDArray[np.float64]:
        a = np.copy(self._values)
        if units == "1":
            return a
        
        try:
            dim, factor = registered_units()[units]
        except KeyError:
            raise ValueError(f"No unit named '{units}'")

        if dim._d != self.dim_coords:
            raise ValueError("Units ... has dimension ... which is incompatible with tensor dimension ...")
        
        return a / factor
    
    @singledispatchmethod
    def transform(self, dim, function) -> Tensor[Any]:
        ...

    @transform.register
    def _(self, dim: DimCoords, function: TensorDataTransformer) -> Tensor[Any]:
        d = D.get_dimension_from_coords(dim)

        return _dimensional_tensor_class_factory(
            d,
            cast(type[Tensor], self.__class__)
        )(
            function(self._values)
        )

    @transform.register
    def _(self, dim: type[Dim], function: TensorDataTransformer) -> Tensor[Any]:
        return _dimensional_tensor_class_factory(
            dim,
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
            return self.transform(
                self.dim_coords * o.dim_coords,
                lambda values: values * o._values
            )
        except:
            raise RuntimeError("...")

    def __rmul__(self, o: object) -> Tensor: # o * self
        o = ensure_tensor(o)
        try:
            return self.transform(
                o.dim_coords * self.dim_coords,
                lambda values: o._values * values
            )
        except:
            raise RuntimeError("...")

    def __truediv__(self, o: object) -> Tensor: # self / o
        o = ensure_tensor(o)
        try:
            return self.transform(
                self.dim_coords / o.dim_coords,
                lambda values: values / o._values
            )
        except:
            raise RuntimeError("...")

    def __rtruediv__(self, o: object) -> Tensor: # o / self
        o = ensure_tensor(o)
        try:
            return self.transform(
                o.dim_coords / self.dim_coords,
                lambda values: o._values / values
            )
        except:
            raise RuntimeError("...")

    def __matmul__(self, o: object) -> Tensor: # self @ o
        o = ensure_tensor(o)
        try:
            return self.transform(
                self.dim_coords * o.dim_coords,
                lambda values: values @ o._values
            )
        except:
            raise RuntimeError("...")

    def __rmatmul__(self, o: object) -> Tensor: # o @ self
        o = ensure_tensor(o)
        try:
            return self.transform(
                o.dim_coords * self.dim_coords,
                lambda values: o._values @ values
            )
        except:
            raise RuntimeError("...")
        
    def __mod__(self, o: object) -> Tensor: # self % o
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        
        try:
            return cast(type[Tensor], self.__class__)(self._values % o._values)
        except:
            raise RuntimeError("...")

    def __rmod__(self, o: object) -> Tensor: # o % self
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
         raise RuntimeError("Tensors dimensions are not compatible!")

        try:
            return cast(type[Tensor], self.__class__)(o._values % self._values)
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
                D.Dimless,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values == o._values)
                )
            )
        except:
            raise RuntimeError("...")
        
    def __neq__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values != o._values)
                )
            )
        except:
            raise RuntimeError("...")
    
    def __lt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values < o._values)
                )
            )
        except:
            raise RuntimeError("...")
    
    def __le__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values <= o._values)
                )
            )
        except:
            raise RuntimeError("...")
    
    def __gt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values > o._values)
                )
            )
        except:
            raise RuntimeError("...")
    
    def __ge__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self.transform(
                D.Dimless._d,
                cast(
                    TensorDataTransformer, 
                    lambda values: np.int8(values >= o._values)
                )
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
    if all(t.dim_coords == t0.dim_coords for t in tensors):
        return True
    
    raise RuntimeError("...")

def _dimensional_tensor_class_factory(
    dim: type[Dim], 
    parent_class: type[Tensor] = Tensor
) -> type[Tensor]:
    base_tensor_class = parent_class if parent_class._is_base_tensor_class() else parent_class._base_tensor_class

    return type(
        "TensorWithDimension",
        (base_tensor_class, ),
        dict(
            _dimension=dim,
            _base_tensor_class=base_tensor_class
        )
    )