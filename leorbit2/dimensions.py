from dataclasses import dataclass
from typing import Any, Never, Self, TypeIs, TypeVar, cast, overload

import numpy as np

@dataclass(frozen=True)
class DimEls:
    length: int = 0
    time: int = 0

    @property
    def dimensionless(self) -> bool:
        return (
            self.length
            == self.time
            == 0
        )
    
    def __mul__(self, o: DimEls) -> DimEls:
        return DimEls(
            length=self.length + o.length,
            time=self.time + o.time
        )
    
    def __truediv__(self, o: DimEls) -> DimEls:
        return DimEls(
            length=self.length - o.length,
            time=self.time - o.time
        )

class Dim:
    d: DimEls | None = None

class DimLess(Dim):
    d = DimEls()

class LengthDim(Dim):
    d = DimEls(length=1)

class TimeDim(Dim):
    d = DimEls(time=1)

class VelocityDim(Dim):
    d = DimEls(length=1, time=-1)

# Registry of known dimensions for lookup
_DIMENSION_REGISTRY: dict[DimEls, type[Dim]] = {
    DimEls(): DimLess,
    DimEls(length=1): LengthDim,
    DimEls(time=1): TimeDim,
    DimEls(length=1, time=-1): VelocityDim,
}

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
SomeDimFull = TypeVar("SomeDimFull", bound=LengthDim | TimeDim | VelocityDim)
TensorData = np.typing.NDArray[np.floating[Any]]

class ProductDim[SomeDim, SomeOtherDim](Dim):
    @classmethod
    def __class_getitem__(cls, params):
        # Handle non-tuple params (e.g., type variables for generic purposes)
        if not isinstance(params, tuple) or len(params) != 2:
            # Return a generic class for type-checking purposes
            return type(f"ProductDim[{params}]", (Dim,), {})

        dim1, dim2 = params
        dim1_els: DimEls | None = getattr(dim1, "d", None)
        dim2_els: DimEls | None = getattr(dim2, "d", None)

        # If dimensions don't have 'd' attribute, return generic class
        if dim1_els is None or dim2_els is None:
            dim1_name = getattr(dim1, '__name__', str(dim1))
            dim2_name = getattr(dim2, '__name__', str(dim2))
            return type(f"ProductDim[{dim1_name}, {dim2_name}]", (Dim,), {})

        # Compute product dimension
        result_els = dim1_els * dim2_els

        # Look up in registry for known dimensions
        result_dim = _DIMENSION_REGISTRY.get(result_els)
        if result_dim is not None:
            return result_dim

        # Create dynamic class for unknown dimension combinations
        return type(
            f"ProductDim[{dim1.__name__}, {dim2.__name__}]",
            (Dim,),
            dict(d=result_els)
        )

class QuotientDim[SomeDim, SomeOtherDim](Dim):
    @classmethod
    def __class_getitem__(cls, params):
        # Handle non-tuple params (e.g., type variables for generic purposes)
        if not isinstance(params, tuple) or len(params) != 2:
            # Return a generic class for type-checking purposes
            return type(f"QuotientDim[{params}]", (Dim,), {})

        dim1, dim2 = params
        dim1_els: DimEls | None = getattr(dim1, "d", None)
        dim2_els: DimEls | None = getattr(dim2, "d", None)

        # If dimensions don't have 'd' attribute, return generic class
        if dim1_els is None or dim2_els is None:
            dim1_name = getattr(dim1, '__name__', str(dim1))
            dim2_name = getattr(dim2, '__name__', str(dim2))
            return type(f"QuotientDim[{dim1_name}, {dim2_name}]", (Dim,), {})

        # Compute quotient dimension
        result_els = dim1_els / dim2_els

        # Look up in registry for known dimensions
        result_dim = _DIMENSION_REGISTRY.get(result_els)
        if result_dim is not None:
            return result_dim

        # Create dynamic class for unknown dimension combinations
        return type(
            f"QuotientDim[{dim1.__name__}, {dim2.__name__}]",
            (Dim,),
            dict(d=result_els)
        )

class Tensor[SomeDim = DimLess]():
    _dimension: DimEls

    def __init__(self, values: Number | TensorData):
        self._values = np.array(values)

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        dim_els: DimEls = getattr(dim, "d")

        return dimensional_tensor_class_factory(
            dim_els,
            cast(type[Tensor], cls)
        )
    
    def ensure_compatible_dimensions(self, o: Tensor) -> TypeIs[Tensor[SomeDim]]:
        return isinstance(o, Tensor) and o._dimension == self._dimension
    
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

def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[SomeDim] | Tensor[DimLess]:
    """ensure tensor. if not a tensor object, creates a dimless tensor"""

    if isinstance(o, Tensor):
        return o
    elif isinstance(o, Number | np.ndarray):
        return Tensor[DimLess](o)

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

class Scalar[SomeDim](Tensor[SomeDim]):
    def __init__(self, value: Number):
        super().__init__(value)

    ### + OPERATOR ###

    @overload
    def __add__(self: Scalar[DimLess], o: Number | Scalar[DimLess]) -> Scalar[DimLess]: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __add__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__add__(o))
    
    @overload
    def __radd__(self: Scalar[DimLess], o: Number) -> Scalar[DimLess]: ...

    @overload
    def __radd__(self: Scalar[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__radd__(o))
    
    ### - OPERATOR ###
    
    @overload
    def __sub__(self: Scalar[DimLess], o: Number | Scalar[DimLess]) -> Scalar[DimLess]: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __sub__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__sub__(o))
    
    @overload
    def __rsub__(self: Scalar[DimLess], o: Number) -> Scalar[DimLess]: ...

    @overload
    def __rsub__(self: Scalar[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: Scalar[DimLess], o: Number | Scalar[DimLess]) -> Scalar[DimLess]: ...

    @overload
    def __mul__(self: Scalar[DimLess], o: Scalar[SomeOtherDim]) -> Scalar[SomeOtherDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Number | Scalar[DimLess]) -> Scalar[SomeDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Scalar[Any]: ...

    def __mul__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Scalar[DimLess], o: Number) -> Scalar[DimLess]: ...

    @overload
    def __rmul__(self: Scalar[SomeDim], o: Number) -> Scalar[SomeDim]: ...

    def __rmul__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Scalar[DimLess], o: Number | Scalar[DimLess]) -> Scalar[DimLess]: ...

    @overload
    def __truediv__(self: Scalar[DimLess], o: Scalar[SomeDim]) -> Scalar[QuotientDim[DimLess, SomeDim]]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Number | Scalar[DimLess]) -> Scalar[SomeDim]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[DimLess]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Scalar[Any]: ...

    def __truediv__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Scalar[DimLess], o: Number) -> Scalar[DimLess]: ...

    @overload
    def __rtruediv__(self: Scalar[SomeDim], o: Number) -> Scalar[QuotientDim[DimLess, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Scalar[Any]:
        return cast(Scalar[Any], super().__rtruediv__(o))
    
    ### @ OPERATOR ###

    def __matmul__(self, o: object) -> Scalar[Any]:
        raise RuntimeError("@ operation not defined for scalar")


a = Scalar[LengthDim](1)
b = Scalar[DimLess](2)
c = Scalar[LengthDim](3)

z = a / c