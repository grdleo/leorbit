from abc import ABC, abstractmethod
from dataclasses import dataclass
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

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
    
class UnitsRegistry:
    """..."""

class Dim:
    _d: ClassVar[DimEls | None] = None

class DimLess(Dim):
    _d = DimEls()

class AngleD(Dim):
    _d = DimEls()

    rad: ClassVar[float] = 1 # base
    deg: ClassVar[float] = np.pi / 180

class AngularAccD(Dim):
    """rad/s**2"""
    _d = DimEls(time=-2)

class AngularJerkD(Dim):
    """rad/s**2"""
    _d = DimEls(time=-3)

class LengthD(Dim):
    """m"""
    _d = DimEls(length=1)

    meter: ClassVar[float] = 1 # base

    milli_meter: ClassVar[float] = 1e-3 * meter
    kilo_meter: ClassVar[float] = 1e3 * meter

class InvLengthD(Dim):
    """1/m"""
    _d = DimEls(length=-1)

class TimeD(Dim):
    _d = DimEls(time=1)

    second: ClassVar[float] = 1 # base

    milli_second: ClassVar[float] = 1e-3 * second
    minute: ClassVar[float] = 60 * second
    hour: ClassVar[float] = 60 * minute
    day: ClassVar[float] = 24 * hour
    month: ClassVar[float] = 30 * day
    year: ClassVar[float] = 365 * day

class VelocityD(Dim):
    _d = DimEls(length=1, time=-1)

    meter_per_second: ClassVar[float] = 1 # base

    kilo_meter_per_hour: ClassVar[float] = meter_per_second / 3.6

# Registry of known dimensions for lookup
_DIMENSION_REGISTRY: dict[DimEls, type[Dim]] = {
    DimEls(): DimLess,
    DimEls(length=1): LengthD,
    DimEls(time=1): TimeD,
    DimEls(length=1, time=-1): VelocityD,
}

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
SomeDimFull = TypeVar("SomeDimFull", bound=LengthD | TimeD | VelocityD)
TensorData = np.typing.NDArray[np.floating[Any]]

class ProductDim[SomeDim, SomeOtherDim](Dim):
    @classmethod
    def __class_getitem__(cls, params):
        # Handle non-tuple params (e.g., type variables for generic purposes)
        if not isinstance(params, tuple) or len(params) != 2:
            # Return a generic class for type-checking purposes
            return type(f"ProductDim[{params}]", (Dim,), {})

        dim1, dim2 = params
        dim1_els: DimEls | None = getattr(dim1, "_d", None)
        dim2_els: DimEls | None = getattr(dim2, "_d", None)

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
        dim1_els: DimEls | None = getattr(dim1, "_d", None)
        dim2_els: DimEls | None = getattr(dim2, "_d", None)

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
    
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

class Tensor[SomeDim = DimLess]():
    _dimension: DimEls

    def __init__(self, values: Number | TensorData):
        """CAREFUL!! User should never instantiante any tensor using `__init__`.
        Always use the classmethod `new`."""
        if not hasattr(self.__class__, "_dimension"):
            raise RuntimeError("!!")
        
        self._values = np.array(values)

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        dim_els: DimEls = getattr(dim, "_d")

        return dimensional_tensor_class_factory(
            dim_els,
            cast(type[Tensor], cls)
        )
    
    def copy(self) -> Self:
        return self.__class__(self._values)

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
        
    def __pos__(self) -> Self:
        return self.__class__(+self._values)
    
    def __neg__(self) -> Self:
        return self.__class__(-self._values)

def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[SomeDim] | Tensor[DimLess]:
    """ensure tensor. if not a tensor object, creates a dimless tensor"""

    if isinstance(o, Tensor):
        return o
    elif isinstance(o, Number | np.ndarray):
        return Tensor[DimLess](o)

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

def scalar_class_factory(dim: DimEls) -> type[Scalar]:
    return cast(
        type[Scalar], 
        dimensional_tensor_class_factory(
            dim,
            Scalar,
        )
    )

def vector3_class_factory(dim: DimEls) -> type[Vector3]:
    return cast(
        type[Vector3], 
        dimensional_tensor_class_factory(
            dim,
            Vector3,
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
        return scalar_class_factory(self._dimension)(self._values[0])
    
    @property
    def z(self) -> Scalar[SomeDim]:
        return scalar_class_factory(self._dimension)(self._values[0])
    
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

    ### + OPERATOR ###

    @overload
    def __add__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __add__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __add__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__add__(o))

    @overload
    def __radd__(self: Vector3[DimLess], o: Vector3[DimLess]) -> Vector3[DimLess]: ...

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
    def __rsub__(self: Vector3[DimLess], o: Vector3[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __rsub__(self: Vector3[SomeDim], o: Vector3[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __rsub__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __rsub__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rsub__(o))

    ### * OPERATOR ###

    @overload
    def __mul__(self: Vector3[DimLess], o: Number | Scalar[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __mul__(self: Vector3[DimLess], o: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Number | Scalar[DimLess]) -> Vector3[SomeDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Vector3[Any]: ...

    def __mul__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Vector3[DimLess], o: Number | Scalar[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __rmul__(self: Vector3[DimLess], o: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __rmul__(self: Vector3[SomeDim], o: Number | Scalar[DimLess]) -> Vector3[SomeDim]: ...

    @overload
    def __rmul__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __rmul__(self, o: Scalar) -> Vector3[Any]: ...

    def __rmul__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Vector3[DimLess], o: Number | Scalar[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __truediv__(self: Vector3[DimLess], o: Scalar[SomeDim]) -> Vector3[QuotientDim[DimLess, SomeDim]]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Number | Scalar[DimLess]) -> Vector3[SomeDim]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeDim]) -> Vector3[DimLess]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Vector3[Any]: ...

    def __truediv__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Vector3[DimLess], o: Number | Scalar[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __rtruediv__(self: Vector3[SomeDim], o: Number | Scalar[DimLess]) -> Vector3[QuotientDim[DimLess, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Vector3[Any]:
        return cast(Vector3[Any], super().__rtruediv__(o))
    
    ### @ OPERATOR (dot product) ###

    @overload
    def __matmul__(self: Vector3[DimLess], o: Vector3[DimLess]) -> Scalar[DimLess]: ...

    @overload
    def __matmul__(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> Scalar[Any]:
        if isinstance(o, Vector3):
            return self.dot(o)
        
        raise TypeError("...")
    
class Matrix33[SomeDim](Tensor[SomeDim]):
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
        dim_els = getattr(dim, "_d", None)
        if dim_els == self._dimension:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @overload
    def inverse(self: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def inverse(self: Matrix33[SomeDim]) -> Matrix33[QuotientDim[DimLess, SomeDim]]: ...
    
    def inverse(self) -> Matrix33:
        raise NotImplementedError()

    ### + OPERATOR ###

    @overload
    def __add__(self: Matrix33[DimLess], o: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __add__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __add__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __add__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__add__(o))

    @overload
    def __radd__(self: Matrix33[DimLess], o: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __radd__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __radd__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __radd__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__radd__(o))

    ### - OPERATOR ###

    @overload
    def __sub__(self: Matrix33[DimLess], o: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __sub__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __sub__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __sub__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__sub__(o))

    @overload
    def __rsub__(self: Matrix33[DimLess], o: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __rsub__(self: Matrix33[SomeDim], o: Matrix33[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __rsub__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    def __rsub__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rsub__(o))

    ### * OPERATOR (element-wise scaling) ###

    @overload
    def __mul__(self: Matrix33[DimLess], o: Number | Scalar[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __mul__(self: Matrix33[DimLess], o: Scalar[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Number | Scalar[DimLess]) -> Matrix33[SomeDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __mul__(self, o: Scalar) -> Matrix33[Any]: ...

    def __mul__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__mul__(o))

    @overload
    def __rmul__(self: Matrix33[DimLess], o: Number | Scalar[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __rmul__(self: Matrix33[DimLess], o: Scalar[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __rmul__(self: Matrix33[SomeDim], o: Number | Scalar[DimLess]) -> Matrix33[SomeDim]: ...

    @overload
    def __rmul__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __rmul__(self, o: Scalar) -> Matrix33[Any]: ...

    def __rmul__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rmul__(o))

    ### / OPERATOR ###

    @overload
    def __truediv__(self: Matrix33[DimLess], o: Number | Scalar[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __truediv__(self: Matrix33[DimLess], o: Scalar[SomeDim]) -> Matrix33[QuotientDim[DimLess, SomeDim]]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Number | Scalar[DimLess]) -> Matrix33[SomeDim]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Scalar[SomeDim]) -> Matrix33[DimLess]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self, o: Scalar) -> Matrix33[Any]: ...

    def __truediv__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__truediv__(o))

    @overload
    def __rtruediv__(self: Matrix33[DimLess], o: Number | Scalar[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __rtruediv__(self: Matrix33[SomeDim], o: Number | Scalar[DimLess]) -> Matrix33[QuotientDim[DimLess, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> Matrix33[Any]:
        return cast(Matrix33[Any], super().__rtruediv__(o))

    ### @ OPERATOR (matrix multiplication) ###

    @overload
    def __matmul__(self: Matrix33[DimLess], o: Matrix33[DimLess]) -> Matrix33[DimLess]: ...

    @overload
    def __matmul__(self: Matrix33[DimLess], o: Matrix33[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Matrix33[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[DimLess], o: Vector3[DimLess]) -> Vector3[DimLess]: ...

    @overload
    def __matmul__(self: Matrix33[DimLess], o: Vector3[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Vector3[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __matmul__(self, o: object) -> Matrix33[Any] | Vector3[Any]:
        if isinstance(o, Matrix33):
            return cast(Matrix33[Any], super().__matmul__(o))
        elif isinstance(o, Vector3):
            return cast(Vector3[Any], super().__matmul__(o))

        raise TypeError("@ operation requires Matrix33 or Vector3")

###########################################

class QuantityMeta(type):
    def __getattr__(cls, name: str) -> Scalar:
        if name == "rad":
            return Scalar[AngleD].new(1)
        elif name == "m":
            return Scalar[LengthD].new(LengthD.meter)
        elif name == "km":
            return Scalar[LengthD].new(LengthD.kilo_meter)
        elif name == "s":
            return Scalar[TimeD].new(TimeD.second)
        elif name == "min":
            return Scalar[TimeD].new(TimeD.minute)
        elif name == "hour":
            return Scalar[TimeD].new(TimeD.hour)
        elif name == "day":
            return Scalar[TimeD].new(TimeD.day)
        elif name == "month":
            return Scalar[TimeD].new(TimeD.month)
        elif name == "year":
            return Scalar[TimeD].new(TimeD.year)
        
        raise ValueError(f"No unit named '{name}'")
        

class Quantity(metaclass=QuantityMeta):
    # ANGLES

    rad: Scalar[AngleD]
    """radians"""

    # DISTANCES

    m: Scalar[LengthD]
    """meter"""

    km: Scalar[LengthD]
    """kilometer"""

    # DURATIONS

    s: Scalar[TimeD]
    """second"""

    min: Scalar[TimeD]
    """minute"""

    hour: Scalar[TimeD]
    """hour"""

    day: Scalar[TimeD]
    """day"""

    month: Scalar[TimeD]
    """month (30 days)"""

    year: Scalar[TimeD]
    """year (365 days)"""


#########################################

### TRANSFORMS ###

T1 = TypeVar("T1")
T2 = TypeVar("T2")

class Transform(Generic[T1, T2], ABC):
    @abstractmethod
    def do(self, tensor: T1) -> T2: ...

    @abstractmethod
    def undo(self, tensor: T2) -> T1: ...

    @abstractmethod
    def copy(self) -> Self: ...

    def reverse(self) -> Transform[T2, T1]:
        t = self.copy()

        do = t.do
        undo = t.undo

        tt = cast(Transform[T2, T1], t)

        tt.do = undo  # type: ignore
        tt.undo = do  # type: ignore

        return tt
    
class TransformIdentify(Generic[T1], Transform[T1, T1]):
    def do(self, tensor: T1) -> T1:
        return tensor
    
    def undo(self, tensor: T1) -> T1:
        return tensor
    
    def copy(self) -> Self:
        t = TransformIdentify[T1]()
        return cast(Self, t)
    
class TransformVector3Linear(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[DimLess]):
        self.matrix = matrix
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix @ tensor
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix.inverse() @ tensor
    
    def copy(self) -> Self:
        t = TransformVector3Linear[SomeDim](
            matrix=self.matrix.copy()
        )
        return cast(Self, t)

class TransformVector3Affine(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[DimLess], translation: Vector3[SomeDim]):
        self.matrix = matrix
        self.translation = translation
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix @ tensor + self.translation
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix.inverse() @ (tensor - self.translation)
    
    def copy(self) -> Self:
        t = TransformVector3Affine[SomeDim](
            matrix=self.matrix.copy(), 
            translation=self.translation.copy()
        )
        return cast(Self, t)
    
class TransformVector3RotationZ(Generic[SomeDim], TransformVector3Linear[SomeDim]):
    def __init__(self, angle_rad: Number | Scalar[AngleD]):
        if isinstance(angle_rad, Scalar):
            angle_rad = angle_rad.magnitude("rad")
        
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        rot_mat = Matrix33[DimLess].new(
            c, -s, 0,
            s, c, 0,
            0, 0, 1
        )

        super().__init__(rot_mat)
    
class TransformChain(Generic[T1, T2], Transform[T1, T2]):
    def __init__(self, *transforms: Transform[Any, Any]):
        """executed in the given order"""
        self.transforms = list(transforms)

    def do(self, tensor: T1) -> T2:
        result: Any = tensor
        for t in self.transforms:
            result = t.do(result)
        return cast(T2, result)

    def undo(self, tensor: T2) -> T1:
        result: Any = tensor
        for t in reversed(self.transforms):
            result = t.undo(result)
        return cast(T1, result)
    
    def copy(self) -> Self:
        t = TransformChain[T1, T2](
            *(t.copy() for t in self.transforms)
        )

        return cast(Self, t)


a = Scalar[LengthD].new(1)
b = Scalar[DimLess].new(2)
c = Scalar[TimeD].new(3)

z = (a / c).cast(VelocityD)