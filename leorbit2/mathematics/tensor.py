from abc import ABC, abstractmethod
import copy
from dataclasses import dataclass
from functools import singledispatchmethod
import inspect
from multiprocessing import Value
from pyclbr import Class
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Generic, Literal, Never, Self, TypeAlias, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np
import numpy.typing as npt

from leorbit2.mathematics.dimensions import Dim, DimCoords, D, registered_units, registered_dimensions

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
TensorData = np.typing.NDArray[np.floating[Any]]

TensorDataTransformer: TypeAlias = Callable[[TensorData], TensorData]

class Tensor[SomeDim = D.Dimless]():
    """Generic n-dimensional tensor carrying a physical dimension.

    Concrete subclasses specialize tensor shape and semantics while reusing
    dimension-aware arithmetic from this base class.
    """

    _dim: type[Dim]
    _base_tensor_class: type[Tensor]

    def __init__(self, values: Number | TensorData):
        """Initialize raw tensor values.

        Notes:
            End-users are expected to construct concrete tensors through
            subclass constructors/helpers such as ``new``.
        """
        if self._is_base_tensor_class():
            raise RuntimeError("...")
        
        self._values = np.array(values)

    @property
    def dim_coords(self) -> DimCoords:
        """Dimension coordinates associated with this tensor."""
        return self._dim._d
    
    @property
    def dim(self) -> type[SomeDim]:
        """Dimension class associated with this tensor."""
        return cast(
            type[SomeDim],
            self._dim
        )
    
    @classmethod
    def _is_base_tensor_class(cls) -> bool:
        """Whether ``cls`` is the undimensionalized root tensor class."""
        return not hasattr(cls, "_base_tensor_class")
    
    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[Tensor]:
        """Create a tensor subclass bound to ``dim_coords``."""
        if not cls._is_base_tensor_class():
            raise RuntimeError("Cannot call `dimensionalize` on a dimensionalized tensor class.")
        
        try:
            dim = D.get_dimension_from_coords(dim_coords)
        except ValueError:
            dim = type(
                f"DynamicDim : {dim_coords.representation}",
                (Dim, ),
                dict(
                    _d=copy.copy(dim_coords)
                )
            )
        
        return cast(
            type[Tensor],
            type(
                "DimensionalizedTensor",
                (cls, ),
                dict(
                    _dim=dim,
                    _base_tensor_class=cls
                )
            )
        )

    def __repr__(self) -> str:
        return f"Tensor[D.{self.dim.__class__.__name__}]({self._values})"

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        """Return a dimensionalized tensor class when given a concrete
        Dimension subclass; otherwise (type-checking / generics) return the
        original class so `Tensor[SomeDim]` works in annotations.
        """
        if inspect.isclass(dim) and issubclass(dim, Dim):
            return _dimensional_tensor_class_factory(
                dim,
                cast(type[Tensor], cls)
            )

        # Allow generic/type-var usage like `Tensor[SomeDim]` in annotations
        # by returning the original class when `dim` is not a concrete
        # Dimension subclass.
        return cls # type: ignore
    
    def cast(self, dim: type[SomeOtherDim]) -> Tensor[SomeOtherDim]:
        """Type-cast to another dimension if coordinates are identical."""
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    def copy(self) -> Self:
        """Return a value copy with the same tensor class and dimension."""
        return self.__class__(self._values)

    def ensure_compatible_dimensions(self: Tensor[SomeDim], o: Tensor[SomeOtherDim]) -> TypeIs[Tensor[SomeOtherDim]]:
        """Return whether two tensors have identical dimension coordinates."""
        return (
            isinstance(o, Tensor) 
            and o.dim_coords == self.dim_coords
        )
    
    def check(self, dim: type[Dim]) -> bool:
        """Return ``True`` if tensor dimension matches ``dim``."""
        return self.dim_coords == dim._d
    
    def get_raw_array(self, units: str = "1") -> npt.NDArray[np.float64]:
        """Return values converted to requested units.

        Args:
            units: Unit symbol registered in ``dimensions.registered_units``.
                Use ``"1"`` to retrieve base-unit values.
        """
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
    
    ### TENSOR OPERATIONS
    # Base operations rules
    # Addition : result is same dimension and same tensor type
    # Multiplication : result is product dimension but same tensor type
    
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
            return self._base_tensor_class.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values * o._values
            )
        except:
            raise RuntimeError("...")

    def __rmul__(self, o: object) -> Tensor: # o * self
        o = ensure_tensor(o)
        try:
            return self._base_tensor_class.dimensionalize(
                o.dim_coords * self.dim_coords
            )(
                o._values * self._values
            )
        except:
            raise RuntimeError("...")

    def __truediv__(self, o: object) -> Tensor: # self / o
        o = ensure_tensor(o)
        try:
            return self._base_tensor_class.dimensionalize(
                self.dim_coords / o.dim_coords
            )(
                self._values / o._values
            )
        except:
            raise RuntimeError("...")

    def __rtruediv__(self, o: object) -> Tensor: # o / self
        o = ensure_tensor(o)
        try:
            return self._base_tensor_class.dimensionalize(
                o.dim_coords / self.dim_coords
            )(
                o._values / self._values
            )
        except:
            raise RuntimeError("...")

    def __matmul__(self, o: object) -> Tensor: # self @ o
        o = ensure_tensor(o)
        try:
            return self._base_tensor_class.dimensionalize(
                self.dim_coords * o.dim_coords
            )(
                self._values @ o._values
            )
        except:
            raise RuntimeError("...")

    def __rmatmul__(self, o: object) -> Tensor: # o @ self
        o = ensure_tensor(o)
        try:
            return self._base_tensor_class.dimensionalize(
                o.dim_coords * self.dim_coords
            )(
                o._values @ self._values
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
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values == o._values)
            )
        except:
            raise RuntimeError("...")

    def __neq__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values != o._values)
            )
        except:
            raise RuntimeError("...")

    def __lt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values < o._values)
            )
        except:
            raise RuntimeError("...")

    def __le__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values <= o._values)
            )
        except:
            raise RuntimeError("...")

    def __gt__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values > o._values)
            )
        except:
            raise RuntimeError("...")

    def __ge__(self, o: object) -> Tensor[D.Dimless]:
        o = ensure_tensor(o)
        ensure_same_dimensions(self, o)

        try:
            return self._base_tensor_class[D.Dimless]( # type: ignore
                np.float64(self._values >= o._values)
            )
        except:
            raise RuntimeError("...")
    
SomeTensor = TypeVar("SomeTensor", bound=Tensor)

def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[SomeDim] | Tensor[D.Dimless]:
    """Return ``o`` as a tensor, wrapping numbers/arrays as dimensionless tensors."""

    if isinstance(o, Tensor):
        return o
    elif isinstance(o, Number | np.ndarray):
        return Tensor[D.Dimless](o)

    raise RuntimeError("...")

def ensure_same_dimensions(*tensors: Tensor[Any]) -> Literal[True]:
    """Validate that all tensors share the same dimension coordinates."""
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
    """Build a runtime tensor subclass bound to ``dim``."""
    base_tensor_class = parent_class if parent_class._is_base_tensor_class() else parent_class._base_tensor_class

    return type(
        "TensorWithDimension",
        (base_tensor_class, ),
        dict(
            _dim=dim,
            _base_tensor_class=base_tensor_class
        )
    )