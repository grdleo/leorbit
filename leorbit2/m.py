from abc import ABC, abstractmethod
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, cached_property
import inspect
from pyclbr import Class
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Generic, Literal, NamedTuple, Never, Self, TypeAlias, TypeGuard, TypeIs, TypeVar, cast, overload
from copy import copy
from unittest import result

import numpy as np
import numpy.typing as npt

class DimCoords:
    """Exponent triplet describing a physical dimension.

    Coordinates are stored as rational exponents on the base axes:
    length (L), time (T), and mass (M).
    """

    _DYNAMIC_DIM_REGISTRY: dict[DimCoords, type[Dim]] = dict()

    def __init__(self,
        length: Fraction | int = 0,
        time: Fraction | int = 0,
        mass: Fraction | int = 0,
    ):
        """Build dimension coordinates from base-axis exponents."""
        self.__length = Fraction(length)
        self.__time = Fraction(time)
        self.__mass = Fraction(mass)

    def __copy__(self) -> DimCoords:
        """Return a shallow copy of this coordinate object."""
        return DimCoords(
            length=self.__length,
            time=self.__time,
            mass=self.__mass
        )
    
    @cached_property
    def representation(self) -> str:
        """Return a stable human-readable representation, e.g. ``L 1 × T -2 × M 0``."""
        l, t, m = self.__length, self.__time, self.__mass

        sl = f"L {l.numerator}" + ("" if l.denominator == 1 else f"/{l.denominator}")
        st = f"T {t.numerator}" + ("" if t.denominator == 1 else f"/{t.denominator}")
        sm = f"M {m.numerator}" + ("" if m.denominator == 1 else f"/{m.denominator}")

        return " × ".join((sl, st, sm))
    
    def __repr__(self) -> str:
        return f"<DimCoords : {self.representation}>"

    def __hash__(self) -> int:
        """Hash based on the three exponent values so instances can be used as
        dictionary keys (registered_dimensions uses DimCoords as keys).
        """
        return hash((self.__length, self.__time, self.__mass))

    @property
    def length(self) -> Fraction:
        """Exponent of length axis ``L``."""
        return self.__length
    
    @property
    def time(self) -> Fraction:
        """Exponent of time axis ``T``."""
        return self.__time
    
    @property
    def mass(self) -> Fraction:
        """Exponent of mass axis ``M``."""
        return self.__mass

    @property
    def dimensionless(self) -> bool:
        """Whether all exponents are zero."""
        return (
            0
            == self.length
            == self.time
            == self.mass
        )
    
    def __eq__(self, o: object) -> bool:
        """Compare two coordinate triplets component-wise."""
        if not isinstance(o, DimCoords):
            return False
        
        return (
            self.length == o.length
            and self.time == o.time
            and self.mass == o.mass
        )
    
    def __mul__(self, o: DimCoords) -> DimCoords:
        """Compose dimensions by multiplying quantities (adds exponents)."""
        return DimCoords(
            length=self.length + o.length,
            time=self.time + o.time,
            mass=self.mass + o.mass
        )
    
    def __truediv__(self, o: DimCoords) -> DimCoords:
        """Compose dimensions by dividing quantities (subtracts exponents)."""
        return DimCoords(
            length=self.length - o.length,
            time=self.time - o.time,
            mass=self.mass - o.mass
        )
    
    def __pow__(self, p: Fraction | int | float) -> DimCoords:
        """Raise a dimension to a scalar power."""
        p = Fraction(p)

        return DimCoords(
            length=self.length * p,
            time=self.time * p,
            mass=self.mass * p
        )
    
    def to_dimension(self) -> type[Dim]:
        """Return the registered dimension class matching these coordinates.
        If no registered class matches, return a dynamic class with these coordinates."""
        try:
            return D.get_dimension_from_coords(self)
        except ValueError:
            pass
        
        class DynamicDim(Dim):
            _d = self

        return self.__class__._DYNAMIC_DIM_REGISTRY.setdefault(self, DynamicDim)

class Dim:
    """Base class for dimensions"""
    _d: ClassVar[DimCoords] = None # type: ignore (all derived dimensions must have this attribute)


class D:
    """Registry namespace for built-in dimensions and canonical unit factors."""

    class Dimless(Dim):
        """1"""
        _d = DimCoords()

    class Angle(Dim):
        """rad"""
        _d = DimCoords()

        rad: ClassVar[float] = 1 # base
        deg: ClassVar[float] = np.pi / 180

    class AngularVelocity(Dim):
        """rad/s"""
        _d = DimCoords(time=-1)

    class AngularAcc(Dim):
        """rad/s**2"""
        _d = DimCoords(time=-2)

    class AngularJerk(Dim):
        """rad/s**2"""
        _d = DimCoords(time=-3)

    class Length(Dim):
        """m"""
        _d = DimCoords(length=1)

        meter: ClassVar[float] = 1 # base

        milli_meter: ClassVar[float] = 1e-3 * meter
        kilo_meter: ClassVar[float] = 1e3 * meter

        radii_earth: ClassVar[float] = 6378135 * meter
        radii_sun: ClassVar[float] = 6.957e8 * meter

    class InvLength(Dim):
        """m**-1"""
        _d = DimCoords(length=-1)

    class Time(Dim):
        """s"""
        _d = DimCoords(time=1)

        second: ClassVar[float] = 1 # base

        milli_second: ClassVar[float] = 1e-3 * second
        minute: ClassVar[float] = 60 * second
        hour: ClassVar[float] = 60 * minute
        day: ClassVar[float] = 24 * hour
        month: ClassVar[float] = 30 * day
        year: ClassVar[float] = 365 * day

    class Mass(Dim):
        """kg"""
        _d = DimCoords(mass=1)

        kilo_gram: ClassVar[float] = 1 # base

        gram: ClassVar[float] = 1e-3 * kilo_gram
        milli_gram: ClassVar[float] = 1e-3 * gram
        metric_ton: ClassVar[float] = 1e3 * kilo_gram

    class Velocity(Dim):
        """m.s**-1"""
        _d = DimCoords(length=1, time=-1)

        meter_per_second: ClassVar[float] = 1 # base

        kilo_meter_per_hour: ClassVar[float] = meter_per_second / 3.6

    class Acceleration(Dim):
        """m.s**-2"""
        _d = DimCoords(length=1, time=-2)

    class Force(Dim):
        """kg.m.s**-2"""
        _d = DimCoords(mass=1, length=1, time=-2)

        newton: ClassVar[float] = 1 # base

    class GrativationnalParam(Dim):
        """m**3.s**-2"""
        _d = DimCoords(length=3, time=-2)

    @staticmethod
    def get_dimension_from_coords(dim_coords: DimCoords) -> type[Dim]:
        """Return the registered dimension class matching ``dim_coords``.

        Raises:
            ValueError: If no registered class matches.
        """
        dim = registered_dimensions().get(dim_coords, None)
        if dim is None:
            raise ValueError(f"No dimension registered with coords '{dim_coords}'")
        
        return dim

@cache
def registered_dimensions() -> dict[DimCoords, type[Dim]]:
    """Return a cached mapping from dimension coordinates to dimension classes."""
    return {
        dim_cls._d: dim_cls
        for dim_cls in D.__dict__.values()
        if isinstance(dim_cls, type) and issubclass(dim_cls, Dim)
    }

@cache
def registered_units() -> dict[str, tuple[type[Dim], Number]]:
    """Return a cached mapping from unit names to ``(dimension, factor)`` tuples."""
    return {
        u: (d, f)
        for d in registered_dimensions().values()
        for u, f in d.__dict__.items()
        if u != "_d" and isinstance(f, (float, int))
    }


Number: TypeAlias = float | int | np.floating
TensorData: TypeAlias = np.typing.NDArray[np.floating[Any]]
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)

class ProductDim[SomeDim, SomeOtherDim](Dim):
    """Type-level marker representing a product of two dimensions."""
    ...

class QuotientDim[SomeDim, SomeOtherDim](Dim):
    """Type-level marker representing a quotient of two dimensions."""
    ...
    
class PowerDim[SomeDim, Numerator, Denominator](Dim):
    """Type-level marker representing a powered dimension."""
    ...

N3: TypeAlias = Literal[-3]
N2: TypeAlias = Literal[-2]
N1: TypeAlias = Literal[-1]
OO: TypeAlias = Literal[0]
P1: TypeAlias = Literal[1]
P2: TypeAlias = Literal[2]
P3: TypeAlias = Literal[3]

class Tensor[SomeDim = D.Dimless]():
    """Generic n-dimensional tensor carrying a physical dimension.

    Concrete subclasses specialize tensor shape and semantics while reusing
    dimension-aware arithmetic from this base class.
    """

    _dim: ClassVar[type[Dim]]
    _base_tensor_class: ClassVar[type[Tensor]]

    def __init__(self, values: Number | TensorData):
        """Initialize raw tensor values.

        NOTE: `Tensor` should never be used as-is
        """
        if not self._is_dimensionalized():
            raise RuntimeError("Class has to be dimensionalized")
        
        self._values = np.asarray(values)

    def __hash__(self) -> int:
        return hash(f"{hash(self._values.data.tobytes())}${hash(self._dim._d)}")

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
    def _is_dimensionalized(cls) -> bool:
        dim = getattr(cls, "_dim", None)
        if dim is None:
            return False
        
        if not issubclass(dim, Dim):
            raise RuntimeError("?")
        
        return True

    def __repr__(self) -> str:
        return f"Tensor[D.{self.dim.__class__.__name__}]({self._values})"

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        """Return a dimensionalized tensor class when given a concrete
        Dimension subclass; otherwise (type-checking / generics) return the
        original class so `Tensor[SomeDim]` works in annotations.
        """
        if inspect.isclass(dim) and issubclass(dim, Dim):
            if cls._is_dimensionalized():
                raise RuntimeError("Cannot subscript a dimensionalized tensor class.")
            
            class DimensionalizedTensor(cls):
                _dim=dim
                # `_base_tensor_class` inherited from parent, that must have it!
            
            return cast(type[Tensor], DimensionalizedTensor)

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
    
    def __pos__(self) -> Self:
        return self.__class__(self._values)
    
    def __neg__(self) -> Self:
        return self.__class__(-self._values)
    
    def __add__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("...")
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["+"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[self.dim]( # type: ignore (is not specialized...)
            self._values + o._values
        )
    
    def __radd__(self, o: object) -> Tensor:
        return ensure_tensor(o) + self
    
    def __sub__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("...")
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["-"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[self.dim]( # type: ignore (is not specialized...)
            self._values - o._values
        )
    
    def __rsub__(self, o: object) -> Tensor:
        return ensure_tensor(o) - self
    
    def __mul__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["*"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[ # type: ignore (is not specialized...)
            (self.dim_coords * o.dim_coords).to_dimension()
        ](self._values * o._values)
    
    def __rmul__(self, o: object) -> Tensor:
        return ensure_tensor(o) * self
    
    def __truediv__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["/"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[ # type: ignore (is not specialized...)
            (self.dim_coords / o.dim_coords).to_dimension()
        ](self._values / o._values)
    
    def __rtruediv__(self, o: object) -> Tensor:
        return ensure_tensor(o) / self
    
    def __matmul__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["@"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[ # type: ignore (is not specialized...)
            (self.dim_coords * o.dim_coords).to_dimension()
        ](self._values @ o._values)
    
    def __rmatmul__(self, o: object) -> Tensor:
        return ensure_tensor(o) @ self
    
    def __mod__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("...")
        
        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["%"]

        if tensor_cls_result is None:
            raise ValueError("...")
        
        return tensor_cls_result[self.dim]( # type: ignore (is not specialized...)
            self._values % o._values
        )
    
    def __rmod__(self, o: object) -> Tensor:
        return ensure_tensor(o) % self
    
    def __eq__(self, o: object) -> bool:
        o = ensure_tensor(o)

        if self._base_tensor_class is not o._base_tensor_class:
            return False

        if not self.ensure_compatible_dimensions(o):
            return False
        
        return np.array_equal(self._values, o._values)

    def __neq__(self, o: object) -> bool:
        return not self.__eq__(o)

def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[SomeDim] | Tensor[D.Dimless]:
    """Return ``o`` as a tensor, wrapping numbers/arrays as dimensionless tensors."""

    if isinstance(o, Tensor):
        return o
    elif isinstance(o, Number):
        return Scalar[D.Dimless](np.asarray(o))
    elif isinstance(o, np.ndarray):
        arr = np.asarray(o)

        if arr.ndim == 0:
            return Scalar[D.Dimless](arr.item())
        if arr.ndim == 1:
            return ScalarArray[D.Dimless](arr)
        if arr.ndim == 2 and arr.shape[0] == 3:
            return Vector3Array[D.Dimless](arr)
        if arr.ndim == 2 and arr.shape == (3, 3):
            return Matrix33[D.Dimless](arr)

        raise RuntimeError("Unsupported array shape for tensor conversion")

    raise RuntimeError("...")

def ensure_same_dimensions(*tensors: Tensor[Any]) -> Literal[True]:
    """Validate that all tensors share the same dimension coordinates."""
    if len(tensors) == 0:
        return True
    
    t0, *_ = tensors
    if all(t.dim_coords == t0.dim_coords for t in tensors):
        return True
    
    raise RuntimeError("...")

from fractions import Fraction
from typing import Any, Generic, Literal, Never, TypeVar, cast, overload

class Tensor_S(Tensor[SomeDim], Generic[SomeDim]):
    """Core scalar tensor type that supports one or many scalar elements."""

    def __init__(self, data: TensorData):
        data = np.asarray(data)

        if data.ndim == 0:
            super().__init__(data.item())
        elif data.ndim == 1:
            super().__init__(data)
        else:
            raise ValueError("Wrong shape")

    def cast(self, dim: type[SomeOtherDim]) -> Tensor_S[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self  # type: ignore
        raise RuntimeError("Cannot cast")
    
    @property
    def size(self) -> int:
        shape = np.asarray(self._values).shape
        if len(shape) == 0:
            return 1
        s, = shape
        return int(s)

    @property
    def base_unit_value(self) -> Number | np.ndarray:
        if np.asarray(self._values).ndim == 0:
            return np.float64(self._values)
        return cast(Number, np.float64(self._values))

    def magnitude(self, units: str = "1") -> np.float64 | npt.NDArray[np.float64]:
        raw = self.get_raw_array(units)
        if np.asarray(raw).ndim == 0:
            return np.float64(raw)
        return cast(npt.NDArray[np.float64], raw)
    
    def _comparison(self, o: object, comparison: Callable[[npt.NDArray, npt.NDArray], npt.NDArray]) -> bool:
        other = ensure_tensor(o)

        if other._base_tensor_class is not Tensor_S:
            raise ValueError("Can only compare scalars")
        
        other = cast(Tensor_S, other)
        if self.size * other.size > 1:
            raise ValueError("Can only compare unidimensionnal scalars")
        
        return bool(comparison(self._values, other._values))

    def __lt__(self, o: object) -> bool:
        return self._comparison(o, np.less)

    def __le__(self, o: object) -> bool:
        return self._comparison(o, np.less_equal)

    def __gt__(self, o: object) -> bool:
        return self._comparison(o, np.greater)

    def __ge__(self, o: object) -> bool:
        return self._comparison(o, np.greater_equal)

Tensor_S._base_tensor_class = Tensor_S

class Tensor_V3(Tensor[SomeDim], Generic[SomeDim]):
    """Core vector3 tensor type that supports one or many elements."""

    def __init__(self, data: TensorData):
        data = np.asarray(data)
        rows, *_ = data.shape
        if data.ndim != 2 or rows != 3:
            raise ValueError("Wrong shape (expected 3 rows)")
        
        super().__init__(data)

    @classmethod
    def from_components(
        cls,
        x: Tensor_S[SomeDim] | Number | npt.NDArray,
        y: Tensor_S[SomeDim] | Number | npt.NDArray,
        z: Tensor_S[SomeDim] | Number | npt.NDArray,
    ) -> Tensor_V3[SomeDim]:
        vec_els = cast(list[object], [x, y, z])
        if all(isinstance(el, (int, float, np.floating)) for el in vec_els):
            return cls(
                np.array(vec_els).reshape((3,1))
            )
        
        if all(isinstance(el, (np.ndarray, list)) for el in vec_els):
            arrs = [np.asarray(el) for el in vec_els]
            if not all(arr.ndim == 1 for arr in arrs):
                raise ValueError("...")
            
            return cls(
                np.stack(arrs)
            )
        
        if not all(isinstance(el, Tensor_S) for el in vec_els):
            raise ValueError("...")
        
        x = cast(Tensor_S[SomeDim], x)
        y = cast(Tensor_S[SomeDim], y)
        z = cast(Tensor_S[SomeDim], z)

        ensure_same_dimensions(x, y, z)
        if cls._dim._d != x.dim_coords:
            raise ValueError("...")

        x_vals = np.asarray(x.base_unit_value).reshape(-1)
        y_vals = np.asarray(y.base_unit_value).reshape(-1)
        z_vals = np.asarray(z.base_unit_value).reshape(-1)

        if not (x_vals.size == y_vals.size == z_vals.size):
            raise ValueError("...")
        
        return cls(
            np.stack([x_vals, y_vals, z_vals])
        )
    
    @staticmethod
    def from_spherical(
        theta: Tensor_S[D.Angle],
        delta: Tensor_S[D.Angle],
        rho: Tensor_S[SomeOtherDim],
    ) -> Tensor_V3[SomeOtherDim]:
        
        cos_delta = cos(delta)
        sin_delta = sin(delta)
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        x = rho * cos_theta * cos_delta
        y = rho * sin_theta * cos_delta
        z = rho * sin_delta

        return Tensor_V3[rho.dim].from_components(
            cast(Tensor_S[SomeOtherDim], x),
            cast(Tensor_S[SomeOtherDim], y),
            cast(Tensor_S[SomeOtherDim], z),
        )
    
    def dot(self, o: Tensor_V3[SomeOtherDim]) -> Tensor_S[ProductDim[SomeDim, SomeOtherDim]]:
        result_dim = (self.dim_coords * o.dim_coords).to_dimension()
        return Tensor_S[result_dim]( # type: ignore
            np.sum(self._values * o._values, axis=0)
        )
    
    def cross(self, o: Tensor_V3[SomeOtherDim]) -> Tensor_V3[ProductDim[SomeDim, SomeOtherDim]]:
        result_dim = (self.dim_coords * o.dim_coords).to_dimension()
        return Tensor_V3[result_dim]( # type: ignore
            np.cross(self._values, o._values, axis=0)
        )

    @property
    def size(self) -> int:
        _, s = self._values.shape
        return int(s)

    def cast(self, dim: type[SomeOtherDim]) -> Tensor_V3[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self  # type: ignore
        raise RuntimeError("Cannot cast")

    @cached_property
    def x(self) -> Tensor_S[SomeDim]:
        return Tensor_S[self.dim](self._values[0, :])

    @cached_property
    def y(self) -> Tensor_S[SomeDim]:
        return Tensor_S[self.dim](self._values[1, :])

    @cached_property
    def z(self) -> Tensor_S[SomeDim]:
        return Tensor_S[self.dim](self._values[2, :])

    @cached_property
    def length(self) -> Tensor_S[SomeDim]:
        return sqrt(self.length_squared) # type: ignore
    
    @cached_property
    def length_squared(self) -> Tensor_S[PowerDim[SomeDim, P2, P1]]:
        return self.dot(self) # type: ignore

    @cached_property
    def theta(self) -> Tensor_S[D.Angle]:
        return cast(Tensor_S[D.Angle], atan2(self.y, self.x))

    @cached_property
    def delta(self) -> Tensor_S[D.Angle]:
        xy = sqrt(square(self.x) + square(self.y))
        return cast(Tensor_S[D.Angle], atan2(self.z, xy))

    def angle(self, o: Tensor_V3[SomeDim]) -> Tensor_S[D.Angle]:
        dot: np.ndarray = np.sum(self._values * o._values, axis=0)
        cos_angle: np.ndarray = (dot / (self.length._values * o.length._values))

        angle = np.acos(cos_angle)
        angle[cos_angle >= 1] = 0
        angle[cos_angle <= -1] = np.pi

        return Tensor_S[D.Angle](angle)

    def normalized(self) -> Tensor_V3[D.Dimless]:
        return cast(Tensor_V3[D.Dimless], self / self.length)

Tensor_V3._base_tensor_class = Tensor_V3

class Tensor_M33(Tensor[SomeDim], Generic[SomeDim]):
    """3×3 matrix carrying a physical dimension."""

    def __init__(self, data: TensorData):
        data = np.asarray(data)
        rows, cols, *_ = data.shape
        if data.ndim != 2 or rows != 3 or cols != 3:
            raise ValueError("Wrong shape (expected 3×3)")
        
        super().__init__(data)
    
    @classmethod
    def from_elements(cls,
        a: Tensor_S[SomeDim] | Number, b: Tensor_S[SomeDim] | Number, c: Tensor_S[SomeDim] | Number,
        d: Tensor_S[SomeDim] | Number, e: Tensor_S[SomeDim] | Number, f: Tensor_S[SomeDim] | Number,
        g: Tensor_S[SomeDim] | Number, h: Tensor_S[SomeDim] | Number, i: Tensor_S[SomeDim] | Number,
    ) -> Tensor_M33[SomeDim]:
        """Create a matrix from row-major coefficients."""
        mat: TensorData
        mat_els = cast(list[object], [a, b, c, d, e, f, g, h, i])

        if all(isinstance(el, (int, float, np.floating)) for el in mat_els):
            mat = np.array(mat_els).reshape((3,3))
        elif all(isinstance(el, Tensor_S) and el.size == 1 for el in mat_els):
            scalars = cast(list[Tensor_S[Any]], mat_els)
            ensure_same_dimensions(*scalars)
            mat = np.array([el.base_unit_value for el in scalars]).reshape((3,3))
        else:
            raise RuntimeError("...")

        return cls(mat)
    
    def cast(self, dim: type[SomeOtherDim]) -> Tensor_M33[SomeOtherDim]:
        """Type-cast to another dimension when coordinates are identical."""
        if dim._d == self.dim_coords:
            return self # type: ignore
        raise RuntimeError("Cannot cast")
    
    @cached_property
    def det(self) -> Number:
        return np.linalg.det(self._values)
    
    def inverse(self) -> Tensor_M33:
        """Return matrix inverse.
        """
        return Tensor_M33[
            (self.dim_coords ** -1).to_dimension()
        ](
            np.linalg.inv(self._values)
        )
    
    def __repr__(self) -> str:
        return f"TensorMatrix33[D.{self.dim.__class__.__name__}]({self._values})"

Tensor_M33._base_tensor_class = Tensor_M33

### convenience classes

class Scalar(Tensor_S[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for one-element scalar tensors."""
    
class ScalarArray(Tensor_S[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for many-element scalar tensors."""

    def __getitem__(self, index: int) -> Scalar[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("...")
        
        return cast(
            Scalar[SomeDim],
            Scalar[self.dim](np.array(self._values)[index])
        )
    
class Vector3(Tensor_V3[SomeDim], Generic[SomeDim]):
    """convenience wrapper for three-element vector tensors."""

class Vector3Array(Tensor_V3[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for many-element vector tensors."""

    def __getitem__(self, index: int) -> Vector3[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("...")
        
        return cast(
            Vector3[SomeDim],
            Vector3[self.dim](self._values[:, index])
        )
    
class Matrix33(Tensor_M33[SomeDim], Generic[SomeDim]):
    """Convience wrapper..."""


######## functions

def abs(tensor: Tensor[Any]) -> Tensor[Any]:
    return tensor._base_tensor_class[tensor.dim]( # type: ignore
        np.abs(tensor._values)
    )

def square(tensor: Tensor[Any]) -> Tensor[Any]:
    return_dim = (tensor.dim_coords ** 2).to_dimension()
    return tensor._base_tensor_class[return_dim]( # type: ignore
        np.square(tensor._values)
    )

def cube(tensor: Tensor[Any]) -> Tensor[Any]:
    return_dim = (tensor.dim_coords ** 3).to_dimension()
    return tensor._base_tensor_class[return_dim]( # type: ignore
        np.power(tensor._values, 3)
    )

def sqrt(tensor: Tensor[Any]) -> Tensor[Any]:
    return_dim = (tensor.dim_coords ** Fraction(1, 2)).to_dimension()
    return tensor._base_tensor_class[return_dim]( # type: ignore
        np.sqrt(tensor._values)
    )

def cbrt(tensor: Tensor[Any]) -> Tensor[Any]:
    return_dim = (tensor.dim_coords ** Fraction(1, 3)).to_dimension()
    return tensor._base_tensor_class[return_dim]( # type: ignore
        np.cbrt(tensor._values)
    )

def cos(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    if not tensor.check(D.Angle):
        raise ValueError("cos() expects an angle-typed tensor")
    
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.cos(tensor._values)
    )

def sin(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    if not tensor.check(D.Angle):
        raise ValueError("sin() expects an angle-typed tensor")
    
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.sin(tensor._values)
    )

def tan(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    if not tensor.check(D.Angle):
        raise ValueError("tan() expects an angle-typed tensor")
    
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.tan(tensor._values)
    )

def acos(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    if not tensor.check(D.Dimless):
        raise ValueError("acos() expects a dimensionless tensor")
    
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.acos(tensor._values)
    )

def asin(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    if not tensor.check(D.Dimless):
        raise ValueError("asin() expects a dimensionless tensor")
    
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.asin(tensor._values)
    )

def atan(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    if not tensor.check(D.Dimless):
        raise ValueError("atan() expects a dimensionless tensor")
    
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.atan(tensor._values)
    )

def atan2(y: Tensor[SomeDim], x: Tensor[SomeDim]) -> Tensor[D.Angle]:
    """Elementwise two-argument arctangent that returns an angle-typed tensor.

    Returns a ``Scalar[D.Angle]`` when both inputs are scalar-like and a
    ``ScalarArray[D.Angle]`` when at least one input has array semantics.
    """
    if not y.ensure_compatible_dimensions(x) or y._base_tensor_class is not x._base_tensor_class:
        raise ValueError("Incompatible dimensions or tensor types")
    
    vals = np.atan2(y._values, x._values)

    return y._base_tensor_class[D.Angle](vals) # type: ignore

def normalize_angle(angle: Tensor[D.Angle]) -> Tensor[D.Angle]:
    """Returns the given angle in its [0, 2π] range."""
    return angle._base_tensor_class[D.Angle](angle._values % (2 * np.pi)) # type: ignore

def normalize_angle_symmetric(angle: Tensor[D.Angle]) -> Tensor[D.Angle]:
    """Returns the given angle in its [-π, π] range."""
    normalized = angle._values % (2 * np.pi)
    normalized[normalized > np.pi] -= 2 * np.pi

    return angle._base_tensor_class[D.Angle](normalized) # type: ignore

def interpolate(t1: Tensor[SomeDim], t2: Tensor[SomeDim], alpha: Number) -> Tensor[SomeDim]:
    """Return the linear interpolation between ``t1`` and ``t2`` with factor ``0 <= alpha <= 1``."""
    if not t1.ensure_compatible_dimensions(t2) or t1._base_tensor_class is not t2._base_tensor_class:
        raise ValueError("Incompatible dimensions or tensor types")
    
    return t1._base_tensor_class[ # type: ignore
        t1.dim_coords.to_dimension()
    ](
        (1 - alpha) * t1._values + alpha * t2._values
    )

########################


OPERATION_RESULT_TYPE: dict[
    tuple[type[Tensor], type[Tensor]],
    dict[str, type[Tensor] | None]
] = { # type: ignore
    (Tensor_S, Tensor_S): {
        "+": Tensor_S,
        "-": Tensor_S,
        "*": Tensor_S,
        "/": Tensor_S,
        "@": None,
        "%": Tensor_S
    },
    (Tensor_S, Tensor_V3): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3,
        "/": None,
        "@": None,
        "%": None
    },
    (Tensor_S, Tensor_M33): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33,
        "/": None,
        "@": None,
        "%": None
    },
    (Tensor_V3, Tensor_S): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3,
        "/": Tensor_V3,
        "@": None,
        "%": Tensor_V3
    },
    (Tensor_V3, Tensor_V3): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3, # term by term product
        "/": Tensor_V3, # term by term div
        "@": None,
        "%": Tensor_V3 # term by term modulo

    },
    (Tensor_V3, Tensor_M33): {
        "+": None,
        "-": None,
        "*": None,
        "/": None,
        "@": None, # wrong side mat product
        "%": None
    },
    (Tensor_M33, Tensor_S): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33,
        "/": Tensor_M33,
        "@": None,
        "%": Tensor_M33
    },
    (Tensor_M33, Tensor_V3): {
        "+": None,
        "-": None,
        "*": None,
        "/": None,
        "@": Tensor_V3, # mat product
        "%": None
    },
    (Tensor_M33, Tensor_M33): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33, # term by term product
        "/": Tensor_M33, # term by term div
        "@": Tensor_M33, # mat product
        "%": Tensor_M33 # term by term modulo
    },
}

#############################

class QuantityMeta(type):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """

    def __getattr__(cls, name: str) -> Scalar:
        """Resolve a unit name into its corresponding scalar quantity."""
        try:
            dim, factor = registered_units()[name]
            return Scalar[dim](np.asarray(factor))
        except KeyError:
            raise ValueError(f"No unit named '{name}'")   

class Quantity(metaclass=QuantityMeta):
    """Convenience namespace for creating unit-scaled scalar values."""

    @classmethod
    def get(cls, value: str) -> Scalar:
        """Return the scalar unit associated with ``value``."""
        return cls.__getattr__(value)

    # ANGLES

    rad: Scalar[D.Angle]
    """radians"""

    deg: Scalar[D.Angle]
    """degrees"""

    # DISTANCES

    meter: Scalar[D.Length]
    """meter"""

    kilo_meter: Scalar[D.Length]
    """kilometer"""

    radii_earth: Scalar[D.Length]
    """Mean radius of planet Earth (R🜨). 

    `R🜨 = 6378135 m`
    """

    radii_sun: Scalar[D.Length]
    """Mean radius of Sun (R☉). 

    `R☉ = 6.957e8 m`
    """

    # DURATIONS

    second: Scalar[D.Time]
    """second"""

    minute: Scalar[D.Time]
    """minute"""

    hour: Scalar[D.Time]
    """hour"""

    day: Scalar[D.Time]
    """day"""

    month: Scalar[D.Time]
    """month (30 days)"""

    year: Scalar[D.Time]
    """year (365 days)"""

    # VELOCITIES

    meter_per_second: Scalar[D.Velocity]
    """meter per second"""

    kilo_meter_per_hour: Scalar[D.Velocity]
    """kilometer per hour"""