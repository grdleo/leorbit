from copy import copy
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from functools import cached_property, wraps
import inspect
from itertools import chain, repeat
from types import EllipsisType
from typing import Any, Callable, ClassVar, Generic, Literal, NamedTuple, Self, Type, TypeAlias, TypeIs, TypeVar, cast

import numpy as np
import numpy.typing as npt


class DimTriplet:
    """Exponent triplet describing a physical dimension.

    Coordinates are stored as rational exponents on the base axes:
    length (L), time (T), and mass (M).
    """

    def __init__(self,
        length: Fraction | int = 0,
        time: Fraction | int = 0,
        mass: Fraction | int = 0,
    ):
        """Build dimension coordinates from base-axis exponents."""
        self.__length = Fraction(length)
        self.__time = Fraction(time)
        self.__mass = Fraction(mass)

    def __copy__(self) -> DimTriplet:
        """Return a shallow copy of this coordinate object."""
        return DimTriplet(
            length=self.__length,
            time=self.__time,
            mass=self.__mass
        )
    
    @cached_property
    def representation(self) -> str:
        """Return a stable human-readable representation, e.g. ``L^1 × T^-2 × M^0``."""
        l, t, m = self.__length, self.__time, self.__mass

        sl = f"L^{l.numerator}" + ("" if l.denominator == 1 else f"/{l.denominator}")
        st = f"T^{t.numerator}" + ("" if t.denominator == 1 else f"/{t.denominator}")
        sm = f"M^{m.numerator}" + ("" if m.denominator == 1 else f"/{m.denominator}")

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
        if not isinstance(o, DimTriplet):
            return False
        
        return (
            self.length == o.length
            and self.time == o.time
            and self.mass == o.mass
        )
    
    def __mul__(self, o: DimTriplet) -> DimTriplet:
        """Compose dimensions by multiplying quantities (adds exponents)."""
        return DimTriplet(
            length=self.length + o.length,
            time=self.time + o.time,
            mass=self.mass + o.mass
        )
    
    def __truediv__(self, o: DimTriplet) -> DimTriplet:
        """Compose dimensions by dividing quantities (subtracts exponents)."""
        return DimTriplet(
            length=self.length - o.length,
            time=self.time - o.time,
            mass=self.mass - o.mass
        )
    
    def __pow__(self, p: Fraction | int | float) -> DimTriplet:
        """Raise a dimension to a scalar power."""
        p = Fraction(p)

        return DimTriplet(
            length=self.length * p,
            time=self.time * p,
            mass=self.mass * p
        )
    
__DYNAMIC_DIM_REGISTRY: dict[DimTriplet, type[Dim]] = dict()
    
class _DimClassAlgebra(type):
    def __mul__(cls, o: type[Dim]) -> type[Dim]:
        """Compose dimensions by multiplying quantities (adds exponents)."""
        cls = cast(type[Dim], cls)
        if not issubclass(o, Dim):
            raise NotImplementedError()
        
        return _dimension_factory(
            cls.triplet() * o.triplet()
        )
    
    def __truediv__(cls, o: type[Dim]) -> type[Dim]:
        """Compose dimensions by dividing quantities (subtracts exponents)."""
        cls = cast(type[Dim], cls)
        if not issubclass(o, Dim):
            raise NotImplementedError()
        
        return _dimension_factory(
            cls.triplet() / o.triplet()
        )
    
    def __rtruediv__(cls, o: Literal[1]) -> type[Dim]:
        """Allows to compute the inverse dimension (1 / Dim)."""
        cls = cast(type[Dim], cls)
        if o != 1:
            raise NotImplementedError()
        
        return _dimension_factory(
            cls.triplet() ** -1
        )
    
    def __pow__(cls, p: Fraction | int | float) -> type[Dim]:
        """Raise a dimension to a scalar power."""
        cls = cast(type[Dim], cls)
        if not isinstance(p, (Fraction, int, float)):
            raise NotImplementedError()
        
        return _dimension_factory(
            cls.triplet() ** p
        )
    
class Dim(metaclass=_DimClassAlgebra):
    __triplet: ClassVar[DimTriplet]

    def __init__(self, *args, **kwargs):
        """Prevent instantiation of dimension classes."""
        raise RuntimeError(f"Cannot instantiate dimension class `{self.__class__.__name__}`.")

    @classmethod
    def triplet(cls) -> DimTriplet:
        """Return the coordinate triplet associated with this dimension."""
        if not hasattr(cls, "__triplet"):
            raise RuntimeError(f"`{cls.__name__}` is not a registered dimension class.")
        
        return getattr(cls, "__triplet")

    def __eq__(self, o: object) -> bool:
        """Compare dimensions by their coordinate triplets."""
        if not isinstance(o, type) or not issubclass(o, Dim):
            return False
        
        return self.triplet() == o.triplet()
    
def _dimension_factory(triplet: DimTriplet) -> type[Dim]:
    """Return a dimension class corresponding to the given coordinate triplet."""
    global __DYNAMIC_DIM_REGISTRY

    return __DYNAMIC_DIM_REGISTRY.setdefault(
        triplet,
        type(
            f"DynamicDim_{triplet.representation}",
            (Dim, ),
            dict(__triplet=triplet)
        )
    )
    ...

_ = Dimless = Angle = _dimension_factory(DimTriplet())
_ = Length = _dimension_factory(DimTriplet(length=1))
_ = Time = _dimension_factory(DimTriplet(time=1))
_ = Mass = _dimension_factory(DimTriplet(mass=1))
_ = Velocity = Length / Time
_ = Acceleration = Velocity / Time
_ = Force = Mass * Acceleration
_ = Frequency = AngularVelocity = 1 / Time
_ = AngularAcceleration = AngularVelocity / Time
_ = AngularJerk = AngularAcceleration / Time
_ = InvLength = 1 / Length
# NOTE: This is a trick discovered accidentally for Pyright to recognize these dimensions as "real types"
# instead of just `type[Dim]` which would be the case if we directly assigned the result of the operations to the variables.

@dataclass(frozen=True)
class U:
    second = s = 1.
    minute = min = 60. * second
    hour = h = 60. * minute
    day = 24. * hour

    meter = m = 1.

    kilogram = kg = 1.

    @classmethod
    def get_factor(cls, units: str) -> float:
        """Return the conversion factor from the given units to SI units."""
        if not hasattr(cls, units):
            raise ValueError(f"Unknown unit '{units}'.")
        
        return getattr(cls, units)

type TensorUnaryOperator = Literal["+", "-"]
type TensorBinaryOperator = Literal["+", "-", "*", "/", "@"]

class TensorKind(Enum):
    SCALAR = "scalar"
    VECTOR3 = "vector3"
    MATRIX33 = "matrix33"


class Tensor:
    """Generic n-dimensional tensor carrying a physical dimension.
    """
    _data: npt.NDArray[np.float64]
    _phy_dimension: type[Dim]

    def __init__(self, data: npt.NDArray[np.float64], dimension: type[Dim] | None = None):
        """Initialize a tensor with the given data and dimension.
        Data units are default SI units corresponding to dimension."""
        if dimension is None:
            dimension = Dimless
        assert dimension.triplet() is not None

        self._data = data
        self._phy_dimension = dimension

    def __hash__(self) -> int:
        return hash(f"{hash(self._data.data.tobytes())}${hash(self._phy_dimension.triplet())}")
    
    def __repr__(self) -> str:
        return f"<Tensor {self._data} [{self._phy_dimension.triplet().representation}]>"
    
    def __copy__(self) -> Tensor:
        """Return a shallow/deep copy of this tensor."""
        return Tensor(
            data=self._data.copy(),
            dimension=self._phy_dimension
        )
    
    @property
    def phy_dimension(self) -> type[Dim]:
        """Return the physical dimension of this tensor."""
        return copy(self._phy_dimension)
    
    @property
    def array_dimensions(self) -> int:
        """Return the number of array dimensions of this tensor."""
        return self._data.ndim
    
    @property
    def size(self) -> int:
        if self.array_dimensions == 1:
            return self._data.size
        elif self.array_dimensions == 2:
            return self._data.size // 3
        elif self.array_dimensions == 3:
            return self._data.size // 9
        
        raise RuntimeError("Unreachable code")
    
    @property
    def kind(self) -> TensorKind:
        if self.array_dimensions == 1:
            return TensorKind.SCALAR
        elif self.array_dimensions == 2:
            return TensorKind.VECTOR3
        elif self.array_dimensions == 3:
            return TensorKind.MATRIX33
        
        raise RuntimeError("Unreachable code")
        
    
    @property
    def shape(self) -> tuple[int, ...]:
        """Return the shape of this tensor."""
        return self._data.shape
    
    def check(self, 
        dimension: type[Dim] | None = None,
        kind: TensorKind | str | None = None,
        size: int | None = None
    ) -> bool:
        """Check if the dimension of this tensor matches the given one."""
        if (dimension == kind == size == None):
            raise ValueError("Nothing to check!")
        
        return (
            (self.phy_dimension == dimension if dimension is not None else True)
            and (self.kind == TensorKind(kind) if kind is not None else True)
            and (self.size == size if size is not None else True)
        )
    
    def secure(self, 
        dimension: type[Dim] | None = None,
        kind: TensorKind | str | None = None,
        size: int | None = None
    ) -> Self:
        """Returns self if the dimension of this tensor matches the given one, otherwise raises a ValueError."""
        if not self.check(dimension, kind, size):
            raise ValueError(f"Tensor does not match prerogatives.")
        
        return self
    
    def raw_data_array(self, units: float | str) -> npt.NDArray[np.float64]:
        """Return the raw data array of this tensor, converted to the given units."""
        if isinstance(units, str):
            factor = U.get_factor(units)
        else:
            factor = units
        return self._data * factor
    
    def __getitem__(self, index: int) -> Tensor:
        """..."""
        slices = chain(
            repeat(slice(None), self.array_dimensions - 1),
            [index]
        )

        tensor = Tensor(
            self._data[*slices],
            self.phy_dimension
        )

        assert tensor.array_dimensions == self.array_dimensions
        assert tensor.size == 1

        return tensor
    
    
    def perform_binary_operation(self, other: Tensor, operator: TensorBinaryOperator) -> Tensor:
        """Perform the given binary operation with another tensor, checking dimension compatibility."""
        if operator in ("+", "-"):
            if self.phy_dimension != other.phy_dimension:
                raise ValueError(f"Cannot perform '{operator}' on tensors with different dimensions: {self.phy_dimension.triplet().representation} and {other.phy_dimension.triplet().representation}.")
            data = self._data + other._data if operator == "+" else self._data - other._data
            dimension = self.phy_dimension
        elif operator in ("*", "@"):
            data = self._data * other._data if operator == "*" else self._data @ other._data
            dimension = self.phy_dimension * other.phy_dimension if operator == "*" else self.phy_dimension / other.phy_dimension
        elif operator == "/":
            data = self._data / other._data
            dimension = self.phy_dimension / other.phy_dimension
        else:
            raise NotImplementedError(f"Unsupported operator '{operator}'.")
        
        return Tensor(data=data, dimension=dimension)

    @property
    def scalar(self) -> float:
        """Returns this tensor as scalar value, if it is one sized."""
        assert self.kind == TensorKind.SCALAR
        assert self.size == 1

        return self._data.item()
    
    class ElementsVector3(NamedTuple):
        x: float
        y: float
        z: float
    
    @property
    def vec3(self) -> ElementsVector3:
        """Returns this tensor as vector3 value, if it is one sized."""
        assert self.kind == TensorKind.VECTOR3
        assert self.size == 1

        return Tensor.ElementsVector3(*self._data.flatten())
    
    class ElementsMatrix33(NamedTuple):
        """(a_ij) i: col; j: row"""
        a11: float
        a21: float
        a31: float
        a12: float
        a22: float
        a32: float
        a13: float
        a23: float
        a33: float
    
    @property
    def mat33(self) -> ElementsMatrix33:
        """Returns this tensor as matrix33 value, if it is one sized."""
        assert self.kind == TensorKind.MATRIX33
        assert self.size == 1

        return Tensor.ElementsMatrix33(*self._data.flatten())

@dataclass
class TensorBound:
    dimension: type[Dim] | None = None
    kind: TensorKind | str | None = None
    size: int | None = None

    def check(self, tensor: Tensor) -> bool:
        return tensor.check(
            dimension=self.dimension,
            kind=self.kind,
            size=self.size
        )

def tensor_inputs(**inputs: TensorBound | type[float]):
    def wrapper(f: Callable) -> Callable:
        @wraps(f)
        def wrapped(*args: Tensor) -> Tensor:
            sig = inspect.signature(f)
            for i, param in enumerate(sig.parameters.values()):
                if param.name not in inputs.keys():
                    raise ValueError(f"Parameter '{param.name}' is not declared in the tensor input specification.")
                
                bound = inputs[param.name]
                value: Any = args[i]
                if bound is float:
                    if not isinstance(value, float | int):
                        raise ValueError(f"Argument '{param.name}' is expected to be a number.")
                elif isinstance(bound, TensorBound):
                    if not bound.check(value):
                        raise ValueError(f"Argument '{param.name}' does not match the expected tensor bound.")
                else:
                    raise ValueError(f"Invalid tensor input specification for parameter '{param.name}'.")
            
            return f(*args)
        
        return wrapped

    return wrapper

def tensor_output(output: TensorBound | type[float]):
    def wrapper(f: Callable) -> Callable:
        @wraps(f)
        def wrapped(*args: Tensor) -> Tensor | float:
            o: Tensor | float = f(*args)

            if output is float:
                if not isinstance(o, float | int):
                    raise ValueError(f"Output is expected to be a number.")
            elif isinstance(output, TensorBound):
                if isinstance(o, float | int):
                    raise ValueError(f"Output is expected to be a tensor, but got a number.")
                if not output.check(o):
                    raise ValueError(f"Output does not match the expected tensor bound.")
            
            return o
        
        return wrapped

    return wrapper