from copy import copy
from dataclasses import dataclass
from decimal import Decimal
from enum import Enum
from fractions import Fraction
from functools import cached_property, wraps
import inspect
from itertools import chain, repeat
from types import EllipsisType
from typing import Annotated, Any, Callable, ClassVar, Generic, Iterable, Literal, NamedTuple, Self, Type, TypeAlias, TypeIs, TypeVar, TypedDict, cast, get_args, get_origin, get_type_hints
import operator
from numbers import Real as RealNumber
from numbers import Rational as RationalNumber

import numpy as np
import numpy.typing as npt

type NumpyFloatArray = npt.NDArray[np.float64]


class DimTriplet:
    """Exponent triplet describing a physical dimension.

    Coordinates are stored as rational exponents on the base axes:
    length (L), time (T), and mass (M).
    """

    def __init__(self,
        length: RationalNumber | float | int = 0,
        time: RationalNumber | float | int = 0,
        mass: RationalNumber | float | int = 0,
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
    
    def __pow__(self, p: RationalNumber | int) -> DimTriplet:
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
    
    def __pow__(cls, p: RationalNumber | int) -> type[Dim]:
        """Raise a dimension to a scalar power."""
        cls = cast(type[Dim], cls)
        if not isinstance(p, (RationalNumber, int)):
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

class TensorBinaryOperator(Enum):
    ADD = "+"
    SUB = "-"
    MUL = "*"
    MATMUL = "@"
    TRUEDIV = "/"
    FLOORDIV = "//"
    MODULO = "%"
    
    @property
    def operator(self) -> Callable[[NumpyFloatArray, NumpyFloatArray], NumpyFloatArray]:
        if self == TensorBinaryOperator.ADD:
            return operator.add
        elif self == TensorBinaryOperator.SUB:
            return operator.sub
        elif self == TensorBinaryOperator.MUL:
            return operator.mul
        elif self == TensorBinaryOperator.MATMUL:
            return operator.matmul
        elif self == TensorBinaryOperator.TRUEDIV:
            return operator.truediv
        elif self == TensorBinaryOperator.FLOORDIV:
            return operator.floordiv
        elif self == TensorBinaryOperator.MODULO:
            return operator.mod
        
        raise NotImplementedError(f"Unsupported operator '{self.value}'.")
    
    def dimension_result(self, left_dim: type[Dim], right_dim: type[Dim]) -> type[Dim] | None:
        """Return the resulting dimension of applying this operator to quantities of the given dimensions."""
        if self in (TensorBinaryOperator.ADD, TensorBinaryOperator.SUB):
            if left_dim != right_dim:
                return None
            return left_dim
        elif self == TensorBinaryOperator.MODULO:
            if left_dim != right_dim and not right_dim.triplet().dimensionless:
                return None
            return left_dim
        elif self in (TensorBinaryOperator.MUL, TensorBinaryOperator.MATMUL):
            return left_dim * right_dim
        elif self in (TensorBinaryOperator.TRUEDIV, TensorBinaryOperator.FLOORDIV):
            return left_dim / right_dim
        
        raise NotImplementedError(f"Unsupported operator '{self.value}' for dimension composition.")

class TensorKind(Enum):
    SCALAR = "scalar"
    VECTOR3 = "vector3"
    MATRIX33 = "matrix33"


class Tensor:
    """Generic n-dimensional tensor carrying a physical dimension.
    """
    _data: NumpyFloatArray
    _phy_dimension: type[Dim]

    def __init__(self, data: NumpyFloatArray | RealNumber, dimension: type[Dim] | None = None):
        """Initialize a tensor with the given data and dimension.
        Data units are default SI units corresponding to dimension."""
        if dimension is None:
            dimension = Dimless
        assert dimension.triplet() is not None

        self._data = np.asarray(data, dtype=np.float64)
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
    
    def raw_data_array(self, units: float | str) -> NumpyFloatArray:
        """Return the raw data array of this tensor, converted to the given units."""
        if isinstance(units, str):
            factor = Quantity.get(units).scalar
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

    def __pos__(self) -> Tensor:
        return Tensor(
            data=+self._data,
            dimension=self.phy_dimension
        )
    
    def __neg__(self) -> Tensor:
        return Tensor(
            data=-self._data,
            dimension=self.phy_dimension
        )
    
    def __pow__(self, exponent: RationalNumber) -> Tensor:
        return Tensor(
            data=self._data ** float(exponent),
            dimension=self.phy_dimension ** exponent
        )

    def __add__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.ADD)
    
    def __radd__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.ADD)
    
    def __sub__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.SUB)
    
    def __rsub__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.SUB)
    
    def __mul__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.MUL)
    
    def __rmul__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.MUL)
    
    def __truediv__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.TRUEDIV)
    
    def __rtruediv__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.TRUEDIV)
    
    def __floordiv__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.FLOORDIV)
    
    def __rfloordiv__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.FLOORDIV)
    
    def __mod__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.MODULO)
    
    def __rmod__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.MODULO)
    
    def __matmul__(self, right: Tensor | RealNumber) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.MATMUL)
    
    def __rmatmul__(self, left: RealNumber) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.MATMUL)
    
    
    def perform_binary_operation(self, other: Tensor | RealNumber, op: TensorBinaryOperator) -> Tensor:
        """Perform the given binary operation with another tensor, checking dimension compatibility."""
        if not isinstance(other, Tensor):
            other = scalar(other)

        dimension = op.dimension_result(self.phy_dimension, other.phy_dimension)
        if dimension is None:
            raise ValueError("Dimensions incompatible for given operator")
        
        return Tensor(
            data=op.operator(self._data, other._data), 
            dimension=dimension
        )

    @property
    def scalar(self) -> float:
        """Returns this tensor as scalar value, if it is one sized."""
        assert self.kind == TensorKind.SCALAR
        assert self.size == 1

        return self._data.item()
    
    
    
    @property
    def vec3(self) -> ElementsVector3:
        """Returns this tensor as vector3 value, if it is one sized."""
        assert self.kind == TensorKind.VECTOR3
        assert self.size == 1

        x, y, z = self._data.flatten()

        return ElementsVector3(x=x, y=y, z=z)
    
    @property
    def mat33(self) -> ElementsMatrix33:
        """Returns this tensor as matrix33 value, if it is one sized."""
        assert self.kind == TensorKind.MATRIX33
        assert self.size == 1

        a11, a21, a31, a12, a22, a32, a13, a23, a33 = self._data.flatten()

        return ElementsMatrix33(
            a11=a11, a21=a21, a31=a31,
            a12=a12, a22=a22, a32=a32,
            a13=a13, a23=a23, a33=a33
        )
    
class ElementsVector3(TypedDict):
    x: float
    y: float
    z: float

class ElementsMatrix33(TypedDict):
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

def tensor_check(f: Callable) -> Callable:
    """Wrapper that checks tensor inputs and outputs of a function based on type annotations."""
    signature = inspect.signature(f)
    hints = get_type_hints(f, include_extras=True)

    def _bound_from_annotation(name: str, annotation: Any, where: str) -> TensorBound:
        if get_origin(annotation) is not Annotated:
            raise TypeError(
                f"{where} '{name}' must be annotated as Annotated[Tensor, TensorBound(...)]."
            )

        args = get_args(annotation)
        if len(args) != 2:
            raise TypeError(
                f"{where} '{name}' must be annotated as Annotated[Tensor, TensorBound(...)]."
            )

        tensor_type, bound = args
        if tensor_type is not Tensor or not isinstance(bound, TensorBound):
            raise TypeError(
                f"{where} '{name}' must be annotated as Annotated[Tensor, TensorBound(...)]."
            )

        return bound

    parameter_bounds: dict[str, TensorBound] = {}
    for name in signature.parameters:
        if name not in hints:
            raise TypeError(
                f"Parameter '{name}' must be annotated as Annotated[Tensor, TensorBound(...)]."
            )
        parameter_bounds[name] = _bound_from_annotation(name, hints[name], "Parameter")

    if "return" not in hints:
        raise TypeError("Return annotation must be Annotated[Tensor, TensorBound(...)].")
    return_bound = _bound_from_annotation("return", hints["return"], "Return annotation")

    @wraps(f)
    def wrapper(*args, **kwargs):
        bound_arguments = signature.bind(*args, **kwargs)
        bound_arguments.apply_defaults()

        for name, value in bound_arguments.arguments.items():
            bound = parameter_bounds[name]
            parameter = signature.parameters[name]

            if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
                values_to_check = value
            elif parameter.kind is inspect.Parameter.VAR_KEYWORD:
                values_to_check = value.values()
            else:
                values_to_check = (value,)

            for checked_value in values_to_check:
                if not isinstance(checked_value, Tensor):
                    raise TypeError(f"Argument '{name}' must be a Tensor.")
                if not bound.check(checked_value):
                    raise ValueError(f"Argument '{name}' does not satisfy declared TensorBound.")

        result = f(*args, **kwargs)

        if not isinstance(result, Tensor):
            raise TypeError("Return value must be a Tensor.")
        if not return_bound.check(result):
            raise ValueError("Return value does not satisfy declared TensorBound.")

        return result

    return wrapper

__UNITS_REGISTRY: dict[str, Tensor] = dict()

def _units_register(units: list[str], dim: type[Dim], base_factor: float):
    global __UNITS_REGISTRY
    __UNITS_REGISTRY |= {
        u: Tensor(
            data=np.asarray(base_factor, dtype=np.float64).reshape((1,)), 
            dimension=dim
        )
        for u in units
    }
class Quantity(type):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """
    @classmethod
    def get(cls, name: str) -> Tensor:
        """Resolve a unit name into its corresponding scalar quantity."""
        global __UNITS_REGISTRY
        qt = __UNITS_REGISTRY.get(name)
        if qt is None:
            raise AttributeError(f"Unknown quantity '{name}'.")
        
        return copy(qt)
    

    dimensionless: ClassVar[Tensor]
    """dimensionless (1)"""
    _units_register(["dimensionless", "dimless"], Dimless, 1.)

    # ANGLES

    radian: ClassVar[Tensor]
    """radian"""
    _units_register(["radian", "rad"], Angle, 1.)

    turn: ClassVar[Tensor]
    """turns (360°)"""
    _units_register(["turn", "rev"], Angle, 2 * np.pi)

    degree: ClassVar[Tensor]
    """degree"""
    _units_register(["degree", "deg"], Angle, np.pi / 180)

    # DISTANCES

    meter: ClassVar[Tensor]
    """meter"""
    _units_register(["meter", "m"], Length, 1.)

    kilo_meter: ClassVar[Tensor]
    """kilometer"""
    _units_register(["kilometer", "km"], Length, 1e3)

    radii_earth: ClassVar[Tensor]
    """Mean radius of planet Earth (R🜨). 

    `R🜨 = 6378135 m`
    """
    _units_register(["radii_earth", "R🜨"], Length, 6378135)

    radii_sun: ClassVar[Tensor]
    """Mean radius of Sun (R☉). 

    `R☉ = 6.957e8 m`
    """
    _units_register(["radii_sun", "R☉"], Length, 6.957e8)


    astronomical_unit: ClassVar[Tensor]
    """Astronomical unit (au)."""
    _units_register(["astronomical_unit", "au"], Length, 149597870700)

    # DURATIONS

    second: ClassVar[Tensor]
    """second"""
    _units_register(["second", "s"], Time, 1.)


    minute: ClassVar[Tensor]
    """minute"""
    _units_register(["minute", "min"], Time, 60.)

    hour: ClassVar[Tensor]
    """hour"""
    _units_register(["hour", "h"], Time, 3600.)

    day: ClassVar[Tensor]
    """day"""
    _units_register(["day", "d"], Time, 86400.)

    month: ClassVar[Tensor]
    """month (30 days)"""
    _units_register(["month", "mo"], Time, 30 * 86400)

    year: ClassVar[Tensor]
    """year (365 days)"""
    _units_register(["year", "y"], Time, 365 * 86400)

    # MASSES

    kilo_gram: ClassVar[Tensor]
    """kilogram"""
    _units_register(["kilogram", "kg"], Mass, 1.)

    gram: ClassVar[Tensor]
    """gram"""
    _units_register(["gram", "g"], Mass, 1e-3)

    metric_ton: ClassVar[Tensor]
    """metric ton (1000 kg)"""
    _units_register(["metric_ton", "ton"], Mass, 1e3)

### CONVENIENCE FACTORY FUNCTIONS ###
### ############################# ###

def scalar(value: RealNumber) -> Tensor:
    """Create a dimensionless scalar tensor with the given value."""
    assert isinstance(value, RealNumber)

    return Tensor(
        data=np.asarray(value, dtype=np.float64).reshape((1,)), 
        dimension=Dimless
    )

def scalar_array(values: Iterable[RealNumber]) -> Tensor:
    """Create a dimensionless scalar tensor with the given values."""
    values = np.asarray(values, dtype=np.float64)
    assert values.ndim == 1

    return Tensor(
        data=values,
        dimension=Dimless
    )

def vec3(**elements: ElementsVector3) -> Tensor:
    """Create a dimensionless vector3 tensor with the given values."""
    return Tensor(
        data=np.asarray((elements["x"], elements["y"], elements["z"]), dtype=np.float64).reshape((3, 1)), 
        dimension=Dimless
    )

def vec3_array(values: Iterable[ElementsVector3]) -> Tensor:
    """Create a dimensionless vector3 tensor with the given values."""
    values = np.asarray(values, dtype=np.float64)
    assert values.ndim == 2 and values.shape[1] == 3

    return Tensor(
        data=values.transpose((1, 0)), 
        dimension=Dimless
    )

def mat33(**elements: ElementsMatrix33) -> Tensor:
    """Create a dimensionless matrix33 tensor with the given values."""
    return Tensor(
        data=np.asarray((
            elements["a11"], elements["a21"], elements["a31"],
            elements["a12"], elements["a22"], elements["a32"],
            elements["a13"], elements["a23"], elements["a33"]
        ), dtype=np.float64).reshape((3, 3, 1)), 
        dimension=Dimless
    )

def mat33_array(values: Iterable[ElementsMatrix33]) -> Tensor:
    """Create a dimensionless matrix33 tensor with the given values."""
    values = np.asarray(values, dtype=np.float64)
    assert values.ndim == 3 and values.shape[1:] == (3, 3)

    return Tensor(
        data=values.transpose((1, 2, 0)), 
        dimension=Dimless
    )