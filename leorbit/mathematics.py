from copy import copy
from dataclasses import dataclass
from decimal import Decimal
from enum import Enum
from fractions import Fraction
from functools import cached_property, wraps
import inspect
from itertools import chain, repeat
from multiprocessing import Value
from re import M
from types import EllipsisType
from typing import Annotated, Any, Callable, ClassVar, Generic, Iterable, Literal, NamedTuple, Self, Type, TypeAlias, TypeIs, TypeVar, TypedDict, cast, get_args, get_origin, get_type_hints
import operator

import numpy as np
import numpy.typing as npt

type NumpyFloatArray = npt.NDArray[np.float64]
RealNumber = float | int | Fraction | Decimal


class DimTriplet:
    """Exponent triplet describing a physical dimension.

    Coordinates are stored as rational exponents on the base axes:
    length (L), time (T), and mass (M).
    """

    def __init__(self,
        length: RealNumber | float | int = 0,
        time: RealNumber | float | int = 0,
        mass: RealNumber | float | int = 0,
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
    
    def __pow__(self, p: RealNumber | int) -> DimTriplet:
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
    
    def __pow__(cls, p: RealNumber | int | float) -> type[Dim]:
        """Raise a dimension to a scalar power."""
        cls = cast(type[Dim], cls)
        if not isinstance(p, RealNumber):
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
        elif self in (TensorBinaryOperator.MUL, ):
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
        if not isinstance(data, np.ndarray):
            data = np.asarray(data, dtype=np.float64).reshape((1,))

        assert dimension.triplet() is not None

        if data.ndim == 0:
            data = data.reshape((1,))

        if data.ndim == 1:
            pass
        elif data.ndim == 2:
            rows, _ = data.shape
            if rows != 3:
                raise ValueError("Vector3 tensor data must have shape (3, N).")
        elif data.ndim == 3:
            rows, cols, _ = data.shape
            if rows != 3 or cols != 3:
                raise ValueError("Matrix33 tensor data must have shape (3, 3, N).")
        else:
            raise ValueError("Tensor data must be 1D, 2D, or 3D.")

        self._data = np.asarray(data, dtype=np.float64)
        self._phy_dimension = dimension

    @property
    def data(self) -> NumpyFloatArray:
        return self._data.copy()

    @property
    def dim_triplet(self) -> DimTriplet:
        return self._phy_dimension.triplet()

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

    def copy(self) -> Tensor:
        return self
    
    @property
    def phy_dimension(self) -> type[Dim]:
        """Return the physical dimension of this tensor."""
        return copy(self._phy_dimension)
    
    @property
    def array_dimensions(self) -> int:
        """Return the number of array dimensions of this tensor."""
        # NOTE: THIS IS APPROUVED. DO NOT TOUCH IT FFS.
        return self._data.ndim 
    
    @property
    def kind(self) -> TensorKind:
        # NOTE: THIS IS APPROUVED. DO NOT TOUCH IT FFS.
        if self.array_dimensions == 1:
            return TensorKind.SCALAR
        elif self.array_dimensions == 2:
            return TensorKind.VECTOR3
        elif self.array_dimensions == 3:
            return TensorKind.MATRIX33
        
        raise RuntimeError("Unreachable code")
    
    @property
    def size(self) -> int:
        # NOTE: THIS IS APPROUVED. DO NOT TOUCH IT FFS.
        if self.kind == TensorKind.SCALAR:
            return self._data.size
        elif self.kind == TensorKind.VECTOR3:
            return self._data.size // 3
        elif self.kind == TensorKind.MATRIX33:
            return self._data.size // 9
        
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
            return True
        
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
    
    def concatenate(self, other: Tensor) -> Tensor:
        """Concatenate this tensor with another one along the first axis, checking dimension compatibility."""
        if self.kind != other.kind:
            raise ValueError("Cannot concatenate tensors of different kinds.")
        if self.phy_dimension != other.phy_dimension:
            raise ValueError("Cannot concatenate tensors of different dimensions.")
        return Tensor(
            np.concatenate((self._data, other._data), axis=self.array_dimensions - 1),
            self.phy_dimension
        )
    
    def raw_data_array(self, units: Tensor | float | str | None = None) -> NumpyFloatArray:
        """Return the raw data array of this tensor, converted to the given units."""
        if isinstance(units, Tensor):
            if units.kind != TensorKind.SCALAR:
                raise ValueError("Units tensor must be a scalar.")
            if units.phy_dimension != self.phy_dimension:
                raise ValueError("Units tensor must have the same physical dimension as the tensor.")
            
            factor = units.scalar.value()
        elif isinstance(units, str):
            factor = Quantity.get(units).scalar.value()
        elif isinstance(units, (int, float)):
            factor = units
        elif units is None:
            factor = 1.0
        return self._data / factor
    
    def __getitem__(self, index: int) -> Tensor:
        """..."""
        slices = chain(
            repeat(slice(None), self.array_dimensions - 1),
            [index]
        )

        selected_data = np.asarray(self._data[*slices])
        if selected_data.ndim < self.array_dimensions:
            selected_data = np.expand_dims(selected_data, axis=-1)

        tensor = Tensor(
            selected_data,
            self.phy_dimension
        )

        assert tensor.array_dimensions == self.array_dimensions
        assert tensor.size == 1

        return tensor

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Tensor):
            return False
        
        return (
            self.phy_dimension == o.phy_dimension 
            and np.array_equal(self._data, o._data)
        )

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
    
    def __abs__(self) -> Tensor:
        return Tensor(
            data=np.abs(self._data),
            dimension=self.phy_dimension
        )
    
    def __pow__(self, exponent: RealNumber | int | float) -> Tensor:
        return Tensor(
            data=self._data ** float(exponent),
            dimension=self.phy_dimension ** exponent
        )

    def __lt__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() < 0)

    def __le__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() <= 0)

    def __gt__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() > 0)

    def __ge__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() >= 0)

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
    
    def __or__(self, right: Tensor) -> Tensor:
        """Concatenation operator, only works for tensors of the same kind and dimension."""
        return self.concatenate(right)
    
    def perform_binary_operation(self, other: Tensor | RealNumber, op: TensorBinaryOperator) -> Tensor:
        """Perform the given binary operation with another tensor, checking dimension compatibility."""
        # NOTE: THIS IS APPROUVED. DO NOT TOUCH IT FFS.
        if not isinstance(other, Tensor):
            other = scalar(other)

        dimension = op.dimension_result(self.phy_dimension, other.phy_dimension)
        if dimension is None:
            raise ValueError("Dimensions incompatible for given operator")
        
        return Tensor(
            data=op.operator(self._data, other._data), 
            dimension=dimension
        )

    @cached_property
    def scalar(self) -> TensorAsScalar:
        """Return the scalar value of this tensor if it is a scalar, otherwise raise a ValueError."""
        if self.kind != TensorKind.SCALAR:
            raise ValueError("Tensor is not a scalar.")
        
        return TensorAsScalar(self._data, self.phy_dimension)
    
    @cached_property
    def vector3(self) -> TensorAsVector3:
        """Return the vector data of this tensor if it is a vector3, otherwise raise a ValueError."""
        if self.kind != TensorKind.VECTOR3:
            raise ValueError("Tensor is not a vector3.")
        
        return TensorAsVector3(self._data, self.phy_dimension)
    
    @cached_property
    def matrix33(self) -> TensorAsMatrix33:
        """Return the matrix data of this tensor if it is a matrix33, otherwise raise a ValueError."""
        if self.kind != TensorKind.MATRIX33:
            raise ValueError("Tensor is not a matrix33.")
        
        return TensorAsMatrix33(self._data, self.phy_dimension)
    
class TensorAsScalar(Tensor):
    def __init__(self, data: NumpyFloatArray | RealNumber, dimension: type[Dim] | None = None):
        super().__init__(data, dimension)

        if self.kind != TensorKind.SCALAR:
            raise ValueError("Tensor is not a scalar.")
        
    @property
    def scalar(self) -> TensorAsScalar:
        return self
        
    def value(self, units: Tensor | str | float | None = None) -> float:
        try:
            return self.raw_data_array(units).item()
        except AttributeError:
            raise ValueError("Tensor data is not a single scalar value.")
        
    def values(self, units: Tensor | str | float | None = None) -> npt.NDArray[np.float64]:
        return self.raw_data_array(units)
    
class TensorAsVector3(Tensor):
    def __init__(self, data: NumpyFloatArray | RealNumber, dimension: type[Dim] | None = None):
        super().__init__(data, dimension)

        if self.kind != TensorKind.VECTOR3:
            raise ValueError("Tensor is not a vector3.")
        
    @property
    def vector3(self) -> TensorAsVector3:
        return self
        
    @property
    def length_squared(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return self.length ** 2
    
    @property
    def length(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.linalg.norm(self._data, axis=0),
            self.phy_dimension
        )
    
    @property
    def x(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[0]).reshape(-1), 
            self.phy_dimension
        )

    @property
    def y(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[1]).reshape(-1), 
            self.phy_dimension
        )

    @property
    def z(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[2]).reshape(-1), 
            self.phy_dimension
        )

    @property
    def theta(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, dimension=Angle)]:
        if self.kind != TensorKind.VECTOR3:
            raise RuntimeError("theta is only available for vector tensors")
        return atan2(self.y, self.x)

    @property
    def delta(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, dimension=Angle)]:
        if self.kind != TensorKind.VECTOR3:
            raise RuntimeError("delta is only available for vector tensors")
        return asin((self.z / self.length))
    
    def dot(self, 
        other: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        if self.kind != other.kind:
            raise RuntimeError("dot product is only available between vectors")
        
        return Tensor(
            np.sum(self._data * other._data, axis=0).reshape(-1),
            self.phy_dimension * other.phy_dimension
        )

    def cross(self, 
        other: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        if self.kind != other.kind:
            raise RuntimeError("cross product is only available between vectors")
        return Tensor(
            np.cross(self._data.T, other._data.T).T,
            self.phy_dimension * other.phy_dimension
        )
    
    def angle(self, 
        other: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, dimension=Angle)]:
        if self.kind != other.kind:
            raise RuntimeError("angle is only available for vector tensors")
        data_len_self = np.sum(self._data ** 2, axis=0).reshape(-1) ** .5
        data_len_other = np.sum(other._data ** 2, axis=0).reshape(-1) ** .5
        data_dot = np.sum(self._data * other._data, axis=0).reshape(-1)

        return acos(ensure_tensor(data_dot / (data_len_self * data_len_other)))

    def normalized(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3, dimension=Dimless)]:
        norm = np.linalg.norm(self._data, axis=0, keepdims=True)
        return Tensor(
            self._data / norm, 
            Dimless
        )
    
    @classmethod
    def from_components(cls, 
        x: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)], 
        y: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)], 
        z: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        tx, ty, tz = ensure_tensor(x), ensure_tensor(y), ensure_tensor(z)
        ensure_same_dimensions(tx, ty, tz)
        return cls(np.vstack((tx._data.reshape(1, -1), ty._data.reshape(1, -1), tz._data.reshape(1, -1))), tx.phy_dimension)

    @classmethod
    def from_spherical(cls, 
        theta: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, dimension=Angle)], 
        delta: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, dimension=Angle)], 
        radius: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        ensure_same_dimensions(theta, delta)
        x = radius * cos(delta) * cos(theta)
        y = radius * cos(delta) * sin(theta)
        z = radius * sin(delta)
        return cls.from_components(x, y, z)

class TensorAsMatrix33(Tensor):
    def __init__(self, data: NumpyFloatArray | RealNumber, dimension: type[Dim] | None = None):
        super().__init__(data, dimension)

        if self.kind != TensorKind.MATRIX33:
            raise ValueError("...")
        
    @property
    def matrix33(self) -> TensorAsMatrix33:
        """Return the matrix data of this tensor if it is a matrix33, otherwise raise a ValueError."""
        return self
        
    def inverse(self) -> Tensor:
        if self.kind != TensorKind.MATRIX33:
            raise RuntimeError("inverse only applies to matrix tensors")
        if self._data.ndim == 2:
            inv = np.linalg.inv(self._data)
        else:
            inv = np.stack([np.linalg.inv(self._data[:, :, i]) for i in range(self._data.shape[2])], axis=2)
        inv = np.where(np.abs(inv) < 1e-15, 0.0, inv)
        return Tensor(inv, 1 / self.phy_dimension)
    
    def matrix_product(self, other: Tensor) -> Tensor:
        if other.kind == TensorKind.SCALAR:
            return self * other
        elif other.kind == TensorKind.VECTOR3:
            return Tensor(
                np.einsum('ijk,jk->ik', self._data, other._data),
                self.phy_dimension * other.phy_dimension
            )
        elif other.kind == TensorKind.MATRIX33:
            return Tensor(
                np.einsum('ijk,jlk->ilk', self._data, other._data),
                self.phy_dimension * other.phy_dimension
            )
        
        raise RuntimeError("Unsupported tensor kind for matrix product.")
    
    def __matmul__(self, other: Tensor) -> Tensor:
        return self.matrix_product(other)
    
    @classmethod
    def from_elements(
        cls,
        a11: RealNumber | Tensor, a21: RealNumber | Tensor, a31: RealNumber | Tensor,
        a12: RealNumber | Tensor, a22: RealNumber | Tensor, a32: RealNumber | Tensor,
        a13: RealNumber | Tensor, a23: RealNumber | Tensor, a33: RealNumber | Tensor,
    ) -> Tensor:
        vals = [ensure_tensor(v) for v in (a11, a21, a31, a12, a22, a32, a13, a23, a33)]
        ensure_same_dimensions(*vals)
        dim = vals[0].phy_dimension
        arr = np.asarray([v.scalar.value() for v in vals], dtype=np.float64).reshape((3, 3, 1))
        return cls(arr, dim)

@dataclass
class TensorBound:
    dimension: type[Dim] | None = None
    kind: TensorKind | str | None = None
    size: int | None = None

    def check(self, tensor: Tensor) -> bool:
        if self.dimension is None and self.kind is None and self.size is None:
            return True
        return tensor.check(
            dimension=self.dimension,
            kind=self.kind,
            size=self.size
        )
    
    def secure(self, tensor: Tensor) -> Tensor:
        return tensor.secure(
            dimension=self.dimension,
            kind=self.kind,
            size=self.size
        )

def tensor_check(f):
    """Wrapper that checks function inputs/outputs based on type annotations.

    Supported runtime checks:
    - ``Annotated[Tensor, TensorBound(...)]``
    - ``tuple[...]`` containing supported element annotations
    - Plain runtime types via ``isinstance``

    String-literal annotations are intentionally ignored.
    """
    signature = inspect.signature(f)
    raw_annotations = getattr(f, "__annotations__", {})
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

    def _check_against_annotation(value: Any, annotation: Any, where: str) -> None:
        origin = get_origin(annotation)

        if origin is Annotated:
            bound = _bound_from_annotation(where, annotation, where)
            if not isinstance(value, Tensor):
                raise TypeError(f"{where} must be a Tensor.")
            if not bound.check(value):
                raise ValueError(f"{where} does not satisfy declared TensorBound.")
            return

        if origin is tuple:
            if not isinstance(value, tuple):
                raise TypeError(f"{where} must be a tuple.")

            tuple_annotations = get_args(annotation)

            if len(tuple_annotations) == 2 and tuple_annotations[1] is Ellipsis:
                element_annotation = tuple_annotations[0]
                for element in value:
                    _check_against_annotation(element, element_annotation, f"{where} tuple element")
                return

            if len(value) != len(tuple_annotations):
                raise TypeError(f"{where} tuple arity does not match its annotation.")

            for element, element_annotation in zip(value, tuple_annotations):
                _check_against_annotation(element, element_annotation, f"{where} tuple element")
            return

        if annotation is Any:
            return

        if isinstance(annotation, type):
            if not isinstance(value, annotation):
                raise TypeError(f"{where} must be of type {annotation.__name__}.")
            return

        # Unsupported typing forms are ignored at runtime.
        return

    parameter_annotations: dict[str, Any] = {}
    for name in signature.parameters:
        raw_annotation = raw_annotations.get(name, None)
        if isinstance(raw_annotation, str):
            continue

        annotation = hints.get(name, None)
        if annotation is None:
            continue

        parameter_annotations[name] = annotation

    raw_return_annotation = raw_annotations.get("return", None)
    return_annotation = None if isinstance(raw_return_annotation, str) else hints.get("return", None)

    @wraps(f)
    def wrapper(*args, **kwargs):
        bound_arguments = signature.bind(*args, **kwargs)
        bound_arguments.apply_defaults()

        for name, value in bound_arguments.arguments.items():
            annotation = parameter_annotations.get(name)
            if annotation is None:
                continue

            parameter = signature.parameters[name]

            if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
                values_to_check = value
            elif parameter.kind is inspect.Parameter.VAR_KEYWORD:
                values_to_check = value.values()
            else:
                values_to_check = (value,)

            for checked_value in values_to_check:
                _check_against_annotation(checked_value, annotation, f"Argument '{name}'")

        result = f(*args, **kwargs)

        if return_annotation is not None:
            _check_against_annotation(result, return_annotation, "Return value")

        return result

    return wrapper

_UNITS_REGISTRY: dict[str, Tensor] = dict()

def _units_register(units: list[str], dim: type[Dim], base_factor: float):
    global _UNITS_REGISTRY
    _UNITS_REGISTRY |= {
        u: Tensor(
            data=np.asarray(base_factor, dtype=np.float64).reshape((1,)), 
            dimension=dim
        )
        for u in units
    }

class QuantityMeta(type):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """

    def __getattr__(cls, name: str) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        """Resolve a unit name into its corresponding scalar quantity."""
        try:
            global _UNITS_REGISTRY
            return _UNITS_REGISTRY[name].copy()
        except KeyError:
            raise ValueError(f"No unit named '{name}'")


class Quantity(metaclass=QuantityMeta):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """
    @classmethod
    def get(cls, value: str) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        """Return the scalar unit associated with ``value``."""
        return cls.__getattr__(value)

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
    _units_register(["kilo_meter"], Length, 1e3)

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
    _units_register(["kilo_meter_per_hour", "kmph"], Velocity, 1000 / 3600)
    _units_register(["radian_per_second", "rad_per_s"], AngularVelocity, 1.0)

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

    return Tensor(data=np.asarray(value, dtype=np.float64).reshape((1,)), dimension=Dimless)

def vector3(x: RealNumber, y: RealNumber, z: RealNumber) -> Tensor:
    """Create a dimensionless vector3 tensor with the given values."""
    return Tensor(
        data=np.asarray((x, y, z), dtype=np.float64).reshape((3, 1)), 
        dimension=Dimless
    )


def matrix33(a11: RealNumber, a12: RealNumber, a13: RealNumber,
          a21: RealNumber, a22: RealNumber, a23: RealNumber,
          a31: RealNumber, a32: RealNumber, a33: RealNumber) -> Tensor:
    """Create a dimensionless matrix33 tensor with the given values.
    Elements given row-wise."""
    return Tensor(
        data=np.asarray((
            a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33
        ), dtype=np.float64).reshape((3, 3, 1)), 
        dimension=Dimless
    )


def ensure_tensor(v: Tensor | RealNumber) -> Tensor:
    if isinstance(v, Tensor):
        return v
    if isinstance(v, np.ndarray):
        arr = np.asarray(v, dtype=np.float64)
        return Tensor(arr.reshape(-1), Dimless)
    return scalar(v)


def ensure_same_dimensions(*tensors: Tensor) -> bool:
    if len(tensors) < 2:
        return True
    d0 = tensors[0].phy_dimension.triplet()
    if not all(t.phy_dimension.triplet() == d0 for t in tensors):
        raise RuntimeError("Dimensions are not compatible")
    return True

@tensor_check
def sin(
    t: Annotated[Tensor, TensorBound(dimension=Angle)]
) -> Annotated[Tensor, TensorBound(dimension=Dimless)]:
    return Tensor(np.sin(t._data), Dimless)


@tensor_check
def cos(
    t: Annotated[Tensor, TensorBound(dimension=Angle)]
) -> Annotated[Tensor, TensorBound(dimension=Dimless)]:
    return Tensor(np.cos(t._data), Dimless)


@tensor_check
def tan(
    t: Annotated[Tensor, TensorBound(dimension=Angle)]
) -> Annotated[Tensor, TensorBound(dimension=Dimless)]:
    return Tensor(np.tan(t._data), Dimless)


@tensor_check
def asin(
    t: Annotated[Tensor, TensorBound(dimension=Dimless)]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    return Tensor(np.arcsin(t._data), Angle)


@tensor_check
def acos(
    t: Annotated[Tensor, TensorBound(dimension=Dimless)]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    return Tensor(np.arccos(t._data), Angle)


@tensor_check
def atan(
    t: Annotated[Tensor, TensorBound(dimension=Dimless)]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    return Tensor(np.arctan(t._data), Angle)


@tensor_check
def atan2(
    y: Annotated[Tensor, TensorBound()],
    x: Annotated[Tensor, TensorBound()]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    ensure_same_dimensions(y, x)
    return Tensor(np.arctan2(y._data, x._data), Angle)


@tensor_check
def normalize_angle(
    angle: Annotated[Tensor, TensorBound(dimension=Angle)]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    return Tensor(np.mod(angle._data, 2 * np.pi), Angle)


@tensor_check
def normalize_angle_symmetric(
    angle: Annotated[Tensor, TensorBound(dimension=Angle)]
) -> Annotated[Tensor, TensorBound(dimension=Angle)]:
    wrapped = np.mod(angle._data + np.pi, 2 * np.pi) - np.pi
    return Tensor(wrapped, Angle)


@tensor_check
def interpolate(
    a: Annotated[Tensor, TensorBound()],
    b: Annotated[Tensor, TensorBound()],
    p: float
) -> Annotated[Tensor, TensorBound()]:
    ensure_same_dimensions(a, b)
    return Tensor(a._data + p * (b._data - a._data), a.phy_dimension)