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

import pint

type NumpyFloatArray = npt.NDArray[np.float64]
RealNumber = float | int | Fraction | Decimal

U = pint.UnitRegistry()
"""`pint.UnitRegistry` instance for unit management."""

class TensorBinaryOperator(Enum):
    """Supported binary operators between tensors and scalar numbers."""

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
    
    def dimension_result(self, left_dim: pint.Unit, right_dim: pint.Unit) -> pint.Unit | None:
        """Return the resulting dimension of applying this operator to quantities of the given dimensions."""
        if self in (TensorBinaryOperator.ADD, TensorBinaryOperator.SUB):
            if not left_dim.is_compatible_with(right_dim):
                return None
            return left_dim
        elif self == TensorBinaryOperator.MODULO:
            if not left_dim.is_compatible_with(right_dim) and not right_dim.dimensionless:
                return None
            return left_dim
        elif self in (TensorBinaryOperator.MUL, ):
            return left_dim * right_dim
        elif self in (TensorBinaryOperator.TRUEDIV, TensorBinaryOperator.FLOORDIV):
            return left_dim / right_dim
        
        raise NotImplementedError(f"Unsupported operator '{self.value}' for dimension composition.")

class TensorKind(Enum):
    """Shape category of a tensor: scalar, 3-vector, or 3x3 matrix."""

    SCALAR = "scalar"
    VECTOR3 = "vector3"
    MATRIX33 = "matrix33"


class Tensor:
    """Generic n-dimensional tensor carrying a physical dimension.
    """
    _data: NumpyFloatArray
    _units: pint.Unit

    def __init__(self, data: NumpyFloatArray | RealNumber, units: pint.Unit | None = None):
        """Initialize a tensor with the given data and dimension.
        Data units are default SI units corresponding to dimension."""
        if units is None:
            units = U.dimensionless
        elif not isinstance(units, pint.Unit):
            units = U.parse_units(str(units))
        if not isinstance(data, np.ndarray):
            data = np.asarray(data, dtype=np.float64).reshape((1,))

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
        self._units = cast(pint.Unit, units)

    @property
    def data(self) -> NumpyFloatArray:
        return self._data.copy()

    def __hash__(self) -> int:
        return hash(f"{hash(self._data.data.tobytes())}${hash(self._units.dimensionality)}")
    
    def __repr__(self) -> str:
        return self.human_repr()
    
    def human_repr(self, units: pint.Unit | str | None = None) -> str:
        """Human-readable inline representation.

        Uses ``<Scalar ...>``, ``<Vector3 ...>``, or ``<Matrix33 ...>``.
        For tensors with multiple samples, only the first sample is shown and
        the output is suffixed with `...`.
        """

        if units is None:
            units = _get_base_units(self._units)

        data = self.raw_data_array(units)
        dim_repr = units if isinstance(units, str) else str(self._units)
        suffix = " ..." if self.size > 1 else ""

        if self.kind == TensorKind.SCALAR:
            shown_data = np.asarray(data[0]).reshape((1,))
            data_repr = np.array2string(
                shown_data,
                separator=", ",
                max_line_width=10_000,
            )
            data_repr = " ".join(data_repr.split())
            return f"<Scalar {data_repr} [{dim_repr}]{suffix}>"
        elif self.kind == TensorKind.VECTOR3:
            shown_data = np.asarray(data[:, 0]).reshape((3,))
            x, y, z = (float(shown_data[0]), float(shown_data[1]), float(shown_data[2]))
            return f"<Vector3 x={x} y={y} z={z} [{dim_repr}]{suffix}>"
        else:
            shown_data = np.asarray(data[:, :, 0]).reshape((3, 3))
            return (
                "<Matrix33 "
                f"a11={float(shown_data[0, 0])} a12={float(shown_data[0, 1])} a13={float(shown_data[0, 2])} "
                f"a21={float(shown_data[1, 0])} a22={float(shown_data[1, 1])} a23={float(shown_data[1, 2])} "
                f"a31={float(shown_data[2, 0])} a32={float(shown_data[2, 1])} a33={float(shown_data[2, 2])} "
                f"[{dim_repr}]{suffix}>"
            )
    
    def __copy__(self) -> Tensor:
        """Return a shallow/deep copy of this tensor."""
        return self.copy()

    def copy(self) -> Tensor:
        return Tensor(
            data=self._data.copy(),
            units=copy(self._units)
        )
    
    @property
    def units(self) -> pint.Unit:
        """Return the physical dimension of this tensor."""
        return copy(self._units)
    
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
        units: pint.Unit | None = None,
        kind: TensorKind | str | None = None,
        size: int | None = None
    ) -> bool:
        """Check if the dimension of this tensor matches the given one."""
        if (units == kind == size == None):
            return True
        
        return (
            (self._units.is_compatible_with(units) if units is not None else True)
            and (self.kind == TensorKind(kind) if kind is not None else True)
            and (self.size == size if size is not None else True)
        )
    
    def secure(self, 
        units: pint.Unit | None = None,
        kind: TensorKind | str | None = None,
        size: int | None = None
    ) -> Self:
        """Returns self if the dimension of this tensor matches the given one, otherwise raises a ValueError."""
        if not self.check(units, kind, size):
            raise ValueError(f"Tensor does not match prerogatives.")
        
        return self
    
    def concatenate(self, other: Tensor) -> Tensor:
        """Concatenate this tensor with another one along the first axis, checking dimension compatibility."""
        if self.kind != other.kind:
            raise ValueError("Cannot concatenate tensors of different kinds.")
        if not self._units.is_compatible_with(other._units):
            raise ValueError("Cannot concatenate tensors of different dimensions.")
        return Tensor(
            np.concatenate((self._data, other._data), axis=self.array_dimensions - 1),
            self._units
        )
    
    def raw_data_array(self, units: pint.Unit | str) -> NumpyFloatArray:
        """Return the raw data array of this tensor, converted to the given units."""
        factor: float | None = None

        if isinstance(units, str):
            units = U.parse_units(units)
        if isinstance(units, pint.Unit):
            if not self._units.is_compatible_with(units):
                raise ValueError("Units must be compatible with the tensor's physical dimension.")
            factor = cast(pint.Quantity, 1.0 * self._units).m_as(units)
        
        if factor is None:
            raise ValueError("Unsupported units specification.")
        
        return self._data * factor
    
    def with_units(self, units: pint.Unit | Tensor | str) -> Tensor:
        """Returns a copy of this `Tensor` object with the same data
        but with the dimension of the provided units"""

        return Tensor(
            data=self._data.copy(),
            units=_retrieve_units(units)
        )
    
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
            copy(self._units)
        )

        assert tensor.array_dimensions == self.array_dimensions
        assert tensor.size == 1

        return tensor

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Tensor):
            return False
        
        return (
            self._units.is_compatible_with(o._units)
            and np.array_equal(self._data, o._data)
        )

    def __pos__(self) -> Tensor:
        return Tensor(
            data=+self._data,
            units=copy(self._units)
        )
    
    def __neg__(self) -> Tensor:
        return Tensor(
            data=-self._data,
            units=copy(self._units)
        )
    
    def __abs__(self) -> Tensor:
        return Tensor(
            data=np.abs(self._data),
            units=copy(self._units)
        )
    
    def __pow__(self, exponent: RealNumber | int | float) -> Tensor:
        return Tensor(
            data=self._data ** float(exponent),
            units=self._units ** exponent # type: ignore (it works...)
        )

    def __lt__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() < 0)

    def __le__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() <= 0)

    def __gt__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() > 0)

    def __ge__(self, other: Tensor) -> bool:
        return bool((self._data - other._data).item() >= 0)

    def __add__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.ADD)
    
    def __radd__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.ADD)
    
    def __sub__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.SUB)
    
    def __rsub__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.SUB)
    
    def __mul__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.MUL)
    
    def __rmul__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.MUL)
    
    def __truediv__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.TRUEDIV)
    
    def __rtruediv__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.TRUEDIV)
    
    def __floordiv__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.FLOORDIV)
    
    def __rfloordiv__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.FLOORDIV)
    
    def __mod__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return self.perform_binary_operation(right, TensorBinaryOperator.MODULO)
    
    def __rmod__(self, left: RealNumber | pint.Quantity | pint.Unit) -> Tensor:
        return scalar(left).perform_binary_operation(self, TensorBinaryOperator.MODULO)
    
    def __or__(self, right: Tensor) -> Tensor:
        """Concatenation operator, only works for tensors of the same kind and dimension."""
        return self.concatenate(right)
    
    def perform_binary_operation(self, other: Tensor | RealNumber | pint.Quantity | pint.Unit, op: TensorBinaryOperator) -> Tensor:
        """Perform the given binary operation with another tensor, checking dimension compatibility."""
        # NOTE: THIS IS APPROUVED. DO NOT TOUCH IT FFS.
        if not isinstance(other, Tensor):
            other = scalar(other)

        units = op.dimension_result(self._units, other._units)
        if units is None:
            raise ValueError("Dimensions incompatible for given operator")
        
        return Tensor(
            data=op.operator(self._data, other._data), 
            units=units
        )

    @cached_property
    def scalar(self) -> TensorAsScalar:
        """Return the scalar value of this tensor if it is a scalar, otherwise raise a ValueError."""
        if self.kind != TensorKind.SCALAR:
            raise ValueError("Tensor is not a scalar.")
        
        return TensorAsScalar(self._data, self._units)
    
    @cached_property
    def vector3(self) -> TensorAsVector3:
        """Return the vector data of this tensor if it is a vector3, otherwise raise a ValueError."""
        if self.kind != TensorKind.VECTOR3:
            raise ValueError("Tensor is not a vector3.")
        
        return TensorAsVector3(self._data, self._units)
    
    @cached_property
    def matrix33(self) -> TensorAsMatrix33:
        """Return the matrix data of this tensor if it is a matrix33, otherwise raise a ValueError."""
        if self.kind != TensorKind.MATRIX33:
            raise ValueError("Tensor is not a matrix33.")
        
        return TensorAsMatrix33(self._data, self._units)
    
class TensorAsScalar(Tensor):
    """Scalar-specialized tensor helper."""

    def __init__(self, data: NumpyFloatArray | RealNumber, units: pint.Unit | None = None):
        super().__init__(data, units)

        if self.kind != TensorKind.SCALAR:
            raise ValueError("Tensor is not a scalar.")
        
    @property
    def scalar(self) -> TensorAsScalar:
        """Return this scalar view itself."""
        return self
        
    def value(self, units: pint.Unit | str) -> float:
        """Return the scalar value converted to optional ``units``."""
        try:
            return self.raw_data_array(units).item()
        except AttributeError:
            raise ValueError("Tensor data is not a single scalar value.")
        
    def values(self, units: pint.Unit | str) -> npt.NDArray[np.float64]:
        """Return scalar data as a NumPy array converted to optional ``units``."""
        return self.raw_data_array(units)
    
class TensorAsVector3(Tensor):
    """Vector3-specialized tensor helper with vector operations."""

    def __init__(self, data: NumpyFloatArray | RealNumber, units: pint.Unit | None = None):
        super().__init__(data, units)

        if self.kind != TensorKind.VECTOR3:
            raise ValueError("Tensor is not a vector3.")
        
    @property
    def vector3(self) -> TensorAsVector3:
        """Return this vector3 view itself."""
        return self
        
    @property
    def length_squared(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return self.length ** 2
    
    @property
    def length(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.linalg.norm(self._data, axis=0),
            self._units
        )
    
    @property
    def x(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[0]).reshape(-1), 
            self._units
        )

    @property
    def y(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[1]).reshape(-1), 
            self._units
        )

    @property
    def z(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]:
        return Tensor(
            np.asarray(self._data[2]).reshape(-1), 
            self._units
        )

    @property
    def theta(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, units=U.radian)]:
        if self.kind != TensorKind.VECTOR3:
            raise RuntimeError("theta is only available for vector tensors")
        return atan2(self.y, self.x)

    @property
    def delta(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, units=U.radian)]:
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
            self._units * other._units
        )

    def cross(self, 
        other: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        if self.kind != other.kind:
            raise RuntimeError("cross product is only available between vectors")
        return Tensor(
            np.cross(self._data.T, other._data.T).T,
            self._units * other._units
        )
    
    def angle(self, 
        other: Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, units=U.radian)]:
        if self.kind != other.kind:
            raise RuntimeError("angle is only available for vector tensors")
        data_len_self = np.sum(self._data ** 2, axis=0).reshape(-1) ** .5
        data_len_other = np.sum(other._data ** 2, axis=0).reshape(-1) ** .5
        data_dot = np.sum(self._data * other._data, axis=0).reshape(-1)

        return acos(ensure_tensor(data_dot / (data_len_self * data_len_other)))

    def normalized(self) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3, units=U.dimensionless)]:
        norm = np.linalg.norm(self._data, axis=0, keepdims=True)
        return Tensor(
            self._data / norm, 
            U.dimensionless
        )
    
    @classmethod
    def from_components(cls, 
        x: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)], 
        y: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)], 
        z: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        tx, ty, tz = ensure_tensor(x), ensure_tensor(y), ensure_tensor(z)
        ensure_same_dimensions(tx, ty, tz)
        return cls(np.vstack((tx._data.reshape(1, -1), ty._data.reshape(1, -1), tz._data.reshape(1, -1))), tx._units)

    @classmethod
    def from_spherical(cls, 
        theta: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, units=U.radian)], 
        delta: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR, units=U.radian)], 
        radius: Annotated[Tensor, TensorBound(kind=TensorKind.SCALAR)]
    ) -> Annotated[Tensor, TensorBound(kind=TensorKind.VECTOR3)]:
        ensure_same_dimensions(theta, delta)
        x = radius * cos(delta) * cos(theta)
        y = radius * cos(delta) * sin(theta)
        z = radius * sin(delta)
        return cls.from_components(x, y, z)

class TensorAsMatrix33(Tensor):
    """Matrix33-specialized tensor helper with linear algebra utilities."""

    def __init__(self, data: NumpyFloatArray | RealNumber, units: pint.Unit | None = None):
        super().__init__(data, units)

        if self.kind != TensorKind.MATRIX33:
            raise ValueError("...")
        
    @property
    def matrix33(self) -> TensorAsMatrix33:
        """Return the matrix data of this tensor if it is a matrix33, otherwise raise a ValueError."""
        return self
        
    def inverse(self) -> Tensor:
        """Return matrix inverse for each matrix sample in the tensor."""
        if self.kind != TensorKind.MATRIX33:
            raise RuntimeError("inverse only applies to matrix tensors")
        if self._data.ndim == 2:
            inv = np.linalg.inv(self._data)
        else:
            inv = np.stack([np.linalg.inv(self._data[:, :, i]) for i in range(self._data.shape[2])], axis=2)
        inv = np.where(np.abs(inv) < 1e-15, 0.0, inv)
        return Tensor(inv, 1 / self._units)
    
    def matrix_product(self, other: Tensor) -> Tensor:
        """Multiply by a scalar, vector3, or matrix33 tensor."""
        if other.kind == TensorKind.SCALAR:
            return self * other
        elif other.kind == TensorKind.VECTOR3:
            return Tensor(
                np.einsum('ijk,jk->ik', self._data, other._data),
                self._units * other._units
            )
        elif other.kind == TensorKind.MATRIX33:
            return Tensor(
                np.einsum('ijk,jlk->ilk', self._data, other._data),
                self._units * other._units
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
        """Build a matrix33 tensor from nine scalar elements."""
        vals = [ensure_tensor(v) for v in (a11, a21, a31, a12, a22, a32, a13, a23, a33)]
        ensure_same_dimensions(*vals)
        dim = vals[0]._units
        arr = np.asarray([v.scalar.value(vals[0].units) for v in vals], dtype=np.float64).reshape((3, 3, 1))
        return cls(arr, dim)

@dataclass
class TensorBound:
    """Runtime constraints describing acceptable tensor dimension/shape/size."""

    units: pint.Unit | str | None = None
    kind: TensorKind | str | None = None
    size: int | None = None

    def check(self, tensor: Tensor) -> bool:
        """Return whether ``tensor`` satisfies this bound."""
        if self.units is None and self.kind is None and self.size is None:
            return True
        
        if isinstance(self.units, str):
            try:
                units = U.parse_units(self.units)
            except pint.UndefinedUnitError:
                raise ValueError(f"Invalid units string '{self.units}' in TensorBound.")
            self.units = units

        return tensor.check(
            units=self.units,
            kind=self.kind,
            size=self.size
        )
    
    def secure(self, tensor: Tensor) -> Tensor:
        """Return ``tensor`` when it satisfies this bound, else raise."""

        if isinstance(self.units, str):
            try:
                units = U.parse_units(self.units)
            except pint.UndefinedUnitError:
                raise ValueError(f"Invalid units string '{self.units}' in TensorBound.")
            self.units = units
        
        return tensor.secure(
            units=self.units,
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

def _retrieve_units(obj: pint.Unit | str | Tensor) -> pint.Unit:
    if isinstance(obj, Tensor):
        return obj._units
    elif isinstance(obj, str):
        return U.parse_units(obj)
    elif isinstance(obj, pint.Unit):
        return obj
    else:
        # Support concrete Pint unit implementations not typed as ``pint.Unit``.
        try:
            return U.parse_units(str(obj))
        except Exception:
            pass
    
    raise TypeError()

def _get_base_units(q: pint.Unit | pint.Quantity) -> pint.Unit:
    if isinstance(q, pint.Unit):
        q = cast(pint.Quantity, 1.0 * q)
    if isinstance(q, pint.Quantity):
        return cast(pint.Unit, q.to_base_units().units)
    
    raise TypeError()

### CONVENIENCE FACTORY FUNCTIONS ###
### ############################# ###

def scalar(value: RealNumber | pint.Unit | pint.Quantity | str) -> Tensor:
    """Create a scalar tensor from a number, unit, or pint quantity."""
    if isinstance(value, str):
        # Parse textual quantities like "12 meter" through this registry.
        value = cast(pint.Quantity, U.Quantity(value))
    
    if isinstance(value, pint.Unit):
        return Tensor(
            data=np.asarray(1.0, dtype=np.float64).reshape((1,)),
            units=value,
        )

    if isinstance(value, pint.Quantity):
        return Tensor(
            data=np.asarray(value.magnitude, dtype=np.float64).reshape((1,)),
            units=_retrieve_units(value.units), # type: ignore ...
        )

    # Accept pint-like quantity objects from other runtimes/registries.
    if hasattr(value, "magnitude") and hasattr(value, "units"):
        q = U.Quantity(getattr(value, "magnitude"), getattr(value, "units"))
        return Tensor(
            data=np.asarray(q.magnitude, dtype=np.float64).reshape((1,)),
            units=_retrieve_units(q.units), # type: ignore ...
        )

    if not isinstance(value, RealNumber):
        raise TypeError(f"Unsupported scalar input type: {type(value)!r}")

    return Tensor(
        data=np.asarray(value, dtype=np.float64).reshape((1,)),
        units=U.dimensionless
    )

def vector3(x: RealNumber, y: RealNumber, z: RealNumber) -> Tensor:
    """Create a dimensionless vector3 tensor with the given values."""
    return Tensor(
        data=np.asarray((x, y, z), dtype=np.float64).reshape((3, 1)), 
        units=U.dimensionless
    )

def matrix33(a11: RealNumber, a12: RealNumber, a13: RealNumber,
          a21: RealNumber, a22: RealNumber, a23: RealNumber,
          a31: RealNumber, a32: RealNumber, a33: RealNumber) -> Tensor:
    """Create a dimensionless matrix33 tensor with the given values.
    Elements given row-wise."""
    return Tensor(
        data=np.asarray(
            (
                a11, a12, a13,
                a21, a22, a23,
                a31, a32, a33
            ), 
            dtype=np.float64
        ).reshape((3, 3, 1)), 
        units=U.dimensionless
    )


def ensure_tensor(v: Tensor | RealNumber | pint.Unit | pint.Quantity) -> Tensor:
    """Return ``v`` as a tensor, wrapping plain numbers/arrays as dimensionless."""
    if isinstance(v, Tensor):
        return v
    if isinstance(v, np.ndarray):
        arr = np.asarray(v, dtype=np.float64)
        return Tensor(
            arr.reshape(-1), 
            units=U.dimensionless
        )
    return scalar(v)


def ensure_same_dimensions(*tensors: Tensor) -> bool:
    """Validate that all tensors share the same physical dimension."""
    if len(tensors) < 2:
        return True
    d0 = tensors[0]._units
    if not all(d0.is_compatible_with(t._units) for t in tensors):
        raise RuntimeError("Dimensions are not compatible")
    return True

@tensor_check
def sin(
    t: Annotated[Tensor, TensorBound(units=U.radian)]
) -> Annotated[Tensor, TensorBound(units=U.dimensionless)]:
    """Element-wise sine on angle tensors."""
    return Tensor(np.sin(t._data), U.dimensionless)


@tensor_check
def cos(
    t: Annotated[Tensor, TensorBound(units=U.radian)]
) -> Annotated[Tensor, TensorBound(units=U.dimensionless)]:
    """Element-wise cosine on angle tensors."""
    return Tensor(np.cos(t._data), U.dimensionless)


@tensor_check
def tan(
    t: Annotated[Tensor, TensorBound(units=U.radian)]
) -> Annotated[Tensor, TensorBound(units=U.dimensionless)]:
    """Element-wise tangent on angle tensors."""
    return Tensor(np.tan(t._data), U.dimensionless)


@tensor_check
def asin(
    t: Annotated[Tensor, TensorBound(units=U.dimensionless)]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Element-wise arcsine returning angle tensors."""
    return Tensor(np.arcsin(t._data), U.radian)


@tensor_check
def acos(
    t: Annotated[Tensor, TensorBound(units=U.dimensionless)]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Element-wise arccosine returning angle tensors."""
    return Tensor(np.arccos(t._data), U.radian)


@tensor_check
def atan(
    t: Annotated[Tensor, TensorBound(units=U.dimensionless)]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Element-wise arctangent returning angle tensors."""
    return Tensor(np.arctan(t._data), U.radian)


@tensor_check
def atan2(
    y: Annotated[Tensor, TensorBound()],
    x: Annotated[Tensor, TensorBound()]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Element-wise two-argument arctangent with dimension checking."""
    ensure_same_dimensions(y, x)
    return Tensor(np.arctan2(y._data, x._data), U.radian)


@tensor_check
def normalize_angle(
    angle: Annotated[Tensor, TensorBound(units=U.radian)]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Wrap angles to the ``[0, 2π)`` interval."""
    return Tensor(np.mod(angle._data, 2 * np.pi), U.radian)


@tensor_check
def normalize_angle_symmetric(
    angle: Annotated[Tensor, TensorBound(units=U.radian)]
) -> Annotated[Tensor, TensorBound(units=U.radian)]:
    """Wrap angles to the ``[-π, π)`` interval."""
    wrapped = np.mod(angle._data + np.pi, 2 * np.pi) - np.pi
    return Tensor(wrapped, U.radian)


@tensor_check
def interpolate(
    a: Annotated[Tensor, TensorBound()],
    b: Annotated[Tensor, TensorBound()],
    p: float
) -> Annotated[Tensor, TensorBound()]:
    """Linear interpolation between tensors ``a`` and ``b`` at ratio ``p``."""
    ensure_same_dimensions(a, b)
    return Tensor(a._data + p * (b._data - a._data), a._units)