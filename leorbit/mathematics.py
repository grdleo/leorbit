from fractions import Fraction
from functools import cached_property
import inspect
from typing import Any, Callable, ClassVar, Generic, Literal, Self, Type, TypeAlias, TypeIs, TypeVar, cast

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

    @classmethod
    def triplet(cls) -> DimTriplet:
        """Return the coordinate triplet associated with this dimension."""
        if not hasattr(cls, "__triplet"):
            raise RuntimeError(f"`{cls.__name__}` is not a registered dimension class.")
        
        return getattr(cls, "__triplet")
    
Number: TypeAlias = float | int | np.floating[Any]
TensorData: TypeAlias = npt.NDArray[np.float64]
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
Numerator = TypeVar("Numerator", bound=int)
Denominator = TypeVar("Denominator", bound=int)
    
class ProductDim(Generic[SomeDim, SomeOtherDim], Dim):
    """Type-level marker representing a product of two dimensions."""
    ...

class QuotientDim(Generic[SomeDim, SomeOtherDim], Dim):
    """Type-level marker representing a quotient of two dimensions."""
    ...
    
class PowerDim(Generic[SomeDim, Numerator, Denominator], Dim):
    """Type-level marker representing a powered dimension."""
    ...

_ = Dimless = Angle = _dimension_factory(DimTriplet())
_ = Length = _dimension_factory(DimTriplet(length=1))
_ = Time = _dimension_factory(DimTriplet(time=1))
_ = Mass = _dimension_factory(DimTriplet(mass=1))
_ = Velocity = Length / Time
_ = Acceleration = Velocity / Time
_ = Force = Mass * Acceleration
_ = Frequency = AngularVelocity = 1 / Time
# NOTE: This is a trick discovered accidentally for Pyright to recognize these dimensions as "real types"
# instead of just `type[Dim]` which would be the case if we directly assigned the result of the operations to the variables.

__UNITS_REGISTRY: dict[type[Dim], dict[str, float]] = {
    Dimless: dict(
        dimensionless=1,
        radian=1,
        turn=1 / (2 * np.pi),
        degree=np.pi / 180
    ),
    Length: dict(
        meter=1,
        kilo_meter=1e3,
        radii_earth=6378135,
        radii_sun=6.957e8,
        astronomical_unit=149597870700,
    ),
    Time: dict(
        second=1,
        minute=60,
        hour=60 * 60,
        day=24 * 60 * 60,
        month=30 * 24 * 60 * 60,
        year=365 * 24 * 60 * 60
    ),
    Mass: dict(
        kilo_gram=1,
        gram=1e-3,
        metric_ton=1e3,
    )
}

__UNITS_FACTORS_DIMENSIONS: dict[str, tuple[type[Dim], float]] = {
    unit: (dim, factor)
    for dim, units in __UNITS_REGISTRY.items()
    for unit, factor in units.items()
}


def _get_registered_unit(unit_name: str) -> tuple[type[Dim], float]:
    """Return the dimension and factor associated with a registered unit name."""
    for dim, units in __UNITS_REGISTRY.items():
        factor = units.get(unit_name, None)
        if factor is not None:
            return dim, factor

    raise ValueError(f"No unit named '{unit_name}'")


class Tensor(Generic[SomeDim]):
    """Generic n-dimensional tensor carrying a physical dimension.

    Concrete subclasses specialize tensor shape and semantics while reusing
    dimension-aware arithmetic from this base class.
    """

    _dim: ClassVar[type[Dim]]
    _base_tensor_class: ClassVar[type[Tensor]]

    def __init__(self, values: Number | TensorData):
        """Initialize raw tensor values.

        NOTE: `Tensor` should never be used as-is.
        """
        if not self._is_dimensionalized():
            raise RuntimeError("Class has to be dimensionalized")

        self._values = np.asarray(values, dtype=np.float64)

    def __hash__(self) -> int:
        return hash(f"{hash(self._values.data.tobytes())}${hash(self.dim.triplet())}")

    @property
    def dim_triplet(self) -> DimTriplet:
        """Dimension coordinates associated with this tensor."""
        return self.dim.triplet()

    @property
    def dim(self) -> type[SomeDim]:
        """Dimension class associated with this tensor."""
        return cast(type[SomeDim], self._dim)

    @classmethod
    def _is_base_tensor_class(cls) -> bool:
        """Whether `cls` is the undimensionalized root tensor class."""
        return not hasattr(cls, "_base_tensor_class")

    @classmethod
    def _is_dimensionalized(cls) -> bool:
        dim = getattr(cls, "_dim", None)
        if dim is None:
            return False

        if not issubclass(dim, Dim):
            raise RuntimeError("Invalid tensor dimension class")

        return True

    def __repr__(self) -> str:
        return f"<Tensor {self._values} [{self.dim_triplet.representation}]>"

    @classmethod
    def __class_getitem__(cls, dim: type[SomeDim]) -> type[Tensor]:
        """Return a dimensionalized tensor class for a concrete dimension."""
        if inspect.isclass(dim) and issubclass(dim, Dim):
            if cls._is_dimensionalized():
                raise RuntimeError("Cannot subscript a dimensionalized tensor class.")

            class DimensionalizedTensor(cls):
                _dim = dim

            return cast(type[Tensor], DimensionalizedTensor)

        return cls  # type: ignore[return-value]

    def cast(self, dim: type[SomeOtherDim]) -> Tensor[SomeOtherDim]:
        """Type-cast to another dimension if coordinates are identical."""
        if dim.triplet() == self.dim_triplet:
            return self  # type: ignore[return-value]

        raise RuntimeError("Cannot cast")

    def copy(self) -> Self:
        """Return a value copy with the same tensor class and dimension."""
        return self.__class__(self._values)

    def ensure_compatible_dimensions(self: Tensor[SomeDim], o: Tensor[SomeOtherDim]) -> TypeIs[Tensor[SomeOtherDim]]:
        """Return whether two tensors have identical dimension coordinates."""
        return isinstance(o, Tensor) and o.dim_triplet == self.dim_triplet

    def check(self, dim: type[Dim]) -> bool:
        """Return `True` if tensor dimension matches `dim`."""
        return self.dim_triplet == dim.triplet()

    def get_raw_array(self, units: str = "1") -> npt.NDArray[np.float64]:
        """Return values converted to requested units."""
        a = np.copy(self._values)
        if units == "1":
            return a

        dim, factor = _get_registered_unit(units)
        if dim.triplet() != self.dim_triplet:
            raise ValueError(
                f"Units '{units}' have dimension '{dim.triplet().representation}' which is incompatible "
                f"with tensor dimension '{self.dim_triplet.representation}'"
            )

        return a / factor

    def __pos__(self) -> Self:
        return self.__class__(self._values)

    def __neg__(self) -> Self:
        return self.__class__(-self._values)

    def __add__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("Incompatible dimensions")

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["+"]

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor addition")

        return tensor_cls_result[self.dim](self._values + o._values)  # type: ignore[index]

    def __radd__(self, o: object) -> Tensor:
        return ensure_tensor(o) + self

    def __sub__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("Incompatible dimensions")

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["-"]

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor subtraction")

        return tensor_cls_result[self.dim](self._values - o._values)  # type: ignore[index]

    def __rsub__(self, o: object) -> Tensor:
        return ensure_tensor(o) - self

    def __mul__(self, o: object) -> Tensor:
        o = ensure_tensor(o)

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["*"]

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor multiplication")

        return tensor_cls_result[self.dim * o.dim](self._values * o._values)  # type: ignore[index]

    def __rmul__(self, o: object) -> Tensor:
        return ensure_tensor(o) * self

    def __truediv__(self, o: object) -> Tensor:
        o = ensure_tensor(o)

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["/"]

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor division")

        return tensor_cls_result[self.dim / o.dim](self._values / o._values)  # type: ignore[index]

    def __rtruediv__(self, o: object) -> Tensor:
        return ensure_tensor(o) / self

    def __matmul__(self, o: object) -> Tensor:
        o = ensure_tensor(o)

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["@"]

        ein_op = "ij,jn->in"
        if (
            self._base_tensor_class is Tensor_M33
            and o._base_tensor_class is Tensor_V3
            and self._values.ndim == 3
        ):
            _, _, nb_matrices = self._values.shape
            _, nb_vecs = o._values.shape
            if nb_matrices == nb_vecs:
                ein_op = "abn,bn->an"
            else:
                raise RuntimeError("Incompatible batched matrix/vector shapes")

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor matrix multiplication")

        return tensor_cls_result[self.dim * o.dim](  # type: ignore[index]
            np.einsum(ein_op, self._values, o._values)
        )

    def __rmatmul__(self, o: object) -> Tensor:
        return ensure_tensor(o) @ self

    def __mod__(self, o: object) -> Tensor:
        o = ensure_tensor(o)
        if not self.ensure_compatible_dimensions(o):
            raise ValueError("Incompatible dimensions")

        tensor_cls_result = OPERATION_RESULT_TYPE[
            self._base_tensor_class, o._base_tensor_class
        ]["%"]

        if tensor_cls_result is None:
            raise ValueError("Unsupported tensor modulo")

        return tensor_cls_result[self.dim](self._values % o._values)  # type: ignore[index]

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


def ensure_tensor(o: Any | Tensor[SomeDim]) -> Tensor[Any]:
    """Return `o` as a tensor, wrapping numbers/arrays as dimensionless tensors."""
    if isinstance(o, Tensor):
        return o
    if isinstance(o, (float, int, np.floating)):
        return Scalar[Dimless](np.asarray(o))
    if isinstance(o, np.ndarray):
        arr = np.asarray(o)

        if arr.ndim == 0:
            return Scalar[Dimless](arr.item())
        if arr.ndim == 1:
            return ScalarArray[Dimless](arr)
        if arr.ndim == 2 and arr.shape[0] == 3:
            return Vector3Array[Dimless](arr)
        if arr.ndim == 2 and arr.shape == (3, 3):
            return Matrix33[Dimless](arr)

        raise RuntimeError("Unsupported array shape for tensor conversion")

    raise RuntimeError("Cannot convert object to tensor")


def ensure_same_dimensions(*tensors: Tensor[Any]) -> Literal[True]:
    """Validate that all tensors share the same dimension coordinates."""
    if len(tensors) == 0:
        return True

    t0, *_ = tensors
    if all(t.dim_triplet == t0.dim_triplet for t in tensors):
        return True

    raise RuntimeError("Tensors must share the same dimensions")


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
        if dim.triplet() == self.dim_triplet:
            return self  # type: ignore[return-value]

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
        return self.magnitude("1")

    def magnitude(self, units: str = "1") -> Number | npt.NDArray[np.float64]:
        raw = self.get_raw_array(units)
        if self.size == 1:
            return raw.item()

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

    def __getitem__(self, index: int) -> Scalar[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("Scalar index out of bounds")

        return cast(Scalar[SomeDim], Scalar[self.dim](np.array(self._values)[index]))


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
        x: Tensor_S[SomeDim] | Number | npt.NDArray[np.float64],
        y: Tensor_S[SomeDim] | Number | npt.NDArray[np.float64],
        z: Tensor_S[SomeDim] | Number | npt.NDArray[np.float64],
    ) -> Tensor_V3[SomeDim]:
        vec_els = cast(list[object], [x, y, z])
        if all(isinstance(el, (int, float, np.floating)) for el in vec_els):
            return cls(np.array(vec_els).reshape((3, 1)))

        if all(isinstance(el, (np.ndarray, list)) for el in vec_els):
            arrs = [np.asarray(el) for el in vec_els]
            if not all(arr.ndim == 1 for arr in arrs):
                raise ValueError("All component arrays must be one-dimensional")

            return cls(np.stack(arrs))

        if not all(isinstance(el, Tensor_S) for el in vec_els):
            raise ValueError("Components must be numbers, arrays, or scalar tensors")

        x = cast(Tensor_S[SomeDim], x)
        y = cast(Tensor_S[SomeDim], y)
        z = cast(Tensor_S[SomeDim], z)

        ensure_same_dimensions(x, y, z)
        if cls._dim.triplet() != x.dim_triplet:
            raise ValueError("Component dimensions are incompatible with vector dimension")

        x_vals = np.asarray(x.base_unit_value).reshape(-1)
        y_vals = np.asarray(y.base_unit_value).reshape(-1)
        z_vals = np.asarray(z.base_unit_value).reshape(-1)

        if not (x_vals.size == y_vals.size == z_vals.size):
            raise ValueError("Component arrays must have the same size")

        return cls(np.stack([x_vals, y_vals, z_vals]))

    @staticmethod
    def from_spherical(
        theta: Tensor_S,
        delta: Tensor_S,
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

    @staticmethod
    def from_vectors(*vectors: Vector3[SomeDim]):
        if len(vectors) == 0:
            raise ValueError("Expected at least one vector")

        if not all(isinstance(v, Tensor_V3) for v in vectors):
            raise ValueError("Expected Vector3-compatible inputs")

        vectors_ = cast(tuple[Tensor_V3[SomeDim], ...], vectors)

        try:
            ensure_same_dimensions(*vectors_)
        except RuntimeError as ex:
            raise ValueError("All vectors must share the same dimensions") from ex

        if not all(v.size == 1 for v in vectors_):
            raise ValueError("from_vectors expects scalar Vector3 values (size == 1)")

        stacked = np.concatenate([v._values for v in vectors_], axis=1)
        first = vectors_[0]
        return Vector3Array[first.dim](stacked)  # type: ignore[index]

    def dot(self, o: Tensor_V3) -> Tensor_S:
        return Tensor_S[self.dim * o.dim](np.sum(self._values * o._values, axis=0))  # type: ignore[index]

    def cross(self, o: Tensor_V3) -> Tensor_V3:
        return Tensor_V3[self.dim * o.dim](np.cross(self._values, o._values, axis=0))  # type: ignore[index]

    @property
    def size(self) -> int:
        _, s = self._values.shape
        return int(s)

    def cast(self, dim: type[SomeOtherDim]) -> Tensor_V3[SomeOtherDim]:
        if dim.triplet() == self.dim_triplet:
            return self  # type: ignore[return-value]

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
    def length(self) -> Tensor_S:
        return sqrt(self.length_squared)  # type: ignore[return-value]

    @cached_property
    def length_squared(self) -> Tensor_S:
        return self.dot(self)

    @cached_property
    def theta(self) -> Tensor_S:
        return cast(Tensor_S, atan2(self.y, self.x))

    @cached_property
    def delta(self) -> Tensor_S:
        xy = sqrt(square(self.x) + square(self.y))
        return cast(Tensor_S, atan2(self.z, xy))

    def angle(self, o: Tensor_V3[SomeDim]) -> Tensor_S:
        dot: np.ndarray = np.sum(self._values * o._values, axis=0)
        cos_angle: np.ndarray = dot / (self.length._values * o.length._values)

        angle = np.acos(cos_angle)
        angle[cos_angle >= 1] = 0
        angle[cos_angle <= -1] = np.pi

        return Tensor_S[Angle](angle)

    def normalized(self) -> Tensor_V3:
        return cast(Tensor_V3, self / self.length)

    def __getitem__(self, index: int) -> Vector3[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("Vector index out of bounds")

        return cast(
            Vector3[SomeDim],
            Vector3[self.dim](self._values[:, index:index + 1])
        )


Tensor_V3._base_tensor_class = Tensor_V3


class Tensor_M33(Tensor[SomeDim], Generic[SomeDim]):
    """3×3 matrix carrying a physical dimension."""

    def __init__(self, data: TensorData):
        data = np.asarray(data)

        if data.ndim not in (2, 3):
            raise ValueError("Wrong shape")

        rows, cols, *_ = data.shape
        if rows != 3 or cols != 3:
            raise ValueError("Wrong shape")

        super().__init__(data)

    @property
    def size(self) -> int:
        if self._values.ndim == 2:
            return 1

        *_, size = self._values.shape
        return int(size)

    @classmethod
    def from_elements(
        cls,
        a: Tensor_S[SomeDim] | Number, b: Tensor_S[SomeDim] | Number, c: Tensor_S[SomeDim] | Number,
        d: Tensor_S[SomeDim] | Number, e: Tensor_S[SomeDim] | Number, f: Tensor_S[SomeDim] | Number,
        g: Tensor_S[SomeDim] | Number, h: Tensor_S[SomeDim] | Number, i: Tensor_S[SomeDim] | Number,
    ) -> Tensor_M33[SomeDim]:
        """Create a matrix from row-major coefficients."""
        mat: TensorData
        mat_els = cast(list[object], [a, b, c, d, e, f, g, h, i])

        if all(isinstance(el, (int, float, np.floating)) for el in mat_els):
            mat = np.array(mat_els).reshape((3, 3))
        elif all(isinstance(el, Tensor_S) and el.size == 1 for el in mat_els):
            scalars = cast(list[Tensor_S[Any]], mat_els)
            ensure_same_dimensions(*scalars)
            mat = np.array([el.base_unit_value for el in scalars]).reshape((3, 3))
        else:
            raise RuntimeError("Matrix elements must be all numbers or all scalar tensors")

        return cls(mat)

    def cast(self, dim: type[SomeOtherDim]) -> Tensor_M33[SomeOtherDim]:
        """Type-cast to another dimension when coordinates are identical."""
        if dim.triplet() == self.dim_triplet:
            return self  # type: ignore[return-value]

        raise RuntimeError("Cannot cast")

    def inverse(self) -> Tensor_M33:
        """Return matrix inverse."""
        vals = np.asarray(self._values)
        if vals.ndim == 2:
            inv_vals = np.linalg.inv(vals)
        elif vals.ndim == 3:
            moved = np.moveaxis(vals, 2, 0)
            inv_moved = np.linalg.inv(moved)
            inv_vals = np.moveaxis(inv_moved, 0, 2)
        else:
            raise ValueError("Unsupported array shape for matrix inverse")

        return Tensor_M33[1 / self.dim](inv_vals)  # type: ignore[index]

    def __repr__(self) -> str:
        return f"TensorMatrix33[{self.dim.__name__}]({self._values})"

    def __getitem__(self, index: int) -> Matrix33[SomeDim]:
        if self.size == 1:
            return cast(Matrix33[SomeDim], self)

        if index < 0 or index >= self.size:
            raise KeyError("Matrix index out of bounds")

        return cast(
            Matrix33[SomeDim],
            Matrix33[self.dim](self._values[:, :, index:index + 1])
        )


Tensor_M33._base_tensor_class = Tensor_M33


class Scalar(Tensor_S[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for one-element scalar tensors."""


class ScalarArray(Tensor_S[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for many-element scalar tensors."""


class Vector3(Tensor_V3[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for three-element vector tensors."""

    O: ClassVar[Tensor_V3]
    X: ClassVar[Tensor_V3]
    Y: ClassVar[Tensor_V3]
    Z: ClassVar[Tensor_V3]


Vector3.O = Vector3[Dimless].from_components(0, 0, 0)
Vector3.X = Vector3[Dimless].from_components(1, 0, 0)
Vector3.Y = Vector3[Dimless].from_components(0, 1, 0)
Vector3.Z = Vector3[Dimless].from_components(0, 0, 1)


class Vector3Array(Tensor_V3[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for many-element vector tensors."""


class Matrix33(Tensor_M33[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for one 3×3 matrix."""

    Id: ClassVar[Tensor_M33]


Matrix33.Id = Matrix33[Dimless].from_elements(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
)


class Matrix33Array(Tensor_M33[SomeDim], Generic[SomeDim]):
    """Convenience wrapper for many 3×3 matrices."""


def abs(tensor: Tensor) -> Tensor:
    return tensor._base_tensor_class[tensor.dim](np.abs(tensor._values))  # type: ignore[index]


def square(tensor: Tensor) -> Tensor:
    return tensor._base_tensor_class[tensor.dim ** 2](np.square(tensor._values))  # type: ignore[index]


def cube(tensor: Tensor) -> Tensor:
    return tensor._base_tensor_class[tensor.dim ** 3](np.power(tensor._values, 3))  # type: ignore[index]


def sqrt(tensor: Tensor) -> Tensor:
    return tensor._base_tensor_class[tensor.dim ** Fraction(1, 2)](np.sqrt(tensor._values))  # type: ignore[index]


def cbrt(tensor: Tensor) -> Tensor:
    return tensor._base_tensor_class[tensor.dim ** Fraction(1, 3)](np.cbrt(tensor._values))  # type: ignore[index]


def cos(tensor: Tensor) -> Tensor:
    if not tensor.check(Angle):
        raise ValueError("cos() expects an angle-typed tensor")

    return tensor._base_tensor_class[Dimless](np.cos(tensor._values))  # type: ignore[index]


def sin(tensor: Tensor) -> Tensor:
    if not tensor.check(Angle):
        raise ValueError("sin() expects an angle-typed tensor")

    return tensor._base_tensor_class[Dimless](np.sin(tensor._values))  # type: ignore[index]


def tan(tensor: Tensor) -> Tensor:
    if not tensor.check(Angle):
        raise ValueError("tan() expects an angle-typed tensor")

    return tensor._base_tensor_class[Dimless](np.tan(tensor._values))  # type: ignore[index]


def acos(tensor: Tensor) -> Tensor:
    if not tensor.check(Dimless):
        raise ValueError("acos() expects a dimensionless tensor")

    return tensor._base_tensor_class[Angle](np.acos(tensor._values))  # type: ignore[index]


def asin(tensor: Tensor) -> Tensor:
    if not tensor.check(Dimless):
        raise ValueError("asin() expects a dimensionless tensor")

    return tensor._base_tensor_class[Angle](np.asin(tensor._values))  # type: ignore[index]


def atan(tensor: Tensor) -> Tensor:
    if not tensor.check(Dimless):
        raise ValueError("atan() expects a dimensionless tensor")

    return tensor._base_tensor_class[Angle](np.atan(tensor._values))  # type: ignore[index]


def atan2(y: Tensor, x: Tensor) -> Tensor:
    """Elementwise two-argument arctangent that returns an angle-typed tensor."""
    if not y.ensure_compatible_dimensions(x) or y._base_tensor_class is not x._base_tensor_class:
        raise ValueError("Incompatible dimensions or tensor types")

    vals = np.atan2(y._values, x._values)
    return y._base_tensor_class[Angle](vals)  # type: ignore[index]


def normalize_angle(angle: Tensor) -> Tensor:
    """Return the given angle in its `[0, 2π]` range."""
    return angle._base_tensor_class[Angle](angle._values % (2 * np.pi))  # type: ignore[index]


def normalize_angle_symmetric(angle: Tensor) -> Tensor:
    """Return the given angle in its `[-π, π]` range."""
    normalized = np.asarray(angle._values) % (2 * np.pi)
    normalized = np.where(normalized > np.pi, normalized - 2 * np.pi, normalized)
    return angle._base_tensor_class[Angle](normalized)  # type: ignore[index]


def interpolate(t1: Tensor[SomeDim], t2: Tensor[SomeDim], alpha: Number) -> Tensor[SomeDim]:
    """Return the linear interpolation between `t1` and `t2` with `0 <= alpha <= 1`."""
    if not t1.ensure_compatible_dimensions(t2) or t1._base_tensor_class is not t2._base_tensor_class:
        raise ValueError("Incompatible dimensions or tensor types")

    return t1._base_tensor_class[t1.dim](  # type: ignore[index]
        (1 - alpha) * t1._values + alpha * t2._values
    )


OPERATION_RESULT_TYPE: dict[
    tuple[type[Tensor], type[Tensor]],
    dict[str, type[Tensor] | None],
] = {
    (Tensor_S, Tensor_S): {
        "+": Tensor_S,
        "-": Tensor_S,
        "*": Tensor_S,
        "/": Tensor_S,
        "@": None,
        "%": Tensor_S,
    },
    (Tensor_S, Tensor_V3): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3,
        "/": None,
        "@": None,
        "%": None,
    },
    (Tensor_S, Tensor_M33): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33,
        "/": None,
        "@": None,
        "%": None,
    },
    (Tensor_V3, Tensor_S): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3,
        "/": Tensor_V3,
        "@": None,
        "%": Tensor_V3,
    },
    (Tensor_V3, Tensor_V3): {
        "+": Tensor_V3,
        "-": Tensor_V3,
        "*": Tensor_V3,
        "/": Tensor_V3,
        "@": None,
        "%": Tensor_V3,
    },
    (Tensor_V3, Tensor_M33): {
        "+": None,
        "-": None,
        "*": None,
        "/": None,
        "@": None,
        "%": None,
    },
    (Tensor_M33, Tensor_S): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33,
        "/": Tensor_M33,
        "@": None,
        "%": Tensor_M33,
    },
    (Tensor_M33, Tensor_V3): {
        "+": None,
        "-": None,
        "*": None,
        "/": None,
        "@": Tensor_V3,
        "%": None,
    },
    (Tensor_M33, Tensor_M33): {
        "+": Tensor_M33,
        "-": Tensor_M33,
        "*": Tensor_M33,
        "/": Tensor_M33,
        "@": Tensor_M33,
        "%": Tensor_M33,
    },
}

class QuantityMeta(type):
    """Metaclass exposing registered units as class attributes.

    Example:
        ``Quantity.km`` returns a ``Scalar[D.Length]`` with value ``1000``.
    """

    def __getattr__(cls, name: str) -> Scalar:
        """Resolve a unit name into its corresponding scalar quantity."""
        global __UNITS_FACTORS_DIMENSIONS
        
        try:
            dim, factor = __UNITS_FACTORS_DIMENSIONS[name]
            return Scalar[dim](np.asarray(factor))
        except KeyError:
            raise ValueError(f"No unit named '{name}'")

class Quantity(metaclass=QuantityMeta):
    """Convenience namespace for creating unit-scaled scalar values."""

    @classmethod
    def get(cls, value: str) -> Scalar:
        """Return the scalar unit associated with ``value``."""
        return cls.__getattr__(value)
    
    # DIMENSIONLESS

    dimensionless: Scalar[Dimless]
    """dimensionless (1)"""

    # ANGLES

    radian: Scalar[Angle]
    """radian"""

    turn: Scalar[Angle]
    """turns (360°)"""

    degree: Scalar[Angle]
    """degree"""

    # DISTANCES

    meter: Scalar[Length] # type: ignore[assignment]
    """meter"""

    kilo_meter: Scalar[Length]
    """kilometer"""

    radii_earth: Scalar[Length]
    """Mean radius of planet Earth (R🜨). 

    `R🜨 = 6378135 m`
    """

    radii_sun: Scalar[Length]
    """Mean radius of Sun (R☉). 

    `R☉ = 6.957e8 m`
    """

    astronomical_unit: Scalar[Length]
    """Astronomical unit (au)."""

    # DURATIONS

    second: Scalar[Time]
    """second"""

    minute: Scalar[Time]
    """minute"""

    hour: Scalar[Time]
    """hour"""

    day: Scalar[Time]
    """day"""

    month: Scalar[Time]
    """month (30 days)"""

    year: Scalar[Time]
    """year (365 days)"""

    # MASSES

    kilo_gram: Scalar[Mass]
    """kilogram"""

    gram: Scalar[Mass]
    """gram"""

    metric_ton: Scalar[Mass]
    """metric ton (1000 kg)"""