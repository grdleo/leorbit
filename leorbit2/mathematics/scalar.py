from __future__ import annotations

from fractions import Fraction
from typing import Any, Generic, Literal, Never, TypeVar, cast, overload

import numpy as np
import numpy.typing as npt

from leorbit2.mathematics.dimensions import D, Dim, DimCoords, PowerDim, ProductDim, QuotientDim, SomeDim, SomeOtherDim
from leorbit2.mathematics.tensor import Tensor

Number = float | int | np.floating
TensorData = np.typing.NDArray[np.floating[Any]]


class TensorScalar(Tensor[SomeDim], Generic[SomeDim]):
    """Core scalar tensor type that supports one or many scalar elements."""

    @classmethod
    def dimensionalize(cls, dim_coords: DimCoords) -> type[TensorScalar]:
        return cast(type[TensorScalar], super().dimensionalize(dim_coords))

    @classmethod
    def new(cls, value: Number | TensorData | list[Number]) -> TensorScalar:
        if isinstance(value, (int, float, np.floating, np.integer)):
            return cls(value)
        if isinstance(value, list):
            return cls(np.array(value).flatten())
        if isinstance(value, np.ndarray):
            if value.ndim == 0:
                return cls(value.item())
            return cls(value.flatten())
        raise TypeError("...")

    def cast(self, dim: type[SomeOtherDim]) -> TensorScalar[SomeOtherDim]:
        if dim._d == self.dim_coords:
            return self  # type: ignore
        raise RuntimeError("Cannot cast")

    @property
    def base_unit_value(self) -> Number | np.ndarray:
        if np.asarray(self._values).ndim == 0:
            return np.float64(self._values)
        return cast(np.ndarray, np.float64(self._values))

    def magnitude(self, units: str = "1") -> np.float64 | npt.NDArray[np.float64]:
        raw = self.get_raw_array(units)
        if np.asarray(raw).ndim == 0:
            return np.float64(raw)
        return cast(npt.NDArray[np.float64], raw)

    @staticmethod
    def _result_class(a: object, b: object) -> type[TensorScalar]:
        if isinstance(a, ScalarArray) or isinstance(b, ScalarArray):
            return ScalarArray
        return Scalar

    @overload
    def __add__(self: TensorScalar[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __add__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: TensorScalar[SomeDim], o: TensorScalar[SomeDim]) -> TensorScalar[SomeDim]: ...

    def __add__(self, o: object) -> TensorScalar[Any]:
        if isinstance(o, Number):
            return cast(TensorScalar[Any], super().__add__(o))
        if not isinstance(o, TensorScalar):
            raise TypeError("...")
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        cls = self._result_class(self, o)
        return cls[self.dim](self._values + o._values)  # type: ignore
    
    @overload
    def __radd__(self: TensorScalar[D.Dimless], o: Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __radd__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    def __radd__(self, o: object) -> TensorScalar[Any]:
        return self.__add__(o)

    @overload
    def __sub__(self: TensorScalar[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __sub__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: TensorScalar[SomeDim], o: TensorScalar[SomeDim]) -> TensorScalar[SomeDim]: ...

    def __sub__(self, o: object) -> TensorScalar[Any]:
        if isinstance(o, Number):
            return cast(TensorScalar[Any], super().__sub__(o))
        if not isinstance(o, TensorScalar):
            raise TypeError("...")
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        cls = self._result_class(self, o)
        return cls[self.dim](self._values - o._values)  # type: ignore

    @overload
    def __rsub__(self: TensorScalar[D.Dimless], o: Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __rsub__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    def __rsub__(self, o: object) -> TensorScalar[Any]:
        return cast(TensorScalar[Any], super().__rsub__(o))

    @overload
    def __mul__(self: TensorScalar[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __mul__(self: TensorScalar[D.Dimless], o: TensorScalar[SomeOtherDim]) -> TensorScalar[SomeOtherDim]: ...

    @overload
    def __mul__(self: TensorScalar[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[SomeDim]: ...

    @overload
    def __mul__(self: TensorScalar[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorScalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    def __mul__(self, o: object) -> TensorScalar[Any]:
        if isinstance(o, Number):
            return cast(TensorScalar[Any], super().__mul__(o))
        if not isinstance(o, TensorScalar):
            raise TypeError("...")
        cls = self._result_class(self, o)
        return cls.dimensionalize(self.dim_coords * o.dim_coords)(self._values * o._values)

    @overload
    def __rmul__(self: TensorScalar[D.Dimless], o: Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __rmul__(self: TensorScalar[SomeDim], o: Number) -> TensorScalar[SomeDim]: ...

    def __rmul__(self, o: object) -> TensorScalar[Any]:
        return cast(TensorScalar[Any], self.__mul__(o))

    @overload
    def __truediv__(self: TensorScalar[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __truediv__(self: TensorScalar[D.Dimless], o: TensorScalar[SomeDim]) -> TensorScalar[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __truediv__(self: TensorScalar[SomeDim], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[SomeDim]: ...

    @overload
    def __truediv__(self: TensorScalar[SomeDim], o: TensorScalar[SomeOtherDim]) -> TensorScalar[QuotientDim[SomeDim, SomeOtherDim]]: ...

    def __truediv__(self, o: object) -> TensorScalar[Any]:
        if isinstance(o, Number):
            return cast(TensorScalar[Any], super().__truediv__(o))
        if not isinstance(o, TensorScalar):
            raise TypeError("...")
        cls = self._result_class(self, o)
        return cls.dimensionalize(self.dim_coords / o.dim_coords)(self._values / o._values)

    @overload
    def __rtruediv__(self: TensorScalar[D.Dimless], o: Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __rtruediv__(self: TensorScalar[SomeDim], o: Number) -> TensorScalar[QuotientDim[D.Dimless, SomeDim]]: ...

    def __rtruediv__(self, o: object) -> TensorScalar[Any]:
        return cast(TensorScalar[Any], super().__rtruediv__(o))

    @overload
    def __mod__(self: TensorScalar[D.Dimless], o: Number | TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    @overload
    def __mod__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __mod__(self: TensorScalar[SomeDim], o: TensorScalar[SomeDim]) -> TensorScalar[SomeDim]: ...

    def __mod__(self, o: object) -> TensorScalar[Any]:
        if isinstance(o, Number):
            return cast(TensorScalar[Any], super().__mod__(o))
        if not isinstance(o, TensorScalar):
            raise TypeError("...")
        if not self.ensure_compatible_dimensions(o):
            raise RuntimeError("Tensors dimensions are not compatible!")
        cls = self._result_class(self, o)
        return cls[self.dim](self._values % o._values)  # type: ignore

    @overload
    def __rmod__(self: TensorScalar[D.Dimless], o: Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __rmod__(self: TensorScalar[SomeDim], o: Number) -> Never: ...

    def __rmod__(self, o: object) -> TensorScalar[Any]:
        return cast(TensorScalar[Any], super().__rmod__(o))

    def __matmul__(self, o: Never) -> Never:
        raise RuntimeError("@ operation not defined for scalar")

    def __eq__(self, o: object) -> np.ndarray:
        return np.asarray(super().__eq__(o)._values)

    def __neq__(self, o: object) -> np.ndarray:
        return np.asarray(super().__neq__(o)._values)

    def __lt__(self, o: object) -> np.ndarray:
        return np.asarray(super().__lt__(o)._values)

    def __le__(self, o: object) -> np.ndarray:
        return np.asarray(super().__le__(o)._values)

    def __gt__(self, o: object) -> np.ndarray:
        return np.asarray(super().__gt__(o)._values)

    def __ge__(self, o: object) -> np.ndarray:
        return np.asarray(super().__ge__(o)._values)

    @overload
    def __pow__(self: TensorScalar[D.Dimless], o: Fraction | Number) -> TensorScalar[D.Dimless]: ...

    @overload
    def __pow__(self: TensorScalar[D.Dimless], o: TensorScalar[D.Dimless]) -> TensorScalar[D.Dimless]: ...

    def __pow__(self, o: object) -> Tensor[Any]:
        if isinstance(o, TensorScalar) and o.check(D.Dimless):
            o = np.asarray(o.base_unit_value).item()

        if not isinstance(o, (Fraction, int, float, np.floating)):
            raise ValueError()

        cls: type[TensorScalar] = ScalarArray if isinstance(self, ScalarArray) else Scalar
        return cls.dimensionalize(self.dim_coords ** Fraction(o))(np.power(self._values, float(o)))


class Scalar(TensorScalar[SomeDim], Generic[SomeDim]):
    """End-user convenience scalar class."""

    @classmethod
    def new(cls, value: Number | TensorData) -> Scalar:
        if isinstance(value, np.ndarray):
            value = value.item()
        return cast(Scalar, super().new(value))

    def magnitude(self, units: str = "1") -> np.float64:
        return np.float64(self.get_raw_array(units))

    @staticmethod
    def _cmp_bool(value: object) -> bool:
        return bool(np.asarray(value).item())

    def __eq__(self, o: object) -> bool:
        return self._cmp_bool(super().__eq__(o))

    def __neq__(self, o: object) -> bool:
        return self._cmp_bool(super().__neq__(o))

    def __lt__(self, o: object) -> bool:
        return self._cmp_bool(super().__lt__(o))

    def __le__(self, o: object) -> bool:
        return self._cmp_bool(super().__le__(o))

    def __gt__(self, o: object) -> bool:
        return self._cmp_bool(super().__gt__(o))

    def __ge__(self, o: object) -> bool:
        return self._cmp_bool(super().__ge__(o))


class ScalarArray(TensorScalar[SomeDim], Generic[SomeDim]):
    """End-user convenience scalar-array class."""

    @classmethod
    def new(cls, values: list[Number] | TensorData) -> ScalarArray:
        return cast(ScalarArray, super().new(np.array(values).flatten()))

    @property
    def size(self) -> int:
        s, = np.array(self._values).shape
        return int(s)

    def __getitem__(self, index: int) -> Scalar[SomeDim]:
        if index < 0 or index >= self.size:
            raise KeyError("...")
        return Scalar[self.dim].new(np.array(self._values)[index])

    def magnitude(self, units: str = "1") -> npt.NDArray[np.float64]:
        return cast(npt.NDArray[np.float64], self.get_raw_array(units))
