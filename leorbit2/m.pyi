from __future__ import annotations

from fractions import Fraction
from typing import Any, ClassVar, Generic, Literal, Never, TypeAlias, TypeVar, overload

import numpy as np
import numpy.typing as npt

Number: TypeAlias = float | int | np.floating[Any]
TensorData: TypeAlias = npt.NDArray[np.floating[Any]]


class DimCoords:
    length: Fraction
    time: Fraction
    mass: Fraction
    dimensionless: bool

    def __mul__(self, o: DimCoords) -> DimCoords: ...
    def __truediv__(self, o: DimCoords) -> DimCoords: ...
    def __pow__(self, p: Fraction | int | float) -> DimCoords: ...
    def to_dimension(self) -> type[Dim]: ...


class Dim:
    _d: ClassVar[DimCoords]


class D:
    class Dimless(Dim): ...
    class Angle(Dim): ...
    class Length(Dim): ...
    class Time(Dim): ...
    class Mass(Dim): ...


SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)


class ProductDim(Dim, Generic[SomeDim, SomeOtherDim]): ...
class QuotientDim(Dim, Generic[SomeDim, SomeOtherDim]): ...
class PowerDim(Dim, Generic[SomeDim]): ...


class Tensor(Generic[SomeDim]):
    _values: TensorData
    _base_tensor_class: ClassVar[type[Tensor[Any]]]

    @property
    def dim_coords(self) -> DimCoords: ...

    @property
    def dim(self) -> type[SomeDim]: ...

    def cast(self, dim: type[SomeOtherDim]) -> Tensor[SomeOtherDim]: ...
    def copy(self) -> Tensor[SomeDim]: ...
    def ensure_compatible_dimensions(self, o: Tensor[Any]) -> bool: ...
    def check(self, dim: type[Dim]) -> bool: ...
    def get_raw_array(self, units: str = "1") -> npt.NDArray[np.float64]: ...


class Tensor_S(Tensor[SomeDim], Generic[SomeDim]):
    @property
    def size(self) -> int: ...

    @property
    def base_unit_value(self) -> Number | npt.NDArray[np.float64]: ...

    def magnitude(self, units: str = "1") -> np.float64 | npt.NDArray[np.float64]: ...


class Tensor_V3(Tensor[SomeDim], Generic[SomeDim]):
    @classmethod
    def from_components(
        cls,
        x: Tensor_S[SomeDim] | Number,
        y: Tensor_S[SomeDim] | Number,
        z: Tensor_S[SomeDim] | Number,
    ) -> Tensor_V3[SomeDim]: ...

    @staticmethod
    def from_spherical(
        theta: Tensor_S[D.Angle],
        delta: Tensor_S[D.Angle],
        rho: Tensor_S[SomeOtherDim],
    ) -> Tensor_V3[SomeOtherDim]: ...

    def dot(self, o: Tensor_V3[SomeOtherDim]) -> Tensor_S[ProductDim[SomeDim, SomeOtherDim]]: ...
    def cross(self, o: Tensor_V3[SomeOtherDim]) -> Tensor_V3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @property
    def size(self) -> int: ...

    @property
    def x(self) -> Tensor_S[SomeDim]: ...

    @property
    def y(self) -> Tensor_S[SomeDim]: ...

    @property
    def z(self) -> Tensor_S[SomeDim]: ...

    @property
    def length(self) -> Tensor_S[SomeDim]: ...

    @property
    def theta(self) -> Tensor_S[D.Angle]: ...

    @property
    def delta(self) -> Tensor_S[D.Angle]: ...

    def angle(self, o: Tensor_V3[SomeDim]) -> Tensor_S[D.Angle]: ...
    def normalized(self) -> Tensor_V3[D.Dimless]: ...


class Tensor_M33(Tensor[SomeDim], Generic[SomeDim]):
    @classmethod
    def from_elements(
        cls,
        a: Tensor_S[SomeDim] | Number,
        b: Tensor_S[SomeDim] | Number,
        c: Tensor_S[SomeDim] | Number,
        d: Tensor_S[SomeDim] | Number,
        e: Tensor_S[SomeDim] | Number,
        f: Tensor_S[SomeDim] | Number,
        g: Tensor_S[SomeDim] | Number,
        h: Tensor_S[SomeDim] | Number,
        i: Tensor_S[SomeDim] | Number,
    ) -> Tensor_M33[SomeDim]: ...

    @property
    def det(self) -> Number: ...

    @overload
    def inverse(self: Tensor_M33[D.Dimless]) -> Tensor_M33[D.Dimless]: ...

    @overload
    def inverse(self: Tensor_M33[SomeDim]) -> Tensor_M33[QuotientDim[D.Dimless, SomeDim]]: ...


class Scalar(Tensor_S[SomeDim], Generic[SomeDim]):
    @overload
    def __add__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    @overload
    def __add__(self: Scalar[SomeDim], o: ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    @overload
    def __sub__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    @overload
    def __sub__(self: Scalar[SomeDim], o: ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    @overload
    def __mul__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __mul__(self: Scalar[D.Dimless], o: Scalar[SomeOtherDim]) -> Scalar[SomeOtherDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Number | Scalar[D.Dimless]) -> Scalar[SomeDim]: ...

    @overload
    def __mul__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Scalar[D.Dimless], o: Number | Scalar[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __truediv__(self: Scalar[D.Dimless], o: Scalar[SomeOtherDim]) -> Scalar[QuotientDim[D.Dimless, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Number | Scalar[D.Dimless]) -> Scalar[SomeDim]: ...

    @overload
    def __truediv__(self: Scalar[SomeDim], o: Scalar[SomeOtherDim]) -> Scalar[QuotientDim[SomeDim, SomeOtherDim]]: ...


class ScalarArray(Tensor_S[SomeDim], Generic[SomeDim]):
    def __getitem__(self, index: int) -> Scalar[SomeDim]: ...

    @overload
    def __add__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __add__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    @overload
    def __add__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    @overload
    def __sub__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __sub__(self: ScalarArray[SomeDim], o: Number) -> Never: ...

    @overload
    def __sub__(self: ScalarArray[SomeDim], o: Scalar[SomeDim] | ScalarArray[SomeDim]) -> ScalarArray[SomeDim]: ...

    @overload
    def __mul__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __mul__(self: ScalarArray[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[SomeOtherDim]: ...

    @overload
    def __mul__(self: ScalarArray[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[SomeDim]: ...

    @overload
    def __mul__(self: ScalarArray[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: ScalarArray[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __truediv__(self: ScalarArray[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[QuotientDim[D.Dimless, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: ScalarArray[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> ScalarArray[SomeDim]: ...

    @overload
    def __truediv__(self: ScalarArray[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> ScalarArray[QuotientDim[SomeDim, SomeOtherDim]]: ...


class Vector3(Tensor_V3[SomeDim], Generic[SomeDim]):
    O: ClassVar[Vector3[D.Dimless]]
    X: ClassVar[Vector3[D.Dimless]]
    Y: ClassVar[Vector3[D.Dimless]]
    Z: ClassVar[Vector3[D.Dimless]]
    ONE: ClassVar[Vector3[D.Dimless]]

    @overload
    @classmethod
    def from_components(
        cls,
        x: Scalar[SomeDim],
        y: Scalar[SomeDim],
        z: Scalar[SomeDim],
    ) -> Vector3[SomeDim]: ...

    @overload
    @classmethod
    def from_components(
        cls,
        x: Number,
        y: Number,
        z: Number,
    ) -> Vector3[D.Dimless]: ...

    @staticmethod
    def from_spherical(
        theta: Scalar[D.Angle],
        delta: Scalar[D.Angle],
        rho: Scalar[SomeOtherDim],
    ) -> Vector3[SomeOtherDim]: ...

    def __add__(self: Vector3[SomeDim], o: Vector3[SomeDim] | Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    def __sub__(self: Vector3[SomeDim], o: Vector3[SomeDim] | Scalar[SomeDim]) -> Vector3[SomeDim]: ...

    @overload
    def __mul__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __mul__(self: Vector3[D.Dimless], o: Scalar[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[SomeDim]: ...

    @overload
    def __mul__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Vector3[D.Dimless], o: Number | Scalar[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __truediv__(self: Vector3[D.Dimless], o: Scalar[SomeOtherDim]) -> Vector3[QuotientDim[D.Dimless, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Number | Scalar[D.Dimless]) -> Vector3[SomeDim]: ...

    @overload
    def __truediv__(self: Vector3[SomeDim], o: Scalar[SomeOtherDim]) -> Vector3[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Vector3[D.Dimless], o: Vector3[D.Dimless]) -> Scalar[D.Dimless]: ...

    @overload
    def __matmul__(self: Vector3[SomeDim], o: Vector3[SomeOtherDim]) -> Scalar[ProductDim[SomeDim, SomeOtherDim]]: ...


class Vector3Array(Tensor_V3[SomeDim], Generic[SomeDim]):
    @overload
    @classmethod
    def from_components(
        cls,
        x: ScalarArray[SomeDim],
        y: ScalarArray[SomeDim],
        z: ScalarArray[SomeDim],
    ) -> Vector3Array[SomeDim]: ...

    @overload
    @classmethod
    def from_components(
        cls,
        x: npt.NDArray[np.floating[Any]],
        y: npt.NDArray[np.floating[Any]],
        z: npt.NDArray[np.floating[Any]],
    ) -> Vector3Array[D.Dimless]: ...

    @staticmethod
    def from_spherical(
        theta: ScalarArray[D.Angle],
        delta: ScalarArray[D.Angle],
        rho: ScalarArray[SomeOtherDim],
    ) -> Vector3Array[SomeOtherDim]: ...

    def __getitem__(self, index: int) -> Vector3[SomeDim]: ...

    def __add__(self: Vector3Array[SomeDim], o: Vector3Array[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    def __sub__(self: Vector3Array[SomeDim], o: Vector3Array[SomeDim] | ScalarArray[SomeDim]) -> Vector3Array[SomeDim]: ...

    @overload
    def __mul__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __mul__(self: Vector3Array[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def __mul__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[SomeDim]: ...

    @overload
    def __mul__(self: Vector3Array[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Vector3Array[D.Dimless], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __truediv__(self: Vector3Array[D.Dimless], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[QuotientDim[D.Dimless, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Vector3Array[SomeDim], o: Number | Scalar[D.Dimless] | ScalarArray[D.Dimless]) -> Vector3Array[SomeDim]: ...

    @overload
    def __truediv__(self: Vector3Array[SomeDim], o: Scalar[SomeOtherDim] | ScalarArray[SomeOtherDim]) -> Vector3Array[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Vector3Array[D.Dimless], o: Vector3[D.Dimless]) -> ScalarArray[D.Dimless]: ...

    @overload
    def __matmul__(self: Vector3Array[SomeDim], o: Vector3[SomeOtherDim]) -> ScalarArray[ProductDim[SomeDim, SomeOtherDim]]: ...


class Matrix33(Tensor_M33[SomeDim], Generic[SomeDim]):
    @overload
    @classmethod
    def from_elements(
        cls,
        a: Scalar[SomeDim], b: Scalar[SomeDim], c: Scalar[SomeDim],
        d: Scalar[SomeDim], e: Scalar[SomeDim], f: Scalar[SomeDim],
        g: Scalar[SomeDim], h: Scalar[SomeDim], i: Scalar[SomeDim],
    ) -> Matrix33[SomeDim]: ...

    @overload
    @classmethod
    def from_elements(
        cls,
        a: Number, b: Number, c: Number,
        d: Number, e: Number, f: Number,
        g: Number, h: Number, i: Number,
    ) -> Matrix33[D.Dimless]: ...

    @overload
    def inverse(self: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def inverse(self: Matrix33[SomeDim]) -> Matrix33[QuotientDim[D.Dimless, SomeDim]]: ...

    @overload
    def __add__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __add__(self: Matrix33[SomeDim], o: Matrix33[SomeDim] | Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __sub__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __sub__(self: Matrix33[SomeDim], o: Matrix33[SomeDim] | Scalar[SomeDim]) -> Matrix33[SomeDim]: ...

    @overload
    def __mul__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __mul__(self: Matrix33[D.Dimless], o: Scalar[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[SomeDim]: ...

    @overload
    def __mul__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Matrix33[D.Dimless], o: Number | Scalar[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __truediv__(self: Matrix33[D.Dimless], o: Scalar[SomeOtherDim]) -> Matrix33[QuotientDim[D.Dimless, SomeOtherDim]]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Number | Scalar[D.Dimless]) -> Matrix33[SomeDim]: ...

    @overload
    def __truediv__(self: Matrix33[SomeDim], o: Scalar[SomeOtherDim]) -> Matrix33[QuotientDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Matrix33[D.Dimless]) -> Matrix33[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Matrix33[SomeOtherDim]) -> Matrix33[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Matrix33[SomeOtherDim]) -> Matrix33[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3[D.Dimless]) -> Vector3[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3[SomeOtherDim]) -> Vector3[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Vector3[SomeOtherDim]) -> Vector3[ProductDim[SomeDim, SomeOtherDim]]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3Array[D.Dimless]) -> Vector3Array[D.Dimless]: ...

    @overload
    def __matmul__(self: Matrix33[D.Dimless], o: Vector3Array[SomeOtherDim]) -> Vector3Array[SomeOtherDim]: ...

    @overload
    def __matmul__(self: Matrix33[SomeDim], o: Vector3Array[SomeOtherDim]) -> Vector3Array[ProductDim[SomeDim, SomeOtherDim]]: ...
