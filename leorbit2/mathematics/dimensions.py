from abc import ABC, abstractmethod
from dataclasses import dataclass
from fractions import Fraction
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

@dataclass(frozen=True)
class DimEls:
    length: Fraction | int = 0
    time: Fraction | int = 0

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
    
    def sqrt(self) -> DimEls:
        return DimEls(
            length=Fraction(self.length) / 2,
            time=Fraction(self.time) / 2
        )

class Dim:
    """Base class for dimensions"""
    _d: ClassVar[DimEls | None] = None

class D:
    """Dimensions registry."""

    class Dimless(Dim):
        """1"""
        _d = DimEls()

    class Angle(Dim):
        """rad"""
        _d = DimEls()

        rad: ClassVar[float] = 1 # base
        deg: ClassVar[float] = np.pi / 180

    class AngularAcc(Dim):
        """rad/s**2"""
        _d = DimEls(time=-2)

    class AngularJerk(Dim):
        """rad/s**2"""
        _d = DimEls(time=-3)

    class Length(Dim):
        """m"""
        _d = DimEls(length=1)

        meter: ClassVar[float] = 1 # base

        milli_meter: ClassVar[float] = 1e-3 * meter
        kilo_meter: ClassVar[float] = 1e3 * meter

    class InvLength(Dim):
        """1/m"""
        _d = DimEls(length=-1)

    class Time(Dim):
        """s"""
        _d = DimEls(time=1)

        second: ClassVar[float] = 1 # base

        milli_second: ClassVar[float] = 1e-3 * second
        minute: ClassVar[float] = 60 * second
        hour: ClassVar[float] = 60 * minute
        day: ClassVar[float] = 24 * hour
        month: ClassVar[float] = 30 * day
        year: ClassVar[float] = 365 * day

    class Velocity(Dim):
        """m/s"""
        _d = DimEls(length=1, time=-1)

        meter_per_second: ClassVar[float] = 1 # base

        kilo_meter_per_hour: ClassVar[float] = meter_per_second / 3.6

# Registry of known dimensions for lookup
_DIMENSION_REGISTRY: dict[DimEls, type[Dim]] = {
    DimEls(): D.Dimless,
    DimEls(length=1): D.Length,
    DimEls(time=1): D.Time,
    DimEls(length=1, time=-1): D.Velocity,
}


Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
SomeDimFull = TypeVar("SomeDimFull", bound=D.Length | D.Time | D.Velocity)
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
    
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)